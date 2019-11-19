from windows import Window
import pyasdf
import numpy as np
from obspy.signal.cross_correlation import correlate, xcorr_max


class Misfit_window(Window):
    def __init__(self, parent_window):
        super().__init__(left=parent_window.left, right=parent_window.right, channel=parent_window.channel,
                         network=parent_window.network, station=parent_window.station, phases=parent_window.phases, gcmtid=parent_window.gcmtid)
        self.net_sta = f"{self.network}.{self.station}"
        self.cc = None
        self.snr = None
        self.deltat = None
        # we should always keep only one channel
        self.component = self.channel[-1]

    def update_first_arrival(self, first_arrival_dict):
        # load from file, just an float
        self.first_arrival = first_arrival_dict[self.gcmtid][self.net_sta]

    def update_snr(self, data_asdf):
        data_wg = data_asdf.waveforms[self.net_sta]
        data_tag = data_wg.get_waveform_tags()[0]
        data_tr = data_wg[data_tag].select(component=self.component).copy()
        # get the noise window
        event_time = data_asdf.events[0].origins[0].time
        noise_st = data_tr.slice(event_time, event_time+self.first_arrival)
        signal_st = data_tr.slice(self.left, self.right)
        # get averaged power ratio
        noise_data = noise_st.data
        signal_data = signal_st.data
        noise_avg_power = (noise_data ** 2) / len(noise_data)
        signal_avg_power = (signal_data**2)/len(signal_data)
        self.snr = 10*np.log10(signal_avg_power/noise_avg_power)

    def update_cc_deltat(self, data_asdf, sync_asdf):
        # we assume the delta and the event_time are the same, but the starttime may be slightly different
        # also we have to make sure the net_sta is existing
        data_wg = data_asdf.waveforms[self.net_sta]
        data_tag = data_wg.get_waveform_tags()[0]
        sync_wg = sync_asdf.waveforms[self.net_sta]
        sync_tag = sync_wg.get_waveform_tags()[0]
        data_tr = data_wg[data_tag].select(component=self.component).copy()
        sync_tr = sync_wg[sync_tag].select(component=self.component).copy()
        # we make the starttime of sync to be the same with data
        sync_tr.stats.starttime = data_tr.stats.starttime
        # cut to the window
        data_win_tr = data_tr.slice(self.left, self.right)
        sync_win_tr = sync_tr.slice(self.left, self.right)
        # use data as the reference, calculate cc and deltat
        cc_all = correlate(data_win_tr, sync_win_tr, None, demean=False)
        self.cc = cc_all[len(cc_all) // 2]
        self.deltat, _ = xcorr_max(cc_all)
