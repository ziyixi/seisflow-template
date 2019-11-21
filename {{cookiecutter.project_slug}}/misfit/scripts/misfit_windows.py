from window import Window
import pyasdf
import numpy as np
from obspy.signal.cross_correlation import correlate, xcorr_max
import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")


class Misfit_window(Window):
    def __init__(self, parent_window):
        super().__init__(left=parent_window.left, right=parent_window.right, channel=parent_window.channel,
                         network=parent_window.network, station=parent_window.station, phases=parent_window.phases.copy(), gcmtid=parent_window.gcmtid)
        self.net_sta = f"{self.network}.{self.station}"
        self.cc = None
        self.similarity = None
        self.snr = None
        self.deltat = None
        # we should always keep only one channel
        self.component = self.channel[-1]

    def update_first_arrival(self, first_arrival_dict):
        # load from file, just an float
        self.first_arrival = first_arrival_dict[self.gcmtid][self.net_sta]

    def update_snr(self, data_asdf):
        if (self.net_sta not in data_asdf.waveforms.list()):
            return
        data_wg = data_asdf.waveforms[self.net_sta]
        data_tag = data_wg.get_waveform_tags()[0]
        if (len(data_wg[data_tag]) != 3):
            return
        data_tr = data_wg[data_tag].select(component=self.component)[0].copy()
        event_time = data_asdf.events[0].origins[0].time
        if (data_tr.stats.starttime >= event_time + self.first_arrival):
            return
        # get the noise window
        signal_st = data_tr.slice(self.left, self.right)
        # get averaged power ratio
        signal_data = signal_st.data
        # noise_avg_power = np.sum(noise_data**2) / len(noise_data)
        noise_avg_power = self.cal_noise_average_energy(
            data_tr, self.first_arrival, event_time)
        signal_avg_power = np.sum(signal_data**2)/len(signal_data)
        self.snr = 10*np.log10(signal_avg_power/noise_avg_power)

    def update_cc_deltat(self, data_asdf, sync_asdf):
        if ((self.net_sta not in data_asdf.waveforms.list()) or (self.net_sta not in sync_asdf.waveforms.list())):
            return
        # we assume the delta and the event_time are the same, but the starttime may be slightly different
        # also we have to make sure the net_sta is existing
        data_wg = data_asdf.waveforms[self.net_sta]
        data_tag = data_wg.get_waveform_tags()[0]
        if (len(data_wg[data_tag]) != 3):
            return
        sync_wg = sync_asdf.waveforms[self.net_sta]
        sync_tag = sync_wg.get_waveform_tags()[0]
        if (len(sync_wg[sync_tag]) != 3):
            return
        data_tr = data_wg[data_tag].select(component=self.component)[0].copy()
        sync_tr = sync_wg[sync_tag].select(component=self.component)[0].copy()
        # we make the starttime of sync to be the same with data
        tolerance_time = 60
        time_difference = np.abs(
            sync_tr.stats.starttime - data_tr.stats.starttime)
        if (time_difference <= data_tr.stats.delta):
            sync_tr.stats.starttime = data_tr.stats.starttime
        elif ((time_difference <= tolerance_time) and (data_tr.stats.starttime <= self.left)):
            pass
        else:
            return
        # cut to the window
        data_win_tr = data_tr.slice(self.left, self.right)
        sync_win_tr = sync_tr.slice(self.left, self.right)
        # use data as the reference, calculate cc and deltat
        cc_all = correlate(data_win_tr, sync_win_tr, None, demean=False)
        self.similarity = cc_all[len(cc_all) // 2]
        self.deltat, self.cc = xcorr_max(cc_all, abs_max=False)
        delta = data_tr.stats.delta
        self.deltat = self.deltat*delta

    def __repr__(self):
        return f"Windows(left={self.left},right={self.right},channel={self.channel},network={self.network},gcmtid={self.gcmtid},station={self.station},phases={self.phases},snr={self.snr},deltat={self.deltat},cc={self.cc},similarity={self.similarity})"

    def cal_noise_average_energy(self, data_tr, first_arrival, event_time):
        """
        Because the data may have a pulse at the beginning, we should remove these parts
        """
        noise_start = None
        if(first_arrival == None):
            # no first arrival, we don't use that trace
            return 1e9
        else:
            if(first_arrival >= 120):
                # the first part of the data may have some problem
                noise_start = 100
            elif(first_arrival >= 70):
                noise_start = 50
            else:
                noise_start = 0
        noise_win_start = event_time + noise_start
        # avoid containning the first arrival
        if(first_arrival > 10):
            noise_win_end = event_time+first_arrival-10
        else:
            noise_win_end = event_time+first_arrival
        tr_noise = data_tr.slice(noise_win_start, noise_win_end)
        noise_average_energy = np.sum(tr_noise.data**2)/len(tr_noise.data)
        return noise_average_energy
