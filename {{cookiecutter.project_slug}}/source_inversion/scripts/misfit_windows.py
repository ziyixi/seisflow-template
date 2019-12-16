from window import Window
import pyasdf
import numpy as np
from obspy.signal.cross_correlation import correlate, xcorr_max
import sys
import warnings
import configparser

if not sys.warnoptions:
    warnings.simplefilter("ignore")


class Misfit_window(Window):
    def __init__(self, parent_window):
        super().__init__(left=parent_window.left, right=parent_window.right, channel=parent_window.channel,
                         network=parent_window.network, station=parent_window.station, phases=parent_window.phases.copy(), gcmtid=parent_window.gcmtid)
        self.net_sta = f"{self.network}.{self.station}"
        self.cc = None
        self.similarity = None
        self.snr_energy = None
        self.snr_amp = None
        self.deltat = None
        # we should always keep only one channel
        self.component = self.channel[-1]

    def update_first_arrival_baz(self, first_arrival_dict, baz_dict):
        # load from file, just an float
        self.first_arrival = first_arrival_dict[self.gcmtid][self.net_sta]
        self.baz = baz_dict[self.gcmtid][self.net_sta]

    def update_cc_deltat(self, data_asdf, sync_dict):
        if (self.net_sta not in data_asdf.waveforms.list()):
            return
        # we assume the delta and the event_time are the same, but the starttime may be slightly different
        # also we have to make sure the net_sta is existing
        data_wg = data_asdf.waveforms[self.net_sta]
        data_tag = data_wg.get_waveform_tags()[0]
        data_tr = data_wg[data_tag].select(component=self.component)[0].copy()
        sync_tr = sync_dict[self.net_sta].select(component=self.component)[0].copy()
        # we make the starttime of sync to be the same with data
        tolerance_time = 60
        time_difference = np.abs(
            sync_tr.stats.starttime - data_tr.stats.starttime)
        if (time_difference <= data_tr.stats.delta):
            sync_tr.stats.starttime = data_tr.stats.starttime
        elif ((time_difference <= tolerance_time) and (data_tr.stats.starttime <= self.left)):
            # ! need to fix here later
            sync_tr.trim(data_tr.stats.starttime, sync_tr.stats.endtime)
            sync_tr.stats.starttime = data_tr.stats.starttime
        else:
            return
        # cut to the window
        data_win_tr = data_tr.slice(self.left, self.right)
        data_win_tr.taper(0.05, type="hann")
        sync_win_tr = sync_tr.slice(self.left, self.right)
        sync_win_tr.taper(0.05, type="hann")
        # use data as the reference, calculate cc and deltat
        cc_all = correlate(data_win_tr, sync_win_tr, None,
                           demean=False, normalize="naive")
        self.similarity = cc_all[len(cc_all) // 2]
        self.deltat, self.cc = xcorr_max(cc_all, abs_max=False)
        delta = data_tr.stats.delta
        self.deltat = self.deltat * delta

    def __repr__(self):
        return f"Windows(left={self.left},right={self.right},channel={self.channel},network={self.network},gcmtid={self.gcmtid},station={self.station},phases={self.phases},snr_energy={self.snr_energy},snr_amp={self.snr_amp},deltat={self.deltat},cc={self.cc},similarity={self.similarity})"
