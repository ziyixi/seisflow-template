"""
Get windows described in tao's 2018 GGG paper.
"""

from window import Window, Windows_collection
import obspy
from os.path import join, basename, dirname, abspath
import configparser
import pickle
import click


def load_configure(config_fname):
    current_file_dir = dirname(abspath(__file__))
    config_path = join(current_file_dir, "..", "configuration", config_fname)
    config = configparser.ConfigParser()
    config.read(config_path)
    # load configs
    traveltime_dir = config["path"]["traveltime"]
    output_dir = config["path"]["output"]
    time_length = config["simulation"].getint("time_length")
    return traveltime_dir, time_length, output_dir


def load_taveltime(traveltime_dir):
    phases = ["S", "sS", "SS", "P", "pP",
              "sP", "PP", "3.3kmps", "4.6kmps", "ScS"]
    traveltime = {}
    for each_phase in phases:
        fpath = join(traveltime_dir, f"traveltime.{each_phase}.pkl")
        with open(fpath, "rb") as handle:
            traveltime[each_phase] = pickle.load(handle)
    return traveltime


def load_eventtime(traveltime_dir):
    fpath = join(traveltime_dir, f"extra.event_time.pkl")
    with open(fpath, "rb") as handle:
        eventtime = pickle.load(handle)
    return eventtime


def generate_windows(traveltime, event_time, time_length):
    # get all the net_sta
    used_gcmtid = list(traveltime["P"].keys())[0]
    all_net_sta = list(traveltime["P"][used_gcmtid].keys())
    # for each net_sta, get all the windows, and merge them accordingly.
    phases_zr = ["P", "pP", "sP", "PP", "S", "sS", "SS"]
    phases_t = ["ScS", "S", "sS", "SS"]
    result = {}
    for each_gcmtid in list(traveltime["P"].keys()):
        result[each_gcmtid] = {}
        for each_net_sta in all_net_sta:
            used_event_time = event_time[each_gcmtid][each_net_sta]
            result[each_gcmtid][each_net_sta] = {
                "zr": Windows_collection(),
                "t": Windows_collection(),
                "surface": Windows_collection()
            }
            # zr
            for each_phase in phases_zr:
                each_phase_traveltime = traveltime[each_phase][each_gcmtid][each_net_sta]
                if (each_phase_traveltime == None):
                    continue
                elif (each_phase_traveltime > time_length):
                    continue
                else:
                    left = used_event_time + each_phase_traveltime - 20
                    if (left < used_event_time):
                        left = used_event_time
                    right = used_event_time + each_phase_traveltime + 50
                    if (right > used_event_time + time_length):
                        right = used_event_time + time_length
                    channel = "ZR"
                    network = each_net_sta.split(".")[0]
                    station = each_net_sta.split(".")[1]
                    phases = [each_phase]
                    result[each_gcmtid][each_net_sta]["zr"].append_window(Window(
                        left=left, right=right, channel=channel, network=network, station=station, phases=phases, gcmtid=each_gcmtid))
            # t
            for each_phase in phases_t:
                each_phase_traveltime = traveltime[each_phase][each_gcmtid][each_net_sta]
                if (each_phase_traveltime == None):
                    continue
                elif (each_phase_traveltime > time_length):
                    continue
                else:
                    left = used_event_time + each_phase_traveltime - 20
                    if (left < used_event_time):
                        left = used_event_time
                    right = used_event_time + each_phase_traveltime + 50
                    if (right > used_event_time + time_length):
                        right = used_event_time + time_length
                    channel = "T"
                    network = each_net_sta.split(".")[0]
                    station = each_net_sta.split(".")[1]
                    phases = [each_phase]
                    result[each_gcmtid][each_net_sta]["t"].append_window(Window(
                        left=left, right=right, channel=channel, network=network, station=station, phases=phases, gcmtid=each_gcmtid))
            # surface
            phase_travel_time_left = traveltime["4.6kmps"][each_gcmtid][each_net_sta]
            phase_travel_time_right = traveltime["3.3kmps"][each_gcmtid][each_net_sta]
            if(phase_travel_time_right < time_length):
                left = used_event_time+phase_travel_time_left-50
                if (left < used_event_time):
                    left = used_event_time
                right = used_event_time + phase_travel_time_right + 50
                if (right > used_event_time + time_length):
                    right = used_event_time+time_length
                channel = "ZRT"
                network = each_net_sta.split(".")[0]
                station = each_net_sta.split(".")[1]
                phases = ["surface"]
                result[each_gcmtid][each_net_sta]["surface"].append_window(Window(
                    left=left, right=right, channel=channel, network=network, station=station, phases=phases, gcmtid=each_gcmtid))

    # at the moment, we have the result dict: gcmtid->net_sta->type
    # merge the windows
    for each_gcmtid in result:
        result_each_gcmtid = result[each_gcmtid]
        for each_net_sta in result_each_gcmtid:
            result_each_net_sta = result_each_gcmtid[each_net_sta]
            result_each_net_sta["zr"].merge_windows()
            result_each_net_sta["t"].merge_windows()

    return result


def save_result(result, output_dir):
    output_fname = join(output_dir, "windows.pkl")
    with open(output_fname, "wb") as handle:
        pickle.dump(result, handle)


@click.command()
@click.option('--conf', required=True, type=str, help="configuration file name in the configuration directory")
def main(conf):
    traveltime_dir, time_length, output_dir = load_configure(conf)
    traveltime = load_taveltime(traveltime_dir)
    eventtime = load_eventtime(traveltime_dir)
    result = generate_windows(traveltime, eventtime, time_length)
    save_result(result, output_dir)


if __name__ == "__main__":
    main()
