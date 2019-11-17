import obspy
from pyasdf import ASDFDataSet
import numpy as np
import multiprocessing
from os.path import join, basename, dirname, abspath
from glob import glob
from collections import namedtuple
from obspy.taup import TauPyModel
import configparser
import click
from functools import partial
import tqdm
import pickle
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

phase_list = ["s", "S", "sS", "SS", "p", "P",
              "pP", "sP", "PP", "3.3kmps", "4.6kmps", "ScS"]
Event_pair = namedtuple('Event_pair', ['gcmtid', 'lat', 'lon', 'dep', 'time'])
model = TauPyModel(model="ak135")


def load_asdf_info(asdf_fname):
    # asdf file
    with ASDFDataSet(asdf_fname, mode="r") as asdf_file:
        lat = asdf_file.events[0].origins[0].latitude
        lon = asdf_file.events[0].origins[0].longitude
        dep = asdf_file.events[0].origins[0].depth
        time = asdf_file.events[0].origins[0].time

    return lat, lon, dep, time


def load_station_info(station_fname):
    # station file
    stations = np.loadtxt(station_fname, dtype=np.str)
    stations[:, 2:] = stations[:, 2:].astype(np.float)
    return stations


def get_event_pairs(asdf_directory):
    # assume all the asdf_directory should be like 201304210322A.preprocessed_10s_to_120s.h5
    rep_events = {}
    all_asdf_files = glob(join(asdf_directory, "*h5"))
    for each_asdf_file in all_asdf_files:
        fname = basename(each_asdf_file)
        gcmtid = fname.split(".")[0]
        rep_events[gcmtid] = each_asdf_file

    all_event_pairs = []
    for gcmtid in rep_events:
        each_rep = rep_events[gcmtid]
        lat, lon, dep, time = load_asdf_info(each_rep)
        all_event_pairs.append(Event_pair(
            gcmtid=gcmtid, lat=lat, lon=lon, dep=dep, time=time))
    return all_event_pairs


def kernel_calculate_travel_time(event_pair, stations=None):
    model = TauPyModel(model="ak135")
    # for each station, we calculate the travel time
    result = {}
    evla = event_pair.lat
    evlo = event_pair.lon
    evdp = event_pair.dep
    event_time = event_pair.time

    for row in stations:
        result_template = {
            "event_time": event_time,
            "evla": evla,
            "evlo": evlo,
            "evdp": evdp,
            "gcarc": None,
            "az": None,
            "baz": None,
            "S": None,
            "sS": None,
            "SS": None,
            "P": None,
            "pP": None,
            "sP": None,
            "PP": None,
            "3.3kmps": None,
            "4.6kmps": None,
            "ScS": None
        }
        station = row[0]
        network = row[1]
        net_sta = f"{network}.{station}"
        stla = row[2]
        stlo = row[3]
        arrivals = model.get_travel_times_geo(
            evdp, evla, evlo, stla, stlo, phase_list)

        gcarc_m, az, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
        gcarc = kilometers2degrees(gcarc_m / 1000)
        result_template["gcarc"] = gcarc
        result_template["az"] = az
        result_template["baz"] = baz

        for each_arrival in arrivals:
            name = each_arrival.name
            time = each_arrival.time
            if ((name in phase_list)):
                if (name == "p"):
                    name = "P"
                if (name == "s"):
                    name = "S"
                if (result_template[name] != None):
                    result_template[name] = time
        result[net_sta] = result_template
    return result


def load_configure(config_fname):
    # current file path
    current_file_dir = dirname(abspath(__file__))
    config_path = join(current_file_dir, "..", "configuration", config_fname)
    config = configparser.ConfigParser()
    config.read(config_path)
    # get values
    asdf_directory = config["path"]["asdf"]
    station_fname = config["path"]["station"]
    output_dir = config["path"]["output"]
    processes = config["parallel"].getint("processes")
    use_tqdm = config["ui"].getboolean("tqdm")
    return asdf_directory, station_fname, output_dir, processes, use_tqdm


@click.command()
@click.option('--conf', required=True, type=str, help="configuration file name in the configuration directory")
def main(conf):
    # load the configuration file
    asdf_directory, station_fname, output_dir, processes, use_tqdm = load_configure(
        conf)
    stations = load_station_info(station_fname)
    all_event_pairs = get_event_pairs(asdf_directory)
    # parallel
    kernel = partial(kernel_calculate_travel_time, stations=stations)
    if(use_tqdm):
        with multiprocessing.Pool(processes) as p:
            travel_result = list(
                tqdm.tqdm(p.imap(kernel, all_event_pairs), total=len(all_event_pairs)))
    else:
        with multiprocessing.Pool(processes) as p:
            travel_result = list(
                p.imap(kernel, all_event_pairs), total=len(all_event_pairs))

    # save to the output directory
    phase_list.remove("p")
    phase_list.remove("s")
    for each_phase in phase_list:
        output_result = {}
        output_fname = join(output_dir, f"traveltime.{each_phase}.pkl")
        for each_pair_result, each_pair in zip(travel_result, all_event_pairs):
            gcmtid = each_pair.gcmtid
            output_result[gcmtid] = {}
            for net_sta in each_pair_result:
                travel_time = each_pair_result[net_sta][each_phase]
                output_result[gcmtid][net_sta] = travel_time
        with open(output_fname, "wb") as handle:
            pickle.dump(output_result, handle)
    # save gcarc,az,baz
    extra_save_list = ["event_time", "gcarc",
                       "az", "baz", "evla", "evlo", "evdp"]
    for each_extra in extra_save_list:
        output_result = {}
        output_fname = join(output_dir, f"extra.{each_extra}.pkl")
        for each_pair_result, each_pair in zip(travel_result, all_event_pairs):
            gcmtid = each_pair.gcmtid
            output_result[gcmtid] = {}
            for net_sta in each_pair_result:
                extra_value = each_pair_result[net_sta][each_extra]
                output_result[gcmtid][net_sta] = extra_value
        with open(output_fname, "wb") as handle:
            pickle.dump(output_result, handle)


if __name__ == "__main__":
    main()
