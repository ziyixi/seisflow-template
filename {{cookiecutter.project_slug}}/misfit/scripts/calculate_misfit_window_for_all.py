"""
Calculate the misfit window for all the events.
"""
import configparser
from get_tao_2018_ggg_misfitwin_each_event import run
from os.path import join, basename, dirname, abspath
import tqdm
import click
import subprocess
import pickle
import multiprocessing


def load_pickle(pickle_path):
    with open(pickle_path, "rb") as f:
        data = pickle.load(f)
    return data


def load_configure(config_fname):
    current_file_dir = dirname(abspath(__file__))
    config_path = join(current_file_dir, "..", "configuration", config_fname)
    config = configparser.ConfigParser()
    config.read(config_path)
    # load configuration
    windows_dir = config["path"]["windows_dir"]
    first_arrival_dir = config["path"]["first_arrival_dir"]
    baz_dir = config["path"]["baz_dir"]
    evdp_dir = config["path"]["evdp_dir"]
    data_asdf_dir = config["path"]["data_asdf_dir"]
    sync_asdf_dir = config["path"]["sync_asdf_dir"]
    output_dir = config["path"]["output_dir"]

    considered_freq_seconds = config["setting"]["considered_freq_seconds"]
    max_time_seconds = config["setting"]["max_time_seconds"]
    use_tqdm = config["setting"].getboolean("use_tqdm")
    processes = config["setting"]["processes"]

    return windows_dir, first_arrival_dir, baz_dir, evdp_dir, data_asdf_dir, sync_asdf_dir, \
        output_dir, considered_freq_seconds, use_tqdm, max_time_seconds, processes


def wrap_single_event(windows_dir, first_arrival_dir, baz_dir, data_asdf_dir, sync_asdf_dir,
                      output_dir, considered_freq_seconds, max_time_seconds, all_gcmtids):
    # mkdir for each considered_freq_seconds
    result = []
    considered_freq_seconds = considered_freq_seconds.split(",")
    for each_considered_freq_second in considered_freq_seconds:
        output_dir_freq_second = join(output_dir, each_considered_freq_second)
        command = f"mkdir -p {output_dir_freq_second}"
        subprocess.call(command, shell=True)
        for used_gcmtid in all_gcmtids:
            data_asdf_body_path = join(
                data_asdf_dir, f"{used_gcmtid}.preprocessed_{each_considered_freq_second}s_to_{max_time_seconds}s.h5")
            sync_asdf_body_path = join(
                sync_asdf_dir, f"{used_gcmtid}.preprocessed_{each_considered_freq_second}s_to_{max_time_seconds}s.h5")
            data_asdf_surface_path = join(
                data_asdf_dir, f"{used_gcmtid}.preprocessed_{each_considered_freq_second}s_to_{max_time_seconds}s.h5")
            sync_asdf_surface_path = join(
                sync_asdf_dir, f"{used_gcmtid}.preprocessed_{each_considered_freq_second}s_to_{max_time_seconds}s.h5")

            result.append((windows_dir, first_arrival_dir, baz_dir, data_asdf_body_path, sync_asdf_body_path,
                           data_asdf_surface_path, sync_asdf_surface_path, output_dir_freq_second, used_gcmtid, True, False))
    return result


def kernel(info):
    run(*info)


@click.command()
@click.option('--conf', required=True, type=str, help="configuration file name in the configuration directory")
def main(conf):
    windows_dir, first_arrival_dir, baz_dir, evdp_dir, data_asdf_dir, sync_asdf_dir, \
        output_dir, considered_freq_seconds, use_tqdm, max_time_seconds, processes = load_configure(
            conf)
    # get all gcmtids, load evdp file
    all_evdps = load_pickle(join(evdp_dir, "extra.evdp.pkl"))
    all_gcmtids = list(all_evdps.keys())
    single_run_iters = wrap_single_event(windows_dir, first_arrival_dir, baz_dir, data_asdf_dir, sync_asdf_dir,
                                         output_dir, considered_freq_seconds, max_time_seconds, all_gcmtids)
    with multiprocessing.Pool(processes) as p:
        unused_list = list(
            tqdm.tqdm(p.imap(kernel, single_run_iters), total=len(single_run_iters)))
