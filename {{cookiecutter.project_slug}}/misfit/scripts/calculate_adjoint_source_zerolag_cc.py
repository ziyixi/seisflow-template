"""
Read in the misfit windows, calculate the weight and weight them.
"""
import obspy
from os.path import join, basename, dirname, abspath
import configparser
import pickle
import pyasdf
import click


def load_configure(config_fname):
    current_file_dir = dirname(abspath(__file__))
    config_path = join(current_file_dir, "..", "configuration", config_fname)
    config = configparser.ConfigParser()
    config.read(config_path)
    # load configs
    misfit_windows_dir = config["path"]["misfit_windows_dir"]
    data_asdf_body_path = config["path"]["data_asdf_body_path"]
    sync_asdf_body_path = config["path"]["sync_asdf_body_path"]
    raw_sync_asdf_body_path = config["path"]["raw_sync_asdf_body_path"]
    data_asdf_surface_path = config["path"]["data_asdf_surface_path"]
    sync_asdf_surface_path = config["path"]["sync_asdf_surface_path"]
    raw_sync_asdf_surface_path = config["path"]["raw_sync_asdf_surface_path"]
    baz_dir = config["path"]["baz_dir"]
    output_dir = config["path"]["output_dir"]
    used_gcmtid = config["setting"]["used_gcmtid"]
    consider_surface = config["setting"].getboolean("consider_surface")
    use_tqdm = config["setting"].getboolean("use_tqdm")
    use_geographical_weight = config["weight"].getboolean(
        "use_geographical_weight")
    use_category_weight = config["weight"].getboolean(
        "use_category_weight")
    snr_weight_min1 = config["weight"].getboolean(
        "snr_weight_min1")
    snr_weight_min2 = config["weight"].getboolean(
        "snr_weight_min2")
    return misfit_windows_dir, data_asdf_body_path, sync_asdf_body_path, raw_sync_asdf_body_path, \
        data_asdf_surface_path, sync_asdf_surface_path, raw_sync_asdf_surface_path, baz_dir, output_dir, \
        used_gcmtid, consider_surface, use_tqdm, \
        use_geographical_weight, use_category_weight, snr_weight_min1, snr_weight_min2


def calculate_adjoint_source_each_window(mistit_window, raw_sync_asdf_trace, sync_asdf_trace, data_asdf_trace):
    # firstly we
