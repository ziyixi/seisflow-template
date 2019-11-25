"""
Read in the misfit windows, calculate the weight and weight them.
"""
import obspy
from os.path import join, basename, dirname, abspath
import configparser
import pickle
import pyasdf
import click
import numpy as np
from collections import namedtuple
from obspy.geodetics import locations2degrees
import tqdm

Weight = namedtuple(
    'Weight', ['snr', 'cc', 'deltat', 'geographical', 'category'])


def load_pickle(pickle_path):
    with open(pickle_path, "rb") as f:
        data = pickle.load(f)
    return data


def save_pickle(pickle_path, to_save):
    with open(pickle_path, "wb") as f:
        pickle.dump(to_save, f, protocol=pickle.HIGHEST_PROTOCOL)


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
    output_dir = config["path"]["output_dir"]
    stations_path = config["path"]["stations_path"]

    used_gcmtid = config["setting"]["used_gcmtid"]
    consider_surface = config["setting"].getboolean("consider_surface")

    use_geographical_weight = config["weight"].getboolean(
        "use_geographical_weight")
    use_category_weight = config["weight"].getboolean(
        "use_category_weight")
    snr_weight = config["weight"]["snr_weight"]
    cc_weight = config["weight"]["cc_weight"]
    deltat_weight = config["weight"]["deltat_weight"]
    return misfit_windows_dir, data_asdf_body_path, sync_asdf_body_path, raw_sync_asdf_body_path, \
        data_asdf_surface_path, sync_asdf_surface_path, raw_sync_asdf_surface_path, output_dir, stations_path,\
        used_gcmtid, consider_surface,  \
        use_geographical_weight, use_category_weight, snr_weight, cc_weight, deltat_weight


def calculate_adjoint_source_each_window(mistit_window, raw_sync_asdf_trace, sync_asdf_trace, data_asdf_trace, mintime, maxtime):
    # firstly we prepare the windowed trace with tapering.
    # all the input traces should be a copy
    similarity = mistit_window.similarity

    # hope they will have the same length, or we have to do some trick
    # actually they don't have the same length, we might do some trick
    windowed_data_trace = data_asdf_trace.slice(
        mistit_window.left, mistit_window.right)
    windowed_sync_trace = sync_asdf_trace.slice(
        mistit_window.left, mistit_window.right)

    true_windowed_sync_time = windowed_sync_trace.stats.starttime
    offset_window_unwindowed = round(
        (true_windowed_sync_time - sync_asdf_trace.stats.starttime) / sync_asdf_trace.stats.delta)
    offset_unwindowed_raw = round(
        (sync_asdf_trace.stats.starttime-raw_sync_asdf_trace.stats.starttime)/(raw_sync_asdf_trace.stats.delta))

    windowed_sync_trace.stats.starttime = windowed_data_trace.stats.starttime
    newright = min(windowed_sync_trace.stats.endtime,
                   windowed_data_trace.stats.endtime)
    windowed_data_trace = windowed_data_trace.slice(windowed_data_trace.stats.starttime, newright
                                                    )
    windowed_sync_trace = windowed_sync_trace.slice(windowed_data_trace.stats.starttime, newright
                                                    )

    # taper
    windowed_data_trace.taper(0.05, type="hann")
    windowed_sync_trace.taper(0.05, type="hann")
    obs_norm = np.sqrt(np.sum(windowed_data_trace.data ** 2))
    syn_norm = np.sqrt(np.sum(windowed_sync_trace.data ** 2))
    Nw = obs_norm * syn_norm
    Aw = similarity * obs_norm / syn_norm
    syn_delta = windowed_sync_trace.stats.delta
    raw_syn_sampling_rate = raw_sync_asdf_trace.stats.sampling_rate

    # just consider one direction
    adjoint_source_windowed = windowed_sync_trace.copy()
    adjoint_source_windowed.data[:] = 0.0
    obs_filt_win = windowed_data_trace.data
    syn_filt_win = windowed_sync_trace.data
    adjoint_source_windowed.data = (
        obs_filt_win - Aw * syn_filt_win) / Nw / syn_delta
    # apply taper and bandpass
    adjoint_source_windowed.taper(0.05, type="hann")
    adjoint_source_windowed.filter(
        "bandpass", freqmin=1/maxtime, freqmax=1/mintime, corners=2, zerophase=True)
    # retrive to unwindowed trace
    len_window = len(adjoint_source_windowed.data)
    adjoint_source_unwindowed = sync_asdf_trace.copy()
    adjoint_source_unwindowed.data[:] = 0.0
    adjoint_source_unwindowed.data[offset_window_unwindowed:
                                   offset_window_unwindowed+len_window] = adjoint_source_windowed.data
    # interpolate to the raw sync
    adjoint_source_unwindowed.interpolate(
        sampling_rate=raw_syn_sampling_rate)
    len_adjoint_source_unwindowed = len(adjoint_source_unwindowed.data)
    adjoint_source_final = raw_sync_asdf_trace.copy()
    adjoint_source_final.data[:] = 0.0
    adjoint_source_final.data[offset_unwindowed_raw:
                              offset_unwindowed_raw+len_adjoint_source_unwindowed] = adjoint_source_unwindowed.data
    return adjoint_source_final


def cal_cos_weight(value, value1, value2):
    result = 0
    if(value < value1):
        result = 0
    elif(value1 <= value < value2):
        result = 0.5+0.5*np.cos(np.pi*(value2-value)/(value2-value1))
    else:
        result = 1
    return result


def cal_cos_weight_deltat(value, value1, value2):
    result = 0
    if (value < value1):
        result = 1
    elif (value1 <= value < value2):
        result = 0.5 - 0.5 * \
            np.cos(np.pi * (value2 - value) / (value2 - value1))
    else:
        result = 0
    return result


def weight_and_write_adjoint_source_asdf(
        result,  output_dir, used_gcmtid, snr_weight, cc_weight, deltat_weight, stations_path, time_offset, event, adjoint_source_length, trace_delta):
    """
    Write the adjoint source trace as asdf format, the waveforms will be the unweighted result, the aux part is the weighted adjoint source.
    """
    output_asdf_path = join(output_dir, f"{used_gcmtid}.h5")
    output_asdf = pyasdf.ASDFDataSet(output_asdf_path, mode="w")
    # we should create a dict storing the weighting info.
    weighting_dict = {}
    snr1, snr2 = map(float, snr_weight.split(","))
    cc1, cc2 = map(float, cc_weight.split(","))
    deltat1, deltat2 = map(float, deltat_weight.split(","))
    # type result[net_sta][category] ->list of (misfit_window,trace)
    # it's better to only consider the result dict as misfit_windows may have some None values
    # * update snr,cc,deltat
    for net_sta in result:
        weighting_dict[net_sta] = {}
        for category in result[net_sta]:
            if(len(result[net_sta][category]) == 0):
                weighting_dict[net_sta][category] = []
                continue
            else:
                weighting_dict[net_sta][category] = []
            for each_item in result[net_sta][category]:
                misfit_window, adjoint_trace = each_item
                wsnr = cal_cos_weight(misfit_window.snr_energy, snr1, snr2)
                wcc = cal_cos_weight(misfit_window.cc, cc1, cc2)
                wdeltat = cal_cos_weight_deltat(
                    np.abs(misfit_window.deltat), deltat1, deltat2)
                # use the weight to add the adjoint source together. (not normalized)
                weighting_dict[net_sta][category].append(
                    Weight(wsnr, wcc, wdeltat, None, None))
    # * update category
    # get the number of measurements for different categories:
    count_categories = {}
    net_sta_sample = list(result.keys())[0]
    all_categories = list(result[net_sta_sample].keys())
    for each_category in all_categories:
        # count number
        count_categories[each_category] = 0
        for net_sta in result:
            if(len(result[net_sta][each_category]) != 0):
                # if exist this category for net_sta
                count_categories[each_category] += 1
    # set the weight
    # we assume that the category counting is not related to snr,cc,deltat
    for net_sta in result:
        for category in result[net_sta]:
            new_weight_list = []
            for each_weight in weighting_dict[net_sta][category]:
                new_weight_list.append(each_weight._replace(
                    category=1 / count_categories[category]))
            weighting_dict[net_sta][category] = new_weight_list
    # * update geographical
    # get the distance matrix for geographical weighting, (consider all used stations)
    # we use data_asdf_waveforms_list to construct the distance matrix
    weighting_dict = update_geographical_weighting(
        weighting_dict, list(result.keys()), stations_path)

    # * and now we can write to the asdf file
    # * result and the weighting_dict have the same structure.
    # firstly we have to get the number of points for the adjoint trace
    weight_normalize_factor = 0
    final_adjoint_source = {}
    for net_sta in result:
        # in the order of E,T,Z
        adjoint_source = np.zeros((3, adjoint_source_length))
        for category in result[net_sta]:
            for each_item, each_weight in zip(result[net_sta][category], weighting_dict[net_sta][category]):
                each_misfit_window, each_adjoint_trace = each_item
                wcc = each_weight.cc
                wdeltat = each_weight.deltat
                wsnr = each_weight.snr
                wcategory = each_weight.category
                wgeographical = each_weight.geographical
                try:
                    weight_normalize_factor += wcc * wdeltat * wsnr * wcategory * wgeographical
                except:
                    print(each_misfit_window)
                    print(each_adjoint_trace)
                    print(each_weight)
                    exit()
                if (each_misfit_window.component == "Z"):
                    adjoint_source[2, :] += each_adjoint_trace.data * \
                        wcc * wdeltat * wsnr * wcategory * wgeographical
                elif (each_misfit_window.component == "R"):
                    r_adjoint_source = each_adjoint_trace.data * \
                        wcc * wdeltat * wsnr * wcategory * wgeographical
                    theta = (each_misfit_window.baz - 180) % 360
                    adjoint_source[0, :] += r_adjoint_source * \
                        np.sin(np.deg2rad(theta))
                    adjoint_source[1, :] += r_adjoint_source * \
                        np.cos(np.deg2rad(theta))
                elif (each_misfit_window.component == "T"):
                    t_adjoint_source = each_adjoint_trace.data*wcc * \
                        wdeltat * wsnr * wcategory * wgeographical
                    theta = (each_misfit_window.baz - 90) % 360
                    adjoint_source[0, :] += t_adjoint_source * \
                        np.sin(np.deg2rad(theta))
                    adjoint_source[1, :] += t_adjoint_source * \
                        np.cos(np.deg2rad(theta))
        final_adjoint_source[net_sta] = adjoint_source
    # normalize the adjoint sources
    for net_sta in final_adjoint_source:
        final_adjoint_source[net_sta] /= weight_normalize_factor
    # write each net_sta
    components = ["MXE", "MXN", "MXZ"]
    for net_sta in final_adjoint_source:
        for index_component in range(3):
            component = components[index_component]
            specfem_adj_source = np.empty((adjoint_source_length, 2))
            specfem_adj_source[:, 0] = np.linspace(
                0, (adjoint_source_length - 1) * trace_delta, adjoint_source_length)
            specfem_adj_source[:, 0] -= time_offset
            specfem_adj_source[:,
                               1] = final_adjoint_source[net_sta][index_component, :]
            tag = net_sta.replace(".", "_") + "_" + components[index_component]
            output_asdf.add_auxiliary_data(
                data=specfem_adj_source, data_type="AdjointSources", path=tag, parameters={})
    del output_asdf
    # save weighting pkl
    save_pickle(join(output_dir, f"weight.{used_gcmtid}.pkl"), weighting_dict)


def write_stations_adjoint(data_asdf_body_waveforms_list, stations_path, output_dir, used_gcmtid):
    stations_loaded = np.loadtxt(stations_path, dtype=np.str)
    output_path = join(output_dir, f"{used_gcmtid}.stations")
    with open(output_path, "w") as f:
        for row in stations_loaded:
            net_sta = f"{row[1]}.{row[0]}"
            if (net_sta in data_asdf_body_waveforms_list):
                f.write(" ".join(row)+"\n")


def update_geographical_weighting(weighting_dict, data_asdf_waveforms_list, stations_path):
    # load stations info
    stations_loc_list = []
    stations_loaded = np.loadtxt(stations_path, dtype=np.str)
    stations_loc_mapper = {}
    for row in stations_loaded:
        net = row[1]
        sta = row[0]
        net_sta = f"{net}.{sta}"
        # lon,lat
        stations_loc_mapper[net_sta] = (float(row[3]), float(row[2]))
    # build up the lon,lat matrix in the order of data_asdf_waveforms_list
    matrix_size = len(data_asdf_waveforms_list)
    lat1 = np.zeros((matrix_size, matrix_size))
    lon1 = np.zeros((matrix_size, matrix_size))
    lat2 = np.zeros((matrix_size, matrix_size))
    lon2 = np.zeros((matrix_size, matrix_size))
    for irow in range(matrix_size):
        for icolumn in range(matrix_size):
            lat1[irow, icolumn] = stations_loc_mapper[data_asdf_waveforms_list[irow]][1]
            lon1[irow, icolumn] = stations_loc_mapper[data_asdf_waveforms_list[irow]][0]
            lat2[irow, icolumn] = stations_loc_mapper[data_asdf_waveforms_list[icolumn]][1]
            lon2[irow, icolumn] = stations_loc_mapper[data_asdf_waveforms_list[icolumn]][0]
    distance_matrix = locations2degrees(lat1, lon1, lat2, lon2)
    # find Delta0 and get geographical weighting for each net_sta
    dref = np.arange(0.5, 8.5, 0.5)
    ratio_list = []
    for iref in range(len(dref)):
        wt = np.zeros(matrix_size)
        for i in range(matrix_size):
            wt[i] = 1.0/np.sum(np.exp(-(distance_matrix[i]/dref[iref])**2))
        ratio = np.max(wt)/np.min(wt)
        ratio_list.append(ratio)
    # we should find the nearest 1/3 max ratio
    ratio_list = np.array(ratio_list)
    max_ratio = np.max(ratio_list)
    candiate_ratio = max_ratio/3
    candiate_id = np.argmin(
        np.abs(ratio_list[:len(ratio_list)//2]-candiate_ratio))
    used_ref = dref[candiate_id]
    used_wt_dict = {}
    for i in range(matrix_size):
        net_sta = data_asdf_waveforms_list[i]
        used_wt_dict[net_sta] = 1.0 / \
            np.sum(np.exp(-(distance_matrix[i]/used_ref)**2))
    # update the weight
    for net_sta in weighting_dict:
        for category in weighting_dict[net_sta]:
            new_weight_list = []
            for each_weight in weighting_dict[net_sta][category]:
                new_weight_list.append(each_weight._replace(
                    geographical=used_wt_dict[net_sta]))
            weighting_dict[net_sta][category] = new_weight_list
    return weighting_dict


def run(misfit_windows_dir, data_asdf_body_path, sync_asdf_body_path, raw_sync_asdf_body_path,
        data_asdf_surface_path, sync_asdf_surface_path, raw_sync_asdf_surface_path, output_dir, stations_path,
        used_gcmtid, consider_surface,
        use_geographical_weight, use_category_weight, snr_weight, cc_weight, deltat_weight
        ):
    # * get some info
    # firstly we load the pickle file of the misfit_windows
    misfit_windows_path = join(misfit_windows_dir, f"{used_gcmtid}.pkl")
    misfit_windows = load_pickle(misfit_windows_path)
    # we can load the asdf files (remember to delete them later)
    if (consider_surface):
        data_asdf_surface = pyasdf.ASDFDataSet(
            data_asdf_surface_path, mode="r")
        sync_asdf_surface = pyasdf.ASDFDataSet(
            sync_asdf_surface_path, mode="r")
        raw_sync_asdf_surface = pyasdf.ASDFDataSet(
            raw_sync_asdf_surface_path, mode="r")
    data_asdf_body = pyasdf.ASDFDataSet(data_asdf_body_path, mode="r")
    sync_asdf_body = pyasdf.ASDFDataSet(sync_asdf_body_path, mode="r")
    raw_sync_asdf_body = pyasdf.ASDFDataSet(raw_sync_asdf_body_path, mode="r")
    # prepare some info
    data_asdf_body_waveforms_list = data_asdf_body.waveforms.list()
    sync_asdf_body_waveforms_list = sync_asdf_body.waveforms.list()
    data_asdf_body_waveforms_list = list(
        set(data_asdf_body_waveforms_list) & set(sync_asdf_body_waveforms_list))
    if(consider_surface):
        data_asdf_surface_waveforms_list = data_asdf_surface.waveforms.list()
        sync_asdf_surface_waveforms_list = sync_asdf_surface.waveforms.list()
    # * calculate the adjoint source
    # loop through all the windows, calculate the adjoint source
    # the result contains the traces for the adjoint sources
    result = {}
    adjoint_source_length = None
    trace_delta = None
    time_offset = None
    for net_sta in tqdm.tqdm(data_asdf_body_waveforms_list):
        result[net_sta] = {}
        for category in misfit_windows[net_sta]:
            result[net_sta][category] = []
            windows = misfit_windows[net_sta][category].windows
            if("surface" in category):
                data_asdf = data_asdf_surface
                sync_asdf = sync_asdf_surface
                raw_sync_asdf = raw_sync_asdf_surface
            else:
                data_asdf = data_asdf_body
                sync_asdf = sync_asdf_body
                raw_sync_asdf = raw_sync_asdf_body

            data_wg = data_asdf.waveforms[net_sta]
            sync_wg = sync_asdf.waveforms[net_sta]
            raw_sync_wg = raw_sync_asdf.waveforms[net_sta]
            data_tag = data_wg.get_waveform_tags()[0]
            sync_tag = sync_wg.get_waveform_tags()[0]
            raw_sync_tag = raw_sync_wg.get_waveform_tags()[0]

            # get mintime and maxtime
            tag_spliter = data_tag.split("_")
            mintime_str = tag_spliter[1][:-1]
            maxtime_str = tag_spliter[3][:-1]
            mintime = float(mintime_str)
            maxtime = float(maxtime_str)

            # loop through all the windows
            for each_window in windows:
                component = each_window.component
                sync_asdf_trace = sync_wg[sync_tag].select(
                    component=component)[0]
                raw_sync_asdf_trace = raw_sync_wg[raw_sync_tag].select(component=component)[
                    0]
                data_asdf_trace = data_wg[data_tag].select(
                    component=component)[0]
                time_offset = raw_sync_asdf_body.events[0].preferred_origin(
                ).time-raw_sync_asdf_trace.stats.starttime
                adjoint_source_trace = calculate_adjoint_source_each_window(
                    each_window, raw_sync_asdf_trace, sync_asdf_trace, data_asdf_trace, mintime, maxtime)
                adjoint_source_length = len(adjoint_source_trace.data)
                trace_delta = adjoint_source_trace.stats.delta
                result[net_sta][category].append(
                    (each_window, adjoint_source_trace))
    event = data_asdf_body.events
    weight_and_write_adjoint_source_asdf(
        result,  output_dir, used_gcmtid, snr_weight, cc_weight, deltat_weight, stations_path, time_offset, event, adjoint_source_length, trace_delta)
    write_stations_adjoint(data_asdf_body_waveforms_list,
                           stations_path, output_dir, used_gcmtid)

    if(consider_surface):
        del data_asdf_surface
        del sync_asdf_surface
        del raw_sync_asdf_surface
    del data_asdf_body
    del sync_asdf_body
    del raw_sync_asdf_body


@click.command()
@click.option('--conf', required=True, type=str, help="configuration file name in the configuration directory")
def main(conf):
    config_path = join("..", "configuration", conf)
    # load configuration
    misfit_windows_dir, data_asdf_body_path, sync_asdf_body_path, raw_sync_asdf_body_path, \
        data_asdf_surface_path, sync_asdf_surface_path, raw_sync_asdf_surface_path, output_dir, stations_path,\
        used_gcmtid, consider_surface,  \
        use_geographical_weight, use_category_weight, snr_weight, cc_weight, deltat_weight = load_configure(
            config_path)
    run(misfit_windows_dir, data_asdf_body_path, sync_asdf_body_path, raw_sync_asdf_body_path,
        data_asdf_surface_path, sync_asdf_surface_path, raw_sync_asdf_surface_path, output_dir, stations_path,
        used_gcmtid, consider_surface,
        use_geographical_weight, use_category_weight, snr_weight, cc_weight, deltat_weight
        )


if __name__ == "__main__":
    main()
