import pyasdf
from os.path import join, basename, dirname, abspath
import obspy
import numpy as np 
import click
import configparser

def load_configure(config_fname):
    current_file_dir = dirname(abspath(__file__))
    config_path = join(current_file_dir, "..", "configuration", config_fname)
    config = configparser.ConfigParser()
    config.read(config_path)
    # load configs
    first_asdf_path=config["path"]["first_asdf_path"]
    second_asdf_path=config["path"]["second_asdf_path"]
    out_asdf_path=config["path"]["out_asdf_path"]
    return first_asdf_path,second_asdf_path,out_asdf_path


def run(first_asdf_path,second_asdf_path,out_asdf_path):
    """
    First,second is defined in tao's paper.
    """
    first_asdf=pyasdf.ASDFDataSet(first_asdf_path,mode="r")
    second_asdf=pyasdf.ASDFDataSet(second_asdf_path,mode="r")
    out_asdf=pyasdf.ASDFDataSet(out_asdf_path,mode="w")
    out_asdf.add_quakeml(first_asdf.events)
    out_event=first_asdf.events[0]

    # we assume the two raw asdf sync has the same structure
    all_net_sta_list=first_asdf.waveforms.list()
    for each_net_sta in all_net_sta_list:
        first_st=first_asdf.waveforms[each_net_sta].synthetic
        second_st=second_asdf.waveforms[each_net_sta].synthetic
        out_st=obspy.Stream()
        for first_tr,second_tr in zip(first_st,second_st):
            out_tr=first_tr.copy()
            out_tr.data=(first_tr.data-second_tr.data)
            out_st+=out_tr
        out_asdf.add_waveforms(out_st,tag="perturbed",event_id=out_event)


    del first_asdf
    del second_asdf
    del out_asdf

@click.command()
@click.option('--conf', required=True, type=str, help="configuration file name in the configuration directory")
def main(conf):
    first_asdf_path,second_asdf_path,out_asdf_path=load_configure(conf)
    run(first_asdf_path,second_asdf_path,out_asdf_path)

if __name__ == "__main__":
    main()