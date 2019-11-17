import argparse
import sys
from glob import glob
from os.path import join, dirname, abspath

from slurmpy import Slurm
import configparser

# N_total = 30
# N_each = 3
# N_iter = 10
# nproc = 441


def get_config(configure_fname):
    config = configparser.ConfigParser()
    config.read(configure_fname)
    # get values
    N_total = config["task"].getint("N_total")
    N_each = config["task"].getint("N_each")
    N_iter = config["task"].getint("N_iter")
    nproc = config["task"].getint("nproc")
    # N_node, ntasks, partition, time, account
    N_node = config["slurm"].getint("N_node")
    ntasks = config["slurm"].getint("ntasks")
    partition = config["slurm"]["partition"]
    time = config["slurm"]["time"]
    account = config["slurm"]["account"]
    return N_total, N_each, N_iter, nproc, N_node, ntasks, partition, time, account


def get_args(args=None):
    parser = argparse.ArgumentParser(
        description='A python script to submit jobs in one sbatch job',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--base',
                        help='the directory to place all the specefem directories',
                        required=True)
    parser.add_argument('--configure',
                        help='the configuration file name in the configuration directory',
                        required=True)
    results = parser.parse_args(args)
    current_file_dir = dirname(abspath(__file__))
    return results.base, join(current_file_dir, "..", "configuration", results.configure)


def get_dirs(base):
    return glob(join(base, "*"))


def get_scripts(thedirs, N_total, N_each, N_iter, nproc, N_node, ntasks, partition, time, account):
    result = "date; "
    result += "module load boost/1.68; "
    result += "module load phdf5/1.8.16; "
    # for xmeshfem3D
    result += f"echo 'start xmeshfem3D'; "
    for iiter in range(N_iter):
        result += f"echo 'start iteration {iiter}'; "
        for ieach in range(N_each):
            # ievent
            ievent = iiter*N_each+ieach
            ievent_dir = thedirs[ievent]
            # cd
            result += f"cd {ievent_dir}; "
            # if N_each-1
            if(ieach == N_each-1):
                inc = ieach*nproc
                result += f"ibrun -n {nproc} -o {inc} ./bin/xmeshfem3D; "
            else:
                inc = ieach*nproc
                result += f"ibrun -n {nproc} -o {inc} ./bin/xmeshfem3D & "
        result += f"wait; "
        result += f"echo 'end iteration {iiter}'; "
        result += f"date; "

    # for xspecfem3D
    result += f"echo 'start xspecfem3D'; "
    for iiter in range(N_iter):
        result += f"echo 'start iteration {iiter}'; "
        for ieach in range(N_each):
            # ievent
            ievent = iiter*N_each+ieach
            ievent_dir = thedirs[ievent]
            # cd
            result += f"cd {ievent_dir}; "
            # if N_each-1
            if(ieach == N_each-1):
                inc = ieach*nproc
                result += f"ibrun -n {nproc} -o {inc} ./bin/xspecfem3D; "
            else:
                inc = ieach*nproc
                result += f"ibrun -n {nproc} -o {inc} ./bin/xspecfem3D & "
        result += f"wait; "
        result += f"echo 'end iteration {iiter}'; "
        result += f"date; "

    return result


def submit_job(thecommand):
    s = Slurm("sync", {"nodes": N_node, "ntasks": ntasks,
                       "partition": partition, "time": time, "account": account})
    s.run(thecommand)


if __name__ == "__main__":
    base, configure_fname = get_args(sys.argv[1:])
    N_total, N_each, N_iter, nproc, N_node, ntasks, partition, time, account = get_config(
        configure_fname)
    thedirs = get_dirs(base)
    thecommand = get_scripts(thedirs, N_total, N_each, N_iter,
                             nproc, N_node, ntasks, partition, time, account)
    submit_job(thecommand)
