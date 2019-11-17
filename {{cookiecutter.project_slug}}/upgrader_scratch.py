import subprocess
from os.path import join, isdir
import os
cwd = os.getcwd()


def upgrader_scratch():
    """
    handle the soft link issue after upgrading
    """
    dirs_in_work = ["cmts", "configuration", "data_info",
                    "misfit", "scripts", "specfem_plugins", "{{cookiecutter.specfem_root}}", "figures"]
    dirs_in_scratch = ["databases_mpis", "datas", "outputs",
                       "simulation_test", "simulations", "specfem_ref", "syncs", "models"]

    scratch_dir = "{{cookiecutter.scratch_path}}"
    # copy files
    for each_scratch_dir in dirs_in_scratch:
        if(isdir(each_scratch_dir)):
            work_dir = join(cwd, each_scratch_dir, ".")
            scratch_dir = join(
                "{{cookiecutter.scratch_path}}", each_scratch_dir)
            command = f"cp -R {work_dir} {scratch_dir}"
            subprocess.call(command, shell=True)
            command = f"rm -rf {work_dir}"
            subprocess.call(command, shell=True)
            command = f"ln -s {scratch_dir} {work_dir}"
            subprocess.call(command, shell=True)


if __name__ == "__main__":
    upgrader_scratch()
