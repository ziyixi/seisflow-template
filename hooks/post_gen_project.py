import subprocess
from os.path import join
import os
cwd = os.getcwd()


def copy_specfem_plugins():
    specfem_root = "{{cookiecutter.specfem_root}}"
    specfem_plugin = "./specfem_plugins"
    command = f"cp -R {specfem_plugin}/* {specfem_root}"
    subprocess.call(command, shell=True)


def setup_GLL_model():
    if("{{cookiecutter.model_type}}" == "GLL"):
        # store the raw bin files
        gll_path = "{{cookiecutter.gll_path}}"
        new_gll_path = "./models/raw"
        command = f"cp -r {gll_path} {new_gll_path}"
        specfem_root = "{{cookiecutter.specfem_root}}"
        command = f"ln -s {new_gll_path} {specfem_root}/DATA/GLL"
        subprocess.call(command, shell=True)
    else:
        pass


def connect_scratch():
    """
    Move some files to the scratch directory and make a soft link
    """
    dirs_in_work = ["cmts", "configuration", "data_info",
                    "misfit", "scripts", "specfem_plugins"]
    dirs_in_scratch = ["databases_mpis", "datas", "outputs",
                       "simulation_test", "simulations", "specfem_ref", "syncs", "models"]

    scratch_dir = "{{cookiecutter.scratch_path}}"
    command = f"mkdir -p {scratch_dir}"
    subprocess.call(command, shell=True)
    # since the scratch will keep the directory structure, it's safe to move the directories.
    for each_scratch_dir in dirs_in_scratch:
        from_dir = join(cwd, each_scratch_dir)
        to_dir = join("{{cookiecutter.scratch_path}}", each_scratch_dir)
        command = f"mv {from_dir} {to_dir}"
        subprocess.call(command, shell=True)
        command = f"ln -s {to_dir} {from_dir}"
        subprocess.call(command, shell=True)
    for each_work_dir in dirs_in_work:
        work_path = join(cwd, each_work_dir)
        scratch_path = join("{{cookiecutter.scratch_path}}", each_work_dir)
        command = f"ln -s {work_path} {scratch_path}"
        subprocess.call(command, shell=True)


copy_specfem_plugins()
setup_GLL_model()
connect_scratch()
