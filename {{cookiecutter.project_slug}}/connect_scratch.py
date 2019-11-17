import subprocess
from os.path import join
import os
cwd = os.getcwd()


def connect_scratch():
    """
    Move some files to the scratch directory and make a soft link
    """
    dirs_in_work = ["cmts", "configuration", "data_info",
                    "misfit", "scripts", "specfem_plugins", "{{cookiecutter.specfem_root}}", "figures"]
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


if __name__ == "__main__":
    connect_scratch()
