import subprocess
from os.path import join
import os
cwd = os.getcwd()


def setup_GLL_model():
    if("{{cookiecutter.model_type}}" == "GLL"):
        # store the raw bin files
        gll_path = "{{cookiecutter.gll_path}}"
        new_gll_path = join(cwd, "models", "raw")
        command = f"cp -r {gll_path} {new_gll_path}"
        specfem_root = "{{cookiecutter.specfem_root}}"
        command = f"ln -s {new_gll_path} {specfem_root}/DATA/GLL"
        subprocess.call(command, shell=True)
    else:
        pass


if __name__ == "__main__":
    setup_GLL_model()
