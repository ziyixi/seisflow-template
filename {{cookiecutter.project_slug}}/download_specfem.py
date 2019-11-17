import subprocess


def download_specfem():
    specfem_root = "{{cookiecutter.specfem_root}}"
    specfem_url = "{{cookiecutter.specfem_url}}"
    command = f"git clone {specfem_url} {specfem_root}"
    subprocess.call(command, shell=True)


def copy_specfem_plugins():
    specfem_root = "{{cookiecutter.specfem_root}}"
    specfem_plugin = "./specfem_plugins"
    command = f"cp -R {specfem_plugin}/* {specfem_root}"
    subprocess.call(command, shell=True)


if __name__ == "__main__":
    download_specfem()
    copy_specfem_plugins()
