import subprocess


def download_specfem():
    specfem_root = "{{cookiecutter.specfem_root}}"
    specfem_url = "{{cookiecutter.specfem_url}}"
    command = f"mkdir -p {specfem_root}"
    # command = f"git clone {specfem_url} {specfem_root}"
    subprocess.call(command, shell=True)


download_specfem()
