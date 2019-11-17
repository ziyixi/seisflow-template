"""
Set the forward simulation structure given a name
"""
import click
import subprocess
import os
from os.path import dirname, abspath, join
import sh

# current file path
current_file_dir = dirname(abspath(__file__))


def get_paths(step_name):
    # get base, cmtfiles, ref, output, database
    base = join(current_file_dir, "..", "simulations", "step_name")
    cmtfiles = join(current_file_dir, "..", "cmts", "step_name")
    ref = join(current_file_dir, "..", "specfem_ref", "step_name")
    output = join(current_file_dir, "..", "outputs", "step_name")
    database = join(current_file_dir, "..", "databases_mpis", "step_name")
    return base, cmtfiles, ref, output, database


@click.command()
@click.option('--name', required=True, type=str, help="forward simulation name")
def main(name):
    base, cmtfiles, ref, output, database = get_paths(name)
    to_run_pyname = join(current_file_dir, "forward_structure_base.py")
    sh.python(to_run_pyname, "--base", base, "--cmtfiles",
              cmtfiles, "--ref", ref, "--output", output, "--database", database)


if __name__ == "__main__":
    main()
