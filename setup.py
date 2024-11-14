import subprocess
import sys

from setuptools import os, setup

from glob import glob
from shutil import move

# make clean every time slows things but is needed because users will mostly install after changing networks
# and a make never seems to full update the new network
# print(f"STAGE: {sys.argv[1]}, full command: {sys.argv}")
if sys.argv[1] in ("editable_wheel"):
    cwd = os.getcwd()
    os.chdir("src/fortran_src/")
    _ = subprocess.call("make clean", shell=True)
    result = subprocess.call("make python", shell=True)
    if result != 0:
        exit(1)
    os.chdir(cwd)

    wrap_file = glob("src/fortran_src/uclchemwrap*.so")[0]
    move(wrap_file, "src/uclchem/uclchemwrap.so")
elif sys.argv[1] in ("bdist_wheel"):
    raise RuntimeError(
        "Cannot install UCLCHEM in the site packages, use `pip install -e .`"
    )
else:
    pass
setup()
