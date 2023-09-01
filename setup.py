import setuptools  # this is the "magic" import
import subprocess
import sys

from setuptools import os, setup

from glob import glob
from shutil import move

#make clean every time slows things but is needed because users will mostly install after changing networks
#and a make never seems to full update the new network
cwd = os.getcwd()
os.chdir('src/fortran_src/')
result = subprocess.call('make clean', shell=True)
# if result != 0:
#     exit(1)
result = subprocess.call('make python', shell=True)
if result != 0:
    exit(1)
os.chdir(cwd)

wrap_file = glob("src/fortran_src/uclchemwrap*.so")[0]
move(wrap_file, "src/uclchem/uclchemwrap.so")

setup()
