import setuptools  # this is the "magic" import
from numpy.distutils.core import setup
from glob import glob
from numpy.distutils import exec_command
from shutil import move
from sys import argv

with open("README.md", "r") as fh:
    long_description = fh.read()


#make clean every time slows things but is needed because users will mostly install after changing networks
#and a make never seems to full update the new network
status, output = exec_command.exec_command(
    "make clean", execute_in="src/fortran_src/", use_shell=True
)

status, output = exec_command.exec_command(
    "make python", execute_in="src/fortran_src/", use_shell=True
)

wrap_file = glob("src/fortran_src/uclchemwrap*.so")[0]
move(wrap_file, "src/uclchem/uclchemwrap.so")


exec(open("src/uclchem/__version__.py").read())
setup(
    name="uclchem",  # Replace with your own username
    version=__version__,
    author="Jonathan Holdship",
    author_email="jonholdship@gmail.com",
    description="A package for chemical modelling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://uclchem.github.io",
    packages=setuptools.find_packages(where="src"),
    package_dir={"uclchem": "src/uclchem"},
    package_data={"uclchem": ["uclchemwrap.so"]},
    data_files=[("uclchem", glob("src/*.csv"))],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=["pandas", "numpy", "pyyaml", "matplotlib", "seaborn"],
)
