import setuptools  # this is the "magic" import
from numpy.distutils.core import setup, Extension
from glob import glob
from numpy.distutils import exec_command
with open("README.md", "r") as fh:
    long_description = fh.read()

try:
    #exec_command.exec_command( "make clean", execute_in='src/fortran_src/', use_shell=True)
    exec_command.exec_command( "make python", execute_in='src/fortran_src/', use_shell=True)
except Exception as e:
    print("Making UCLCHEM failed")
    print("You likely do not have gfortran and/or cmake installed")
    print("Printing python exception...\n")
    print(e)

# wrap_ext=Extension('uclchemwrap', ['src/fortran_src/uclchemwrap.so'])

print(setuptools.find_packages(where='src'))

exec(open('src/uclchem/version.py').read())
setup(
    name="uclchem", # Replace with your own username
    version=__version__,
    author="Jonathan Holdship",
    author_email="jonholdship@gmail.com",
    description="A package for chemical modelling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://uclchem.github.io",
    packages=setuptools.find_packages(where='src'),
    package_dir={'uclchem': 'src/uclchem'},
    package_data={'uclchem': ['uclchemwrap.so']},
    data_files=[('uclchem',glob('src/*.csv'))],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['pandas','numpy']
)