
project('uclchemwrap',
        ['c', 'fortran'],
        version : run_command('sh', '-c',
  'git describe --tags --always 2>/dev/null | ' +
  'python3 -c "import sys, re; ' +
  'v = sys.stdin.read().strip(); ' +
  'm = re.match(r\'v?(\\d+\\.\\d+\\.\\d+)(?:-(\\d+)-g([0-9a-f]+))?\', v); ' +
  'print(f\'{m.group(1)}.post{m.group(2)}+{m.group(3)}\' if m and m.group(2) else (m.group(1) if m else \'0.0.0dev+nogit\'))"'
).stdout().strip(),
        meson_version: '>= 1.7.0',
        default_options : [
                            'warning_level=1',
                            'buildtype=release'
                          ])

# Check if this is a makerates-only build
makerates_only = get_option('makerates_only')

fc = meson.get_compiler('fortran')

fsource_dir = meson.current_source_dir() / 'src' / 'fortran_src'

py = import('python').find_installation('''python''', pure: false)
py_dep = py.dependency()


version_file = vcs_tag(
  input : 'src/uclchem/version.py.in',
  output : 'version.py',
  install_dir: py.get_install_dir() / 'uclchem',
  install : true,
  # Literal copy paste since we cannot define a variable before the project() call
  command : run_command('sh', '-c',
  'git describe --tags --always 2>/dev/null | ' +
  'python3 -c "import sys, re; ' +
  'v = sys.stdin.read().strip(); ' +
  'm = re.match(r\'v?(\\d+\\.\\d+\\.\\d+)(?:-(\\d+)-g([0-9a-f]+))?\', v); ' +
  'print(f\'{m.group(1)}.post{m.group(2)}+{m.group(3)}\' if m and m.group(2) else (m.group(1) if m else \'0.0.0dev+nogit\'))"'
).stdout().strip(),
)

# Only compile Fortran extension if not makerates_only build
if not makerates_only

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

inc_np = include_directories(incdir_numpy, incdir_f2py)
# gh-25000
quadmath_dep = fc.find_library('quadmath', required: false)

dvode_sources = files(fsource_dir / 'dvode.f90')

dvode_lib = static_library('dvode', dvode_sources,
  # include_directories: inc_np,
  dependencies: [quadmath_dep],
)

uclchemwrap_source = files(fsource_dir / 'constants.f90',
                     fsource_dir / 'defaultparameters.f90',
                     fsource_dir / 'f2py_constants.f90',
                     fsource_dir / 'physics-core.f90',
                     fsource_dir / 'network.f90',
                     fsource_dir / 'cloud.f90',
                     fsource_dir / 'hotcore.f90',
                     fsource_dir / 'postprocess.f90',
                     fsource_dir / 'surfacereactions.f90',
                     fsource_dir / 'sputtering.f90',
                     fsource_dir / 'cshock.f90',
                     fsource_dir / 'jshock.f90',
                     fsource_dir / 'collapse.f90',
                     fsource_dir / 'photoreactions.f90',
                     fsource_dir / 'rates.f90',
                     fsource_dir / 'odes.f90',
                     fsource_dir / 'chemistry.f90',
                     fsource_dir / 'io.f90',
                     fsource_dir / 'wrap.f90'
                    )

# Custom target to generate the f2py wrappers
uclchem_wrapper = custom_target('uclchemwrapmodule.c',
  input : uclchemwrap_source,
  output : ['uclchemwrapmodule.c', 'uclchemwrap-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'uclchemwrap', '--lower']
)

py.extension_module('uclchemwrap',
                     [
                      uclchemwrap_source,
                      uclchem_wrapper,
                      fortranobject_c
                     ],
                     include_directories: [
                     inc_np,
                     ],
                     dependencies : [
                     py_dep,
                     quadmath_dep,
                     ],
                     link_with : dvode_lib,
                     fortran_args: ['-O3', '-fPIC', '-ffree-line-length-0'],
                     install : true)

endif  # End of Fortran compilation section

usource = 'src' / 'uclchem'

# Always install makerates/ subpackage (core functionality)
py.install_sources([
  usource / 'makerates' / '__init__.py',
  usource / 'makerates' / 'io_functions.py',
  usource / 'makerates' / 'makerates.py',
  usource / 'makerates' / 'network.py',
  usource / 'makerates' / 'reaction.py',
  usource / 'makerates' / 'species.py'
], subdir: 'uclchem' / 'makerates')

# Always install nested f2py_constants.f90 template
py.install_sources([
  usource / 'makerates' / 'fortran_templates' / 'f2py_constants.f90'
], subdir: 'uclchem' / 'makerates' / 'fortran_templates')

# Full installation with all Python modules
py.install_sources([
  usource / '__init__.py',
  # usource / 'version.py', # version is automatically generated by vcs_tag and placed in the site-packages directory
  usource / 'analysis.py',
  usource / 'constants.py',
  usource / 'debug.py',
  usource / 'model.py',
  usource / 'trialmodel.py',
#  usource / 'legacy_model.py',
  usource / 'plot.py',
  usource / 'reactions.csv',
  usource / 'species.csv',
  usource / 'tests.py',
  usource / 'utils.py',
], subdir:'uclchem')

# Install makerates/ subpackage
py.install_sources([
  usource / 'makerates' / '__init__.py',
  usource / 'makerates' / 'io_functions.py',
  usource / 'makerates' / 'makerates.py',
  usource / 'makerates' / 'network.py',
  usource / 'makerates' / 'reaction.py',
  usource / 'makerates' / 'species.py'
], subdir: 'uclchem' / 'makerates')

# Install nested f2py_constants.f90 template
py.install_sources([
  usource / 'makerates' / 'fortran_templates' / 'f2py_constants.f90'
], subdir: 'uclchem' / 'makerates' / 'fortran_templates')
