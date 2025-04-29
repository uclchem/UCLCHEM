project('uclchemwrap',
        ['c', 'fortran'],
        version : '0.1',
        meson_version: '>= 1.1.0',
        default_options : [
                            'warning_level=1',
                            'buildtype=release'
                          ])

fc = meson.get_compiler('fortran')

fsource_dir = meson.current_source_dir() / 'src' / 'fortran_src'

py = import('python').find_installation('''python''', pure: false)
py_dep = py.dependency()

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

uclchem_wrapper = custom_target('uclchemwrapmodule.c',
  input : uclchemwrap_source,  # .f so no F90 wrappers
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
                     # dvode_lib,
                     ],
                     link_with : dvode_lib,
                     install : true)
