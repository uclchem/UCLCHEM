# This is a backup of the version.py file in the uclchem package.
# If everything goes right, meson will actually use version.py.in to automatically
# generate the file in the build directory from the git tags and commits.
# If this fails for some reason, reenable the meson build line including the version.py
# explicitly and run the meson build command again. It will then use this very file.
__version__ = "0.0.0dev"
