# layouteditor-wrapper

A wrapper for LayoutScript, the Python module for Juspertor LayoutEditor. This package is not developed or endorsed in any way by Juspertor.

## Contents

The name of the installed package is `layouteditor_wrapper`.
It includes the following modules:
- `wrapper.py`, which contains wrapper classes for the LayoutScript objects;
- `path.py`, which contains classes and functions useful for drawing co-planar waveguide components; 
- `components.py`, which contains a few example functions that create useful components.

There is also a template script `app/interactive.py` that starts LayoutEditor with `Layout` and `Drawing` objects in the namespace:
```
$ /path/to/layouteditor-wrapper$ python -i app/interactive.py
Variable 'l' is <class 'layouteditor_wrapper.wrapper.Layout'>
Variable 'd' is <class 'layouteditor_wrapper.wrapper.Drawing'>
>>> d.cells
OrderedDict([('noname', <layouteditor_wrapper.wrapper.Cell object at 0x00000262A2FAC518>)])
```

This script can be used as a template to create layouts entirely in code.

## Installation

The only dependencies are `setuptools` for installation and `numpy`.

### Linux

Re-organized package not yet installed on Linux.

### macOS

If LayoutEditor has been installed in `/Applications`, then `LayoutScript` is in `/Applications/layout.app/Contents/python`.

### Windows

Download the .zip Windows installer (not the .msi) and unzip it to the desired location.
The .zip installer includes a Python 3.7 distribution, and LayoutScript seems to rely on the specific DLLs that come with this distribution. 
I did not succeed in running LayoutScript in a conda environment.
Change directory to the LayoutEditor Python 3.7 distribution root directory, then run the simple installer script `install.py` in the `layouteditor-wrapper` root directory:
```
C:\path\to\layout-yyyymmdd-win-64bit\layout\python\python37> python C:\path\to\layouteditor-wrapper\install.py
``` 
On Windows, this script does the following:
- Bootstrap the packaging tools `pip` and `setuptools`, which are not included with the LayoutEditor distribution.
- Use pip to install the packages in `requirements.txt`.
- Use pip to perform an editable install of the `layouteditor_wrapper` package.
- Write a batch file `windows-environment.bat` in the package root directory.
The cleanest way to run is to create the suggested link, which opens a command window and runs the batch file. 
The batch file appends the directories to `PATH` that are necessary for the Python 3.7 executable (and pip) to run.
It also creates an environment variable `LAYOUTSCRIPT_PATH` that `layouteditor_wrapper` uses to import the LayoutScript components.
This method avoids polluting the `PATH`.
With the packaging tools installed, you can customize the environment with any packages needed to draw your layouts.

## Credits

Earlier generations of this package were developed in the [Experimental Cosmology Group at Columbia University](https://github.com/ColumbiaCMB) along with Glenn Jones and Heather McCarrick.