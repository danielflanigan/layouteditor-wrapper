layouteditor-wrapper
====================

[![docs](https://readthedocs.org/projects/layouteditor-wrapper/badge/?version=latest)](https://layouteditor-wrapper.readthedocs.io/en/latest/)

A Python wrapper for [LayoutScript](https://layouteditor.org/layoutscript), the scripting interface of [Juspertor LayoutEditor](https://layouteditor.com).
This package is not developed or endorsed in any way by Juspertor.

Introduction
------------

The main goal of this package is to provide an interface to LayoutScript that uses familiar Python objects and is harmonized with the LayoutEditor GUI (especially in the occasional cases where LayoutScript uses different conventions). For example, instead of using the LayoutScript ``point`` and ``pointArray`` classes, the interface uses numpy arrays with size 2 containing (x, y) coordinates. This allows points to be treated as vectors for mathematical operations such as calculating distances or angles between them. The interface by default uses the float user unit, which is often more convenient than using the integer database units.

Contents
--------

The name of the installed package is ``layouteditor_wrapper``. It includes the following modules::

- ``wrapper.py``, the core of the package, which contains stateless wrapper classes for the LayoutScript objects;
- ``transmission_line.py``, which contains classes and functions for drawing transmission lines; 
- ``cpw.py``, which draws co-planar waveguide structures, especially those useful for drawing superconducting resonators; 
- ``examples.py``, which contains a few example functions that create simple components;
- ``__main__.py``, a script for interactive use.

Install
-------

This package is compatible with versions of LayoutEditor that include LayoutScript instead of the old pylayout interface. (Earlier versions of this package were developed for pylayout.) The only dependencies are ``setuptools`` for installation and ``numpy``.

**Linux**

Re-organized package not yet installed on Linux.

**macOS**

If LayoutEditor has been installed in ``/Applications``, then the LayoutScript files are in ``/Applications/layout.app/Contents/python``.

**Windows**

Download the .zip Windows installer (not the .msi) and unzip it to the desired location. The .zip installer includes a Python 3.7 distribution, and LayoutScript seems to rely on the specific DLLs that come with this distribution. The recommended use of this package is to install it in the distribution included with LayoutEditor. (I have not succeeded in creating a ``conda`` environment in which ``LayoutScript`` will run.)

Change directory to the LayoutEditor Python 3.7 distribution root directory, then run the simple installer script ``install.py`` in the ``layouteditor-wrapper`` root directory::

    C:\path\to\layout-yyyymmdd-win-64bit\layout\python\python37> python C:\path\to\layouteditor-wrapper\install.py

On Windows, this script does the following:
- Bootstrap the packaging tools ``pip`` and ``setuptools``, which are not included with the LayoutEditor distribution;
- Use ``pip`` to install the necessary packages listed in ``setup.py``;
- Use ``pip`` to perform an editable install of the ``layouteditor_wrapper`` package;
- Write a batch file ``windows-environment.bat`` in the package root directory.

The recommended way to run is to create a shortcut with the target suggested by the script, which opens a command window and runs the batch file to setup the environment. The batch file appends the directories to ``PATH`` that are necessary for the Python 3.7 executable (and ``pip``) to run; it also creates an environment variable ``LAYOUTSCRIPT_PATH`` that ``layouteditor_wrapper`` uses to import the LayoutScript components. With the packaging tools installed, you can customize the environment with any packages needed to draw your layouts.

Credits
-------

Earlier generations of this package were developed in the [Experimental Cosmology Group at Columbia University](https://github.com/ColumbiaCMB) along with Glenn Jones and Heather McCarrick.
