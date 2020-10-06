import os
import platform
import subprocess
import sys


def windows_install():
    layouteditor_wrapper_path = os.path.dirname(os.path.abspath(__file__))
    print(subprocess.check_output([sys.executable, '-m', 'ensurepip', '--default-pip']).decode())
    print(subprocess.check_output([sys.executable, '-m', 'pip', 'install', '-e', layouteditor_wrapper_path]).decode())
    # Strip the last two directories from C:\path\to\layout-20200801-win-64bit\layout\python\python37
    layout_path = os.path.split(os.path.split(sys.exec_prefix)[0])[0]
    lines = [
        r"@echo off",
        r"set layout_path={}".format(layout_path),
        r"set LAYOUTSCRIPT_PATH=%layout_path%\python",
        r"set PATH=%PATH%;%layout_path%\bin;%layout_path%\python\python37;%layout_path%\python\python37\Scripts",
        r"set layout_path=",
        r"echo Set LAYOUTSCRIPT_PATH and appended to PATH."
    ]
    environment_filename = os.path.join(layouteditor_wrapper_path, 'windows-environment.bat')
    with open(environment_filename, 'w') as f:
        f.writelines(["{}{}".format(line, os.linesep) for line in lines])
        print("Wrote {}".format(environment_filename))
        print("To start a command window with the layout environment, create a link with the following target:")
        print(r'%windir%\System32\cmd.exe "/K" {}'.format(environment_filename))


if __name__ == '__main__':
    operating_system = platform.system()
    if operating_system == 'Windows':
        windows_install()
    elif operating_system == 'Linux':
        print("Linux installation not yet suppported.")
        sys.exit(1)
    elif operating_system == 'Darwin':
        print("macOS installation not yet suppported.")
        sys.exit(1)
    else:
        print("platform.system() returned {}".format(operating_system))
        sys.exit(1)
