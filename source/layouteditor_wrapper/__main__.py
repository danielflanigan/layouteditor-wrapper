"""This module is intended to be used interactively; see ``__init__.py``."""
import sys

from layouteditor_wrapper import wrapper


def main(filename=None):
    layout = wrapper.Layout()
    if filename is not None:
        layout.load(filename)
    drawing = layout.drawing()
    return layout, drawing


if __name__ == '__main__':
    try:
        l, d = main(filename=sys.argv[1])
    except IndexError:
        l, d = main()
    print("Variable 'l' is {}".format(type(l)))
    print("Variable 'd' is {}".format(type(d)))
