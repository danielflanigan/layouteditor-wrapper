"""Usage: > python -i -m layouteditor_wrapper <optional filename>

This module is intended to be used interactively. Assuming that the layouteditor_wrapper package is installed:
> python -i -m layouteditor_wrapper
Variable 'l' is <class 'layouteditor_wrapper.wrapper.Layout'>
Variable 'd' is <class 'layouteditor_wrapper.wrapper.Drawing'>
>>> d.cells
OrderedDict([('noname', <layouteditor_wrapper.wrapper.Cell object at 0x000001CDFF756B38>)])
"""
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
