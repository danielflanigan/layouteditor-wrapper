import os
import sys
sys.path.append('C:\\Program Files (x86)\\LayoutEditor\\python')
import LayoutScript

from layouteditorwrapper import wrapper


def main():
    layout = wrapper.Layout()
    drawing = layout.drawing()
    return layout, drawing

if __name__ == '__main__':
    layout, drawing = main()
