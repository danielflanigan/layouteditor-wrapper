import os
import sys
sys.path.append(r'C:\Users\dflaniga\code\layout-20200504-win-64bit\layout\python')
import LayoutScript as pylayout

from layouteditorwrapper import wrapper


def main():
    layout = wrapper.Layout()
    drawing = layout.drawing()
    return layout, drawing

if __name__ == '__main__':
    layout, drawing = main()
