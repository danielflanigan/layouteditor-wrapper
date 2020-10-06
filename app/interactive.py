from layouteditor_wrapper import wrapper


def main():
    layout = wrapper.Layout()
    drawing = layout.drawing()
    return layout, drawing


if __name__ == '__main__':
    l, d = main()
    print("Variable 'l' is {}".format(type(l)))
    print("Variable 'd' is {}".format(type(d)))
