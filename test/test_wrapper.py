from __future__ import absolute_import, division, print_function

import numpy as np
import pytest

from layouteditor_wrapper import wrapper


@pytest.fixture
def layout_and_drawing():
    layout = wrapper.Layout()
    drawing = layout.drawing()
    return layout, drawing


# ToDo: finish creating elements and comparing them to... what?
def test_elements(layout_and_drawing):
    layout, drawing = layout_and_drawing
    noname = drawing.cells['noname']
    new_cell = drawing.add_cell('new_cell')
    cellref = noname.add_cell(cell=new_cell, origin=(1, 1), angle=90)



