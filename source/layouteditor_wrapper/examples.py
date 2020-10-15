"""
This module contains a few examples that show how to draw simple shapes using the classes in the
:mod:`wrapper.py <layouteditor_wrapper.wrapper>` module.
"""

from __future__ import absolute_import, division, print_function

import numpy as np


def interdigitated_capacitor(drawing, space, length, width, base, offset, turns, layer, cell_name=None):
    """Create and return a new :class:`~layouteditor_wrapper.wrapper.Cell` containing an interdigitated capacitor (IDC).

    This function is a copy of ``combdrive.layout`` from the LayoutEditor shapes library.

    :param drawing: the object to which the new cell should be added.
    :type drawing: Drawing
    :param space: the spacing between tines.
    :type space: float or int
    :param length: the length of each tine from base to end.
    :type length: float or int
    :param width: the width of each tine.
    :type width: float or int
    :param base: the size of the base connecting the tines on each end.
    :type base: float or int
    :param offset: the distance from the base to the tines of the opposite group.
    :type offset: float or int
    :param int turns: the number of pairs of tines; there is an extra tine on the bottom.
    :param int layer: the layer on which to create the IDC.
    :param str cell_name: the name of the cell; the default includes all the parameters.
    :return: the new cell containing the IDC.
    :rtype: Cell
    """
    if cell_name is None:
        cell_name = 'IDC_{:.3f}_{:.3f}_{:.3f}_{:.3f}_{:.3f}_{:.0f}_{:.0f}'.format(space, length, width, base, offset,
                                                                                  turns, layer)
    cell = drawing.add_cell(cell_name)
    for turn in range(turns):
        left_lower = (width + space) * 2 * turn
        left_upper = (width + space) * (2 * turn + 1)
        cell.add_box(left_lower, base, width, length, layer)
        cell.add_box(left_upper, base + offset, width, length, layer)
    total_width = 2 * turns * (width + space) + width
    cell.add_box(total_width - width, base, width, length, layer)  # rightmost tine
    cell.add_box(0, 0, total_width, base, layer)  # lower base
    cell.add_box(0, base + length + offset, total_width, base, layer)  # upper base
    return cell


def meander(drawing, length, spacing, width, turns, layer, cell_name=None):
    """Create and return a new :class:`~layouteditor_wrapper.wrapper.Cell` containing a meandered inductor.

    The lower left corner of the meander is at ``(0, 0)``, and the center of the first trace is at
    ``(width / 2, width / 2)``. The upper left corner is at ``(0, length)``. The lower right corner is at
    ``(2 * turns * width + (2 * turns - 1) * spacing, 0)`` because the final turn has no connecting piece to the right.

    :param drawing: the drawing to which the new cell should be added.
    :type drawing: Drawing
    :param length: the length of each turn, from outer edge to outer edge.
    :type length: float or int
    :param spacing: the edge-to-edge spacing between traces.
    :type spacing: float or int
    :param width: the width of the trace.
    :type width: float or int
    :param int turns: the number of out-and-back turns.
    :param int layer: the layer on which to create the meander.
    :param str cell_name: the name of the cell; the default includes all the parameters.
    :return: a new cell containing the meander.
    :rtype: Cell
    """
    if cell_name is None:
        cell_name = 'meander_{:.3f}_{:.3f}_{:.3f}_{:.0f}_{:.0f}'.format(length, spacing, width, turns, layer)

    cell = drawing.add_cell(cell_name)
    points = [np.array([width / 2, width / 2])]
    for turn in range(turns):
        points.append(points[-1] + np.array([0, length - width]))
        points.append(points[-1] + np.array([spacing + width, 0]))
        points.append(points[-1] + np.array([0, -(length - width)]))
        points.append(points[-1] + np.array([spacing + width, 0]))
    points.pop()
    cell.add_path(points, int(layer), width=width, cap=2)
    return cell


def double_meander(drawing, length, spacing, width, turns, layer, cell_name=None):
    """Create and return a new :class:`~layouteditor_wrapper.wrapper.Cell` containing a double-wound meandered inductor.

    The lower left corner of the meander is at ``(0, 0)``, and the center of the first trace is at
    ``(width /2, width / 2)``. The upper left corner is at ``(0, length)``. The lower right corner is at
    ``(2 * turns * width + (2 * turns - 1) * spacing, 0)`` because the final turn has no connecting piece to the right.

    :param drawing: the drawing to which the new cell should be added.
    :type drawing: Drawing
    :param length: the length of each turn, from outer edge to outer edge.
    :type length: float or int
    :param spacing: the edge-to-edge spacing between traces.
    :type spacing: float or int
    :param width: the width of the trace.
    :type width: float or int
    :param int turns: the number of turns, where each turn is a pair of traces.
    :param int layer: the layer on which to create the meander.
    :param str cell_name: the name of the cell; the default includes all the parameters.
    :return: a new cell containing the meander.
    :rtype: Cell
    """
    if cell_name is None:
        cell_name = 'double_meander_{:.3f}_{:.3f}_{:.3f}_{:.0f}_{:.0f}'.format(length, spacing, width, turns, layer)

    cell = drawing.add_cell(cell_name)
    out = [np.array([width / 2, width / 2 + width + spacing])]
    back = [np.array([width / 2 + width + spacing, width / 2])]
    horizontal = width + spacing
    vertical = (length - width) - horizontal
    for turn in range(turns):
        if turn % 2:
            out.append(out[-1] + np.array([0, -vertical]))
            back.append(back[-1] + np.array([0, -vertical]))
            out.append(out[-1] + np.array([horizontal, 0]))
            back.append(back[-1] + np.array([3 * horizontal, 0]))
        else:
            out.append(out[-1] + np.array([0, vertical]))
            back.append(back[-1] + np.array([0, vertical]))
            out.append(out[-1] + np.array([3 * horizontal, 0]))
            back.append(back[-1] + np.array([horizontal, 0]))
    # Fix the first point.
    out[0] = np.array([width / 2, width / 2])
    # Delete the last horizontal connection and close the loop.
    out.pop()
    back.pop()
    if turns % 2:
        back[-1] = np.array([back[-1][0], out[-1][1]])
    else:
        out[-1] = np.array([out[-1][0], back[-1][1]])
    out.extend(reversed(back))
    cell.add_path(out, int(layer), width=width, cap=2)
    return cell
