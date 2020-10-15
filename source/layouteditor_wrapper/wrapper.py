"""This module is a wrapper for LayoutScript, the Python module of LayoutEditor.

This module includes several types of objects that each wrap a LayoutScript object. The naming conventions for
LayoutScript objects are the following: in methods and functions, LayoutScript objects have a ``ls_`` prefix; as
attributes of wrapper classes, the attribute ``.ls`` is the wrapped object. For example, the constructor for Cell
includes ``ls_cell``, which is a ``LayoutScript.cell`` object; also, if cell is a
:class:`~layouteditor_wrapper.wrapper.Cell` object, then ``cell.ls`` is the wrapped ``LayoutScript.cell`` object.

Where possible, the classes use the same conventions and notation as the LayoutEditor GUI, which sometimes differ
from from those used by LayoutScript. For example, this interface by default expresses lengths in terms of the user
unit, like the GUI, but can also use the integer database units, like LayoutScript.

The classes representing the objects in a drawing are nearly-stateless wrappers for the underlying LayoutScript objects.
The wrapper class methods and properties do not have any attributes of their own, but just provide an interface to the
LayoutScript objects. Since the wrappers store no state, there is no need for them to be unique.
For example, every call to :meth:`Drawing.cells <layouteditor_wrapper.wrapper.Drawing.cells>` generates a new list of :class:`~layouteditor_wrapper.wrapper.Cell` objects.  These can be cached for speed, if necessary.

In method docstrings the word *point* refers to a point with two coordinates. Internally, the classes use numpy arrays
with shape ``(2,)``, but their methods accept anything that allows ``point[0]`` and ``point[1]`` to be indexed, such as
a tuple.
"""
from __future__ import absolute_import, division, print_function
from collections import OrderedDict
import os
import sys
import warnings

import numpy as np

try:
    import LayoutScript as ls
except ImportError:
    try:
        layoutscript_path = os.environ['LAYOUTSCRIPT_PATH']
        sys.path.insert(0, layoutscript_path)
        import LayoutScript as ls
        sys.path.remove(layoutscript_path)
    except KeyError:
        warnings.warn("LayoutScript is not available on the Python path, and LAYOUTSCRIPT_PATH is not set.")


# The two following simple functions are available to code that uses (lists of) numpy arrays as points.
# This makes it easy for methods to accept lists of tuples, for example.

def to_point(indexable):
    """Return a numpy ndarray in the two-dimensional point format used by this module.

    :param indexable: an indexable object with integer indices 0 and 1, such as a two-element tuple.
    :type indexable: object
    :return: an array with shape (2,) containing the values at these two indices.
    :rtype: numpy.ndarray
    """
    return np.array([indexable[0], indexable[1]])


def to_point_list(iterable):
    """Return a list of numpy arrays in the two-dimensional point format used by this module.

    :param iterable iterable: an iterable of indexable objects that all have integer indices 0 and 1.
    :return: a list of arrays with shape (2,) containing the values at these two indices.
    :rtype: list[numpy.ndarray]
    """
    return [to_point(point) for point in iterable]


def instantiate_element(ls_element, drawing):
    """Instantiate the appropriate wrapper class for the given LayoutScript element type.

    :param ls_element: a LayoutScript element object.
    :type ls_element: LayoutScript.element
    :param drawing: a Drawing object.
    :type drawing: :class:
    :return: a wrapper instance for the LayoutScript element.
    :rtype: :class:`~layouteditor_wrapper.wrapper.Element`
    """
    elements = (Box, Cellref, CellrefArray, Circle, Path, Polygon, Text)
    for element in elements:
        if getattr(ls_element, 'is' + element.__name__)():
            return element(ls_element, drawing)
    raise ValueError("Unknown LayoutScript element.")


class Layout(object):
    """Wrap a LayoutScript.layout object."""

    def __init__(self):
        """Instantiate a new Layout, which creates a new LayoutScript layout."""
        self.ls = ls.project.newLayout()

    def drawing(self, use_user_unit=True, auto_number=False):
        """Return a wrapper.Drawing object that wraps the current LayoutScript.drawing object."""
        return Drawing(self.ls.drawing, use_user_unit=use_user_unit, auto_number=auto_number)

    def load(self, filename):
        """Load the file from disk with the given filename."""
        self.ls.open(filename)

    def save(self):
        """Save the current layout to the filename from which it was loaded.

        :raises: RuntimeError if the filename attribute of the layout is empty.
        """
        if self.filename:
            self.ls.save()
        else:
            raise RuntimeError("Filename is empty, probably because the layout was not loaded from disk.")

    def save_as(self, filename):
        """Save the current layout to the given filename."""
        self.ls.saveAs(filename)

    # ToDo: can this be set?
    @property
    def filename(self):
        """The filename used to load the current layout from disk, or an empty string if it was never loaded."""
        return self.ls.filename


class Drawing(object):
    """Wrap a LayoutScript.drawingField object."""

    def __init__(self, ls_drawing, use_user_unit=True, auto_number=False, validate_cell_names=True):
        """Create a new Drawing object from the given ``LayoutScript.drawingField`` object.

        All of the LayoutScript classes expect and return integer database units. However, this class interface can
        instead use the float user units, which are often more convenient. For example, if the database units are nm
        and the user units are µm, then if use_use_unit is True then all class methods will expect and return length
        values in µm. If it is False, then all class methods will expect and return length values in nm. The
        underlying data are always saved as integers in database units, so these must be small enough for the desired
        precision.

        :param ls_drawing: the drawing instance to wrap.
        :type ls_drawing: LayoutScript.drawingField
        :param bool use_user_unit: if True, all values input to and returned from this instance are in user units
                                   (float); if False, they are in database units (int).
        :param bool auto_number: if True, the `add_cell` method ensures unique cell names by appending an integer to the
                                 names of all cells it creates in this drawing.
        :param bool validate_cell_names: if True, check that the new cell name (after auto-numbering, if applicable) is
                                         not already used to avoid corruption; use False to add large numbers of cells
                                         more rapidly.
        :return: the new Drawing object.
        :rtype: Drawing
        """
        self.ls = ls_drawing
        self.use_user_unit = use_user_unit
        self.auto_number = auto_number
        if auto_number:
            self._cell_number = 0
        self.validate_cell_names = validate_cell_names

    @property
    def database_unit(self):
        """The physical length of the database unit in meters (``float``); the default seems to be 1e-9, or 1 nm."""
        return self.ls.databaseunits

    @database_unit.setter
    def database_unit(self, unit):
        self.ls.databaseunits = unit

    @property
    def user_unit(self):
        """The ratio of the database unit to the user unit (``float``).

        All points are saved as integer values in database units, so this number is the data resolution in user units.
        For example, if the database units are nm and the user units are µm, then this value is 0.001 = 1 nm / 1 µm.
        """
        return self.ls.userunits

    @user_unit.setter
    def user_unit(self, unit):
        self.ls.userunits = unit

    def to_database_units(self, value_or_array):
        """Return the given value or array scaled to the database units.

        The behavior of this function depends on the value of the `__init__` argument use_user_unit in the following way: if
        this is True, then this function expects values in user units, which are scaled appropriately and rounded to
        the nearest integer; if False, this function expects values in database units, which are simply rounded to
        the nearest integer. The return value is always an int or a numpy ndarray of dtype either np.int32 or
        possibly np.int64.

        :param value_or_array: a value or array to be converted to integer database units.
        :return: the converted int or int array.
        """
        try:
            # Ensure that zero-length arrays are treated as scalars and converted to ints
            assert value_or_array.shape
            if self.use_user_unit:
                return np.round(value_or_array / self.user_unit).astype(np.int)
            else:
                return np.round(value_or_array).astype(np.int)
        # The LayoutScript methods do not accept zero-length numpy arrays where they expect scalars
        # AttributeError: not an array-like object; AssertionError: zero-length array.
        except (AttributeError, AssertionError):
            if self.use_user_unit:
                return int(round(value_or_array / self.user_unit))
            else:
                return int(value_or_array)

    def from_database_units(self, value_or_array):
        """Return the given value or array scaled from the database units.

        :return: float or numpy.ndarray of float
        """
        try:
            if self.use_user_unit:
                return (value_or_array * self.user_unit).astype(np.float)
            else:
                return value_or_array.astype(np.int)
        except AttributeError:  # not an array-like object
            if self.use_user_unit:
                return float(value_or_array * self.user_unit)
            else:
                return int(value_or_array)

    @property
    def cells(self):
        """An OrderedDict that maps all the cell names (``str``) in the layout to :class:`Cell` instances (read-only).

        Modifications to the returned OrderedDict are **not** propagated back to the layout. To add a cell to the
        drawing use :meth:`~Drawing.add_cell`, to delete a cell from the drawing use :meth:`~Drawing.delete_cell`, and
        to change the name of a cell assign the new name to :attr:`~Cell.name`. Modifications to the layout are **not**
        propagated to the returned OrderedDict, and this property must be accessed again to obtain an updated mapping.

        LayoutScript uses a linked list to store its cell instances. (``drawing.thisCell`` -> ``cellList``,
        ``cellList.thisCell`` -> ``cell``, ``cellList.nextCell`` -> ``cellList`` or None.) New cells are prepended to
        the internal linked list, so the index of a cell changes as new cells are added.

        Note that LayoutEditor forbids changing the name of a cell to an existing name but LayoutScript does not. This
        module attempts to enforce unique cell names, so the mapping returned by this method should always be valid.
        """
        n_cells = 0
        cell_dict = OrderedDict()
        current = self.ls.firstCell
        while current is not None:
            cell = Cell(current.thisCell, self)
            cell_dict[cell.name] = cell
            n_cells += 1
            current = current.nextCell
        if n_cells != len(cell_dict):
            raise RuntimeError("Duplicate cell name; use method cell_list to diagnose and repair.")
        return cell_dict

    @property
    def cell_list(self):
        """A list of all the :class:`Cell` instances in the layout (read-only).

        Modifications to the returned list are **not** propagated back to the layout. To add a cell to the drawing
        use :meth:`~Drawing.add_cell`, and to delete a cell from the drawing use :meth:`~Drawing.delete_cell`.
        Modifications to the layout are **not** propagated to the returned list, and this property must be accessed
        again to obtain an updated list.

        Note that new cells added to the drawing are prepended to the internal linked list, so the index of a cell in
        the returned list will change as cells are added or deleted.
        """
        cl = list()
        current = self.ls.firstCell
        while current is not None:
            cl.append(Cell(current.thisCell, self))
            current = current.nextCell
        return cl

    def _np_to_layoutscript(self, array):
        """Create and return a point object (without adding it to the drawing) with the coordinates scaled to integer
        database units.

        :param array: a two-element array containing the x- and y-coordinates of the point in either user units or
        database units; see documentation for `use_user_unit` in ` __init__`.
        :type array: numpy.ndarray
        :return: LayoutScript.point
        """
        return ls.point(int(self.to_database_units(array[0])), int(self.to_database_units(array[1])))

    def _layoutscript_to_np(self, point):
        """Create and return a numpy.ndarray with shape ``(2,)``, with the coordinates scaled from integer database
        units according to `use_user_unit`.

        :param point: a point containing x- and y-coordinates.
        :type array: LayoutScript.point
        :return: numpy.ndarray
        """
        return self.from_database_units(np.array([point.x(), point.y()]))

    def _to_point_array(self, list_of_numpy_arrays):
        """Return a ``LayoutScript.pointArray`` of ``LayoutScript.point`` objects converted from the given list of
        numpy ndarrays with shape ``(2,)``.

        :return: LayoutScript.pointArray
        """
        pa = ls.pointArray(len(list_of_numpy_arrays))
        for i, array in enumerate(list_of_numpy_arrays):
            pa.setPoint(i, self._np_to_layoutscript(array))
        return pa

    def _to_list_of_numpy_arrays(self, point_array):
        array_list = []
        for i in range(point_array.size()):
            array_list.append(self._layoutscript_to_np(point_array.point(i)))
        return array_list

    def add_cell(self, name):
        """Add an empty cell with the given name to the drawing and return a new :class:`Cell` instance that wraps a
        ``LayoutScript.cell`` instance.

        If ``auto_number`` is True, the string ``_x`` is appended to the given cell name, where x is the number of cells
        that have been created so far by this class. This allows code to repeatedly call this method with the same name
        string without producing an error, since the appended number usually ensures that their names will be unique.

        If ``validate_cell_names`` is True, check before adding the cell that the given name is not already used in the
        drawing. This avoids corruption of the cell list at the cost of recreating and traversing the ``cells``
        OrderedDict. Set this to False to skip this check.

        :param str name: the name of the new cell.
        :return: the new, empty cell.
        :rtype: Cell
        :raises ValueError: if a cell with the given name already exists in the drawing.
        """
        if self.auto_number:
            name = '{}_{}'.format(name, self._cell_number)
        if self.validate_cell_names and name in self.cells:
            raise ValueError("A Cell with name '{}' already exists.".format(name))
        ls_cell = self.ls.addCell().thisCell
        ls_cell.cellName = name
        self._cell_number += 1
        # Adding a cell does not update currentCell. Without the line below, boolean operations (and probably others)
        #  will operate on currentCell instead of the cell created by this method.
        # ToDo: this may no longer be necessary
        self.ls.currentCell = ls_cell
        return Cell(ls_cell, self)

    def delete_cell(self, cell):
        """Delete the given cell from the drawing.

        The value of the current cell is unspecified after the deletion.

        :param cell: the cell to delete.
        :type cell: Cell
        :return: None
        :raises KeyError: if the given Cell is not in the drawing.
        """
        self.delete_cell_name(cell.name)

    def delete_cell_name(self, name):
        """Delete a cell with the given name.

        The value of the current cell is unspecified after the deletion.

        :param str name: the name of the cell to delete.
        :return: None
        :raises KeyError: if no cell in the drawing has the given name.
        """
        self.current_cell_name = name
        self.ls.deleteCurrentCell()

    @property
    def current_cell(self):
        """The Cell instance that wraps the ``currentCell`` attribute of the drawing."""
        return Cell(ls_cell=self.ls.currentCell, drawing=self)

    @current_cell.setter
    def current_cell(self, cell):
        self.ls.currentCell = cell.ls

    @property
    def current_cell_name(self):
        """The name (``str``) of the current cell."""
        return self.current_cell.name

    @current_cell_name.setter
    def current_cell_name(self, name):
        self.current_cell = self.cells[name]


class Cell(object):
    """Wrap a ``LayoutScript.cell`` object."""

    def __init__(self, ls_cell, drawing):
        """
        :param ls_cell: the cell object to wrap.
        :type ls_cell: LayoutScript.cell
        :param drawing: the drawing object that contains this cell.
        :type drawing: LayoutScript.drawingField
        """
        self.ls = ls_cell
        self.drawing = drawing

    @property
    def name(self):
        """The name of this cell (``str``), which should be unique.

        :raises ValueError: if the name is changed to a name already used in the drawing.
        """
        try:
            return self.ls.cellName.toAscii().data()
        except AttributeError:
            return self.ls.cellName

    @name.setter
    def name(self, new_name):
        if new_name != self.name and new_name in self.drawing.cells:
            raise ValueError("A Cell with name '{}' already exists.".format(new_name))
        self.ls.cellName = str(new_name)

    @property
    def elements(self):
        """Generate and return a list of all elements in the cell.

        Since LayoutScript elements do not have an internal name, unlike cells, there is no obvious way to create
        dictionary keys for them. Note that new elements are prepended to the internal linked list as they are added
        to a cell, so the index of each element is not constant.
        """
        element_list = []
        current = self.ls.firstElement
        while current is not None:
            element_list.append(instantiate_element(current.thisElement, self.drawing))
            current = current.nextElement
        return element_list

    def __str__(self):
        return 'Cell {}: {}'.format(self.name, [str(e) for e in self.elements])

    def subtract(self, positive_layer, negative_layer, result_layer, delete=True):
        """Perform the boolean operation ``positive - negative = result`` on the given layers.

        :param int positive_layer: Structures on this layer remain unless subtracted.
        :param int negative_layer: Structures on this layer are subtracted from the positive layer.
        :param int result_layer: The structures resulting from the subtraction are created on this layer.
        :param bool delete: if True, delete all structures on the positive and negative layers after the subtraction.
        :return: None
        """
        self.drawing.ls.setCell(self.ls)
        bh = ls.booleanHandler(self.drawing.ls)
        # bh.boolOnLayer(positive_layer, negative_layer, result_layer, ls.string('A-B'), 0, 0, 0)
        bh.boolOnLayer(positive_layer, negative_layer, result_layer, 'A-B', 0, 0, 0)
        if delete:
            self.drawing.ls.deleteLayer(positive_layer)
            self.drawing.ls.deleteLayer(negative_layer)

    def add_cell(self, cell, origin, angle=0):
        """Add a single cell to this cell.

        The origin coordinates are expected to be floats in user units if use_user_unit is True, and ints in database
        units if it is False.

        :param cell: the cell object to add to this cell.
        :type cell: Cell
        :param origin: a point containing the origin x- and y-coordinates.
        :type origin: numpy.ndarray
        :param float angle: the cell orientation in degrees.
        :return: a Cellref object with a reference to the given Cell.
        :rtype: Cellref
        """
        ls_cell = self.ls.addCellref(cell.ls, self.drawing._np_to_layoutscript(to_point(origin)))
        cell = Cellref(ls_cell, self.drawing)
        cell.angle = angle
        return cell

    def add_cell_array(self, cell, origin=(0, 0), step_x=(0, 0), step_y=(0, 0), repeat_x=1, repeat_y=1, angle=0):
        """Return a CellrefArray object produced by adding to this cell an array of the given cells.

        The coordinates of the origin and step poins are expected to be floats in user units if use_user_unit is True,
        and ints in database units if it is False.

        :param cell: the Cell object to use to create the CellrefArray in this cell.
        :type cell: wrapper.Cell
        :param origin: a point containing the origin x- and y-coordinates.
        :type origin: numpy.ndarray
        :param step_x: a point containing the x- and y-increment for all cells in each row.
        :type step_x: numpy.ndarray
        :param step_y: a point containing the x- and y-increment for all cells in each column.
        :type step_y: numpy.ndarray
        :param int repeat_x: the number of columns.
        :params int repeat_y: the number of rows.
        :param float angle: the cell orientation in degrees.
        :return: a CellrefArray object with a reference to the given Cell.
        :rtype: CellrefArray
        """
        repeat_x = int(repeat_x)
        repeat_y = int(repeat_y)
        # Strange but true: the constructor for this object expects three points that are different from both the
        # points returned by getPoints() and the GUI interface points.
        ls_origin = to_point(origin)
        ls_total_x = repeat_x * to_point(step_x) + ls_origin
        ls_total_y = repeat_y * to_point(step_y) + ls_origin
        point_array = self.drawing._to_point_array([ls_origin, ls_total_x, ls_total_y])
        ls_cell_array = self.ls.addCellrefArray(cell.ls, point_array, repeat_x, repeat_y)
        cell_array = CellrefArray(ls_cell_array, self.drawing)
        cell_array.angle = angle
        return cell_array

    def add_box(self, x, y, width, height, layer):
        """Add a rectangular box to this cell and return the corresponding object.

        All length values are expected to be floats in user units if use_user_unit is True, and ints in database units
        if it is False.

        :param x: the x-coordinate of the origin.
        :type x: float or int
        :param y: the y-coordinate of the origin.
        :type y: float or int
        :param width: the horizontal width of the box, positive if the box is to extend to the right from the origin.
        :type width: float or int
        :param height: the vertical height of the box, positive if the box is to extend upward from the origin.
        :type height: float or int
        :param int layer: the layer on which the box is created.
        :return: the new box.
        :rtype: Box
        """
        ls_box = self.ls.addBox(self.drawing.to_database_units(x),
                                self.drawing.to_database_units(y),
                                self.drawing.to_database_units(width),
                                self.drawing.to_database_units(height),
                                int(layer))
        return Box(ls_box, self.drawing)

    # ToDo: addChamferedBox
    # ToDo: addCircleBox
    # ToDo: addEllipse
    # ToDo: addRoundedBox

    def add_circle(self, origin, radius, layer, number_of_points=0):
        """Add a circular polygon to this cell and return the corresponding object.

        All length values are expected to be floats in user units if use_user_unit is True, and ints in database
        units if it is False. Note that LayoutScript considers any regular polygon with 8 or more points to be a
        circle, and once created a circle has no special properties.

        :param origin: a point containing the (x, y) coordinates of the circle center.
        :type origin: numpy.ndarray
        :param radius: the circle radius.
        :type radius: float or int
        :param int layer: the layer on which the circle is created.
        :param int number_of_points: the number of unique points to use in creating the circle; the default of 0 uses
                                     the current LayoutScript default.
        :type number_of_points: int
        :return: the new circle.
        :rtype: Circle
        """
        ls_circle = self.ls.addCircle(int(layer), self.drawing._np_to_layoutscript(to_point(origin)),
                                      self.drawing.to_database_units(radius), int(number_of_points))
        return Circle(ls_circle, self.drawing)

    def add_polygon(self, points, layer):
        """Add a polygon to this cell and return the corresponding object.

        All length values are expected to be floats in user units if use_user_unit is True, and ints in database
        units if it is False. If the given list of points does not close, LayoutScript will automatically add the
        first point to the end of the point list in order to close it.

        :param points: an iterable of points that are the vertices of the polygon.
        :type points: iterable
        :param int layer: the layer on which the polygon is created.
        :return: the new polygon.
        :rtype: Polygon
        """
        ls_polygon = self.ls.addPolygon(self.drawing._to_point_array(points), int(layer))
        return Polygon(ls_polygon, self.drawing)

    def add_polygon_arc(self, center, inner_radius, outer_radius, layer, start_angle=0, stop_angle=0):
        """Add a polygon in the shape of a full or partial annulus to this cell and return the corresponding object.

        All length values are expected to be floats in user units if use_user_unit is True, and ints in database
        units if it is False. The default start and stop angles create an arc that touches itself, forming a full
        annulus. The angles are taken mod 360, so it is not possible to create a polygon that overlaps itself.

        :param center: a point containing the (x, y) coordinates of the circle center.
        :type center: numpy.ndarray
        :param inner_radius: the inner radius of the arc.
        :type inner_radius: float or int
        :param outer_radius: the outer radius of the arc.
        :type outer_radius: float or int
        :param int layer: the layer on which the arc is created.
        :param float start_angle: the start angle in degrees, measured counterclockwise from the x-axis.
        :param float stop_angle: the stop angle in degrees, measured counterclockwise from the x-axis.
        :return: the new polygon arc.
        :rtype: Polygon
        """
        ls_polygon = self.ls.addPolygonArc(self.drawing._np_to_layoutscript(to_point(center)),
                                           self.drawing.to_database_units(inner_radius),
                                           self.drawing.to_database_units(outer_radius),
                                           float(start_angle), float(stop_angle), int(layer))
        return Polygon(ls_polygon, self.drawing)

    def add_path(self, points, layer, width=None, cap=None):
        """Add a path to this cell and return the corresponding object.

        All length values are expected to be floats in user units if use_user_unit is True, and ints in database
        units if it is False. A path may be closed or open, and may have one of three styles of end cap.

        :param points: an iterable of points that are the vertices of the path.
        :type points: numpy.ndarray
        :param int layer: the layer on which the path is created.
        :param width: the width of the path; if None, use the current default, which may be zero (not recommended for
                      drawing physical structures).
        :type width: float or int
        :param int cap: the cap style of the path, where 0 means no cap, 1 means a round cap with radius equal to half
                        the path width, and 2 means a rectangular cap that extends past the end of the path by an amount
                        equal to half its width; the default None creates a path with the current default cap style.
        :return: the new path.
        :rtype: Path
        """
        ls_path = self.ls.addPath(self.drawing._to_point_array(points), int(layer))
        path = Path(ls_path, self.drawing)
        if width is not None:
            path.width = float(width)
        if cap is not None:
            path.cap = int(cap)
        return path

    def add_text(self, origin, text, layer, height=None):
        """Add text to this cell and return the corresponding object.

        :param origin: a point representing the origin of the text object, which appears to be to the upper left of
                       where the text begins.
        :type origin: numpy.ndarray
        :param str text: the text to be displayed.
        :param int layer: the layer on which to create the text.
        :param height: the height of the text; positive values are interpreted as user units, negative values create
                       a fixed height in pixels, and the default uses the current default.
        :type height: float or int
        :return: the new text.
        :rtype: Text
        """
        ls_text = self.ls.addText(int(layer), self.drawing._np_to_layoutscript(to_point(origin)), str(text))
        text_ = Text(ls_text, self.drawing)
        if height is not None:
            text_.height = height
        return text_


class Element(object):
    """Wrap a generic ``LayoutScript.element`` object.

    LayoutScript does not define subclasses for separate types of elements. Instead, the ``LayoutScript.element`` class
    has methods like ``isPolygon`` and ``isCellref`` that can be used to determine the element type.

    This class is not used directly, but its subclasses correspond to the specific types of elements that appear in
    cells in a layout, such as polygons and references to other cells. This class implements only the methods that
    are relevant for all types of ``LayoutScript.element`` objects.
    """

    def __init__(self, ls_element, drawing):
        """
        :param ls_element: the element to wrap.
        :type ls_element: LayoutScript.element
        :param drawing: the drawing in which the element exists.
        :type drawing: Drawing
        """
        self.ls = ls_element
        self.drawing = drawing

    @property
    def points(self):
        """The list of (x, y) points (``numpy.ndarray``).

        Every Element has an array of one or more points, with a size and meaning that depend on the type of element.
        For example, for a :class:`Polygon`, the array contains the points of the polygon, with identical first and
        last points. For a :class:`Cellref`, which represents a reference in one cell to another cell,
        the array contains a single point that is the location of the referenced cell.
        """
        return self.drawing._to_list_of_numpy_arrays(self.ls.getPoints())

    @points.setter
    def points(self, points):
        self.ls.setPoints(self.drawing._to_point_array(points))

    # ToDo: figure out what this means...
    @property
    def data_type(self):
        """Every Element has an integer datatype attribute, which typically equals 0."""
        return self.ls.getDatatype()

    @data_type.setter
    def data_type(self, data_type):
        self.ls.setDatatype(int(data_type))

    def __str__(self):
        """Return a str representation of the Element."""
        return '{} {}'.format(self.__class__.__name__, self.points)

    @property
    def angle(self):
        """The angle of the element in degrees, measured counterclockwise from the x-axis."""
        return self.ls.getTrans().getAngle()

    @angle.setter
    def angle(self, angle):
        transformation = self.ls.getTrans()
        # Bizarrely, strans.rotate(angle) rotates by -angle; these lines rotate to zero then to the desired angle.
        transformation.rotate(transformation.getAngle())
        transformation.rotate(-angle)
        self.ls.setTrans(transformation)

    # ToDo: explain better the rotation
    @property
    def scale(self):
        """The scale of the element, with a default of 1.0, used for expandable elements.

        The returned value is always positive. However, setting a negative scale will produce a rotation by 180
        degrees along with a scaling by the absolute value of the given scale.
        """
        return self.ls.getTrans().getScale()

    # ToDo: is this calculation correct?
    @scale.setter
    def scale(self, scale):
        transformation = self.ls.getTrans()
        transformation.scale(scale / transformation.getScale())
        self.ls.setTrans(transformation)

    @property
    def mirror_x(self):
        """A bool that controls whether the Element is reflected across the x-axis."""
        return self.ls.getTrans().getMirror_x()

    @mirror_x.setter
    def mirror_x(self, mirror):
        transformation = self.ls.getTrans()
        if bool(mirror) ^ transformation.getMirror_x():
            transformation.toggleMirror_x()
        self.ls.setTrans(transformation)

    # ToDo: verify
    def reset_transformation(self):
        """Reset the Element transformation: reset the scale to 1.0, the angle to 0.0, and mirror_x to False."""
        transformation = self.ls.getTrans()
        transformation.reset()
        self.ls.setTrans(transformation)

    # ToDo: accelerate this by using array stacking instead of iterating over the points
    def __eq__(self, other):
        return (len(self.points) == len(other.points)
                and np.all([np.all(s == o) for s, o in zip(self.points, other.points)])
                and self.data_type == other.data_type
                and self.scale == other.scale
                and self.angle == other.angle
                and self.mirror_x == other.mirror_x)


class LayerElement(Element):
    """An Element that exists on a single layer.

    This class is identical to :class:`Element` except that it adds a layer property for those that exist on a single
    layer, which is all of them except for the :class:`Cellref` and :class:`CellrefArray` classes.
    """

    @property
    def layer(self):
        """The layer number (``int``)."""
        return self.ls.layerNum

    @layer.setter
    def layer(self, layer):
        self.ls.layerNum = int(layer)


# ToDo: figure out ls_element type
class CellElement(Element):
    """This class adds a cell attribute to Element."""

    def __init__(self, ls_element, drawing):
        """
        :param ls_element:
        :type ls_element: LayoutScript.what?
        :param drawing: the drawing to which the
        :type drawing: Drawing
        """
        super(CellElement, self).__init__(ls_element, drawing)
        self.cell = Cell(ls_element.depend(), drawing)


# ToDo: finish adding comments below
class Cellref(CellElement):

    def __str__(self):
        return 'Cell {} at ({:.3f}, {:.3f})'.format(self.cell.name, self.origin[0], self.origin[1])

    @property
    def origin(self):
        return self.points[0]

    @origin.setter
    def origin(self, origin):
        self.points = [to_point(origin)]


class CellrefArray(CellElement):

    def __str__(self):
        return 'Cell {}: {} {} by {}'.format(self.cell.name, self.points, self.repeat_x, self.repeat_y)

    @staticmethod
    def _to_pylayout(points):
        origin, step_x, step_y = points
        return [origin, step_x + origin, step_y + origin]

    @staticmethod
    def _from_pylayout(points):
        origin, ls_x, ls_y = points
        return [origin, ls_x - origin, ls_y - origin]

    @property
    def points(self):
        return self._from_pylayout(self.drawing._to_list_of_numpy_arrays(self.ls.getPoints()))

    @points.setter
    def points(self, points):
        self.ls.setPoints(self.drawing._to_point_array(self._to_pylayout(points)))

    @property
    def origin(self):
        return self.points[0]

    @origin.setter
    def origin(self, origin):
        self.points = [to_point(origin), self.points[1], self.points[2]]

    @property
    def step_x(self):
        return self.points[1]

    @step_x.setter
    def step_x(self, step_x):
        self.points = [self.points[0], to_point(step_x), self.points[2]]

    @property
    def step_y(self):
        return self.points[2]

    @step_y.setter
    def step_y(self, step_y):
        self.points = [self.points[0], self.points[1], to_point(step_y)]

    @property
    def repeat_x(self):
        return self.ls.getNx()

    @repeat_x.setter
    def repeat_x(self, repeat):
        self.ls.setNx(int(repeat))

    @property
    def repeat_y(self):
        return self.ls.getNy()

    @repeat_y.setter
    def repeat_y(self, repeat):
        self.ls.setNy(int(repeat))


class Box(LayerElement):

    @property
    def _points(self):
        (x_upper_left, y_upper_left), (x_lower_right, y_lower_right) = self.points
        x = x_upper_left
        y = y_lower_right
        width = x_lower_right - x_upper_left
        height = y_upper_left - y_lower_right
        return x, y, width, height

    @property
    def x(self):
        """The x-coordinate of the origin, currently read-only."""
        return self._points[0]

    @property
    def y(self):
        """The y-coordinate of the origin, currently read-only."""
        return self._points[1]

    @property
    def width(self):
        """The horizontal extent measured from the origin, currently read-only."""
        return self._points[2]

    @property
    def height(self):
        """The vertical extent measured from the origin, currently read-only."""
        return self._points[3]

    @property
    def perimeter(self):
        """The perimeter of the box."""
        return 2 * self.width + 2 * self.height


class Circle(LayerElement):
    """LayoutEditor considers any regular polygon with more than 8 points to be a circle."""

    @property
    def center(self):
        """The center of the circle (``numpy.ndarray``, read-only)."""
        # The last point is always the same as the first.
        x = np.mean([p[0] for p in self.points[:-1]])
        y = np.mean([p[1] for p in self.points[:-1]])
        return self.drawing.from_database_units(self.drawing.to_database_units(np.array([x, y])))

    @property
    def radius(self):
        """The radius of the circle (float, read-only)."""
        return np.sqrt(np.sum((self.points[0] - self.center) ** 2))

    @property
    def perimeter(self):
        """The perimeter of the circle (float, read-only)."""
        x, y = np.vstack(self.points).T
        return np.sum(np.hypot(np.diff(x), np.diff(y)))


class Path(LayerElement):

    # ToDo: does this need to be converted to user units?
    @property
    def width(self):
        """The width of the path (float or int); note that zero-width paths are allowed."""
        return self.ls.getWidth()

    @width.setter
    def width(self, width):
        self.ls.setWidth(self.drawing.to_database_units(width))

    @property
    def cap(self):
        return self.ls.getCap()

    @cap.setter
    def cap(self, cap):
        self.ls.setCap(int(cap))

    @property
    def length(self):
        """Return the length of the path in user units, not including the caps.

        :return: the path length
        """
        x, y = np.vstack(self.points).T
        return np.sum(np.hypot(np.diff(x), np.diff(y)))


class Polygon(LayerElement):

    @property
    def perimeter(self):
        """The perimeter of the polygon.

        For a Polygon the first and last point are always the same.
        """
        x, y = np.vstack(self.points).T
        return np.sum(np.hypot(np.diff(x), np.diff(y)))


class Text(LayerElement):

    def __str__(self):
        return 'Text "{}" at ({:.3f}, {:.3f})'.format(self.text, self.origin[0], self.origin[1])

    @property
    def text(self):
        try:
            return self.ls.getName().toAscii().data()  # pylayout
        except AttributeError:
            return self.ls.getName()  # LayoutScript

    @text.setter
    def text(self, text):
        self.ls.setName(text)

    @property
    def height(self):
        return self.drawing.from_database_units(self.ls.getWidth())

    @height.setter
    def height(self, height):
        self.ls.setWidth(self.drawing.to_database_units(height))

    @property
    def origin(self):
        return self.points[0]

    @origin.setter
    def origin(self, origin):
        self.points = [to_point(origin)]
