"""This module contains classes and functions are for drawing transmission lines.

Points passed to these functions and methods are converted to package format using :func:`to_point`, and sequences of
points are converted to package format using :func:`to_point_list`.
The phrase 'package format' refers to the two-element (x, y) numpy ndarray point format used everywhere in the package.
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from . import wrapper

DEFAULT_POINTS_PER_RADIAN = 60


def from_increments(increments, origin=(0, 0)):
    """Return a list of points starting from the given origin and separated by the given increments.

    This function exists because it is often easier to specify paths in terms of the differences between points than in
    terms of the absolute values. Example::

        >>> from_increments(increments=[(200, 0), (0, 300)], origin=(100, 0))
        [array([100, 0]), array([300, 0]), array([300, 300])]

    :param iterable[indexable] increments: a list of points in the package format; the differences between consecutive
                                           returned points.
    :param indexable origin: the starting point of the list.
    :return: a list of points in the package format.
    :rtype: list[numpy.ndarray]
    """
    points = [wrapper.to_point(origin)]
    for increment in [wrapper.to_point(point) for point in increments]:
        points.append(points[-1] + increment)
    return points


# ToDo: warn if consecutive points are too close together to bend properly.
def smooth(points, radius, points_per_radian=DEFAULT_POINTS_PER_RADIAN):
    """Return a list of smoothed points constructed by adding points to change the given corners into arcs.

    At each corner, the original point is replaced by points that form circular arcs that are tangent to the original
    straight sections. The returned path will not contain any of the given points except for the starting and ending
    points, because all of the interior points will be replaced.

    If the given radius is too large there is no way to make this work, and the results will be ugly. **No warning is
    currently given when this happens, so you need to inspect your design.** The given radius should be smaller than
    about half the length of the shorted straight section. If several points lie on the same line, the redundant ones
    are removed.

    :param iterable[indexable] points: a list of points in package format.
    :param float radius: the radius of the circular arcs used to connect the straight segments.
    :param int points_per_radian: the number of points per radian of arc; the default of 60 (about 1 per degree) is
                                  usually enough.
    :return: a list of smoothed points.
    :rtype: list[numpy.ndarray]
    """
    bends = []
    angles = []
    corners = []
    offsets = []
    for before, current, after in zip(points[:-2], points[1:-1], points[2:]):
        before_to_current = current - before
        current_to_after = after - current
        # The angle at which the path bends at the current point, in (-pi, pi)
        bend_angle = np.angle(np.inner(before_to_current, current_to_after) +
                              1j * np.cross(before_to_current, current_to_after))
        if np.abs(bend_angle) > 0:  # If the three points are co-linear then drop the current point
            # The distance from the corner point to the arc center
            h = radius / np.cos(bend_angle / 2)
            # The absolute angle of the arc center point, in (-pi, pi)
            theta = (np.arctan2(before_to_current[1], before_to_current[0]) +
                     bend_angle / 2 + np.sign(bend_angle) * np.pi / 2)
            # The offset of the arc center relative to the corner
            offset = h * np.array([np.cos(theta), np.sin(theta)])
            # The absolute angles of the new points (at least two), using the absolute center as origin
            arc_angles = (theta + np.pi + np.linspace(-bend_angle / 2, bend_angle / 2,
                                                      int(np.ceil(np.abs(bend_angle) * points_per_radian) + 1)))
            bend = [current + offset + radius * np.array([np.cos(phi), np.sin(phi)]) for phi in arc_angles]
            bends.append(bend)
            angles.append(bend_angle)
            corners.append(current)
            offsets.append(offset)
    return bends, angles, corners, offsets


class SegmentList(list):
    """A list subclass to hold Segments that are joined sequentially to form a path."""

    def draw(self, cell, origin, *args, **kwargs):
        """Draw all of the segments contained in this instance into the given cell.

        The segments are drawn so that the start point of the first segment is the origin, and the start of each
        subsequent segment is the end of the previous segment. If each segment has (0, 0) as its start point, they will
        be continuously connected.

        :param Cell cell: The cell into which the result is drawn.
        :param indexable origin: The point to use for the origin of the first segment.
        :param args: arguments passed to :meth:`Segment.draw`.
        :param kwargs: keyword arguments passed to :meth:`Segment.draw`.
        :return: None.
        """
        # It's crucial to avoiding input modification that this also makes a copy.
        point = wrapper.to_point(origin)
        for segment in self:
            segment.draw(cell, point, *args, **kwargs)
            # NB: using += produces an error when casting int to float.
            point = point + segment.end

    @property
    def start(self):
        """The start point of the SegmentList."""
        return self[0].start

    @property
    def end(self):
        """The end point of the SegmentList."""
        return np.sum(np.vstack([element.end for element in self]), axis=0)

    @property
    def span(self):
        """The difference between start and end points: span = end - start, in the vector sense."""
        return self.end - self.start

    @property
    def length(self):
        """The sum of the lengths of the Segments in this SegmentList."""
        return np.sum([element.length for element in self])


class Segment(object):
    """An element of a :class:`SegmentList` that can draw itself into a cell."""

    def __init__(self, points, round_to=None):
        """The given points are saved as :attr:`_points`, and should generally not be modified.

        If the first point is (0, 0) then this segment will be connected to the previous segment when drawn.

        :param iterable[indexable] points: the points that form the Segment.
        :param round_to: if not None, the coordinates of each point are rounded to this value; useful for ensuring that
                         all the points in a design lie on a grid (larger than the database unit size).
        :type round_to: float, int, or None
        """
        points = wrapper.to_point_list(points)
        if round_to is not None:
            points = [round_to * np.round(p / round_to) for p in points]
        self._points = points

    @property
    def points(self):
        """The points, rounded to :attr:`round_to` (:class:`list`[:class:`numpy.ndarray`], read-only)."""
        return self._points

    @property
    def start(self):
        """The start point (:class:`numpy.ndarray`, read-only)."""
        return self._points[0]

    @property
    def end(self):
        """The end point (:class:`numpy.ndarray`, read-only)."""
        return self._points[-1]

    @property
    def x(self):
        """The x-coordinates of all points (:class:`numpy.ndarray`, read-only)."""
        return np.array([point[0] for point in self.points])

    @property
    def y(self):
        """The y-coordinates of all points (:class:`numpy.ndarray`, read-only)."""
        return np.array([point[1] for point in self.points])

    @property
    def length(self):
        """The segment length calculated by adding the lengths of straight lines connecting the points."""
        return np.sum(np.hypot(np.diff(self.x), np.diff(self.y)))

    def draw(self, cell, origin, **kwargs):
        """Draw this instance in the given cell.

        Subclasses implement this method to draw themselves.

        :param Cell cell: the cell into which this segment will be drawn.
        :param indexable origin: draw the segment relative to this point; the meaning depends on the segment type.
        :return: None.
        """
        pass


class SmoothedSegment(Segment):
    """An element of a :class:`SegmentList` that can draw itself into a cell.

     The corners of the given outline points are smoothed by replacing them with circular arcs using :func:`smooth`.
     """

    def __init__(self, outline, radius, points_per_radian, round_to=None):
        """The given outline points are passed to :func:`smooth` and the result is stored in the instance attributes
        :attr:`bends`, :attr:`angles`, :attr:`corners`, and :attr:`offsets`.

        :param iterable[indexable] outline: the outline points, before smoothing.
        :param float radius: the radius of the circular arcs used to connect the straight segments; see :func:`smooth`.
        :param int points_per_radian: the number of points per radian of arc; see :func:`smooth`.
        :param round_to: if not None, the coordinates of the given outline points are rounded to this value before
                         smoothing, which is useful for ensuring that key points in a design lie on a grid that is
                         larger than the database unit; the smoothed points are **not** rounded and have the precision
                         of the database unit.
        :type round_to: float or None
        """
        super(SmoothedSegment, self).__init__(points=outline, round_to=round_to)
        self.radius = radius
        self.points_per_radian = points_per_radian
        self.bends, self.angles, self.corners, self.offsets = smooth(self._points, radius, points_per_radian)

    @property
    def points(self):
        """The smoothed points (:class:`list`[:class:`numpy.ndarray`]`, read-only); the outline points are
        :attr:`_points`.
        """
        p = [self.start]
        for bend in self.bends:
            p.extend(bend)
        p.append(self.end)
        return p


class Trace(SmoothedSegment):
    """A single positive wire that could be used as microstrip trace.

    It can be drawn to overlap at either end with the adjacent segments, and the overlap lengths are not counted when
    calculating the total length. This is useful, for example, when drawing a trace where the electrical connection is
    formed by an overlap.
    """

    def __init__(self, outline, width, start_overlap=0, end_overlap=0, radius=None,
                 points_per_radian=DEFAULT_POINTS_PER_RADIAN, round_to=None):
        """
        :param iterable[indexable] outline: the outline points, before smoothing.
        :param float width: the width of the trace.
        :param float start_overlap: the overlap length at the start.
        :param float end_overlap: the overlap length at the end.
        :param float radius: see :func:`smooth`.
        :param int points_per_radian: see :func:`smooth`.
        :param float round_to: see :class:`SmoothedSegment`.
        """
        self.width = width
        self.start_overlap = start_overlap
        self.end_overlap = end_overlap
        if radius is None:
            radius = 2 * width
        super(Trace, self).__init__(outline=outline, radius=radius, points_per_radian=points_per_radian,
                                    round_to=round_to)

    def draw(self, cell, origin, layer):
        """Draw this trace into the given cell.

        :param Cell cell: the cell into which to draw the trace.
        :param iterable[indexable] origin: the point at which to place the start of the trace.
        :param int layer: the layer on which to draw the trace.
        :return: None.
        """
        origin = wrapper.to_point(origin)
        points = [origin + point for point in self.points]
        cell.add_path(points=points, layer=layer, width=self.width)
        # Note that the overlap points are not stored or counted in the calculation of the length.
        if self.start_overlap > 0:
            v_start = points[0] - points[1]
            phi_start = np.arctan2(v_start[1], v_start[0])
            start_points = [points[0],
                            points[0] + self.start_overlap * np.array([np.cos(phi_start), np.sin(phi_start)])]
            cell.add_path(points=start_points, layer=layer, width=self.width)
        if self.end_overlap > 0:
            v_end = points[-1] - points[-2]
            phi_end = np.arctan2(v_end[1], v_end[0])
            end_points = [points[-1],
                          points[-1] + self.end_overlap * np.array([np.cos(phi_end), np.sin(phi_end)])]
            cell.add_path(points=end_points, layer=layer, width=self.width)
