"""
This module contains classes and functions that are useful for drawing co-planar waveguide components, especially
superconducting resonators.

The phrase 'package format' refers to the two-element (x, y) numpy ndarray point format used everywhere in the package.
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from . import wrapper


def from_increments(increments, origin=(0, 0)):
    """Return a list of points starting from the given origin and separated by the given increments.

    This function exists because it is often easier to specify paths in terms of the differences between points than in
    terms of the absolute values. Example::

        >>> from_increments(increments=[(200, 0), (0, 300)], origin=(100, 0))
        [np.array([100, 0]), np.array([300, 0]), np.array([300, 300])]

    :param increments: a list of points in the package format; the differences between consecutive returned points.
    :type increments: list[indexable]
    :param origin: the starting point of the list.
    :type origin: indexable
    :return: a list of points in the package format.
    :rtype: list[numpy.ndarray]
    """
    points = [wrapper.to_point(origin)]
    for increment in [wrapper.to_point(point) for point in increments]:
        points.append(points[-1] + increment)
    return points


# ToDo: warn if consecutive points are too close together to bend properly.
def smooth_path(points, radius, points_per_radian):
    """Return a list of smoothed points constructed by adding points to change the given corners into arcs.

    At each corner, points are added so that straight sections are connected by circular arcs that are tangent to the
    straight sections. If the given radius is too large there is no way to make this work, and the results will be ugly.
    The given radius should be smaller than about half the length of the shorted straight section. If several points lie
    on the same line, the redundant ones are removed. **Note that the returned path will not contain any of the given
    points except for the starting and ending points, because all of the interior points will be replaced.**

    :param points: a list of points in package format.
    :type points: list[indexable]
    :param float radius: the radius of the circular arcs used to connect the straight segments.
    :param int points_per_radian: the number of points per radian of arc radius; usually 60 (about 1 per degree) works.
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
    """A list subclass for PathElements that are joined sequentially to form a path."""

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """Draw all of the PathElements contained in this SegmentList into the given cell.

        The PathElements are drawn so that the origin of each element after the first is the end of the previous element.

        :param cell: The cell into which the result is drawn.
        :type cell: Cell
        :param origin: The point to use for the origin of the first Segment.
        :type origin: indexable
        :param int positive_layer: The positive layer for the boolean operations.
        :param int negative_layer: The negative layer for the boolean operations.
        :param int result_layer: The layer on which the final result is drawn.
        :return: None.
        """
        # It's crucial to avoiding input modification that this also makes a copy.
        point = wrapper.to_point(origin)
        for element in self:
            element.draw(cell, point, positive_layer, negative_layer, result_layer)
            # NB: using += produces an error when casting int to float.
            point = point + element.end

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
        """The sum of the lengths of the PathElements in this SegmentList."""
        return np.sum([element.length for element in self])


class Segment(object):
    """An element of a SegmentList."""

    def __init__(self, points, round_to=None):
        """

        :param points:
        :param round_to:
        """
        points = wrapper.to_point_list(points)
        if round_to is not None:
            points = [round_to * np.round(p / round_to) for p in points]
        self._points = points

    @property
    def points(self):
        """The points (``list[numpy.ndarray]``) in the element, rounded to ``round_to`` (read-only)."""
        return self._points

    @property
    def start(self):
        """The start point (``numpy.ndarray``) of the Segment (read-only)."""
        return self._points[0]

    @property
    def end(self):
        """The end point (``numpy.ndarray``) of the Segment (read-only)."""
        return self._points[-1]

    @property
    def x(self):
        """A ``numpy.ndarray`` containing the x-coordinates of all points (read-only)."""
        return np.array([point[0] for point in self.points])

    @property
    def y(self):
        """A ``numpy.ndarray`` containing the y-coordinates of all points (read-only)."""
        return np.array([point[1] for point in self.points])

    @property
    def length(self):
        """The length of the Segment, calculating by adding the lengths of straight lines connecting the points."""
        return np.sum(np.hypot(np.diff(self.x), np.diff(self.y)))

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """Draw this Segment in the given cell.

        Subclasses implement this method to draw themselves.

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :return: None
        """
        pass


class SmoothedSegment(Segment):
    """
    document me!
    """

    def __init__(self, outline, radius, points_per_radian, round_to=None):
        """

        :param outline:
        :param radius:
        :param points_per_radian:
        :param round_to:
        """
        super(SmoothedSegment, self).__init__(points=outline, round_to=round_to)
        self.radius = radius
        self.points_per_radian = points_per_radian
        self.bends, self.angles, self.corners, self.offsets = smooth_path(self._points, radius, points_per_radian)

    @property
    def points(self):
        """

        :return:
        """
        p = [self.start]
        for bend in self.bends:
            p.extend(bend)
        p.append(self.end)
        return p


class Trace(SmoothedSegment):
    """A single positive trace that could be used for microstrip.

    It can be drawn to overlap at either end with the adjacent elements, and the overlap lengths are not counted when
    calculating the total length. This is useful, for example, when drawing a trace where the electrical connection is
    formed by an overlap.
    """

    def __init__(self, outline, width, start_overlap=0, end_overlap=0, radius=None, points_per_radian=60,
                 round_to=None):
        """

        :param outline:
        :param width:
        :param start_overlap:
        :param end_overlap:
        :param radius:
        :param points_per_radian:
        :param round_to:
        """
        self.width = width
        self.start_overlap = start_overlap
        self.end_overlap = end_overlap
        if radius is None:
            radius = 2 * width
        super(Trace, self).__init__(outline=outline, radius=radius, points_per_radian=points_per_radian,
                                    round_to=round_to)

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :return:
        """
        origin = wrapper.to_point(origin)
        points = [origin + point for point in self.points]
        cell.add_path(points=points, layer=result_layer, width=self.width)
        # Note that the overlap points are not stored or counted in the calculation of the length.
        if self.start_overlap > 0:
            v_start = points[0] - points[1]
            phi_start = np.arctan2(v_start[1], v_start[0])
            start_points = [points[0],
                            points[0] + self.start_overlap * np.array([np.cos(phi_start), np.sin(phi_start)])]
            cell.add_path(points=start_points, layer=result_layer, width=self.width)
        if self.end_overlap > 0:
            v_end = points[-1] - points[-2]
            phi_end = np.arctan2(v_end[1], v_end[0])
            end_points = [points[-1],
                          points[-1] + self.end_overlap * np.array([np.cos(phi_end), np.sin(phi_end)])]
            cell.add_path(points=end_points, layer=result_layer, width=self.width)


class CPW(SmoothedSegment):
    """Negative co-planar waveguide."""

    def __init__(self, outline, width, gap, radius=None, points_per_radian=60, round_to=None):
        """

        :param outline:
        :param width:
        :param gap:
        :param radius:
        :param points_per_radian:
        :param round_to:
        """
        self.width = width
        self.gap = gap
        if radius is None:
            radius = width / 2 + gap
        super(CPW, self).__init__(outline=outline, radius=radius, points_per_radian=points_per_radian,
                                  round_to=round_to)

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :return:
        """
        points = [wrapper.to_point(origin) + point for point in self.points]
        cell.add_path(points, negative_layer, self.width)
        cell.add_path(points, positive_layer, self.width + 2 * self.gap)
        cell.subtract(positive_layer=positive_layer, negative_layer=negative_layer, result_layer=result_layer)


class CPWBlank(SmoothedSegment):
    """Negative co-planar waveguide with the center trace missing, used when the center trace is a separate layer."""

    def __init__(self, outline, width, gap, radius=None, points_per_radian=60, round_to=None):
        """

        :param outline:
        :param width:
        :param gap:
        :param radius:
        :param points_per_radian:
        :param round_to:
        """
        self.width = width
        self.gap = gap
        if radius is None:
            radius = width / 2 + gap
        super(CPWBlank, self).__init__(outline=outline, radius=radius, points_per_radian=points_per_radian,
                                       round_to=round_to)

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        points = [wrapper.to_point(origin) + point for point in self.points]
        cell.add_path(points, positive_layer, self.width + 2 * self.gap)
        cell.subtract(positive_layer=positive_layer, negative_layer=negative_layer, result_layer=result_layer)


class CPWElbowCoupler(SmoothedSegment):
    """Negative co-planar waveguide elbow coupler."""

    def __init__(self, tip_point, elbow_point, joint_point, width, gap, radius=None, points_per_radian=60,
                 round_to=None):
        """

        :param tip_point:
        :param elbow_point:
        :param joint_point:
        :param width:
        :param gap:
        :param radius:
        :param points_per_radian:
        :param round_to:
        """
        self.width = width
        self.gap = gap
        if radius is None:
            radius = width / 2 + gap
        super(CPWElbowCoupler, self).__init__(outline=[tip_point, elbow_point, joint_point], radius=radius,
                                              points_per_radian=points_per_radian, round_to=round_to)

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer, round_tip=True):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :param round_tip:
        :return:
        """
        points = [wrapper.to_point(origin) + point for point in self.points]
        cell.add_path(points, negative_layer, self.width)
        cell.add_path(points, positive_layer, self.width + 2 * self.gap)
        if round_tip:
            v = points[0] - points[1]
            theta = np.degrees(np.arctan2(v[1], v[0]))
            cell.add_polygon_arc(points[0], self.width / 2, self.width / 2 + self.gap, result_layer,
                                 theta - 90, theta + 90)
        else:
            raise NotImplementedError("Need to code this up.")
        cell.subtract(positive_layer=positive_layer, negative_layer=negative_layer, result_layer=result_layer)


class CPWElbowCouplerBlank(SmoothedSegment):
    """Negative co-planar wavaeguide elbow coupler with the center trace missing, used when the center trace is a
    separate layer.
    """

    def __init__(self, tip_point, elbow_point, joint_point, width, gap, radius=None, points_per_radian=60,
                 round_to=None):
        """

        :param tip_point:
        :param elbow_point:
        :param joint_point:
        :param width:
        :param gap:
        :param radius:
        :param points_per_radian:
        :param round_to:
        """
        self.width = width
        self.gap = gap
        if radius is None:
            radius = width / 2 + gap
        super(CPWElbowCouplerBlank, self).__init__(outline=[tip_point, elbow_point, joint_point], radius=radius,
                                                   points_per_radian=points_per_radian, round_to=round_to)

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer, round_tip=True):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :param round_tip:
        :return:
        """
        points = [wrapper.to_point(origin) + point for point in self.points]
        cell.add_path(points, result_layer, self.width + 2 * self.gap)
        if round_tip:
            v = points[0] - points[1]
            theta = np.degrees(np.arctan2(v[1], v[0]))
            cell.add_polygon_arc(points[0], 0, self.width / 2 + self.gap, result_layer, theta - 90, theta + 90)
        else:
            raise NotImplementedError("Need to code this up.")


class CPWTransition(Segment):
    """Negative transition between two sections of co-planar waveguide."""

    def __init__(self, start_point, end_point, start_width, end_width, start_gap, end_gap, round_to=None):
        """

        :param start_point:
        :param end_point:
        :param start_width:
        :param end_width:
        :param start_gap:
        :param end_gap:
        :param round_to:
        """
        super(CPWTransition, self).__init__(points=[start_point, end_point], round_to=round_to)
        self.start_width = start_width
        self.start_gap = start_gap
        self.end_width = end_width
        self.end_gap = end_gap

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :return:
        """
        v = self.end - self.start
        phi = np.arctan2(v[1], v[0])
        rotation = np.array([[np.cos(phi), -np.sin(phi)],
                             [np.sin(phi), np.cos(phi)]])
        upper = [(0, self.start_width / 2),
                 (0, self.start_width / 2 + self.start_gap),
                 (self.length, self.end_width / 2 + self.end_gap),
                 (self.length, self.end_width / 2)]
        lower = [(x, -y) for x, y in upper]
        upper_rotated = [np.dot(rotation, wrapper.to_point(p).T).T for p in upper]
        lower_rotated = [np.dot(rotation, wrapper.to_point(p).T).T for p in lower]
        upper_rotated_shifted = [wrapper.to_point(origin) + self.start + p for p in upper_rotated]
        lower_rotated_shifted = [wrapper.to_point(origin) + self.start + p for p in lower_rotated]
        cell.add_polygon(upper_rotated_shifted, result_layer)
        cell.add_polygon(lower_rotated_shifted, result_layer)


class CPWTransitionBlank(Segment):
    """Negative transition between two sections of co-planar waveguide, used when the center trace is a separate
    layer."""

    def __init__(self, start_point, end_point, start_width, end_width, start_gap, end_gap, round_to=None):
        """

        :param start_point:
        :param end_point:
        :param start_width:
        :param end_width:
        :param start_gap:
        :param end_gap:
        :param round_to:
        """
        super(CPWTransitionBlank, self).__init__(points=[start_point, end_point], round_to=round_to)
        self.start_width = start_width
        self.start_gap = start_gap
        self.end_width = end_width
        self.end_gap = end_gap

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :return:
        """
        v = self.end - self.start
        phi = np.arctan2(v[1], v[0])
        rotation = np.array([[np.cos(phi), -np.sin(phi)],
                             [np.sin(phi), np.cos(phi)]])
        poly = [(0, self.start_width / 2 + self.start_gap),
                (self.length, self.end_width / 2 + self.end_gap),
                (self.length, -self.end_width / 2 - self.end_gap),
                (0, -self.start_width / 2 - self.start_gap)]
        poly_rotated = [np.dot(rotation, wrapper.to_point(p).T).T for p in poly]
        poly_rotated_shifted = [wrapper.to_point(origin) + self.start + p for p in poly_rotated]
        cell.add_polygon(poly_rotated_shifted, result_layer)


# ToDo: split this into separate classes
# ToDo: refactor this so that it's clear that how to mix it with Segment subclasses
class Mesh(object):
    """A mix-in class that allows Segment subclasses that have the same outlines to share mesh code."""

    def path_mesh(self):
        """Return the centers of the mesh elements calculated using the Segment's parameters.

        :return: the centers of the mesh elements.
        :rtype: list[numpy.ndarray]
        """
        mesh_centers = []
        center_to_first_row = self.width / 2 + self.gap + self.mesh_border
        # Mesh the straight sections
        starts = [self.start] + [bend[-1] for bend in self.bends]
        ends = [bend[0] for bend in self.bends] + [self.end]
        for start, end in zip(starts, ends):
            v = end - start
            length = np.linalg.norm(v)
            phi = np.arctan2(v[1], v[0])
            R = np.array([[np.cos(phi), -np.sin(phi)],
                          [np.sin(phi), np.cos(phi)]])
            num_mesh_columns = np.floor(length / self.mesh_spacing)
            if num_mesh_columns == 0:
                continue
            elif num_mesh_columns == 1:
                x = np.array([length / 2])
            else:
                x = np.linspace(self.mesh_spacing / 2, length - self.mesh_spacing / 2, num_mesh_columns)
            y = center_to_first_row + self.mesh_spacing * np.arange(self.num_mesh_rows)
            xx, yy = np.meshgrid(np.concatenate((x, x)), np.concatenate((y, -y)))
            Rxy = np.dot(R, np.vstack((xx.flatten(), yy.flatten())))
            mesh_centers.extend(zip(start[0] + Rxy[0, :], start[1] + Rxy[1, :]))
        # Mesh the curved sections
        for row in range(self.num_mesh_rows):
            center_to_row = center_to_first_row + row * self.mesh_spacing
            for radius in [self.radius - center_to_row, self.radius + center_to_row]:
                if radius < self.mesh_spacing / 2:
                    continue
                for angle, corner, offset in zip(self.angles, self.corners, self.offsets):
                    num_points = np.round(radius * np.abs(angle) / self.mesh_spacing)
                    if num_points == 1:
                        max_angle = 0
                    else:
                        max_angle = (1 - 1 / num_points) * angle / 2
                    mesh_centers.extend([corner + offset +
                                         radius * np.array([np.cos(phi), np.sin(phi)]) for phi in
                                         (np.arctan2(-offset[1], -offset[0]) +
                                          np.linspace(-max_angle, max_angle, num_points))])
        return mesh_centers

    def trapezoid_mesh(self):
        """

        :return:
        """
        mesh_centers = []
        v = self.end - self.start
        length = np.linalg.norm(v)
        phi = np.arctan2(v[1], v[0])
        R = np.array([[np.cos(phi), -np.sin(phi)],
                      [np.sin(phi), np.cos(phi)]])
        start_to_first_row = self.start_width / 2 + self.start_gap + self.start_mesh_border
        difference_to_first_row = self.end_width / 2 + self.end_gap + self.end_mesh_border - start_to_first_row
        num_mesh_columns = np.floor(length / self.mesh_spacing)
        if num_mesh_columns == 0:
            return []
        elif num_mesh_columns == 1:
            x = np.array([length / 2])
        else:
            x = np.linspace(self.mesh_spacing / 2, length - self.mesh_spacing / 2, num_mesh_columns)
        y = self.mesh_spacing * np.arange(self.num_mesh_rows, dtype=np.float)
        xxp, yyp = np.meshgrid(x, y)  # These correspond to the positive y-values
        y_shift = start_to_first_row + difference_to_first_row * x / length
        yyp += y_shift
        xx = np.concatenate((xxp, xxp))
        yy = np.concatenate((yyp, -yyp))  # The negative y-values are reflected
        Rxy = np.dot(R, np.vstack((xx.flatten(), yy.flatten())))
        mesh_centers.extend(zip(self.start[0] + Rxy[0, :], self.start[1] + Rxy[1, :]))
        return mesh_centers


class CPWMesh(CPW, Mesh):
    """doc me!"""

    def __init__(self, outline, width, gap, mesh_spacing, mesh_border, mesh_radius, num_circle_points, num_mesh_rows,
                 radius=None, points_per_radian=60, round_to=None):
        """

        :param outline:
        :param width:
        :param gap:
        :param mesh_spacing:
        :param mesh_border:
        :param mesh_radius:
        :param num_circle_points:
        :param num_mesh_rows:
        :param radius:
        :param points_per_radian:
        :param round_to:
        """
        super(CPWMesh, self).__init__(outline=outline, width=width, gap=gap, radius=radius,
                                      points_per_radian=points_per_radian, round_to=round_to)
        self.mesh_spacing = mesh_spacing
        self.mesh_radius = mesh_radius
        self.mesh_border = mesh_border
        self.num_circle_points = num_circle_points
        self.num_mesh_rows = num_mesh_rows
        self.mesh_centers = self.path_mesh()

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :return:
        """
        super(CPWMesh, self).draw(cell=cell, origin=origin, positive_layer=positive_layer,
                                  negative_layer=negative_layer, result_layer=result_layer)
        for mesh_center in self.mesh_centers:
            cell.add_circle(origin=origin + mesh_center, radius=self.mesh_radius, layer=result_layer,
                            number_of_points=self.num_circle_points)


class CPWBlankMesh(CPWBlank, Mesh):
    """doc me!"""

    def __init__(self, outline, width, gap, mesh_spacing, mesh_border, mesh_radius, num_circle_points, num_mesh_rows,
                 radius=None, points_per_radian=60, round_to=None):
        """

        :param outline:
        :param width:
        :param gap:
        :param mesh_spacing:
        :param mesh_border:
        :param mesh_radius:
        :param num_circle_points:
        :param num_mesh_rows:
        :param radius:
        :param points_per_radian:
        :param round_to:
        """
        super(CPWBlankMesh, self).__init__(outline=outline, width=width, gap=gap, radius=radius,
                                           points_per_radian=points_per_radian, round_to=round_to)
        self.mesh_spacing = mesh_spacing
        self.mesh_radius = mesh_radius
        self.mesh_border = mesh_border
        self.num_circle_points = num_circle_points
        self.num_mesh_rows = num_mesh_rows
        self.mesh_centers = self.path_mesh()

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :return:
        """
        super(CPWBlankMesh, self).draw(cell=cell, origin=origin, positive_layer=positive_layer,
                                       negative_layer=negative_layer, result_layer=result_layer)
        for mesh_center in self.mesh_centers:
            cell.add_circle(origin=origin + mesh_center, radius=self.mesh_radius, layer=result_layer,
                            number_of_points=self.num_circle_points)


class CPWElbowCouplerMesh(CPWElbowCoupler):
    """doc me!"""
    pass


class CPWElbowCouplerBlankMesh(CPWElbowCouplerBlank):
    """doc me!"""
    pass


class CPWTransitionMesh(CPWTransition, Mesh):
    """doc me!"""

    def __init__(self, start_point, end_point, start_width, end_width, start_gap, end_gap, mesh_spacing,
                 start_mesh_border, end_mesh_border, mesh_radius, num_circle_points, num_mesh_rows, round_to=None):
        """

        :param start_point:
        :param end_point:
        :param start_width:
        :param end_width:
        :param start_gap:
        :param end_gap:
        :param mesh_spacing:
        :param start_mesh_border:
        :param end_mesh_border:
        :param mesh_radius:
        :param num_circle_points:
        :param num_mesh_rows:
        :param round_to:
        """
        super(CPWTransitionMesh, self).__init__(start_point=start_point, end_point=end_point, start_width=start_width,
                                                end_width=end_width, start_gap=start_gap, end_gap=end_gap,
                                                round_to=round_to)
        self.mesh_spacing = mesh_spacing
        self.start_mesh_border = start_mesh_border
        self.end_mesh_border = end_mesh_border
        self.mesh_radius = mesh_radius
        self.num_circle_points = num_circle_points
        self.num_mesh_rows = num_mesh_rows
        self.mesh_centers = self.trapezoid_mesh()

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :return:
        """
        super(CPWTransitionMesh, self).draw(cell=cell, origin=origin, positive_layer=positive_layer,
                                            negative_layer=negative_layer, result_layer=result_layer)
        for mesh_center in self.mesh_centers:
            cell.add_circle(origin=origin + mesh_center, radius=self.mesh_radius, layer=result_layer,
                            number_of_points=self.num_circle_points)


class CPWTransitionBlankMesh(CPWTransitionBlank, Mesh):
    """doc me!"""

    def __init__(self, start_point, end_point, start_width, end_width, start_gap, end_gap, mesh_spacing,
                 start_mesh_border, end_mesh_border, mesh_radius, num_circle_points, num_mesh_rows, round_to=None):
        """

        :param start_point:
        :param end_point:
        :param start_width:
        :param end_width:
        :param start_gap:
        :param end_gap:
        :param mesh_spacing:
        :param start_mesh_border:
        :param end_mesh_border:
        :param mesh_radius:
        :param num_circle_points:
        :param num_mesh_rows:
        :param round_to:
        """
        super(CPWTransitionBlankMesh, self).__init__(start_point=start_point, end_point=end_point,
                                                     start_width=start_width, end_width=end_width,
                                                     start_gap=start_gap, end_gap=end_gap, round_to=round_to)
        self.mesh_spacing = mesh_spacing
        self.start_mesh_border = start_mesh_border
        self.end_mesh_border = end_mesh_border
        self.mesh_radius = mesh_radius
        self.num_circle_points = num_circle_points
        self.num_mesh_rows = num_mesh_rows
        self.mesh_centers = self.trapezoid_mesh()

    def draw(self, cell, origin, positive_layer, negative_layer, result_layer):
        """

        :param cell:
        :param origin:
        :param positive_layer:
        :param negative_layer:
        :param result_layer:
        :return:
        """
        super(CPWTransitionBlankMesh, self).draw(cell=cell, origin=origin, positive_layer=positive_layer,
                                                 negative_layer=negative_layer, result_layer=result_layer)
        for mesh_center in self.mesh_centers:
            cell.add_circle(origin=origin + mesh_center, radius=self.mesh_radius, layer=result_layer,
                            number_of_points=self.num_circle_points)
