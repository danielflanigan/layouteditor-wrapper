"""
This module contains classes for drawing co-planar waveguide transmission lines.

All of the classes draw structures with **negative polarity**, meaning that a structure represents the absence of metal,
because this avoids the problem of having determining the extent of the ground plane.
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from layouteditor_wrapper.wrapper import to_point
from layouteditor_wrapper.transmission_line import Segment, SmoothedSegment, DEFAULT_POINTS_PER_RADIAN


class CPW(SmoothedSegment):
    """The negative space of a segment of co-planar waveguide."""

    def __init__(self, outline, width, gap, radius=None, points_per_radian=DEFAULT_POINTS_PER_RADIAN, round_to=None):
        """
        :param outline:
        :param float width: the width of the center trace.
        :param float gap: the gaps on each side of the center trace between it and the ground planes.
        :param float radius: see :func:`~wrapper.smooth`.
        :param int points_per_radian: see :func:`~wrapper.smooth`.
        :param float round_to: see :class:`SmoothedSegment`.
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
        points = [to_point(origin) + point for point in self.points]
        cell.add_path(points, negative_layer, self.width)
        cell.add_path(points, positive_layer, self.width + 2 * self.gap)
        cell.subtract(positive_layer=positive_layer, negative_layer=negative_layer, result_layer=result_layer)


class CPWBlank(SmoothedSegment):
    """Negative co-planar waveguide with the center trace missing, used when the center trace is a separate layer."""

    def __init__(self, outline, width, gap, radius=None, points_per_radian=DEFAULT_POINTS_PER_RADIAN, round_to=None):
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
        points = [to_point(origin) + point for point in self.points]
        cell.add_path(points, positive_layer, self.width + 2 * self.gap)
        cell.subtract(positive_layer=positive_layer, negative_layer=negative_layer, result_layer=result_layer)


class CPWElbowCoupler(SmoothedSegment):
    """Negative co-planar waveguide elbow coupler."""

    def __init__(self, tip_point, elbow_point, joint_point, width, gap, radius=None, points_per_radian=DEFAULT_POINTS_PER_RADIAN,
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
        points = [to_point(origin) + point for point in self.points]
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

    def __init__(self, tip_point, elbow_point, joint_point, width, gap, radius=None, points_per_radian=DEFAULT_POINTS_PER_RADIAN,
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
        points = [to_point(origin) + point for point in self.points]
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
        upper_rotated = [np.dot(rotation, to_point(p).T).T for p in upper]
        lower_rotated = [np.dot(rotation, to_point(p).T).T for p in lower]
        upper_rotated_shifted = [to_point(origin) + self.start + p for p in upper_rotated]
        lower_rotated_shifted = [to_point(origin) + self.start + p for p in lower_rotated]
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
        poly_rotated = [np.dot(rotation, to_point(p).T).T for p in poly]
        poly_rotated_shifted = [to_point(origin) + self.start + p for p in poly_rotated]
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
                 radius=None, points_per_radian=DEFAULT_POINTS_PER_RADIAN, round_to=None):
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
                 radius=None, points_per_radian=DEFAULT_POINTS_PER_RADIAN, round_to=None):
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
