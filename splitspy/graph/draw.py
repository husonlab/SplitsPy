# draw.py
"""Draws a phylogenetic outline

Draws a phylogenetic outline, using John Zelle's graphics library

See: Huson et al (2021)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""
import math

from splitspy.graph.graph import Graph
from .graphics import *

__author__ = 'Daniel H. Huson'


def draw(graph: Graph, label_angles: [float] = None, fit: float = -1.0, title: str = "Outline",
         width: int = 1000, height: int = 1000,
         m_left: int = 150, m_right: int = 150, m_top: int = 150, m_bot: int = 150,
         font_size: int = 12) -> None:

    bw = width - m_left - m_right
    bh = height - m_top - m_bot

    if bw > bh:
        m_left += (bw - bh)/2
        m_right += (bw - bh)/2
        bw = bh
    elif bh > bw:
        m_top += (bh - bw) / 2
        m_bot += (bh - bw) / 2
        bh = bw

    x_min, x_max, y_min, y_max = graph.bbox()

    x = lambda a: bw / (x_max - x_min) * (a - x_min) + m_left

    y = lambda a: bh / (y_max - y_min) * (a - y_min) + m_top

    win = GraphWin(title, width, height)
    win.setBackground('white')

    if fit != -1:
        Text(Point(40, 10), "Fit: " + ("{:.2f}".format(fit))).draw(win)

    center = Point(0.5 * width, 0.5 * height)

    points = {}
    for v in graph.nodes():
        points[v] = Point(x(v.pos[0]), y(v.pos[1]))

    for e in graph.edges():
        ln = Line(points[e.src()], points[e.tar()])
        ln.draw(win)

    boxes = []
    i = 1
    for v in graph.nodes():
        # Circle(points[v], 0.5).draw(win)
        if v.label is not None:
            if v.label == "Root":
                angle = 90
            else:
                angle = label_angles[i] if label_angles is not None else __angle(center,points[v])
            i += 1
            pos = __label_pos(v.label, font_size, angle, points[v], boxes)
            label = Text(pos, v.label)
            label.setSize(font_size)
            label.draw(win)

    win.wait_window()


def __label_pos(label: str, font_size: int, angle: float, pt: Point, boxes: [[Point, float, float]]) -> Point:
    direct = __translate(Point(0, 0), angle)

    lw = 0.6 * font_size * len(label)
    lh = font_size
    delta = 0.2 * font_size

    if -direct.x <= direct.y <= direct.x:  # right
        loc = Point(pt.x + delta + 0.5 * lw, pt.y)
    elif direct.y >= direct.x and direct.y >= -direct.x:  # bottom
        loc = Point(pt.x, pt.y + lh + delta)
    elif direct.x <= direct.y <= -direct.x:  # left
        loc = Point(pt.x - 0.5 * lw - delta, pt.y)
    elif direct.y <= direct.x and direct.y <= -direct.x:  # top
        loc = Point(pt.x, pt.y - lh - delta)
    else:
        loc = Point(0, 0)

    box = [loc, lw, lh]
    changed = True
    while changed:
        changed = False
        for other in boxes:
            while __intersects(box, other):
                box[0] = __translate(box[0], angle, 5)
                changed = True
            if changed:
                break

    boxes.append(box)

    return box[0]


def __angle(a: Point, b: Point) -> float:
    p = Point(b.x-a.x, b.y-a.y)
    if p.x != 0:
        a = math.atan(math.fabs(p.x)/math.fabs(p.x))

        if p.x >0:
            if p.y >0:
                return a/math.pi*180
            else:
                return (2*math.pi-a)/math.pi*180
        else:
            if p.y >0:
                return (math.pi-a)/math.pi*180
            else:
                return (math.pi+a)/math.pi*180
    elif p.y > 0:
        return (0.5*math.pi) / math.pi * 180
    else:
        return (-0.5 * math.pi) / math.pi * 180


def __translate (pt: Point, angle: float, dist: float = 5.0) -> Point:
    dx = dist * math.cos(math.pi / 180.0 * angle)
    dy = dist * math.sin(math.pi / 180.0 * angle)

    if math.fabs(dx) < 0.000001:
        dx = 0.0
    if math.fabs(dy) < 0.000001:
        dy = 0.0
    return Point(pt.x + dx, pt.y + dy)


def __intersects(box: [Point,float, float], other: [Point, float, float]) -> bool:

    if box[0].x+box[1] < other[0].x or other[0].x+other[1] < box[0].x:
        return False
    elif box[0].y+box[2] < other[0].y or other[0].y+other[2] < box[0].y:
        return False
    else:
        return True

