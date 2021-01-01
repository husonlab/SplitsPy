# draw.py
"""Draws a phylogenetic outline

Draws a phylogenetic outline

See: Huson et al (2021)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""
import math
from typing import Tuple
from PIL import Image, ImageDraw, ImageFont
from splitspy.graph.graph import Graph

__author__ = 'Daniel H. Huson'


def draw(outfile: str, graph: Graph, label_angles: [float] = None, fit: float = -1.0,
         width: int = 1000, height: int = 1000,
         m_left: int = 150, m_right: int = 150, m_top: int = 150, m_bot: int = 150,
         font_size: int = 12, scale_factor: int =5) -> None:

    width *= scale_factor
    height *= scale_factor
    m_left *= scale_factor
    m_right *= scale_factor
    m_top *= scale_factor
    m_bot *= scale_factor
    font_size *= scale_factor
    line_width = scale_factor

    bw = width - m_left - m_right
    bh = height - m_top - m_bot

    if bw > bh:
        m_left += (bw - bh) / 2
        m_right += (bw - bh) / 2
        bw = bh
    elif bh > bw:
        m_top += (bh - bw) / 2
        m_bot += (bh - bw) / 2
        bh = bw

    x_min, x_max, y_min, y_max = graph.bbox()

    x = lambda a: bw / (x_max - x_min) * (a - x_min) + m_left

    y = lambda a: bh / (y_max - y_min) * (a - y_min) + m_top

    im = Image.new("RGB", (width, height), (255, 255, 255))

    im_draw = ImageDraw.Draw(im)

    font = ImageFont.truetype("Arial",size=font_size)
    black = (0, 0, 0)

    if fit != -1:
        im_draw.text((40*scale_factor, 10*scale_factor), "Fit: " + ("{:.2f}".format(fit)), font=ImageFont.truetype("Arial",size=10*scale_factor), fill=black)

    center = (0.5 * width, 0.5 * height)

    points = {}
    for v in graph.nodes():
        points[v] = (x(v.pos[0]), y(v.pos[1]))

    for e in graph.edges():
        im_draw.line([points[e.src()], points[e.tar()]], width=line_width, fill=black)

    boxes = []
    i = 1
    for v in graph.nodes():
        if v.label is not None:
            if v.label == "Root":
                angle = 90
            else:
                angle = label_angles[i] if label_angles is not None else __angle(center, points[v])
            i += 1
            pos = __label_pos(v.label, font_size, angle, points[v], boxes)
            im_draw.text(pos, v.label, font=font, fill=black)

    if outfile is None or outfile == "":
        im.show()
    else:
        im.save(outfile)


def __label_pos(label: str, font_size: int, angle: float, pt: Tuple[float,float], boxes: [[Tuple[float,float], float, float]]) -> Tuple[float,float]:
    direct = __translate((0, 0), angle)

    lw = 0.6 * font_size * len(label)
    lh = font_size
    delta = 0.2 * font_size

    if -direct[0] <= direct[1] <= direct[0]:  # right
        loc = (pt[0] + delta, pt[1] - 0.5*lh)
    elif direct[1] >= direct[0] and direct[1] >= -direct[0]:  # bottom
        loc = (pt[0] - 0.5*lw, pt[1] + 0.5*lh + delta)
    elif direct[0] <= direct[1] <= -direct[0]:  # left
        loc = (pt[0] - lw - delta, pt[1] - 0.5*lh)
    elif direct[1] <= direct[0] and direct[1] <= -direct[0]:  # top
        loc = (pt[0] - 0.5*lw, pt[1] - 1.5*lh - delta)
    else:
        loc = (0, 0)

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


def __angle(a: Tuple[float, float], b: Tuple[float, float]) -> float:
    p = (b[0] - a[0], b[1] - a[1])
    if p[0] != 0:
        a = math.atan(math.fabs(p[0]) / math.fabs(p[0]))

        if p[0] > 0:
            if p[1] > 0:
                return a / math.pi * 180
            else:
                return (2 * math.pi - a) / math.pi * 180
        else:
            if p[1] > 0:
                return (math.pi - a) / math.pi * 180
            else:
                return (math.pi + a) / math.pi * 180
    elif p[1] > 0:
        return (0.5 * math.pi) / math.pi * 180
    else:
        return (-0.5 * math.pi) / math.pi * 180


def __translate(pt: Tuple[float, float], angle: float, dist: float = 5.0) -> Tuple[float, float]:
    dx = dist * math.cos(math.pi / 180.0 * angle)
    dy = dist * math.sin(math.pi / 180.0 * angle)

    if math.fabs(dx) < 0.000001:
        dx = 0.0
    if math.fabs(dy) < 0.000001:
        dy = 0.0
    return (pt[0] + dx, pt[1] + dy)


def __intersects(box: [Tuple[float,float], float, float], other: [Tuple[float,float], float, float]) -> bool:
    if box[0][0] + box[1] < other[0][0] or other[0][0] + other[1] < box[0][0]:
        return False
    elif box[0][1] + box[2] < other[0][1] or other[0][1] + other[2] < box[0][1]:
        return False
    else:
        return True
