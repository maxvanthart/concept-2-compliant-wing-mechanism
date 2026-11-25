"""Write coloured AP214 surface references for CAD, in millimetres."""
from math import pi

import numpy as np

from .fabrication import fabrication_geometry, rotate


PIVOT_RADIUS_MM = 2.0
SHEET_DEPTH_MM = 10.0
ORIENTATION_SPAN_MM = 8.0
STEP_NAMES = {'L9': 'green', 'L4': 'blue', 'L11': 'coupler', 'L0': 'ground'}


def write_step(path, lines, pivots, orientations, link_colors, support=None):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.Interface import Interface_Static
    from OCP.Quantity import Quantity_Color, Quantity_TOC_sRGB
    from OCP.STEPCAFControl import STEPCAFControl_Writer
    from OCP.STEPControl import STEPControl_AsIs
    from OCP.TCollection import TCollection_ExtendedString
    from OCP.TDataStd import TDataStd_Name
    from OCP.TDocStd import TDocStd_Document
    from OCP.XCAFDoc import XCAFDoc_ColorSurf, XCAFDoc_DocumentTool
    from OCP.gp import gp_Ax3, gp_Cylinder, gp_Dir, gp_Pnt

    doc = TDocStd_Document(TCollection_ExtendedString('Concept 2 STEP'))
    shapes = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())
    colors = XCAFDoc_DocumentTool.ColorTool_s(doc.Main())
    half_depth = SHEET_DEPTH_MM / 2

    def add(name, face, color):
        label = shapes.AddShape(face, False)
        TDataStd_Name.Set_s(label, TCollection_ExtendedString(name))
        rgb = [int(color[i:i + 2], 16) / 255 for i in (1, 3, 5)]
        colors.SetColor(label, Quantity_Color(*rgb, Quantity_TOC_sRGB), XCAFDoc_ColorSurf)

    def polygon(vertices):
        wire = BRepBuilderAPI_MakePolygon()
        for x, y, z in vertices:
            wire.Add(gp_Pnt(float(x), float(y), float(z)))
        wire.Close()
        return BRepBuilderAPI_MakeFace(wire.Wire(), True).Face()

    def sheet(start, end):
        return polygon([(*start, -half_depth), (*end, -half_depth),
                        (*end, half_depth), (*start, half_depth)])

    for name, (start, end) in lines.items():
        add(f'Link_{name}', sheet(start, end), link_colors[int(name[1:])])
    for name, (x, y) in pivots.items():
        axis = gp_Ax3(gp_Pnt(float(x), float(y), 0), gp_Dir(0, 0, 1))
        cylinder = gp_Cylinder(axis, PIVOT_RADIUS_MM)
        face = BRepBuilderAPI_MakeFace(cylinder, 0, 2 * pi, -half_depth, half_depth)
        add(f'Pivot_{name}', face.Face(), link_colors[0])
    # B's ground-side frames can coincide. Export one surface with both labels.
    unique = {}
    for name, centre, angle, owner in orientations:
        end = centre + rotate(np.array([ORIENTATION_SPAN_MM, 0]), angle)
        key = tuple(np.round(np.r_[centre, end], 8))
        if key not in unique:
            unique[key] = ([], centre, end, owner)
        unique[key][0].append(name)
    for names, start, end, owner in unique.values():
        add('Mount_' + '__'.join(names), sheet(start, end), link_colors[owner])
    if support is not None:
        add('Ground_support_profile', polygon([(*xy, 0) for xy in support['web_outline_mm']]),
            link_colors[0])
        start, end = support['floor_contact_mm']
        add('Ground_floor_contact', sheet(start, end), link_colors[0])

    writer = STEPCAFControl_Writer()
    if not Interface_Static.SetCVal_s('write.step.schema', 'AP214IS'):
        raise OSError('Could not select the AP214 STEP schema.')
    writer.SetColorMode(True)
    if not writer.Transfer(doc, STEPControl_AsIs) or writer.Write(str(path)) != IFSelect_RetDone:
        raise OSError(f'Could not write STEP surfaces to {path}')


def save_step_exports(out, xy, pos, mounts, endpoints, poses, link_colors):
    """Export the four relaxed fabrication groups."""
    _, groups = fabrication_geometry(xy, pos, mounts, endpoints, poses)
    out.mkdir(parents=True, exist_ok=True)
    for group in groups:
        root = int(group['root_body'][1:])
        pivots = {pin['pivot']: pin['center_mm'] for pin in group['pin_axis_centres']}
        orientations = []
        for joint, cfp in group['cfps'].items():
            centre = np.asarray(cfp['pivot_mm'])
            pivots[mounts[joint][0]] = centre
            orientations.append((joint + '_relaxed', centre,
                                 np.deg2rad(cfp['clamp_angle_deg']), root))
        for joint, face in group['mount_faces'].items():
            centre = np.asarray(face['pivot_mm'])
            pivots[mounts[joint][0]] = centre
            orientations.append((joint + '_' + face['body'], centre,
                                 np.deg2rad(face['face_angle_deg']), int(face['body'][1:])))
        name = STEP_NAMES[group['root_body']]
        write_step(out / f'concept2_{name}_relaxed.step', group['centerlines_mm'],
                   pivots, orientations, link_colors, group.get('support'))
