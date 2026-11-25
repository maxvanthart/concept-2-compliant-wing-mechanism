"""Clamp frames and separated, body-local fabrication drawings."""
import json

import numpy as np

from .config import (CFP_BLADE_WIDTH_MM, CFP_CROSSING_DEG, CFP_LENGTH_MM,
                     CFP_THICKNESS_MM)
from .prbm import FLEXURE_PAIRS


CFP_WIDTH_MM = CFP_LENGTH_MM * np.sin(np.deg2rad(CFP_CROSSING_DEG / 2))
CFP_HEIGHT_MM = CFP_LENGTH_MM * np.cos(np.deg2rad(CFP_CROSSING_DEG / 2))
NEUTRAL_CLAMPS = np.array([
    [[-CFP_WIDTH_MM / 2, -CFP_HEIGHT_MM / 2], [CFP_WIDTH_MM / 2, -CFP_HEIGHT_MM / 2]],
    [[-CFP_WIDTH_MM / 2, CFP_HEIGHT_MM / 2], [CFP_WIDTH_MM / 2, CFP_HEIGHT_MM / 2]],
])
BODY_NAMES = {0: 'Ground', 4: 'Triangle 1', 9: 'Triangle 2', 11: 'Wing coupler',
              1: 'AP rocker', 5: 'BQ link', 10: 'BD link', 6: 'AC crank'}
BODY_LINKS = {body: (body,) for body in BODY_NAMES}
BODY_LINKS.update({4: (2, 3, 4), 9: (7, 8, 9)})
REVOLUTE_INTERFACES = (('A', 0, 1), ('A', 0, 6), ('C', 6, 9))
PRINT_GROUPS = (
    (9, 'Green assembly', ('D',), ('B_BD',)),
    (4, 'Blue assembly', ('P', 'Q'), ('B_BQ',)),
    (11, 'Red wing coupler', (), ('T1', 'T2')),
    (0, 'Rigid ground', (), ()),
)


def rotate(points, angle):
    c, s = np.cos(angle), np.sin(angle)
    return points @ np.array([[c, s], [-s, c]])


def clamp_offsets(joint, neutral_angle):
    """Return parent and child clamp-line offsets from their directed body axes."""
    if joint in ('T1', 'T2'):
        # The moving clamp centres face into the coupler from its opposite ends.
        child_offset = np.deg2rad(-90 if joint == 'T1' else 90)
        return child_offset + neutral_angle, child_offset
    # The parent clamp is on the inboard side of the parent bar's end.
    parent_offset = -np.pi / 2
    return parent_offset, parent_offset - neutral_angle


def relaxed_frames(root, internal_joints, parts, mounts):
    """Place the group's rigid bodies by releasing its internal CFPs."""
    frames = {root: (np.zeros(2), 0.0)}
    for joint in internal_joints:
        parent, child = FLEXURE_PAIRS[joint]
        known, new = (parent, child) if parent in frames else (child, parent)
        _, _, neutral = mounts[joint]
        shift, angle = frames[known]
        new_angle = angle + (neutral if new == child else -neutral)
        known_pivot = np.array(parts[f'L{known}']['cfp_mounts'][joint]['pivot_mm'])
        new_pivot = np.array(parts[f'L{new}']['cfp_mounts'][joint]['pivot_mm'])
        new_shift = shift + rotate(known_pivot, angle) - rotate(new_pivot, new_angle)
        frames[new] = (new_shift, new_angle)
    return frames


def body_local_geometry(xy, pos, mounts, link_endpoints):
    parts = {}
    for body, name in BODY_NAMES.items():
        origin = xy[link_endpoints[body][0]]
        angle = pos[f'theta{body}']

        def local(points):
            return rotate(points - origin, -angle)

        part = {'name': name, 'origin_pivot': link_endpoints[body][0],
                'centerlines': {}, 'cfp_mounts': {}, 'pin_axis_centres': []}
        parts[f'L{body}'] = part
        for link in BODY_LINKS[body]:
            endpoints = local(np.array([xy[p] for p in link_endpoints[link]]))
            part['centerlines'][f'L{link}'] = endpoints.tolist()
        for pivot, parent, child in REVOLUTE_INTERFACES:
            if body in (parent, child):
                centre = local(xy[pivot])
                mate = child if body == parent else parent
                part['pin_axis_centres'].append({'pivot': pivot, 'mate': f'L{mate}',
                                                 'center_mm': centre.tolist()})
        for joint, (pivot, offsets, _) in mounts.items():
            pair = FLEXURE_PAIRS[joint]
            if body not in pair:
                continue
            side = pair.index(body)
            centre = local(xy[pivot])
            clamp = rotate(NEUTRAL_CLAMPS[side], offsets[side]) + centre
            part['cfp_mounts'][joint] = {
                'side': ('parent', 'child')[side], 'pivot_mm': centre.tolist(),
                'clamp_angle_deg': (float(np.rad2deg(offsets[side])) + 180) % 360 - 180,
                'blade_endpoints_mm': clamp.tolist(),
            }
    return parts


def print_group(root, name, internal, carried, parts, mounts, *, root_angle=0.0):
    frames = relaxed_frames(root, internal, parts, mounts)
    group = {'name': name, 'root_body': f'L{root}', 'body_transforms': {},
             'centerlines_mm': {}, 'pin_axis_centres': [], 'cfps': {}, 'mount_faces': {}}
    for body, (shift, angle) in frames.items():
        shift = rotate(shift, root_angle)
        angle += root_angle
        group['body_transforms'][f'L{body}'] = {
            'origin_mm': shift.tolist(), 'angle_deg': float(np.rad2deg(angle)),
        }
        part = parts[f'L{body}']
        for link, ends in part['centerlines'].items():
            group['centerlines_mm'][link] = (rotate(ends, angle) + shift).tolist()
        for pin in part['pin_axis_centres']:
            centre = rotate(pin['center_mm'], angle) + shift
            group['pin_axis_centres'].append({**pin, 'center_mm': centre.tolist(), 'body': f'L{body}'})
        for joint, face in part['cfp_mounts'].items():
            pivot = rotate(face['pivot_mm'], angle) + shift
            clamp_angle = angle + np.deg2rad(face['clamp_angle_deg'])
            side = ('parent', 'child').index(face['side'])
            if joint in internal or joint in carried:
                # The child and parent define the same frame at zero twist.
                group['cfps'][joint] = {
                    'attachment': 'internal' if joint in internal else face['side'],
                    'pivot_mm': pivot.tolist(),
                    'clamp_angle_deg': float((np.rad2deg(clamp_angle) + 180) % 360 - 180),
                    'clamps_mm': (rotate(NEUTRAL_CLAMPS, clamp_angle) + pivot).tolist(),
                }
            else:
                ends = rotate(face['blade_endpoints_mm'], angle) + shift
                normal = rotate(np.array([0, 1 - 2 * side]), clamp_angle)
                group['mount_faces'][joint] = {
                    'body': f'L{body}', 'pivot_mm': pivot.tolist(),
                    'face_line_mm': ends.tolist(), 'towards_cfp_xy': normal.tolist(),
                    'face_angle_deg': float((np.rad2deg(clamp_angle) + 180) % 360 - 180),
                }
    return group


def ground_support(group, poses, mounts):
    """Define a side-view support concept with its foot below the nominal sweep."""
    a, b = np.asarray(group['centerlines_mm']['L0'])
    lowest = min((pose[point][1] - pose['A'][1]) * 1000
                 for pose in poses for point in ('A', 'B', 'C', 'D', 'P', 'Q', 'T1', 'T2'))
    # Bound the nominal in-plane CFP footprint about every virtual pivot.
    footprint_radius = np.hypot(CFP_WIDTH_MM, CFP_HEIGHT_MM) / 2
    lowest_mount = min((pose[pivot][1] - pose['A'][1]) * 1000 - footprint_radius
                       for pose in poses for pivot, _, _ in mounts.values())
    foot_top = min(lowest, lowest_mount) - 10
    floor = foot_top - 6
    left, right = sorted((a[0], b[0]))
    upper = sorted((a, b), key=lambda point: point[0])
    return {
        'view': 'Side elevation projected onto linkage plane',
        'layout': 'Support web behind the moving links',
        'web_outline_mm': [[left - 5, foot_top], [right + 5, foot_top],
                           (upper[1] + [5, 5]).tolist(), (upper[0] + [-5, 5]).tolist()],
        'foot_outline_mm': [[left - 12, floor], [right + 12, floor],
                            [right + 12, foot_top], [left - 12, foot_top]],
        'floor_contact_mm': [[left - 12, floor], [right + 12, floor]],
    }


def draw_print_group(axis, group, color):
    axis.set_title(group['name'], color=color, loc='left', fontweight='bold')
    if 'support' in group:
        support = group['support']
        for outline, alpha in (('web_outline_mm', .18), ('foot_outline_mm', .30)):
            axis.fill(*np.asarray(support[outline]).T, color=color, alpha=alpha)
        contact = np.asarray(support['floor_contact_mm'])
        axis.plot(*contact.T, color=color, lw=2)
        axis.annotate('Flat foot', contact.mean(axis=0), xytext=(0, -12),
                      textcoords='offset points', ha='center', fontsize=9)
        slope = (group['body_transforms']['L0']['angle_deg'] + 180) % 360 - 180
        axis.annotate(f'AB {slope:+.2f}°', np.mean(group['centerlines_mm']['L0'], axis=0),
                      xytext=(0, 7), textcoords='offset points', ha='center', fontsize=9)
    for ends in group['centerlines_mm'].values():
        axis.plot(*np.asarray(ends).T, color=color, lw=2, alpha=.48)
    pin_labels = {pin['pivot']: pin['center_mm'] for pin in group['pin_axis_centres']}
    for pivot, centre in pin_labels.items():
        axis.plot(*centre, 'o', color=color, mfc='white', ms=7)
        axis.annotate(f'{pivot} pin', centre, xytext=(-8, 8), textcoords='offset points',
                      ha='right', fontsize=9)
    for joint, cfp in group['cfps'].items():
        clamps = np.asarray(cfp['clamps_mm'])
        for i in (0, 1):
            blade = np.array([clamps[0, i], clamps[1, 1 - i]])
            axis.plot(*blade.T, color=color, lw=1.35)
        for ends in clamps:
            axis.plot(*ends.T, color=color, lw=2.4)
        axis.annotate(joint.replace('_', '–'), clamps.mean(axis=(0, 1)),
                      xytext=(8, 9), textcoords='offset points', fontsize=9)
    for joint, face in group['mount_faces'].items():
        ends = np.asarray(face['face_line_mm'])
        axis.plot(*ends.T, color=color, ls='--', marker='s', mfc='white', lw=2, ms=4)
        centre = ends.mean(axis=0)
        axis.annotate('', centre + 6 * np.asarray(face['towards_cfp_xy']), xytext=centre,
                      arrowprops={'arrowstyle': '->', 'color': color})
        if joint != 'B_BD':
            label = 'B · two mount faces' if joint == 'B_BQ' else f'{joint} mount face'
            axis.annotate(label, centre, xytext=(8, -13), textcoords='offset points', fontsize=9)


def fabrication_geometry(xy, pos, mounts, link_endpoints, poses):
    parts = body_local_geometry(xy, pos, mounts, link_endpoints)
    groups = [print_group(*spec, parts, mounts, root_angle=pos['theta0'] if spec[0] == 0 else 0)
              for spec in PRINT_GROUPS]
    ground = groups[-1]
    ground['support'] = ground_support(ground, poses, mounts)
    return parts, groups


def save_fabrication_layout(out, xy, pos, mounts, link_endpoints, colors, poses):
    """Draw four printable groups with all included blades at their free shape."""
    import matplotlib.pyplot as plt

    parts, groups = fabrication_geometry(xy, pos, mounts, link_endpoints, poses)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained',
                            gridspec_kw={'height_ratios': [2.7, 1]})
    fig.suptitle('Concept 2 | Four relaxed fabrication groups\n'
                 'Solid crosses: neutral CFPs   ·   Dashed lines and arrows: mating faces',
                 fontsize=15, fontweight='bold')
    panel_bounds = []
    for axis, group, (root, *_rest) in zip(axes.flat, groups, PRINT_GROUPS):
        draw_print_group(axis, group, colors[root])
        limits = axis.dataLim.get_points()
        panel_bounds.append(limits)
    bounds = np.stack(panel_bounds)
    spans = bounds[:, 1] - bounds[:, 0]
    half_width = (spans[:, 0].max() + 64) / 2
    row_heights = [max(spans[:2, 1]) + 60, max(spans[2:, 1]) + 60]
    axes[0, 0].get_gridspec().set_height_ratios(row_heights)
    for i, (axis, limits) in enumerate(zip(axes.flat, bounds)):
        x, y = limits.mean(axis=0)
        half_height = row_heights[i // 2] / 2
        axis.set(xlim=(x - half_width, x + half_width), ylim=(y - half_height, y + half_height),
                 aspect='equal')
        axis.set_box_aspect(half_height / half_width)
        if i % 2 == 0:
            axis.set_ylabel('y [mm]')
        if i >= 2:
            axis.set_xlabel('x [mm]')
        axis.grid(color='#e4eaf0', zorder=0)
    fig.savefig(out / 'prbm_fabrication_layout.svg')
    fig.savefig(out / 'prbm_fabrication_layout.png', dpi=160)
    plt.close(fig)
    data = {'units': 'mm', 'body_frames': 'Origin at named pivot, +x along directed L axis',
            'pin_axes': 'Perpendicular to the drawing plane',
            'cfp': {'free_blade_length_mm': CFP_LENGTH_MM, 'crossing_deg': CFP_CROSSING_DEG,
                    'thickness_mm': CFP_THICKNESS_MM, 'blade_depth_mm': CFP_BLADE_WIDTH_MM,
                    'neutral_clamps_parent_then_child_mm': NEUTRAL_CLAMPS.tolist()},
            'bodies': parts, 'fabrication_groups': groups,
            'separate_input_crank': 'L6'}
    (out / 'prbm_fabrication_frames.json').write_text(json.dumps(data, indent=2) + '\n',
                                                     encoding='utf-8')
