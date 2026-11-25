"""Draw the nominal linkage, revolute joints and candidate torsion springs."""
import csv

import numpy as np

from .config import (CFP_BLADE_WIDTH_MM, CFP_CROSSING_DEG, CFP_LENGTH_MM,
                     CFP_THICKNESS_MM, ELEMENT_WIDTH_M)
from .fabrication import (CFP_HEIGHT_MM, CFP_WIDTH_MM, NEUTRAL_CLAMPS,
                          clamp_offsets, rotate, save_fabrication_layout)
from .geometry import sweep_driven_positions
from .prbm import FLEXURE_PAIRS, relative_flexure_angles


LINK_ENDPOINTS = (
    ('A', 'B'), ('A', 'P'), ('P', 'T1'), ('Q', 'T1'), ('P', 'Q'), ('B', 'Q'),
    ('A', 'C'), ('C', 'T2'), ('D', 'T2'), ('C', 'D'), ('B', 'D'), ('T1', 'T2'),
)
PIVOT_INTERFACES = (
    ('A', 'Ground to AP and AC', None),
    ('B', 'Ground to BQ', 'B_BQ'),
    ('B', 'Ground to BD', 'B_BD'),
    ('P', 'AP to triangle 1', 'P'),
    ('Q', 'BQ to triangle 1', 'Q'),
    ('C', 'AC to triangle 2', None),
    ('D', 'BD to triangle 2', 'D'),
    ('T1', 'Triangle 1 to L11', 'T1'),
    ('T2', 'Triangle 2 to L11', 'T2'),
)
REVOLUTE_PIVOTS = frozenset({'A', 'C'})
NEUTRAL_AXIS_COLOR = '#c21870'
LINK_COLORS = ('#526174',) + ('#23649a',) * 5 + ('#177965',) * 5 + ('#a8402c',)


def neutral_angles(ranges):
    for pivot, _, spring in PIVOT_INTERFACES:
        if spring not in ranges:
            continue
        parent, child = FLEXURE_PAIRS[spring]
        low, high = ranges[spring]
        unwrapped = (low + high) / 2
        wrapped = (unwrapped + 180) % 360 - 180
        label = f'{pivot} ({spring[2:]})' if pivot == 'B' else pivot
        yield pivot, label, parent, child, wrapped, unwrapped


def overlay_neutral_angles(ax, xy, pos, ranges, length_mm):
    """Place neutral-angle arcs at the virtual pivots, using assembly parent axes."""
    label_offsets = {
        'B (BQ)': (47, -1), 'B (BD)': (47, -21),
        'P': (-19, -26), 'Q': (26, 3), 'D': (25, -19),
        'T1': (-21, 27), 'T2': (27, -22),
    }
    color = '#5b287e'
    for pivot, label, parent, _, wrapped, _ in neutral_angles(ranges):
        origin = xy[pivot]
        start = pos[f'theta{parent}']
        end = start + np.deg2rad(wrapped)
        radius = length_mm * (0.102 if label == 'B (BQ)' else 0.065)
        angles = np.linspace(start, end, 50)
        arc = origin + radius * np.column_stack((np.cos(angles), np.sin(angles)))
        ax.plot(*arc.T, color=color, linewidth=1.35, zorder=4)
        for angle, style, ray_color, width in (
            (start, ':', color, 1.0),
            (end, '--', NEUTRAL_AXIS_COLOR, 2.0),
        ):
            ray_end = origin + 1.28 * radius * np.array([np.cos(angle), np.sin(angle)])
            ax.plot(*np.stack((origin, ray_end)).T, style, color=ray_color, linewidth=width,
                    zorder=3)
        offset = label_offsets[label]
        caption = f'{wrapped:+.2f}°'
        if pivot == 'B':
            caption = label.removeprefix('B (').removesuffix(')') + ' ' + caption
        ax.annotate(caption, xy=arc[len(arc) // 2], xytext=offset,
                    textcoords='offset points', xycoords='data',
                    ha='left' if offset[0] > 0 else 'right', va='center',
                    color=color, fontsize=8.5, fontweight='bold', zorder=6,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1.3},
                    arrowprops={'arrowstyle': '-', 'color': color, 'lw': 0.7})


def draw_neutral_angles(info, ranges):
    """Draw each neutral orientation in its own parent body frame."""
    info.clear()
    info.set_axis_off()
    if not ranges:
        info.text(0, 0.9, 'Neutral angles require a closed AC revolution.', fontsize=12)
        return

    parent_color, child_color = '#23649a', NEUTRAL_AXIS_COLOR
    rows = list(neutral_angles(ranges))

    info.text(0, 0.965, 'Relative neutral angles', fontsize=13, fontweight='bold')
    info.text(0, 0.925, 'Positive angles rotate counterclockwise from parent to child.',
              fontsize=9, color='#526174')
    table = info.table(
        cellText=[(label, f'L{parent} → L{child}', f'{wrapped:+.2f}°', f'{unwrapped:+.2f}°')
                  for _, label, parent, child, wrapped, unwrapped in rows],
        colLabels=['Joint', 'Body axes', 'CAD angle', 'Sweep angle'],
        cellLoc='left', colLoc='left', bbox=[0, 0.595, 1, 0.305],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor('white')
        cell.set_facecolor('#edf1f7' if row % 2 else '#f7f9fc')
        if row == 0:
            cell.set_facecolor('#dce5ef')
            cell.set_text_props(fontweight='bold')

    info.text(0, 0.565, 'Neutral orientation in each joint’s local frame',
              fontsize=11, fontweight='bold')
    for i, (_, label, _, _, wrapped, _) in enumerate(rows):
        row, col = divmod(i, 4)
        icon = info.inset_axes([0.25 * col, 0.30 - 0.22 * row, 0.24, 0.20])
        angle = np.deg2rad(wrapped)
        icon.annotate('', xy=(1, 0), xytext=(0, 0),
                      arrowprops={'arrowstyle': '-|>', 'color': parent_color, 'lw': 1.4})
        icon.annotate('', xy=(0.77 * np.cos(angle), 0.77 * np.sin(angle)), xytext=(0, 0),
                      arrowprops={'arrowstyle': '-|>', 'color': child_color, 'lw': 1.4})
        arc = np.linspace(0, angle, 45)
        icon.plot(0.36 * np.cos(arc), 0.36 * np.sin(arc), color='#5b287e', lw=1)
        icon.set(xlim=(-1.13, 1.13), ylim=(-1.05, 1.05), aspect='equal',
                 title=f'{label}   {wrapped:+.1f}°')
        icon.title.set_fontsize(8.5)
        icon.set_axis_off()
    info.text(0.75, 0.20, 'Parent axis (0°)', fontsize=8.5, color=parent_color)
    info.text(0.75, 0.17, 'Child axis', fontsize=8.5, color=child_color)
    info.text(0.75, 0.14, 'Neutral angle', fontsize=8.5, color='#5b287e')
    info.text(0, 0.012, 'Each diagram shows a separate neutral pose. The left view is the assembly.\n'
              'CAD angles are wrapped to ±180°. Sweep angles follow the unwrapped model.',
              fontsize=8.5, color='#526174', linespacing=1.35)


def draw_cfp_placement(ax, info, xy, pos, title, mounts):
    """Show both body-fixed clamp frames at each nominal virtual pivot."""
    ax.clear()
    info.clear()
    info.set_axis_off()

    def draw_mount(axis, joint, center, *, mark_ground=False):
        _, offsets, _ = mounts[joint]
        parent, child = FLEXURE_PAIRS[joint]
        color = LINK_COLORS[parent]
        clamps = [rotate(points, pos[f'theta{body}'] + offset) + center
                  for points, body, offset in zip(NEUTRAL_CLAMPS, (parent, child), offsets)]
        for i in (0, 1):
            blade = np.array([clamps[0][i], clamps[1][1-i]])
            axis.plot(*blade.T, color=color, lw=1, ls='--', zorder=4)
        axis.plot(*clamps[0].T, color=color, lw=2, zorder=5,
                  marker='s' if mark_ground else None, mfc='white', ms=4)
        axis.plot(*clamps[1].T, color=LINK_COLORS[child], lw=2, zorder=5)
        axis.plot(*center, 'o', color=color, mfc='white', ms=3, zorder=6)

    for (a, b), color in zip(LINK_ENDPOINTS, LINK_COLORS):
        ax.plot(*np.stack((xy[a], xy[b])).T, color=color, lw=2.5, alpha=.57, zorder=2)

    offsets = {'A': (-10, 8), 'B': (23, -16), 'C': (-12, -10), 'D': (18, -11),
               'P': (-16, 14), 'Q': (13, 14), 'T1': (-19, 15), 'T2': (17, 14)}
    for pivot, _, spring in PIVOT_INTERFACES:
        if spring == 'B_BD':
            # Both mounts have the same centre and footprint in the planar view.
            continue
        label = pivot
        if spring is None:
            ax.plot(*xy[pivot], 'o', color='#34495e', markerfacecolor='white', markersize=6, zorder=5)
        elif spring in mounts:
            draw_mount(ax, spring, xy[pivot], mark_ground=pivot == 'B')
            if pivot == 'B':
                label = 'B · 2 mounts'
        ax.annotate(label, xy[pivot], xytext=offsets[pivot], textcoords='offset points',
                    ha='center', va='center', fontsize=9, fontweight='bold', zorder=7,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1})
    ax.set(xlabel='x [mm]', ylabel='y [mm]', aspect='equal')
    ax.set_title(title, loc='left', fontsize=11)
    ax.margins(x=.23, y=.16)
    ax.grid(color='#e4eaf0', zorder=0)
    ax.spines[['top', 'right']].set_visible(False)

    info.text(0, .965, 'One reference CFP design', fontsize=13, fontweight='bold')
    info.text(0, .910, f'{CFP_LENGTH_MM:g} mm blade length  ·  {CFP_CROSSING_DEG:g}° crossing',
              fontsize=10)
    info.text(0, .865, f'{CFP_WIDTH_MM:.2f} × {CFP_HEIGHT_MM:.2f} mm neutral clamp footprint',
              fontsize=10)
    info.text(0, .820, f'{CFP_THICKNESS_MM:g} mm blade thickness  ·  {CFP_BLADE_WIDTH_MM:g} mm depth per blade',
              fontsize=10)
    info.text(0, .755, 'Solid clamp lines use their body colours. Dashed diagonals',
              fontsize=9, color='#526174')
    info.text(0, .722, 'connect blade ends schematically in the parent colour.',
              fontsize=9, color='#526174')
    if not mounts:
        info.text(0, .60, 'A closed full sweep is required to set neutral clamp offsets.',
                  fontsize=9, color='#526174')
        return
    rows = []
    for joint, (_, offsets, neutral) in mounts.items():
        parent, child = FLEXURE_PAIRS[joint]
        twist_deg = np.rad2deg(pos[f'theta{child}'] - pos[f'theta{parent}']
                              + offsets[1] - offsets[0])
        rows.append((joint.replace('_', '–'), f'{np.rad2deg(neutral):+.2f}°',
                     f'{(twist_deg + 180) % 360 - 180:+.2f}°'))
    table = info.table(cellText=rows, colLabels=['Joint', 'Neutral body angle', 'Assembly twist'],
                       cellLoc='left', colLoc='left', bbox=[0, .285, 1, .375])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor('white')
        cell.set_facecolor('#dce5ef' if row == 0 else '#edf1f7')
    for i, joint in enumerate(('B_BQ', 'B_BD')):
        icon = info.inset_axes([.08 + i * .45, .027, .34, .135])
        draw_mount(icon, joint, np.zeros(2), mark_ground=True)
        icon.set(xlim=(-18, 18), ylim=(-18, 18), aspect='equal')
        icon.set_title(f'B → L{FLEXURE_PAIRS[joint][1]}', fontsize=9)
        icon.set_axis_off()
    info.text(0, .225, 'B shares one pivot in plan. Squares show the separate ground mount.',
              fontsize=8.5, color='#526174')


def schematic_geometry(cfg, pos, reference_poses=None):
    """Share the closed sweep and neutral mounting frames across all exports."""
    xy = {name: pos[name] * 1000 for name, _, _ in PIVOT_INTERFACES}
    poses = reference_poses
    if poses is None:
        angles = pos['theta6'] - np.linspace(0, 2 * np.pi, 721)
        try:
            poses = list(sweep_driven_positions(cfg, angles, theta1_guess_rad=pos['theta1']))
        except ValueError:
            poses = []
    if poses and not np.allclose(poses[0]['P'], poses[-1]['P'], rtol=0, atol=1e-8):
        poses = []
    ranges = {}
    for name, (a, b) in FLEXURE_PAIRS.items() if poses else ():
        relative_deg = np.rad2deg(relative_flexure_angles(poses, a, b))
        ranges[name] = (float(relative_deg.min()), float(relative_deg.max()))
    mounts = {}
    for pivot, _, joint in PIVOT_INTERFACES:
        if joint in ranges:
            low, high = ranges[joint]
            neutral = np.deg2rad(((low + high) / 2 + 180) % 360 - 180)
            mounts[joint] = (pivot, clamp_offsets(joint, neutral), neutral)
    return xy, poses, ranges, mounts


def save_prbm_schematic(out, cfg, pos, *, reference_poses=None, geometry=None):
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnnotationBbox, DrawingArea
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, Polygon

    if geometry is None:
        geometry = schematic_geometry(cfg, pos, reference_poses)
    xy, poses, ranges, mounts = geometry
    lengths_mm = cfg.L * 1000
    width_mm = ELEMENT_WIDTH_M * 1000
    spring_color = '#5b287e'
    revolute_color = '#38495a'
    colors = LINK_COLORS
    joint_offsets = {
        'A': (-12, 13), 'B': (-14, -13), 'P': (-6, -14), 'Q': (11, -12),
        'C': (-13, -17), 'D': (10, -14), 'T1': (-16, 13), 'T2': (10, 11),
    }
    link_offsets = [(5, -10), (-11, 2), (-12, 3), (12, 3), (3, 9), (11, 1),
                    (14, -1), (12, -1), (-12, -2), (1, -10), (-12, -5), (8, 11)]

    fig, (ax, info) = plt.subplots(
        1, 2, figsize=(12, 7.4), gridspec_kw={'width_ratios': [1.37, 1]},
        layout='constrained',
    )
    fig.get_layout_engine().set(h_pad=0.12, w_pad=0.10)
    fig.suptitle('Concept 2 | Lumped-flexure reference geometry', x=0.03, ha='left',
                 fontsize=19, fontweight='bold')
    title = f'AC angle {np.rad2deg(pos["theta6"]):+.2f}° from +x'
    if reference_poses is not None:
        title = 'Minimum-energy assembly pose\n' + title
    ax.set_title(title, loc='left', fontsize=11)
    for index, (start, end) in enumerate(LINK_ENDPOINTS):
        coordinate_pair = np.array([xy[start], xy[end]])
        ax.plot(*coordinate_pair.T, color=colors[index], linewidth=2.5, zorder=2)
        ax.annotate(f'L{index}', coordinate_pair.mean(axis=0), xytext=link_offsets[index],
                    textcoords='offset points', ha='center', va='center', fontsize=8,
                    color=colors[index], bbox={'fc': 'white', 'ec': 'none', 'pad': 0.6})

    theta = np.linspace(0, 3.6 * np.pi, 80)
    radius = np.linspace(0.6, 5.4, len(theta))
    for name in ('A', 'B'):
        ground = DrawingArea(30, 30)
        shift = -15 if name == 'A' else 0
        ground.add_artist(Polygon([(15, 29), (8 + shift, 11), (22 + shift, 11)],
                                  facecolor='white', edgecolor=revolute_color, linewidth=1.1))
        ground.add_artist(Line2D([2 + shift, 28 + shift], [11, 11],
                                 color=revolute_color, linewidth=1.1))
        for x in range(4, 27, 5):
            ground.add_artist(Line2D([x + shift, x - 4 + shift], [11, 5],
                                     color=revolute_color, linewidth=0.8))
        ax.add_artist(AnnotationBbox(ground, xy[name], box_alignment=(0.5, 1),
                                     frameon=False, pad=0, zorder=3))
    for name, point in xy.items():
        symbol = DrawingArea(15, 15)
        symbol.add_artist(Circle((7.5, 7.5), 7.2, color='white'))
        if name in REVOLUTE_PIVOTS:
            joint_color = revolute_color
            symbol.add_artist(Circle((7.5, 7.5), 4.7, fill=False, ec=joint_color, lw=1.25))
            symbol.add_artist(Circle((7.5, 7.5), 1.4, color=joint_color))
        else:
            joint_color = spring_color
            symbol.add_artist(Line2D(7.5 + radius * np.cos(theta), 7.5 + radius * np.sin(theta),
                                     color=joint_color, linewidth=1.05))
        ax.add_artist(AnnotationBbox(symbol, point, frameon=False, pad=0, zorder=4))
        ax.annotate(name, point, xytext=joint_offsets[name], textcoords='offset points',
                    ha='center', va='center', color=joint_color, fontweight='bold', zorder=5)
    wing = pos['wing_position_m'] * 1000
    ax.plot(*wing, marker='D', color=colors[11], markersize=5, zorder=5)
    ax.annotate('Wing midpoint', wing, xytext=(53, 32), textcoords='offset points',
                color=colors[11], fontsize=8, ha='left', va='center', zorder=5,
                arrowprops={'arrowstyle': '-', 'color': colors[11], 'lw': 0.65})
    overlay_neutral_angles(ax, xy, pos, ranges, lengths_mm[4])
    ax.set(xlabel='x [mm]', ylabel='y [mm]', aspect='equal')
    ax.margins(x=0.24, y=0.16)
    ax.grid(color='#e4eaf0', zorder=0)
    ax.spines[['top', 'right']].set_visible(False)

    info.set_axis_off()
    info.text(0, 0.965, 'Pivot-to-pivot dimensions', fontsize=12, fontweight='bold')
    info.text(0.99, 0.965, f'Width: {width_mm:g} mm', ha='right', fontsize=9, color='#526174')
    rows = [(f'L{i}', f'{a} – {b}', f'{length:.3f}')
            for i, ((a, b), length) in enumerate(zip(LINK_ENDPOINTS, lengths_mm))]
    table = info.table(cellText=rows, colLabels=['Link', 'Endpoints', 'Length [mm]'],
                       cellLoc='left', colLoc='left', bbox=[0, 0.433, 1, 0.505])
    joint_rows = []
    for pivot, bodies, spring in PIVOT_INTERFACES:
        if spring is None:
            joint_rows.append((pivot, 'Revolute', bodies, '—'))
        else:
            if spring in ranges:
                minimum, maximum = ranges[spring]
                twist = f'±{(maximum-minimum)/2:.1f}°'
            else:
                twist = 'n/a'
            joint_rows.append((pivot, 'Spring candidate', bodies, twist))
    interfaces = info.table(cellText=joint_rows,
                            colLabels=['Pivot', 'Joint type', 'Connected bodies', 'Twist'],
                            colWidths=[0.10, 0.29, 0.43, 0.18], cellLoc='left', colLoc='left',
                            bbox=[0, 0.089, 1, 0.292])
    for section in (table, interfaces):
        section.auto_set_font_size(False)
        section.set_fontsize(9)
        for (row, _), cell in section.get_celld().items():
            cell.set_edgecolor('white')
            cell.set_facecolor('#edf1f7' if row % 2 else '#f7f9fc')
            if row == 0:
                cell.set_facecolor('#dce5ef')
                cell.set_text_props(fontweight='bold')
    info.text(0, 0.404, 'Joint types and torsion spring travel',
              fontsize=11, fontweight='bold')
    travel_note = ('Twist is relative to the midpoint of each spring’s angular range over 360°.'
                   if poses else 'Twist is unavailable for this assembly branch over 360°.')
    info.text(0, 0.013, 'Hatched triangle: ground. Ring + dot: revolute joint.\n'
              'Purple arc: neutral angle from dotted parent to magenta dashed child.\n' + travel_note,
              fontsize=8.5, color='#526174', linespacing=1.35)
    fig.savefig(out / 'prbm_schematic.svg')
    fig.savefig(out / 'prbm_schematic.png', dpi=160)
    draw_neutral_angles(info, ranges)
    fig.suptitle('Concept 2 | Neutral joint angles', x=0.03, ha='left',
                 fontsize=19, fontweight='bold')
    fig.savefig(out / 'prbm_neutral_angles.svg')
    fig.savefig(out / 'prbm_neutral_angles.png', dpi=160)
    draw_cfp_placement(ax, info, xy, pos, title, mounts)
    fig.suptitle('Concept 2 | Crossed-flexure placement', x=0.03, ha='left',
                 fontsize=19, fontweight='bold')
    fig.savefig(out / 'prbm_cfp_placement.svg')
    fig.savefig(out / 'prbm_cfp_placement.png', dpi=160)
    plt.close(fig)
    if mounts:
        save_fabrication_layout(out, xy, pos, mounts, LINK_ENDPOINTS, LINK_COLORS, poses)
    with (out / 'prbm_link_geometry.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['link', 'pivot_a', 'pivot_b', 'xa_mm', 'ya_mm', 'xb_mm', 'yb_mm',
                         'nominal_length_mm', 'element_width_mm'])
        for i, (a, b) in enumerate(LINK_ENDPOINTS):
            dimensions = (*xy[a], *xy[b], lengths_mm[i], width_mm)
            writer.writerow([f'L{i}', a, b, *(f'{value:.9f}' for value in dimensions)])
    return bool(mounts)
