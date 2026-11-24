"""Draw the nominal linkage, revolute joints and candidate torsion springs."""
import csv
import json
from pathlib import Path

import numpy as np

from .config import ELEMENT_WIDTH_M
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


def save_flexure_table(out, ranges, width_mm):
    joints = []
    for pivot, _, spring in PIVOT_INTERFACES:
        if spring in ranges:
            minimum, maximum = ranges[spring]
            joints.append({
                'id': spring,
                'label': f'{pivot} · Ground–{spring[2:]}' if pivot == 'B' else pivot,
                'minDeg': minimum,
                'maxDeg': maximum,
            })
    template = Path(__file__).with_name('flexure_table.html').read_text(encoding='utf-8')
    payload = json.dumps({'widthMm': width_mm, 'joints': joints})
    (out / 'flexure_table.html').write_text(template.replace('__GEOMETRY_JSON__', payload),
                                           encoding='utf-8')


def save_prbm_schematic(out, cfg, pos):
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnnotationBbox, DrawingArea
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, Polygon

    xy = {name: pos[name] * 1000 for name, _, _ in PIVOT_INTERFACES}
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
    lengths_mm = cfg.L * 1000
    width_mm = ELEMENT_WIDTH_M * 1000
    spring_color = '#5b287e'
    revolute_color = '#38495a'
    colors = ['#526174'] + ['#23649a'] * 5 + ['#177965'] * 5 + ['#a8402c']
    joint_offsets = {
        'A': (-12, 13), 'B': (17, -13), 'P': (-13, 13), 'Q': (11, 11),
        'C': (-13, -17), 'D': (10, -14), 'T1': (-16, 13), 'T2': (10, 11),
    }
    link_offsets = [(5, -10), (-11, 2), (-12, 3), (12, 3), (3, 9), (11, 1),
                    (-27, 1), (12, -1), (-12, -2), (1, -10), (-12, -5), (8, 11)]

    fig, (ax, info) = plt.subplots(
        1, 2, figsize=(12, 7.4), gridspec_kw={'width_ratios': [1.37, 1]},
        layout='constrained',
    )
    fig.get_layout_engine().set(h_pad=0.12, w_pad=0.10)
    fig.suptitle('Concept 2 | Lumped-flexure reference geometry', x=0.03, ha='left',
                 fontsize=19, fontweight='bold')
    ax.set_title(f'AC angle {np.rad2deg(pos["theta6"]):+.2f}° from +x', loc='left', fontsize=11)
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
        ground.add_artist(Polygon([(15, 29), (8, 11), (22, 11)],
                                  facecolor='white', edgecolor=revolute_color, linewidth=1.1))
        ground.add_artist(Line2D([2, 28], [11, 11], color=revolute_color, linewidth=1.1))
        for x in range(4, 27, 5):
            ground.add_artist(Line2D([x, x-4], [11, 5], color=revolute_color, linewidth=0.8))
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
    ax.annotate('Wing midpoint', wing, xytext=(37, 7), textcoords='offset points',
                color=colors[11], fontsize=8, ha='left', va='center', zorder=5,
                arrowprops={'arrowstyle': '-', 'color': colors[11], 'lw': 0.65})
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
              'Purple spiral: torsion spring candidate.\n' + travel_note,
              fontsize=8.5, color='#526174', linespacing=1.35)
    fig.savefig(out / 'prbm_schematic.svg')
    fig.savefig(out / 'prbm_schematic.png', dpi=160)
    plt.close(fig)
    save_flexure_table(out, ranges, width_mm)

    with (out / 'prbm_link_geometry.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['link', 'pivot_a', 'pivot_b', 'xa_mm', 'ya_mm', 'xb_mm', 'yb_mm',
                         'nominal_length_mm', 'element_width_mm'])
        for i, (a, b) in enumerate(LINK_ENDPOINTS):
            dimensions = (*xy[a], *xy[b], lengths_mm[i], width_mm)
            writer.writerow([f'L{i}', a, b, *(f'{value:.9f}' for value in dimensions)])
