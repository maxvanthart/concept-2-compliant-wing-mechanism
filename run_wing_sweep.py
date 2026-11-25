"""Plot the Concept 2 wing path and quasi-static input torque contributions."""
import argparse
import csv
from pathlib import Path

import numpy as np

from concept2.config import Concept2Config, MOTIONGEN_THETA1_RAD, MOTIONGEN_THETA6_RAD
from concept2.geometry import sweep_driven_positions
from concept2.prbm import FLEXURE_PAIRS, relative_flexure_angles
from concept2.velocity import solve_driven_velocities
from dimensioning.io import load_config


BODY_ANGLES = {f'L{i}': i for pair in FLEXURE_PAIRS.values() for i in pair}
WING_CHORD_M = 0.04
WING_SPAN_M = 1.05
WING_AREA_M2 = WING_CHORD_M * WING_SPAN_M


def sweep(direction):
    cfg = Concept2Config()
    q = np.linspace(0, 2 * np.pi, 721)
    positions = list(sweep_driven_positions(cfg, MOTIONGEN_THETA6_RAD + direction * q,
                                           theta1_guess_rad=MOTIONGEN_THETA1_RAD))
    velocities = [solve_driven_velocities(cfg, pos, direction) for pos in positions]
    return q, positions, velocities


def hinge_motion(a, b, positions, velocities):
    angle = relative_flexure_angles(positions, a, b)

    def rate(link, state):
        return state['omega11_geometric' if link == 11 else f'omega{link}']

    sensitivity = np.array([rate(b, v) - rate(a, v) for v in velocities])
    return angle, sensitivity


def wing_incidence(phi, motion):
    """Acute angle to the flow and projected area fraction for a thin plate."""
    motion_angle = np.arctan2(motion[:, 1], motion[:, 0])
    moving = np.any(motion != 0, axis=1)
    projected_fraction = np.where(moving, np.abs(np.sin(phi - motion_angle)), 0.0)
    angle_deg = np.where(moving, np.rad2deg(np.arcsin(projected_fraction)), np.nan)
    return angle_deg, projected_fraction


def calculate(hinges=None, rpm=None, drag_factor=None, *, direction=-1, common_k=None):
    q, positions, velocities = sweep(direction)
    xy = np.array([pos['wing_position_m'] for pos in positions])
    phi = np.unwrap([pos['wing_orientation_rad'] for pos in positions])
    dr_dq = np.array([vel['wing_velocity_m_s'] for vel in velocities])
    dphi_dq = np.array([vel['wing_omega_rad_s'] for vel in velocities])
    speed_per_rad_s = np.linalg.norm(dr_dq, axis=1)
    incidence_deg, projected_fraction = wing_incidence(phi, dr_dq)
    data = {
        'q_rad': q, 'ac_travel_deg': np.rad2deg(q),
        'ac_angle_rad': MOTIONGEN_THETA6_RAD + direction * q,
        'wing_x_m': xy[:, 0], 'wing_y_m': xy[:, 1],
        'wing_orientation_rad': phi, 'wing_dx_dq_m_rad': dr_dq[:, 0],
        'wing_dy_dq_m_rad': dr_dq[:, 1], 'wing_dphi_dq_rad_rad': dphi_dq,
        'wing_speed_per_input_rad_s': speed_per_rad_s,
        'wing_incidence_deg': incidence_deg, 'wing_projected_area_fraction': projected_fraction,
        'wing_projected_area_m2': WING_AREA_M2 * projected_fraction,
    }
    flex = np.zeros_like(q)
    if hinges is None:
        for name, (a, b) in FLEXURE_PAIRS.items():
            angle, sensitivity = hinge_motion(a, b, positions, velocities)
            rest = 0.5 * (angle.min() + angle.max())
            flex += (angle - rest) * sensitivity
            data[f'hinge_{name}_rad'] = angle
        flex_key = 'flexure_torque_per_common_k'
        if common_k is not None:
            flex *= common_k
            flex_key = 'flexure_torque_nm'
    else:
        for hinge in hinges:
            a, b = (BODY_ANGLES.get(link) for link in (hinge.link_a, hinge.link_b))
            if not any({a, b} == set(pair) for pair in FLEXURE_PAIRS.values()):
                raise ValueError(f'{hinge.name}: use a candidate flexure body pair (A and C are revolute)')
            angle, sensitivity = hinge_motion(a, b, positions, velocities)
            flex += hinge.k_nm_rad * (angle - hinge.rest_angle_rad) * sensitivity
            data[f'hinge_{hinge.name}_rad'] = angle
        flex_key = 'flexure_torque_nm'
    data[flex_key] = flex

    if rpm is not None:
        speed = (rpm * 2 * np.pi / 60) * dr_dq
        data['wing_speed_m_s'] = np.linalg.norm(speed, axis=1)
    if drag_factor is None:
        drag_key = 'wing_drag_torque_per_drag_factor_and_speed_squared'
        data[drag_key] = projected_fraction * speed_per_rad_s**3
    else:
        # Still air. The drag factor is based on the broadside planform area.
        force = -(drag_factor * projected_fraction * data['wing_speed_m_s'])[:, None] * speed
        data.update(wing_fx_n=force[:, 0], wing_fy_n=force[:, 1])
        drag_key = 'wing_drag_torque_nm'
        data[drag_key] = -np.einsum('ij,ij->i', force, dr_dq)
    if flex_key == 'flexure_torque_nm' and drag_factor is not None:
        data['total_quasistatic_torque_nm'] = data[flex_key] + data[drag_key]
    return data, flex_key, drag_key


def save_sweep(out, data, flex_key, drag_key):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    out.mkdir(parents=True, exist_ok=True)
    with (out / 'wing_sweep.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(data)
        writer.writerows(zip(*data.values()))

    q = data['ac_travel_deg']
    sense = 'Clockwise' if data['ac_angle_rad'][-1] < data['ac_angle_rad'][0] else 'Counterclockwise'
    travel_label = f'{sense} AC travel from start [deg]'
    x = data['wing_x_m'] * 1000
    y = data['wing_y_m'] * 1000
    blue, green, violet, orange = '#23649a', '#177965', '#743f8c', '#c95715'
    fig, layout = plt.subplot_mosaic(
        [['path', 'position'], ['path', 'rotation']],
        figsize=(12, 6), width_ratios=[1.2, 1], layout='constrained',
    )
    fig.suptitle('Concept 2 | Wing midpoint and rotation', fontsize=15, fontweight='bold')
    path, position, rotation = (layout[name] for name in ('path', 'position', 'rotation'))
    path.plot(x, y, color=blue)
    start = path.scatter(x[0], y[0], color=blue, label='Start / end', zorder=4)
    distance = np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
    stations = np.linspace(0, distance[-1], 18, endpoint=False) + distance[-1] / 36
    samples = np.searchsorted(distance, stations)
    phi = data['wing_orientation_rad']
    arrow_length_mm = 8
    path.quiver(x[samples], y[samples],
                arrow_length_mm * np.cos(phi[samples]), arrow_length_mm * np.sin(phi[samples]),
                angles='xy', scale_units='xy', scale=1, units='dots', width=1.7,
                color=violet, zorder=3)
    first, last = np.searchsorted(q, [188, 203])
    path.annotate('', xy=(x[last], y[last]), xytext=(x[first], y[first]),
                  arrowprops={'arrowstyle': '-|>', 'color': orange, 'lw': 2.1,
                              'mutation_scale': 19, 'shrinkA': 0, 'shrinkB': 0},
                  zorder=5)
    path.set(xlabel='x [mm]', ylabel='y [mm]', title='Midpoint trajectory', aspect='equal')
    path.margins(x=0.09, y=0.25)
    path.legend(handles=[
        start, Line2D([], [], marker='>', linestyle='none', color=orange, label='Motion direction'),
        Line2D([], [], marker='>', linestyle='none', color=violet, label='L11 direction (T1 to T2)'),
    ], loc='upper center', bbox_to_anchor=(0.5, -0.31), ncols=3, frameon=False)
    position.plot(q, x, color=blue, label='x')
    position.plot(q, y, color=green, label='y')
    position.set(ylabel='Position [mm]', title='Midpoint coordinates')
    position.legend()
    rotation.plot(q, np.rad2deg(data['wing_orientation_rad']), color=violet)
    rotation.set(ylabel='T1 to T2 angle [deg]', title='Wing orientation')
    for ax in (position, rotation):
        ax.set(xlabel=travel_label, xlim=(0, 360), xticks=range(0, 361, 90))
    for ax in (path, position, rotation):
        ax.grid(alpha=0.24)
    fig.savefig(out / 'wing_kinematics.png', dpi=165)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(9, 8.3), sharex=True, layout='constrained')
    fig.suptitle('Concept 2 | Input torque contributions', fontsize=15, fontweight='bold')
    flex_units = ('Input torque [N·m]' if flex_key.endswith('_nm')
                  else 'Input torque / common stiffness [rad]')
    drag_units = ('Input torque [N·m]' if drag_key.endswith('_nm')
                  else 'Input torque / (c · input speed²) [m³]')
    flex_title = ('Flexures' if flex_key.endswith('_nm')
                  else 'Flexures | normalized, common stiffness and centered rest angles')
    drag_title = ('Wing drag | projected area, still air' if drag_key.endswith('_nm')
                  else 'Wing drag | normalized, projected area, still air')
    for ax, key, color, units, title in zip(
        axes[:2], (flex_key, drag_key), (blue, green), (flex_units, drag_units), (flex_title, drag_title)
    ):
        ax.plot(q, data[key], color=color)
        ax.axhline(0, color='#657284', linewidth=0.7)
        ax.set(ylabel=units, title=title)
        ax.grid(alpha=0.24)
    axes[-1].plot(q, data['wing_incidence_deg'], color=violet)
    axes[-1].set(ylabel='Acute angle [deg]', title='L11 angle to the relative airflow', ylim=(0, 90))
    axes[-1].grid(alpha=0.24)
    axes[-1].set(xlabel=travel_label, xlim=(0, 360), xticks=range(0, 361, 45))
    fig.savefig(out / 'input_torque.png', dpi=165)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hinges', type=Path, help='Physical flexure JSON with stiffness and rest angles')
    cfp = parser.add_argument_group('Identical crossed-flexure pivots (dimensions per blade)')
    cfp.add_argument('--common-k', type=float, help='Shared fitted stiffness in N m/rad')
    cfp.add_argument('--cfp-length-mm', type=float, help='Free diagonal blade length')
    cfp.add_argument('--cfp-width-mm', type=float, help='Width of one blade')
    cfp.add_argument('--cfp-thickness-mm', type=float, help='Thickness of one blade')
    cfp.add_argument('--cfp-e-gpa', type=float, help="Material Young's modulus")
    parser.add_argument('--rpm', type=float, help='Non-negative AC speed in revolutions per minute')
    parser.add_argument('--counterclockwise', action='store_true', help='Reverse the default clockwise input rotation')
    parser.add_argument('--drag-factor', type=float, help='Broadside c = rho * C_D * A / 2 [kg/m]')
    parser.add_argument('--air-density', type=float, help='Air density [kg/m³] for the rectangular wing')
    parser.add_argument('--drag-coefficient', type=float, help='Broadside C_D of the rectangular wing')
    parser.add_argument('--out', type=Path, default=Path('results/C2_WING_SWEEP'))
    args = parser.parse_args()
    if args.rpm is not None and (not np.isfinite(args.rpm) or args.rpm < 0):
        parser.error('--rpm must be finite and non-negative; use --counterclockwise to reverse')
    common_k = args.common_k
    if common_k is not None and (args.hinges is not None or not np.isfinite(common_k) or common_k <= 0):
        parser.error("--common-k must be positive and finite and cannot be used with --hinges.")
    blade = (args.cfp_length_mm, args.cfp_width_mm, args.cfp_thickness_mm, args.cfp_e_gpa)
    if any(value is not None for value in blade):
        if args.hinges is not None or common_k is not None:
            parser.error('Use one of the shared CFP dimensions, --common-k, or --hinges.')
        if not all(value is not None and np.isfinite(value) and value > 0 for value in blade):
            parser.error('Supply positive finite length, width, thickness and modulus for the CFP.')
        length, width, thickness = (value / 1000 for value in blade[:3])
        modulus = blade[3] * 1e9
        common_k = 2 * modulus * width * thickness**3 / (12 * length)
    if common_k is not None:
        print(f'Shared CFP stiffness: {common_k:.6g} N m/rad per pivot '
              f'({len(FLEXURE_PAIRS)} pivots)')
    drag_factor = args.drag_factor
    if args.air_density is not None or args.drag_coefficient is not None:
        if drag_factor is not None:
            parser.error('Use --drag-factor or supply --air-density and --drag-coefficient')
        if args.air_density is None or args.drag_coefficient is None:
            parser.error('--air-density and --drag-coefficient must be supplied together')
        if not np.isfinite(args.air_density) or args.air_density <= 0:
            parser.error('--air-density must be finite and positive')
        if not np.isfinite(args.drag_coefficient) or args.drag_coefficient < 0:
            parser.error('--drag-coefficient must be finite and non-negative')
        drag_factor = 0.5 * args.air_density * args.drag_coefficient * WING_AREA_M2
    if drag_factor is not None:
        if args.rpm is None:
            parser.error('Computing dimensional drag requires --rpm')
        if not np.isfinite(drag_factor) or drag_factor < 0:
            parser.error('--drag-factor must be finite and non-negative')
    try:
        hinges = load_config(args.hinges) if args.hinges else None
        direction = 1 if args.counterclockwise else -1
        data, flex_key, drag_key = calculate(hinges, args.rpm, drag_factor, direction=direction, common_k=common_k)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.error(str(error))
    save_sweep(args.out, data, flex_key, drag_key)
    print('Results saved to', args.out.resolve())


if __name__ == '__main__':
    main()
