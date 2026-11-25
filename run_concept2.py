"""Run Concept 2 checks and, when possible, individual-loop kinematics."""
import argparse
from dataclasses import replace
from pathlib import Path
import numpy as np
from concept2.assembly import minimum_energy_pose
from concept2.config import (Concept2Config, CAD_ANGLES_RAD, LOOP_REFERENCE_ANGLES_RAD,
                             MOTIONGEN_THETA1_RAD, MOTIONGEN_THETA6_RAD)
from concept2.geometry import (feasibility_checks, loop_residuals, solve_driven_positions,
                               solve_positions, unit)
from concept2.schematic import LINK_COLORS, LINK_ENDPOINTS, save_prbm_schematic, schematic_geometry
from concept2.velocity import solve_driven_velocities, solve_velocities


def run(cfg, *, make_plot=True, require_closed=False, theta6_rad=None, omega6_rad_s=-1.0,
        use_minimum_energy=False, export_step=False):
    check = feasibility_checks(cfg)
    print('Corrected lengths [cm], L0..L11:', np.round(cfg.L*100, 4))
    print('Ratios relative to L0:', np.round(cfg.L/cfg.scale_m, 4))
    for label, angles in [('Raw CAD', CAD_ANGLES_RAD),
                          ('Directed CAD', LOOP_REFERENCE_ANGLES_RAD)]:
        residuals = loop_residuals([length * unit(theta) for length, theta in zip(cfg.L, angles)])
        print(f'{label} table closure norms [mm]:', np.round(np.linalg.norm(residuals, axis=1)*1000, 5))
    if check['loop2_always_impossible']:
        print('GEOMETRY ERROR: loop 2 is impossible for ANY angles.')
        print(f"  L6 + L9 = {check['loop2_moving_max_m']:.6g} m")
        print(f"  |L0 - L10| = {check['loop2_ground_min_m']:.6g} m")
        print('  Confirm length ratios/labels against the original CAD model.')
        return False
    reference_poses = None
    ac_driven = use_minimum_energy or theta6_rad is not None
    try:
        if not ac_driven:
            print('Fixed-phase diagnostic: AP and AC forced to turn together')
            pos = solve_positions(cfg)
            vel = solve_velocities(cfg, pos)
            per_rad = solve_velocities(replace(cfg, omega1_rad_s=1.0), pos)
        else:
            print('AC-driven model: AP follows from the fixed length of L11')
            if use_minimum_energy:
                pos, reference_poses = minimum_energy_pose(cfg)
                travel_deg = np.rad2deg(MOTIONGEN_THETA6_RAD - pos['theta6'])
                print(f'Minimum-energy assembly pose: {travel_deg:.4f} deg clockwise from MotionGen start')
            else:
                # Use the shortest path from the saved orientation.
                travel = (theta6_rad - MOTIONGEN_THETA6_RAD + np.pi) % (2*np.pi) - np.pi
                segments = max(1, int(np.ceil(abs(travel) / np.deg2rad(5))))
                guess = MOTIONGEN_THETA1_RAD
                for offset in np.linspace(0, travel, segments + 1):
                    step_angle = theta6_rad + offset - travel
                    pos = solve_driven_positions(cfg, step_angle, theta1_guess_rad=guess)
                    guess = pos['theta1']
            vel = solve_driven_velocities(cfg, pos, omega6_rad_s)
            per_rad = solve_driven_velocities(cfg, pos, 1.0)
    except (ValueError, np.linalg.LinAlgError) as exc:
        print('Cannot solve at this input angle:', exc)
        return False
    for loop in (1,2,3):
        print(f'Loop {loop} position closure norm [m]: '
              f"{np.linalg.norm(pos[f'loop{loop}_residual']):.7g}")
        print(f'Loop {loop} velocity closure norm [m/s]: '
              f"{np.linalg.norm(vel[f'loop{loop}_velocity_residual']):.7g}")
    print(f"Required L11 [m]: {pos['required_L11']:.7g}")
    print(f"Specified L11 [m]: {cfg.L[11]:.7g}")
    print(f"Loop 3 gap [m]: {pos['closure_gap_m']:+.7g}")
    print(f"Crank angles AP, AC from +x [deg]: {np.rad2deg(pos['theta1']):+.7g}, "
          f"{np.rad2deg(pos['theta6']):+.7g}")
    wing_x, wing_y = pos['wing_position_m']
    phi = pos['wing_orientation_rad']
    print(f'Geometric wing midpoint (x, y) [m]: ({wing_x:.7g}, {wing_y:.7g})')
    print(f'Wing angle T1 to T2 from +x [rad, deg]: {phi:+.7g}, {np.rad2deg(phi):+.7g}')
    vx, vy = vel['wing_velocity_m_s']
    print(f'Geometric wing midpoint velocity (vx, vy) [m/s]: ({vx:+.7g}, {vy:+.7g})')
    print(f"Geometric wing angular velocity [rad/s]: {vel['wing_omega_rad_s']:+.7g}")
    print('Angular speeds of each four-bar (rad/s):')
    for link in (1,4,5,6,9,10):
        print(f'  omega{link}: {vel[f"omega{link}"]:+.6g}')
    print(f"Required link-11 length rate [m/s]: {vel['required_L11_rate_m_s']:+.6g}")
    length_gradient = per_rad['required_L11_rate_m_s']
    print(f'd(required L11)/dq [m/rad]: {length_gradient:+.6g}')
    tol = max(1e-8, cfg.scale_m*1e-6)
    position_closed = all(np.linalg.norm(pos[f'loop{i}_residual']) <= tol for i in (1, 2, 3))
    input_compatible = all(np.linalg.norm(per_rad[f'loop{i}_velocity_residual']) <= tol
                           for i in (1, 2, 3))
    if not position_closed:
        print('WARNING: The pose exceeds the numerical closure tolerance.')
    if not input_compatible:
        print('WARNING: The prescribed input motion violates rigid-link velocity closure.')
    out = Path(__file__).resolve().parent / 'results' / 'C2'
    geometry = None
    if ac_driven and position_closed and input_compatible and (make_plot or export_step):
        geometry = schematic_geometry(cfg, pos, reference_poses)
    if export_step:
        if geometry is None or not geometry[3]:
            print('STEP export requires a closed AC-driven revolution to determine neutral angles.')
            return False
        from concept2.step_export import save_step_exports
        xy, poses, _, mounts = geometry
        try:
            save_step_exports(out / 'step', xy, pos, mounts, LINK_ENDPOINTS, poses, LINK_COLORS)
        except ImportError as error:
            print('STEP export needs cadquery-ocp-novtk. Install it with '
                  'python -m pip install cadquery-ocp-novtk.')
            print(error)
            return False
        except OSError as error:
            print(error)
            return False
        print('Saved coloured AP214 STEP reference surfaces under results/C2/step/')
    if make_plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator

        out.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 6))
        for nodes, color, label in [(['A','P','Q','B'], 'tab:blue','Loop 1 and triangle'),
                                    (['A','C','D','B'], 'tab:orange','Loop 2 and triangle'),
                                    (['P','T1','Q'],'tab:blue',None),
                                    (['C','T2','D'],'tab:orange',None),
                                    (['A','B'],'0.35','Ground link L0')]:
            xy = np.array([pos[k] for k in nodes])
            ax.plot(xy[:,0], xy[:,1], '-o', color=color, label=label)
        ax.plot([pos['T1'][0],pos['T2'][0]], [pos['T1'][1],pos['T2'][1]],
                'r--', label='Required span T1–T2')
        ax.plot(wing_x, wing_y, 'D', color='tab:green', label='Wing midpoint')
        for key in ('A','B','P','Q','C','D','T1','T2'):
            x,y = pos[key]; ax.annotate(key,(x,y),xytext=(5,5),textcoords='offset points')
        ax.set_aspect('equal', adjustable='box'); ax.grid(True)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.set_xlabel('x [m]'); ax.set_ylabel('y [m]')
        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
        if use_minimum_energy:
            ax.set_title(f'Minimum-energy assembly pose | AC {np.rad2deg(pos["theta6"]):+.2f}°')
        else:
            ax.set_title(f"Concept 2 geometry check — L11 mismatch = {pos['closure_gap_m']:+.4g} m")
        fig.tight_layout(); fig.savefig(out/'geometry_check.png', dpi=160); plt.close(fig)
        print('Saved results/C2/geometry_check.png')
        if geometry is not None:
            fabrication_saved = save_prbm_schematic(out, cfg, pos, reference_poses=reference_poses,
                                                     geometry=geometry)
            print('Saved results/C2/prbm_schematic.svg, prbm_schematic.png and prbm_link_geometry.csv')
            print('Saved results/C2/prbm_neutral_angles.svg and prbm_neutral_angles.png')
            print('Saved results/C2/prbm_cfp_placement.svg and prbm_cfp_placement.png')
            if fabrication_saved:
                print('Saved results/C2/prbm_fabrication_layout.svg, prbm_fabrication_layout.png '
                      'and prbm_fabrication_frames.json')
    return not require_closed or (position_closed and input_compatible)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    defaults = Concept2Config()
    p.add_argument('--L0-m', type=float, default=defaults.scale_m, help='Physical ground link length [m]')
    p.add_argument('--theta0-deg', type=float, default=np.rad2deg(CAD_ANGLES_RAD[0]))
    p.add_argument('--theta6-deg', type=float, help='Override the minimum-energy pose with an AC angle from global +x [deg]')
    p.add_argument('--omega6', type=float, help='Signed AC input speed [rad/s], default -1 (clockwise)')
    p.add_argument('--fixed-phase', action='store_true', help='Run the old equal-speed diagnostic')
    p.add_argument('--theta1-deg', type=float, help='AP input angle for the fixed-phase diagnostic')
    p.add_argument('--theta6-offset-deg', type=float, help='AC–AP offset for the fixed-phase diagnostic')
    p.add_argument('--omega1', type=float, help='AP input speed for the fixed-phase diagnostic [rad/s]')
    p.add_argument('--branch1', type=int, choices=(-1,1), default=1)
    p.add_argument('--branch2', type=int, choices=(-1,1), default=1)
    p.add_argument('--triangle1-side', type=int, choices=(-1,1), default=1)
    p.add_argument('--triangle2-side', type=int, choices=(-1,1), default=1)
    p.add_argument('--no-plot', action='store_true')
    p.add_argument('--step', action='store_true', help='Export four relaxed groups as coloured AP214 STEP surfaces')
    p.add_argument('--require-closed', action='store_true', help='Fail unless position and input velocity close at this pose')
    args = p.parse_args()
    fixed = args.fixed_phase or any(value is not None
                                    for value in (args.theta1_deg, args.theta6_offset_deg, args.omega1))
    if fixed and (args.theta6_deg is not None or args.omega6 is not None):
        p.error('AC driver options cannot be combined with fixed-phase diagnostic options')
    cfg = replace(defaults, scale_m=args.L0_m, theta0_rad=np.deg2rad(args.theta0_deg),
          theta1_rad=defaults.theta1_rad if args.theta1_deg is None else np.deg2rad(args.theta1_deg),
          theta6_minus_theta1_rad=(defaults.theta6_minus_theta1_rad if args.theta6_offset_deg is None
                                  else np.deg2rad(args.theta6_offset_deg)),
          omega1_rad_s=defaults.omega1_rad_s if args.omega1 is None else args.omega1,
          branch_loop1=args.branch1, branch_loop2=args.branch2,
          triangle1_side=args.triangle1_side, triangle2_side=args.triangle2_side)
    theta6 = None if args.theta6_deg is None else np.deg2rad(args.theta6_deg)
    success = run(cfg, make_plot=not args.no_plot, require_closed=args.require_closed,
                  theta6_rad=theta6, use_minimum_energy=not fixed and theta6 is None,
                  omega6_rad_s=-1.0 if args.omega6 is None else args.omega6,
                  export_step=args.step)
    raise SystemExit(0 if success else 1)
