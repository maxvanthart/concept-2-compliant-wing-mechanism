"""Differentiate loops 1 and 2: solve 2x2 Jacobian systems.

At a rigid-body compatible configuration this also gives apex velocities.
Loop 3 is a consistency check; without closure, omega11 is NOT a
physically valid full-mechanism angular speed.
"""
import numpy as np
from .geometry import rotate90, unit, cross2, loop_residuals


def solve_velocities(cfg, pos, *, omega1_rad_s=None, omega6_rad_s=None):
    L = cfg.L
    w1 = cfg.omega1_rad_s if omega1_rad_s is None else omega1_rad_s
    w6 = w1 if omega6_rad_s is None else omega6_rad_s
    # omega1*R1_perp + omega4*R4_perp - omega5*R5_perp = 0
    j1 = np.column_stack((L[4]*rotate90(unit(pos['theta4'])),
                          -L[5]*rotate90(unit(pos['theta5']))))
    # omega6*R6_perp + omega9*R9_perp - omega10*R10_perp = 0
    j2 = np.column_stack((L[9]*rotate90(unit(pos['theta9'])),
                          -L[10]*rotate90(unit(pos['theta10']))))
    if any(np.linalg.matrix_rank(jacobian) < 2 for jacobian in (j1, j2)):
        raise ValueError('Four-bar velocity is indeterminate at a toggle.')
    w4, w5 = np.linalg.solve(j1, -w1*L[1]*rotate90(unit(pos['theta1'])))
    w9, w10 = np.linalg.solve(j2, -w6*L[6]*rotate90(unit(pos['theta6'])))
    # Apex velocities from rigid triangle bodies.
    vp = w1*rotate90(pos['P']-pos['A'])
    vc = w6*rotate90(pos['C']-pos['A'])
    vt1 = vp+w4*rotate90(pos['T1']-pos['P'])
    vt2 = vc+w9*rotate90(pos['T2']-pos['C'])
    d = pos['T2']-pos['T1']
    vrel = vt2-vt1
    denom = float(d@d)
    if denom < 1e-16:
        omega11_geom = np.nan
        radial_rate = np.nan
    else:
        omega11_geom = cross2(d, vrel)/denom
        radial_rate = float((d@vrel)/np.sqrt(denom))
    speeds = (0., w1, w4, w4, w4, w5, w6, w9, w9, w9, w10, omega11_geom)
    derivatives = [w * length * rotate90(unit(pos[f'theta{i}']))
                   for i, (w, length) in enumerate(zip(speeds, L))]
    r1, r2, r3 = loop_residuals(derivatives)
    return dict(omega0=0., omega1=w1, omega6=w6, omega4=float(w4),omega5=float(w5),
                omega2=float(w4), omega3=float(w4), omega7=float(w9),
                omega8=float(w9), omega9=float(w9), omega10=float(w10),
                omega11_geometric=float(omega11_geom),
                wing_velocity_m_s=(vt1 + vt2) / 2,
                wing_omega_rad_s=float(omega11_geom),
                required_L11_rate_m_s=radial_rate,
                loop1_velocity_residual=r1, loop2_velocity_residual=r2,
                loop3_velocity_residual=r3,
                velocity_T1_m_s=vt1, velocity_T2_m_s=vt2,
                jacobian_loop1=j1, jacobian_loop2=j2)


def solve_driven_velocities(cfg, pos, omega6_rad_s):
    """Differentiate the L11 constraint with AC as the sole driven crank."""
    if not abs(pos['closure_gap_m']) <= max(1e-10, cfg.scale_m * 1e-9):
        raise ValueError('AC-driven velocity requires a closed position')
    ap = solve_velocities(cfg, pos, omega1_rad_s=1.0, omega6_rad_s=0.0)
    ac = solve_velocities(cfg, pos, omega1_rad_s=0.0, omega6_rad_s=1.0)
    radial_per_ap = ap['required_L11_rate_m_s']
    if not np.isfinite(radial_per_ap) or abs(radial_per_ap) <= cfg.scale_m * 1e-10:
        raise ValueError('Loop 3 is singular for an AC-driven input at this pose')
    ratio = -ac['required_L11_rate_m_s'] / radial_per_ap
    return solve_velocities(cfg, pos, omega1_rad_s=ratio * omega6_rad_s, omega6_rad_s=omega6_rad_s)
