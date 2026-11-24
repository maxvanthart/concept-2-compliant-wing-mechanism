"""Two analytical four-bar solutions plus the third-loop compatibility check.

The loops use EXACTLY the PDF's vector order:
  (1) R1 + R4 - R0 - R5 = 0
  (2) R6 + R9 - R0 - R10 = 0
  (3) R1 + R2 + R11 - R8 - R10 - R0 = 0
The triangular rigid bodies (2,3,4) and (7,8,9) add orientation constraints.
"""
import numpy as np
from scipy.optimize import root_scalar
from .config import Concept2Config


def unit(theta):
    return np.array([np.cos(theta), np.sin(theta)], dtype=float)


def cross2(a, b):
    return float(a[0]*b[1] - a[1]*b[0])


def rotate90(v):
    return np.array([-v[1], v[0]], dtype=float)


def angle(v):
    return float(np.arctan2(v[1], v[0]))


def circle_intersections(p, rp, q, rq, *, tol=1e-10):
    """Return both points X with |X-p|=rp and |X-q|=rq."""
    dvec = q-p
    d = np.linalg.norm(dvec)
    if d < tol:
        raise ValueError('Circle centres coincide; geometry is singular.')
    if d > rp+rq+tol:
        raise ValueError(f'Cannot close linkage: centre distance {d:.6g} > sum of links {rp+rq:.6g}.')
    if d < abs(rp-rq)-tol:
        raise ValueError(f'Cannot close linkage: centre distance {d:.6g} < difference of links {abs(rp-rq):.6g}.')
    direction = dvec / d
    a = (rp*rp - rq*rq + d*d)/(2*d)
    h2 = rp*rp-a*a
    if h2 < -tol:
        raise ValueError('No real circle intersection.')
    h = np.sqrt(max(0., h2))
    midpoint = p+a*direction
    perpendicular = rotate90(direction)
    return midpoint+h*perpendicular, midpoint-h*perpendicular


def triangle_vertex(start, end, side_from_start, side_from_end, orientation=1):
    points = circle_intersections(start, side_from_start, end, side_from_end)
    return points[0] if orientation == 1 else points[1]


def loop_residuals(vectors):
    """Evaluate all three PDF closures on directed link vectors indexed 0–11."""
    r0, r1, r2, _, r4, r5, r6, _, r8, r9, r10, r11 = vectors
    return (r1 + r4 - r0 - r5,
            r6 + r9 - r0 - r10,
            r1 + r2 + r11 - r8 - r10 - r0)


def solve_positions(cfg: Concept2Config, *, theta1_rad=None, theta6_rad=None):
    L = cfg.L
    A = np.zeros(2)
    B = L[0]*unit(cfg.theta0_rad)
    t1 = cfg.theta1_rad if theta1_rad is None else theta1_rad
    t6 = t1 + cfg.theta6_minus_theta1_rad if theta6_rad is None else theta6_rad
    P = A + L[1]*unit(t1)
    C = A + L[6]*unit(t6)

    # Loop 1: P + L4*u(theta4) = B + L5*u(theta5)
    Q = circle_intersections(P, L[4], B, L[5])[0 if cfg.branch_loop1 == 1 else 1]
    t4, t5 = angle(Q-P), angle(Q-B)

    # Loop 2: C + L9*u(theta9) = B + L10*u(theta10)
    D = circle_intersections(C, L[9], B, L[10])[0 if cfg.branch_loop2 == 1 else 1]
    t9, t10 = angle(D-C), angle(D-B)

    # Rigid triangle 2,3,4: link 2 points from P to top T1.
    T1 = triangle_vertex(P, Q, L[2], L[3], cfg.triangle1_side)
    t2 = angle(T1-P)
    t3 = angle(T1-Q)  # geometrical orientation from Q to apex

    # Rigid triangle 7,8,9: link 7 points from C to top T2,
    # and link 8 points from D to T2.
    T2 = triangle_vertex(C, D, L[7], L[8], cfg.triangle2_side)
    t7, t8 = angle(T2-C), angle(T2-D)

    link11_vector = T2-T1  # R11 required from T1 to T2
    required_L11 = float(np.linalg.norm(link11_vector))
    t11 = angle(link11_vector) if required_L11 > 1e-14 else float('nan')
    wing_position = (T1 + T2) / 2
    gap = required_L11-L[11]

    # Closed-loop vector residuals (loop 3 uses required R11 direction, desired L11).
    angles = (cfg.theta0_rad, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11)
    r1, r2, r3 = loop_residuals([length * unit(theta) for length, theta in zip(L, angles)])
    return dict(A=A, B=B, P=P, Q=Q, C=C, D=D, T1=T1, T2=T2,
                theta0=cfg.theta0_rad, theta1=t1, theta6=t6, theta4=t4, theta5=t5,
                theta2=t2, theta3=t3, theta7=t7, theta8=t8,
                theta9=t9, theta10=t10, theta11=t11,
                wing_position_m=wing_position, wing_orientation_rad=t11,
                required_L11=required_L11, closure_gap_m=gap,
                loop1_residual=r1, loop2_residual=r2, loop3_residual=r3)


def solve_driven_positions(cfg, theta6_rad, *, theta1_guess_rad=None):
    """Drive AC and determine AP from the fixed-length L11 constraint.

    Use the previous pose's AP angle when following a sweep on one branch.
    """
    guess = cfg.theta1_rad if theta1_guess_rad is None else theta1_guess_rad

    def positions(theta1):
        return solve_positions(cfg, theta1_rad=theta1, theta6_rad=theta6_rad)

    answer = root_scalar(lambda theta1: positions(theta1)['closure_gap_m'],
                         x0=guess, x1=guess + 1e-5, xtol=1e-12)
    if not answer.converged:
        raise ValueError('Could not close loop 3 on the selected assembly branch')
    pos = positions(answer.root)
    if not abs(pos['closure_gap_m']) <= max(1e-10, cfg.scale_m * 1e-9):
        raise ValueError('Loop 3 solution exceeds the numerical closure tolerance')
    return pos


def sweep_driven_positions(cfg, theta6_angles_rad, *, theta1_guess_rad=None):
    """Follow one assembly branch by seeding each pose with its predecessor."""
    guess = theta1_guess_rad
    for angle in theta6_angles_rad:
        pos = solve_driven_positions(cfg, angle, theta1_guess_rad=guess)
        guess = pos['theta1']
        yield pos


def feasibility_checks(cfg: Concept2Config):
    """Input-independent necessary geometric condition for loop 2."""
    L = cfg.L
    lower = abs(L[0]-L[10])
    upper = L[0]+L[10]
    moving_max = L[6]+L[9]
    moving_min = abs(L[6]-L[9])
    return {
        'loop2_always_impossible': moving_max < lower or moving_min > upper,
        'loop2_moving_max_m': moving_max,
        'loop2_ground_min_m': lower,
    }
