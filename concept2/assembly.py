"""Select the unloaded assembly pose for identical linear flexures."""
import numpy as np
from scipy.optimize import minimize_scalar

from .config import MOTIONGEN_THETA1_RAD, MOTIONGEN_THETA6_RAD
from .geometry import solve_driven_positions, sweep_driven_positions
from .prbm import FLEXURE_PAIRS, relative_flexure_angles


def minimum_energy_pose(cfg):
    """Return the lowest-energy pose and the full clockwise reference sweep.

    The shared positive stiffness cancels from the minimisation. Each hinge's
    rest angle is the midpoint of its unwrapped relative rotation range.
    """
    travel = np.linspace(0, 2 * np.pi, 721)
    poses = list(sweep_driven_positions(cfg, MOTIONGEN_THETA6_RAD - travel,
                                      theta1_guess_rad=MOTIONGEN_THETA1_RAD))
    if not np.allclose(poses[0]['P'], poses[-1]['P'], rtol=0, atol=1e-8):
        raise ValueError('The selected assembly branch does not repeat after one AC revolution.')

    pairs = tuple(FLEXURE_PAIRS.values())
    relative = np.column_stack([relative_flexure_angles(poses, a, b) for a, b in pairs])
    rest = (relative.min(axis=0) + relative.max(axis=0)) / 2
    # The first and last poses coincide. Use one of them in the periodic search.
    energy = 0.5 * np.sum((relative[:-1] - rest)**2, axis=1)
    candidates = np.flatnonzero((energy <= np.roll(energy, 1)) & (energy <= np.roll(energy, -1)))
    step = travel[1]

    def pose_energy(q):
        phase = q % (2 * np.pi)
        nearest = int(round(phase / step))
        pose = solve_driven_positions(cfg, MOTIONGEN_THETA6_RAD - phase,
                                      theta1_guess_rad=poses[nearest]['theta1'])
        angles = np.array([pose[f'theta{b}'] - pose[f'theta{a}'] for a, b in pairs])
        angles += 2 * np.pi * np.round((relative[nearest] - angles) / (2 * np.pi))
        deflection = angles - rest
        return pose, 0.5 * float(deflection @ deflection)

    best = None
    for index in candidates:
        result = minimize_scalar(lambda q: pose_energy(q)[1],
                                 bounds=(travel[index] - step, travel[index] + step),
                                 method='bounded', options={'xatol': 1e-12})
        if not result.success:
            raise ValueError('Could not refine the minimum-energy assembly pose.')
        if best is None or result.fun < best.fun:
            best = result
    if best is None:
        raise ValueError('No finite energy minimum was found over the AC revolution.')
    return pose_energy(best.x)[0], poses
