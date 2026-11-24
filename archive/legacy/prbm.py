"""Generic PRBM helper functions. Units are SI."""
import numpy as np
from scipy.optimize import root


def rotational_spring_torque(theta_rad, theta_rest_rad, k_nm_rad):
    return k_nm_rad * (np.asarray(theta_rad) - theta_rest_rad)


def rotational_spring_energy(theta_rad, theta_rest_rad, k_nm_rad):
    return 0.5 * k_nm_rad * (np.asarray(theta_rad) - theta_rest_rad)**2


def solve_static(model, u, q0):
    """Solve equilibrium for independent q. Extend for constrained systems."""
    answer = root(lambda q: model.generalized_force(q, u), np.asarray(q0, dtype=float))
    if not answer.success or not np.all(np.isfinite(answer.x)):
        raise RuntimeError(f"Equilibrium solve failed at u={u}: {answer.message}")
    return answer.x


def output_velocity(model, q, qdot):
    """Velocity relationship: rdot_out = J(q) @ qdot."""
    return model.jacobian(np.asarray(q)) @ np.asarray(qdot)
