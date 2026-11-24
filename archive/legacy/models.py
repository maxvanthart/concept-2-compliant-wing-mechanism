"""Kinematics / mechanics interface for one concept.

This example uses *one independent rotational coordinate*. Adapt the arrays and
solver for multiple DOFs and holonomic constraints when equations are available.
"""
import numpy as np


class ConceptModel:
    def __init__(self, config):
        self.config = config

    def position(self, q: np.ndarray) -> np.ndarray:
        """Return output x,y position [m] from generalized coordinates q [rad]."""
        raise NotImplementedError("Add your position equations: r_out(q).")

    def jacobian(self, q: np.ndarray) -> np.ndarray:
        """Return 2 x n Jacobian dr_out/dq [m/rad]."""
        raise NotImplementedError("Add the analytical/geometric Jacobian.")

    def generalized_force(self, q: np.ndarray, u: float) -> np.ndarray:
        """Equilibrium residual: internal generalized forces - applied forces.

        At static equilibrium: residual(q,u)=0.
        If using a closed linkage, include constraints in an augmented solver.
        """
        raise NotImplementedError("Add PRBM equilibrium equations.")

    def mass_matrix(self, q: np.ndarray) -> np.ndarray:
        """Optional: M(q) for future dynamic equations M*qdd + h = Q."""
        raise NotImplementedError("Add the dynamics only if needed.")


class DemoOneLinkModel(ConceptModel):
    """ONLY A DEMONSTRATION, NOT A MODEL OF YOUR TEAM'S CONCEPT."""
    def __init__(self, config):
        super().__init__(config)
        self.length_m = 0.10
        self.k_nm_rad = 0.1

    def position(self, q):
        theta = float(q[0])
        return np.array([self.length_m * np.cos(theta),
                         self.length_m * np.sin(theta)])

    def jacobian(self, q):
        theta = float(q[0])
        return np.array([[-self.length_m * np.sin(theta)],
                         [self.length_m * np.cos(theta)]])

    def generalized_force(self, q, u):
        return np.array([self.k_nm_rad * (q[0] - self.config.rest_angle_rad) - u])
