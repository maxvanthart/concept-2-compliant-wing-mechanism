"""Per-concept input parameters. Replace all None entries with measured/derived values."""
from dataclasses import dataclass
import numpy as np

@dataclass
class ConceptConfig:
    name: str
    # Dimensions and beam design data (SI units)
    beam_length_m: float | None = None
    beam_width_m: float | None = None
    beam_thickness_m: float | None = None
    youngs_modulus_pa: float | None = None
    # PRBM (values are model specific; do not guess coefficients)
    characteristic_radius_factor: float | None = None
    torsional_stiffness_nm_rad: float | None = None
    rest_angle_rad: float = 0.0
    # Solver / actuation
    initial_q_rad: float = 0.0
    input_torque_min_nm: float = -0.01
    input_torque_max_nm: float = 0.01
    samples: int = 101

    def second_moment_m4(self):
        if self.beam_width_m is None or self.beam_thickness_m is None:
            return None
        return self.beam_width_m * self.beam_thickness_m**3 / 12.0


def get_config(name: str) -> ConceptConfig:
    return ConceptConfig(name=name.upper())
