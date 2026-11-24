"""Concept 2 geometry from the annotated CAD table supplied 17 Nov 2026.

Link indices L0..L11 are explicitly mapped; all dimensions are in metres.
Raw CAD angles and the inferred vector directions are stored separately.
"""
from dataclasses import dataclass
import numpy as np

# Row order in the screenshot is NOT index order. Values below are correctly
# reordered by link label, rather than by row position.
LENGTHS_CM = np.array([
    12.637, 8.553, 23.892, 29.981,
    26.115, 11.394, 3.060, 27.323,
    29.558, 16.344, 6.884, 4.972,
], dtype=float)

# These angles were copied as displayed. A CAD package may report different
# orientation conventions than the signed vectors in our loop equations.
CAD_ANGLES_RAD = np.array([
    5.814, 2.405, 1.144, 5.274,
    6.144, 3.901, 4.988, 1.227,
    4.946, 2.955, 0.958, 0.237,
], dtype=float)

# Reversing these vectors matches the PDF loop and triangle directions to
# within the precision of the rounded CAD table. R0 fixes the reference frame.
LOOP_REFERENCE_ANGLES_RAD = CAD_ANGLES_RAD.copy()
LOOP_REFERENCE_ANGLES_RAD[[3, 5, 6, 8, 9]] += np.pi
LOOP_REFERENCE_ANGLES_RAD %= 2 * np.pi

# Saved joint directions in 2-fourbars.motiongen. The lengths used by the
# project still come from the supplied rounded CAD table.
MOTIONGEN_THETA1_RAD = np.arctan2(4.484785 - (-0.963279), -11.203188 - (-4.610042))
MOTIONGEN_THETA6_RAD = np.arctan2(-3.695073 - (-0.963279), -3.231620 - (-4.610042))

# Full mechanism is uniformly rescaled to L4 = 20 cm.
TARGET_L4_CM = 20.0
SCALE_FACTOR = TARGET_L4_CM / LENGTHS_CM[4]
SCALED_LENGTHS_CM = LENGTHS_CM * SCALE_FACTOR
# Width of the rigid elements; NOT automatically a validated flexure cross-section.
ELEMENT_WIDTH_M = 0.01
RATIOS_TO_L0 = LENGTHS_CM / LENGTHS_CM[0]

@dataclass
class Concept2Config:
    # A convenient normalization (L0=1). Set scale_m=physical L0 in metres later.
    scale_m: float = SCALED_LENGTHS_CM[0] / 100.0
    # Ground link angle is written as 5.814 rad in the page-2 table.
    theta0_rad: float = CAD_ANGLES_RAD[0]
    # Rounded CAD table reference pose.
    theta1_rad: float = CAD_ANGLES_RAD[1]
    # Constant phase for the legacy equal-speed PDF diagnostic.
    theta6_minus_theta1_rad: float = LOOP_REFERENCE_ANGLES_RAD[6] - LOOP_REFERENCE_ANGLES_RAD[1]
    omega1_rad_s: float = 1.0  # example input angular velocity
    branch_loop1: int = 1  # +1 or -1 for the two circle-intersection branches
    branch_loop2: int = 1
    triangle1_side: int = 1
    triangle2_side: int = 1

    @property
    def L(self) -> np.ndarray:
        return RATIOS_TO_L0 * self.scale_m
