import unittest
from dataclasses import replace
import numpy as np
from scipy.optimize import brentq
from concept2.config import (
    Concept2Config, LENGTHS_CM, SCALED_LENGTHS_CM, ELEMENT_WIDTH_M,
    LOOP_REFERENCE_ANGLES_RAD, MOTIONGEN_THETA1_RAD, MOTIONGEN_THETA6_RAD,
)
from concept2.geometry import loop_residuals, solve_driven_positions, solve_positions, unit
from concept2.velocity import solve_driven_velocities, solve_velocities

class TestScaledGeometry(unittest.TestCase):
    def test_length_l4(self):
        self.assertAlmostEqual(SCALED_LENGTHS_CM[4], 20.0)
    def test_all_ratios_preserved(self):
        np.testing.assert_allclose(SCALED_LENGTHS_CM / SCALED_LENGTHS_CM[4], LENGTHS_CM / LENGTHS_CM[4])
    def test_config_matches_scaled_lengths(self):
        np.testing.assert_allclose(Concept2Config().L * 100, SCALED_LENGTHS_CM)
    def test_width(self):
        self.assertEqual(ELEMENT_WIDTH_M, 0.01)

    def test_directed_cad_table_agrees_with_pdf_loops(self):
        cfg = Concept2Config()
        vectors = [length * unit(theta)
                   for length, theta in zip(cfg.L, LOOP_REFERENCE_ANGLES_RAD)]
        self.assertTrue(np.all(np.linalg.norm(loop_residuals(vectors), axis=1) < 1e-4))
        self.assertLess(np.linalg.norm(vectors[3] + vectors[4] - vectors[2]), 1e-4)
        self.assertLess(np.linalg.norm(vectors[7] - vectors[9] - vectors[8]), 1e-4)

    def test_pose_closure_does_not_imply_rigid_motion(self):
        cfg = Concept2Config()

        def positions_at(q):
            return solve_positions(replace(cfg, theta1_rad=q))

        q = brentq(lambda q: positions_at(q)['closure_gap_m'], 2.4, 2.41)
        cfg_at_root = replace(cfg, theta1_rad=q)
        pos = positions_at(q)
        vel = solve_velocities(cfg_at_root, pos)
        for i in (1, 2, 3):
            self.assertLess(np.linalg.norm(pos[f'loop{i}_residual']), 1e-10)
        for i in (1, 2):
            self.assertLess(np.linalg.norm(vel[f'loop{i}_velocity_residual']), 1e-10)
        self.assertGreater(np.linalg.norm(vel['loop3_velocity_residual']), 0.06)
        step = 1e-5
        lo, hi = positions_at(q-step), positions_at(q+step)
        for key in ('T1', 'T2'):
            np.testing.assert_allclose((hi[key]-lo[key])/(2*step),
                                       vel[f'velocity_{key}_m_s'], atol=1e-9)
        self.assertAlmostEqual((hi['required_L11']-lo['required_L11'])/(2*step),
                               vel['required_L11_rate_m_s'], places=8)

    def test_wing_output_is_midpoint_and_directed_bar_angle(self):
        cfg = Concept2Config()
        pos = solve_positions(cfg)
        vel = solve_velocities(cfg, pos)
        np.testing.assert_allclose(pos['wing_position_m'], (pos['T1'] + pos['T2']) / 2)
        span = pos['T2'] - pos['T1']
        np.testing.assert_allclose(unit(pos['wing_orientation_rad']), span / np.linalg.norm(span))
        step = 1e-5
        lo = solve_positions(replace(cfg, theta1_rad=cfg.theta1_rad - step))
        hi = solve_positions(replace(cfg, theta1_rad=cfg.theta1_rad + step))
        np.testing.assert_allclose((hi['wing_position_m'] - lo['wing_position_m']) / (2 * step),
                                   vel['wing_velocity_m_s'], atol=1e-9)
        self.assertAlmostEqual((hi['wing_orientation_rad'] - lo['wing_orientation_rad']) / (2 * step),
                               vel['wing_omega_rad_s'], places=8)

    def test_motiongen_driver_closes_through_a_full_turn(self):
        cfg = Concept2Config()
        guess = MOTIONGEN_THETA1_RAD
        for q in np.linspace(MOTIONGEN_THETA6_RAD, MOTIONGEN_THETA6_RAD + 2 * np.pi, 361):
            pos = solve_driven_positions(cfg, q, theta1_guess_rad=guess)
            vel = solve_driven_velocities(cfg, pos, 1.0)
            for loop in (1, 2, 3):
                self.assertLess(np.linalg.norm(pos[f'loop{loop}_residual']), 1e-10)
                self.assertLess(np.linalg.norm(vel[f'loop{loop}_velocity_residual']), 1e-10)
            guess = pos['theta1']
        first = solve_driven_positions(cfg, MOTIONGEN_THETA6_RAD, theta1_guess_rad=MOTIONGEN_THETA1_RAD)
        np.testing.assert_allclose(pos['wing_position_m'], first['wing_position_m'], atol=1e-10)
        self.assertAlmostEqual(pos['wing_orientation_rad'], first['wing_orientation_rad'], places=10)
        self.assertAlmostEqual(guess, first['theta1'], places=10)
        self.assertLess(solve_driven_velocities(cfg, first, 1.)['omega1'], 0.)

    def test_ac_driven_wing_velocity_matches_differenced_positions(self):
        cfg = Concept2Config()
        q = MOTIONGEN_THETA6_RAD
        step = 1e-5
        pos, lo, hi = (solve_driven_positions(cfg, angle, theta1_guess_rad=MOTIONGEN_THETA1_RAD)
                       for angle in (q, q-step, q+step))
        vel = solve_driven_velocities(cfg, pos, 1.)
        np.testing.assert_allclose((hi['wing_position_m'] - lo['wing_position_m']) / (2*step),
                                   vel['wing_velocity_m_s'], atol=1e-9)
        self.assertAlmostEqual((hi['wing_orientation_rad'] - lo['wing_orientation_rad']) / (2*step),
                               vel['wing_omega_rad_s'], places=8)
