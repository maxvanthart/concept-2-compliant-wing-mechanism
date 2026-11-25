import unittest
from dataclasses import replace
import numpy as np
from scipy.optimize import brentq
from concept2.assembly import minimum_energy_pose
from concept2.config import (
    Concept2Config, LENGTHS_CM, SCALED_LENGTHS_CM, ELEMENT_WIDTH_M,
    LOOP_REFERENCE_ANGLES_RAD, MOTIONGEN_THETA1_RAD, MOTIONGEN_THETA6_RAD,
)
from concept2.fabrication import NEUTRAL_CLAMPS, clamp_offsets, rotate
from concept2.geometry import loop_residuals, solve_driven_positions, solve_positions, unit
from concept2.prbm import FLEXURE_PAIRS, relative_flexure_angles
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

    def test_cfp_clamps_align_at_neutral_and_face_inboard_at_coupler(self):
        for joint in FLEXURE_PAIRS:
            neutral = np.deg2rad(-23)
            parent_offset, child_offset = clamp_offsets(joint, neutral)
            self.assertAlmostEqual(neutral + child_offset - parent_offset, 0)
        centre = NEUTRAL_CLAMPS[1].mean(axis=0)
        for joint, direction in (('T1', 1), ('T2', -1)):
            _, child_offset = clamp_offsets(joint, neutral)
            x, y = rotate(centre, child_offset)
            self.assertGreater(direction * x, 0)
            self.assertAlmostEqual(y, 0)

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

    def test_minimum_energy_pose_is_a_constrained_stationary_point(self):
        cfg = Concept2Config()
        pos, poses = minimum_energy_pose(cfg)
        q_deg = np.rad2deg(MOTIONGEN_THETA6_RAD - pos['theta6'])
        self.assertAlmostEqual(q_deg, 21.39775, places=4)
        pairs = tuple(FLEXURE_PAIRS.values())
        relative = np.column_stack([relative_flexure_angles(poses, a, b) for a, b in pairs])
        rest = (relative.min(axis=0) + relative.max(axis=0)) / 2
        angles = np.array([pos[f'theta{b}'] - pos[f'theta{a}'] for a, b in pairs])
        delta = (angles - rest + np.pi) % (2 * np.pi) - np.pi
        vel = solve_driven_velocities(cfg, pos, -1.0)
        rates = {i: vel[f'omega{i}'] for i in (0, 1, 4, 5, 9, 10)}
        rates[11] = vel['omega11_geometric']
        derivative = np.array([rates[b] - rates[a] for a, b in pairs])
        self.assertLess(abs(delta @ derivative), 1e-6)
        self.assertLess(np.sum(delta**2), np.min(np.sum((relative - rest)**2, axis=1)))
        for i in (1, 2, 3):
            self.assertLess(np.linalg.norm(pos[f'loop{i}_residual']), 1e-10)
