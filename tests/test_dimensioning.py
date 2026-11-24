import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from concept2.prbm import reference_geometry
from dimensioning.core import Material, Hinge, analyze, checks
from dimensioning.io import load_kinematics

class TestPRBM(unittest.TestCase):
    def setUp(self):
        self.mat=Material('test',100e9,allowable_stress_pa=100e6)
        self.h=Hinge('H','a','b',0.01,0.01,0.001,self.mat)
        self.q=np.linspace(-0.1,0.1,101)
    def test_stiffness(self):
        self.assertAlmostEqual(self.h.k_nm_rad,100e9*0.01*(0.001**3)/12/0.01)
    def test_spring_torque(self):
        q=self.q; z=np.zeros_like(q)
        r=analyze(q,[self.h],{'H':q},z,z,z,z)
        np.testing.assert_allclose(r['flexure_torque_nm'],self.h.k_nm_rad*q,atol=1e-12)
    def test_external_force(self):
        q=self.q; z=np.zeros_like(q)
        r=analyze(q,[],{},q,z,np.ones_like(q)*3,z)
        np.testing.assert_allclose(r['wing_torque_nm'],-3,atol=1e-12)
    def test_checks(self):
        q=self.q;z=np.zeros_like(q)
        r=analyze(q,[self.h],{'H':q},z,z,z,z)
        c=checks([self.h],r)[0]
        self.assertAlmostEqual(c['peak_strain'],0.005)
        self.assertIsNone(c['fatigue_amplitude_ratio'])
    def test_requires_orientation_for_external_moment(self):
        q=self.q;z=np.zeros_like(q)
        with self.assertRaises(ValueError): analyze(q,[],{},z,z,z,z,wing_external_moment_nm=np.ones_like(q))

    def test_target_thickness_and_static_limit(self):
        self.assertAlmostEqual(self.h.thickness_for_stiffness_m(self.h.k_nm_rad),
                               self.h.thickness_m)
        q = self.q
        z = np.zeros_like(q)
        result = analyze(q, [self.h], {'H': q}, z, z, z, z)
        self.assertAlmostEqual(checks([self.h], result)[0]['max_thickness_static_m'], 0.0002)

    def test_wrapped_link_angles_from_csv(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'sweep.csv'
            path.write_text(
                'q_rad,wing_x_m,wing_y_m,wing_fx_n,wing_fy_n,theta_a_rad,theta_b_rad\n'
                f'0,0,0,0,0,3.0,3.05\n'
                f'1,0,0,0,0,3.1,{3.25-2*np.pi}\n'
                f'2,0,0,0,0,{3.2-2*np.pi},{3.45-2*np.pi}\n',
                encoding='utf-8',
            )
            data = load_kinematics(path, [self.h])
        np.testing.assert_allclose(data['hinge_angles_rad']['H'], [0.05, 0.15, 0.25])

    def test_virtual_pivot_and_effective_link_length(self):
        springs = {'H': {'start_xy_m': [0, 0], 'end_xy_m': [0.02, 0],
                         'pivot_fraction': 0.25}}
        layout = reference_geometry(springs, {'link': ['base', 'H']}, {'base': [0, 0]})
        self.assertEqual(layout['virtual_spring_xy_m']['H'], [0.005, 0.0])
        self.assertAlmostEqual(layout['effective_link_length_m']['link'], 0.005)

if __name__=='__main__': unittest.main()
