"""Quasi-static torsion-spring PRBM, loads, and preliminary beam stress/fatigue.

Angles: rad; lengths: m; forces: N; torque: N*m; E and stress: Pa.
No claim that this model predicts all local notch stresses or 3D printed fatigue.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Material:
    name: str
    youngs_modulus_pa: float
    allowable_stress_pa: float | None = None
    allowable_strain: float | None = None
    fatigue_amplitude_allowable_pa: float | None = None
    # Fatigue limit above is user supplied at the stated number of cycles,
    # mean stress, finish and manufacturing condition; not a universal limit.

    def __post_init__(self):
        if not np.isfinite(self.youngs_modulus_pa) or self.youngs_modulus_pa <= 0:
            raise ValueError('E must be finite and positive')
        for name in ('allowable_stress_pa','allowable_strain','fatigue_amplitude_allowable_pa'):
            value = getattr(self,name)
            if value is not None and (not np.isfinite(value) or value <= 0):
                raise ValueError(f'{name} must be finite and positive')


def thickness_for_stiffness_m(target_k_nm_rad, length_m, width_m, material, stiffness_factor):
    """Invert the spring equation for a uniform rectangular flexure."""
    values = np.array([target_k_nm_rad, length_m, width_m, stiffness_factor])
    if not np.all(np.isfinite(values) & (values > 0)):
        raise ValueError('Target stiffness, length, width and stiffness factor must be finite and positive')
    return (12 * length_m * target_k_nm_rad /
            (stiffness_factor * material.youngs_modulus_pa * width_m)) ** (1/3)


@dataclass(frozen=True)
class Hinge:
    name: str
    link_a: str
    link_b: str
    length_m: float
    width_m: float
    thickness_m: float
    material: Material
    stiffness_factor: float = 1.0
    # Relative rest angle is established from the assembly's unforced geometry.
    rest_angle_rad: float = 0.0
    target_stiffness_nm_rad: float | None = None

    def __post_init__(self):
        values = np.array([self.length_m, self.width_m, self.thickness_m, self.stiffness_factor])
        if not np.all(np.isfinite(values) & (values > 0)):
            raise ValueError(f'{self.name}: dimensions and stiffness_factor must be finite and positive')
        if not np.isfinite(self.rest_angle_rad):
            raise ValueError(f'{self.name}: rest angle must be finite')
        if self.target_stiffness_nm_rad is not None:
            self.thickness_for_stiffness_m(self.target_stiffness_nm_rad)

    @property
    def second_moment_m4(self): return self.width_m * self.thickness_m**3 / 12

    @property
    def k_nm_rad(self):
        return (self.stiffness_factor * self.material.youngs_modulus_pa *
                self.second_moment_m4 / self.length_m)

    def stress_pa(self, delta_theta_rad):
        """Straight prismatic beam, uniform pure bending; NOT notch peak stress."""
        # Moment per beam in this simplified model: E I delta_theta/L.
        return self.material.youngs_modulus_pa * self.thickness_m * np.asarray(delta_theta_rad) / (2*self.length_m)

    def strain(self, delta_theta_rad):
        return self.thickness_m * np.asarray(delta_theta_rad) / (2*self.length_m)

    def thickness_for_stiffness_m(self, target_k_nm_rad):
        """Invert the spring equation at the supplied flexure length and width."""
        return thickness_for_stiffness_m(target_k_nm_rad, self.length_m, self.width_m,
                                         self.material, self.stiffness_factor)


def _as_1d(values, n=None, name='vector'):
    a = np.asarray(values, float)
    if a.ndim != 1 or (n is not None and len(a) != n) or not np.all(np.isfinite(a)):
        raise ValueError(f'{name} must be finite 1D array of length {n}')
    return a


def relative_rotation_rad(theta_a_rad, theta_b_rad):
    """Continuous theta_b - theta_a; adjacent samples must resolve each rotation."""
    a = _as_1d(theta_a_rad, name='link A rotation')
    b = _as_1d(theta_b_rad, len(a), 'link B rotation')
    return np.unwrap(b) - np.unwrap(a)


def analyze(q_rad, hinges, hinge_angles_rad, wing_x_m, wing_y_m,
            wing_force_x_n, wing_force_y_n, wing_external_moment_nm=None,
            wing_orientation_rad=None):
    """Find actuator torque resisting spring and applied wing load (quasi-static).

    hinge_angles_rad: mapping hinge name -> relative continuous unwrapped angle array.
    q_rad: monotonically increasing generalized input rotation.
    External forces point in actual direction they act ON the wing.
    Wing positions/orientation and hinge angles must be from a valid closed mechanism.
    Returns T_total = dU/dq - F.dot(dr/dq) - M*(dphi/dq).
    """
    q = _as_1d(q_rad, name='q')
    if len(q) < 3 or not np.all(np.diff(q) > 0):
        raise ValueError('q must contain >=3 strictly increasing input angles [rad]')
    n = len(q)
    x = _as_1d(wing_x_m,n,'wing x')
    y = _as_1d(wing_y_m,n,'wing y')
    fx = _as_1d(wing_force_x_n,n,'wing Fx')
    fy = _as_1d(wing_force_y_n,n,'wing Fy')
    moment = np.zeros(n) if wing_external_moment_nm is None else _as_1d(wing_external_moment_nm,n,'wing moment')
    if np.any(moment) and wing_orientation_rad is None:
        raise ValueError('Wing orientation required when external moments are nonzero')
    phi = np.zeros(n) if wing_orientation_rad is None else _as_1d(wing_orientation_rad,n,'wing orientation')
    edge_order = 2
    dx_dq = np.gradient(x,q,edge_order=edge_order)
    dy_dq = np.gradient(y,q,edge_order=edge_order)
    dphi_dq = np.gradient(phi,q,edge_order=edge_order)
    t_hinges = {}
    stresses = {}
    strains = {}
    relative_angles = {}
    t_flex = np.zeros(n)
    energy = np.zeros(n)
    for h in hinges:
        if h.name in t_hinges:
            raise ValueError(f'Duplicate hinge {h.name}')
        if h.name not in hinge_angles_rad: raise ValueError(f'Missing relative angle for hinge {h.name}')
        theta = _as_1d(hinge_angles_rad[h.name],n,h.name)
        relative_angles[h.name] = theta
        delta = theta - h.rest_angle_rad
        dtheta_dq = np.gradient(theta,q,edge_order=edge_order)
        t = h.k_nm_rad * delta * dtheta_dq
        t_hinges[h.name] = t
        t_flex += t
        energy += 0.5*h.k_nm_rad*delta**2
        stresses[h.name] = h.stress_pa(delta)
        strains[h.name] = h.strain(delta)
    t_wing = -(fx*dx_dq + fy*dy_dq + moment*dphi_dq)
    return {'q_rad':q,'hinge_torque_nm':t_hinges,'flexure_torque_nm':t_flex,
            'wing_torque_nm':t_wing,'total_torque_nm':t_flex+t_wing,
            'elastic_energy_j':energy,'hinge_stress_pa':stresses,
            'hinge_strain':strains,'hinge_relative_angle_rad':relative_angles}


def checks(hinges, result):
    """Screen max/min stress, strain, and stress amplitude; no fatigue life inference."""
    records = []
    for h in hinges:
        s = result['hinge_stress_pa'][h.name]
        e = result['hinge_strain'][h.name]
        sig_peak = float(np.max(np.abs(s)))
        eps_peak = float(np.max(np.abs(e)))
        sig_a = float((max(s)-min(s))/2)
        sig_m = float((max(s)+min(s))/2)
        m = h.material
        stress_sf = None if m.allowable_stress_pa is None else (m.allowable_stress_pa / sig_peak if sig_peak else np.inf)
        strain_sf = None if m.allowable_strain is None else (m.allowable_strain / eps_peak if eps_peak else np.inf)
        fatigue_ratio = None if m.fatigue_amplitude_allowable_pa is None else sig_a/m.fatigue_amplitude_allowable_pa
        max_thickness = min((factor * h.thickness_m for factor in (stress_sf, strain_sf)
                             if factor is not None and np.isfinite(factor)), default=None)
        target_thickness = (None if h.target_stiffness_nm_rad is None
                            else h.thickness_for_stiffness_m(h.target_stiffness_nm_rad))
        records.append({'hinge':h.name,'material':m.name,'k_nm_rad':h.k_nm_rad,
                        'peak_stress_mpa':sig_peak/1e6,'peak_strain':eps_peak,
                        'mean_stress_mpa':sig_m/1e6,'stress_amplitude_mpa':sig_a/1e6,
                        'static_stress_safety_factor':stress_sf,
                        'strain_safety_factor':strain_sf,
                        'max_thickness_static_m':max_thickness,
                        'thickness_for_target_stiffness_m':target_thickness,
                        'fatigue_amplitude_ratio':fatigue_ratio})
    return records
