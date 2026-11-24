"""Strict real-data CSV and JSON configuration IO. Never silently use demo loads."""
import csv
import json
from pathlib import Path
import numpy as np
from concept2.prbm import reference_geometry
from .core import Material, Hinge, relative_rotation_rad, thickness_for_stiffness_m


def load_config(path):
    raw = json.loads(Path(path).read_text(encoding='utf-8'))
    for name, props in raw['materials'].items():
        if props.get('youngs_modulus_pa') is None:
            raise ValueError(f'{name}: supply youngs_modulus_pa')
    materials = {name:Material(name=name,**props) for name,props in raw['materials'].items()}
    hinges = []
    for item in raw['hinges']:
        data = dict(item)
        for key in ('length_m', 'width_m', 'rest_angle_rad', 'stiffness_factor'):
            if key in data and data[key] is None:
                raise ValueError(f"{data['name']}: supply {key}")
        mat = data.pop('material')
        if mat not in materials: raise ValueError(f'Unknown material {mat}')
        if data.get('thickness_m') is None:
            target = data.get('target_stiffness_nm_rad')
            if target is None:
                raise ValueError(f"{data['name']}: supply thickness_m or target_stiffness_nm_rad")
            data['thickness_m'] = thickness_for_stiffness_m(
                target, data['length_m'], data['width_m'], materials[mat],
                data.get('stiffness_factor', 1.0))
        hinges.append(Hinge(material=materials[mat],**data))
    if len({h.name for h in hinges}) != len(hinges): raise ValueError('Hinge names must be unique')
    if not hinges: raise ValueError('At least one hinge required')
    return hinges


def load_kinematics(path,hinges):
    with Path(path).open(newline='',encoding='utf-8-sig') as fp:
        rows = list(csv.DictReader(fp))
    if len(rows)<3: raise ValueError('At least 3 rows of validated mechanism kinematics required')
    required = ['q_rad','wing_x_m','wing_y_m','wing_fx_n','wing_fy_n']
    for k in required:
        if k not in rows[0]: raise ValueError(f'Missing column {k}; see inputs/kinematics_TEMPLATE.csv')
    def col(name): return np.array([float(row[name]) for row in rows])
    angles = {}
    for hinge in hinges:
        hinge_column = f'hinge_{hinge.name}_rad'
        a = f'theta_{hinge.link_a}_rad'
        b = f'theta_{hinge.link_b}_rad'
        if hinge_column in rows[0]:
            angles[hinge.name] = col(hinge_column)
        elif a in rows[0] and b in rows[0]:
            angles[hinge.name] = relative_rotation_rad(col(a), col(b))
        else:
            raise ValueError(f'{hinge.name}: supply {hinge_column} or both {a} and {b}')
    out = {'q_rad':col('q_rad'), 'wing_x_m':col('wing_x_m'), 'wing_y_m':col('wing_y_m'),
           'wing_fx_n':col('wing_fx_n'), 'wing_fy_n':col('wing_fy_n'),
           'hinge_angles_rad':angles}
    for name,key in [('wing_external_moment_nm','wing_external_moment_nm'),('wing_orientation_rad','wing_orientation_rad')]:
        if name in rows[0]: out[key]=col(name)
    return out


def load_layout(path, hinges):
    raw = json.loads(Path(path).read_text(encoding='utf-8'))
    if set(raw['virtual_springs']) != {hinge.name for hinge in hinges}:
        raise ValueError('Layout virtual spring names must match the configured hinges')
    return reference_geometry(raw['virtual_springs'], raw['rigid_links'], raw.get('fixed_points'))


def save_outputs(out_dir,result,records,layout=None,*,demo=False):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True,exist_ok=True)
    keys = ['q_rad','flexure_torque_nm','wing_torque_nm','total_torque_nm','elastic_energy_j']
    with (out_dir/'torque_sweep.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f); writer.writerow(keys)
        for i in range(len(result['q_rad'])):writer.writerow([result[key][i] for key in keys])
    with (out_dir/'hinge_checks.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    with (out_dir/'hinge_sweep.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['hinge', 'q_rad', 'relative_angle_rad', 'input_torque_nm', 'stress_pa', 'strain'])
        for hinge, angles in result['hinge_relative_angle_rad'].items():
            values = zip(result['q_rad'], angles, result['hinge_torque_nm'][hinge],
                         result['hinge_stress_pa'][hinge], result['hinge_strain'][hinge])
            writer.writerows((hinge, *row) for row in values)
    if layout is not None:
        (out_dir/'prbm_geometry.json').write_text(json.dumps(layout, indent=2) + '\n', encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,5))
    qdeg=np.rad2deg(result['q_rad'])
    for key,label in [('flexure_torque_nm','Flexures'),('wing_torque_nm','Wing loads'),('total_torque_nm','Total')]:
        ax.plot(qdeg,result[key],label=label)
    ax.set(xlabel='Input rotation [deg]',ylabel='Required actuator torque [N m]')
    if demo:
        ax.set_title('Synthetic PRBM demo')
    ax.legend();ax.grid(True);fig.tight_layout();fig.savefig(out_dir/'torque_plot.png',dpi=160);plt.close(fig)
