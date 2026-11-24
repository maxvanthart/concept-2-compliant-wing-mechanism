"""Run: python run_dimensioning.py --demo OR --config ... --kinematics ..."""
import argparse
from pathlib import Path
import numpy as np
from dimensioning.core import analyze,checks,Material,Hinge
from dimensioning.io import load_config,load_kinematics,load_layout,save_outputs


def demonstration():
    """Invented ONE-hinge mechanism for software testing, NOT actual Concept 2."""
    steel=Material('STEEL_DEMO',youngs_modulus_pa=200e9,allowable_stress_pa=250e6)
    hinges=[Hinge('H_DEMO','ground','link1',0.015,0.010,0.0003,steel)]
    q=np.linspace(-0.15,0.15,121)
    return hinges,dict(q_rad=q,hinge_angles_rad={'H_DEMO':q},
        wing_x_m=0.1*np.cos(q),wing_y_m=0.1*np.sin(q),
        wing_fx_n=np.zeros_like(q),wing_fy_n=-2*np.ones_like(q))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--demo',action='store_true',help='Generate illustrative test results; NOT real Concept 2 results')
    parser.add_argument('--config',default='inputs/hinges_TEMPLATE.json')
    parser.add_argument('--kinematics',default='inputs/kinematics_TEMPLATE.csv')
    parser.add_argument('--out',default=None)
    parser.add_argument('--layout', help='JSON with reference flexure endpoints and rigid link connectivity')
    args=parser.parse_args()
    if args.demo:
        hinges,data=demonstration()
        folder=args.out or 'results/PRBM_DEMO_NOT_CONCEPT2'
        print('DEMO: synthetic one-hinge geometry and forces. NOT validated Concept 2.')
    else:
        try:
            hinges=load_config(args.config)
            data=load_kinematics(args.kinematics,hinges)
        except (ValueError, OSError, KeyError, TypeError) as exc:
            parser.error(str(exc))
        folder=args.out or 'results/PRBM_CONCEPT2_UNVALIDATED'
        print('INPUT DATA: this program cannot verify physical closure or material test data.')
    result=analyze(data['q_rad'],hinges,data['hinge_angles_rad'],data['wing_x_m'],data['wing_y_m'],
        data['wing_fx_n'],data['wing_fy_n'],data.get('wing_external_moment_nm'),data.get('wing_orientation_rad'))
    records=checks(hinges,result)
    layout = load_layout(args.layout, hinges) if args.layout else None
    save_outputs(folder,result,records,layout,demo=args.demo)
    print(f"Peak |total torque|: {np.max(np.abs(result['total_torque_nm'])):.6g} N m")
    for row in records: print(row)
    print('Results:',Path(folder).resolve())


if __name__=='__main__': main()
