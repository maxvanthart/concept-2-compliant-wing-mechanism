# Concept 2 · Compliant wing mechanism

Run the following commands from the project directory with Python 3.10 or later.

```powershell
python -m pip install -r requirements.txt
python run_concept2.py --require-closed
python run_wing_sweep.py
python -B -m unittest discover -s tests -v
```

| Task | Entry point and output |
| --- | --- |
| Check loop closure and generate the PRBM schematic | `run_concept2.py` → `results/C2/` |
| Compare flexure designs interactively | Open `results/C2/flexure_table.html` in a browser |
| Calculate the wing path, rotation, and input torque | `run_wing_sweep.py` → `results/C2_WING_SWEEP/` |
| Size flexures using design inputs | `run_dimensioning.py`, with templates in `inputs/` |
| Find the original single-link demo | `archive/legacy/` |

## Reference geometry and loop closure

*Berekeningen voor Compliant.pdf* gives the first two loops and the rigid triangles on p. 1, followed by the third loop and the CAD table on p. 2. The ratios on p. 3 agree with the CAD series in link order. The CAD table on p. 2 lists the links in this order: 0, 10, 6, 1, 5, 9, 7, 8, 4, 2, 3, 11. The model scales all link lengths by `20 / 26.115`. L4 is 200 mm, L11 is 38.0777 mm, and the element width is 10 mm.

The `2-fourbars.motiongen` file has one motor at A. AB is the fixed reference and AC (L6) is the driven crank. AP (L1) follows from the fixed length of L11. `run_concept2.py` uses the saved MotionGen pose as the starting point for following the assembly branch. `--theta6-deg` sets the AC angle in degrees from the global +x axis. `--omega6` sets the signed angular speed in rad/s and defaults to −1 rad/s (clockwise).

| Loop | Vector equation |
| --- | --- |
| 1 | R1 + R4 − R0 − R5 = 0 |
| 2 | R6 + R9 − R0 − R10 = 0 |
| 3 | R1 + R2 + R11 − R8 − R10 − R0 = 0 |

With R0 as the reference, the CAD directions of R3, R5, R6, R8, and R9 were reversed to match the vector loops and triangles. `concept2/config.py` also records the original CAD angles. The fixed phase `theta6 - theta1 = -0.55859265 rad` belongs to the separate `--fixed-phase` diagnostic. For that diagnostic, the default branches `(+1,+1,+1,+1)` reproduce the directed CAD angles to within about 0.00081 rad.

| Position residual | Loop 1 | Loop 2 | Loop 3 |
| --- | ---: | ---: | ---: |
| Raw CAD angles, mm | 174.432 | 233.267 | 452.657 |
| Directed CAD angles, mm | 0.08841 | 0.04527 | 0.08276 |
| Fixed-phase diagnostic at θ1 = 2.405 rad, mm | < 1e-12 | < 1e-12 | 0.01945 |

The fixed-phase diagnostic imposes equal angular speeds on AP and AC. At θ1 = 2.405308738589 rad, its positions close, while the required L11 length changes by +63.008 mm per radian of input. With AC as the sole driver, AP rocks between approximately 99.4° and 190.4° over one revolution. The tests check position and velocity closure at 361 poses and verify that the wing returns to its initial position and orientation.

`run_concept2.py` reports the position and velocity closure residuals and produces a geometry plot. `--require-closed` returns a nonzero exit code if the pose or motion fails the numerical closure check. The motion check uses displacement per radian, including when the specified motor speed is zero.

## PRBM schematic and flexure selection

The pseudo-rigid-body model (PRBM) represents each candidate flexure as one torsion spring between two rigid bodies. Each triangle is a single rigid body. The virtual pivots lie at the MotionGen joints to preserve the nominal kinematics, so the effective lengths are the pivot-to-pivot distances. The free bending length of a physical flexure is a separate input.

| Pivot | Rotational interfaces | PRBM |
| --- | --- | --- |
| A | Ground–AP and ground–AC | Revolute joints |
| B | Ground–BQ and ground–BD | Two candidate torsion springs |
| P, Q | AP–triangle P–Q–T1 at P and BQ–the same triangle at Q | Candidate torsion springs |
| C | AC–triangle C–D–T2 | Revolute joint |
| D | BD–triangle C–D–T2 | Candidate torsion spring |
| T1, T2 | Each triangle–L11 at its apex | Candidate torsion springs |

AC makes a full revolution relative to the frame at A and relative to the second triangle at C. AP rocks about its own connection to the frame at A. The revolute joints at A and C carry zero elastic spring torque in this model. Ground–BQ and ground–BD have separate relative rotations at B.

For a closed, AC-driven pose, `python run_concept2.py --require-closed` generates the SVG and PNG schematics under `results/C2/`. The **Twist** column shows the relative angular travel of each candidate torsion spring over a full AC revolution, centred between its two extreme angles. The physical mounting determines the actual rest angle. `prbm_link_geometry.csv` lists the pivot coordinates, nominal link lengths, and element width in mm for CAD and ANSYS. A detailed flexure model also needs the bending-zone geometry, transition radii, material properties, bearings, and loads.

Open `results/C2/flexure_table.html` locally in a browser. Run `python run_concept2.py --require-closed` after a geometry change to update the angle limits from 721 poses. For each candidate flexure, enter the material, joint type, strip dimensions, strain limit, and any rest-angle offset. The width starts at 10 mm. Enter Young's modulus to calculate the linear torsional stiffness of a straight strip. Enter the strain limit as a percentage.

For a straight strip, the table assumes uniform pure bending. For a single-stage crossed-flexure pivot (CFP) with straight, directly clamped blades, the calculation gives a lower bound on peak bending strain per blade. The change in end-tangent angles is the integral of curvature, so the maximum absolute curvature is at least as large as its absolute mean. Enter the peak strain from ANSYS or a measurement over the full travel to check a CFP or compound joint against the chosen limit. Pivot drift and fatigue life require separate checks.

## Wing trajectory and torque sweep

The wing is mounted at the midpoint of T1–T2 (L11). Its location is `(T1 + T2) / 2`. Its angle is `atan2(T2_y - T1_y, T2_x - T1_x)`, positive counterclockwise from the global +x axis. `run_concept2.py` also reports the midpoint velocity and the angular velocity of L11.

`python run_wing_sweep.py` calculates 721 poses over one clockwise revolution of AC. Use `--counterclockwise` to reverse it. The output is `wing_kinematics.png`, `input_torque.png`, and `wing_sweep.csv` in `results/C2_WING_SWEEP/`. Zero travel is the saved MotionGen pose. In the CSV, `q_rad` runs from 0 to 2π in the chosen direction, and `ac_angle_rad` is the absolute AC angle measured counterclockwise from the global +x axis.

Without a flexure configuration, the upper torque curve shows `T_flex/k` for a common torsional stiffness. The rest angle of each of the seven candidate springs is centred in its own angular range. Use `--hinges inputs/hinges_C2.json` to calculate flexure torque in N·m with physical flexure inputs. Use L0 for ground, L4 and L9 for the two triangles, and L1, L5, L10, and L11 for the remaining bodies. The order of the two body labels determines the sign of the relative angle. Give the rest angle on the same unwrapped angle branch as the exported values.

The drag model assumes still air and a thin rectangular wing with a **40 mm chord along L11** and a **1.05 m total span**. The full planform area is **0.0420 m²**. The broadside factor is `c = rho*C_D*A/2`. The acute angle `alpha` between L11 and the relative airflow at the midpoint determines the projected area `A*abs(sin(alpha))`. The force on the midpoint is `F = -c*abs(sin(alpha))*|v|*v`, with `v = (dr/dq)*dq/dt`. The CSV includes the incidence angle, area fraction, and `wing_projected_area_m2`. The input torque required to overcome drag is `-F·dr/dq`. Positive torque acts in the chosen direction of rotation. [NASA Glenn, drag equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/)

The default drag curve is normalised by `c*(dq/dt)^2`. For drag torque in N·m, supply `--rpm` and either `--air-density` with `--drag-coefficient` or a direct `--drag-factor`. With both aerodynamic and flexure inputs, the CSV also contains the total quasi-static torque. This drag model covers resistance due to translation and projected area. For lift, profile drag, pitching moment, or a varying aerodynamic load, supply forces and moments for each pose to `run_dimensioning.py`.

## Flexure dimensioning

1. Copy `inputs/hinges_TEMPLATE.json`. For each flexure, enter its two adjacent bodies, free length, thickness, rest angle, justified stiffness factor, and material. The default width is 0.010 m. If the length is known and `target_stiffness_nm_rad` is specified, leave `thickness_m` blank to calculate a candidate thickness.
2. Provide at least three increasing values of `q_rad` with validated motion and wing loads. Supply a continuous relative angle in `hinge_<name>_rad` or absolute body angles in `theta_<link_a>_rad` and `theta_<link_b>_rad`. The tool unwraps jumps at ±π and calculates `theta_b - theta_a` for the absolute-angle input. Sample often enough that each body rotates by less than π between adjacent rows. Specify the rest angle on the same unwrapped branch as the first row.
3. Supply the midpoint in `wing_x_m` and `wing_y_m`, and the resultant force in `wing_fx_n` and `wing_fy_n`. For a resultant moment about the midpoint, add `wing_external_moment_nm` and the continuous L11 angle in `wing_orientation_rad`. Both columns can be omitted when the resultant moment is zero.
4. To calculate equivalent pivot locations and link lengths from a layout, pass `--layout path/to/layout.json`. Each entry in `virtual_springs` uses the flexure name and provides `start_xy_m`, `end_xy_m`, and the `pivot_fraction` from the selected PRBM. The coordinates share one planar reference frame, and the pivot is interpolated along the specified segment. In `rigid_links`, give each link two point names from the virtual springs or `fixed_points`, the latter specified as `[x, y]` pairs. Supply the kinematics of this geometry separately for the sweep.

```powershell
python run_dimensioning.py --config inputs/hinges_C2.json --kinematics inputs/kinematics_C2.csv --out results/C2_PRBM
```

| Output | Contents |
| --- | --- |
| `torque_sweep.csv`, `torque_plot.png` | Flexure, wing, and total input torque versus input angle |
| `hinge_sweep.csv` | Relative angle, input torque contribution, stress, and strain per flexure and pose |
| `hinge_checks.csv` | Stiffness, peak stress and strain, static safety factors, and nominal stress amplitude. Optional outputs include thickness for the target stiffness, static thickness limit at the same length and rotation, and the ratio to the supplied fatigue stress amplitude. |
| `prbm_geometry.json` | Virtual pivot positions and equivalent link lengths when a layout is supplied |

The quasi-static model uses `k = factor * E * b * t^3 / (12 * Lf)` and `U = sum(0.5*k*delta_theta^2)`. The actuator torque is `dU/dq - Fx*dx/dq - Fy*dy/dq - M*dphi/dq`. For a straight flexure in uniform pure bending, nominal surface strain is `t*delta_theta/(2*Lf)`, and linear-elastic stress is `E*strain`. The stiffness factor applies to the spring law. Notch stresses, combined loads, and nonlinear or time-dependent behaviour require separate models or measurements. [MIT Mechanics of Materials, pure bending](https://ocw.mit.edu/courses/2-001-mechanics-materials-i-fall-2006/1340543ced026dab5c8444a53fe92153_lec20.pdf)

A design-level result requires the physical flexure locations and body pairs, the selected flexure models, rest angles, strip dimensions or target stiffnesses, and the resultant wing force and moment throughout the stroke. Recheck loop closure if the virtual pivots change the effective lengths. Spring steel, PLA, and TPU require material data for the chosen production process and operating conditions, including allowable stress, strain, and fatigue data for the relevant mean stress and design life. The optional fatigue check compares nominal amplitude with the limit you supply.

`python run_dimensioning.py --demo` uses synthetic data and writes to `results/PRBM_DEMO_NOT_CONCEPT2/` by default. Run the original single-link demo from `archive/legacy/` if needed.

## Source files

- `concept2/config.py` records the CAD lengths, angles, and reference pose.
- `concept2/geometry.py` solves position closure. `concept2/velocity.py` differentiates the loop constraints.
- `concept2/prbm.py` defines the candidate flexure pairs, relative angles, and reference pivot geometry.
- `concept2/schematic.py` generates the PRBM drawing, link geometry CSV, and interactive table from `concept2/flexure_table.html`.
- `dimensioning/core.py` calculates relative rotations, thickness from target stiffness, nominal stress and strain, and input torque. `dimensioning/io.py` reads the inputs and exports the results.
- `inputs/hinges_TEMPLATE.json` leaves unknown design values blank. `inputs/kinematics_TEMPLATE.csv` gives the expected CSV layout.
- `tests/` checks position and velocity closure, angle unwrapping, virtual pivot geometry, and inverse sizing.

