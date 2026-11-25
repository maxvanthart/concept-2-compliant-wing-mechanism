# Concept 2 · Compliant wing mechanism

This TU Delft project turns a rigid MotionGen wing linkage into a printable
mechanism with seven crossed-flexure pivots (CFPs). The driven crank AC
makes one full turn. The wing is mounted at the midpoint of T1–T2 (L11),
with its chord directed from T1 to T2.

## Run the project

Run these commands from the project directory with Python 3.10 or newer.

```powershell
python -m pip install -r requirements.txt
python run_concept2.py --require-closed
python run_wing_sweep.py
python -m unittest discover -s tests
```

| Tool | Use it to |
| --- | --- |
| [`run_concept2.py`](run_concept2.py) | Check position and velocity closure, find the default assembly pose, and regenerate the schematics and fabrication frames in `results/C2/`. `--theta6-deg` selects an absolute AC angle from global +x. `--omega6` defaults to −1 rad/s (clockwise). |
| [`run_flexure_dashboard.py`](run_flexure_dashboard.py) | Configure one CFP and preview its motion. **Check configuration** plots peak stress and reports cyclic stress amplitude. **Fit stiffness** returns a linear fit over the selected angle range. **Save CAD PNG** exports a dimensioned neutral sketch. |
| [`run_crossed_flexure.py`](run_crossed_flexure.py) | Run the isolated nonlinear CFP model and export stress, strain, moment and pivot drift to `results/CFP/`. |
| [`run_wing_sweep.py`](run_wing_sweep.py) | Plot the wing path, orientation and input torque over 721 poses. `--counterclockwise` reverses the sweep. The default outputs are `wing_kinematics.png`, `input_torque.png` and `wing_sweep.csv` in `results/C2_WING_SWEEP/`. |
| [`run_dimensioning.py`](run_dimensioning.py) | Size flexures and evaluate supplied wing loads. Fill the [hinge JSON](inputs/hinges_TEMPLATE.json) and [kinematics CSV](inputs/kinematics_TEMPLATE.csv), then pass them with `--config` and `--kinematics`. Outputs include torque, hinge stress and strain, sizing and optional fatigue comparisons. `--layout` accepts a virtual-pivot layout. `--demo` runs a synthetic example. |
| `run_concept2.py --step` | Export the four relaxed fabrication groups as coloured AP214 surfaces. See [CAD handoff](#5-relaxed-fabrication-and-cad-handoff). |

## Design process

### 1. Establish the rigid kinematics

The calculation notes define the CAD dimensions and three vector loops
[1, pp. 1–3]. The MotionGen file supplies the starting pose and AC driver
[2]. Scale every link by `20 / 26.115` to make L4 **200 mm**. The resulting
L11 length is **38.078 mm** and the reference element width is **10 mm**.

Solve all three loops with AC prescribed and AP free to rock. The fixed
length of L11 determines AP’s angle. Differentiating the loop equations
checks that the velocities close as well. This preserves the original
motion before hinge stiffness is introduced. The [PRBM schematic](results/C2/prbm_schematic.svg)
shows the closed geometry, and the [link CSV](results/C2/prbm_link_geometry.csv)
exports the pivot coordinates and lengths. The legacy `--fixed-phase`
diagnostic illustrates why imposing equal speeds on AP and AC is
incompatible with the closed mechanism.

Keep **A–P**, **A–C** and the joint at **C** as pin-and-hole interfaces.
A–P and A–C share the A axis and can be stacked. AC turns continuously
relative to its neighbours, and keeping AP on the same pin simplifies
assembly.

### 2. Choose the flexures and their rest angles

Use the lumped pseudo-rigid-body model (PRBM) [3, Ch. 5]. Each CFP becomes
one torsion spring at the original pivot, so the nominal rigid lengths
remain pivot-to-pivot distances. The crossed blades fit a longer bending
length into a compact joint and guide motion near the original pivot.
Using one common geometry simplifies fabrication and gives the seven
joints a shared moment–angle curve [4, Fig. 2].

For each parent–child pair, unwrap the signed relative angle
$\delta_i=\theta_b-\theta_a$ over one revolution. Choose the neutral
angle halfway between its extrema. This minimises the largest absolute
twist of that joint for symmetric travel limits. The table shows CAD
angles wrapped to ±180°. The [neutral-angle schematic](results/C2/prbm_neutral_angles.svg)
also lists the unwrapped values for torque calculations.

| Joint | Parent → child axis | Neutral angle | Required travel |
| --- | --- | ---: | ---: |
| B–BQ | L0 → L5 | +64.35° | ±21.83° |
| B–BD | L0 → L10 | +56.64° | ±43.50° |
| P | L1 → L4 | −151.26° | ±26.71° |
| Q | L5 → L4 | −46.59° | ±5.85° |
| D | L10 → L9 | −40.22° | ±32.16° |
| T1 | L4 → L11 | −7.06° | ±34.07° |
| T2 | L9 → L11 | −7.91° | ±38.35° |

Each blade sees the full clamp-to-clamp rotation. The pure-bending
relation gives $\varepsilon=t\delta/(2L)$ for uniform curvature [5].
The isolated nonlinear two-blade model resolves the nonuniform strain,
holding moment and pivot drift. The **reference PETG design** uses a
30 mm free blade length, 0.8 mm thickness, 4.5 mm depth per blade and a
90° crossing. B–BD sets the largest required range, **±43.5°**.
Use the neutral-angle schematic for the required travel and the CFP
dashboard to evaluate a hinge geometry against that travel.

### 3. Orient the mounts and choose the assembly pose

Place the neutral CFP crossing at its virtual pivot. If $\gamma_a$
and $\gamma_b$ are the clamp angles measured from their respective body
axes, choose them so $\gamma_b-\gamma_a=-\delta_{0,i}$. This makes the
clamps parallel at the joint’s neutral angle. At T1 and T2 the coupler
clamps face inward from opposite ends, perpendicular to L11, so they fit
between the two pivots. The [placement schematic](results/C2/prbm_cfp_placement.svg)
shows the actual clamp directions in the assembly pose.

The hinge rest angles are independent. To choose one assembly pose, let
$q$ be clockwise AC travel and $e_i(q)=\delta_i(q)-\delta_{0,i}$. With
$\mathcal L=T-U$, the static Lagrange equation gives $U'(q)=Q_q$ [6].
For a common linear stiffness $k>0$ and zero applied input moment,

$$
U(q)=\tfrac12 k\sum_i e_i(q)^2,\qquad
U'(q)=k\sum_i e_i(q)\,[\theta_b'(q)-\theta_a'(q)]=0.
$$

The differentiated loop constraints determine the angle derivatives and the common $k$
cancels from the stationarity condition. Equivalently,
$\nabla U+J_g^\mathsf T\lambda=0$ for the loop constraints $g=0$.
Projecting onto the admissible tangent $J_g\theta'=0$ gives $U'(q)=0$.
The solver refines the minima found over a full revolution and chooses
the lowest-energy pose.

The default pose is **21.40° clockwise from the MotionGen start**
(**−84.62° absolute AC**). For the reference linear PETG stiffness it
stores **3.28 mJ**. The largest remaining twist is **18.27° at B–BQ**.
[Pose and energy plot](results/C2/least_loaded_pose.png)

### 4. Calculate wing motion and actuator torque

Virtual work maps each holding moment $M_i$ to the actuator [6].
With $q$ increasing in the direction of the sweep,

$$
T_\mathrm{flex}=\sum_i M_i(e_i)\,\frac{\mathrm d\delta_i}{\mathrm dq},
\qquad
T_\mathrm{drag}=-\mathbf F\cdot\frac{\mathrm d\mathbf r}{\mathrm dq}.
$$

The model uses the wing midpoint $\mathbf r=(T1+T2)/2$ and the orientation
of L11. It treats the wing as a **40 mm × 1.05 m** flat plate in still air.
The drag estimate uses $\tfrac12\rho C_D A$ and the projected area factor
$|\sin\alpha|$, where $\alpha$ is the angle between L11 and the relative
flow. The drag force opposes the midpoint velocity and scales with its
squared speed [7]. Give `--rpm` plus `--air-density` and `--drag-coefficient`
(or `--drag-factor`) for drag and combined quasi-static torque in N·m.
With no drag inputs, the plot shows a normalized drag term. Positive
actuator torque acts in the direction of the sweep.

For the reference PETG modulus [8], the linear shared-stiffness model gives
**−0.0308 to +0.0279 N·m**. The nonlinear model gives **−0.0319 to
+0.0290 N·m**. Both are elastic torque estimates on the nominal linkage
kinematics. The isolated beam model applies zero net external
force to the moving block, and the sweep excludes inertia and bearing
friction. [Torque comparison](results/C2_WING_SWEEP_PETG_NONLINEAR/input_torque.png)
· [CSV](results/C2_WING_SWEEP_PETG_NONLINEAR/wing_sweep.csv)

```powershell
python run_wing_sweep.py --nonlinear-cfp --cfp-length-mm 30 --cfp-width-mm 4.5 --cfp-thickness-mm 0.8 --cfp-e-gpa 2.1168 --out results/C2_WING_SWEEP_PETG_NONLINEAR
```

The simple pin-and-spring and drift-aware four-bar PRBMs are described in
[9]. For our separate crossed blades, $K=2EI/L$ gives
**0.02710 N·m/rad** near zero [10, Sec. 5.1]. At +43.5°, the dashboard
shows **20.57 N·mm** from this linear model and **22.37 N·mm** from the
nonlinear solver. Use the nonlinear moment curve for the full-stroke
torque calculation.

### 5. Relaxed fabrication and CAD handoff

Print four groups with their **own CFPs in the relaxed shape**. Rotate the
links within each group to the appropriate neutral angles, then align the
groups using their pin and clamp interfaces during assembly. Keeping the
ground rigid makes it a stable reference and lets the B flexures mount
separately. The course’s hybrid fabrication example [11, slide 16] and
the STAGE distinction between functional and stress-free fabrication
geometry [4, Sec. 2], [12, Sec. 2.2] inform this approach.

| Relaxed group | Contains | Assembly interface |
| --- | --- | --- |
| [Green](results/C2/step/concept2_green_relaxed.step) | Triangle C–D–T2, BD link, D and B–BD CFPs | C pin axis and T2 mounting direction |
| [Blue](results/C2/step/concept2_blue_relaxed.step) | Triangle P–Q–T1, AP and BQ links, P, Q and B–BQ CFPs | A pin axis and T1 mounting direction |
| [Red coupler](results/C2/step/concept2_coupler_relaxed.step) | L11 and the T1 and T2 CFPs | Free clamps mate to the triangle faces |
| [Rigid ground](results/C2/step/concept2_ground_relaxed.step) | Inclined AB (−26.88°) and a horizontal foot | Shared A pin axis and two B mounting interfaces |

The [fabrication drawing](results/C2/prbm_fabrication_layout.svg) shows
all four groups. The [frame JSON](results/C2/prbm_fabrication_frames.json)
gives the body-local and relaxed coordinates. The separate AC crank
connects A to C during assembly.

```powershell
python -m pip install cadquery-ocp-novtk==7.9.3.1.1
python run_concept2.py --require-closed --step --no-plot
```

Each AP214 file contains coloured link-centreline sheets, cylindrical
pivot datums along +Z, and planes through the pivot axes in the neutral
clamp directions. The datums use millimetres and stay fixed when blade
dimensions change. In Onshape, create explicit Mate connectors at the
axis centres (Z=0) and use the radial edge of the mount plane to orient
each connector [13]. Derive the parametric hinge at that connector and
build the solids from its clamp-face references. For a hinge in a
separate document, update the Derived reference to its new version
after changing the hinge geometry [14].

## Validation and design limits

- **Kinematics.** Tests cover a full revolution, position and velocity
  closure, relative angles and nominal output motion. The STEP files were
  re-imported in Open CASCADE to verify the surfaces, datums and colours.
- **Materials.** Dashboard presets come from manufacturer data for PLA
  [15], PLA Pro [16], PETG [8], conditioned CoPA [17], TPU95 [18]
  and Alleima 12R11 steel [19]. Filament tensile values and steel
  proof-stress values are screening references. The Steel 1500 preset
  also includes published reversed-bending fatigue references for 0.50
  and 0.75 mm strip [19, p. 4]. The check compares stress amplitude with
  the selected reference. Final sizing needs data for the actual material,
  print orientation and cycle life.
- **Prototype.** Choose the solid profiles, bores, fasteners and stacking
  heights in CAD. Check clearances through the full stroke, then compare
  strain, pivot drift and input torque with ANSYS or bench measurements.
  The ground support outline is a side-view concept. The model treats the
  hinge crossing as a nominal pivot. Full assembly loads, blade contact
  and out-of-plane deflections require a detailed assembly model.

## Acknowledgements

OpenAI Codex assisted with developing and checking Python code, deriving and
explaining the Concept 2 calculations, generating plots and schematics,
and editing this README and its references.

## References

[1] “Berekeningen voor Compliant,” unpublished Concept 2 calculation notes
and CAD dimensions, Nov. 17, 2026. Project file: `Berekeningen voor Compliant.pdf`.

[2] “2-fourbars,” unpublished Concept 2 MotionGen linkage model. Project file:
`2-fourbars.motiongen`.

[3] L. L. Howell, *Compliant Mechanisms*. New York, NY, USA: Wiley, 2001.
[Publisher record](https://www.wiley.com/en-us/-p-9780471384786).

[4] J. Rommers, V. van der Wijk, A. M. Aragón, and J. L. Herder, “The STAGE
method for simultaneous design of the stress and geometry of flexure
mechanisms,” *Precis. Eng.*, vol. 89, pp. 103–112, Aug. 2024,
doi: [10.1016/j.precisioneng.2024.05.021](https://doi.org/10.1016/j.precisioneng.2024.05.021).
[TU Delft repository copy](https://repository.tudelft.nl/record/uuid:a7900e35-5308-4ad1-960a-931a43ec39c1).

[5] C. Livermore, “Beam in pure bending,” lecture notes for 2.001, Mechanics
and Materials I, Massachusetts Institute of Technology, Cambridge, MA, USA,
Nov. 27, 2006. Accessed: Nov. 25, 2026. [Online]. Available:
[MIT OpenCourseWare](https://ocw.mit.edu/courses/2-001-mechanics-materials-i-fall-2006/1340543ced026dab5c8444a53fe92153_lec20.pdf).

[6] J. K. Vandiver, “An introduction to Lagrange equations,” lecture notes for
2.003SC, Engineering Dynamics, Massachusetts Institute of Technology,
Cambridge, MA, USA, Fall 2011. Accessed: Nov. 25, 2026. [Online]. Available:
[MIT OpenCourseWare](https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/d1984062464dcc3df68b5741fd37192b_MIT2_003SCF11_Lagrange.pdf).

[7] NASA Glenn Research Center, “Drag equation,” *Beginner’s Guide to
Aeronautics*. Accessed: Nov. 25, 2026. [Online]. Available:
[NASA](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/).

[8] Polymaker, “PolyLite PETG,” mechanical properties. Accessed: Nov. 25,
2026. [Online]. Available:
[Polymaker](https://polymaker.com/product/polylite-petg).

[9] B. D. Jensen and L. L. Howell, "The modeling of cross-axis flexural pivots,"
*Mech. Mach. Theory*, vol. 37, no. 5, pp. 461-476, May 2002,
doi: [10.1016/S0094-114X(02)00007-1](https://doi.org/10.1016/S0094-114X(02)00007-1).

[10] F. Cosandier, S. Henein, M. Richard, and L. Rubbert,
*The Art of Flexure Mechanism Design*. Lausanne, Switzerland: EPFL Press,
2017, ch. 5. [Open-access edition](https://www.epflpress.org/open_access_download/2229/1121).

[11] D. Farhadi Machekposhti, “Fabrication and experimental measurement of
compliant mechanisms,” Lecture 10, ME46115 Compliant Mechanisms, Delft
University of Technology, Delft, The Netherlands, 2026. Course slides:
`Lecture_Slides_CM_Week_8_Fabrication_Measurment.pdf`.

[12] J. Rommers, "Design strategies for large range flexure mechanisms,"
Ph.D. dissertation, Delft University of Technology, Delft, The Netherlands,
2022, doi: [10.4233/uuid:4eaf8bdb-53e2-49e5-9730-2550945b267f](https://doi.org/10.4233/uuid:4eaf8bdb-53e2-49e5-9730-2550945b267f).
[Full text](https://pure.tudelft.nl/ws/portalfiles/portal/122952921/JRommers_dissertation_7_6_2022.pdf).

[13] Onshape, “Mate connector,” *Onshape Help*. Accessed: Nov. 25, 2026.
[Online]. Available: [Onshape Help](https://cad.onshape.com/help/Content/mateconnector_a.htm).

[14] Onshape, “Derived,” *Onshape Help*. Accessed: Nov. 25, 2026.
[Online]. Available: [Onshape Help](https://cad.onshape.com/help/Content/PartStudio/derived.htm).

[15] Polymaker, “PolyLite PLA,” technical data sheet. Accessed: Nov. 25,
2026. [Online]. Available:
[Polymaker Wiki](https://wiki.polymaker.com/polymaker-products/more-about-our-products/documents/technical-data-sheets/pla/polylite-tm-pla).

[16] Polymaker, “PolyLite PLA Pro,” technical data sheet. Accessed: Nov. 25,
2026. [Online]. Available:
[Polymaker Wiki](https://wiki.polymaker.com/polymaker-products/more-about-our-products/documents/technical-data-sheets/pla/polylite-tm-pla-pro).

[17] Polymaker, “PolyMide CoPA,” technical data sheet. Accessed: Nov. 25,
2026. [Online]. Available:
[Polymaker Wiki](https://wiki.polymaker.com/polymaker-products/more-about-our-products/documents/technical-data-sheets/nylon/polymide-tm-copa).

[18] Polymaker, “PolyFlex TPU95,” technical data sheet. Accessed: Nov. 25,
2026. [Online]. Available:
[Polymaker Wiki](https://wiki.polymaker.com/polymaker-products/more-about-our-products/documents/technical-data-sheets/tpu/polyflex-tm-tpu95).

[19] Alleima, “Alleima 12R11,” strip steel datasheet, May 9, 2025. Accessed:
Nov. 25, 2026. [Online]. Available:
[Alleima](https://www.alleima.com/contentassets/56dacf45f3484ed7ad6a391aae380fd2/datasheet-alleima-12r11-en-v2025-05-09-0637-version-1.pdf/download).
