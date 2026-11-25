"""Nonlinear planar crossed-flexure pivot, based on Rommers et al. (2024)."""

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_bvp


def sweep(lx, ly, thickness, width, modulus, angles, preload=None, tolerance=1e-7):
    """Solve two extensible Euler–Bernoulli rods with a free-translating plate.

    Coordinates are scaled by blade length L, forces by EI/L² and moments
    by EI/L. States per rod are x, y, tangent angle and bending moment.
    """
    length = np.hypot(lx, ly)
    half_span, height = lx / (2 * length), ly / length
    starts = np.array([-half_span, half_span])
    arms = -starts
    beta = np.arctan2(height, 2 * arms)
    slenderness = thickness**2 / (12 * length**2)
    stress_scale = modulus * thickness / (2 * length)
    ei = modulus * width * thickness**3 / 12
    mesh = np.linspace(0, 1, 41)
    sample = np.linspace(0, 1, 401)
    neutral = np.zeros((2, 2)) if preload is None else np.array([preload, -preload])
    slope = neutral[:, 1] - neutral[:, 0]
    neutral_fx = slope / np.sin(beta)
    natural_stretch = 1 - slenderness * neutral_fx * np.cos(beta)

    def ode(s, state, params):
        rods = state.reshape(2, 4, -1)
        fx, fy = params[:4].reshape(2, 2).T[:, :, None]
        c, sn = np.cos(rods[:, 2]), np.sin(rods[:, 2])
        stretch = natural_stretch[:, None] + slenderness * (fx * c + fy * sn)
        curvature = rods[:, 3] - neutral[:, :1] - slope[:, None] * s
        return np.stack((stretch * c, stretch * sn, curvature,
                         stretch * (sn * fx - c * fy)), axis=1).reshape(8, -1)

    def solve(angle, state, params):
        arm_vectors = arms[:, None] * np.array([np.cos(angle), np.sin(angle)])

        def boundary(left, right, p):
            left = left.reshape(2, 4)
            right = right.reshape(2, 4)
            force = p[:4].reshape(2, 2)
            return np.r_[left[:, 0] - starts, left[:, 1], left[:, 2] - beta,
                         (right[:, :2] - p[4:] - arm_vectors).ravel(),
                         right[:, 2] - beta - angle, force.sum(axis=0)]

        solution = solve_bvp(ode, boundary, mesh, state, p=params, tol=tolerance)
        if not solution.success:
            raise RuntimeError(f"Equilibrium failed at {np.degrees(angle):.2f} deg: {solution.message}")
        rods = solution.sol(sample).reshape(2, 4, -1)
        forces = solution.p[:4].reshape(2, 2)
        axial = (forces[:, :1] * np.cos(rods[:, 2]) +
                 forces[:, 1:] * np.sin(rods[:, 2]))
        axial_stress = modulus * slenderness * axial
        bending_stress = stress_scale * rods[:, 3]
        surface_stress = np.stack((axial_stress + bending_stress,
                                   axial_stress - bending_stress), axis=1)
        torque = (rods[:, 3, -1] + arm_vectors[:, 0] * forces[:, 1] -
                  arm_vectors[:, 1] * forces[:, 0]).sum() * ei / length
        center = solution.p[4:]
        pivot_drift = center - [0, height] - height / 2 * np.array(
            [-np.sin(angle), np.cos(angle) - 1])
        return solution, (rods[:, :2] * length * 1000, surface_stress / 1e6,
                          torque, pivot_drift * length * 1000, rods[:, 3])

    initial = np.stack((starts[:, None] + 2 * arms[:, None] * mesh,
                        np.broadcast_to(height * mesh, (2, mesh.size)),
                        np.broadcast_to(beta[:, None], (2, mesh.size)),
                        neutral[:, :1] + slope[:, None] * mesh), axis=1).reshape(8, -1)
    parameters = np.ravel(np.column_stack((neutral_fx, np.zeros(2))))
    parameters = np.r_[parameters, 0, height]
    poses = {}
    for direction in (-1, 1):
        state, params = initial, parameters
        for angle in angles[len(angles) // 2::direction]:
            solution, poses[angle] = solve(angle, state, params)
            state, params = solution.sol(mesh), solution.p
    return [poses[angle] for angle in angles]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize

    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in (("lx-mm", 100), ("ly-mm", 100), ("thickness-mm", .5),
                          ("width-mm", 15), ("e-gpa", 114), ("angle-deg", 45)):
        parser.add_argument("--" + name, type=float, default=default)
    parser.add_argument("--steps", type=int, default=91)
    parser.add_argument("--out", type=Path, default=Path("results/CFP"))
    args = parser.parse_args()
    values = (args.lx_mm, args.ly_mm, args.thickness_mm, args.width_mm, args.e_gpa)
    if (not all(np.isfinite(v) and v > 0 for v in values)
            or not 0 < args.angle_deg <= 45 or args.steps < 3 or args.steps % 2 == 0):
        parser.error("Use positive dimensions and E, 0 < angle <= 45, and an odd step count >= 3.")

    angles = np.radians(np.linspace(-args.angle_deg, args.angle_deg, args.steps))
    geometry = (*(v / 1000 for v in values[:4]), args.e_gpa * 1e9)
    straight = sweep(*geometry, angles)
    clamp_moments = np.array([pose[4][0, [0, -1]] for pose in straight])
    # STAGE: center the bending-moment extremes at each clamp, then realize
    # that neutral moment field with a stress-free curvature and axial strain.
    preload = -(clamp_moments.min(axis=0) + clamp_moments.max(axis=0)) / 2
    stage = sweep(*geometry, angles, preload)
    designs = {"Straight blades": straight, "STAGE-inspired": stage}
    args.out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
    sections = [("Peak surface stress (MPa)", axes[0, 0]),
                ("Reaction moment (N·m)", axes[0, 1]),
                ("Virtual pivot drift (mm)", axes[0, 2])]
    peak = max(np.abs(pose[1]).max() for poses in designs.values() for pose in poses)
    norm = Normalize(0, peak)
    with (args.out / "crossed_flexure.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["design", "angle_deg", "peak_stress_mpa", "peak_abs_strain_pct",
                         "torque_nm", "pivot_dx_mm", "pivot_dy_mm"])
        for column, (label, poses) in enumerate(designs.items()):
            stresses = [np.abs(pose[1]).max() for pose in poses]
            torques = [pose[2] for pose in poses]
            drifts = np.array([pose[3] for pose in poses])
            for angle, stress, torque, drift in zip(np.degrees(angles), stresses, torques, drifts):
                writer.writerow([label, angle, stress, 100 * stress / (args.e_gpa * 1000),
                                 torque, *drift])
            for series, (_, axis) in zip((stresses, torques, np.linalg.norm(drifts, axis=1)), sections):
                axis.plot(np.degrees(angles), series, label=label)
            axis = axes[1, column]
            for points, surface, *_ in (poses[0], poses[-1]):
                for blade, stress in zip(points, surface):
                    xy = blade.T
                    segments = np.stack((xy[:-1], xy[1:]), axis=1)
                    colors = np.abs(stress).max(axis=0)
                    axis.add_collection(LineCollection(segments, cmap="viridis", norm=norm,
                                                        array=(colors[:-1]+colors[1:])/2, linewidth=2))
                axis.plot(*points[:, :, -1].T, color=".3", linewidth=2)
            axis.plot(*points[:, :, 0].T, color=".3", linewidth=3)
            axis.autoscale()
            axis.set(aspect="equal", title=label + " at both limits", xlabel="x (mm)", ylabel="y (mm)")
            print(f"{label}: {max(stresses):.2f} MPa peak, "
                  f"{max(np.abs(torques)):.3f} N m peak absolute moment, "
                  f"{np.linalg.norm(drifts, axis=1).max():.3f} mm peak drift")
    for ylabel, axis in sections:
        axis.set(xlabel="Angle (deg)", ylabel=ylabel)
        axis.grid(alpha=.25)
        axis.legend()
    axis = axes[1, 2]
    length = np.hypot(*geometry[:2])
    s = np.linspace(0, 1, 401)
    neutral_shape = straight[len(angles) // 2][0]
    for blade, moments in zip(neutral_shape, (preload, -preload)):
        delta = blade[:, -1] - blade[:, 0]
        beta = np.arctan2(delta[1], delta[0])
        slope = moments[1] - moments[0]
        tangent = beta - moments[0] * s - slope * s**2 / 2
        stretch = 1 - geometry[2]**2 / (12 * length**2) * slope / np.tan(beta)
        released = blade[:, :1] + length * 1000 * stretch * cumulative_trapezoid(
            [np.cos(tangent), np.sin(tangent)], s, axis=1, initial=0)
        axis.plot(*blade, "--", color=".55")
        axis.plot(*released, color="tab:orange")
    axis.set(aspect="equal", title="Released blades (orange)", xlabel="x (mm)", ylabel="y (mm)")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"), ax=axes[1, :2],
                 label="Peak absolute surface stress (MPa)", shrink=.8)
    fig.savefig(args.out / "crossed_flexure.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
