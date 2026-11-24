"""Run quasi-static torque sweep; plots and CSV for implemented concepts."""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from models import ConceptModel, DemoOneLinkModel
from config import get_config
from prbm import solve_static


def run(concept: str, demo: bool):
    config = get_config(concept)
    model = DemoOneLinkModel(config) if demo else ConceptModel(config)
    if not demo:
        print(f"[{concept}] Equations not implemented yet. Edit models.py to start actual analysis.")
        return
    torques = np.linspace(config.input_torque_min_nm, config.input_torque_max_nm, config.samples)
    q = np.array([config.initial_q_rad], dtype=float)
    rows = []
    for u in torques:
        q = solve_static(model, u, q)
        x, y = model.position(q)
        rows.append((u, q[0], x, y))
    out = Path("results") / concept.upper()
    out.mkdir(parents=True, exist_ok=True)
    with (out / "sweep.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["input_torque_Nm", "joint_angle_rad", "output_x_m", "output_y_m"])
        writer.writerows(rows)
    data = np.asarray(rows)
    fig, ax = plt.subplots()
    ax.plot(data[:,2]*1e3, data[:,3]*1e3)
    ax.set(xlabel="Output x [mm]", ylabel="Output y [mm]", title=f"{concept.upper()} — DEMO ONLY: output trajectory")
    ax.axis("equal")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(out / "trajectory_demo.png", dpi=160)
    plt.close(fig)
    print(f"DEMO ONLY results written to {out.resolve()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--concept", default="C2")
    parser.add_argument("--demo", action="store_true", help="Run artificial 1-link example")
    args = parser.parse_args()
    run(args.concept, args.demo)
