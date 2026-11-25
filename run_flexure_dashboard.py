"""Explore an identical crossed-flexure design for the Concept 2 pivots."""

from functools import lru_cache

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, Button, TextBox
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches

from run_crossed_flexure import sweep

# ==========================================
# 1. MATERIAALDATABASE
# ==========================================
materialen = {
    "PLA": {"E": 3426.9e6, "stress_reference": 52.3e6,
            "naam": "PolyLite PLA", "basis": "XY tensile strength"},
    "PLA+": {"E": 2932.2e6, "stress_reference": 49.8e6,
             "naam": "PolyLite PLA Pro", "basis": "XY tensile strength"},
    "PETG": {"E": 2116.8e6, "stress_reference": 50.8e6,
             "naam": "PolyLite PETG", "basis": "XY tensile strength"},
    "PA": {"E": 724e6, "stress_reference": 34.3e6,
           "naam": "PolyMide CoPA (conditioned)", "basis": "XY tensile strength"},
    "TPU": {"E": 33.6e6, "stress_reference": 11.1e6,
            "naam": "PolyFlex TPU95", "basis": "XY stress at 100% strain"},
    "Steel": {"E": 185e9, "stress_reference": 1150e6,
              "naam": "Alleima 12R11 (cold rolled, Rm 1300 MPa)",
              "basis": "Rp0.2 - 0.2% offset proof stress"},
    "Steel 1500": {"E": 185e9, "stress_reference": 1375e6,
                   "naam": "Alleima 12R11 (cold rolled, Rm 1500 MPa)",
                   "basis": "Rp0.2 - 0.2% offset proof stress",
                   "fatigue_by_thickness": {0.5: 555, 0.75: 525}},
}

# Startwaarden
state = {
    "mat_key": "PETG",
    "L_mm": 30.0,
    "t_mm": 0.8,
    "b_mm": 4.5,
    "kruis_hoek": 90.0,
    "rot_deg": 0.0,
    "auto_play": True,
    "phase": 0.0
}

# ==========================================
# 2. BEREKENINGEN & KINEMATICA
# ==========================================
def bereken_alles():
    mat = materialen[state["mat_key"]]
    L = state["L_mm"] * 1e-3
    t = state["t_mm"] * 1e-3
    b = state["b_mm"] * 1e-3
    alpha = np.radians(state["kruis_hoek"] / 2.0)

    I = (b * t**3) / 12.0
    K_Nm_rad = 2.0 * mat["E"] * I / L
    K_Nmm_deg = K_Nm_rad * 1000.0 * (np.pi / 180.0)

    W_mm = state["L_mm"] * np.sin(alpha)
    H_mm = state["L_mm"] * np.cos(alpha)
    basis_hoek = 90.0 - (state["kruis_hoek"] / 2.0)
    fillet_R = max(1.0, round(state["t_mm"] * 1.5, 1))

    return {
        "mat": mat, "alpha": alpha, "K_Nm_rad": K_Nm_rad,
        "K_Nmm_deg": K_Nmm_deg,
        "W_mm": W_mm, "H_mm": H_mm, "basis_hoek": basis_hoek, "fillet_R": fillet_R
    }

MOTION_ANGLES = np.linspace(-45, 45, 91)


@lru_cache(maxsize=1)
def motion_poses(width_mm, height_mm, thickness_mm, blade_width_mm, modulus):
    return sweep(width_mm / 1000, height_mm / 1000, thickness_mm / 1000,
                 blade_width_mm / 1000, modulus, np.radians(MOTION_ANGLES))


# ==========================================
# 3. DASHBOARD OPZETTEN
# ==========================================
fig = plt.figure(figsize=(13.5, 7.8))
fig.canvas.manager.set_window_title("Crossed Flexure Pivot - Interactive Simulator & CAD")

# Twee hoofdplots: Links = Animatie, Rechts = CAD Schets
ax_sim = fig.add_axes([0.05, 0.34, 0.40, 0.58])
ax_cad = fig.add_axes([0.52, 0.34, 0.44, 0.58])

def teken_dashboard():
    res = bereken_alles()
    L_mm, t_mm, b_mm = state["L_mm"], state["t_mm"], state["b_mm"]
    W_mm, H_mm = res["W_mm"], res["H_mm"]
    th_deg = state["rot_deg"]

    # --- LINKER PLOT: SIMULATIE ---
    ax_sim.clear()
    lim = max(35.0, L_mm * 0.85)
    ax_sim.set_xlim(-lim, lim)
    ax_sim.set_ylim(-lim, lim)
    ax_sim.set_aspect('equal')
    ax_sim.grid(True, linestyle='--', alpha=0.5)
    ax_sim.set_title(f"Nonlinear 2D motion ({res['mat']['naam']})", fontweight='bold')
    ax_sim.set_xlabel("X [mm]")
    ax_sim.set_ylabel("Y [mm]")

    A1, B1_0 = np.array([-W_mm/2, -H_mm/2]), np.array([W_mm/2, H_mm/2])
    A2, B2_0 = np.array([W_mm/2, -H_mm/2]), np.array([-W_mm/2, H_mm/2])

    try:
        poses = motion_poses(W_mm, H_mm, t_mm, b_mm, res["mat"]["E"])
    except RuntimeError:
        ax_sim.text(.5, .5, "2D motion could not be solved for this geometry",
                    transform=ax_sim.transAxes, ha="center")
    else:
        index = np.clip(np.searchsorted(MOTION_ANGLES, th_deg) - 1, 0, len(poses) - 2)
        share = (th_deg - MOTION_ANGLES[index]) / (MOTION_ANGLES[index + 1] - MOTION_ANGLES[index])
        points = (1 - share) * poses[index][0] + share * poses[index + 1][0]
        pivot = (1 - share) * poses[index][3] + share * poses[index + 1][3]
        # The solver uses the fixed clamps as y = 0. Center the drawing on the neutral pivot.
        points[:, 1] -= H_mm / 2
        ax_sim.plot(*points[:, :, 0].T, "k-", lw=6, label="Fixed block")
        ax_sim.plot(*points[:, :, -1].T, "b-", lw=6, label="Moving block")
        ax_sim.plot(*points[0], "r-", lw=max(1.5, t_mm * 2.5), label="Flexure blades")
        ax_sim.plot(*points[1], "r-", lw=max(1.5, t_mm * 2.5))
        ax_sim.plot(*pivot, "go", markersize=8, label="Moving pivot reference")
        ax_sim.legend(loc="lower right", fontsize=8)

    moment_Nmm = res["K_Nmm_deg"] * abs(th_deg)
    info = (
        f"Rotation: {th_deg:+5.1f}\u00b0\n"
        f"Stiffness K: {res['K_Nmm_deg']:.4f} Nmm/° ({res['K_Nm_rad']:.4f} Nm/rad)\n"
        f"Linear torque: {moment_Nmm:5.2f} Nmm"
    )
    ax_sim.text(0.03, 0.96, info, transform=ax_sim.transAxes, va='top', fontsize=9,
                color="black",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.9))

    # --- RECHTER PLOT: 2D CAD SCHETS ---
    ax_cad.clear()
    blok_h = max(4.0, L_mm * 0.16)
    blok_w = W_mm + max(8.0, L_mm * 0.35)

    ax_cad.add_patch(patches.Rectangle((-blok_w/2, -H_mm/2 - blok_h), blok_w, blok_h,
                                       lw=1.5, ec='black', fc='#d9d9d9'))
    ax_cad.add_patch(patches.Rectangle((-blok_w/2, H_mm/2), blok_w, blok_h,
                                       lw=1.5, ec='black', fc='#e6f0ff'))

    ax_cad.plot(*np.stack((A1, B1_0)).T, color='#2b5c8f', lw=6, alpha=0.45)
    ax_cad.plot(*np.stack((A2, B2_0)).T, color='#d9534f', lw=6, alpha=0.45)
    ax_cad.plot(*np.stack((A1, B1_0)).T, 'k-.', lw=1)
    ax_cad.plot(*np.stack((A2, B2_0)).T, 'k-.', lw=1)
    ax_cad.plot(0, 0, 'go', markersize=6)

    # Maatlijnen
    y_dim = H_mm/2 + blok_h + L_mm * 0.10
    ax_cad.annotate("", xy=(-W_mm/2, y_dim), xytext=(W_mm/2, y_dim), arrowprops=dict(arrowstyle="<->"))
    ax_cad.plot([-W_mm/2, -W_mm/2], [H_mm/2, y_dim], 'k:', lw=0.8)
    ax_cad.plot([W_mm/2, W_mm/2], [H_mm/2, y_dim], 'k:', lw=0.8)
    ax_cad.text(0, y_dim + 0.6, f"W = {W_mm:.2f} mm", ha='center', va='bottom', fontweight='bold', fontsize=9)

    x_dim = blok_w/2 + L_mm * 0.10
    ax_cad.annotate("", xy=(x_dim, -H_mm/2), xytext=(x_dim, H_mm/2), arrowprops=dict(arrowstyle="<->"))
    ax_cad.plot([W_mm/2, x_dim], [-H_mm/2, -H_mm/2], 'k:', lw=0.8)
    ax_cad.plot([W_mm/2, x_dim], [H_mm/2, H_mm/2], 'k:', lw=0.8)
    ax_cad.text(x_dim + 0.8, 0, f"H = {H_mm:.2f} mm", ha='left', va='center', rotation=90, fontweight='bold', fontsize=9)

    ax_cad.text(0, H_mm*0.18, f"{state['kruis_hoek']:.0f}°", ha='center', fontweight='bold', fontsize=9)
    ax_cad.text(-W_mm/2 + 3.5, -H_mm/2 + 0.6, f"{res['basis_hoek']:.1f}°", ha='left', fontsize=8)

    cad_txt = (
        f"CAD DIMENSIONS: L = {L_mm:.1f} mm | t = {t_mm:.2f} mm | Width per blade = {b_mm:.1f} mm\n"
        f"Assembly depth = {2*b_mm:.1f} mm + clearance | Fillet R = {res['fillet_R']:.1f} mm"
    )
    ax_cad.text(0.02, 0.02, cad_txt, transform=ax_cad.transAxes, fontsize=8.5, va='bottom',
                bbox=dict(boxstyle="round", facecolor="#fffbe6", edgecolor="#ccaa00"))

    ax_cad.set_title("2D CAD drawing (centerlines)", fontweight='bold')
    ax_cad.set_aspect('equal')
    margin = max(14.0, L_mm * 0.45)
    ax_cad.set_xlim(-blok_w/2 - margin*0.5, blok_w/2 + margin)
    ax_cad.set_ylim(-H_mm/2 - blok_h - margin*0.8, H_mm/2 + blok_h + margin*0.8)
    ax_cad.grid(True, linestyle=':', alpha=0.4)
    ax_cad.set_xlabel("X [mm]")

# ==========================================
# 4. INTERACTIEVE WIDGETS (SLIDERS & KNOPPEN)
# ==========================================
# Materiaal RadioButtons linksonder
ax_radio = fig.add_axes([0.05, 0.04, 0.14, 0.23], facecolor='#f4f4f4')
ax_radio.set_title("Material", fontsize=10, fontweight='bold')
radio = RadioButtons(ax_radio, list(materialen.keys()), active=2)

# Sliders in het middenonder
ax_L    = fig.add_axes([0.29, 0.22, 0.45, 0.03])
ax_t    = fig.add_axes([0.29, 0.17, 0.45, 0.03])
ax_b    = fig.add_axes([0.29, 0.12, 0.45, 0.03])
ax_hoek = fig.add_axes([0.29, 0.07, 0.45, 0.03])
ax_rot  = fig.add_axes([0.29, 0.02, 0.45, 0.03])

s_L    = Slider(ax_L,    "Length L [mm] ", 10.0, 60.0, valinit=state["L_mm"], valstep=0.5)
s_t    = Slider(ax_t,    "Thickness t [mm]  ", 0.3,  2.5,  valinit=state["t_mm"], valstep=0.05)
s_b    = Slider(ax_b,    "Blade width [mm]", 2.0,  30.0, valinit=state["b_mm"], valstep=0.5)
s_hoek = Slider(ax_hoek, "Crossing angle [°] ", 40.0, 140.0, valinit=state["kruis_hoek"], valstep=1.0)
s_rot  = Slider(ax_rot,  "Manual angle θ [°]", -45.0, 45.0, valinit=0.0, valstep=0.5)

# Knoppen rechtsonder
ax_play = fig.add_axes([0.80, 0.23, 0.15, 0.05])
ax_save = fig.add_axes([0.80, 0.09, 0.15, 0.05])
ax_fit = fig.add_axes([0.80, 0.02, 0.15, 0.05])
btn_play = Button(ax_play, "Pause / Manual", color="#e0efff", hovercolor="#cce0ff")
btn_save = Button(ax_save, "Save CAD PNG", color="#e6ffe6", hovercolor="#ccffcc")

btn_fit = Button(ax_fit, "Fit stiffness")
ax_check = fig.add_axes([0.80, 0.16, 0.15, 0.05])
btn_check = Button(ax_check, "Check configuration")
ax_fatigue = fig.add_axes([0.43, 0.946, 0.18, 0.032])
fatigue_input = TextBox(ax_fatigue, "Fatigue amplitude limit [MPa] ", initial="")
fig.text(0.62, 0.96, "Leave blank to use matching preset data", fontsize=9, va="center")

# Event handlers
def on_slider_change(val):
    state["L_mm"] = s_L.val
    state["t_mm"] = s_t.val
    state["b_mm"] = s_b.val
    state["kruis_hoek"] = s_hoek.val
    if not state["auto_play"]:
        state["rot_deg"] = s_rot.val
    teken_dashboard()
    fig.canvas.draw_idle()

def on_rot_manual(val):
    state["auto_play"] = False
    btn_play.label.set_text("Start animation")
    state["rot_deg"] = s_rot.val
    teken_dashboard()
    fig.canvas.draw_idle()

def on_mat_change(label):
    state["mat_key"] = label
    teken_dashboard()
    fig.canvas.draw_idle()

def toggle_play(event):
    state["auto_play"] = not state["auto_play"]
    btn_play.label.set_text("Pause / Manual" if state["auto_play"] else "Start animation")

def save_cad(event):
    extent = ax_cad.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    filename = f"flexure_cad_{state['mat_key']}_L{state['L_mm']:.0f}_t{state['t_mm']:.2f}.png"
    fig.savefig(filename, bbox_inches=extent.expanded(1.18, 1.18), dpi=200)
    print(f"CAD drawing saved to: {filename}")

for s in [s_L, s_t, s_b, s_hoek]:
    s.on_changed(on_slider_change)
s_rot.on_changed(on_rot_manual)
radio.on_clicked(on_mat_change)
btn_play.on_clicked(toggle_play)
btn_save.on_clicked(save_cad)

def analyze_configuration(event, *, show_fit):
    res = bereken_alles()
    fatigue_limit = None
    fatigue_label = "No matching fatigue data"
    entered = fatigue_input.text.strip()
    if entered:
        try:
            fatigue_limit = float(entered)
            if not np.isfinite(fatigue_limit) or fatigue_limit <= 0:
                raise ValueError
        except ValueError:
            message = "Enter a positive fatigue amplitude limit in MPa, or leave it blank."
            print(message)
            plt.figure(figsize=(8, 2)).text(.05, .5, message, wrap=True)
            plt.gcf().show()
            return
        fatigue_label = "user-supplied limit"
    else:
        thickness = round(state["t_mm"], 8)
        fatigue_limit = res["mat"].get("fatigue_by_thickness", {}).get(thickness)
        if fatigue_limit is not None:
            fatigue_label = "2 million cycles, reversed bending, median survival"
    limit = abs(s_rot.val)
    if limit == 0:
        limit = 25.0 if show_fit else 45.0
    angles = np.radians(np.linspace(-limit, limit, 51))
    try:
        poses = sweep(res["W_mm"] / 1000, res["H_mm"] / 1000,
                      state["t_mm"] / 1000, state["b_mm"] / 1000,
                      res["mat"]["E"], angles)
    except RuntimeError as error:
        message = f"Check incomplete: {error}"
        print(message)
        plt.figure(figsize=(8, 2)).text(.05, .5, message, wrap=True)
        plt.gcf().show()
        return
    moments = np.array([pose[2] for pose in poses])
    fitted_k = np.dot(angles, moments) / np.dot(angles, angles)
    fig_fit = plt.figure(figsize=(8, 7 if show_fit else 5), layout="constrained")
    ax_stress = fig_fit.add_subplot(2, 1, 2) if show_fit else fig_fit.add_subplot()
    if show_fit:
        ax_fit_plot = fig_fit.add_subplot(2, 1, 1)
        ax_fit_plot.plot(np.degrees(angles), moments, label="Nonlinear 2D model")
        ax_fit_plot.plot(np.degrees(angles), fitted_k * angles, "--",
                         label=f"Linear fit: {fitted_k:.6g} N m/rad")
        ax_fit_plot.set(xlabel="Rotation (\u00b0)", ylabel="Reaction moment (N m)",
                        title=f"Shared stiffness over \u00b1{limit:g} degrees")
        ax_fit_plot.grid(alpha=.25)
        ax_fit_plot.legend()
    surface_stresses = np.array([pose[1] for pose in poses])
    stress = np.abs(surface_stresses).max(axis=(1, 2, 3))
    local_low = surface_stresses.min(axis=0)
    local_high = surface_stresses.max(axis=0)
    amplitude = ((local_high - local_low) / 2).max()
    largest_mean = np.abs((local_high + local_low) / 2).max()
    fatigue_summary = f"Fatigue: peak alternating stress {amplitude:.1f} MPa"
    if fatigue_limit is None:
        fatigue_summary += " | No matching fatigue data"
    else:
        comparison = "below" if amplitude < fatigue_limit else "at or above"
        fatigue_summary += f" | {comparison} {fatigue_limit:g} MPa ({fatigue_label})"
    mean_note = (f"Max absolute local mean stress: {largest_mean:.1f} MPa. "
                 "Amplitude comparison only; mean stress and clearance are unchecked.")
    threshold = res["mat"]["stress_reference"] / 1e6
    ax_stress.plot(np.degrees(angles), stress, label="Peak blade stress")
    ax_stress.axhline(threshold, color="tab:red", linestyle="--",
                      label=f"{res['mat']['basis']} ({threshold:g} MPa)")
    crossed = np.flatnonzero(stress >= threshold)
    if crossed.size:
        middle = len(angles) // 2
        envelope = np.maximum(stress[middle:], stress[middle::-1])
        i = np.flatnonzero(envelope >= threshold)[0]
        onset = np.interp(threshold, envelope[i - 1:i + 1],
                          np.degrees(angles[middle + i - 1:middle + i + 1]))
        status = f"Reference stress reached near \u00b1{onset:.1f}\u00b0 in the 2D model"
    else:
        status = f"Reference stress not reached through \u00b1{limit:g}\u00b0 in the 2D model"
    ax_stress.set(xlabel="Rotation (\u00b0)", ylabel="Peak stress (MPa)", title=status)
    ax_stress.grid(alpha=.25)
    ax_stress.legend()
    angle_range = f"Tested range: \u00b1{limit:g}\u00b0"
    result = "AT OR ABOVE REFERENCE" if crossed.size else "BELOW REFERENCE"
    summary = f"{result} | {angle_range} | {stress.max():.1f} / {threshold:g} MPa"
    fig_fit.suptitle(f"{summary}\n{res['mat']['basis']} = {threshold:g} MPa")
    fig_fit.supxlabel(f"{fatigue_summary}\n{mean_note}", fontsize=9)
    ax_stress.set_xlabel("Rotation (degrees)")
    print(f"{res['mat']['naam']} | E = {res['mat']['E'] / 1e9:g} GPa | {res['mat']['basis']}")
    print(summary)
    print(status)
    print(fatigue_summary)
    print(mean_note)
    if show_fit:
        print(f"python run_wing_sweep.py --common-k {fitted_k:.9g} --out results/C2_SHARED_CFP")
    fig_fit.show()
    return fitted_k


btn_fit.on_clicked(lambda event: analyze_configuration(event, show_fit=True))
btn_check.on_clicked(lambda event: analyze_configuration(event, show_fit=False))


def animate(frame):
    if state["auto_play"]:
        state["phase"] += 0.12
        state["rot_deg"] = 25 * np.sin(state["phase"])
        teken_dashboard()


teken_dashboard()
ani = FuncAnimation(fig, animate, interval=80, cache_frame_data=False)
if __name__ == "__main__":
    plt.show()
