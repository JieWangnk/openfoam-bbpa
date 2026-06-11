#!/usr/bin/env python3
"""
Build a 3-panel phase animation (instantaneous | IPA | BBPA) for the
stenosis-pipe tutorial, on a longitudinal mid-plane slice of |U|.

For each phase point k = 0 .. nPhases-1 (phase theta = k/nPhases of the cycle):
  - instantaneous : |U| from the last solved cycle, at time t = (Ncyc-1)*T + k*T/n
  - IPA           : |U_IPA_phase{k}| from the final time directory (cumulative
                    pointwise phase average)
  - BBPA          : |U_PA| from the phase-aligned bin directory t = k*T/n
                    (cumulative bin-integrated phase average)

Usage:
    python3 make_phase_animation.py <caseDir> [outGif]
Requires: pyvista, imageio, numpy (and a headless GL stack -- start_xvfb is used).
"""
import os, sys, glob
import numpy as np
import pyvista as pv
import imageio.v2 as imageio

CASE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
OUT  = sys.argv[2] if len(sys.argv) > 2 else os.path.join(CASE, "stenosis_phase_animation.gif")
T, NPH, NCYC, FPS = 1.0, 50, 3, 12        # period, phase points, cycles, gif fps
VMAX = 2.6                                # |U| colour-scale top [m/s]

try:
    pv.start_xvfb()
except Exception:
    pass
pv.global_theme.font.size = 14

foam = os.path.join(CASE, "case.foam")
if not os.path.exists(foam):
    open(foam, "w").close()
reader = pv.OpenFOAMReader(foam)
reader.cell_to_point_creation = True
times = np.array(reader.time_values)
def nearest(t): return float(times[np.argmin(np.abs(times - t))])

def slice_mag(t, field):
    reader.set_active_time_value(nearest(t))
    mb = reader.read()
    blk = mb["internalMesh"] if "internalMesh" in mb.keys() else mb[0]
    if field not in blk.array_names:
        return None, None
    sl = blk.slice(normal="z")
    V = sl[field]
    sl["mag"] = np.linalg.norm(V, axis=1) if V.ndim == 2 else np.abs(V)
    return sl, sl["mag"]

frames = []
for k in range(NPH):
    th = k / NPH
    tphase = (NCYC - 1) * T + k * T / NPH      # last-cycle time at phase k
    panels = [
        ("instantaneous",       tphase,     "U"),
        (f"IPA  (N={NCYC-1})",   times[-1],  f"U_IPA_phase{k}"),   # cumulative, final dir
        (f"BBPA (N={NCYC-1})",   tphase,     "U_PA"),              # companion bin, last cycle
    ]
    p = pv.Plotter(off_screen=True, shape=(3, 1), window_size=(1280, 470))
    for j, (title, t, fld) in enumerate(panels):
        p.subplot(j, 0)
        sl, mag = slice_mag(t, fld)
        if sl is not None:
            p.add_mesh(sl, scalars="mag", cmap="turbo", clim=[0, VMAX],
                       show_scalar_bar=(j == 2),
                       scalar_bar_args=dict(title="|U| [m/s]", vertical=True))
        p.add_text(f"{title}  (phase {th:.2f})", font_size=11)
        p.view_xy()
        try:
            p.camera.tight(padding=0.06, view="xy", adjust_render_window=False)
        except Exception:
            p.camera.zoom(3.5)
    img = p.screenshot(return_img=True); p.close()
    frames.append(img)
    print(f"  frame {k+1}/{NPH} (phase {th:.2f})", flush=True)

imageio.mimsave(OUT, frames, fps=FPS, loop=0)
print(f"[wrote] {OUT}  ({len(frames)} frames)")
