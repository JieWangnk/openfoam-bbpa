"""
Render a 3-panel wall-shear-stress comparison on the stenosis/aorta wall.

Each panel colours the wall_aorta (or equivalent) patch by |WSS|*rho in Pa,
using a Black-Body Radiation LUT clamped to [0, 30] Pa for visual
consistency across panels.

Panels (default):
    A: instantaneous WSS                         -> wallShearStress
    B: running-window mean,  w = 0.001 s (I=500) -> wallShearStressMean_w_500
    C: running-window mean,  w = 0.01  s (I=50)  -> wallShearStressMean_w_50

Purpose: visually demonstrate the "filter-width" argument in the paper's
Table II / \S4.1 - widening the averaging window progressively low-passes
the near-wall WSS content. Directly supports the -28.2% WSS spread claim.

Usage (pvbatch):

    pvbatch screenshot_WSS_BPM120.py \\
        --case /path/to/reconstructed/case \\
        --time 16.07088 \\
        --fields wallShearStress wallShearStressMean_w_500 wallShearStressMean_w_50 \\
        --labels "Instantaneous" "w=0.001s (I=500)" "w=0.01s (I=50)" \\
        --out /path/to/figures/fig_wss_filter_width.png

For a local dry-run on the BBPA tutorial (only wallShearStress available):

    pvbatch screenshot_WSS_BPM120.py \\
        --case /home/mchi4jw4/OpenFOAM/mchi4jw4-12/tutorials/stenosisPipeBBPA \\
        --time 1.48 --times 0.5 1.0 1.48 \\
        --fields wallShearStress wallShearStress wallShearStress \\
        --labels "t=0.5s" "t=1.0s" "t=1.48s" \\
        --out /tmp/wss_test.png
"""

import argparse
import os
import sys
from pathlib import Path

from paraview.simple import *  # noqa: F401, F403  (ParaView globals)


# -----------------------------------------------------------------------------
#  Defaults - tuned for the aortic coarctation geometry; override via CLI
# -----------------------------------------------------------------------------
DEFAULT_RHO = 1060.0               # kg/m^3 blood density
DEFAULT_WSS_RANGE = (0.0, 30.0)    # Pa, colour-map clamp
DEFAULT_PATCH = "wall_aorta"       # boundary patch name

# Camera aligned with previous screenshot trace; override via --camera
DEFAULT_CAMERA = {
    "position": [0.18820, -0.13075, 0.00777],
    "focal_point": [-0.00992, -0.00887, -0.01984],
    "view_up": [-0.12383, 0.02344, 0.99203],
    "parallel_scale": 0.06063,
}

DEFAULT_IMAGE_SIZE = (1600, 900)


# -----------------------------------------------------------------------------
#  Rendering helpers
# -----------------------------------------------------------------------------
def open_case(case_dir: str):
    """Return the ParaView OpenFOAMReader source for case_dir."""
    foam_stub = os.path.join(case_dir, "case.foam")
    if not os.path.exists(foam_stub):
        Path(foam_stub).touch()
    reader = OpenFOAMReader(FileName=foam_stub)
    # Settings are best-effort: newer ParaView versions rename these props;
    # we silently skip unknowns so the script stays portable.
    for attr, val in [
        ("SkipZeroTime", 0),
        ("CaseType", "Reconstructed Case"),
    ]:
        try:
            setattr(reader, attr, val)
        except (AttributeError, Exception):
            pass
    reader.UpdatePipeline()
    return reader


def open_vtk(vtk_path: str):
    """Open a single legacy VTK file (e.g. from foamToVTK -patches ...)."""
    reader = LegacyVTKReader(FileNames=[vtk_path])
    reader.UpdatePipeline()
    return reader


def restrict_to_patch(reader, patch_name: str):
    """Return a source that contains only the named wall patch block."""
    # Default Mesh Parts in ParaView's OpenFOAMReader:
    #   internalMesh, patch/<name>, patch/<other> ...
    # We select only the requested boundary patch.
    try:
        reader.MeshRegions = ["patch/" + patch_name]
    except Exception:
        # Older ParaView may expose as just the patch name
        reader.MeshRegions = [patch_name]
    reader.UpdatePipeline()
    return reader


def pa_calculator(source, field: str, rho: float):
    """Return a Calculator that outputs |field|*rho in Pa under name 'WSS_Pa'.

    Inserts CellDatatoPointData so the vector array is visible to the
    function parser on both cell- and point-centred inputs.
    """
    c2p = CellDatatoPointData(Input=source)
    c2p.ProcessAllArrays = 1
    c2p.UpdatePipeline()
    calc = Calculator(Input=c2p)
    calc.AttributeType = "Point Data"
    calc.ResultArrayName = "WSS_Pa"
    calc.Function = f"{rho}*mag({field})"
    calc.UpdatePipeline()
    return calc


def make_lut(view, wss_range):
    """Configure a shared Black-Body Radiation LUT for WSS_Pa."""
    lut = GetColorTransferFunction("WSS_Pa")
    lut.ApplyPreset("Black-Body Radiation", True)
    lut.RescaleTransferFunction(wss_range[0], wss_range[1])
    pwf = GetOpacityTransferFunction("WSS_Pa")
    pwf.RescaleTransferFunction(wss_range[0], wss_range[1])
    return lut


def place_camera(view, cam=None):
    view.OrientationAxesVisibility = 0
    if cam is None:
        view.ResetCamera(False)
        return
    view.CameraPosition = cam["position"]
    view.CameraFocalPoint = cam["focal_point"]
    view.CameraViewUp = cam["view_up"]
    view.CameraParallelScale = cam["parallel_scale"]


def style_colorbar(view, lut, title):
    bar = GetScalarBar(lut, view)
    bar.AutoOrient = 0
    bar.Orientation = "Horizontal"
    bar.Title = title
    bar.ComponentTitle = ""
    bar.TitleBold = 1
    bar.TitleFontSize = 18
    bar.LabelBold = 1
    bar.LabelFontSize = 16
    bar.WindowLocation = "Any Location"
    bar.Position = [0.22, 0.06]
    bar.ScalarBarLength = 0.55
    bar.AddRangeLabels = 0


def render_panel(
    case_dir: str,
    time_value: float,
    field: str,
    rho: float,
    patch: str,
    wss_range,
    out_png: str,
    label: str,
    image_size,
    camera,
    vtk_path: str = None,
):
    """Render a single panel to out_png.

    Opens EITHER an OpenFOAM case (case_dir) with optional time selection
    and patch restriction, OR a single legacy VTK file (vtk_path).
    The VTK path wins if both are supplied.
    """
    if vtk_path:
        reader = open_vtk(vtk_path)
    else:
        reader = open_case(case_dir)
        restrict_to_patch(reader, patch)

    view = CreateView("RenderView")
    view.ViewSize = list(image_size)
    view.Background = [1, 1, 1]
    view.UseColorPaletteForBackground = 0

    if not vtk_path:
        scene = GetAnimationScene()
        scene.UpdateAnimationUsingDataTimeSteps()
        scene.AnimationTime = time_value

    calc = pa_calculator(reader, field, rho)
    disp = Show(calc, view)
    disp.Representation = "Surface"
    ColorBy(disp, ("POINTS", "WSS_Pa"))

    lut = make_lut(view, wss_range)
    disp.LookupTable = lut
    disp.SetScalarBarVisibility(view, True)

    place_camera(view, camera)
    style_colorbar(view, lut, f"{label} — WSS (Pa)")

    view.Update()
    SaveScreenshot(
        out_png,
        view,
        ImageResolution=list(image_size),
        FontScaling="Do not scale fonts",
        TransparentBackground=0,
    )
    Delete(disp)
    Delete(calc)
    Delete(reader)
    Delete(view)
    print(f"  wrote {out_png}")


# -----------------------------------------------------------------------------
#  CLI
# -----------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", type=str, default=None,
                   help="OpenFOAM case directory (reconstructed). Alternative to --vtk.")
    p.add_argument("--vtk", type=str, default=None,
                   help="Single legacy VTK file (e.g. from foamToVTK -patches); alternative to --case.")
    p.add_argument("--time", type=float, default=None,
                   help="Single time value used for all panels (if --times not given)")
    p.add_argument("--times", type=float, nargs="+", default=None,
                   help="Per-panel time values (overrides --time)")
    p.add_argument("--fields", nargs="+", required=True,
                   help="Per-panel OpenFOAM field name (vector; magnitude is taken)")
    p.add_argument("--labels", nargs="+", default=None,
                   help="Per-panel caption label; defaults to field names")
    p.add_argument("--patch", type=str, default=DEFAULT_PATCH,
                   help=f"Wall patch name (default: {DEFAULT_PATCH})")
    p.add_argument("--rho", type=float, default=DEFAULT_RHO,
                   help=f"Blood density kg/m^3 (default: {DEFAULT_RHO})")
    p.add_argument("--wss-min", type=float, default=DEFAULT_WSS_RANGE[0])
    p.add_argument("--wss-max", type=float, default=DEFAULT_WSS_RANGE[1])
    p.add_argument("--out", required=True, type=str,
                   help="Output PNG path (per-panel images get A/B/C suffix)")
    p.add_argument("--size", type=int, nargs=2, default=list(DEFAULT_IMAGE_SIZE))
    p.add_argument("--use-default-camera", action="store_true",
                   help="Use the hard-coded aorta camera (default: auto-fit to data)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    if not args.case and not args.vtk:
        sys.exit("Must supply either --case or --vtk")

    n = len(args.fields)
    times = args.times if args.times else [args.time if args.time else 0.0] * n
    if len(times) != n:
        sys.exit(f"--times and --fields length mismatch ({len(times)} vs {n})")
    labels = args.labels if args.labels else args.fields

    out_base = Path(args.out)
    out_base.parent.mkdir(parents=True, exist_ok=True)

    for i, (t, field, label) in enumerate(zip(times, args.fields, labels)):
        suffix = "ABC"[i] if n <= 3 else str(i)
        panel_png = out_base.parent / f"{out_base.stem}_{suffix}{out_base.suffix}"
        print(f"[panel {suffix}] t={t} field={field} label={label!r}")
        render_panel(
            case_dir=args.case,
            time_value=t,
            field=field,
            rho=args.rho,
            patch=args.patch,
            wss_range=(args.wss_min, args.wss_max),
            out_png=str(panel_png),
            label=label,
            image_size=tuple(args.size),
            camera=DEFAULT_CAMERA if args.use_default_camera else None,
            vtk_path=args.vtk,
        )
    print(f"\nDone. {n} panels written under {out_base.parent}/")
    print("Compose into a single figure with ImageMagick, for example:")
    print(f"  montage {out_base.stem}_A{out_base.suffix} "
          f"{out_base.stem}_B{out_base.suffix} {out_base.stem}_C{out_base.suffix} "
          f"-geometry +10+10 -tile 3x1 {out_base}")


if __name__ == "__main__":
    main()
