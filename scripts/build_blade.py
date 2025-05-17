"""
build_blade.py  ─ Fusion 360 add-in
-----------------------------------
Builds a wind-turbine blade from a JSON specification, then exports a local .f3d
file and quits Fusion (so the calling process can continue).

Environment variables (required / optional):
    FUSION_BLADE_JSON   absolute path to blade_spec.json          (required)
    FUSION_AIRFOIL_DIR  folder containing all .dat airfoil files  (required)
    FUSION_BLADE_SAVE   absolute output path  (defaults to <json>.f3d)

Notes
-----
* Every blade section **must** have a valid "coordinate_file" (e.g. "2412.dat").
* All .dat files live in FUSION_AIRFOIL_DIR.
* Units in JSON are metres; script converts to the document’s native units.
"""

import adsk.core, adsk.fusion, traceback, os, json, math
from pathlib import Path
from aeropy.airfoils.shapes import naca, cst




# # ------------------------------------------- #
# # ------------- INPUT JSON NAME ------------- #
# # ------------------------------------------- #
# JSON_NAME = "blade a.json"
# # ------------------------------------------- #
# # ------------------------------------------- #
# # ------------------------------------------- #
#
# # ── 0  Input / output paths ─────────────────────────────────────────────
# DATA_DIR = Path(__file__).resolve().parent.parent / "naca_data"
#
# JSON_PATH = DATA_DIR / "blade_profiles" / JSON_NAME
# AIRFOIL_DIR = DATA_DIR / "airfoil_profiles"
#
#
# if not JSON_PATH or not os.path.isfile(JSON_PATH):
#     raise RuntimeError("Set env var FUSION_BLADE_JSON to a valid file.")
# if not AIRFOIL_DIR or not os.path.isdir(AIRFOIL_DIR):
#     raise RuntimeError("Set env var FUSION_AIRFOIL_DIR to a valid folder.")
#
# SAVE_PATH = DATA_DIR / "blade_models" / (JSON_NAME.replace('.json', '.f3d'))
# # ─────────────────────────────────────────────────────────────────────────


JSON_PATH   = os.getenv("FUSION_BLADE_JSON")    # required
AIRFOIL_DIR = os.getenv("FUSION_AIRFOIL_DIR")   # required
if not JSON_PATH or not os.path.isfile(JSON_PATH):
    raise RuntimeError("Set env var FUSION_BLADE_JSON to a valid file.")
if not AIRFOIL_DIR or not os.path.isdir(AIRFOIL_DIR):
    raise RuntimeError("Set env var FUSION_AIRFOIL_DIR to a valid folder.")

SAVE_PATH = os.getenv("FUSION_BLADE_SAVE") or os.path.splitext(JSON_PATH)[0] + ".f3d"

# ── helpers ──────────────────────────────────────────────────────────────
def read_dat(path):
    pts = []
    with open(path, "r") as fp:
        for ln in fp:
            try:
                x, y = [float(v) for v in ln.strip().split()[:2]]
                pts.append((x, y))
            except ValueError:
                continue
    if pts and pts[0] != pts[-1]:
        pts.append(pts[0])
    return pts


def get_naca_coordinates(naca_code):
    foil = naca(naca_code, 200)
    return foil


def to_object_collection(points, chord_len, twist_rad):
    oc = adsk.core.ObjectCollection.create()
    cos_t, sin_t = math.cos(twist_rad), math.sin(twist_rad)
    for x, y in points:
        # shift TE (x=1) to origin, scale to chord
        px = (x - 1.0) * chord_len
        py = y * chord_len
        rx = px * cos_t - py * sin_t
        ry = px * sin_t + py * cos_t
        oc.add(adsk.core.Point3D.create(rx, ry, 0))
    return oc
# ─────────────────────────────────────────────────────────────────────────


def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui  = app.userInterface

        doc    = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
        design = adsk.fusion.Design.cast(app.activeProduct)
        root   = design.rootComponent
        um     = design.unitsManager
        du     = um.defaultLengthUnits     # e.g. "mm"
        to_du  = lambda metres: um.convert(metres, "m", du)

        with open(JSON_PATH, "r") as f:
            spec = json.load(f)
        sections   = sorted(spec["blades"], key=lambda s: s["radial_pos_m"])
        n_blades   = spec.get("num_blades", 1)

        base_plane = root.xYConstructionPlane
        profiles   = []

        # build every section ------------------------------------------------
        for sec in sections:
            z_off  = to_du(sec["radial_pos_m"])
            chord  = to_du(sec["chord_len_m"])
            twist  = sec["twist_rad"]

            dat_file = sec.get("coordinate_file")

            foil = None
            if dat_file:
                foil = read_dat(dat_file)
            else:
                foil = get_naca_coordinates(sec.get("naca"))

            # plane offset along +Z
            p_in = root.constructionPlanes.createInput()
            p_in.setByOffset(base_plane, adsk.core.ValueInput.createByReal(z_off))
            plane = root.constructionPlanes.add(p_in)

            # sketch
            sk = root.sketches.add(plane)
            pts2d = foil
            oc    = to_object_collection(pts2d, chord, twist)
            sk.sketchCurves.sketchFittedSplines.add(oc)
            profiles.append(sk.profiles.item(0))

        # loft blade ---------------------------------------------------------
        lofts   = root.features.loftFeatures
        loft_in = lofts.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        for pr in profiles:
            loft_in.loftSections.add(pr)
        blade_body = lofts.add(loft_in).bodies.item(0)

        # circular pattern ---------------------------------------------------
        if n_blades > 1:
            coll = adsk.core.ObjectCollection.create()
            coll.add(blade_body)
            pat  = root.features.circularPatternFeatures
            p_in = pat.createInput(coll, root.zConstructionAxis)
            p_in.quantity   = adsk.core.ValueInput.createByReal(n_blades)
            p_in.totalAngle = adsk.core.ValueInput.createByString("360 deg")
            pat.add(p_in)

        # export local .f3d ---------------------------------------------------
        exp_mgr = design.exportManager
        opts    = exp_mgr.createFusionArchiveExportOptions(SAVE_PATH)
        exp_mgr.execute(opts)

        app.quit()

    except BaseException as err:
        if ui:
            ui.messageBox("build_blade.py failed:\n" + traceback.format_exc())
        raise err