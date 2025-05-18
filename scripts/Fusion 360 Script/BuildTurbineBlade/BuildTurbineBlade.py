import os, json, math, traceback
from pathlib import Path
import adsk.core, adsk.fusion

HUB_RADIUS = 0.1
TIP_RADIUS = 0.3
CHORD_DISTRIBUTION = lambda r:  3.0 * (1 - 0.7 * r**1.5)
TWIST_DISTRIBUTION = lambda r: math.radians(20 * (1 - r)**1.2)


def generate_le_guide(root, to_du):
    le_sketch = root.sketches.add(root.xYConstructionPlane)
    le_obj_collection = adsk.core.ObjectCollection.create()

    r = linspace(0.0, 1.0, 100)

    r_arr = [to_du(curr_r * (TIP_RADIUS - HUB_RADIUS) + HUB_RADIUS) for curr_r in r]
    chord_arr = [to_du(CHORD_DISTRIBUTION(curr_r)) for curr_r in r]
    twist_arr = [TWIST_DISTRIBUTION(curr_r) for curr_r in r]

    pts = []
    for (chord, twist, r) in zip(chord_arr, twist_arr, r_arr):
        pts.append(transform([(1.0, 0.0)], chord, twist, r))

    for pt in pts:
        le_obj_collection.add(le_sketch.modelToSketchSpace(pt))

    le_spline = le_sketch.sketchCurves.sketchFittedSplines.add(le_obj_collection)
    return root.features.createPath(le_spline)

def generate_te_guide(root, to_du):
    rail_sk   = root.sketches.add(root.yZConstructionPlane)
    z_min = to_du(HUB_RADIUS)
    z_max = to_du(TIP_RADIUS)
    # z_min = to_du(sections[0]["radial_pos_m"])
    # z_max = to_du(sections[-1]["radial_pos_m"])
    p0 = rail_sk.modelToSketchSpace(adsk.core.Point3D.create(0, 0, z_min))
    p1 = rail_sk.modelToSketchSpace(adsk.core.Point3D.create(0, 0, z_max))
    rail_line = rail_sk.sketchCurves.sketchLines.addByTwoPoints(p0, p1)
    return root.features.createPath(rail_line)















def linspace(start, stop, np):
    """
    Emulate Matlab linspace
    """
    return [start + (stop - start) * i / (np - 1) for i in range(np)]
def interpolate(xa, ya, queryPoints):
    """
    A cubic spline interpolation on a given set of points (x,y)
    Recalculates everything on every call which is far from efficient but does the job for now
    should eventually be replaced by an external helper class
    """

    # PreCompute() from Paint Mono which in turn adapted:
    # NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
    # ISBN 0-521-43108-5, page 113, section 3.3.
    # http://paint-mono.googlecode.com/svn/trunk/src/PdnLib/SplineInterpolator.cs

    # number of points
    n = len(xa)
    u, y2 = [0] * n, [0] * n

    for i in range(1, n - 1):
        # This is the decomposition loop of the tridiagonal algorithm.
        # y2 and u are used for temporary storage of the decomposed factors.

        wx = xa[i + 1] - xa[i - 1]
        sig = (xa[i] - xa[i - 1]) / wx
        p = sig * y2[i - 1] + 2.0

        y2[i] = (sig - 1.0) / p

        ddydx = (ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]) - (ya[i] - ya[i - 1]) / (xa[i] - xa[i - 1])

        u[i] = (6.0 * ddydx / wx - sig * u[i - 1]) / p

    y2[n - 1] = 0

    # This is the backsubstitution loop of the tridiagonal algorithm
    # ((int i = n - 2; i >= 0; --i):
    for i in range(n - 2, -1, -1):
        y2[i] = y2[i] * y2[i + 1] + u[i]

    # interpolate() adapted from Paint Mono which in turn adapted:
    # NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
    # ISBN 0-521-43108-5, page 113, section 3.3.
    # http://paint-mono.googlecode.com/svn/trunk/src/PdnLib/SplineInterpolator.cs

    results = [0] * n

    # loop over all query points
    for i in range(len(queryPoints)):
        # bisection. This is optimal if sequential calls to this
        # routine are at random values of x. If sequential calls
        # are in order, and closely spaced, one would do better
        # to store previous values of klo and khi and test if

        klo = 0
        khi = n - 1

        while khi - klo > 1:
            k = (khi + klo) >> 1
            if xa[k] > queryPoints[i]:
                khi = k
            else:
                klo = k

        h = xa[khi] - xa[klo]
        a = (xa[khi] - queryPoints[i]) / h
        b = (queryPoints[i] - xa[klo]) / h

        # Cubic spline polynomial is now evaluated.
        results[i] = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2[klo] + (b * b * b - b) * y2[khi]) * (h * h) / 6.0

    return results
def naca4(number, n, finite_te=False, half_cosine_spacing=False):
    """
    Returns 2*n+1 points in [0 1] for the given 4 digit NACA number string
    """

    m = float(number[0]) / 100.0
    p = float(number[1]) / 10.0
    t = float(number[2:]) / 100.0

    a0 = +0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = +0.2843

    if finite_te:
        a4 = -0.1015  # For finite thick TE
    else:
        a4 = -0.1036  # For zero thick TE

    if half_cosine_spacing:
        beta = linspace(0.0, math.pi, n + 1)
        x = [(0.5 * (1.0 - math.cos(xx))) for xx in beta]  # Half cosine based spacing
    else:
        x = linspace(0.0, 1.0, n + 1)

    yt = [5 * t * (a0 * math.sqrt(xx) + a1 * xx + a2 * pow(xx, 2) + a3 * pow(xx, 3) + a4 * pow(xx, 4)) for xx in x]

    xc1 = [xx for xx in x if xx <= p]
    xc2 = [xx for xx in x if xx > p]

    if p == 0:
        xu = x
        yu = yt

        xl = x
        yl = [-xx for xx in yt]

        #  xc = xc1 + xc2
        #  zc = [0] * len(xc)
    else:
        yc1 = [m / pow(p, 2) * xx * (2 * p - xx) for xx in xc1]
        yc2 = [m / pow(1 - p, 2) * (1 - 2 * p + xx) * (1 - xx) for xx in xc2]
        zc = yc1 + yc2

        dyc1_dx = [m / pow(p, 2) * (2 * p - 2 * xx) for xx in xc1]
        dyc2_dx = [m / pow(1 - p, 2) * (2 * p - 2 * xx) for xx in xc2]
        dyc_dx = dyc1_dx + dyc2_dx

        theta = [math.atan(xx) for xx in dyc_dx]

        xu = [xx - yy * math.sin(zz) for xx, yy, zz in zip(x, yt, theta)]
        yu = [xx + yy * math.cos(zz) for xx, yy, zz in zip(zc, yt, theta)]

        xl = [xx + yy * math.sin(zz) for xx, yy, zz in zip(x, yt, theta)]
        yl = [xx - yy * math.cos(zz) for xx, yy, zz in zip(zc, yt, theta)]

    X = xu[::-1] + xl[1:]
    Z = yu[::-1] + yl[1:]

    return X, Z
def naca5(number, n, finite_te=False, half_cosine_spacing=False):
    """
    Returns 2*n+1 points in [0 1] for the given 5 digit NACA number string
    """

    naca1 = int(number[0])
    naca23 = int(number[1:3])
    naca45 = int(number[3:])

    cld = naca1 * (3.0 / 2.0) / 10.0
    p = 0.5 * naca23 / 100.0
    t = naca45 / 100.0

    a0 = +0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = +0.2843

    if finite_te:
        a4 = -0.1015  # For finite thickness trailing edge
    else:
        a4 = -0.1036  # For zero thickness trailing edge

    if half_cosine_spacing:
        beta = linspace(0.0, math.pi, n + 1)
        x = [(0.5 * (1.0 - math.cos(x))) for x in beta]  # Half cosine based spacing
    else:
        x = linspace(0.0, 1.0, n + 1)

    yt = [5 * t * (a0 * math.sqrt(xx) + a1 * xx + a2 * pow(xx, 2) + a3 * pow(xx, 3) + a4 * pow(xx, 4)) for xx in x]

    P = [0.05, 0.1, 0.15, 0.2, 0.25]
    M = [0.0580, 0.1260, 0.2025, 0.2900, 0.3910]
    K = [361.4, 51.64, 15.957, 6.643, 3.230]

    m = interpolate(P, M, [p])[0]
    k1 = interpolate(M, K, [m])[0]

    xc1 = [xx for xx in x if xx <= p]
    xc2 = [xx for xx in x if xx > p]
    #  xc = xc1 + xc2

    if p == 0:
        xu = x
        yu = yt

        xl = x
        yl = [-x for x in yt]

        #  zc = [0] * len(xc)
    else:
        yc1 = [k1 / 6.0 * (pow(xx, 3) - 3 * m * pow(xx, 2) + pow(m, 2) * (3 - m) * xx) for xx in xc1]
        yc2 = [k1 / 6.0 * pow(m, 3) * (1 - xx) for xx in xc2]
        zc = [cld / 0.3 * xx for xx in yc1 + yc2]

        dyc1_dx = [cld / 0.3 * (1.0 / 6.0) * k1 * (3 * pow(xx, 2) - 6 * m * xx + pow(m, 2) * (3 - m)) for xx in xc1]
        dyc2_dx = [cld / 0.3 * (1.0 / 6.0) * k1 * pow(m, 3)] * len(xc2)

        dyc_dx = dyc1_dx + dyc2_dx
        theta = [math.atan(xx) for xx in dyc_dx]

        xu = [xx - yy * math.sin(zz) for xx, yy, zz in zip(x, yt, theta)]
        yu = [xx + yy * math.cos(zz) for xx, yy, zz in zip(zc, yt, theta)]

        xl = [xx + yy * math.sin(zz) for xx, yy, zz in zip(x, yt, theta)]
        yl = [xx - yy * math.cos(zz) for xx, yy, zz in zip(zc, yt, theta)]

    X = xu[::-1] + xl[1:]
    Z = yu[::-1] + yl[1:]

    return X, Z
def naca(number, n, finite_te=False, half_cosine_spacing=False):
    if len(number) == 4:
        return naca4(number, n, finite_te, half_cosine_spacing)
    elif len(number) == 5:
        return naca5(number, n, finite_te, half_cosine_spacing)
    else:
        raise Exception

def read_dat(path):
    x = []
    y = []
    with open(path, "r") as fp:
        for line in fp:
            try:
                values = line.strip().split()
                x.append(float(values[0]))
                y.append(float(values[1]))
            except ValueError:
                continue
        if len(x) > 0 and (x[0] != x[-1] or y[0] != y[-1]):
            x.append(x[0])
            y.append(y[0])
    return (x, y)




def get_naca_coordinates(naca_code):
    pts = naca(naca_code, 200)
    return pts
def reorder_to_te_clockwise(points):
    if isinstance(points[0], tuple):
        pts = points
    else:
        pts = list(zip(*points))

    max_x = max(p[0] for p in pts)
    te_candidates = [i for i, p in enumerate(pts) if abs(p[0] - max_x) < 1e-6]
    start = te_candidates[0]

    ordered = pts[start:] + pts[1:start+1]

    area = 0.0
    for (x0, y0), (x1, y1) in zip(ordered, ordered[1:]):
        area += (x1 - x0) * (y1 + y0)
    if area > 0:
        ordered = [ordered[0]] + ordered[:0:-1]

    return ordered


def transform(points, chord_len, twist_rad, z_off = None):
    pts = []
    cos_t, sin_t = math.cos(twist_rad), math.sin(twist_rad)
    for x, y in points:
        px = (x - 1.0) * chord_len
        py = y * chord_len
        rx = px * cos_t - py * sin_t
        ry = px * sin_t + py * cos_t
        pts.append(adsk.core.Point3D.create(rx, ry, 0.0 if z_off is None else z_off))
    return pts
def to_object_collection(points, chord_len, twist_rad):
    oc = adsk.core.ObjectCollection.create()
    pts = transform(points, chord_len, twist_rad, z_off = None)
    for pt in pts:
        oc.add(pt)
    return oc

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui  = app.userInterface

        fd = ui.createFileDialog()
        fd.isMultiSelectEnabled = False
        fd.title  = "Select blade JSON specification"
        fd.filter = "JSON (*.json)"
        if fd.showOpen() != adsk.core.DialogResults.DialogOK:
            return
        json_path = Path(fd.filename)



        ROOT_DIR = Path(__file__).resolve().parent.parent
        airfoil_dir = ROOT_DIR / "naca_data" / "airfoil_profiles"
        save_path = ROOT_DIR / "naca_data" / "blade_models" / f"{Path(json_path).stem}.f3d"




        #
        #
        # fol = ui.createFolderDialog()
        # fol.title = "Select folder containing .dat air-foil files"
        # if fol.showDialog() != adsk.core.DialogResults.DialogOK:
        #     return
        # airfoil_dir = Path(fol.folder)
        #
        # sf = ui.createFileDialog()
        # sf.title  = "Save blade as Fusion archive"
        # sf.filter = "Fusion Archive (*.f3d)"
        # sf.initialFilename = json_path.stem + ".f3d"
        # if sf.showSave() != adsk.core.DialogResults.DialogOK:
        #     return
        # save_path = Path(sf.filename)

        doc    = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
        design = adsk.fusion.Design.cast(app.activeProduct)
        root   = design.rootComponent
        um     = design.unitsManager
        du     = um.defaultLengthUnits
        to_du  = lambda m: um.convert(m, "m", du)

        with open(json_path, "r") as fp:
            spec = json.load(fp)

        sections = sorted(spec["blades"], key=lambda s: s["radial_pos_m"])
        n_blades = spec.get("num_blades", 1)

        base_plane = root.xYConstructionPlane
        profiles   = []

        for sec in sections:
            z_off  = to_du(sec["radial_pos_m"])
            chord  = to_du(sec["chord_len_m"])
            twist  = sec["twist_rad"]

            if sec.get("coordinate_file"):
                foil = airfoil_dir / sec["coordinate_file"].with_suffix(".dat")
                pts2d = reorder_to_te_clockwise(read_dat(foil))
            else:
                pts2d = reorder_to_te_clockwise(get_naca_coordinates(sec.get("naca")))

            p_in = root.constructionPlanes.createInput()
            p_in.setByOffset(base_plane, adsk.core.ValueInput.createByReal(z_off))
            plane = root.constructionPlanes.add(p_in)

            sk = root.sketches.add(plane)
            sk.sketchCurves.sketchFittedSplines.add(
                to_object_collection(pts2d, chord, twist)
            )
            profiles.append(sk.profiles.item(0))



        lofts  = root.features.loftFeatures
        lf_in  = lofts.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)

        for pr in profiles:
            lf_in.loftSections.add(pr)

        le_path = generate_le_guide(root, to_du)
        rail_path = generate_te_guide(root, to_du)

        lofts       = root.features.loftFeatures
        blade_body  = None

        lf_in = lofts.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation if blade_body is None else adsk.fusion.FeatureOperations.JoinFeatureOperation)
        for pr in profiles:
            lf_in.loftSections.add(pr)

        lf_in.centerLineOrRails.addCenterLine(le_path)
        lf_in.centerLineOrRails.addRail(rail_path)

        # for i in range(1, len(profiles)):
        #     lf_in = lofts.createInput(
        #         adsk.fusion.FeatureOperations.NewBodyFeatureOperation
        #         if blade_body is None
        #         else adsk.fusion.FeatureOperations.JoinFeatureOperation)
        #
        #     lf_in.loftSections.add(profiles[i - 1])
        #     lf_in.loftSections.add(profiles[i])
        #     lf_in.centerLineOrRails.addCenterLine(le_path)
        #     lf_in.centerLineOrRails.addRail(rail_path)
        #
        #     loft = lofts.add(lf_in)
        #
        #     if blade_body is None:
        #         blade_body = loft.bodies.item(0)


        # exp_mgr = design.exportManager
        # opts    = exp_mgr.createFusionArchiveExportOptions(str(save_path))
        # exp_mgr.execute(opts)

        ui.messageBox(f"Blade created and saved to:\n{save_path}")

    except Exception:
        if ui:
            ui.messageBox("Blade add-in failed:\n" + traceback.format_exc())
        raise


def stop(context):
    pass
