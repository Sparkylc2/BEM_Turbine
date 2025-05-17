#!/usr/bin/env python3
"""
run_fusion_build.py – launches Fusion 360 head-lessly,
drives build_blade.py, and waits for completion.

Usage
-----
    # Use all defaults set in script:
    python run_fusion_build.py

    # Specify a JSON explicitly (optional):
    python run_fusion_build.py  path/to/my_blade.json

Environment overrides (all optional)
---------------------
    FUSION360_EXE       – full path to Fusion360.exe (if auto-detection fails)
    FUSION_BLADE_JSON   – JSON path (overrides CLI / default)
    FUSION_AIRFOIL_DIR  – folder of .dat files (overrides default)
    FUSION_BLADE_SAVE   – output .f3d (overrides default)
"""

import os, sys, glob, subprocess, pathlib, platform

# ─── EDIT THESE DEFAULTS ───────────────────────────────────────────────
PROJECT_ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_JSON = PROJECT_ROOT / "naca_data" / "blade_profiles" / "blade a.json"
DEFAULT_AIRFOIL_DIR = PROJECT_ROOT / "naca_data" / "airfoil_profiles"
DEFAULT_SAVE_DIR = PROJECT_ROOT / "output"
BUILD_SCRIPT = pathlib.Path(__file__).with_name("build_blade.py")
# ─────────────────────────────────────────────────────────────────────────


def find_fusion_exe() -> str:
    env = os.getenv("FUSION360_EXE")
    print("os.getenv(FUSION360_EXE): " + env)
    if env and os.path.isfile(env):
        return env

    if platform.system() == "Darwin":  # macOS
        search_paths = [
            "/Applications/Autodesk Fusion 360.app/Contents/MacOS/Fusion360",
            f"{os.path.expanduser('~')}/Applications/Autodesk Fusion 360.app/Contents/MacOS/Fusion360",
        ]

        webdeploy = pathlib.Path(f"{os.path.expanduser('~')}/Library/Application Support/Autodesk/webdeploy")
        if webdeploy.exists():
            production = webdeploy / "production"
            if production.exists():
                for app_name in ["Autodesk Fusion 360.app", "Autodesk Fusion.app"]:
                    search_paths.append(str(production / app_name / "Contents/MacOS/Fusion360"))

                    symlink_path = production / app_name
                    if symlink_path.exists() and symlink_path.is_symlink():
                        search_paths.append(str(symlink_path / "Contents/MacOS/Fusion360"))

                    for subdir in production.glob("*"):
                        if subdir.is_dir() and not subdir.name.startswith('.'):
                            search_paths.append(str(subdir / app_name / "Contents/MacOS/Fusion360"))

    elif platform.system() == "Windows":
        search_paths = []

        local_appdata = os.getenv("LOCALAPPDATA", "")
        if local_appdata:
            root = pathlib.Path(local_appdata) / "Autodesk" / "webdeploy" / "production"
            if root.is_dir():
                for subdir in root.glob("*"):
                    if subdir.is_dir():
                        exe_path = subdir / "Fusion360.exe"
                        if exe_path.is_file():
                            search_paths.append(str(exe_path))


        program_files = [
            os.getenv("ProgramFiles", "C:\\Program Files"),
            os.getenv("ProgramFiles(x86)", "C:\\Program Files (x86)")
        ]
        for pf in program_files:
            search_paths.append(f"{pf}\\Autodesk\\Fusion 360\\Fusion360.exe")

    else:
        search_paths = []


    for path in search_paths:
        if os.path.isfile(path):
            return path


    raise FileNotFoundError(
        "Fusion 360 executable not found. Either:\n"
        "1. Install Fusion 360 to a standard location\n"
        "2. Set the environment variable FUSION360_EXE to the full path"
    )

def main():
    if len(sys.argv) > 2:
        sys.exit("Usage: run_fusion_build.py  [blade.json]")


    json_path = os.getenv("FUSION_BLADE_JSON") or \
                (sys.argv[1] if len(sys.argv) == 2 else DEFAULT_JSON)
    json_path = pathlib.Path(json_path).expanduser().resolve()
    if not json_path.is_file():
        sys.exit(f"JSON spec not found: {json_path}")


    foil_dir = os.getenv("FUSION_AIRFOIL_DIR") or DEFAULT_AIRFOIL_DIR
    foil_dir = pathlib.Path(foil_dir).expanduser().resolve()
    if not foil_dir.is_dir():

        try:
            foil_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created airfoil directory: {foil_dir}")
        except:
            sys.exit(f"Could not create airfoil directory: {foil_dir}")


    if os.getenv("FUSION_BLADE_SAVE"):
        save_path = pathlib.Path(os.getenv("FUSION_BLADE_SAVE")).expanduser().resolve()
    else:

        save_path = DEFAULT_SAVE_DIR / f"{json_path.stem}.f3d"
        save_path = save_path.expanduser().resolve()


    save_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        fusion_exe = find_fusion_exe()
    except FileNotFoundError as e:
        sys.exit(str(e))


    env = os.environ.copy()
    env["FUSION_BLADE_JSON"]  = str(json_path)
    env["FUSION_AIRFOIL_DIR"] = str(foil_dir)
    env["FUSION_BLADE_SAVE"]  = str(save_path)

    cmd = [fusion_exe, "--script", str(BUILD_SCRIPT)]

    print("Launching Fusion 360")
    print("   ", " ".join(cmd))
    print(f"   JSON: {json_path}")
    print(f"   Airfoil Directory: {foil_dir}")
    print(f"   Save Path: {save_path}")

    try:
        subprocess.run(cmd, env=env, check=True)
        print(f"\nBlade built and saved to {save_path}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Fusion exited with code {e.returncode}")


if __name__ == "__main__":
    main()