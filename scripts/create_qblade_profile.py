import math
from datetime import datetime
import numpy as np
import json
from pathlib import Path
ROOT_DIR = Path(__file__).resolve().parent.parent
BASE_PATH = ROOT_DIR / "naca_data" / "blade_profiles"
FILENAME = "blade_profile_test.json"


json_path = BASE_PATH / FILENAME
def write_qblade_bld(
        filename: str,
        blade_name: str,
        rotor_type: str,
        num_blades: int,
        pos, chord, twist,
        offset_x=None, offset_y=None,
        pitch_axis=0.0,
        polar_file=""
):
    n = len(pos)
    offset_x = offset_x if offset_x is not None else [0.0] * n
    offset_y = offset_y if offset_y is not None else [0.0] * n



    now = datetime.now()

    with open(filename, "w") as f:
        f.write("-" * 40 + "QBlade Blade Definition File" + "-" * 40 + "\n")
        f.write("Generated with : QBlade CE v2.0.8.6_beta windows\n")
        f.write("Archive Format: 310035\n")
        f.write(f"Time : {now.strftime('%H:%M:%S')}\n")
        f.write(f"Date : {now.strftime('%d.%m.%Y')}\n\n")

        f.write("-" * 40 + "Object Name" + "-" * 65 + "\n")
        f.write(f"{blade_name:<50} OBJECTNAME         - the name of the blade object\n\n")

        f.write("-" * 40 + "Parameters" + "-" * 65 + "\n")
        f.write(f"{rotor_type:<50} ROTORTYPE          - the rotor type\n")
        f.write(f"{'false':<50} INVERTEDFOILS      - invert the airfoils? (only VAWT) [bool]\n")
        f.write(f"{num_blades:<50} NUMBLADES          - number of blades\n\n")

        f.write("-" * 40 + "Blade Data" + "-" * 65 + "\n")
        f.write("POS_[m]             CHORD_[m]           TWIST_[deg]         OFFSET_X_[m]        OFFSET_Y_[m]        P_AXIS [-]          POLAR_FILE          \n")

        for i in range(n):
            f.write(f"{pos[i]:<20.5f}{chord[i]:<20.5f}{twist[i]:<20.2f}"
                    f"{offset_x:<20.5f}{offset_y:<20.5f}"
                    f"{pitch_axis:<20.5f}{polar_file}\n")

    print(f"[+] QBlade .bld file saved as: {filename}")




with open(json_path, "r") as fp:
    spec = json.load(fp)

blade_name = "Blade"
rotor_type = "HAWT"
num_blades = spec["num_blades"]
off_x, off_y, off_z = tuple(spec["shift"])

blades = spec["blades"]
pos = [b["radial_pos_m"] for b in blades]
chord = [b["chord_len_m"] for b in blades]
twist = [b["twist_rad"] for b in blades]

write_qblade_bld(BASE_PATH/"blade.bld", blade_name, rotor_type, num_blades, pos, chord, np.degrees(np.array(twist)), off_x, off_y)


