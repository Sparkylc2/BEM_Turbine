from pathlib import Path
import numpy as np
from dataclasses import dataclass
from typing import Optional, List
import json
from json import JSONEncoder
from dataclasses import asdict
# ---------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------- VARS ----------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# goes up one level to the project root, acts as prefix for the json file
ROOT_DIR = Path(__file__).resolve().parent.parent
BASE_PATH = ROOT_DIR / "naca_data" / "blade_profiles"
# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------------------------------------------- #
# --------------------------------------------------- CLASSES ---------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
@dataclass
class Blade:
    radial_pos_m: float
    chord_len_m: float
    twist_rad: float
    naca: Optional[str] = None
    coordinate_file: Optional[str] = None

@dataclass
class BladeProfile:
    blades: List[Blade]
    rotor_hub_radius_m: float
    rotor_radius_m: float
    num_blades: int
    tip_speed_ratio: float
    wind_speed_mps: float


class BladeEncoder(JSONEncoder):
    def default(self, obj):
        if hasattr(obj, '__dict__'):
            return asdict(obj)
        return super().default(obj)

# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------- FUNCS ----------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #

def generate_even_profiles(initial_profiles, profile_resolution):
    # sorts input by radial pos
    sorted_profiles = sorted(initial_profiles, key=lambda p: p["radial_pos"])

    # gen evenly spaced positions
    radial_positions = np.linspace(0.0, 1.0, profile_resolution)

    # build output
    result = []
    current_idx = 0
    for r in radial_positions:
        # find the most recent profile
        while (current_idx + 1 < len(sorted_profiles) and
               sorted_profiles[current_idx + 1]["radial_pos"] <= r):
            current_idx += 1

        profile_copy = sorted_profiles[current_idx].copy()
        profile_copy["radial_pos"] = round(float(r), 5)
        result.append(profile_copy)
    return result


def save_blade_profile_json(blade_profile: BladeProfile, name: str):
    # make sure we don't forget the file extension
    if not name.endswith('.json'):
        name += '.json'

    # create the new path
    file_path = BASE_PATH / name

    # create the new file, dump the object based on the encoder into the file
    with open(file_path, 'w') as f:
        json.dump(blade_profile, f, cls=BladeEncoder, indent=4)


