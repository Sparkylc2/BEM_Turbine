import math

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from create_blade_profile import *



def get_dat_coordinates(filename):
    filename = ROOT_DIR / "naca_data" / "airfoil_profiles" / filename

    with open(filename, 'r') as f:
        lines = f.readlines()

    coordinates = []
    for line in lines:
        line = line.strip()

        try:
            if line and not line.startswith('#'):
                parts = line.split()
                x = float(parts[0])
                y = float(parts[1])
                coordinates.append((x, y))
        except ValueError:
            continue

    return np.array(coordinates)



def transform(points, chord_len, twist_rad, z_off, shift):
    x_coord, y_coord, z_coord = [], [], []
    shift_x, shift_y, shift_z = shift
    cos_t, sin_t = math.cos(twist_rad), math.sin(twist_rad)
    for x, y in points:
        ox = shift_x + x
        oy = shift_y + y
        oz = shift_z + z_off

        px = (ox - 1.0) * chord_len
        py = oy * chord_len
        rx = px * cos_t - py * sin_t
        ry = px * sin_t + py * cos_t
        x_coord.append(rx)
        y_coord.append(ry)
        z_coord.append(oz)
    return np.array(x_coord), np.array(y_coord), np.array(z_coord)

def plot_airfoils(data):
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    for section in data:
        filename = ROOT_DIR / "naca_data" / "airfoil_profiles"/ section["coordinate_file"]
        z_off = section["radial_pos_m"]
        chord_len = section["chord_len_m"]
        twist_rad = section["twist_rad"]

        airfoil_coords = get_dat_coordinates(filename)
        x, y, z = transform(airfoil_coords, chord_len, twist_rad, z_off, (1 - 0.321, 0.0, 0.0))
        x = x/TIP_RADIUS
        y = y/TIP_RADIUS
        z = z/TIP_RADIUS
        ax.plot(np.append(x, x[0]), np.append(y, y[0]), np.append(z, z[0]), label=f'z={z_off:.2f} m')

    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    # ax.set_zlim(0, 1)
    ax.set_title('3D Airfoil Shapes')
    # ax.legend()
    plt.show()



with open(ROOT_DIR / "naca_data" / "blade_profiles" / "blade_profile_test.json", 'r') as f:
    data = json.load(f)

blades = data["blades"]

plot_airfoils(blades)

r = np.linspace(HUB_RADIUS, TIP_RADIUS, 100)
chord = np.array([CHORD_DISTRIBUTION(curr_r) for curr_r in r])
beta = np.array([TWIST_DISTRIBUTION(curr_r) for curr_r in r])

print(max(beta))




plt.figure(figsize=(12, 8))
plt.plot()
plt.figure(figsize=(10, 6))
plt.plot(r / TIP_RADIUS, chord / TIP_RADIUS, label="$c = \\frac{8 \\pi r}{B C_{l, des}} (1 - \\cos(\\phi))$", color='blue')
plt.xlabel('$\\frac{r}{R}$', fontsize=16)
plt.ylabel('$\\frac{c}{R}$', fontsize=16)
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.title('Blade Chord Length Distribution')
# Let matplotlib handle the legend
plt.legend(fontsize=16)
plt.grid()
plt.show()
# plt.savefig('blade_chord_length_distribution.svg', dpi=300)



plt.figure(figsize=(10, 6))
plt.plot(r / TIP_RADIUS, math.degrees(beta), label='$\\frac{2}{3}\\tan^{-1}(\\frac{1}{\\lambda_r \\frac{r}{R_{tip}}})$', color='red')
plt.xlabel('$\\frac{r}{R}$', fontsize=16)
plt.ylabel('$\\beta$  (rad)', fontsize=16)
plt.xlim(0, 1)
# plt.ylim(-0.5, 0.5)
plt.title('Blade Twist Angle Distribution')
plt.legend()
plt.grid()
# plt.savefig('blade_twist_angle_distribution.svg', dpi=300)
plt.show()






