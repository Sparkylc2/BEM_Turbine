# config.py
import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import FixedLocator
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.ticker import ScalarFormatter
from scipy.optimize import fsolve

from create_blade_profile import *

def configure_plotting():
    fm.fontManager.addfont("/Library/Fonts/SF-Pro-Text-Light.otf")
    fm.fontManager.addfont("/Library/Fonts/SF-Pro-Text-Regular.otf")
    fm.fontManager.addfont("/Library/Fonts/SF-Pro-Text-Bold.otf")

    sf_pro_font = fm.FontProperties(fname="/Library/Fonts/SF-Pro-Text-Ultralight.otf")
    sf_pro_bold_font = fm.FontProperties(fname="/Library/Fonts/SF-Pro-Text-Bold.otf")
    sf_pro_italic_font = fm.FontProperties(fname="/Library/Fonts/SF-Pro-Text-RegularItalic.otf")

    plt.rcParams["font.family"] = sf_pro_font.get_name()
    plt.rcParams["font.weight"] = "light"
    plt.rcParams["mathtext.fontset"] = "custom"

    plt.rcParams["mathtext.rm"] = sf_pro_italic_font.get_name() + ":weight=light"+":style=normal"
    plt.rcParams["mathtext.it"] = sf_pro_italic_font.get_name() + ":weight=regular"+":style=italic"
    plt.rcParams["mathtext.bf"] = sf_pro_italic_font.get_name() + ":weight=regular"+":style=normal"
    plt.rcParams["mathtext.sf"] = sf_pro_italic_font.get_name() + ":weight=light"+":style=italic"
    plt.rcParams["lines.antialiased"] = True

    plt.rcParams["lines.linewidth"] = 2
    plt.rcParams["lines.markersize"] = 8
    plt.rcParams["errorbar.capsize"] = 5
    plt.rcParams["lines.markeredgewidth"] = 1.25

    plt.rcParams["legend.frameon"] = True
    plt.rcParams["legend.framealpha"] = 1.0
    plt.rcParams["legend.edgecolor"] = 'black'
    plt.rcParams["legend.fontsize"] = 16
    plt.rcParams["legend.handlelength"] = 3

    plt.rcParams["axes.grid"] = True
    plt.rcParams["grid.alpha"] = 0.8
    plt.rcParams["grid.linestyle"] = '-'
    plt.rcParams["grid.linewidth"] = 0.7

    plt.rcParams["axes.labelsize"] = 20
    plt.rcParams["xtick.labelsize"] = 16
    plt.rcParams["ytick.labelsize"] = 16
    plt.rcParams["axes.spines.right"] = True
    plt.rcParams["axes.spines.top"] = True
    plt.rcParams["axes.linewidth"] = 1

    plt.rcParams["xtick.direction"] = 'out'
    plt.rcParams["ytick.direction"] = 'out'
    plt.rcParams["xtick.major.pad"] = 10
    plt.rcParams["ytick.major.pad"] = 10

    plt.rcParams["axes.facecolor"] = '#f9f9f9'

    plt.rcParams["figure.figsize"] = (10, 6.5)
    plt.rcParams["figure.subplot.left"] = 0.05
    plt.rcParams["figure.subplot.right"] = 0.95
    plt.rcParams["figure.subplot.top"] = 0.95
    plt.rcParams["figure.subplot.bottom"] = 0.05

    plt.rcParams["lines.markersize"] = 8
    plt.rcParams["lines.markeredgewidth"] = 1
    plt.rcParams["errorbar.capsize"] = 5


def apply_grid_styling(ax):
    ax.grid(which='major', linestyle='-', linewidth=0.7, alpha=0.8)
    ax.grid(which='minor', linestyle='--', linewidth=0.5, alpha=0.5)

    for spine in ax.spines.values():
        spine.set_linewidth(1)



def setup_3d_plot(ax, elev=25, azim=300):
    """Configure a 3D axes with consistent styling"""
    # Set orthographic projection
    ax.set_proj_type('ortho')

    # Set default viewing angle
    ax.view_init(elev=elev, azim=azim)

    # Configure tick parameters
    for axis in ['x', 'y', 'z']:
        ax.tick_params(axis=axis, pad=10, labelsize=16)

    # Configure grid
    ax.grid(which='major', linestyle=':', linewidth=0.7, alpha=0.8)
    ax.grid(which='minor', linestyle='--', linewidth=0.5, alpha=0.5)

    # Ensure consistent spine styling
    for spine in ax.spines.values():
        spine.set_linewidth(1)

    return ax


def apply_3d_labels(ax, x_label="X-axis", y_label="Y-axis", z_label="Z-axis", fontsize=20):
    label_padding = 1.1

    xlim, ylim, zlim = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()

    ax.text(1.04 * (xlim[0] + xlim[1]) / 2, ylim[1], 1.0025 * label_padding * zlim[1],
            x_label, ha='center', va='center', fontsize=fontsize)
    ax.text(xlim[0], (ylim[0] + ylim[1]) / 2, 1.025 * label_padding * zlim[1],
            y_label, ha='center', va='center', fontsize=fontsize)
    ax.text(label_padding * 0.9 * xlim[0], ylim[0] * label_padding * 0.999,
            1.1 * (zlim[0] + zlim[1]) / 2, z_label, ha='center', va='center', fontsize=fontsize)


def save_plot(fig, filename, dpi=300):
    fig.savefig(f"{SAVE_DIR}/{filename}.png", format="png", dpi=dpi)
    plt.close(fig)

configure_plotting()


red_color = (1, 191 / 255, 191 / 255)
blue_color = (191 / 255, 223 / 255, 1)
green_color = (191 / 255, 1, 191 / 255)
yellow_color = (1, 223 / 255, 191 / 255)

ROOT_DIR = Path(__file__).resolve().parent.parent
BASE_DIR = ROOT_DIR / "naca_data" / "blade_profiles"
SAVE_DIR = ROOT_DIR / "scripts" / "plots"



def get_dat_coordinates(filename):
    filename = BASE_DIR / filename

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


def plot_chord_distribution_with_airfoils(data, save_filename=None, elev=10, azim=300):
    data, R, CHORD = data
    blades = data["blades"]
    shift = tuple(data["shift"])
    shift = (0, 0, 0)

    fig = plt.figure(figsize=(10, 6.5))
    ax = fig.add_subplot(111, projection='3d')

    setup_3d_plot(ax, elev=elev, azim=azim)

    start_z = 10.0
    end_z = -10.0

    min_y = 10.0
    max_y = -10.0

    for section in blades:
        filename = ROOT_DIR / "naca_data" / "airfoil_profiles" / section["coordinate_file"]
        z_off = section["radial_pos_m"]
        chord_len = section["chord_len_m"]

        airfoil_coords = get_dat_coordinates(filename)
        x, y, z = transform(airfoil_coords, chord_len, 0.0, z_off, shift)
        x = x
        y = y
        z = z

        if min(z) < start_z:
            start_z = min(z)
        if max(z) > end_z:
            end_z = max(z)

        if min(y) < min_y:
            min_y = min(y)
        if max(y) > max_y:
            max_y = max(y)



        ax.plot(np.append(z, z[0]), np.append(y, y[0]), -np.append(x, x[0]),
                color='black',
                alpha=0.8,
                linestyle='--',
                linewidth=0.8,
                markeredgewidth=1,
                solid_capstyle='round',
                )


    ax.plot(R, 0, CHORD,
            label='Chord Distribution',
            color=red_color,
            linestyle='-',
            linewidth=2,
            markeredgewidth=1,
            solid_capstyle='round',
            )

    ax = plt.gca()
    y_min = min_y * 0.9 if min_y > 0 else min_y * 1.1
    y_max = max_y * 0.9 if max_y > 0 else max_y * 1.1

    ax.set_xlim(0, 0.227)
    ax.set_zlim(0, 0.15)
    ax.set_ylim(0, 0.15)

    ax.legend(loc='upper right', fontsize=16, frameon=True,
              edgecolor='black', framealpha=1.0,
              bbox_to_anchor=(1.35, 0.8))

    ax.set_title('3D Airfoil Shapes', fontsize=20, pad=20)

    fig.tight_layout()

    if save_filename:
        save_plot(fig, save_filename, dpi=2000)
    else:
        plt.show()

    return fig, ax

R = np.linspace(HUB_RADIUS, TIP_RADIUS, 100)
CHORD = np.array([CHORD_DISTRIBUTION(r) for r in R])
BETA = np.array([TWIST_DISTRIBUTION(r) for r in R])


with open(ROOT_DIR / "naca_data" / "blade_profiles" / "blade_profile_test.json", 'r') as f:
    data = json.load(f)

plot_chord_distribution_with_airfoils((data, R, CHORD))