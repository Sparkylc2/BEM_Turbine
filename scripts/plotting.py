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


def plot_chord_distribution(data):

    chord_eqn = 'c(r) = $\\frac{\\mathbf{8} \\pi r(\\mathbf{1} - \\mathbf{cos\\mathbf{φ}})} {BC_{L,des}}$'
    fig, ax = plt.subplots(figsize=(10, 10))
    data, R, CHORD = data
    ax.plot(R*100, CHORD*100,
            label=chord_eqn,
            color=red_color,
            linestyle='-',
            linewidth=3,  # Increased thickness
            solid_capstyle='round',
            )

    ax.set_xlabel('Radial Position (mm)', fontsize=20)
    ax.set_ylabel('Chord Length (mm)', fontsize=20, labelpad = 20)
    ax.set_title('Chord Distribution', fontsize=20)
    ax.legend(loc='upper right',
              frameon=True,
              framealpha=0.7,
              edgecolor='black',
              fontsize=20)
    ax.set_xlim(min(R)*100, max(R)*100)
    ax.set_ylim(min(R)*100, max(R)*100)
    apply_grid_styling(ax)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    plt.tight_layout()
    plt.show()

def plot_twist_distribution(data):
    twist_eqn = '$\\mathbf{β(r)} = \\mathbf{\\frac{\\mathbf{2}}{\\mathbf{3}} arctan\\;\\frac{R_{tip}}{λ_r r} - α_{des}}$'
    fig, ax = plt.subplots()
    data, R, TWIST = data

    ax.plot(R * 100, TWIST * 180/math.pi,
            label=twist_eqn,
            color=blue_color,
            linestyle='-',
            linewidth=3,
            solid_capstyle='round',
            )

    ax.set_xlabel('Radial Position (mm)', fontsize=20)
    ax.set_ylabel('Twist  (deg)', fontsize=20, labelpad = 20)
    ax.set_title('Twist Distribution', fontsize=20)
    ax.legend(loc='upper right',
              frameon=True,
              framealpha=0.7,
              edgecolor='black',
              fontsize=20)
    # ax.set_xlim(0, 0.227)
    # ax.set_ylim(0, 0.15)
    ax.set_xlim(min(R) * 100, max(R) * 100)
    ax.set_ylim(min(TWIST * 180/math.pi), max(TWIST * 180/math.pi))
    apply_grid_styling(ax)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    plt.tight_layout()
    plt.show()

def plot_chord_distribution_with_airfoils(data, save_filename=None, elev=10, azim=300):
    data, R, CHORD = data
    blades = data["blades"]
    shift = tuple(data["shift"])
    shift = (0, 0, 0)

    fig = plt.figure(figsize=(10, 6.5))
    ax = fig.add_subplot(111, projection='3d')


    ax.view_init(elev=20, azim=300)  # Adjust these values for best view

    start_z = 10.0
    end_z = -10.0
    min_y = 10.0
    max_y = -10.0

    # Plot airfoil profiles
    for section in blades:
        filename = ROOT_DIR / "naca_data" / "airfoil_profiles" / section["coordinate_file"]
        z_off = section["radial_pos_m"]
        chord_len = section["chord_len_m"]

        airfoil_coords = get_dat_coordinates(filename)
        x, y, z = transform(airfoil_coords, chord_len, section["twist_rad"], z_off, shift)

        if min(z) < start_z:
            start_z = min(z)
        if max(z) > end_z:
            end_z = max(z)

        if min(y) < min_y:
            min_y = min(y)
        if max(y) > max_y:
            max_y = max(y)

        # Plot the airfoil profiles
        ax.plot(np.append(z, z[0]), np.append(y, y[0]), -np.append(x, x[0]),
                color='black',
                alpha=0.8,
                linestyle='-',  # Changed from dashed to solid
                linewidth=0.8,
                solid_capstyle='round',
                )

    # Add reference plane for chord distribution
    z_min, z_max = ax.get_zlim()
    x_min, x_max = start_z, end_z
    # Optional - create a slightly transparent reference plane
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, 2), np.linspace(z_min, z_max, 2))
    ax.plot_surface(xx, np.zeros_like(xx), yy, alpha=0.1, color='gray')

    # Plot chord distribution with thicker, more prominent line
    ax.plot(R, 0, CHORD,
            label='Chord Distribution',
            color=red_color,
            linestyle='-',
            linewidth=3,  # Increased thickness
            solid_capstyle='round',
            zorder=100,  # Ensure it's on top
            )


    ax.plot([0, 0.227], [0, 0], [0, 0], 'k-', alpha=0.3, linewidth=0.5)

    # Set tight limits
    ax.set_xlim(0, 0.227)
    ax.set_zlim(0, 0.15)
    ax.set_ylim(0, 0.15)

    # Orthographic projection for cleaner display
    ax.set_proj_type('ortho')

    # Create a custom legend without a box
    from matplotlib.lines import Line2D
    custom_line = Line2D([0], [0], color=red_color, lw=3)
    ax.legend([custom_line], ['Chord Distribution'],
              loc='upper right',
              frameon=True,
              framealpha=0.7,
              edgecolor='none')

    # Remove title for cleaner look
    # ax.set_title('Blade Geometry', fontsize=20, pad=20)

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
plot_chord_distribution((data, R, CHORD))
plot_twist_distribution((data, R, BETA))