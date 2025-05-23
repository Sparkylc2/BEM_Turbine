import re

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.ticker as mticker
import pandas as pd
from scipy.interpolate import splprep, splev

from create_blade_profile import *


ROOT_DIR = Path(__file__).resolve().parent.parent
BASE_DIR = ROOT_DIR / "naca_data" / "blade_profiles"
SAVE_DIR = ROOT_DIR / "scripts" / "plots"


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
    ax.set_proj_type('ortho')

    ax.view_init(elev=elev, azim=azim)

    for axis in ['x', 'y', 'z']:
        ax.tick_params(axis=axis, pad=10, labelsize=16)

    ax.grid(which='major', linestyle=':', linewidth=0.7, alpha=0.8)
    ax.grid(which='minor', linestyle='--', linewidth=0.5, alpha=0.5)

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

def parse_qblade_file(filename):
    filepath = ROOT_DIR / "naca_data" / "bem_data" / "NREL_DATA_QBLADE"/ filename
    x = []
    y = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('New Turbine') or line.startswith('Windspeed'):
                continue
            parts = re.split(r'\s+|\t', line.strip())
            if len(parts) >= 2 and parts[0] != 'END':
                try:
                    curr_x = float(parts[0])
                    curr_y = float(parts[1])
                    x.append(curr_x)
                    y.append(curr_y)
                except ValueError:
                    continue

    return x, y

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

def get_bem_data(filename):
    filename = BASE_DIR / filename
    df = pd.read_csv(filename)
    return df.to_numpy()


def plot_chord_distribution(data):

    chord_eqn = 'c(r) = $\\frac{\\mathbf{8} \\pi r(\\mathbf{1} - \\mathbf{cos\\mathbf{φ}})} {BC_{L,des}}$'
    fig, ax = plt.subplots(figsize=(10, 10))
    data, R, CHORD = data
    ax.plot(R*1000, CHORD*1000,
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
    ax.set_xlim(min(R)*1000, max(R)*1000)
    ax.set_ylim(0, max(R)*1000)
    apply_grid_styling(ax)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    plt.tight_layout()
    plt.show()
def plot_twist_distribution(data):
    twist_eqn = '$\\mathbf{β(r)} = \\mathbf{\\frac{\\mathbf{2}}{\\mathbf{3}} arctan\\;\\frac{R_{tip}}{λ_r r} - α_{des}}$'
    fig, ax = plt.subplots()
    data, R, TWIST = data

    ax.plot(R * 1000, TWIST * 180/math.pi,
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
    apply_grid_styling(ax)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())

    ax.set_xlim(min(R) * 1000, max(R) * 1000)
    ax.set_ylim(0, max(TWIST * 180/math.pi))
    # ax.set_ylim(min(TWIST * 180/math.pi), max(TWIST * 180/math.pi))
    plt.tight_layout()
    plt.show()



def plot_blade_with_distribution_info(data, save_filename=None, elev=10, azim=300):
    data, R, CHORD = data
    blades = data["blades"]
    shift = tuple(data["shift"])

    fig = plt.figure(figsize=(10, 6.5))
    ax = fig.add_subplot(111, projection='3d')


    ax.view_init(elev=30, azim=140)


    airfoil_zs = []
    te_pts = []
    le_pts = []

    for section in blades:
        filename = ROOT_DIR / "naca_data" / "airfoil_profiles" / section["coordinate_file"]
        z_off = section["radial_pos_m"]
        chord_len = section["chord_len_m"]

        airfoil_zs.append(z_off)

        airfoil_coords = get_dat_coordinates(filename)
        x, y, z = transform(airfoil_coords, chord_len, section["twist_rad"], z_off, shift)
        le_idx = np.argmin(airfoil_coords[:, 0])
        le_point = airfoil_coords[le_idx]
        le_idx = np.argmax(airfoil_coords[:, 0])
        te_point = airfoil_coords[le_idx]

        le_x, le_y, le_z = transform(le_point[None, :], chord_len, section["twist_rad"], z_off, shift)
        te_x, te_y, te_z = transform(te_point[None, :], chord_len, section["twist_rad"], z_off, shift)

        te_pts.append([te_x[0], te_y[0], te_z[0]])
        le_pts.append([le_x[0], le_y[0], le_z[0]])

        ax.plot(np.append(z, z[0])*1000, np.append(y, y[0])*1000, -np.append(x, x[0])*1000,
                color='black',
                alpha=0.8,
                linestyle='--',
                linewidth=0.8,
                solid_capstyle='round',
                )

    center_line_x = np.array([HUB_RADIUS, TIP_RADIUS * 1000])
    center_line_yz = np.array([0, 0])

    ax.plot(center_line_x, center_line_yz, center_line_yz,
            color=red_color,
            linestyle='--',
            linewidth=2,
            solid_capstyle='round',
            label='Rotation Axis',
            alpha=1.0
            )


    te_pts = np.array(te_pts)
    le_pts = np.array(le_pts)

    te_x = te_pts[:, 0]
    te_y = te_pts[:, 1]
    te_z = te_pts[:, 2]

    le_x = le_pts[:, 0]
    le_y = le_pts[:, 1]
    le_z = le_pts[:, 2]

    tck, u = splprep([le_y, le_z], s=0)
    u_fine = np.linspace(0, 1, 100)
    y_spline, z_spline = splev(u_fine, tck)

    x_spline = np.interp(u_fine, np.linspace(0, 1, len(le_x)), le_x)

    ax.plot(z_spline * 1000, y_spline * 1000, -x_spline * 1000, color="black", linewidth=0.8, label='Chord Line', alpha = 0.8)

    tck, u = splprep([te_y, te_z], s=0)
    y_spline, z_spline = splev(u_fine, tck)

    x_spline = np.interp(u_fine, np.linspace(0, 1, len(te_x)), te_x)

    ax.plot(z_spline * 1000, y_spline * 1000, -x_spline * 1000, color="black", linewidth=0.8, alpha = 0.8)


    ax.scatter(np.array(airfoil_zs) * 1000, np.zeros(len(airfoil_zs)), np.zeros(len(airfoil_zs)),
               marker='o', color=green_color, edgecolor="black", label='Max Thickness on Camber Line')
    ax.set_xlim(0.0, 227)
    ax.set_proj_type('ortho')
    ax.legend(loc='upper right', frameon=True, framealpha=0.7, edgecolor='black', fontsize=10)

    fig.tight_layout()

    ax.set_aspect('equal')
    if save_filename:
        save_plot(fig, save_filename, dpi=2000)
    else:
        plt.show()

    return fig, ax


def plot_bem_verification():

    df = pd.read_csv(ROOT_DIR / "naca_data" / "bem_data" / "NREL_DATA" / "NREL_Reference_5MW_126.csv")

    nrel_wind_speed = df["Wind Speed [m/s]"]
    nrel_power = df["Power [kW]"]
    nrel_cp = df["Cp [-]"]
    nrel_thrust = df["Thrust [kN]"]

    # bem_data = pd.read_csv(ROOT_DIR / "naca_data" / "bem_data" / "NREL_BEM.csv")
    all_corr_bem_data = pd.read_csv(ROOT_DIR / "naca_data" / "bem_data" / "NREL_BEM_ALL_CORRECTIONS.csv")
    all_corr_bem_wind_speed = all_corr_bem_data["WIND SPEED (m/s)"]
    all_corr_bem_thrust = all_corr_bem_data[" THRUST (N)"]
    all_corr_bem_power = all_corr_bem_data[" PRODUCED POWER (W)"]
    all_corr_bem_cp = all_corr_bem_data[" C_P"]

    mach_corr_bem_data = pd.read_csv(ROOT_DIR / "naca_data" / "bem_data" / "NREL_BEM_MACH_CORRECTION.csv")
    mach_corr_bem_wind_speed = mach_corr_bem_data["WIND SPEED (m/s)"]
    mach_corr_bem_thrust = mach_corr_bem_data[" THRUST (N)"]
    mach_corr_bem_power = mach_corr_bem_data[" PRODUCED POWER (W)"]
    mach_corr_bem_cp = mach_corr_bem_data[" C_P"]

    stall_corr_bem_data = pd.read_csv(ROOT_DIR / "naca_data" / "bem_data" / "NREL_BEM_STALL_DELAY_CORR.csv")
    stall_corr_bem_wind_speed = stall_corr_bem_data["WIND SPEED (m/s)"]
    stall_corr_bem_thrust = stall_corr_bem_data[" THRUST (N)"]
    stall_corr_bem_power = stall_corr_bem_data[" PRODUCED POWER (W)"]
    stall_corr_bem_cp = stall_corr_bem_data[" C_P"]

    tiproot_corr_bem_data = pd.read_csv(ROOT_DIR / "naca_data" / "bem_data" / "NREL_BEM_TIPROOT_CORRECTION.csv")
    tiproot_corr_bem_wind_speed = tiproot_corr_bem_data["WIND SPEED (m/s)"]
    tiproot_corr_bem_thrust = tiproot_corr_bem_data[" THRUST (N)"]
    tiproot_corr_bem_power = tiproot_corr_bem_data[" PRODUCED POWER (W)"]
    tiproot_corr_bem_cp = tiproot_corr_bem_data[" C_P"]

    no_corr_bem_data = pd.read_csv(ROOT_DIR / "naca_data" / "bem_data" / "NREL_BEM_NO_CORRECTIONS.csv")
    no_corr_bem_wind_speed = no_corr_bem_data["WIND SPEED (m/s)"]
    no_corr_bem_thrust = no_corr_bem_data[" THRUST (N)"]
    no_corr_bem_power = no_corr_bem_data[" PRODUCED POWER (W)"]
    no_corr_bem_cp = no_corr_bem_data[" C_P"]




    qblade_wind_speed, qblade_cp = parse_qblade_file("cp_v_windspeed.txt")
    qblade_wind_speed, qblade_power = parse_qblade_file("power_v_windspeed.txt")
    qblade_wind_speed, qblade_thrust = parse_qblade_file("thrust_v_windspeed.txt")
    qblade_wind_speed = np.array(qblade_wind_speed)

    qblade_cp = np.array(qblade_cp)
    qblade_power = np.array(qblade_power)
    qblade_thrust = np.array(qblade_thrust)

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(nrel_wind_speed, nrel_cp, label='NREL $\\text{C}_{\\mathbf{P}}$', color=red_color, linestyle='-', linewidth=3)
    ax.plot(all_corr_bem_wind_speed, all_corr_bem_cp, label='BEM $\\text{C}_{\\mathbf{P}}$', color=blue_color, linestyle='-', linewidth=3)
    ax.plot(qblade_wind_speed, qblade_cp, label='QBlade $\\text{C}_{\\mathbf{P}}$', color=green_color, linestyle='-', linewidth=3)
    ax.plot(no_corr_bem_wind_speed, no_corr_bem_cp, label='BEM No Corrections $\\text{C}_{\\mathbf{P}}$', color=yellow_color, linestyle='--', linewidth=1)
    ax.plot(mach_corr_bem_wind_speed, mach_corr_bem_cp, label='BEM Mach Correction $\\text{C}_{\\mathbf{P}}$', color="blue", linestyle='--', linewidth=1)
    ax.plot(stall_corr_bem_wind_speed, stall_corr_bem_cp, label='BEM Stall Delay $\\text{C}_{\\mathbf{P}}$', color="purple", linestyle='--', linewidth=1)
    ax.plot(tiproot_corr_bem_wind_speed, tiproot_corr_bem_cp, label='BEM Tip/Root $\\text{C}_{\\mathbf{P}}$', color="green", linestyle='--', linewidth=1)

    # ax.plot(nrel_wind_speed, np.interp(nrel_wind_speed, bem_wind_speed, bem_cp),
    #     linestyle='None', marker='o', color=blue_color, markeredgecolor='black', markersize=5)
    ax.plot(nrel_wind_speed, np.interp(nrel_wind_speed, qblade_wind_speed, qblade_cp),
            linestyle='None', marker='o', color=green_color, markeredgecolor='black', markersize=5)
    ax.plot(nrel_wind_speed, nrel_cp, linestyle='None', marker='o', color=red_color, markeredgecolor='black', markersize=5)


    ax.set_xlim(3.0, 25.0)
    ax.set_xlabel('Wind Speed (m/s)', fontsize=20)
    ax.set_ylabel('$\\mathbf{C}_{\\mathbf{P}}$', fontsize=20)
    ax.set_title('$\\mathbf{C}_{\\mathbf{P}}$ vs Wind Speed', fontsize=20)
    ax.legend(loc='lower left', frameon=True, framealpha=0.7, edgecolor='black', fontsize=16)

    ax.set_xlim(min(nrel_wind_speed), max(nrel_wind_speed))
    ax.set_ylim(0, 0.59)
    apply_grid_styling(ax)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    plt.tight_layout()
    plt.savefig("Cp_vs_WindSpeed.png", dpi=300)
    plt.show()


    fig, ax = plt.subplots(figsize=(10, 6.5))
    bem_cp_interp = np.interp(nrel_wind_speed, bem_wind_speed, bem_cp)
    qblade_cp_interp = np.interp(nrel_wind_speed, qblade_wind_speed, qblade_cp)
    ax.plot(nrel_wind_speed, 100 * np.abs(bem_cp_interp - nrel_cp), label='BEM $\\text{C}_{\\mathbf{P}}$ - NREL $\\text{C}_{\\mathbf{P}}$', color=blue_color, linestyle='-', linewidth=3)
    ax.plot(nrel_wind_speed, 100 * np.abs(qblade_cp_interp - nrel_cp), label='QBlade $\\text{C}_{\\mathbf{P}}$ - NREL $\\text{C}_{\\mathbf{P}}$', color=green_color, linestyle='-', linewidth=3)
    ax.set_xlim(7.5, 11.25)
    ax.set_xlabel('Wind Speed (m/s)', fontsize=20)
    ax.set_ylabel("$\\mathbf{C}_{\\mathbf{P}}$ Residual ($\\mathbf{\\%}$)", fontsize=20)
    ax.set_title('$\\mathbf{C}_{\\mathbf{P}}$ Error vs Wind Speed', fontsize=20)
    ax.legend(loc='upper right', frameon=True, framealpha=0.7, edgecolor='black', fontsize=20)
    ax.set_ylim(0, 10)
    apply_grid_styling(ax)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    plt.tight_layout()
    plt.show()






    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(nrel_wind_speed, nrel_power, label='NREL Power', color=red_color, linestyle='-', linewidth=3)
    ax.plot(bem_wind_speed, bem_power/1000, label='BEM Power', color=blue_color, linestyle='-', linewidth=3)
    ax.plot(qblade_wind_speed, qblade_power/1000, label='QBlade Power', color=green_color, linestyle='-', linewidth=3)

    ax.plot(nrel_wind_speed, np.interp(nrel_wind_speed, bem_wind_speed, bem_power/1000),
            linestyle='None', marker='o', color=blue_color, markeredgecolor='black', markersize=5)
    ax.plot(nrel_wind_speed, np.interp(nrel_wind_speed, qblade_wind_speed, qblade_power/1000),
            linestyle='None', marker='o', color=green_color, markeredgecolor='black', markersize=5)
    ax.plot(nrel_wind_speed, nrel_power, linestyle='None', marker='o', color=red_color, markeredgecolor='black', markersize=5)

    ax.set_xlim(min(nrel_wind_speed), max(nrel_wind_speed))
    ax.set_xlabel('Wind Speed (m/s)', fontsize=20)
    ax.set_ylabel('Power (kW)', fontsize=20)
    ax.set_title('Power vs Wind Speed', fontsize=20)
    ax.legend(loc='upper right', frameon=True, framealpha=0.7, edgecolor='black', fontsize=20)
    # ax.set_xlim(min(nrel_wind_speed), max(nrel_wind_speed))
    # ax.set_ylim(0, 1.2 * max(nrel_power))
    apply_grid_styling(ax)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    plt.tight_layout()
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(nrel_wind_speed, nrel_thrust, label='NREL Thrust', color=red_color, linestyle='-', linewidth=3)
    ax.plot(bem_wind_speed, bem_thrust/1000, label='BEM Thrust', color=blue_color, linestyle='-', linewidth=3)
    ax.plot(qblade_wind_speed, qblade_thrust/1000, label='QBlade Thrust', color=green_color, linestyle='-', linewidth=3)

    ax.plot(nrel_wind_speed, np.interp(nrel_wind_speed, bem_wind_speed, bem_thrust/1000),
            linestyle='None', marker='o', color=blue_color, markeredgecolor='black', markersize=5)
    ax.plot(nrel_wind_speed, np.interp(nrel_wind_speed, qblade_wind_speed, qblade_thrust/1000),
            linestyle='None', marker='o', color=green_color, markeredgecolor='black', markersize=5)
    ax.plot(nrel_wind_speed, nrel_thrust, linestyle='None', marker='o', color=red_color, markeredgecolor='black', markersize=5)



    ax.set_xlabel('Wind Speed (m/s)', fontsize=20)
    ax.set_ylabel('Thrust (kN)', fontsize=20)
    ax.set_title('Thrust vs Wind Speed', fontsize=20)
    ax.legend(loc='upper right', frameon=True, framealpha=0.7, edgecolor='black', fontsize=20)
    # ax.set_xlim(min(nrel_wind_speed), max(nrel_wind_speed))
    # ax.set_ylim(0, 1.2 * max(nrel_power))
    apply_grid_styling(ax)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    plt.tight_layout()
    plt.show()





# bem_data = get_bem_data("blade_profile_test.json")

# def plot_





R = np.linspace(HUB_RADIUS, TIP_RADIUS, 100)
CHORD = np.array([CHORD_DISTRIBUTION(r) for r in R])
BETA = np.array([TWIST_DISTRIBUTION(r) for r in R])


with open(ROOT_DIR / "naca_data" / "blade_profiles" / "blade_profile_test.json", 'r') as f:
    data = json.load(f)


plot_blade_with_distribution_info((data, R, CHORD))
plot_chord_distribution((data, R, CHORD))
plot_twist_distribution((data, R, BETA))

plot_bem_verification()