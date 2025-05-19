import math

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt



hub = 0.02
tip = 0.227
tsr = 3.5

B = 2 # number of blades
a_des = math.radians(2.25) # design angle of attack
cl_des = 1.012 # design lift coefficient


r = np.linspace(hub, tip, 100) # radial position

phi = 2.0 / 3.0 * np.arctan(1.0 / (tsr * r/tip)) # inflow angle
beta = phi - a_des # twist angle

chord = (8 * np.pi * r / (B * cl_des) * (1 - np.cos(phi)))/tip # chord length


plt.figure(figsize=(10, 6))
plt.plot(r/tip, chord, label='Chord Length', color='blue')
# plt.plot(r, beta, label='Twist Angle', color='red')
plt.xlabel('$\\frac{r}{R}$', fontsize=16)
plt.ylabel('$\\frac{c}{R}$', fontsize=16)
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.title('Blade Chord Length Distribution')
plt.legend()
plt.grid()
plt.show()
plt.savefig('blade_chord_length_distribution.svg', dpi=300)


plt.figure(figsize=(10, 6))
plt.plot(r/tip, beta, label='Twist Angle', color='red')
plt.xlabel('$\\frac{r}{R}$', fontsize=16)
plt.ylabel('$\\beta$  (rad)', fontsize=16)
plt.xlim(0, 1)
# plt.ylim(-0.5, 0.5)
plt.title('Blade Twist Angle Distribution')
plt.legend()
plt.grid()
plt.savefig('blade_twist_angle_distribution.svg', dpi=300)
plt.show()






