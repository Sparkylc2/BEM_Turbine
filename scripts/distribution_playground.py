import numpy as np
from pathlib import Path
import matplotlib as plt



hub = 0.02
tip = 0.227
tsr = 3.5

B = 2 # number of blades
a_des = 2.25 # design angle of attack
cl_des = 1.012 # design lift coefficient


r = np.linspace(hub, tip, 100) # radial position

phi = 2.0 / 3.0 * np.arctan(1.0 / (tsr * r/tip)) # inflow angle
beta = phi - a_des # twist angle

chord = 8 * np.pi * r / (B * cl_des) * (1 - np.cos(phi)) # chord length


# Plotting the results
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.plot(r/tip, chord, label='Chord Length', color='blue')
# plt.plot(r, beta, label='Twist Angle', color='red')
plt.xlabel('r/R (non-dimensional)')
plt.ylabel('Chord Length (m)')
plt.xlim(0, 1)
plt.title('Blade Chord Length Distribution')
plt.legend()
plt.grid()
plt.show()





