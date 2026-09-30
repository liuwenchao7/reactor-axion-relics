import numpy as np
from scipy.integrate import quad

I0 = 70.7  # m-2s-1sr-1

n = 3.01
E0 = 4.29
epsilon = 854.
Rd = 174.
theta_max = 0.616327

E_low = 0.105658
E_up = 1000.

def thetaSpectrum(theta):
    Dtheta = np.sqrt(Rd * Rd * np.cos(theta) * np.cos(theta) + 2 * Rd + 1) - Rd * np.cos(theta)
    Dtheta_n = np.power(Dtheta, -(n - 1)) * np.sin(theta)
    return Dtheta_n

def energySpectrum(energy):
    IE = np.power(E0 + energy, -n) / (1 + energy / epsilon)
    return IE

N = 1 / quad(energySpectrum, E_low, E_up)[0]

def muon_rate():
    return I0 * N * quad(lambda theta: thetaSpectrum(theta), 0, np.pi / 2)[0] * quad(energySpectrum, E_low, E_up)[0]
