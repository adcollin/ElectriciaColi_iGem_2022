# loading python packages
import numpy as np
from scipy import integrate as sp
import matplotlib.pyplot as plt

# General carateristics
time_span = 100
time_points = 300

# Constantes definition TO BE CHANGED FOLLOWING BIBLIOGRAPHY
# Intermediates
k_cat = 2.68 #[µM/s]
Km = 45.0 #[µM]
#E0 = 10.0

# Mediator reduction
k_2 = 0.722 #[1/s]

# Mediator oxidation
k_3 = 0.15  #[1/s]
k_4 = 28.95 #[mA]

# ODE function

def ODE(x, t):
    # Stoichiometric matrix
    #                NADH,H+   NAD+    M_ox    M_red      I
    SM = np.matrix([
                   [1,      -1,         0,       0,        0 ],  # intermediates
                   [-1,      1,        -1,       1,        0 ],  # mediator reduction
                   [0,      0,          1,      -1,        0 ],  # mediator oxidation
                   [0,      0,          0,       0,        1 ],  # current conversion
                   [0,      0,          0,       0,       -1 ],  # already existing electic flux
                   ])

    # Rate matrix
    RM = np.matrix([
                   [k_cat*(x[1]/(Km+x[1]))],  # intermediates raye
                   [k_2*x[0]*x[2]],  # mediator reduction rate
                   [k_3*x[3]],  # mediator oxidation rate
                   [k_4*x[3]], # current conversion from new electronic input
                   [x[4]] # already existing electic flux
                   ])

    d_dt = np.transpose(SM)*RM
    return([d_dt[0, 0], d_dt[1, 0], d_dt[2, 0], d_dt[3, 0], d_dt[4, 0]])


# ODE solution using scipy package
# Time distribution to simulate
t = np.linspace(0, time_span, time_points)

# Initial conditions (state variable values at t=0)
NADH_0 = 0.0
NAD_0 = 50.0
MB_ox_0 = 0.0
MB_red_0 = 50.0
I_0 = 0.0

x0 = [NADH_0, NAD_0, MB_ox_0, MB_red_0, I_0]
out = sp.odeint(ODE, x0, t, args=())

# Data visualization
plt.plot(t, out[:, 0], label = '[NADH] (M)')
plt.plot(t, out[:, 1], 'tab:green', label ='[NAD](M)')
plt.plot(t, out[:, 2], 'tab:orange', label = '[MB_ox](M)')
plt.plot(t, out[:, 3], 'tab:pink', label = '[MB_red] (M)')
plt.plot(t, out[:, 4], 'tab:red', label = 'I (A)')
plt.xlabel('Time (s)')
plt.legend()
# fig, axs = plt.subplots(5, 1)
# axs[0].plot(t, out[:, 0])
# axs[0].set_xlabel('Time (s)')
# axs[0].set_ylabel('[NADH] (M)')
# axs[1].plot(t, out[:, 1], 'tab:green')
# axs[1].set_xlabel('Time (s)')
# axs[1].set_ylabel('[NAD](M)')
# axs[2].plot(t, out[:, 2], 'tab:orange')
# axs[2].set_xlabel('Time (s)')
# axs[2].set_ylabel('[MB_ox](M)')
# axs[3].plot(t, out[:, 3], 'tab:pink')
# axs[3].set_xlabel('Time (s)')
# axs[3].set_ylabel('[MB_red] (M)')
# axs[4].plot(t, out[:, 4], 'tab:red')
# axs[4].set_xlabel('Time (s)')
# axs[4].set_ylabel('I (A)')
plt.show()
