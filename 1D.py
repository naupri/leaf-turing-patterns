import numpy as np
import random
import matplotlib.pyplot as plt

MAXTIME = 1000 # timesteps to simulate
SYSLENGTH = 50 # length of the system
PERTURB = 0.05 # width of mirkanarf random perturbation Gaussian

# create initial arrays
hdzip3 = np.zeros((MAXTIME, SYSLENGTH))
mirkanarf = np.zeros((MAXTIME, SYSLENGTH))
astasiarf = np.zeros((MAXTIME, SYSLENGTH))
space = np.arange(0,SYSLENGTH,1) # auxiliary array for plotting

# introduce random perturbations in mirkanarf
random.seed(2)
for i in range(SYSLENGTH):
    mirkanarf[0][i] = random.gauss(mu=0, sigma=PERTURB)

# visualize the system
fig, ax = plt.subplots(1,2)
for i in [0, 1]:
    timestep_to_plot = i if i==0 else MAXTIME - 1
    ax[i].plot(space, hdzip3[timestep_to_plot], color='red')
    ax[i].plot(space, mirkanarf[timestep_to_plot], color='blue')
    ax[i].plot(space, astasiarf[timestep_to_plot], color='green')
    
    ax[i].set_xlim(0, 50)
    ax[i].set_ylim(-0.25,0.25)
    ax[i].set_title(f'System state at timepoint {timestep_to_plot}')



