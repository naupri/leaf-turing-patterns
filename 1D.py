import math, random
import matplotlib.pyplot as plt
import numpy as np

### Simulation parameters ###
MAXTIME = 5 # run the simulation until this time
TIMESTEP = 0.1 # size of the timestep
SYSLENGTH = 50 # length of the system to simulate
SPACESTEP = 1 # size of the spatial step

PERTURB = 0.05 # sd of the Gaussian from which to sample the initial perturbation of mirkanarf
PDEparams = [0, -0.5, -1, -0.5, -1, -1, 1, 0, 0] # k-parameters of the PDE
D_mirkanarf = 7.5 # diffusion coefficient for miRKANARF
D_astasiarf = 75 # diffusion coefficient for AStasiARF

SEED = 2 # Seed for the random number generator
SAMPLE_TIMES = [i for i in range(5)] # the timesteps to plot
TIMECOUNT = math.ceil(MAXTIME / TIMESTEP) + 1 # number of timesteps to compute
SPACECOUNT = math.ceil(SYSLENGTH / SPACESTEP) + 1 # number of spatial coordinates to compute


def main():
    # create initial arrays
    hdzip3 = np.zeros((TIMECOUNT, SPACECOUNT))
    mirkanarf = np.copy(hdzip3)
    astasiarf = np.copy(hdzip3)
    space = np.arange(0, SYSLENGTH + SPACESTEP, SPACESTEP) # auxiliary array for plotting
    
    # introduce random perturbations in mirkanarf
    random.seed(SEED)
    for i in range(SPACECOUNT):
        mirkanarf[0][i] = random.gauss(mu=0, sigma=PERTURB)
    
    # run the simulation
    for t in range(1, TIMECOUNT):
        for x in range(SPACECOUNT):
            hdzip3[t][x] = hdzip3[t-1][x] + TIMESTEP * hdzip3_step(t-1, x, hdzip3, mirkanarf, astasiarf)
            mirkanarf[t][x] = mirkanarf[t-1][x] + TIMESTEP * mirkanarf_step(t-1, x, hdzip3, mirkanarf, astasiarf)
            astasiarf[t][x] = astasiarf[t-1][x] + TIMESTEP * astasiarf_step(t-1, x, hdzip3, mirkanarf, astasiarf)
            
        # generate plots at selected timepoints
        if t * TIMESTEP in SAMPLE_TIMES:
            fig, ax = plt.subplots(1, 2)
            for i in [0, 1]:
                timestep_to_plot = 0 if i == 0 else t
                
                ax[i].plot(space, mirkanarf[timestep_to_plot], color='blue')
                ax[i].plot(space, astasiarf[timestep_to_plot], color='green')
                ax[i].plot(space, hdzip3[timestep_to_plot], color='red')
                
                ax[i].set_xlim(0, SYSLENGTH + SPACESTEP)
                ax[i].set_title(f'System state at timepoint {timestep_to_plot * TIMESTEP}')

def hdzip3_step(t, x, hdzip3, mirkanarf, astasiarf):
    '''Compute the change in HDZIPIII given the current state of the system'''
    return PDEparams[0] * hdzip3[t][x] + PDEparams[1] * mirkanarf[t][x] + PDEparams[2] * astasiarf[t][x]

def mirkanarf_step(t, x, hdzip3, mirkanarf, astasiarf):
    '''Compute the change in miRKANARF given the current state of the system'''
    return PDEparams[3] * hdzip3[t][x] + PDEparams[4] * mirkanarf[t][x] + PDEparams[5] * astasiarf[t][x] + D_mirkanarf * laplace(mirkanarf[t], x)

def astasiarf_step(t, x, hdzip3, mirkanarf, astasiarf):
    '''Compute the change in AStasiARF given the current state of the system'''
    return PDEparams[6] * hdzip3[t][x] + PDEparams[7] * mirkanarf[t][x] + PDEparams[8] * astasiarf[t][x] + D_astasiarf * laplace(astasiarf[t], x)

def laplace(array, x):
    '''Compute 1D Laplacian at position x with periodic boundaries. Function created with the help of MistralAI'''
    x_prev = (x - 1) % SPACECOUNT  # Wrap around for periodic BC
    x_next = (x + 1) % SPACECOUNT
    return (array[x_next] - 2 * array[x] + array[x_prev]) / (SPACESTEP ** 2)
   
main()
   
