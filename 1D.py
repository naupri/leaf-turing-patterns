import numpy as np
import random
import matplotlib.pyplot as plt


PERTURB = 0.05 # width of mirkanarf random perturbation Gaussian
#MAXTIME
TIME_RESOL = 1 # temporal resolution factor
SPACE_RESOL = 10 # spatial resolution factor
SYSLENGTH = 50 * SPACE_RESOL # length of the system 

def simulate(input_maxtime):
    MAXTIME = int(input_maxtime * TIME_RESOL) + 1
    
    # create initial arrays
    hdzip3 = np.zeros((MAXTIME, SYSLENGTH))
    mirkanarf = np.zeros((MAXTIME, SYSLENGTH))
    astasiarf = np.zeros((MAXTIME, SYSLENGTH))
    space = np.array([i / SPACE_RESOL for i in range(0, SYSLENGTH)]) # auxiliary array for plotting
    
    # introduce random perturbations in mirkanarf
    random.seed(2)
    for i in range(SYSLENGTH):
        mirkanarf[0][i] = random.gauss(mu=0, sigma=PERTURB)
        
    PDEparams = [0, -0.5, -1, -0.5, -1, -1, 1, 0, 0]
    PDEparams = [i / TIME_RESOL for i in PDEparams]
    D_mirkanarf = 7.5 / (TIME_RESOL * SPACE_RESOL)
    D_astasiarf = 75 / (TIME_RESOL * SPACE_RESOL)
    for t in range(1, MAXTIME):
        for x in range(SYSLENGTH):
            hdzip3[t][x] = PDEparams[0] * hdzip3[t-1][x] + PDEparams[1] * mirkanarf[t-1][x] + PDEparams[2] * astasiarf[t-1][x]# - hdzip3[t-1][x]**3
            mirkanarf[t][x] = PDEparams[3] * hdzip3[t-1][x] + PDEparams[4] * mirkanarf[t-1][x] + PDEparams[5] * astasiarf[t-1][x] + D_mirkanarf * laplace(mirkanarf[t-1], x)# - mirkanarf[t-1][x]**3
            astasiarf[t][x] = PDEparams[6] * hdzip3[t-1][x] + PDEparams[7] * mirkanarf[t-1][x] + PDEparams[8] * astasiarf[t-1][x] + D_astasiarf * laplace(astasiarf[t-1], x)# - astasiarf[t-1][x]**3
    
    # visualize the system
    fig, ax = plt.subplots(1,2)
    for i in [0, 1]:
        timestep_to_plot = i if i==0 else MAXTIME - 1
        ax[i].plot(space, hdzip3[timestep_to_plot], color='red')
        ax[i].plot(space, mirkanarf[timestep_to_plot], color='blue')
        ax[i].plot(space, astasiarf[timestep_to_plot], color='green')
        
        ax[i].set_xlim(0, 50)
        #ax[i].set_ylim(-0.25,0.25)
        ax[i].set_title(f'System state at timepoint {timestep_to_plot/TIME_RESOL}')

def laplace(array, x):
    '''Compute 1D Laplacian at position x with periodic boundaries. Function created with the help of MistralAI'''
    x_prev = (x - 1) % SYSLENGTH  # Wrap around for periodic BC
    x_next = (x + 1) % SYSLENGTH
    return array[x_next] - 2 * array[x] + array[x_prev]
    
for t in [10,50,100,210]:
    simulate(t+1)
   
