import numpy as np
import random
import matplotlib.pyplot as plt

MAXTIME = 1000 # timesteps to simulate
SYSLENGTH = 50 # length of the system

# create initial arrays
hdzip3 = np.zeros((MAXTIME, SYSLENGTH))
mirkanarf = np.zeros((MAXTIME, SYSLENGTH))
astasiarf = np.zeros((MAXTIME, SYSLENGTH))

# introduce random perturbations in mirkanarf
random.seed(2)
for i in range(SYSLENGTH):
    mirkanarf[0][i] = random.gauss(mu=0, sigma=0.05)
