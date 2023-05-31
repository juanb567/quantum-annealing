import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh
import functions as func


############################# VARIABLES ##################################


J = 1
h = 0.0001
step = 0.01
n = 3
m = 3

directory = './Graphs'  # Directory where the plots will be saved.                 

assert n*m <= 16   # The problem is too big.

start = time.time()

p = m*n
param_lambda = np.arange(0, 1+step, step)
amplitudes_values = np.zeros((len(param_lambda), 2**p))
magnetization_values = np.zeros(len(param_lambda))

configurations = func.possible_configurations(n, m)
sites, spins = func.lattice(n, m)
neighbors = func.neighbors_classification(sites)


H_0 = func.H_0(n, m, configurations)
H_1 = func.H_1(n, m, spins, configurations, neighbors, J, h)

print('Matrices created')


############################# SIMULATION ##################################


l = 0

for t in param_lambda:
    
    mag = 0.0
    H = (1-t)*H_0 + t*H_1
    energy, ground_state = eigsh(H, 1, which='SA')
    
    for i in range(2**p):
        amplitudes_values[l,i] = abs(ground_state[i,0])**2
        ind_mag = func.magnetization(configurations[i])
        mag += abs((ground_state[i][0].real))**2*ind_mag
    magnetization_values[l] = mag
    
    l += 1

b = list(abs(ground_state[:,0])).index(max(abs(ground_state[:,0])))    # Position of the list 'configurations' where the solution is

print('Energy of the configuration:' , energy)    # Minimum eigenvalue of H

final = time.time()


############################# GRAPHICS ##################################

def save_plot(plot_function, filename, directory=None):
    
    """This module saves the plot in the desired directory."""
    
    if directory is not None:
        os.makedirs(directory, exist_ok=True)
        filepath = os.path.join(directory, filename)
    else:
        filepath = filename

    plot_function()
    plt.savefig(filepath)
    plt.close()


def plot_configuration():
    
    """This module plots the configuration of the ground state
    at the end of the process."""
    
    plt.figure()
    colores = {1: "green", -1: "red"}
    for i in range(len(configurations[b])):
        x , y = sites[i]
        plt.quiver(x, y, 0, configurations[b][i], pivot="middle", color=colores[configurations[b][i]],
                   linewidths=0.3, scale=10)
    plt.xticks(range(-1, n+1))
    plt.yticks(range(-1, m+1))
    plt.gca().set_aspect("equal")
    plt.grid(visible=False)


def plot_magnetization():
    
    """This module plots the magnetization of the 
    system during the process."""
    
    plt.figure()
    plt.plot(param_lambda, magnetization_values)
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Magnetization$")
    plt.xlim([0, 1])
    plt.grid()

def plot_amplitudes():
    
    """This module plots the amplitude of the each configuration
    during the process."""
    
    plt.figure()
    for i in range(2**p):
        plt.plot(param_lambda, amplitudes_values[:,i])
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Amplitudes$")
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.grid()
    
def plot_amplitudes_log():
    
    """This module plots the amplitude of the each configuration
    during the process in a logarithmic scale."""
    
    plt.figure()
    plt.yscale(value='log')
    for i in range(2**p):
        plt.plot(param_lambda, amplitudes_values[:,i])
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Amplitudes$")
    plt.xlim([0, 1])
    plt.ylim([1e-5, 1])
    plt.grid()

########################################################################

save_plot(plot_configuration, "configuration.png", directory)
save_plot(plot_magnetization, "magnetization.png", directory)
save_plot(plot_amplitudes, "amplitudes.png", directory)
save_plot(plot_amplitudes_log, "amplitudes_log.png", directory)

print('Spent time for the simulation: ' , final - start , ' seconds')
