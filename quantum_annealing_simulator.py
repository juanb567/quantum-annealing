import itertools
import random
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh

start = time.time()

#System variables

J = 1
h = 0.001

n = 3
m = 3                   #Number of spins m*n
p = n*m
neighbors = {}
spins = {}
sites = []
configurations = []

# List with coordinates of each spin:

def lattice(n, m, sites, spins):
    for i,j in itertools.product(range(n), range(m)):
        sites.append((i,j))
    for spin in sites:
        spins[spin]=random.choice([-1,1])
    return spins
	
# Value-assignment to each spin:	

def neighbors_classification(spins, neighbors):
    for site in sites:
    	i, j = site
    	neighbors[site]=[]
    	for neighbor_spin in sites:
    		k, l = neighbor_spin
    		if i+1==k and j==l:
    			neighbors[site].append(neighbor_spin)
    		if i-1==k and j==l:
    			neighbors[site].append(neighbor_spin)
    		if j+1==l and i==k:
    			neighbors[site].append(neighbor_spin)
    		if j-1==l and i==k:
    			neighbors[site].append(neighbor_spin)
            
# Creation of all possible configurations of spins:

def possible_configurations(n, m, configurations):
    for k in itertools.product([1,-1] , repeat=n*m):
    	configurations.append(list(k))
    

lattice(n, m, sites, spins)
neighbors_classification(spins, neighbors)
possible_configurations(n, m, configurations)


# Hamiltonians

    # Disperse matrix H_0
   
def H_0(n, m, configurations):
    p = m*n
    H_0 = lil_matrix((2**p, 2**p))
    for i in range(2**p):
        for j in range(1, p+1):
            idx = int(abs(i+1+configurations[i][j-1]*2**(p-j)))-1
            H_0[idx, i] = -1
    return H_0

    # Disperse matrix H_1

def H_1(n, m, spins, configurations, J, h):
    p = m*n
    H_1 = lil_matrix((2**p, 2**p))
    for i in range(2**p):
        j = -1
        for spin in spins:
            spins[spin] = configurations[i][j+1]
            j += 1
        for spin in spins:
            H_1[i, i] += -h*spins[spin]
            for neighbor in neighbors[spin]:
                H_1[i, i] += -0.5*J*spins[spin]*spins[neighbor]
    return H_1

H_0 = H_0(n, m, configurations)
H_1 = H_1(n, m, spins, configurations, J, h)

print('Matrices created')

def magnetization(configuration):
    mag = np.sum(configuration)/len(configuration)
    return mag

# Adiabatic process

step = 0.01
param_lambda = np.arange(0, 1+step, step)

amplitudes_values = np.zeros((len(param_lambda), 2**p))
magnetization_values = np.zeros(len(param_lambda))

l = 0

for t in param_lambda:
    
    print(t*100 , '%')
    
    mag = 0.0
    H = (1-t)*H_0 + t*H_1
    energy, ground_state = eigsh(H, 1, which='SA')
    
    for i in range(2**p):
        amplitudes_values[l,i] = abs(ground_state[i,0])**2
        ind_mag = magnetization(configurations[i])
        mag += abs((ground_state[i][0].real))**2*ind_mag
    magnetization_values[l] = mag
    
    l += 1

b = list(abs(ground_state[:,0])).index(max(abs(ground_state[:,0])))    # Posición de la lista 'configuraciones' donde se encuentra la solución al problema

############################# GRAPHICS ##################################

def plot_configuration():
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
    plt.show()

def plot_magnetization():
    plt.figure()
    plt.plot(param_lambda, magnetization_values)
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Magnetization$")
    plt.xlim([0, 1])
    plt.grid()
    plt.show()    

def plot_amplitudes():
    plt.figure()
    for i in range(2**p):
        plt.plot(param_lambda, amplitudes_values[:,i])
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Amplitudes$")
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.grid()
    plt.show()
    
def plot_amplitudes_log():
    plt.figure()
    plt.yscale(value='log')
    for i in range(2**p):
        plt.plot(param_lambda, amplitudes_values[:,i])
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Amplitudes$")
    plt.xlim([0, 1])
    plt.ylim([1e-5, 1])
    plt.grid()
    plt.show()

########################################################################

plot_configuration()
plot_magnetization()
plot_amplitudes()
plot_amplitudes_log()


final = time.time()

print('Spent time: ' , final - start , ' seconds')