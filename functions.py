import time
import itertools
import random
from scipy.sparse import lil_matrix
import numpy as np
from scipy.sparse.linalg import eigsh


def lattice(n, m):
    
    """This method generates the lattice of the system
        
    Parameters:
        n: number of particles in the horizontal direction.
        m: number of particles in the vertical direction.
            
    Returns:
        sites: coordinates of the particles.
        spins: dictionary with the spin-value of each coordinate."""
        
    sites = []
    spins = {}
    for i,j in itertools.product(range(n), range(m)):
        sites.append((i,j))
        spins[(i,j)] = random.choice([-1,1])
        
    return sites, spins
	
# Value-assignment to each spin:	

def neighbors_classification(sites):
    
    """This method classifies the neighbors of each particle
        
    Parameter:
        sites: coordinates of the particles.
            
    Returns:
        neighbors: dictionary with the neighbors of each particle."""
        
    neighbors = {}
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
    return neighbors
            
# Creation of all possible configurations of spins:

def possible_configurations(n, m):
    
    """This method generates all the possible configurations of the spins.
        
    Parameters:
        n: number of particles in the horizontal direction.
        m: number of particles in the vertical direction.
            
    Returns:
        configurations: list with all the configurations of the lattice."""
        
    configurations = []
    for k in itertools.product([1,-1] , repeat=n*m):
    	configurations.append(list(k))
    return configurations

# Hamiltonians

    # Disperse matrix H_0
   
def H_0(n, m, configurations):
    
    """This method generates the starting hamiltonian H_0.
        
    Parameters:
        n: number of particles in the horizontal direction.
        m: number of particles in the vertical direction.
        configurations: list with all the configurations of the lattice.
            
    Returns:
        H_0: starting hamiltonian."""
        
    p = m*n
    H_0 = lil_matrix((2**p, 2**p))
    for i in range(2**p):
        for j in range(1, p+1):
            idx = int(abs(i+1+configurations[i][j-1]*2**(p-j)))-1
            H_0[idx, i] = -1
    return H_0

    # Disperse matrix H_1

def H_1(n, m, spins, configurations, neighbors, J, h):
    
    """This method generates the final hamiltonian H_1 according to the
    parameters of the system.
        
    Parameters:
        n: number of particles in the horizontal direction.
        m: number of particles in the vertical direction.
        configurations: list with all the configurations of the lattice.
        neighbors: dictionary with the neighbors of each particle.
        J : coupling strength
        h : magnetic field
    Returns:
        H_1: starting hamiltonian."""
        
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


def magnetization(configuration):
    
    """This method calculates the magnetization of a given configuration
    of the lattice.
        
    Parameter:
        configurations: list with all the configurations of the lattice.

    Returns:
        mag: magnetization of the configuration."""
        
    mag = np.sum(configuration)/len(configuration)
    return mag

def simulation(J, h, step, n, m):
    
    # INITIALIZATION OF THE SYSTEM
    
    p = m*n
    param_lambda = np.arange(0, 1+step, step)
    amplitudes_values = np.zeros((len(param_lambda), 2**p))
    magnetization_values = np.zeros(len(param_lambda))

    configurations = possible_configurations(n, m)
    sites, spins = lattice(n, m)
    neighbors = neighbors_classification(sites)

    H0 = H_0(n, m, configurations)
    H1 = H_1(n, m, spins, configurations, neighbors, J, h)
    
    start = time.time()
    
    # SIMULATION
    
    l = 0

    for t in param_lambda:
        
        mag = 0.0
        H = (1-t)*H0 + t*H1
        energy, ground_state = eigsh(H, 1, which='SA')
        
        for i in range(2**p):
            amplitudes_values[l,i] = abs(ground_state[i,0])**2
            ind_mag = magnetization(configurations[i])
            mag += abs((ground_state[i][0].real))**2*ind_mag
        magnetization_values[l] = mag
        
        l += 1

    b = list(abs(ground_state[:,0])).index(max(abs(ground_state[:,0])))   # Position of the list 'configurations' where the solution is
    
    final = time.time()

    print('Spent time for the simulation: ' , final - start , ' seconds')
    
    print('Energy of the configuration:' , energy)    # Minimum eigenvalue of H
    return configurations, sites, b, param_lambda, magnetization_values, amplitudes_values
