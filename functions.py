import itertools
import random
from scipy.sparse import lil_matrix
import numpy as np


# List with coordinates of each spin:

def lattice(n, m):
    """This method calculates the energy of a given configuration, given the exchange constant J = 1.
        
    Parameters:
        config: state of the configuration created by initialstate(N,M).
            
    Returns:
        The calculated energy of the modified configuration of the lattice."""
    sites = []
    spins = {}
    for i,j in itertools.product(range(n), range(m)):
        sites.append((i,j))
        spins[(i,j)] = random.choice([-1,1])
        
    return sites, spins
	
# Value-assignment to each spin:	

def neighbors_classification(sites, spins):
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
    configurations = []
    for k in itertools.product([1,-1] , repeat=n*m):
    	configurations.append(list(k))
    return configurations

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

def H_1(n, m, spins, configurations, neighbors, J, h):
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
    mag = np.sum(configuration)/len(configuration)
    return mag