import functions as func
import numpy as np
from hypothesis import strategies as st
from hypothesis import settings
from hypothesis import given

@given(n = st.integers(1,4), m = st.integers(1,3))
@settings(max_examples = 1)

def test_lattice(n, m):
    
    """ Test to verify that the generated lattice consists of spins with values
        of +1 or -1 with appropiate dimensions."""
        
    sites, spins = func.lattice(n, m)
    assert len(spins) == n*m
    abs_spins = np.abs(list(spins.values()))
    assert abs_spins.all() == 1

@given(n = st.integers(1,4), m = st.integers(1,3))
@settings(max_examples = 1)

def test_configurations(n, m):
    
    """ Test to verify that the quantity of configurations
        is the same as the base of states."""
        
    configurations = func.possible_configurations(n, m)
    assert len(configurations) == 2**(n*m)
       
@given(n = st.integers(1,4), m = st.integers(1,3))
@settings(max_examples = 1)

def test_matrixH_0(n, m):
    
    """ Test to verify that the hamiltonian H_0
        affects only to exactly the particles of the lattice."""
        
    configurations = func.possible_configurations(n, m)
    H = func.H_0(n, m, configurations)
    for i in range(2**(m*n)):
        assert np.sum(H[i]) == -n*m
        
@given(n = st.integers(1, 4), m = st.integers(1, 3), J = st.floats(-1, 1), h = st.floats(-1, 1))
@settings(max_examples = 1)

def test_matrixH_1(n, m, J, h):
    sites, spins = func.lattice(n, m)
    configurations = func.possible_configurations(n, m)
    neighbors = func.neighbors_classification(sites)
    H = func.H_1(n, m, spins, configurations, neighbors, J, h)
    for i in range(2**(m*n)):
        for j in range(2**(m*n)):
            assert H[i, j] == 0
        
@given(n = st.integers(1,4), m = st.integers(1,3))
@settings(max_examples = 1)      
        
def test_magnetization(n, m):
    configurations = func.possible_configurations(n, m)
    for configuration in configurations:
        mag = func.magnetization(configuration)
        assert -1 <= mag <= 1
