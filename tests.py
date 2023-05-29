import quantum_annealing_simulator as aqs
import numpy as np
from hypothesis import strategies as st
from hypothesis import settings
from hypothesis import given

@given(n = st.integers(1,5), m = st.integers(1,5))
@settings(max_examples = 1)

def test_lattice(n, m):
    spins = {}
    sites = []
    spins = aqs.lattice(n, m, sites, spins)
    assert len(spins) == n*m
    abs_spins = np.abs(list(spins.values()))
    assert abs_spins.all() == 1

@given(n = st.integers(1,5), m = st.integers(1,5))
@settings(max_examples = 1)

def test_configurations(n, m):
    configurations = []
    aqs.possible_configurations(n, m, configurations)
    assert len(configurations) == 2**(n*m)