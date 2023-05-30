# Quantum annealing simulator for the 2D Ising model

Quantum annealing is a computational approach that leverages principles of quantum physics to solve optimization problems. It involves using a quantum annealer, a specialized quantum computing device, to explore the energy landscape of a problem and find the optimal solution.

In quantum annealing, the system is initialized in a quantum superposition of states, representing possible solutions. By gradually reducing quantum fluctuations, the system is driven towards the state with the lowest energy, corresponding to the optimal solution. This process is known as annealing.

Quantum annealing is particularly well-suited for solving combinatorial optimization problems, where finding the best configuration among a large number of possibilities is required. It offers the potential for improved efficiency compared to classical optimization algorithms for certain types of problems.

In this project, we are going to work with the Ising model in 2D. The model consists on a lattice of particles which can have the spin pointing upwards or downwards. The system is described by its energy, which depends on the interaction between the particles and an applied magnetic field, due to the magnetic moment of a particle with spin. In this way, the hamiltonian of the system is defined as:

$$ H = -J\sum_{\langle i j \rangle}\sigma_i \sigma_j -h\sum_{i} \sigma_i$$

Where, ignoring some physical constants due to their irrelevance for the combinatorial problem, the terms of the equation represent:

- $J$: the interaction strength or coupling constant between neighboring spins.
- $h$: external magnetic field applied to the system.
- $\sigma_i$: value of the spin of the i-th particle of the lattice.

Here, $\langle i j \rangle$ indicates that the interaction happens only between spins which are next to each other.

Therefore, we want to find the configuration of particles of the lattice which minimizes the energy of the system. There are several algorithms to solve this kind of problem. However, we are going to simulate the process of quantum annealing that a quantum computer would implement to solve this combinatorial problem.

In this model, we evolve a state of a mixture of all possible configurations to the solution of the problem. For that, we introduce the evolution hamiltonian $H(\lambda)$

$$ H(\lambda) = \lambda H_1 + (1-\lambda) H_0 $$

Where $H_1$ is the hamiltonian of the Ising model and:

$$ H_0 = -\sum_{i}S_{i}^{x} $$

$S_{i}^{x}$ being the operator of spin in the x direction. As this operator acts on states of the wavefunction of the spin making the amplitude of proabilities of the states $\ket{\uparrow}$ and $\ket{\downarrow}$ equal, we evolve the wavefuntion from a system with a superposition of all states to a system with a high probability of being in the state of the minimum energy of $H_1$.

# Simulation of the process


