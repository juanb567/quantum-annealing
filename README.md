# Quantum annealing simulator for the 2D Ising model

This project consists on a program to simulate the resolution of the 2D Ising model with quantum annealing.

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

Where $H_1$ is the hamiltonian of the Ising model:

$$H_1 = -J\sum_{\langle i j \rangle}S_{i}^{z} S_{j}^{z} -h\sum_{i} S_{i}^{z}$$

and:

$$ H_0 = -\sum_{i}S_{i}^{x} $$

$S_{i}^{x}$ being the operator of spin in the x direction. As this operator acts on states of the wavefunction of the spin making the amplitude of proabilities of the states $\ket{\uparrow}$ and $\ket{\downarrow}$ equal, we evolve the wavefuntion from a system with a superposition of all states to a system with a high probability of being in the state of the minimum energy of $H_1$.

# Simulation of the process

The code simulates the process of quantum annealing that a quantum computer would implement in an abstract way. For that, the code computes the matrix of the hamiltonian $H(\lambda)$ in the base of the states of $S_{i}^{z}$ and, for each step of $\lambda$, diagonalizes it, obtaining the eigenvalue with the minimum value with its corresponding eigenvector: the ground state and the solution of the system.

For that, the base of the vectors is created: all the possible configurations of the lattice, where $\ket{\uparrow}=\ket{+1}$ and $\ket{\downarrow}=\ket{-1}$. The base has a dimension of $2^N$, where $N$ is the number of particles in the lattice.

<<<<<<< HEAD
The ground state of the system is calculated as a linear combination of those states. The objective is to start from the state where all configurations are equally probable to the state where the minimum energy of $H_1$ has a probabilty of 100%. For that, the code plots the coeficients squared of each configuration of the linear combination of the ground state of $H(\lambda)$
=======
The ground state of the system is calculated as a linear combination of those states. The objective is to start from the state where all configurations are equally probable to the state where the minimum energy of $H_1$ has a probabilty of 100%. For that, the code plots the coeficients squared of each configuration of the linear combination of the ground state of $H(\lambda)$, simulating the evolution of the states of the system.

In addition, it calculates the magnetization of the system, defined as the sum of all spin-values divided by the number of particles:

$$m = \frac{\sum_i\sigma_i}{N}$$

# Contents of the project

There are 4 main contents in the project:

- In the file [quantum_annealing_simulator](https://github.com/juanb567/quantum-annealing/blob/master/quantum_annealing_simulator.py) the code starts with the definition of the parameters of the problem, where it can be chosen the values:
  - J : coupling strength
  - h : magnetic field 
  - step : increment of $\lambda$
  - n : number of particles in horizontal direction
  - m : number of particles in vertical direction

Next to it, the code contains the part where the simulation is done. After that, there are some functions to plot the evolution of the system.
- The file [functions](https://github.com/juanb567/quantum-annealing/blob/master/functions.py) contains functions to create the elements to do the simulation, such as generating the lattice of the system, the possible configurations and creating the matrices $H_0$ and $H_1$.
- The file [plots](https://github.com/juanb567/quantum-annealing/blob/master/plots.py) contains functions to plot the solutions of the simulation: the configuration of the solution of the system and the graphs of the evolution of the magnetization and the amplitudes of the system in a linear and logarithmic scale during the process.
- The file [configuration](https://github.com/juanb567/quantum-annealing/blob/master/congfiguration.txt) contains the parameters of the system, as well as the directory where the plots will be saved. To change them, the user will have to modify the default values.
- The file [tests](https://github.com/juanb567/quantum-annealing/blob/master/tests.py) contains tests implemented for the functions in the file [functions](https://github.com/juanb567/quantum-annealing/blob/master/functions.py).
- The directory [graphs](https://github.com/juanb567/quantum-annealing/blob/master/Graphs) contains the results obtained for a configuration of J = 1, h = 0.0001, step = 0.01, m = 3 and n = 3.

To execute the program, the user needs to clone the repository and type 'python quantum_annealing_simulator.py'. A directory with the images of the simulation will be saved.
>>>>>>> refs/remotes/origin/master
