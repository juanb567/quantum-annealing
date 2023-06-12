import os
import matplotlib.pyplot as plt


def plot_configuration(n, m, configurations, sites, b, directory):
    
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
    filepath = os.path.join(directory, 'configuration.png')
    os.makedirs(directory, exist_ok=True)
    plt.savefig(filepath)


def plot_magnetization(param_lambda, magnetization_values, directory):
    
    """This module plots the magnetization of the 
    system during the process."""
    
    plt.figure()
    plt.plot(param_lambda, magnetization_values)
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Magnetization$")
    plt.xlim([0, 1])
    plt.grid()
    filepath = os.path.join(directory, 'magnetization.png')
    os.makedirs(directory, exist_ok=True)
    plt.savefig(filepath)

def plot_amplitudes(n, m, param_lambda, amplitudes_values, directory):
    
    """This module plots the amplitude of the each configuration
    during the process."""
    
    plt.figure()
    for i in range(2**(n*m)):
        plt.plot(param_lambda, amplitudes_values[:,i])
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Amplitudes$")
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.grid()
    filepath = os.path.join(directory, 'amplitudes.png')
    os.makedirs(directory, exist_ok=True)
    plt.savefig(filepath)
    
def plot_amplitudes_log(n, m, param_lambda, amplitudes_values, directory):
    
    """This module plots the amplitude of the each configuration
    during the process in a logarithmic scale."""
    
    plt.figure()
    plt.yscale(value='log')
    for i in range(2**(n*m)):
        plt.plot(param_lambda, amplitudes_values[:,i])
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Amplitudes$")
    plt.xlim([0, 1])
    plt.ylim([1e-5, 1])
    plt.grid()
    filepath = os.path.join(directory, 'amplitudes_log.png')
    os.makedirs(directory, exist_ok=True)
    plt.savefig(filepath)

########################################################################
