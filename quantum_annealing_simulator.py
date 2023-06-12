import functions as func
import plots
import configparser

# Create a ConfigParser object
config = configparser.ConfigParser()

# Read the variables from the text file
config.read('configuration.txt')

# Get the values of the variables
J = int(config['J']['value'])
h = float(config['h']['value'])
step = float(config['step']['value'])
n = int(config['n']['value'])
m = int(config['m']['value'])
directory = config['directory']['value']

assert n*m <= 16   # The problem is too big.

configurations, sites, b, param_lambda, magnetization_values, amplitudes_values = func.simulation(J, h, step, n, m)


############################# GRAPHICS ##################################

plots.plot_configuration(n, m, configurations, sites, b, directory)
plots.plot_magnetization(param_lambda, magnetization_values, directory)
plots.plot_amplitudes(n, m, param_lambda, amplitudes_values, directory)
plots.plot_amplitudes_log(n, m, param_lambda, amplitudes_values, directory)
