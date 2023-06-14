import functions as func
import plots
import configparser

config = configparser.ConfigParser()

config.read('configuration.txt')

J = int(config['J']['value'])
h = float(config['h']['value'])
step = float(config['step']['value'])
n = int(config['n']['value'])
m = int(config['m']['value'])
directory = config['directory']['value']

assert n*m <= 16   # The problem is too big.

############################# SIMULATION ##################################

configurations, sites, b, param_lambda, magnetization_values, amplitudes_values = func.simulation(J, h, step, n, m)

################################ PLOTS #####################################

plots.plot_configuration(n, m, configurations, sites, b, directory)
plots.plot_magnetization(param_lambda, magnetization_values, directory)
plots.plot_amplitudes(n, m, param_lambda, amplitudes_values, directory)
plots.plot_amplitudes_log(n, m, param_lambda, amplitudes_values, directory)
