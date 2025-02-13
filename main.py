import matplotlib.pyplot as plt
import numpy as np

n_array = np.arange(1, 2, 0.001) #Creates list of refractive indices

def michelson_calc (T, theta_i, theta_f, wavelength, N):
    wavelength *= 10**(-6) #Nanometers to millimeters
    theta_i = np.deg2rad(theta_i)
    theta_f = np.deg2rad(theta_f) #Degrees to radians

    value_array = []

    # Calculates both functions for all values of n
    for n in n_array:
        sqrt_i = np.sqrt(1-((np.sin(theta_i)**2)/n**2))
        left_i =  -(2 * T / sqrt_i) * np.cos(theta_i - np.arccos(sqrt_i))
        right_i = (2 * n * T) / sqrt_i
        complete_i = left_i + right_i

        sqrt_f = np.sqrt(1 - ((np.sin(theta_f) ** 2) / n ** 2))
        left_f = -(2 * T / sqrt_f) * np.cos(theta_f - np.arccos(sqrt_f))
        right_f = (2 * n * T) / sqrt_f
        complete_f = left_f + right_f

        # Subtracts the measured number of fringes from the theoretical number
        value = ((complete_f - complete_i) / wavelength) - N
        value_array.append(value)

    # Finds the refractive index used in the previous equation where the answer
    # was closest to zero.
    index = np.interp(0, value_array, n_array)
    print('Index: ' + str(index))

    # Plots a graph that visualise the calculations
    fig, ax = plt.subplots()
    ax.plot(n_array,value_array)
    ax.axis
    ax.axhline(c='grey', lw=1)

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{textgreek}')

    plt.ylabel(r'$\dfrac{f (θ_f) - f (θ_i) }{λ}$')
    plt.xlabel('Refractive index')
    plt.subplots_adjust(left=0.15, bottom=0.1)
    plt.savefig('michelsonGraph.png', )
    plt.show()