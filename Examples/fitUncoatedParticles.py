import random
import numpy as np
import ACMHystAreaBrown
import sys

# Generate random numbers (similar to generateRandomDoubleb)
def generate_random_double(first, last):
    return random.uniform(first, last)

# Read data from the file (similar to readAreas)
def read_areas(filename):
    mu0 = 4 * np.pi * 1e-7 * 1e9  # Magnetic permeability constant
    fb = []
    areas = []

    with open(filename, 'r') as infile:
        for line in infile:
            # Ignore empty lines or comments
            if not line.strip() or line.startswith('#'):
                continue

            # Parse the line
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Error reading line: {line.strip()}")
            freq, field, area = map(float, parts)

            # Create field parameters and adjust area
            fp = ACMHystAreaBrown.FieldParameters(amplitude=field * 1e-6 * mu0, frequency=freq * 1e-6)
            fb.append(fp)
            areas.append(area / field)

    return fb, areas

# Generate initial guesses (similar to generateInitialGuesses)
def generate_initial_guesses(n_temperatures):
    initial_guesses = []
    for _ in range(n_temperatures):
        rc = generate_random_double(5, 50)
        msat = generate_random_double(0.00005, 0.0005)
        num_particles = generate_random_double(0.01, 10)

        guess = {
            "coreRadius": rc,
            "msat": msat,
            "numParticles": num_particles,
        }
        initial_guesses.append(guess)
    return initial_guesses

# Configure parameters
def main():
    # ParallelTempering configuration
    mc_temperatures = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]
    jump_size       = [0.1] * len(mc_temperatures)
    initial_guesses = generate_initial_guesses(len(mc_temperatures))

    pt_par = ACMHystAreaBrown.PTParameters()
    pt_par.temperatures   = mc_temperatures
    pt_par.numStepsFinish = 2000
    pt_par.numStepsSwap   = 100
    pt_par.maxIterations  = 2000
    pt_par.jumpSize       = jump_size
    pt_par.tolerance      = 0.00001
    pt_par.printSteps     = 10000

    # GaussNewton configuration
    gn_par = ACMHystAreaBrown.GNParameters()
    gn_par.maxIterations  = 500
    gn_par.tolerance      = 1e-4
    gn_par.regularization = 1e-6
    gn_par.printSteps     = 500

    # Read areas from a file
    
    filename  = "Areas_rc_15.data"
    fb, areas = read_areas(filename)
    
    # Additional parameters configuration
    extra_parameters = {
        "kBT": 0.0041124,
        "viscosity": 0.0009,
        "coatingWidth": 0.0,
    }

    target_parameters = {"coreRadius": 15,
                         "msat":0.000150,
                         "numParticles":1};

    # Fit the areas using the Python wrapper
    fit_result = ACMHystAreaBrown.fitAreas(
        fb, areas, pt_par, gn_par, initial_guesses, extra_parameters
    )

    # Print adjusted parameters
    print("Adjusted Parameters:")
    for key, value in fit_result[0].items():  # fit_result[0] contains the adjusted parameters
        print(f"{key}: {value}")

    print("\nTarget Parameters:")
    for key, value in target_parameters.items():
        print(f"{key}: {value}")

    # Optionally, print the adjusted errors
    print("\nStandard Errors:")
    for key, value in fit_result[1].items():  # fit_result[1] contains the adjusted errors
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
