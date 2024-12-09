import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d


# Constants
PI = np.pi

def findStartOfPeriods(time, omega):
    """
    Find the indices in time where the field crosses zero from non-positive to positive.
    """
    field = np.sin(omega * time)
    zero_crossing_indices = np.where((field[:-1] <= 0) & (field[1:] > 0))[0]
    return zero_crossing_indices

def calculateArea_MT(time, magnetization, b0, omega, initialPeriod=0):
    startOfPeriods  = findStartOfPeriods(time,omega)
    numberOfPeriods = len(startOfPeriods)-initialPeriod
    t_int           = time[startOfPeriods[initialPeriod]:]
    mz_int          = magnetization[startOfPeriods[initialPeriod]:]
    integrand       = b0 * omega * mz_int * np.cos(omega * t_int)
    area            = np.trapz(integrand, t_int)/numberOfPeriods
    return np.abs(area)


def findStartOfCycles(x, y, initial_point=0):
    upward_crossings   = np.where((x[:-1] <= initial_point) & (x[1:] > initial_point))[0]
    downward_crossings = np.where((x[:-1] >= initial_point) & (x[1:] < initial_point))[0]
    if upward_crossings[0]>downward_crossings[0]:
        startOfCycles = downward_crossings
    else:
        startOfCycles = upward_crossings

    return startOfCycles

def calculateArea_MB(field, magnetization, initialCycle=0):
    """
    Calculate the enclosed area for magnetization (magnetization vs field) starting from a given cycle.
    """

    # Find the starting indices of each cycle
    startOfCycles  = findStartOfCycles(field, magnetization)
    numberOfCycles = len(startOfCycles)-initialCycle

    # Determine the range of field and magnetization values to consider
    startIndex = startOfCycles[initialCycle]
    field_int  = field[startIndex:]
    mz_int     = magnetization[startIndex:]

    # Ensure the cycle is closed
    if field_int[0] != field_int[-1] or mz_int[0] != mz_int[-1]:
        field_int = np.append(field_int, field_int[0])
        mz_int = np.append(mz_int, mz_int[0])

    # Calculate the enclosed area using the shoelace formula
    area = 0.5 * np.abs(np.sum(field_int[:-1] * mz_int[1:] - field_int[1:] * mz_int[:-1]))

    return area/numberOfCycles
