import numpy as np
from scipy.integrate import odeint
from .integrateAreaCycles import calculateArea_MT

# Constants
PI = np.pi

def computeLegendreCoefficientsDerivative(legendreCoeffs, time, h0, omega, brownTime):
    """
    Compute the time derivatives of Legendre coefficients.
    """
    numCoeffs    = len(legendreCoeffs)
    h            = h0 * np.sin(time * omega)
    dLegendre    = np.zeros(numCoeffs)
    dLegendre[0] = 0  # Zeroth coefficient does not evolve
    dr           = 0.5 / brownTime  # Rotational diffusion coefficient
    for n in range(1, numCoeffs - 1):
        dLegendre[n] = dr * n * (n + 1) * (
            -legendreCoeffs[n]
            + h * (
                legendreCoeffs[n - 1] / (2 * n - 1)
                - legendreCoeffs[n + 1] / (2 * n + 3)
            )
        )
    return dLegendre

def integrateLegendreCoefficients(h0, omega, brownTime, numCoeffs, time, initialConditions=None):
    """
    Integrate Legendre coefficients over time using the Fokker-Planck equation.
    """
    if initialConditions is None:
        initialConditions    = np.zeros(numCoeffs)
        wadim                = brownTime * omega
        initialConditions[0] = 0.5
        initialConditions[1] = -(h0 * wadim / (2 * (wadim**2 - 1)))

    legendreCoeffsT = odeint(computeLegendreCoefficientsDerivative,
                             initialConditions,
                             time,
                             args=(h0, omega, brownTime)
                             )
    return legendreCoeffsT


def computeTimeDependentLegendreCoeffs(h0, frequency, brownTime, numCoeffs=50,
                                       pointsPerCycle=1000, tolerance=1e-4, initialConditions=None):
    """
    Compute the magnetization of the system by integrating the Legendre coefficients until an equilibrium cycle is reached.
    """
    numPoints = pointsPerCycle + 1
    time      = np.linspace(0, 1.0 / frequency, numPoints)
    omega     = 2 * PI * frequency
    error     = np.inf

    # Refine until the error is within tolerance
    #TODO: Maximum number tries
    while error > tolerance:
        legendreCoeffsT   = integrateLegendreCoefficients(h0, omega, brownTime, numCoeffs, time, initialConditions)
        magnetization     = 2 * legendreCoeffsT[:, 1] / 3.0
        error             = np.abs(np.max(magnetization) + np.min(magnetization)) / np.max(magnetization)
        initialConditions = legendreCoeffsT[-1, :]
        
    return time, legendreCoeffsT

def computeMagnetizationCycle(amplitude, frequency, kBT, viscosity, hydroRadius, magneticMoment,
                              numCoeffsStart=50, coeffsStep=50, pointsPerCycle=1000,
                              toleranceCycle=1e-4, toleranceArea=1e-5):
    """
    Compute the magnetization cycle by iteratively refining the Legendre coefficients
    until the area under the curve converges within a given tolerance.
    """
    # Initial parameters
    brownTime = kBT / (8 * PI * viscosity * hydroRadius**3)
    h0        = magneticMoment * amplitude / kBT
    omega     = 2 * PI * frequency

    def computeCycleAndArea(numCoeffs, initialConditions=None):
        """
        Computation the magnetization cycle and area for a given number of coefficients.
        """
        time, legendreCoeffsT = computeTimeDependentLegendreCoeffs(
            h0, frequency, brownTime, numCoeffs, pointsPerCycle, toleranceCycle, initialConditions
        )
        field         = amplitude * np.sin(omega * time)
        magnetization = 2 * legendreCoeffsT[:, 1] / 3.0
        area          = calculateArea_MT(time, magnetization, amplitude, omega)
        return field, magnetization, legendreCoeffsT, area

    # Initial computation
    field, magnetization, legendreCoeffsT, areaOld = computeCycleAndArea(numCoeffsStart)

    # Iterate to refine the area
    error = np.inf
    numCoeffs = numCoeffsStart
    while error > toleranceArea:
        numCoeffs += coeffsStep
        initialConditions = np.zeros(numCoeffs)
        initialConditions[:numCoeffs - coeffsStep] = legendreCoeffsT[-1, :]
        field, magnetization, legendreCoeffsT, areaNew = computeCycleAndArea(numCoeffs,
                                                                             initialConditions)
        error   = np.abs(areaNew - areaOld) / areaNew
        areaOld = areaNew

    return field, magnetization
