from scipy.integrate import odeint
from numpy import pi, zeros, linspace, sin
from calculateArea import getArea
import numpy as np
from scipy.signal import find_peaks

#----------------------------------------utils----------------------------------------------------#
def getLastCycle(field, magnetization, pointsPerCycle=-1, ncycles=-1):
     if (pointsPerCycle>0 and ncycles>0):
          if (pointsPerCycle*ncycles + 1 == len(field)):
               return field[-pointsPerCycle-1:], magnetization[-pointsPerCycle-1:]
          else:
               return field[pointsPerCycle*ncycles-1:], magnetization[pointsPerCycle*ncycles-1:]

     else:
          lastField = field[-1]
          distanceToLastField = abs(field[:-2]-lastField)
          minimumDistances, _ = find_peaks(-distanceToLastField)
          return field[minimumDistances[-2]:], magnetization[minimumDistances[-2]:]
     
#-------------------------------Brownian relaxation-----------------------------------------------#

def compute_dadt(a,t,dr,alpha0,w):
    nindex = len(a)
    field = alpha0*sin(t*w)
    dan = zeros(nindex)
    dan[0] = 0
    for n in range(1,nindex-1):
        dan[n]=dr*n*(n+1)*(-a[n]+field*(a[n-1]/(2*n-1)-a[n+1]/(2*n+3)))
    return dan

def integrateFP(Dr,alpha0, w, nindex, time, an_0 = None):
     wadim = 0.5*w/Dr

     if an_0 is None:
          an_0 = zeros(nindex)
          an_0[0] = 0.5
          an_0[1] = -(alpha0*wadim/(2*(wadim**2-1)))
     an_t = odeint(compute_dadt, an_0, time, args = (Dr,alpha0,w))
     return an_t

def computeAn_t(T, vis, rh, m0, f, b0,
                nindex = 50,
                pointsPerCycle = 1000,
                tolerance = 1e-4,
                an_0 = None):
     
     Dr           = T/(8*pi*vis*rh**3)
     alpha0       = m0*b0/T
     dt           = 1./(f*pointsPerCycle)
     
     
     numberPoints  = pointsPerCycle + 1
     time          = linspace(0,1.0/f, numberPoints)
     an_t          = integrateFP(Dr, alpha0, 2*pi*f, nindex, time, an_0)
     magnet        = 2*an_t[:,1]/3.0
     error         = abs(max(magnet)+min(magnet))/max(magnet)
     while error > tolerance:
          an_t          = integrateFP(Dr, alpha0, 2*pi*f, nindex, time, an_t[-1,:])
          magnet        = 2*an_t[:,1]/3.0
          error         = abs(max(magnet)+min(magnet))/max(magnet)
     return time, an_t

def computeCycle(T, vis, rh, m0, f, b0,
                  nindex_0=50, d_nindex=50,
                  pointsPerCycle=1000,
                  toleranceCycle=1e-4,
                  toleranceArea=1e-5):
    
    def compute_single_cycle(T, vis, rh, m0, f, b0, nindex, pointsPerCycle, tolerance, an_0=None):
        time, an_t = computeAn_t(T, vis, rh, m0, f, b0,
                                 nindex=nindex,
                                 pointsPerCycle=pointsPerCycle,
                                 tolerance=tolerance,
                                 an_0=an_0)
        field = b0 * np.sin(2 * np.pi * f * time)
        magnet = 2 * an_t[:, 1] / 3.0
        return field, magnet, an_t

    def update_area(cycle):
        return getArea(cycle)[0]

    # Primera iteración
    nindex = nindex_0
    field, magnet, an_t = compute_single_cycle(T, vis, rh, m0, f, b0, nindex, pointsPerCycle, toleranceCycle)
    
    cycle = np.column_stack((field, magnet))
    area_old = update_area(cycle)

    # Siguiente iteración inicial
    nindex += d_nindex
    an_0 = np.zeros(nindex)
    an_0[:nindex - d_nindex] = an_t[-1, :]

    field, magnet, an_t = compute_single_cycle(T, vis, rh, m0, f, b0, nindex, pointsPerCycle, toleranceCycle, an_0)
    cycle = np.column_stack((field, magnet))
    area_new = update_area(cycle)

    # Ciclo de ajuste por tolerancia
    error = abs(area_new - area_old) / area_new
    while error > toleranceArea:
        area_old = area_new
        nindex += d_nindex
        an_0 = np.zeros(nindex)
        an_0[:nindex - d_nindex] = an_t[-1, :]

        field, magnet, an_t = compute_single_cycle(T, vis, rh, m0, f, b0, nindex, pointsPerCycle, toleranceCycle, an_0)
        cycle = np.column_stack((field, magnet))
        area_new = update_area(cycle)
        error = abs(area_new - area_old) / area_new

    return field, magnet

def computeMagnetization(T, vis, rh, m0, f, b0,
                         nindex = 30, ncycles = None,
                         pointsPerCycle = 500,
                         tolerance = None):

     #Check that only tolerance or ncycles is not none
     if tolerance is not None and ncycles is not None:
          raise ValueError("Specify only one: tolerance or ncycles.")
     
     if tolerance is None and ncycles is None:
          raise ValueError("Specify one: tolerance or ncycles.")
     
     
     Dr           = T/(8*pi*vis*rh**3)
     alpha0       = m0*b0/T
     dt           = 1./(f*pointsPerCycle)

     if ncycles is None:
          numberPoints  = pointsPerCycle + 1
          time          = linspace(0,1.0/f, numberPoints)
          an_t          = integrateFP(Dr, alpha0, 2*pi*f, nindex, time)
          magnet        = 2*an_t[:,1]/3.0
          error         = abs(max(magnet)+min(magnet))/max(magnet)
          while error > tolerance:
               print(error)
               an_t          = integrateFP(Dr, alpha0, 2*pi*f, nindex, time, an_t[-1,:])
               magnet        = 2*an_t[:,1]/3.0
               error         = abs(max(magnet)+min(magnet))/max(magnet)
     else:
          numberPoints = int(ncycles*pointsPerCycle)+1
          an_t         = integrateFP(Dr, alpha0, 2*pi*f, nindex, time)
          magnet       = 2*an_t[:,1]/3.0
          
     field         = b0*sin(2*pi*f*time)
     return time, field, magnet

#--------------------------------------Neel relaxation-----------------------------------------#

def compute_dadt_neel(a,t,tneel,alpha0,alphak,w):
    nindex = len(a)
    field = alpha0*sin(t*w)
    dan = zeros(nindex)
    dan[0] = 0
    for n in range(1,nindex-2):
        termField = (-a[n]+field*(a[n-1]/(2*n-1)-a[n+1]/(2*n+3)))
        termAnisotropy = alphak*(n*a[n]/((2*n-1)*(2*n+1))-(n+1)*a[n]/((2*n+1)*(2*n+3))-(n+2)*a[n+2]/((2*n+3)*(2*n+5)))
        if n>1:
            termAnisotropy += alphak*(n-1)*a[n-2]/((2*n-3)*(2*n-1))
        dan[n]=0.5/tneel*n*(n+1)*(termField+termAnisotropy)
    return dan


def integrateFP_neel(tn, alpha0, alphak, w, nindex, time):
    an_0 = zeros(nindex)
    an_0[0] = 0.5
    an_t = odeint(compute_dadt_neel, an_0, time, args = (tn,alpha0, alphak, w))
    return an_t


def computeMagnetization_neel(kbT, rc, K, damp, m0, gamma0, f, b0,
                              nindex = 30, ncycles = 10, pointsPerCycle = 500, lastCycle = None):

    t_neel = m0*(1+damp*damp)/(2*kbT*gamma0*damp)
    alpha0 = m0*b0/kbT
    Vc = 4./3.*np.pi*rc*rc*rc
    alphak = 2*K*Vc/kbT
    time = linspace(0,ncycles/f, pointsPerCycle*ncycles)
    if (lastCycle is not None):
         lastCycle+=time[-1]
         time = np.concatenate((time,lastCycle))
    an_t = integrateFP_neel(t_neel, alpha0, alphak, 2*pi*f, nindex, time)
    magnetization = 2*an_t[:,1]/3.0
    field = b0*sin(2*pi*f*time)
    return time, field, magnetization


def computeCycle_neel(kbT, rc, K, damp, m0, gamma0, f, b0, nindex = 30,
                      ncycles = 10, pointsPerCycle = 500, lastCycle = None):
     time, field, magnetization = computeMagnetization_neel(kbT, rc, K, damp, m0, gamma0, f, b0,
                                                            nindex, ncycles, pointsPerCycle,
                                                            lastCycle)
     
     field, magnetization = getLastCycle(field, magnetization, pointsPerCycle, ncycles)
     return field, magnetization       



