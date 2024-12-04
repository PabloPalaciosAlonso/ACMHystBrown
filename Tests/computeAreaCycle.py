import ACMHystAreaBrown as acmb
import numpy as np
import scientificPlot as splt
import numpy as np
from scipy.integrate import simps

b0    = 2
f     = 8
omega = 2*np.pi*f

time          = np.linspace(0, 2/f, 50000)
field         = b0*np.sin(omega*time)
magnetization = b0/(3*(omega**2+1))*(np.sin(omega*time) + omega*np.cos(omega*time))

A1 = acmb.calculateArea_MT(time,  magnetization, b0, omega)
A2 = acmb.calculateArea_MB(field, magnetization)

print(A1, A2, np.pi*b0**2*omega/(3*(1+omega**2)))

