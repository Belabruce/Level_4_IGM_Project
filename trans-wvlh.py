from __future__ import division
import numpy
import matplotlib.pyplot as pyplot
import scipy
import scipy.constants as const
from scipy import integrate



c = const.c                     #speed of light in a vacuum (m s^-1)  
k = const.k                     #20836612000.0               #boltzmann constant (Hz K^-1)
mh = const.m_p + const.m_e      #mass of neutral hydrogen (Kg)  
pi = numpy.pi
l0 = 1216e-10                   #Lyman-Alpha absorption line wavelength (m)

T = 10000                      #Temperature (K)
sigma_0 = 6.3e-18
n_dens = 1e-3
l_path = 1e+21
n_sigma = 5

def delta_l(l0,T):
    del_l = (l0/c)*((2*k*T)/mh)**0.5
    return del_l
        
def phi_l(l, l0, delta_l):
    phi = numpy.e**((-(l-l0)**2)/delta_l**2)    #* (1/(pi*delta_l**2)**0.5)
    return phi



del_l = delta_l(l0,T)

range_low = l0 - (n_sigma*del_l)
range_high = l0 + (n_sigma*del_l)
lambdas = numpy.linspace(range_low,range_high,1000)
lambdas = numpy.array(lambdas)

phis = []
for i in range(len(lambdas)):
    phis.append(phi_l(lambdas[i], l0, del_l))
phis = numpy.array(phis)

#integral = scipy.integrate.simps(phis, lambdas)
#phis = phis/integral
#integral = scipy.integrate.simps(phis, lambdas)
#print integral

sig = phis * sigma_0
tau = sig*n_dens*l_path
Tr = numpy.e**(-tau)

vel =  c*((lambdas/l0)-1)/1000
    
pyplot.figure()
pyplot.plot(vel, tau)
pyplot.show()