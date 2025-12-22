"""
This python file contains all the parameters used in the code.
"""

import numpy as np

R = 8.3 # J mol-1 K-1
F = 95600 # C mol-1
T = 273 + 32 # Recording temperature

# K+ concentrations :
Ki = 0.13 # M, inside
Krest = 0.0035 # Resting Ko. M., outside

# Na+ concentrations
Ni=0.002 # M, inside
Nrest= 0.150 # Resting M., outside


# Myelin electrical properties
Gm = 60*1e-9 # Membrane total conductance S
El = -0.086# V resting potential
Gg = 540*1e-12 #conductance of the connection with oligo and astrocytes
Gp= 6.3*1e-9# Conductance of the processes
Ga= 30*1e-9 # Conductance of astrocytes membrane S

# Kir parameters
Gk = 9.286*1e-9 #S

#HCN parameters
Gh=0.484*1e-9#S
Eh=-0.03 #V
a = 0.63 #mV-1
b = 0.063  #mV-1
c = 0.079  #mV-1
V_half = -100.0 #mV
gamma_h= 2.6

#ATPase NA+/K+
KmK = 0.005 # M the affinity of the ATPase Na+/K+ for K+
KmNa= 0.01 #M/L
Imax =3.15025237e-11 #A

#Hemichannels leak
Gl=5.65744323e-07
gamma_l= 1.3


#  Stimulation (axonal activity characteristics)
dt = 1e-6 # s, 10 milliseconds of delta t
APefflux = 0.002 # M change in K(submyelin) with AP.
APtimes = np.arange(0,0.2,0.01) # Roughly 20 APs at 10ms intervals. # this mimicks the stimulation in the paper
times = np.arange(0,4,1e-6) #time-window studied


# Axon-myelin properties
Cm = 24*1e-12 # farad #capacitance of the myelin
tau=25.24*1e-3
n_internode =50
Vint = 6e-15 # litres. Volume of the internode peri-axonal space.
l_internode= 5e-5 #m 50um
r_axon= 0.5*1e-6 #m 0.5um
e_myelin = 250*1e-9 #m 250 nm
V_myelin = (((e_myelin+r_axon)**2*np.pi-r_axon**2*np.pi)*l_internode)*1000#en L
e_submyelin = 2*1e-8 #m
V_submyelin = (((e_submyelin+r_axon)**2*np.pi-r_axon**2*np.pi)*l_internode)*1000#en L (bonne estimation par Vint)
l_paranode = 1.9 *1e-6
S_section_submyelin= (e_submyelin+r_axon)**2*np.pi-r_axon**2*np.pi

D= 1.85*1e-9 #diffusion coefficent for K+