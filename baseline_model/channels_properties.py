"""
This python file plots I-V curves of the channels used in the model, based on the established parameters
"""


import numpy as np
from core.Current_functions import *
from core.Parameters import *


V = np.arange(-0.120, 0.040, 0.01)
Ih_curve = np.full(len(V), 0.)
I_kcurve = np.full(len(V), 0.)
Il_curve = np.full(len(V), 0.)
# gamma_h=Fgamma_h(El,Krest,Ki,Nrest,Ni)
gamma_l = Fgamma_l(Krest, Ki, Nrest, Ni)
for i in range(len(V)):
    Ih_curve[i] = sum(I_H(P_h(a, b, c, V[i], V_half), Gh, V[i], gamma_h, Krest, Ki, Nrest, Ni))
    I_kcurve[i] = I_KIR(R, T, F, Gk, V[i], Krest, Ki)
    Il_curve[i] = sum(gap_j_leak(V[i], gamma_l, Gl, Krest, Ki, Nrest, Ni))

fig, ax = mpl.subplots()
mpl.plot(V * 1000, I_kcurve * 1e12, color='red')
ax.set_xlabel("V(mV)")
ax.set_ylabel("I(pA)")
mpl.title("Ik-V curve at resting K+ concentration")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
mpl.axhline(y=0, color='k', linestyle='--')
mpl.axvline(x=0, color='k', linestyle='--')
mpl.axvline(x=El * 1000, color='grey', label='Resting potential')
# mpl.xlim(-120, 40)
# mpl.ylim(-500, 20)
mpl.legend()

fig, ax = mpl.subplots()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
mpl.plot(V * 1000, Ih_curve * 1e12, color='red')
ax.set_xlabel("V (mV)")
ax.set_ylabel("I(pA)")
mpl.title(" Ih-V curve at resting K+ concentration")
mpl.axhline(y=0, color='k', linestyle='--')
mpl.axvline(x=0, color='k', linestyle='--')
mpl.axvline(x=El * 1000, color='grey', label='Resting potential')
mpl.xlim(-120, 40)
mpl.ylim(-160, 20)
mpl.legend()

fig, ax = mpl.subplots()
mpl.plot(V * 1000, Il_curve, color='red')
ax.set_xlabel("V(mV)")
ax.set_ylabel("I(pA)")
mpl.title("IL-V curve at resting K+ concentration")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
mpl.axhline(y=0, color='k', linestyle='--')
mpl.axvline(x=0, color='k', linestyle='--')
mpl.axvline(x=El * 1000, color='grey', label='Resting potential')
# mpl.xlim(-120, 40)
# mpl.ylim(-500, 20)
mpl.legend()

mpl.show()