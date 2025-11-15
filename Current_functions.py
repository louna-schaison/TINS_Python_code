"""
This python file contains the functions used in the code to calculate current through the HCN, KIR, leak channels
It also display the I-V curves of KIR and HCN channel when ran as main.
"""

import Parameters
from Parameters import *
import numpy as np
import matplotlib.pyplot as mpl

def ATPase (Vmax,Km, c_ext):
    return (Vmax*c_ext)/(Km+c_ext)


def I_KIR(R,T,F,Gk, V, c_ext,c_int):
    """

    :param R:
    :param T:
    :param F:
    :param Gk: conductance of the KIR channels of the internode (S)
    :param V: membrane potential (V)
    :param c_ext: K+ concentration of the external compartment (M)
    :param c_int:  K+ concentration of the internal compartment (M)
    :return:
    """
    EK = E(R,T,F,c_ext,c_int) #calculating the Equilibrium potential for K+ with Nernst Law
    DV=(V - EK)
    return Gk/(1 + np.exp((DV+0.015)/0.007))*np.sqrt(c_ext)* DV

def E( R,T,F,c_ext,c_int):
    """

    :param R:
    :param T:
    :param F:
    :param c_ext: Ion concentration of the external compartment (M)
    :param c_int: Ion concentration of the internal compartment (M)
    :return: Equilibrium potential of the ion depending on the concentrations
    """
    return (R * T / F) * np.log(c_ext / c_int)


def I_leak(Gl, V, El,ENa, c_ext,c_int):
    Ek = E(c_ext, c_int) #Calcul of the potential of equilibrium of EK
    ratio_Gl_Na_K = (El - Ek) / (ENa - El)
    GlK = 1 / (1 + ratio_Gl_Na_K) #Potassium conductance ratio compared to Na+
    return GlK * Gl * (V - El)

def P_openH(a, b, c, V, V_half):
    """

    :param a:(mV-1)
    :param b:(mV-1)
    :param c:(mV-1)
    :param V: membrane potential (V)
    :param V_half: potential of half-activation (V)
    :return: Probability of HCN channels to be open at a given V
    """
    alpha = a * np.exp(-b * (V * 1000 - V_half))
    beta = a * np.exp(c * (V * 1000 - V_half))
    return alpha/(alpha + beta)

def tau_h(a,b,c,V,V_half):
    alpha = a * np.exp(-b * (V * 1000 - V_half))
    beta = a * np.exp(c * (V * 1000 - V_half))
    return 1/(alpha + beta)

def I_H(Ph, Gh,V,gamma,c_ext,c_int,Ni=0.001,Nrest= 0.140 ):
    """
    :param Ph: probability of HCN channels to be open at a given V
    :param Gh: conductance of the HCN channel in the internode (S)
    :param gamma: permeability ratio between K+ and Na+ (Pk+/PNa+)
    :param V: membrane potential (V)
    :param c_ext: K+ concentration of the external compartment (M)
    :param c_int: K+ concentration of the internal compartment (M)
    :param Ni: Na+ concentration of the internaal compartment (M)
    :param Nrest:Na+ concentration of the external compartment (M)
    :return: Current (in A) based on H-H model
    """
    Eh= (R*T/F)*np.log((gamma*c_ext+Nrest)/(gamma*c_int+Ni))
    return Gh * Ph  * (V - Eh)

def I_H_k (Ih,Ena,Ek,V,gamma):
    return (V-Ek)*Ih / ((V - Ek)+(1/gamma)*(V-Ena))


if __name__ == "__main__":
    """ Tracing the I-V curves of both channels
    """
    V = np.arange(-0.120, 0.040, 0.01)
    Ih_curve = np.full(len(V), 0.)
    I_kcurve = np.full(len(V), 0.)
    Ih_K_curve = np.full(len(V), 0.)
    for i in range(len(V)):
        Ph = P_openH(a, b, c, V[i], V_half)
        Ih_curve[i] = I_H( Ph, Gh, V[i],gamma,Krest,Ki )
        I_kcurve[i] = I_KIR(R, T, F, Gk, V[i], Krest, Ki)
        Ih_K_curve[i] = I_H_k(Ih_curve[i],E(R,T,F,Nrest,Ni), E(R,T,F,Krest,Ki),V[i],gamma)
    fig,ax=mpl.subplots()
    mpl.plot(V*1000, I_kcurve*1e12, color='red')
    ax.set_xlabel("V(mV)")
    ax.set_ylabel("I(pA)")
    mpl.title("Ik-V curve at resting K+ concentration")
    mpl.axhline(y=0, color='k',linestyle='--')
    mpl.axvline(x=0, color='k',linestyle='--')
    mpl.axvline(x=-72, color='grey',label='Resting potential')
    mpl.xlim(-120, 40)
    mpl.ylim(-1e-1, 2e-2)
    mpl.legend()
    fig, ax = mpl.subplots()
    mpl.plot(V * 1000, Ih_curve*1e12, color='red')
    ax.set_xlabel("V (mV)")
    ax.set_ylabel("I(pA)")
    mpl.title(" Ih-V curve at resting K+ concentration")
    mpl.axhline(y=0, color='k', linestyle='--')
    mpl.axvline(x=0, color='k', linestyle='--')
    mpl.axvline(x=-72, color='grey', label='Resting potential')
    mpl.xlim(-120, 40)
    mpl.ylim(-4 , 1)
    mpl.legend()
    mpl.show()