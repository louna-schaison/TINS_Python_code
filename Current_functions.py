"""
This python file contains the functions used in the code to calculate current through the HCN, KIR, leak channels
It also display the I-V curves of KIR and HCN channel when ran as main.
"""

import Parameters
from Parameters import *
import numpy as np
import matplotlib.pyplot as mpl
from scipy.optimize import root


def ATPase (Imax, KmK, KmNa, Ke, Ni=Ni):
    F= (1+KmK/Ke)**(-2)+ (1+KmNa/Ni)**(-3)
    Ina= 3*Imax*F
    Ik= -2*Imax*F
    return(Ik, Ina)


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
    return Gk*P_k(R,T,F,V,c_ext, c_int)*np.sqrt(c_ext)* DV

def P_k(R,T,F,V,c_ext,c_int):
    EK = E(R,T,F,c_ext,c_int) #calculating the Equilibrium potential for K+ with Nernst Law
    DV=(V - EK)
    return 1/(1+np.exp((DV+0.015)/0.007))


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



def P_h(a, b, c, V, V_half):
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

def Fgamma_h(Eh,c_ext,c_int,Nrest,Ni):
    EK=E(R,T,F,c_ext,c_int)
    ENa= E(R,T,F,Nrest,Ni)
    return (ENa/Eh-1)/(1-EK/Eh)


def Fgamma_l(c_ext,c_int,Nrest,Ni):
    EK=E(R,T,F,c_ext,c_int)
    ENa= E(R,T,F,Nrest,Ni)
    return -ENa/EK

def Gi(Gm, Gp, Gg,Ga, n_internode):
    return 1/((n_internode/Gm)-1/Gp)-1/((1/Ga)+(1/Gg))

def gap_j_leak(V,gamma_l,Gl,c_ext,c_int,Nrest, Ni ):
    EK=E(R,T,F,c_ext,c_int)
    ENa= E(R,T,F,Nrest,Ni)
    IlK= Gl * (gamma_l/(1+gamma_l))* (V - EK)
    IlNa = Gl * (1 / (1 + gamma_l))  * (V - ENa)
    return IlK, IlNa


def gap_j(R,T,F,myelinK,V,Gg,El,Ki=Ki):
    return Gg*(V-El)/(np.exp(F*(V-El)/(R*T))-1)*(myelinK-Ki*np.exp(F*(V-El)/(R*T)))

def I_H(Ph, Gh,V,gamma_h,c_ext,c_int,Ni=0.001,Nrest= 0.140 ):
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
    EK=E(R,T,F,c_ext,c_int)
    ENa= E(R,T,F,Nrest,Ni)
    IhK= Gh * (gamma_h/(1+gamma_h))* Ph  * (V - EK)
    IhNa = Gh * (1 / (1 + gamma_h)) * Ph * (V - ENa)
    return IhK, IhNa #pour l'instant return le courant sortant de K (positif) entraint de Na+ (negatif)

def diffK(D,l_paranode, S,V_submyelin,submyelinK, extK):
    return -(S/V_submyelin)*(D/l_paranode)*(submyelinK-extK)

def equations(vars):
    Imax, Gk, Gh , Gl = vars
    eq1 =  1e10*(I_KIR(R,T,F,Gk,El,Krest,Ki)+sum(I_H(P_openH(a,b,c,El,V_half),Gh,El,gamma,Krest,Ki,Ni,Nrest))+sum(ATPase(Imax, KmK,KmNa,Krest,Ni)))+Gl*El
    eq2 = dt*(I_KIR(R,T,F,Gk,El,Krest,Ki)+I_H(P_openH(a,b,c,El,V_half),Gh,El,gamma,Krest,Ki,Ni,Nrest)[0]+ATPase(Imax, KmK,KmNa,Krest,Ni)[1])/(V_myelin*F)
    eq3= Gi(Gm,Gp,Gg,Ga,n_internode) - Gk*(1/(1 + np.exp((El-E(R,T,F,Krest,Ki)+0.015)/0.007)))*np.sqrt(Krest) -Gh*P_openH(a,b,c,El,V_half)-Gl
    return [eq1, eq2, eq3]

initial_guess = [1e-10,200*1e-10, 10*1e-10]
#solution = root(equations, initial_guess,method='lm'  # Levenberg-Marquardt )
#Imax, Gk, Gh = solution.x
#print(solution.x)


if __name__ == "__main__":
    """ Tracing the I-V curves of both channels
    """
    V = np.arange(-0.120, 0.040, 0.01)
    Ih_curve = np.full(len(V), 0.)
    I_kcurve = np.full(len(V), 0.)
    Il_curve = np.full(len(V), 0.)
    gamma_h=Fgamma_h(El,Krest,Ki,Nrest,Ni)
    gamma_l = Fgamma_l( Krest, Ki, Nrest, Ni)
    for i in range(len(V)):
        Ph = P_openH(a, b, c, V[i], V_half)
        Ih_curve[i] = sum(I_H( Ph, Gh, V[i],gamma_h,Krest,Ki ))
        I_kcurve[i] = I_KIR(R, T, F, Gk, V[i], Krest, Ki)
        Il_curve[i]= sum(gap_j_leak(V[i],gamma_l,Gl,Krest,Ki,Nrest,Ni))

    fig,ax=mpl.subplots()
    mpl.plot(V*1000, I_kcurve*1e12, color='red')
    ax.set_xlabel("V(mV)")
    ax.set_ylabel("I(pA)")
    mpl.title("Ik-V curve at resting K+ concentration")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    mpl.axhline(y=0, color='k',linestyle='--')
    mpl.axvline(x=0, color='k',linestyle='--')
    mpl.axvline(x=El*1000, color='grey',label='Resting potential')
    #mpl.xlim(-120, 40)
    #mpl.ylim(-500, 20)
    mpl.legend()

    fig, ax = mpl.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    mpl.plot(V * 1000, Ih_curve*1e12, color='red')
    ax.set_xlabel("V (mV)")
    ax.set_ylabel("I(pA)")
    mpl.title(" Ih-V curve at resting K+ concentration")
    mpl.axhline(y=0, color='k', linestyle='--')
    mpl.axvline(x=0, color='k', linestyle='--')
    mpl.axvline(x=El*1000, color='grey', label='Resting potential')
    mpl.xlim(-120, 40)
    mpl.ylim(-160 , 20)
    mpl.legend()

    fig, ax = mpl.subplots()
    mpl.plot(V * 1000, Il_curve , color='red')
    ax.set_xlabel("V(mV)")
    ax.set_ylabel("I(pA)")
    mpl.title("IL-V curve at resting K+ concentration")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    mpl.axhline(y=0, color='k', linestyle='--')
    mpl.axvline(x=0, color='k', linestyle='--')
    mpl.axvline(x=El * 1000, color='grey', label='Resting potential')
    #mpl.xlim(-120, 40)
    #mpl.ylim(-500, 20)
    mpl.legend()


    mpl.show()