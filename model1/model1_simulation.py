import matplotlib.pyplot as mpl
from core.Current_functions import *
from core.Parameters import *

Gk = 8.802316894091538e-07
Gh= 1.5564219967693693e-09
Imax =1.9692280100250965e-11
Gl=1.22407553e-10



submyelinK = np.full(len(times), Krest)
myelinK = np.full(len(times), Ki)
Ik_s=np.full(len(times), 0.)
Ik_e=np.full(len(times), 0.)
Ih_s=np.full(len(times), 0.)
Ih_e=np.full(len(times), 0.)
Il_s=np.full(len(times), 0.)
Il_e=np.full(len(times), 0.)
V=np.full(len(times),El)


for i in range(len(times) - 1):

    submyelinK[i + 1] = submyelinK[i]
    myelinK[i + 1] = myelinK[i]
    if sum(APtimes < (times[i] - dt)) < sum(APtimes < times[i]):  # meaning that the t+dt is reaching the timing of a new AP
        submyelinK[i + 1] += APefflux

    Ia_s = ATPase(Imax / 2, KmK, KmNa, submyelinK[i], Ni)
    Ia_e = ATPase(Imax / 2, KmK, KmNa, Krest, Ni)

    Ik_s[i] = I_KIR(R, T, F, Gk / 2, V[i], submyelinK[i], myelinK[i])
    Ik_e[i] = I_KIR(R, T, F, Gk / 2, V[i], Krest, myelinK[i])

    IH_s = I_H(P_h(a, b, c, V[i], V_half), Gh / 2, V[i], gamma_h, submyelinK[i], myelinK[i], Nrest, Ni)
    IH_e= I_H(P_h(a, b, c, V[i], V_half), Gh / 2, V[i], gamma_h, Krest, myelinK[i], Nrest, Ni)
    Ih_e[i]=sum(IH_e)
    Ih_s[i]=sum(IH_s)

    IL_s=gap_leak(V[i],0,gamma_l,Gl/2,submyelinK[i],myelinK[i],Nrest,Ni)
    IL_e = gap_leak(V[i], 0, gamma_l, Gl/2, Krest, myelinK[i], Nrest, Ni)
    Il_e[i]=sum(IL_e)
    Il_s[i]=sum(IL_s)


    myelinK[i + 1] -= flux(Ia_e[0] + Ia_s[0] + Ik_s[i] + Ik_e[i] + IH_s[0] + IH_e[0]+IL_e[0]+IL_s[0], V_myelin, F, dt)
    submyelinK[i + 1] += flux(Ia_s[0] + Ik_s[i] + IH_s[0]+IL_s[0], V_submyelin, F, dt)

    V[i + 1] = V_t_dt_1(V[i], Cm, dt, Ik_s[i]+ Ik_e[i]+ Ih_s[i]+ Ih_e[i]+sum(Ia_s)+ sum(Ia_e)+ Il_e[i]+Il_s[i])



fig, ax = mpl.subplots()
mpl.plot(times, V*1e3, color='green')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Vm (mV)")
mpl.ylim(-100,0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig, ax1 = mpl.subplots()

line1 = ax1.plot(times, submyelinK, color='blue', label="Submyelin")
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("[K+] (M)",color='blue')
ax2 = ax1.twinx()
line2 = ax2.plot(times, myelinK, color='darkblue', label="Intracellular")
ax2.set_ylabel("[K+] (M)",color='darkblue')
ax2.set_ylim(0.130,0.140)

lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='lower right')

mpl.title("Evolution of the intracellular [K+] and submyelin [K+] overtime")


fig,ax=mpl.subplots()
mpl.plot(times, Ik_e*1e12,color='red',label="external interface")
mpl.plot(times, Ik_s*1e12,color='darkred', label='submyelin interface')
ax.set_xlabel("Time (s)")
ax.set_ylabel("I (pA)")
mpl.title("KIR current  over time ")
mpl.legend()


fig,ax=mpl.subplots()
mpl.plot(times, Ih_e*1e12,color='red',label="external interface")
mpl.plot(times, Ih_s*1e12,color='darkred', label='submyelin interface')
ax.set_xlabel("Time (s)")
ax.set_ylabel("I (pA)")
mpl.title("HCN current  over time ")
mpl.legend()

fig,ax=mpl.subplots()
mpl.plot(times, Il_e*1e12,color='red',label="external interface")
mpl.plot(times, Il_s*1e12,color='darkred', label='submyelin interface')
ax.set_xlabel("Time (s)")
ax.set_ylabel("I (pA)")
mpl.title("Leak current  over time ")
mpl.legend()


mpl.show()