import matplotlib.pyplot as mpl
from core.Current_functions import *
from core.Parameters import *

Gk =4.858637921371028e-07
Gh= 1.5564219967693693e-09
Imax =1.5744292397264454e-11
Gl=1.05601553e-10
Gi=1.5*1e-9

submyelinK = np.full(len(times), Krest)
myelinK = np.full(len(times), Ki)
Ik_s=np.full(len(times), 0.)
Ik_e=np.full(len(times), 0.)
Ih_s=np.full(len(times), 0.)
Ih_e=np.full(len(times), 0.)
Il_s=np.full(len(times), 0.)
Il_e=np.full(len(times), 0.)
Ig=np.full(len(times), 0.)
Io=np.full(len(times),0.)

V=np.full(len(times),El)

EK= E(R,T,F,Krest,Ki)

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

    Io[i]=1/((1/Gp)+(1/Gi+1/Gp)/49)*(V[i]-El)

    IG=gap_leak(V[i],EK,1,Gg,Ki,myelinK[i],Ni,Ni)
    Ig[i]=sum(IG)



    myelinK[i + 1] -= flux(Ia_e[0] + Ia_s[0] + Ik_s[i] + Ik_e[i] + IH_s[0] + IH_e[0]+IL_e[0]+IL_s[0]+IG[0], V_myelin, F, dt)
    submyelinK[i + 1] += flux(Ia_s[0] + Ik_s[i] + IH_s[0]+IL_s[0], V_submyelin, F, dt)

    V[i + 1] = V_t_dt_1(V[i], Cm, dt, Ik_s[i]+Ik_e[i]+ Ih_s[i]+Ih_e[i]+ sum(Ia_s)+ sum(Ia_e)+ Il_e[i]+Il_s[i]+Ig[i]+Io[i])



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



mpl.show()