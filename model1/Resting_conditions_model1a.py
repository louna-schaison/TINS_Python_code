
import matplotlib.pyplot as mpl
from core.Current_functions import *
from core.Parameters import *


Gk = 6.44127490e-08
Gh=  2.32040177e-10
Imax =1.51027973e-12


submyelinK = np.full(len(times), Krest)
myelinK = np.full(len(times), Ki)
V1=np.full(len(times),-0.100)
V2=np.full(len(times),-0.06)

for i in range (len(times)-1):
    Ia = ATPase(Imax, KmK, KmNa, submyelinK[i], Ni)
    Ik = I_KIR(R, T, F, Gk, V1[i], submyelinK[i], myelinK[i])
    Ih = I_H(P_h(a, b, c, V1[i], V_half), Gh, V1[i], gamma_h, submyelinK[i], myelinK[i], Nrest, Ni)

    myelinK[i+1] = myelinK[i] - flux(Ia[0]+Ik+Ih[0], V_myelin, F, dt)
    submyelinK[i+1] = submyelinK[i] + flux(Ia[0]+Ik+Ih[0], V_submyelin, F, dt)

    V1[i+1] =V_t_dt_1(V1[i],Cm,dt,Ik,0,sum(Ih),0,sum(Ia),0,0)



submyelinK = np.full(len(times), Krest)
myelinK = np.full(len(times), Ki)

for i in range (len(times)-1):
    Ia = ATPase(Imax, KmK, KmNa, submyelinK[i], Ni)
    Ik = I_KIR(R, T, F, Gk, V2[i], submyelinK[i], myelinK[i])
    Ih = I_H(P_h(a, b, c, V2[i], V_half), Gh, V2[i], gamma_h, submyelinK[i], myelinK[i], Nrest, Ni)

    myelinK[i+1] = myelinK[i] - flux(Ia[0]+Ik+Ih[0], V_myelin, F, dt)
    submyelinK[i+1] = submyelinK[i] + flux(Ia[0]+Ik+Ih[0], V_submyelin, F, dt)

    V2[i+1] =V_t_dt_1(V2[i],Cm,dt,Ik,0,sum(Ih),0,sum(Ia),0,0)

fig, ax = mpl.subplots()

ax.plot(times, V1 * 1e3, color='lightgreen', label='V0<Vrest')
ax.plot(times, V2 * 1e3, color='green', label='V0>Vrest')

ax.set_xlabel("Time (s)")
ax.set_ylabel("Voltage (mV)")
ax.legend()
mpl.show()