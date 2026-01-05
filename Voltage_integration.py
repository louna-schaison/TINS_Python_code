from core.Current_functions import *
import matplotlib.pyplot as mpl


# Intialisation

submyelinK = np.full(len(times), Krest)
myelinK = np.full(len(times), Ki)
extK = np.full(len(times),Krest)
atpase_current_K=np.full(len(times),0.)
IKIR_s=np.full(len(times),0.0)
IHCN_s=np.full(len(times),0.0)
IATP_s=np.full(len(times),0.0)

IKIR_e=np.full(len(times),0.0)
IHCN_e=np.full(len(times),0.0)
IATP_e=np.full(len(times),0.0)

EK_s=np.full(len(times),0.)
EK_e=np.full(len(times),0.)

V=np.full(len(times),(R*T/F)*np.log((gamma*Krest+Nrest)/(gamma*Ki+Ni)))

# Over time modelisation

for i in range(len(times)-1):

    submyelinK[i+1]=submyelinK[i] #on aligne la concentration au pas de temps précedent
    myelinK[i+1]=myelinK[i]

    #AP efflux

    if sum(APtimes < (times[i] - dt)) < sum(APtimes < times[i]):  #meaning that the t+dt is reaching the timing of a new AP
        submyelinK[i + 1] += APefflux  # addition of K+ efflux from the neuronal activity


    #interface submyelin

     #Les courants
    atpase_current=ATPase_2(Imax/20,KmK,KmNa,submyelinK[i],Ni)
    IATP_s[i]=sum(atpase_current)

   # ikir= I_KIR(R,T,F,Gk/20,V[i],submyelinK[i],myelinK[i])
   # IKIR_s[i]=ikir
    ihcn= I_H(P_openH(a,b,c,V[i],V_half),Gh/20,V[i],gamma,submyelinK[i],myelinK[i])
    IHCN_s[i]=sum(ihcn)


   #Les concentrations K+
    atpase_current_K=atpase_current[1]
    myelinK[i+1] -=  atpase_current_K*dt/(F*V_myelin)
    submyelinK[i+1]+=atpase_current_K*dt/(F*V_submyelin)

   # myelinK[i+1]-= ikir * dt / (V_myelin * F)
   # submyelinK[i+1]+= ikir * dt / (V_submyelin * F)

    ihcn_k= ihcn[0]
    myelinK[i+1]-= ihcn_k* dt / (V_myelin * F)
    submyelinK[i+1]+= ihcn_k* dt / (V_submyelin * F)

    #interface ext

     #Les courants
    atpase_current=ATPase_2(Imax/20,KmK,KmNa,Krest,Ni)
    IATP_e[i]=sum(atpase_current)

   # ikir= I_KIR(R,T,F,Gk/20,V[i],Krest,myelinK[i])
   # IKIR_e[i]=ikir
    ihcn= I_H(P_openH(a,b,c,V[i],V_half),Gh/20,V[i],gamma,Krest,myelinK[i])
    IHCN_e[i]=sum(ihcn)


   #Les concentrations K+
    atpase_current_K=atpase_current[1]
    myelinK[i+1] -=  atpase_current_K*dt/(F*V_myelin)

    #myelinK[i+1]-= ikir * dt / (V_myelin * F)

    ihcn_k= ihcn[0]
    myelinK[i+1]-= ihcn_k* dt / (V_myelin * F)


    V[i + 1] = V[i] + dt / Cm * ( - IHCN_s[i] -IATP_s[i]-IHCN_e[i] -IATP_e[i])


#To avoid the last value of the I to be null (since were going until  times -1)

IHCN_e[-1]=IHCN_e[-2]
IHCN_s[-1]=IHCN_s[-2]
IKIR_e[-1] = IKIR_e[-2]
IKIR_s[-1] = IKIR_s[-2]



#Plotting the results

#Plotting Vm

fig,ax=mpl.subplots()
mpl.plot(times, V*1000,color='darkgreen')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Vm (mV)")
mpl.title("Evolution of the myelin membrane potential overtime")

# Plotting the [K+]

fig, ax1 = mpl.subplots()
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
line1 = ax1.plot(times, submyelinK, color='blue', label="Submyelin")
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("[K+] (M)",color='blue')
ax2 = ax1.twinx()
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
line2 = ax2.plot(times, myelinK, color='darkblue', label="Intracellular")
ax2.set_ylabel("[K+] (M)",color='darkblue')
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='lower right')
mpl.title("Evolution of the intracellular [K+] and submyelin [K+] overtime")

# Plotting the currents

#fig,ax=mpl.subplots()
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
#mpl.plot(times, 1e12*IKIR_e,color='red',label="External interface")
#mpl.plot(times, 1e12*IKIR_s,color='darkred', label='Submyelin interface')
#ax.set_xlabel("Time (s)")
#ax.set_ylabel("I (pA)")
#mpl.title("KIR current over time ")
#mpl.legend()

#fig,ax=mpl.subplots()
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
#mpl.plot(times, 1e12*IHCN_e,color='red',label="External interface")
#mpl.plot(times,1e12* IHCN_s,color='darkred', label='Submyelin interface')
#ax.set_xlabel("Time (s)")
#ax.set_ylabel("I (pA)")
#mpl.title("HCN current over time ")
#mpl.legend()


mpl.show()
