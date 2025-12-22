from matplotlib.pyplot import subplots

import Parameters
from Parameters import *
import numpy as np
import matplotlib.pyplot as mpl
import pandas as pd

D=pd.read_csv('Traces-control_Fig2B.csv')
D=D.iloc[:, :3]
D.iloc[614,0]=5.196359
D.iloc[1067,0]=7.552666
D=D.to_numpy()
D=D.astype(float)
D=D[D[:,0]>=2.8]
D=D[D[:,0]<=7]
D=D.T
D[1:2]=((-D[1:2]+0.2516)/-0.004037) *1e-3

fig,ax= mpl.subplots()
mpl.plot(D[0],D[1], color = 'lightgreen')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("Time (s)")
ax.set_ylabel("V (V)")
mpl.title('Estimated V from fluorescence traces')
mpl.show()

