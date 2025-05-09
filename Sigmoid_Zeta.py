import numpy as np
import matplotlib.pyplot as plt

# Defining variables
K = 20    #max k
T = 20     #max time
Nk = 10000    
Nt = 10000
zetamax = 0.03  #Max zeta
dk = K / Nk  
dt = T / Nt        
E = np.zeros((Nt, Nk))
k = np.linspace(0.5, Lk, Nk)  
k_j = 2     # Jean's wavenumber
zetaf = 0   # final zeta at k_j integration (ensure continuity)

# initial conditions
count = 0
for ki in k:
    if ki < 1:
        E[0,count] = ki**4
    elif ki<10:
        E[0,count] = ki**(-1)
    elif ki>= 10:
        E[0,count] =  (4.64)*ki**(-5/3)
    count+=1

# integration loop
for ti in range(0, Nt-1):
    for ki in range(0, Nk):
        
        # define zeta
        if k[ki] < k_j:
            zeta = np.tanh(k[ki]) / (np.pi/2) * zetamax
            zetaf  = zeta 
        if k[ki] > k_j:
            zeta = (1 - np.tanh(1*(k[ki]-k_j))) * zetaf


        # integrate
        if ki == 0: #Left boundary
            E[ti+1, ki] = E[ti, ki] - dt * zeta*((E[ti, ki+1] - E[ti, ki]) / (0.5*dk) - (4 / k[ki]) * E[ti, ki]) 
        if 0<ki<Nk-1: 
            E[ti+1, ki] = E[ti, ki] - dt * zeta*((E[ti, ki+1] - E[ti, ki-1]) / dk - (4 / k[ki]) * E[ti, ki])
        if ki == Nk-1: # Right boundary
            E[ti+1, ki] = E[ti, ki] - dt * zeta*((E[ti, ki] - E[ti, ki-1]) / (0.5*dk) - (4 / k[ki]) * E[ti, ki])

#plot
times = [0,5,10,15]
timesp = ["t = "+ str(int(time/dt)) for time in times]
plt.figure(figsize=(20, 12))
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 40
plt.rcParams['mathtext.fontset'] = 'stix'
colors = ["indigo","blueviolet", "violet","palevioletred"]
index = 0
style ="solid"
for time in times:
    plt.loglog(k,E[int(time/dt),:],linestyle = style,color=colors[index],linewidth=3)
    style = "dashed"
    index+=1
plt.legend(timesp) 
plt.xlabel(r'$k$')
plt.ylabel(r'$E_M(k)$')
plt.tick_params(which='minor', length=9)
plt.tick_params(which='major', length=9,width=3)
plt.show()


