import numpy as np
import matplotlib.pyplot as plt

# Initial Parameters
K = 20    # Maximum wavenumber 
T = 20  # Maximum time
Nk = 10000   
Nt = 10000
zeta = 0.01 
dk = K / Nk  
dt = T / Nt       

E = np.zeros((Nt, Nk))  
k = np.linspace(dk + 0.5, K, Nk)  #avoids zero


# Initial conditions
count = 0
for ki in k:
    if ki < 1:
        E[0,count] = ki**4
    elif ki<10:
        E[0,count] = ki**(-1)
    elif ki>= 10:
        E[0,count] =  (4.64)*ki**(-5/3)
    count+=1
    

# Integration loop
for ti in range(0, Nt-1):
    for ki in range(0, Nk):
        if ki == 0:
            E[ti+1, ki] = E[ti, ki] - dt * zeta*((E[ti, ki+1] - E[ti, ki]) / (0.5*dk) - (4 / k[ki]) * E[ti, ki])
        elif 0<ki<Nk-1:
            E[ti+1, ki] = E[ti, ki] - dt * zeta*((E[ti, ki+1] - E[ti, ki-1]) / dk - (4 / k[ki]) * E[ti, ki])
        elif ki == Nk-1:
            E[ti+1, ki] = E[ti, ki] - dt * zeta*((E[ti, ki] - E[ti, ki-1]) / (0.5*dk) - (4 / k[ki]) * E[ti, ki])




# Plot Parameters
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


