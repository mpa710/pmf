import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

# Initial Parameters
K = 20
T = 20
Nk = 10000
Nt = 10000
dk = K / Nk
dt = T / Nt
zeta_0 = 0.01

k = np.linspace(dk + 0.5, Lk, Nk)  
E = np.zeros((Nt, Nk))

# Initial Conditions
count = 0
for ki in k:
    if ki < 1:
        E[0,count] = ki**4
    elif ki<10:
        E[0,count] = ki**(-1)
    elif ki>= 10:
        E[0,count] =  (4.64)*ki**(-5/3)
    count+=1


main = np.zeros_like(k) # Main diagonal of banded matrix
lower = np.zeros_like(k) # Lower Diagnoal of banded matrix

# Updated main and lower diagonals
for i in range(len(k)):
    lambda_i = zeta_0 * dt / dk
    main[i] = 1 + lambda_i - 4 * zeta_0 * dt / k[i]
    lower[i] = - lambda_i
    if i ==  len(k):
        lower[i] = 0

# Banded Matrix
bandedmatrix = np.zeros((2,Nk)) 
bandedmatrix[0, :] = lower 
bandedmatrix[1, :] = main        

# Solve
for n in range(Nt-1):
    E[n+1,:] = solve_banded((1, 0), np.vstack([main, lower]), E[n,:]) # one non zero lower diagonal, zero non zero upper diagonal

# Plot parameters
plt.figure(figsize=(20, 12))
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 40
plt.rcParams['mathtext.fontset'] = 'stix'
colors = ["indigo","blueviolet", "violet","palevioletred"]
index = 0
style ="solid"
times = [0,5,10,15]
timesp = ["t = "+ str(int(time/dt)) for time in times]
for time in times:
    plt.loglog(k,E[int(time/dt),:],linestyle = style,color=colors[index],linewidth=3)
    style = "dashed"
    index+=1
plt.legend(timesp) 
plt.xlabel(r'$k$')
plt.ylabel(r'$E_M(k)$')
plt.ylim(1e-4,2)
plt.tick_params(which='minor', length=9)
plt.tick_params(which='major', length=9,width=3)
plt.show()

