import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#Negative Feedback Model

#Parameters
#Let p denote p53 and m denote Mdm2

ks = 2
k1 = 5
k2 = 2
K1 = 0.1
K2 = 1
dp = 1
dm = 1
n = 2
kp = 5
Kp = 1
i = 2

#Initial conditions
s0 = [0.2,0.1]

#Equations
def dpdt(t,p,m):
    return ks + kp*p**i / (Kp**i + p**i) - k1*m*p / (K1+p) - dp*p

#Here i is used in place of m from the paper to reduce confusion

def dmdt(t,p,m):
    return k2*p**n / (K2**n + p**n) - dm*m

def system(t, y):
    p0 = y[0]
    m0 = y[1]
    return [dpdt(t,p0,m0), dmdt(t,p0,m0)]

#Time scale
t = np.linspace(0,15,1000)   
 
#Solve
sols = odeint(system, y0=s0, t=t, tfirst=True)
p = sols[:,0]
m = sols[:,1]


#Plot
plt.plot(t,p, label = "p53")
plt.plot(t,m, label = "Mdm2")


plt.ylabel("Concentration")
plt.xlabel("Time")
#plt.title("Negative-positive Feedback Model")
plt.legend()
plt.show()





