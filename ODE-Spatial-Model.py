import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#BIOL 382 Project - ODE Spatial Model

#Parameters
kdph1 = 78 #min^-1
Kdph1 = 25 #microMol

kph1 = 3
Kph1 = 0.1

k1 = 10
K1 = 1.01

p0 = 0.083 #min^-1
delt0 = 0.2
Vr = 10 #volume ratio
p1 = 0.04
delt1 = 0.16 #min^-1
kSm = 0.005 # microMol min^-1
kSpm = 1
KSpm = 0.1
p2 = 0.083
delt2 = 0.0001
ktm = 1
ks = 0.015
p5 = 0.04
delt5 = 0.2
kSw = 0.003
kSpw = 1
KSpw = 0.1
p6 = 0.083
delt6 = 0.001
ktw = 1
kdph2 = 96
Kdph2 = 26
kph2 = 1
Kph2 = 0.1
E = 0.1
ATMt = 1.3

# =============================================================================
# #Possibly new permeability parameters
# p0 = 10
# p1 = 0.36
# p3 = 0
# p4 = 0
# p5 = 10
# p6 = 0.36
# 
# =============================================================================

#Initial conditions
s0 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]

#Equations
def du0(t,u0,v0,u1,u3,u4,u5):
    return kdph1*u5*u3/(Kdph1+u3) - k1*u1*u0/(K1+u0) -  kph1*u4*u0/(Kph1+u0) - p0*Vr*(u0-v0)

def du1(t,u1,v1):
    return -p1*Vr*(u1-v1) - delt1*u1

def du2(t,u2,u3):
    return kSm + kSpm*(u3**4)/(KSpm**4 +u3**4) - p2*Vr*u2 - delt2*u2

def du3(t,u0,u3,u4,u5):
    return kph1*u4*u0/(Kph1+u0) - kdph1*u5*u3/(Kdph1+u3)

def du4(t,u4,u5):
    return kph2*E*(ATMt-u4)/(Kph2+0.5*(ATMt-u4)) - 2*kdph2*u5*(u4**2)/(Kdph2+u4**2)

def du5(t,u5,v5):
    return p5*Vr*v5 - delt5*u5

def du6(t,u3,u6):
    return kSw + kSpw*(u3**4)/(KSpw**4 +u3**4) - p6*Vr*u6 - delt6*u6

def dv0(t,u0,v0,v1):
    return ks - k1*v1*v0/(K1+v0) - p0*(v0-u0) - delt0*v0

def dv1(t,u1,v1,v2):
    return ktm*v2 - p1*(v1-u1) - delt1*v1

def dv2(t,u2,v2):
    return p2*u2 - ktm*v2 - delt2*v2

def dv3(t):
    return 0

def dv4(t):
    return 0

def dv5(t,v5,v6):
    return ktw*v6 - p5*v5 - delt5*v5

def dv6(t,u6,v6):
    return p6*u6 - ktw*v6 - delt6*v6


def system(t, y):
    u0 = y[0]
    u1 = y[1]
    u2 = y[2]
    u3 = y[3]
    u4 = y[4]
    u5 = y[5]
    u6 = y[6]
    v0 = y[7]
    v1 = y[8]
    v2 = y[9]
    v3 = y[10]
    v4 = y[11]
    v5 = y[12]
    v6 = y[13]
    
    return [du0(t,u0,v0,u1,u3,u4,u5),du1(t,u1,v1),du2(t,u2,u3),\
            du3(t,u0,u3,u4,u5),du4(t,u4,u5),du5(t,u5,v5),du6(t,u3,u6),\
                dv0(t,u0,v0,v1),dv1(t,u1,v1,v2),dv2(t,u2,v2),\
                    dv3(t),dv4(t),dv5(t,v5,v6),dv6(t,u6,v6)]

#Time scale
t = np.linspace(0,300,1000)   
 
#Solve
sols = odeint(system, y0=s0, t=t, tfirst=True)

u0 = sols[:,0] #p53(n)
u1 = sols[:,1] #Mdm2(n)
u2 = sols[:,2] #Mdm2 mRNA (n)
u3 = sols[:,3] #p53-p (n)
u4 = sols[:,4] #ATM-p (n)
u5 = sols[:,5] #Wip1 (n)
u6 = sols[:,6] #Wip1mRNA (n)
v0 = sols[:,7] #p53(n)
v1 = sols[:,8] #Mdm2(n)
v2 = sols[:,9] #Mdm2 mRNA (n)
v3 = sols[:,10] #p53-p (n)
v4 = sols[:,11] #ATM-p (n)
v5 = sols[:,12] #Wip1 (n)
v6 = sols[:,13] #Wip1mRNA (n)

#Plot
plt.plot(t,u0, label = "p53 (n)")
plt.plot(t,u1, label = "Mdm2 (n)")
#plt.plot(t,u2, label = "Mdm2-mRNA (n)")
plt.plot(t,u3, label = "p53-p (n)")
plt.plot(t,u4, label = "ATM-p (n)")
plt.plot(t,u5, label = "Wip1 (n)")
#plt.plot(t,u6, label = "Wip1-mRNA (n)")

plt.plot(t,v0, label = "p53 (c)")
plt.plot(t,v1, label = "Mdm2 (c)")
#plt.plot(t,v2, label = "Mdm2-mRNA (c)")
plt.plot(t,v3, label = "p53-p (c)")
plt.plot(t,v4, label = "ATM-p (c)")
plt.plot(t,v5, label = "Wip1 (c)")
#plt.plot(t,v6, label = "Wip1-mRNA (c)")

plt.ylabel("Concentration (dimensionless)")
plt.xlabel("Time (dimensionless)")
#plt.title("Negative Feedback Model")
plt.legend()
plt.show()


plt.plot(t,u0+u3, label = "Total p53")
plt.plot(t,u1, label = "Mdm2")
#plt.plot(t,u3, label = "p53-p")
plt.plot(t,u4, label = "ATM-p")
plt.plot(t,u5, label = "Wip1")
plt.ylabel("Concentration (dimensionless)")
plt.xlabel("Time (dimensionless)")
plt.title("Nucleus")
plt.legend()
plt.show()

plt.plot(t,v0, label = "p53 (inactive)")
plt.plot(t,v1, label = "Mdm2")
#plt.plot(t,v3, label = "p53-p")
plt.plot(t,v4, label = "ATM-p")
plt.plot(t,v5, label = "Wip1")
plt.ylabel("Concentration (dimensionless)")
plt.xlabel("Time (dimensionless)")
plt.title("Cytoplasm")
plt.legend()
plt.show()








