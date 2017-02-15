from numpy import *
from scipy.integrate import quad
from matplotlib.pyplot import *
#parameters
n_bands=4
hbar=1.05 #e-34Js
vf=1.0 #e6 m/s
B=10.0 #T
e=1.6 #e-19C
epsilon0=8.85 #e-12
lb=sqrt(hbar/(e*B)*10) #e-8m
Ef=1.28
    #(E1=1.29 E2=1.833 E3 = 2.245 E4=2.59 for B=5)(E1=1.83 E2=2.59 E3=3.17 E4=3.66 for B=10)
    #(E1=0.58 E2=0.82  E3=1.00 E4= 1.16 for B=1)
kf0=Ef/hbar #e8
#kf1=sqrt((Ef/hbar)**2-2/(lb)**2) #e8

k=0.00138#e-20

def getbands(n,kz):
    if n==0: return -hbar*vf*kz
    elif n<0: return -hbar*vf*sqrt(kz**2-2*n/(lb**2))
    else: return hbar*vf*sqrt(kz**2+2*n/(lb**2))


def getdensity(E,Ef,T):
    return 1/(exp((E-Ef)/k/T)+1)

def getfermiL(Ef,l):
    return [Ef/0.016]*l



kz = arange(-5,5,0.01)
for i in range(-n_bands,n_bands+1):
    plot(kz, getbands(i,kz)/0.016)
E = arange(-5,5,0.01)
plot(getdensity(E,Ef,100),E/0.016)
plot(getdensity(E,Ef,200),E/0.016)
plot(getdensity(E,Ef,300),E/0.016)
#Ef_array = [Ef/0.016]*len(kz)
#plot(kz, Ef_array,"--")
#Ef_array1 = getfermiL(0.8,len(kz))
#Ef_array2 = getfermiL(1.28, len(kz))
#Ef_array3 = getfermiL(1.92, len(kz))
#plot(kz,Ef_array1,"--",kz,Ef_array2,"--",kz,Ef_array3,"--") 

show()
