
from numpy import *
from scipy.integrate import quad
from matplotlib.pyplot import *
#parameters
points1=100
ecut=1.0
hbar=1.05 #e-34Js
vf=1.0 #e6 m/s
B=10.0 #T
e=1.6 #e-19C
epsilon0=8.85 #e-12
lb=sqrt(hbar/(e*B)*10) #e-8m
Ef=2.0
    #(E1=1.29 E2=1.833 E3 = 2.245 E4=2.59 for B=5)(E1=1.83 E2=2.59 E3=3.17 E4=3.66 for B=10)
    #(E1=0.58 E2=0.82  E3=1.00 E4= 1.16 for B=1)
kf0=Ef/hbar #e8
#kf1=sqrt((Ef/hbar)**2-2/(lb)**2) #e8
hb_gamma=0.01 #e-20J
T=300.0#Kelvin
k=0.00138#e-20

#suseptbility

#def integrand (kz,hw):
#    CC=1/sqrt(2)
#    En=hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
#    Em=-hbar*vf*kz #e-34+6+8=e-20
#    deltaE=En-Em
#    nu=1.0
#    mu=sqrt((1+kz/sqrt(2/lb**2+kz**2))/2.0)
#   return 200*e**2*hbar**2*vf**2*(CC*nu*mu)**2/(8*pi**3*lb**2*epsilon0*deltaE**2) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)

#def integrand2 (kz,hw):
#    CC=0.5
#    En=hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
#    Em=-hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
#    deltaE=En-Em
#    nu=sqrt((1+kz/sqrt(2/lb**2+kz**2))/2.0)
#    mu=sqrt((1+kz/sqrt(4/lb**2+kz**2))/2.0)
#    return 200*e**2*hbar**2*vf**2*(CC*nu*mu)**2/(8*pi**3*lb**2*epsilon0*deltaE**2) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)

#conductivity

def integrand (kz,hw):# 0 to 1
    CC=1/sqrt(2)
    En=hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=-hbar*vf*kz #e-34+6+8=e-20
    deltaE=En-Em
    nu=1.0
    mu=sqrt((1+kz/sqrt(2/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def Nintegrand (kz,hw):# -1 to 0
    CC=1/sqrt(2)
    Em=-hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    En=-hbar*vf*kz #e-34+6+8=e-20
    deltaE=En-Em
    nu=1.0
    mu=sqrt((1-kz/sqrt(2/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def integrand2 (kz,hw): # -1 to 2
    CC=0.5
    En=hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=-hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1+kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1+kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def Nintegrand2 (kz,hw): # -2 to 1
    CC=0.5
    Em=-hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    En=hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1-kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1-kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def integrand3 (kz, hw): # 1 to 2 
    CC= 0.5
    En=hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1-kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1+kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def Nintegrand3 (kz, hw): # -2 to -1 
    CC= 0.5
    Em=-hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    En=-hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1+kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1-kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def integrand4 (kz, hw): # 2 to 3 
    CC= 0.5
    En=hbar*vf*sqrt(6/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1-kz/sqrt(4/lb**2+kz**2))/2.0)
    mu=sqrt((1+kz/sqrt(6/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def Nintegrand4 (kz, hw): # -3 to -2 
    CC= 0.5
    Em=-hbar*vf*sqrt(6/(lb)**2+kz**2) #e-34+6+8=e-20
    En=-hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1+kz/sqrt(4/lb**2+kz**2))/2.0)
    mu=sqrt((1-kz/sqrt(6/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def integrand5 (kz, hw): # -2 to 3 
    CC= 0.5
    En=hbar*vf*sqrt(6/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=-hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1+kz/sqrt(4/lb**2+kz**2))/2.0)
    mu=sqrt((1+kz/sqrt(6/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def Nintegrand5 (kz, hw): # -3 to 2 
    CC= 0.5
    Em=-hbar*vf*sqrt(6/(lb)**2+kz**2) #e-34+6+8=e-20
    En=hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1-kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1-kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * hb_gamma/((deltaE-hw)**2+hb_gamma**2)*FDdistri
    
result10=[]
result0_1=[]
#result12=[]
homega=[]
ab_L=[]
ab_R=[]
for i in range (1,5*points1+1):
    
    #result10.append(quad(integrand, -kf0,-kf1,args=(i/float(points1)))[0])
    #result10r=(quad(integrand, kf1,5*ecut,args=(i/float(points1)))[0])
    #result21=(quad(integrand3, -kf1,kf1, args=(i/float(points1)))[0])
    #result2_1=quad(integrand2, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    #result10[-1]+=(result10r+result21+result2_1)
    
    result10.append(quad(integrand, -5*ecut,5*ecut,args=(i/float(points1)))[0])#(-kf0,5*ecut for T=0)
    result2_1=quad(integrand2, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    result21=quad(integrand3,-5*ecut,5*ecut,args=(i/float(points1)))[0]
    result32=quad(integrand4, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    result3_2=quad(integrand5, -5*ecut, 5*ecut, args=(i/float(points1)))[0]
    result10[-1]+=(result2_1+result21+result32+result3_2)
    
    #result0_1.append(quad(integrand, kf0,5*ecut,args=(i/float(points1)))[0])
    #result1_2l=quad(integrand2, -5*ecut,-kf1,args=(i/float(points1)))[0]
    #result1_2r=quad(integrand2, kf1, 5*ecut, args=(i/float(points1)))[0]
    #result0_1[-1]+=(result1_2l+result1_2r)
    

    result0_1.append(quad(Nintegrand, -5*ecut,5*ecut,args=(i/float(points1)))[0])
    result1_2=quad(Nintegrand2, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    result_1_2=quad(Nintegrand3,-5*ecut,5*ecut,args=(i/float(points1)))[0]
    result2_3=quad(Nintegrand5, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    result_2_3=quad(Nintegrand4,-5*ecut,5*ecut,args=(i/float(points1)))[0]
    result0_1[-1]+=(result1_2+result_1_2+result2_3+result_2_3)
    
    homega.append(62.42*i/points1)
#    ab_L.append(result10[-1]*4*pi*i/float(points1)/3/hbar)#e-20+34-8
#    ab_R.append(result0_1[-1]*4*pi*i/float(points1)/3/hbar)#e+6
plot( homega, result10,homega, result0_1)
ylim([0,0.4])
show()

