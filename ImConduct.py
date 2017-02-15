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
Ef=1
kf0=Ef/hbar #e8
#kf1=sqrt((Ef/hbar)**2-2/(lb)**2) #e8
hb_gamma=0.008 #e-20J
eb=6
d=20 #e-6
T=10
0.0#Kelvin0.8
k=0.00138#e=-20


def Imintegrand (kz,hw):# 0 to 1
    CC=1/sqrt(2)
    En=hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=-hbar*vf*kz #e-34+6+8=e-20
    deltaE=En-Em
    nu=1.0
    mu=sqrt((1+kz/sqrt(2/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * (hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def NImintegrand (kz,hw):# -1 to 0
    CC=1/sqrt(2)
    Em=-hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    En=-hbar*vf*kz #e-34+6+8=e-20
    deltaE=En-Em
    nu=1.0
    mu=sqrt((1-kz/sqrt(2/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) *(hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri


def Imintegrand2 (kz,hw): # -1 to 2
    CC=0.5
    En=hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=-hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1+kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1+kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * (hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri
def NImintegrand2 (kz,hw): # -2 to 1
    CC=0.5
    Em=-hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    En=hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1-kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1-kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * (hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def Imintegrand3 (kz, hw): # 1 to 2 
    CC= 0.5
    En=hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1-kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1+kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * (hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri#e4 conductivity
def NImintegrand3 (kz, hw): # -2 to -1 
    CC= 0.5
    Em=-hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    En=-hbar*vf*sqrt(2/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1+kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1-kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * (hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def Imintegrand4 (kz, hw): # 2 to 3 
    CC= 0.5
    En=hbar*vf*sqrt(6/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1-kz/sqrt(4/lb**2+kz**2))/2.0)
    mu=sqrt((1+kz/sqrt(6/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * (hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def NImintegrand4 (kz, hw): # -3 to -2 
    CC= 0.5
    Em=-hbar*vf*sqrt(6/(lb)**2+kz**2) #e-34+6+8=e-20
    En=-hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1+kz/sqrt(4/lb**2+kz**2))/2.0)
    mu=sqrt((1-kz/sqrt(6/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * (hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def Imintegrand5 (kz, hw): # -2 to 3 
    CC= 0.5
    En=hbar*vf*sqrt(6/(lb)**2+kz**2) #e-34+6+8=e-20
    Em=-hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1+kz/sqrt(4/lb**2+kz**2))/2.0)
    mu=sqrt((1+kz/sqrt(6/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * (hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri

def NImintegrand5 (kz, hw): # -3 to 2 
    CC= 0.5
    Em=-hbar*vf*sqrt(6/(lb)**2+kz**2) #e-34+6+8=e-20
    En=hbar*vf*sqrt(4/(lb)**2+kz**2) #e-34+6+8=e-20
    deltaE=En-Em
    nu=sqrt((1-kz/sqrt(2/lb**2+kz**2))/2.0)
    mu=sqrt((1-kz/sqrt(4/lb**2+kz**2))/2.0)
    FDdistri=1/(exp((Em-Ef)/k/T)+1)-1/(exp((En-Ef)/k/T)+1)
    return 2*e**2*hbar*vf**2*(CC*nu*mu)**2/(2*pi**2*lb**2*deltaE) * (hw-deltaE)/((deltaE-hw)**2+hb_gamma**2)*FDdistri



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
Imresult10=[]
result0_1=[]
Imresult0_1=[]
#result12=[]
homega=[]
eta=[]
totab=[]

for i in range (1,5*points1+1):
        #result10.append(quad(integrand, -kf0,-kf1,args=(i/float(points1)))[0])
    #result10r=(quad(integrand, kf1,5*ecut,args=(i/float(points1)))[0])
    #result21=(quad(integrand3, -kf1,kf1, args=(i/float(points1)))[0])
    #result2_1=quad(integrand2, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    #result10[-1]+=(result10r+result21+result2_1)
    
    result10.append(quad(integrand, -5*ecut,5*ecut,args=(i/float(points1)))[0])
    result2_1=quad(integrand2, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    result21=quad(integrand3,-5*ecut,5*ecut,args=(i/float(points1)))[0]
    result32=quad(integrand4, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    result3_2=quad(integrand5, -5*ecut, 5*ecut, args=(i/float(points1)))[0]
    result10[-1]+=(result2_1+result21+result32+result3_2)

    Imresult10.append(quad(Imintegrand, -5*ecut,5*ecut,args=(i/float(points1)))[0])
    Imresult2_1=quad(Imintegrand2, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    Imresult21=quad(Imintegrand3,-5*ecut,5*ecut,args=(i/float(points1)))[0]
    Imresult32=quad(Imintegrand4, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    Imresult3_2=quad(Imintegrand5, -5*ecut, 5*ecut, args=(i/float(points1)))[0]
    Imresult10[-1]+=(Imresult2_1+Imresult21+Imresult32+Imresult3_2)

    
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

    
    Imresult0_1.append(quad(NImintegrand, -5*ecut,5*ecut,args=(i/float(points1)))[0])
    Imresult1_2=quad(NImintegrand2, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    Imresult_1_2=quad(NImintegrand3,-5*ecut,5*ecut,args=(i/float(points1)))[0]
    Imresult2_3=quad(NImintegrand5, -5*ecut,5*ecut,args=(i/float(points1)))[0]
    Imresult_2_3=quad(NImintegrand4,-5*ecut,5*ecut,args=(i/float(points1)))[0]
    Imresult0_1[-1]+=(Imresult1_2+Imresult_1_2+Imresult2_3+Imresult_2_3)

    homega.append(i/points1)
#    ab_L.append(result10[-1]*4*pi*i/float(points1)/3/hbar)#e-20+34-8
#    ab_R.append(result0_1[-1]*4*pi*i/float(points1)/3/hbar)#e+6

#convert form conductivity to dielectric constand and index of refraction
indexn10=[]
indexn0_1=[]
ab_L=[]
ab_R=[]
ReDiel10=[]
ImDiel10=[]
ReDiel0_1=[]
ImDiel0_1=[]
trans10=[]
trans0_1=[]
ReConduct10=[]
ReConduc0_1=[]
for i in range(len(result10)):
    ReDiel10.append(eb-(100*hbar/(homega[i]*epsilon0))*Imresult10[i])
    ReDiel0_1.append(eb-(100*hbar/(homega[i]*epsilon0))*Imresult0_1[i])
    ImDiel10.append((100*hbar/(homega[i]*epsilon0))*result10[i])
    ImDiel0_1.append((100*hbar/(homega[i]*epsilon0))*result0_1[i])


    
    indexn10.append(sqrt(0.5*(sqrt(ImDiel10[i]**2+ReDiel10[i]**2)+ReDiel10[i])))
    indexn0_1.append(sqrt(0.5*(sqrt(ImDiel0_1[i]**2+ReDiel0_1[i]**2)+ReDiel0_1[i])))
    ab_L.append(result10[i]/epsilon0/3/indexn10[i]*100)  #e6
    ab_R.append(result0_1[i]/epsilon0/3/indexn0_1[i]*100) #e6
    trans10.append(exp(-2*d*ab_L[i]))
    trans0_1.append(exp(-2*d*ab_R[i]))
    eta.append((trans0_1[i]-trans10[i])/(trans10[i]+trans0_1[i]))

    homega[i]=62.42*homega[i]
    totab.append(ab_L[i]+ab_R[i])
    
plot(homega,result10,homega,result0_1)
ylim([0,0.35])
show()
plot(homega,Imresult10,homega,Imresult0_1)
show()
#plot(homega,ReDiel10,homega,ReDiel0_1)
#ylim([-10,20])
#show()
#plot(homega,ImDiel10,homega,ImDiel0_1)
#show()
#plot(homega, indexn10, homega, indexn0_1)
#ylim([0,5])
#show()
#plot(homega, ab_L, homega, ab_R)
#ylim([0,1])
#show()
#plot(homega,trans10,homega,trans0_1)
#ylim([0,1])
#show()
#plot(homega,eta)
#ylim([0,1])
#show()

