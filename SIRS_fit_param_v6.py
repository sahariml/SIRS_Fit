import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.linear_model import LinearRegression
import csv
from scipy.optimize import minimize
from scipy.signal import savgol_filter
###########
pas=1
nbday=1000
deviation_Iw=0.0001
deviation_S=0.001
deviation_I=0.001
deviation_R=0.001

window=120 #180 #120
def somme_par_pas(vecteur,pas):
    sommes_partielles = []  # Vecteur pour stocker les sommes partielles
    somme_totale = 0
    for i in range(0, len(vecteur), pas):
        somme_partielle = sum(vecteur[i:i+pas])
        sommes_partielles.append(somme_partielle)  # Stocke la somme partielle dans le vecteur
    return sommes_partielles
def orbite_SIRS(params,A,mu,S0,I0,R0,window):
    ####Paramètres
    r,beta,lambd=params
    #h=0.1
    x=S0
    y=I0
    z=R0    
    tmax=window
    S= []
    I= []
    R= []
    S_new=[]
    I_new=[]
    R_new=[]
    for k in range((tmax*pas)-1):
        if np.linalg.norm([x, y, z]) > 10e5:
            print("Divergence")
            break
        else:            
            S.append(x)
            I.append(y)
            R.append(z)
            #x1=x
            #y1=y
            #z1=z            
            S_new.append(h * (A - mu* x - lambd * x * y + beta * z))
            I_new.append(h * (lambd * x * y - (mu+ r) * y))
            R_new.append( h * (r * y - mu* z - beta * z))
            x, y, z = rec(x, y, z, h, beta, A, mu, r, lambd)
            #S_new.append(x-x1)
            #I_new.append(y-y1)
            #R_new.append(z-z1)
    #resultat=somme_par_pas(S_new,pas),somme_par_pas(I_new,pas),somme_par_pas(R_new,pas),somme_par_pas(S,pas),somme_par_pas(I,pas),somme_par_pas(R,pas)
    resultat=somme_par_pas(I_new,pas),somme_par_pas(S,pas),somme_par_pas(I,pas),somme_par_pas(R,pas)
    return resultat


def SIRS(params,A,mu,data_I_new,data_S,data_I,data_R,S0,I0,R0,window):
    # Paramètres
    r,beta,lambd=params
    m = 0
    #tmax = 12000    
    tmax = window
    x=S0
    y=I0
    z=R0    
    S= []
    I= []
    R= []
    I_new=[]
    I_new=orbite_SIRS(params,A,mu,S0,I0,R0,window)[0]
    S=orbite_SIRS(params,A,mu,S0,I0,R0,window)[1]
    I=orbite_SIRS(params,A,mu,S0,I0,R0,window)[2]
    R=orbite_SIRS(params,A,mu,S0,I0,R0,window)[3]
    #resultat=np.linalg.norm((somme_par_pas(I_new,pas) - data), ord=2)            
    #resultat=np.linalg.norm((I_new - data[:-1]), ord=2)
    #resultat=np.linalg.norm((I_new - data_I_new[:-1]), ord=2)            
    #resultat=np.linalg.norm((I_new - data_I_new[:-1]), ord=2)+np.linalg.norm((S - data_S[:-1]), ord=2)+np.linalg.norm((I - data_I[:-1]), ord=2)+np.linalg.norm((R - data_R[:-1]), ord=2)            
    #resultat=np.linalg.norm((I - data_I[:-1]), ord=2)
    #print("lenn=",len(I),len(data_I))
    resultat=np.linalg.norm((I - data_I[:len(I)]), ord=2)
    return resultat
def rec(x, y, z, h, beta, A, mu, r, lambd):
    xp = f(x, y, z, h, beta, A, mu, r, lambd)
    yp = g(x, y, z, h, beta, A, mu, r, lambd)
    zp = hh(x, y, z, h, beta, A, mu, r, lambd)
    return xp, yp, zp

def f(x, y, z, h, beta, A, mu, r, lambd):
    e = x + h * (A - mu* x - lambd * x * y + beta * z)
    return e

def g(x, y, z, h, beta, A, mu, r, lambd):
    e = y + h * (lambd * x * y - (mu+ r) * y)
    return e

def hh(x, y, z, h, beta, A, mu, r, lambd):
    e = z + h * (r * y - mu* z - beta * z)
    return e

def add_noise(data, std):
  noise = np.random.normal(0, std, size=data.shape)
  return data + noise
#############################################################
firstligne = 1
row_Cummuled_Infected=1
row_Recovered=3
row_Deaths=4
row_Infected=5
##########
### us: USA; dz: algérie; ca=canada, ....
country="pk"
with open(str(country)+'-covid.observer.csv', 'r') as file:
    reader = csv.reader(file)
    data = list(reader)
############################################
Cummuled_Infected = np.array([int(row[row_Cummuled_Infected]) for row in data[firstligne:]])
Cummuled_Recovered = np.array([int(row[row_Recovered]) for row in data[firstligne:]])    
Cummuled_Death = np.array([int(row[row_Deaths]) for row in data[firstligne:]])
Infected = np.array([int(row[row_Infected]) for row in data[firstligne:]])    
######################@
Infected_fitted=savgol_filter(Infected, 51, 3)
print("len=",len(Infected),len(Infected_fitted))
#plt.show()
#kkk=mm
#################################
data=np.zeros(nbday)
Pop=220309834
#####
h = 0.10 #fixe 
#####
#r = 0.32 #recovery rate 
#beta = 0.0083   #loss of immunity
#lambd = 0.65 #incidence bilineary
#####
A=8939/Pop 
mu=1/(365*67.7) 
#params=r,beta,lambd
#################################
####
Sw=Infected_fitted[:window]/Pop
Iw=Infected_fitted[:window]/Pop
Rw=Infected_fitted[:window]/Pop
Iw_new=Infected_fitted[:window]/Pop
#################@
S0=Sw[0]
I0=Iw[0]
R0=Rw[0]
###########Minimizing
bnds = ((0, None), (0, None), (0, None))
methods=['Nelder-Mead','Powell','CG','BFGS','L-BFGS-B','TNC','COBYLA','SLSQP']
meth=3
optimal = minimize(SIRS, [0.1, 0.1, 0.1], args=(A,mu,Iw_new,Sw,Iw,Rw,S0,I0,R0,window), method=methods[meth], bounds=bnds,options={'gtol': 1e-8, 'disp': True})
print(optimal.x)
params=optimal.x
#########################
plt.figure(1)
plt.plot(Iw, 'ko', mfc='none',markersize=3)
plt.plot(savgol_filter(orbite_SIRS(params,A,mu,S0,I0,R0,window)[2], 51, 3),'g', linewidth=1)
plt.legend(['Reported cases','Model fit'])
plt.xlabel("t (Days)")
plt.ylabel("Covid-19 Incidence case in Pakistan")
#############
plt.figure(2)
plt.plot(Iw,'g', linewidth=1)
plt.legend(['Iw'])
plt.legend(['Reported cases'])
plt.xlabel("t (Days)")
##############
plt.figure(3)
plt.plot(Infected/Pop,'c', linewidth=1)
plt.plot(Infected_fitted/Pop,'r', linewidth=1)
plt.legend(['Infected','Infected_fitted'])
plt.xlabel("t(Days)")
plt.show()

