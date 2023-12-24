""" This is the code in which we are equilibrating the lattices before run the montecarlo steps on it.
    We have defined a lattice randomly and then randomly assigning the spins as +1 and -1.
    """



import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt



def lattice(N, init='random'):
    if init == 'random':
        lat = np.random.random((N,N)) 
        lat[lat >= 0.5] = 1
        lat[lat < 0.5] = -1
    elif init == 'up':
        lat = np.ones((N,N))  
    return lat

def calcMag(lattice):
    mag = np.sum(lattice)
    return mag


def calcEnergy(lattice):
    energy = 0
    for a in range(len(lattice)):
        for b in range(len(lattice)):
            s = lattice[(a+1)%len(lattice),b] + lattice[(a-1)%len(lattice),b] + lattice[a,(b+1)%len(lattice)] + lattice[a,(b-1)%len(lattice)]
            energy += -s*lattice[a,b]
    return energy/2



def monteCarloSteps(lattice, beta):
    for a in range(len(lattice)):
        for b in range(len(lattice)):
            s =  lattice[a, b]
            nb = lattice[(a+1)%len(lattice),b] + lattice[a,(b+1)%len(lattice)] + lattice[(a-1)%len(lattice),b] + lattice[a,(b-1)%len(lattice)]
            cost = 2*s*nb
            if cost < 0:
                s *= -1
            elif rand() < np.exp(-cost*beta):
                s *= -1
            lattice[a, b] = s
    return lattice


n1 = 2
n2= 4
n3 = 8
n4 = 16
n5 = 50
n6 = 100
nt = 45
T       = np.linspace(5, 0.5, nt); 
E1,E2,E3,E4,E5,E6, M1,M2,M3,M4,M5,M6  = np.zeros(nt), np.zeros(nt),np.zeros(nt),np.zeros(nt),np.zeros(nt),np.zeros(nt),np.zeros(nt),np.zeros(nt),np.zeros(nt),np.zeros(nt),np.zeros(nt),np.zeros(nt)

eqsteps = 2**8
mcsteps = 2**10

for t in range(nt):
    a1 = lattice(n1)
    a2 = lattice(n2)
    a3 = lattice(n3)
    a4 = lattice(n4)
    a5 = lattice(n5)
    a6 = lattice(n6)
    
    iT = 1/T[t]
    Ene1  = Ene2 = Ene3 = Ene4 = Ene5 = Ene6 = Mag1 = Mag2 = Mag3 = Mag4 = Mag5 = Mag6  = 0
    
    for i in range(eqsteps):         # equilibrate
        a1 = monteCarloSteps(a1, iT)           # Monte Carlo moves
        a2 = monteCarloSteps(a2, iT)           # Monte Carlo moves
        a3 = monteCarloSteps(a3, iT)           # Monte Carlo moves
        a4 = monteCarloSteps(a4, iT)           # Monte Carlo moves
        a5 = monteCarloSteps(a5, iT)           # Monte Carlo moves
        a6 = monteCarloSteps(a6, iT)           # Monte Carlo moves

    for i in range(mcsteps):
        a1 = monteCarloSteps(a1, iT)           # Monte Carlo moves
        a2 = monteCarloSteps(a2, iT)           # Monte Carlo moves
        a3 = monteCarloSteps(a3, iT)           # Monte Carlo moves
        a4 = monteCarloSteps(a4, iT) 
        a5 = monteCarloSteps(a5,iT)           
        a6 = monteCarloSteps(a6,iT)      
        
        ene1 = calcEnergy(a1)
        ene2 = calcEnergy(a2)
        ene3 = calcEnergy(a3)
        ene4 = calcEnergy(a4)
        ene5 = calcEnergy(a5)
        ene6 = calcEnergy(a6)

        
        
        mag1 = calcMag(a1)        
        mag2 = calcMag(a2)        
        mag3 = calcMag(a3)        
        mag4 = calcMag(a4)        
        mag5 = calcMag(a5)        
        mag6 = calcMag(a6)        
       

        Ene1 += ene1
        Ene2 += ene2
        Ene3 += ene3
        Ene4 += ene4
        Ene5 += ene5
        Ene6 += ene6
        
        
        Mag1 += mag1
        Mag2 += mag2
        Mag3 += mag3
        Mag4 += mag4
        Mag5 += mag5
        Mag6 += mag6
        
   
    
    E1[t] = Ene1/(mcsteps*n1*n1)
    E2[t] = Ene2/(mcsteps*n2*n2)
    E3[t] = Ene3/(mcsteps*n3*n3)
    E4[t] = Ene4/(mcsteps*n4*n4)
    E5[t] = Ene5/(mcsteps*n5*n5)
    E6[t] = Ene6/(mcsteps*n6*n6)
    
    M1[t] = Mag1/(mcsteps*n1*n1) 
    M2[t] = Mag2/(mcsteps*n2*n2) 
    M3[t] = Mag3/(mcsteps*n3*n3) 
    M4[t] = Mag4/(mcsteps*n4*n4) 
    M5[t] = Mag5/(mcsteps*n5*n5) 
    M6[t] = Mag6/(mcsteps*n6*n6) 

     
    
f = plt.figure(figsize=(18, 20)); #  


sp =  f.add_subplot(2, 1, 1 );
plt.plot(T, E1,   marker='o', color='IndianRed',label = "L=2")
plt.plot(T, E2,  marker='o', color='RoyalBlue', label = "L=4")
plt.plot(T, E3,  marker='o', color='Magenta', label = "L=8")
plt.plot(T, E4,  marker='o', color='Orange', label = "L=16")
plt.plot(T, E5,  marker='o', color='yellow', label = "L=50")
plt.plot(T, E6,  marker='o', color='green', label = "L=100")
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Energy ", fontsize=20);         plt.axis('tight');
plt.legend()
# plt.show()


sp =  f.add_subplot(2, 1, 2 );
plt.plot(T, abs(M1),  marker='o', color='RoyalBlue',label = "L=2")
plt.plot(T, abs(M2),  marker='o', color='IndianRed',label = "L=4")
plt.plot(T, abs(M3),  marker='o', color='Magenta',label = "L=8")
plt.plot(T, abs(M4),  marker='o', color='cyan',label = "L=16")
plt.plot(T, abs(M5),  marker='o', color='yellow',label = "L=50")
plt.plot(T, abs(M6),  marker='o', color='green',label = "L=100")

plt.xlabel("Temperature (T)", fontsize=20); 
plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');
plt.legend()

plt.show()
