
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt



def lattice(N, init='random'):
    """This function creates a lattice of size NxN with random values in the range [-1,1].

    Args:
        N (Positive Integer): It is the size of the lattice which you want to initialize.
        init (str, optional): It is the configuration with which you want to initialize the lattice. Defaults to 'random'. If you want to initialize the lattice with all sites as 1, then you can use 'up'.

    Returns:
        2D Numpy array: size NxN.
    """
    if init == 'random':
        lat = np.random.random((N,N)) 
        lat[lat >= 0.5] = 1
        lat[lat < 0.5] = -1
    elif init == 'up':
        lat = np.ones((N,N))  
    return lat


def calcMag(lattice):
    """This function is used to calculate the magnetization value by summing up all the spin values in the lattice. This calculation neglects the constant value that is multiplied while actual magnetization calculation.

    Args:
        lattice (2D Array): The lattice array which contains the spin values, for which you want to calculate the magnetization.

    Returns:
        Integer: Magnetization value.
    """
    mag = np.sum(lattice)
    return mag
 

def calcEnergy(lattice, mode='hv'):
    """Calculate lattice energy.

    Args:
        lattice (2D array): The spin lattice.
        mode (str): Interaction mode, can be:
            'hv' - horizontal and vertical (default)
            'd' - diagonal only
            'hvd' - horizontal, vertical, and diagonal
    """
    energy = 0
    if mode == 'hv':
        for i in range(len(lattice)):
            for j in range(len(lattice)):
                energy += -lattice[i,j] * (lattice[i, (j+1)%len(lattice)] + lattice[i, (j-1)%len(lattice)] + lattice[(i+1)%len(lattice), j] + lattice[(i-1)%len(lattice), j]) 
    elif mode == 'd':
        for i in range(len(lattice)):
            for j in range(len(lattice)):
                energy += -lattice[i,j] * (lattice[(i+1)%len(lattice), (j+1)%len(lattice)] + lattice[(i-1)%len(lattice), (j-1)%len(lattice)] + lattice[(i-1)%len(lattice), (j+1)%len(lattice)] + lattice[(i+1)%len(lattice), (j-1)%len(lattice)])
    elif mode == 'hvd':
        # horizontal/vertical
        energy += calcEnergy(lattice, 'hv')  
        # diagonal
        energy += calcEnergy(lattice, 'd')
    return energy / 2




def monteCarloSteps(lattice, beta):
    """Performs Monte Carlo steps on a lattice at a given temperature.
    
    This function carries out Monte Carlo steps on a lattice to move it towards 
    thermal equilibrium at a given inverse temperature beta. At each site, a 
    spin is randomly flipped with probability determined by the energy cost of 
    flipping that spin.
    
    Args:
        lattice (2D numpy array): The spin lattice to simulate
        beta (float): The inverse temperature 1/T
    
    Returns:
        lattice (2D numpy array): The lattice after Monte Carlo steps are performed
    
    """
    

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




def simulation(lattice_size, eqsteps=1000, mcsteps=10000):
    """This function is used for running the simulation. 
    Simulation here refers to plotting the graphs of observables like magnetization, energy etc. against temperature.
    

    Args:
        lattice_size (Integer): Lattice initialized for the simulation.
        eqsteps (int, optional): steps for which we want out lattice to go through monteCarloSteps function in order to reach an equilibrium state at that particular temperature. Defaults to 1000.
        mcsteps (int, optional): steps for which we want out lattice to go through monteCarloSteps function before calculating the observable values. Defaults to 10000.
    """
    nt = 45
    E1 = np.zeros(nt)
    M1 = np.zeros(nt)
    Temp = np.linspace(0.5, 5, 45)
    T = 1/Temp
    for t in range(len(T)):
        a1 = lattice(lattice_size)
        for i in range(eqsteps):
            a1 = monteCarloSteps(a1, T[t])
            
        for i in range(mcsteps):
            a1 = monteCarloSteps(a1, T[t]) 
            E1[t] += calcEnergy(a1)
            M1[t] += calcMag(a1)
            
    E1 = E1/(mcsteps*lattice_size*lattice_size)
    M1 = M1/(mcsteps*lattice_size*lattice_size)
    plt.plot(Temp, E1,  marker='o', color='IndianRed', label = "L=2")
    plt.xlabel("Temperature (T)", fontsize=20)
    plt.ylabel("Energy ", fontsize=20)
    plt.legend()
    plt.show()
    plt.plot(Temp, abs(M1),  marker='o', color='RoyalBlue',label = "L=2")
    plt.xlabel("Temperature (T)", fontsize=20); 
    plt.ylabel("Magnetization ", fontsize=20)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    simulation(8, 2**8, 2**9)
