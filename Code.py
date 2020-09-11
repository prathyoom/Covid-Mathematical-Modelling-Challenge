# Model 2 code

# Importing libraries
from matplotlib import pyplot as plt

# Defining constants
Beta = 2.0 # Infection constant
Alpha = 0.0005 # Susceptible constant
GammaR = 0.9 # Recovery Constant
GammaD = 0.1 # Death Constant
N = 100 # Total number of people

# Setting Initial conditions
simulation = 0
Limit = 100

Susceptible = 99 # Number of susceptible people
Infected = 1 # Number of infected people
Recovery = 0 # Number of people temporarily immmune to covid    
Death = 0 # Number of people deceased by covid

# Creating lists to obtain data points
S = []
I = []
R = []
D = []
T = [] # Time

# Performing the discrete simulation
while(simulation <= Limit):
    # Adding data points to plot
    S.append(Susceptible)
    I.append(Infected)
    R.append(Recovery)
    D.append(Death)
    T.append(simulation)

    # Change in number of susceptible people
    dS = ((-1 * Beta * Infected * Susceptible) / N) + (Alpha * Recovery)
    # Change in number of Infected people
    dI = ((Beta * Infected * Susceptible)/N) - (GammaR * Infected) - (GammaD * Infected)
    # Change in number of Recovered people
    dR = (GammaR * Infected) - (Alpha* Recovery)
    # Change in number of deceased people
    dD = (GammaD * Infected)

    # Updating the values
    Susceptible += dS
    Infected += dI 
    Recovery += dR
    Death += dD

    simulation+= 1

"""for i in range(len(S)):
    print(S[i], end=' ')
    print(I[i], end=' ')
    print(R[i], end=' ')
    print(D[i], end=' ')
    print(S[i]+I[i]+R[i]+D[i])"""

plt.plot(T,S)
plt.plot(T,I)
plt.plot(T,R)
plt.plot(T,D)   

plt.legend(["Susceptible", "Infected", "Recovered", "Deceased"])
plt.figtext( .65, 0.8, ["Beta = ", Beta])
plt.figtext( .65, .75, ["Alpha = ", Alpha])
plt.figtext( .65, .7, ["GammaR = ", GammaR])
plt.figtext( .65, .65, ["GammaD = ", GammaD])

plt.show()