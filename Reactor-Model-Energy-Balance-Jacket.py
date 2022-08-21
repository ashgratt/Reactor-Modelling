import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

################## DATA FILE #################

# Reactor Parameters
V = 1.36         # Reactor Volume

# Reaction Parameters
Ea = 69795          # Activation Energy
A = 19666666.67     # Arrhenius Pre-exponential Factor
H_Rx = 69795        # Heat of Reaction (should be negative but doesn't work if it is)

# Feed Properties
rho = 801      # Feed Density
cP = 3140      # Feed Specific Heat Capacity

# Feed Variables
F0 = 3*10**-4   # Feed Flow Rate In
cA0 = 8007.46   # Reactor Feed Concentration
T0 = 294.4      # Reactor Feed Temperature


# Jacket Values
UJ = 851.7      # Overall Heat Transfer Coefficient of Jacket / Coil
AJ = 23.23      # Heat Transfer Area
rhoJ = 998      # Cooling Water Density
cPJ = 4187      # Cooling Water Heat Capacity
VJ = 0.1090     # Jacket / Coil Volume
FJ0 = 0.0003925  # Cooling Water Feed Rate
TJ0 = 294.4     # Cooling Water Inlet Temperature

# Universal Values
R = 8.314       # Universal Gas Constant


####### Differential Equation Function #######
def y_prime(Y, t):
    
    # Means to Inject Disturbances #
    if t>tf/2:
        FJ = FJ0*1.5 # Example of Increasing Cooling Water Flow Rate at Half Time
        T0 = 298
    else:
        FJ = FJ0
        T0 = 294.4 # Not very elegant but can do it this way, or rename variables
   
    # Main Equations #
    cA, T, TJ = Y
    
    k = ((A * np.exp((-Ea)/(R*T))))
    
    c_prime = (F0/V)*(cA0 - cA) - k * cA
    
    T_prime = (F0/V)*(T0 - T) + ((H_Rx*k*cA)/(rho*cP) - ((UJ*AJ*(T-TJ))/(rho*cP*V)))
    
    TJ_prime = (FJ/VJ)*(TJ0 - TJ) + ((UJ*AJ*(T-TJ))/(rhoJ*cPJ*VJ))
    
    return [c_prime, T_prime, TJ_prime]
##############################################

# Define the Time Domain #
t0 = 0          # Initial Time
tf = (24)*3600    # Final Time

t = []          # Empty Time Array

# Loop to Generate a Time Array #
for i in range(t0,tf):
    t.append(i)

# Generate Initial Conditions Matrix
y0 = [cA0, T0, TJ0]

# Start in-built ODE solver #
y = odeint(y_prime, y0, t)

# Extract Columns from Array #
def column(matrix, i):
    return [row[i] for row in matrix]

cA = column(y,0)
T = column(y, 1)
TJ = column(y,2)

# Calculate Conversion
X = []

for i in range(0, len(cA)):
    conv = 1 - (cA[i]/cA0)
    X.append(conv)

# Results Processing

cA_final = cA[-1]
conversion = round(100*(1-(cA_final/cA0)), 2)

print('Final concentration = ' + str(cA_final)) 
print('Conversion = ' + str(conversion) + ' %')   

# Plotting Results

fig, ax = plt.subplots(nrows=4,ncols=1)
plt.subplot(4,1,1)
plt.plot(cA, color = 'green', label = 'Concentration')
plt.legend()

plt.subplot(4,1,2)
plt.plot(T, color = 'red', label = 'Reactor Temperature')
plt.legend()

plt.subplot(4,1,3)
plt.plot(TJ, color = 'blue', label = 'Jacket Temperature')
plt.legend()

plt.subplot(4,1,4)
plt.plot(X, color = 'black', label = 'Conversion')
plt.legend()
