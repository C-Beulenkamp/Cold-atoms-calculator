#Atomic data of 39K
#
# Sources:
# Tiecke, Properties of Potassium
# J. E. Sansonetti, J. Phys. Chem. Ref. Data, Vol. 37, No. 1, 2008
#
#

import scipy.constants as sc
import numpy as np

a0 = sc.hbar/(sc.m_e * sc.c * sc.alpha) # bohr radius

m  = 39.96399848*1.66053906660e-27 # kg, atomic mass
I =   4.0 # Nuclear spin
gI =  0.000176490                 #Nuclear g-factor
states = np.array(["3p64s,J=1/2","3p64p,J=1/2","3p64p,J=3/2","3p65s,J=1/2","3p63d,J=5/2","3p63d,J=3/2","3p65p,J=1/2","3p65p,J=3/2","3p64d,J=5/2","3p64d,J=3/2","3p66s,J=1/2","3p64f,J=5/2","3p64f,J=7/2","3p66p,J=1/2","3p66p,J=3/2"])
statenumber = len(states)

polarizability_plottable_states = [0,1,2]

def index(statename):
	return np.argwhere(states==statename)[0,0]

stateJ     = np.zeros(statenumber,dtype=float)
hyperfineA = np.zeros(statenumber)
hyperfineB = np.zeros(statenumber)
hyperfineC = np.zeros(statenumber)
state_energies = np.zeros(statenumber)
gfactors       = np.zeros(statenumber)
decayrates     = np.zeros((statenumber,statenumber))
transition_energies = np.zeros((statenumber,statenumber))
oscillatorstrengths = np.zeros((statenumber,statenumber))
dipolemoments = np.zeros((statenumber,statenumber))
transition_omegas = np.zeros((statenumber,statenumber))

#Ground "3p64s,J=1/2"
hyperfineA[index("3p64s,J=1/2")] = sc.h * (-285.731e6)
stateJ[index("3p64s,J=1/2")]     = 0.5
gfactors[index("3p64s,J=1/2")]   = 2.002295 
state_energies[index("3p64s,J=1/2")] = sc.h*sc.c*0.0

#D1 "3p64p,J=1/2"
hyperfineA[index("3p64p,J=1/2")] = sc.h * (-34.523e6)
stateJ[index("3p64p,J=1/2")]     = 0.5
gfactors[index("3p64p,J=1/2")]   = 2.0/3.0
state_energies[index("3p64p,J=1/2")] = sc.h*sc.c*12985.189386e2

decayrates[index("3p64p,J=1/2"),index("3p64s,J=1/2")] = 3.734e+07 
decayrates[index("3p64s,J=1/2"),index("3p64p,J=1/2")] = 3.734e+07 

#D2 "3p64p,J=3/2"
hyperfineA[index("3p64p,J=3/2")] = sc.h *   (-7.585e6)
hyperfineB[index("3p64p,J=3/2")] = sc.h *   (-3.44e6)
hyperfineC[index("3p64p,J=3/2")] = sc.h *   0.0  #### ??
stateJ[index("3p64p,J=3/2")]     = 1.5
gfactors[index("3p64p,J=3/2")]   = 4.0/3.0
state_energies[index("3p64p,J=3/2")] = sc.h*sc.c*13042.899700e2

decayrates[index("3p64p,J=3/2"),index("3p64s,J=1/2")] = 3.779e+07 
decayrates[index("3p64s,J=1/2"),index("3p64p,J=3/2")] = 3.779e+07 

#"3p65s,J=1/2"
hyperfineA[index("3p65s,J=1/2")] = sc.h * 0.0  #### ??
hyperfineB[index("3p65s,J=1/2")] = sc.h * 0.0  #### ??
hyperfineC[index("3p65s,J=1/2")] = sc.h * 0.0  #### ??
stateJ[index("3p65s,J=1/2")]     = 0.5
gfactors[index("3p65s,J=1/2")]   = 0.0  #### ??
state_energies[index("3p65s,J=1/2")] = sc.h*sc.c*21026.551e2  # Not isotope specific

decayrates[index("3p64p,J=1/2"),index("3p65s,J=1/2")] = 7.951e+06 
decayrates[index("3p65s,J=1/2"),index("3p64p,J=1/2")] = 7.951e+06 
decayrates[index("3p64p,J=3/2"),index("3p65s,J=1/2")] = 1.582e+07 
decayrates[index("3p65s,J=1/2"),index("3p64p,J=3/2")] = 1.582e+07 

# "3p63d,J=5/2"
hyperfineA[index("3p63d,J=5/2")] = sc.h * 0.0  #### ??
hyperfineB[index("3p63d,J=5/2")] = sc.h * 0.0  #### ??
hyperfineC[index("3p63d,J=5/2")] = sc.h * 0.0  #### ??
stateJ[index("3p63d,J=5/2")]     = 2.5
gfactors[index("3p63d,J=5/2")]   = 0.0  #### ??
state_energies[index("3p63d,J=5/2")] = sc.h*sc.c*21534.680e2  # Not isotope specific


decayrates[index("3p63d,J=5/2"),index("3p64s,J=1/2")] = 1.54e+02 
decayrates[index("3p64s,J=1/2"),index("3p63d,J=5/2")] = 1.54e+02 
decayrates[index("3p63d,J=5/2"),index("3p64p,J=3/2")] = 2.383e+07 
decayrates[index("3p64p,J=3/2"),index("3p63d,J=5/2")] = 2.383e+07 


# "3p63d,J=3/2"
hyperfineA[index("3p63d,J=3/2")] = sc.h * 0.0  #### ??
hyperfineB[index("3p63d,J=3/2")] = sc.h * 0.0  #### ??
hyperfineC[index("3p63d,J=3/2")] = sc.h * 0.0  #### ??
stateJ[index("3p63d,J=3/2")]     = 1.5
gfactors[index("3p63d,J=3/2")]   = 0.0  #### ??
state_energies[index("3p63d,J=3/2")] = sc.h*sc.c*21536.988e2  # Not isotope specific

decayrates[index("3p63d,J=3/2"),index("3p64s,J=1/2")] = 1.54e+02 
decayrates[index("3p64s,J=1/2"),index("3p63d,J=3/2")] = 1.54e+02 
decayrates[index("3p63d,J=3/2"),index("3p64p,J=3/2")] = 3.974e+06 
decayrates[index("3p64p,J=3/2"),index("3p63d,J=3/2")] = 3.974e+06 
decayrates[index("3p63d,J=3/2"),index("3p64p,J=1/2")] = 2.017e+07 
decayrates[index("3p64p,J=1/2"),index("3p63d,J=3/2")] = 2.017e+07 

# "3p65p,J=1/2"
hyperfineA[index("3p65p,J=1/2")] = sc.h * 0.0  #### ??
stateJ[index("3p65p,J=1/2")]     = 0.5
gfactors[index("3p65p,J=1/2")]   = 0.665    
state_energies[index("3p65p,J=1/2")] = sc.h*sc.c*24701.382e2  # Not isotope specific

decayrates[index("3p65p,J=1/2"),index("3p63d,J=3/2")] = 1.65e+06 
decayrates[index("3p63d,J=3/2"),index("3p65p,J=1/2")] = 1.65e+06 
decayrates[index("3p65p,J=1/2"),index("3p65s,J=1/2")] = 4.53e+06 
decayrates[index("3p65s,J=1/2"),index("3p65p,J=1/2")] = 4.53e+06 
decayrates[index("3p65p,J=1/2"),index("3p64s,J=1/2")] = 1.07e+06 
decayrates[index("3p64s,J=1/2"),index("3p65p,J=1/2")] = 1.07e+06 


# "3p65p,J=3/2"
hyperfineA[index("3p65p,J=3/2")] = sc.h * 0.0  #### ??
hyperfineB[index("3p65p,J=3/2")] = sc.h * 0.0  #### ??
hyperfineC[index("3p65p,J=3/2")] = sc.h * 0.0  #### ??
stateJ[index("3p65p,J=3/2")]     = 1.5
gfactors[index("3p65p,J=3/2")]   = 1.34 
state_energies[index("3p65p,J=3/2")] = sc.h*sc.c*24720.139e2  # Not isotope specific

decayrates[index("3p65p,J=3/2"),index("3p63d,J=3/2")] = 1.66e+05 
decayrates[index("3p63d,J=3/2"),index("3p65p,J=3/2")] = 1.66e+05 
decayrates[index("3p65p,J=3/2"),index("3p63d,J=5/2")] = 1.50e+06 
decayrates[index("3p63d,J=5/2"),index("3p65p,J=3/2")] = 1.50e+06 
decayrates[index("3p65p,J=3/2"),index("3p65s,J=1/2")] = 4.58e+06 
decayrates[index("3p65s,J=1/2"),index("3p65p,J=3/2")] = 4.58e+06 
decayrates[index("3p65p,J=3/2"),index("3p64s,J=1/2")] = 1.15e+06 
decayrates[index("3p64s,J=1/2"),index("3p65p,J=3/2")] = 1.15e+06 

# "3p64d,J=5/2"
hyperfineA[index("3p64d,J=5/2")] = sc.h * 0.0  #### ??
hyperfineB[index("3p64d,J=5/2")] = sc.h * 0.0  #### ??
hyperfineC[index("3p64d,J=5/2")] = sc.h * 0.0  #### ??
stateJ[index("3p64d,J=5/2")]     = 2.5
gfactors[index("3p64d,J=5/2")]   = 0.0  #### ??
state_energies[index("3p64d,J=5/2")] = sc.h*sc.c*27397.077e2  # Not isotope specific


decayrates[index("3p64d,J=5/2"),index("3p65p,J=3/2")] = 3.406e+06 
decayrates[index("3p65p,J=3/2"),index("3p64d,J=5/2")] = 3.406e+06 
decayrates[index("3p64d,J=5/2"),index("3p64p,J=3/2")] = 1.37e+04 
decayrates[index("3p64p,J=3/2"),index("3p64d,J=5/2")] = 1.37e+04 


# "3p64d,J=3/2"
hyperfineA[index("3p64d,J=3/2")] = sc.h * 0.0  #### ??
hyperfineB[index("3p64d,J=3/2")] = sc.h * 0.0  #### ??
hyperfineC[index("3p64d,J=3/2")] = sc.h * 0.0  #### ??
stateJ[index("3p64d,J=3/2")]     = 1.5
gfactors[index("3p64d,J=3/2")]   = 0.0  #### ??
state_energies[index("3p64d,J=3/2")] = sc.h*sc.c*27398.147e2  # Not isotope specific

decayrates[index("3p64d,J=3/2"),index("3p65p,J=3/2")] = 5.68e+05 
decayrates[index("3p65p,J=3/2"),index("3p64d,J=3/2")] = 5.68e+05 
decayrates[index("3p64d,J=3/2"),index("3p65p,J=1/2")] = 2.885e+06 
decayrates[index("3p65p,J=1/2"),index("3p64d,J=3/2")] = 2.885e+06 
decayrates[index("3p64d,J=3/2"),index("3p64p,J=3/2")] = 2.4e+03 
decayrates[index("3p64p,J=3/2"),index("3p64d,J=3/2")] = 2.4e+03 
decayrates[index("3p64d,J=3/2"),index("3p64p,J=1/2")] = 1.9e+04 
decayrates[index("3p64p,J=1/2"),index("3p64d,J=3/2")] = 1.9e+04 


# "3p66s,J=1/2"
hyperfineA[index("3p66s,J=1/2")] = sc.h * 0.0  #### ??
stateJ[index("3p66s,J=1/2")]     = 0.5
gfactors[index("3p66s,J=1/2")]   = 0.0  #### ??
state_energies[index("3p66s,J=1/2")] = sc.h*sc.c*27450.7104e2  # Not isotope specific

decayrates[index("3p66s,J=1/2"),index("3p65p,J=3/2")] = 4.955e+06 
decayrates[index("3p65p,J=3/2"),index("3p66s,J=1/2")] = 4.955e+06 
decayrates[index("3p66s,J=1/2"),index("3p65p,J=1/2")] = 1.627e+06 
decayrates[index("3p65p,J=1/2"),index("3p66s,J=1/2")] = 1.627e+06 
decayrates[index("3p66s,J=1/2"),index("3p64p,J=3/2")] = 4.956e+06 
decayrates[index("3p64p,J=3/2"),index("3p66s,J=1/2")] = 4.956e+06 
decayrates[index("3p66s,J=1/2"),index("3p64p,J=1/2")] = 2.500e+06 
decayrates[index("3p64p,J=1/2"),index("3p66s,J=1/2")] = 2.500e+06 


# "3p64f,J=5/2"
hyperfineA[index("3p64f,J=5/2")] = sc.h * 0.0  #### ??
hyperfineB[index("3p64f,J=5/2")] = sc.h * 0.0  #### ??
hyperfineC[index("3p64f,J=5/2")] = sc.h * 0.0  #### ??
stateJ[index("3p64f,J=5/2")]     = 2.5
gfactors[index("3p64f,J=5/2")]   = 0.0  #### ??
state_energies[index("3p64f,J=5/2")] = sc.h*sc.c*28127.85e2  # Not isotope specific

decayrates[index("3p64f,J=5/2"),index("3p64d,J=3/2")] = 8.3e+04 
decayrates[index("3p64d,J=3/2"),index("3p64f,J=5/2")] = 8.3e+04 
decayrates[index("3p64f,J=5/2"),index("3p64d,J=5/2")] = 5.9e+03 
decayrates[index("3p64d,J=5/2"),index("3p64f,J=5/2")] = 5.9e+03  
decayrates[index("3p64f,J=5/2"),index("3p63d,J=3/2")] = 1.463e+07  
decayrates[index("3p63d,J=3/2"),index("3p64f,J=5/2")] = 1.463e+07 
decayrates[index("3p64f,J=5/2"),index("3p63d,J=5/2")] = 1.035e+06 
decayrates[index("3p63d,J=5/2"),index("3p64f,J=5/2")] = 1.035e+06 


# "3p64f,J=7/2"
hyperfineA[index("3p64f,J=7/2")] = sc.h * 0.0  #### ??
hyperfineB[index("3p64f,J=7/2")] = sc.h * 0.0  #### ??
hyperfineC[index("3p64f,J=7/2")] = sc.h * 0.0  #### ??
stateJ[index("3p64f,J=7/2")]     = 3.5
gfactors[index("3p64f,J=7/2")]   = 0.0  #### ??
state_energies[index("3p64f,J=7/2")] = sc.h*sc.c*28127.85e2  # Not isotope specific

decayrates[index("3p64f,J=7/2"),index("3p64d,J=5/2")] = 8.9e+04 
decayrates[index("3p64d,J=5/2"),index("3p64f,J=7/2")] = 8.9e+04 
decayrates[index("3p64f,J=7/2"),index("3p63d,J=5/2")] = 1.547e+07 
decayrates[index("3p63d,J=5/2"),index("3p64f,J=7/2")] = 1.547e+07 

# "3p66p,J=1/2"
hyperfineA[index("3p66p,J=1/2")] = sc.h * 0.0  #### ??
stateJ[index("3p66p,J=1/2")]     = 0.5
gfactors[index("3p66p,J=1/2")]   = 0.6663   
state_energies[index("3p66p,J=1/2")] = sc.h*sc.c*28999.27e2  # Not isotope specific

decayrates[index("3p66p,J=1/2"),index("3p64d,J=3/2")] = 8.3e+05 
decayrates[index("3p64d,J=3/2"),index("3p66p,J=1/2")] = 8.3e+05 
decayrates[index("3p66p,J=1/2"),index("3p63d,J=3/2")] = 4.47e+05 
decayrates[index("3p63d,J=3/2"),index("3p66p,J=1/2")] = 4.47e+05 
decayrates[index("3p66p,J=1/2"),index("3p65s,J=1/2")] = 4.16e+05 
decayrates[index("3p65s,J=1/2"),index("3p66p,J=1/2")] = 4.16e+05 
decayrates[index("3p66p,J=1/2"),index("3p64s,J=1/2")] = 1.45e+05 
decayrates[index("3p64s,J=1/2"),index("3p66p,J=1/2")] = 1.45e+05  

# "3p66p,J=3/2"
hyperfineA[index("3p66p,J=3/2")] = sc.h * 0.0  #### ??
hyperfineB[index("3p66p,J=3/2")] = sc.h * 0.0  #### ??
hyperfineC[index("3p66p,J=3/2")] = sc.h * 0.0  #### ??
stateJ[index("3p66p,J=3/2")]     = 1.5
gfactors[index("3p66p,J=3/2")]   = 1.3337
state_energies[index("3p66p,J=3/2")] = sc.h*sc.c*29007.71e2  # Not isotope specific

decayrates[index("3p66p,J=3/2"),index("3p64d,J=3/2")] = 8.7e+04  
decayrates[index("3p64d,J=3/2"),index("3p66p,J=3/2")] = 8.7e+04 
decayrates[index("3p66p,J=3/2"),index("3p64d,J=5/2")] = 7.8e+05   
decayrates[index("3p64d,J=5/2"),index("3p66p,J=3/2")] = 7.8e+05 
decayrates[index("3p66p,J=3/2"),index("3p63d,J=3/2")] = 4.55e+04 
decayrates[index("3p63d,J=3/2"),index("3p66p,J=3/2")] = 4.55e+04 
decayrates[index("3p66p,J=3/2"),index("3p63d,J=5/2")] = 4.08e+05 
decayrates[index("3p63d,J=5/2"),index("3p66p,J=3/2")] = 4.08e+05 
decayrates[index("3p66p,J=3/2"),index("3p65s,J=1/2")] = 4.35e+05  
decayrates[index("3p65s,J=1/2"),index("3p66p,J=3/2")] = 4.35e+05 
decayrates[index("3p66p,J=3/2"),index("3p64s,J=1/2")] = 1.65e+05 
decayrates[index("3p64s,J=1/2"),index("3p66p,J=3/2")] = 1.65e+05  


#########


for i in np.arange(statenumber):
    for j in np.arange(statenumber):
            transition_energies[i,j] = np.absolute(state_energies[i]-state_energies[j])
            transition_omegas[i,j]   = transition_energies[i,j]/sc.hbar + 0.0001
            
            if i < j:
                    J = stateJ[i]
                    Jprime = stateJ[j]

                    oscillatorstrengths[i,j]        = decayrates[i,j]*( (2.0*Jprime+1)*2.0*np.pi*sc.epsilon_0*sc.c*sc.c*sc.c*sc.m_e  )/(  (2.0*J+1)* transition_omegas[i,j]*transition_omegas[i,j] *sc.e*sc.e )
                    dipolemoments[i,j]              = np.sqrt(3.0*sc.hbar*sc.e*sc.e*(2.0*J+1.0)*oscillatorstrengths[i,j]/(2.0*sc.m_e*transition_omegas[i,j]))
                    oscillatorstrengths[j,i]        = oscillatorstrengths[i,j]
                    dipolemoments[j,i]              = dipolemoments[i,j]
            
transition_f = transition_omegas/(2.0*np.pi)
transition_wavelengths = sc.c/transition_f




#Efimov parameters
#from PRL 111, 053202
s0       = 1.00624
# 1,-1, 33.6 G
aminus   = -830.0
etaminus = 0.204
aplus    = -aminus/4.9 #Theoretical prediction
etaplus  = 0.25 ##Guess
# 1,-1, 162.3 G
aminus   = -730.0
etaminus = 0.26
aplus    = -aminus/4.9#Theoretical prediction
etaplus  = 0.25 ##Guess