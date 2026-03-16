#Atomic data of 6Li
#
# Sources:
# Gehm, Properties of 6Li
# Measurement of isotope shifts, fine and hyperfine structure splittings of the lithium D lines, Eur. Phys. J. D 22, 159–162 (2003)
#

import scipy.constants as sc
import numpy as np

a0 = sc.hbar/(sc.m_e * sc.c * sc.alpha) # bohr radius

m  = 6.0151214*1.66053906660e-27 # kg, atomic mass
I = 1.0 # Nuclear spin
gI = -0.0004476540                   #Nuclear g-factor

states = np.array(["1s22s,J=1/2","1s22p,J=1/2","1s22p,J=3/2","1s23s,J=1/2","1s23p,J=1/2","1s23p,J=3/2"])
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

#Ground "1s22s,J=1/2"
hyperfineA[index("1s22s,J=1/2")] = sc.h * 152.1368407e6
stateJ[index("1s22s,J=1/2")]     = 0.5
gfactors[index("1s22s,J=1/2")]   =  2.0023010
state_energies[index("1s22s,J=1/2")] = sc.h*sc.c*0.0

#D1 "1s22p,J=1/2"
hyperfineA[index("1s22p,J=1/2")] = sc.h * 17.386e6
stateJ[index("1s22p,J=1/2")]     = 0.5
gfactors[index("1s22p,J=1/2")]   = 0.6668
state_energies[index("1s22p,J=1/2")] = sc.h*sc.c*14903.66e2

decayrates[index("1s22p,J=1/2"),index("1s22s,J=1/2")] = 3.6890e+07 
decayrates[index("1s22s,J=1/2"),index("1s22p,J=1/2")] = 3.6890e+07 

#D2 "1s22p,J=3/2"
hyperfineA[index("1s22p,J=3/2")] = -sc.h * 1.155e6 
hyperfineB[index("1s22p,J=3/2")] = -sc.h * 0.106 
hyperfineC[index("1s22p,J=3/2")] =  sc.h * 0.0  #### ??
stateJ[index("1s22p,J=3/2")]     = 1.5
gfactors[index("1s22p,J=3/2")]   = 1.335
state_energies[index("1s22p,J=3/2")] = sc.h*sc.c*14904.00e2

decayrates[index("1s22p,J=3/2"),index("1s22s,J=1/2")] = 3.6891e+07 
decayrates[index("1s22s,J=1/2"),index("1s22p,J=3/2")] = 3.6891e+07 

# "1s23s,J=1/2"
hyperfineA[index("1s23s,J=1/2")] = sc.h *   0.0  #### ??
hyperfineB[index("1s23s,J=1/2")] = sc.h *   0.0  #### ??
hyperfineC[index("1s23s,J=1/2")] = sc.h *   0.0  #### ??
stateJ[index("1s23s,J=1/2")]     = 0.5
gfactors[index("1s23s,J=1/2")]   = 0.0  #### ??
state_energies[index("1s23s,J=1/2")] = sc.h*sc.c*27206.12e2

decayrates[index("1s23s,J=1/2"),index("1s22p,J=1/2")] = 1.1156e+07 
decayrates[index("1s22p,J=1/2"),index("1s23s,J=1/2")] = 1.1156e+07 
decayrates[index("1s23s,J=1/2"),index("1s22p,J=3/2")] = 2.2309e+07 
decayrates[index("1s22p,J=3/2"),index("1s23s,J=1/2")] = 2.2309e+07 

# "1s23p,J=1/2"
hyperfineA[index("1s23p,J=1/2")] = sc.h *   0.0  #### ??
stateJ[index("1s23p,J=1/2")]     = 0.5  
gfactors[index("1s23p,J=1/2")]   = 0.0  #### ??
state_energies[index("1s23p,J=1/2")] = sc.h*sc.c*30925.38e2

decayrates[index("1s23p,J=1/2"),index("1s23s,J=1/2")] = 3.738e+06 
decayrates[index("1s23s,J=1/2"),index("1s23p,J=1/2")] = 3.738e+06 
decayrates[index("1s23p,J=1/2"),index("1s22s,J=1/2")] = 1.002e+06 
decayrates[index("1s22s,J=1/2"),index("1s23p,J=1/2")] = 1.002e+06 

# "1s23p,J=3/2"
hyperfineA[index("1s23p,J=3/2")] = sc.h *   0.0  #### ??
hyperfineB[index("1s23p,J=3/2")] = sc.h *   0.0  #### ??
hyperfineC[index("1s23p,J=3/2")] = sc.h *   0.0  #### ??
stateJ[index("1s23p,J=3/2")]     = 1.5
gfactors[index("1s23p,J=3/2")]   = 0.0  #### ??
state_energies[index("1s23p,J=3/2")] = sc.h*sc.c*30925.38e2

decayrates[index("1s23p,J=3/2"),index("1s23s,J=1/2")] = 3.737e+06 
decayrates[index("1s23s,J=1/2"),index("1s23p,J=3/2")] = 3.737e+06 
decayrates[index("1s23p,J=3/2"),index("1s22s,J=1/2")] = 1.002e+06 
decayrates[index("1s22s,J=1/2"),index("1s23p,J=3/2")] = 1.002e+06 


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
s0       = 1.00624
aminus   = 1
etaminus = 1
aplus    = 1 
etaplus  = 1