#Atomic data of 87 Cs
#
# Sources:
# Steck, Rubidium 87 D Line Data
# J. E. Sansonetti, J. Phys. Chem. Ref. Data 35, 301 (2006)
#
#

import scipy.constants as sc
import numpy as np

a0 = sc.hbar/(sc.m_e * sc.c * sc.alpha)

m  = 86.909180520*1.66053906660e-27 # kg, atomic mass
I = 3.0/2.0 # Nuclear spin
gI = -0.0009951414             #Nuclear g-factor

states = np.array(["4p65s,J=1/2","4p65p,J=1/2","4p65p,J=3/2","4p66p,J=1/2","4p66p,J=3/2"])
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

#Ground "4p65s,J=1/2"
hyperfineA[index("4p65s,J=1/2")] = sc.h*sc.c* 0.113990236053642*100.0
stateJ[index("4p65s,J=1/2")]     = 0.5
gfactors[index("4p65s,J=1/2")]   = 2.00254032
state_energies[index("4p65s,J=1/2")] = sc.h*sc.c*0.0e2

#D1 ,"4p65p,J=1/2",
hyperfineA[index("4p65p,J=1/2",)] = sc.h*sc.c* 0.01365e2
stateJ[index("4p65p,J=1/2",)]     = 0.5
gfactors[index("4p65p,J=1/2",)]   = 0.66590
state_energies[index("4p65p,J=1/2",)] = sc.h*sc.c*12578.950e2

decayrates[index("4p65p,J=1/2"),index("4p65s,J=1/2")] = 3.61e7
decayrates[index("4p65s,J=1/2"),index("4p65p,J=1/2")] = 3.61e7

#D2 "4p65p,J=3/2"
hyperfineA[index("4p65p,J=3/2")] = sc.h*sc.c* 0.00282590e2
hyperfineB[index("4p65p,J=3/2")] = sc.h*sc.c* 0.00041684e2
hyperfineC[index("4p65p,J=3/2")] = sc.h * 0.0  #### ??
stateJ[index("4p65p,J=3/2")]     = 1.5
gfactors[index("4p65p,J=3/2")]   = 1.3340
state_energies[index("4p65p,J=3/2")] = sc.h*sc.c* 12816.545e2

decayrates[index("4p65p,J=3/2"),index("4p65s,J=1/2")] = 3.81e7
decayrates[index("4p65s,J=1/2"),index("4p65p,J=3/2")] = 3.81e7


# "4p66p,J=1/2"
hyperfineA[index("4p66p,J=1/2")] = sc.h*sc.c*0.0044217e2
stateJ[index("4p66p,J=1/2")]     = 0.5
gfactors[index("4p66p,J=1/2")]   = 0.0#### ??
state_energies[index("4p66p,J=1/2")] = sc.h*sc.c*23715.081e2

decayrates[index("4p66p,J=1/2"),index("4p65s,J=1/2")] = 1.50e6
decayrates[index("4p65s,J=1/2"),index("4p66p,J=1/2")] = 1.50e6

# "4p66p,J=3/2"
hyperfineA[index("4p66p,J=3/2")] = sc.h*sc.c*0.0009240e2
hyperfineB[index("4p66p,J=3/2")] = sc.h*sc.c*0.000132e2
hyperfineC[index("4p66p,J=3/2")] = sc.h * 0.0  #### ??
stateJ[index("4p66p,J=3/2")]     = 1.5
gfactors[index("4p66p,J=3/2")]   = 0.0 #### ??
state_energies[index("4p66p,J=3/2")] = sc.h*sc.c*23792.591e2

decayrates[index("4p66p,J=1/2"),index("4p65s,J=1/2")] = 1.77e6
decayrates[index("4p65s,J=1/2"),index("4p66p,J=1/2")] = 1.77e6

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