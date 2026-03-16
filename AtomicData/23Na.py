#Atomic data of 23Na
#
# Sources:
# Steck, Sodium D Line Data
# #xxxx#J. E. Sansonetti, J. Phys. Chem. Ref. Data, Vol. 37, No. 1, 2008
#
#

import scipy.constants as sc
import numpy as np

a0 = sc.hbar/(sc.m_e * sc.c * sc.alpha) # bohr radius

m  = 22.989769280*1.66053906660e-27 # kg, atomic mass
I = 3.0/2.0 # Nuclear spin
gI = -0.00080461080                     #Nuclear g-factor

states = np.array(["2p63s,J=1/2","2p63p,J=1/2","2p63p,J=3/2","2p64s,J=1/2","2p63d,J=5/2","2p63d,J=3/2","2p64p,J=1/2","2p64p,J=3/2"])
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

#Ground "2p63s,J=1/2"
hyperfineA[index("2p63s,J=1/2")] = sc.h * 885.8130644e6
stateJ[index("2p63s,J=1/2")]     = 0.5
gfactors[index("2p63s,J=1/2")]   = 2.00229600
state_energies[index("2p63s,J=1/2")] = sc.h*sc.c*0.0e2

#D1 ,"2p63p,J=1/2",
hyperfineA[index("2p63p,J=1/2",)] = sc.h * 94.465e6
stateJ[index("2p63p,J=1/2",)]     = 0.5
gfactors[index("2p63p,J=1/2",)]   = 0.66581
state_energies[index("2p63p,J=1/2",)] = sc.h*sc.c*16956.17025e2

decayrates[index("2p63p,J=1/2"),index("2p63s,J=1/2")] = 6.14e+07 
decayrates[index("2p63s,J=1/2"),index("2p63p,J=1/2")] = 6.14e+07 

#D2 "2p63p,J=3/2"
hyperfineA[index("2p63p,J=3/2")] = sc.h * 18.534e6
hyperfineB[index("2p63p,J=3/2")] = sc.h * 2.72e6
hyperfineC[index("2p63p,J=3/2")] = sc.h * 0.0  #### ??
stateJ[index("2p63p,J=3/2")]     = 1.5
gfactors[index("2p63p,J=3/2")]   = 1.33420 
state_energies[index("2p63p,J=3/2")] = sc.h*sc.c*16973.36619e2

decayrates[index("2p63p,J=3/2"),index("2p63s,J=1/2")] = 6.16e+07 
decayrates[index("2p63s,J=1/2"),index("2p63p,J=3/2")] = 6.16e+07 

# "2p64s,J=1/2"
hyperfineA[index("2p64s,J=1/2")] = sc.h * 203.6e6
stateJ[index("2p64s,J=1/2")]     = 0.5
gfactors[index("2p64s,J=1/2")]   = 0.0#### ??
state_energies[index("2p64s,J=1/2")] = sc.h*sc.c*25739.999e2

decayrates[index("2p64s,J=1/2"),index("2p63p,J=3/2")] = 1.76e+07 
decayrates[index("2p63p,J=3/2"),index("2p64s,J=1/2")] = 1.76e+07 
decayrates[index("2p64s,J=1/2"),index("2p63p,J=1/2")] = 8.80e+06 
decayrates[index("2p63p,J=1/2"),index("2p64s,J=1/2")] = 8.80e+06 
decayrates[index("2p64s,J=1/2"),index("2p63s,J=1/2")] = 6.95e-04 
decayrates[index("2p63s,J=1/2"),index("2p64s,J=1/2")] = 6.95e-04 

# "2p63d,J=5/2"
hyperfineA[index("2p63d,J=5/2")] = sc.h * 0.0#### ??
hyperfineB[index("2p63d,J=5/2")] = sc.h * 0.0 #### ??
hyperfineC[index("2p63d,J=5/2")] = sc.h * 0.0  #### ??
stateJ[index("2p63d,J=5/2")]     = 2.5
gfactors[index("2p63d,J=5/2")]   = 0.0#### ??
state_energies[index("2p63d,J=5/2")] = sc.h*sc.c* 29172.837e2

decayrates[index("2p63d,J=5/2"),index("2p63p,J=3/2")] = 5.14e+07 
decayrates[index("2p63p,J=3/2"),index("2p63d,J=5/2")] = 5.14e+07 

# "2p63d,J=3/2"
hyperfineA[index("2p63d,J=3/2")] = sc.h * 0.0#### ??
hyperfineB[index("2p63d,J=3/2")] = sc.h * 0.0 #### ??
hyperfineC[index("2p63d,J=3/2")] = sc.h * 0.0  #### ??
stateJ[index("2p63d,J=3/2")]     = 1.5
gfactors[index("2p63d,J=3/2")]   = 0.0#### ??
state_energies[index("2p63d,J=3/2")] = sc.h*sc.c*29172.887e2

decayrates[index("2p63d,J=3/2"),index("2p63p,J=3/2")] = 8.57e+06 
decayrates[index("2p63p,J=3/2"),index("2p63d,J=3/2")] = 8.57e+06 
decayrates[index("2p63d,J=3/2"),index("2p63p,J=1/2")] = 4.29e+07 
decayrates[index("2p63p,J=1/2"),index("2p63d,J=3/2")] = 4.29e+07 

# "2p64p,J=1/2"
hyperfineA[index("2p64p,J=1/2")] = sc.h * 30.6e6
stateJ[index("2p64p,J=1/2")]     = 0.5
gfactors[index("2p64p,J=1/2")]   = 0.0#### ??
state_energies[index("2p64p,J=1/2")] = sc.h*sc.c*30266.99e2

decayrates[index("2p64p,J=1/2"),index("2p63d,J=3/2")] = 1.57e+05 
decayrates[index("2p63d,J=3/2"),index("2p64p,J=1/2")] = 1.57e+05 
decayrates[index("2p64p,J=1/2"),index("2p64s,J=1/2")] = 6.62e+06 
decayrates[index("2p64s,J=1/2"),index("2p64p,J=1/2")] = 6.62e+06  
decayrates[index("2p64p,J=1/2"),index("2p63s,J=1/2")] = 2.73e+06 
decayrates[index("2p63s,J=1/2"),index("2p64p,J=1/2")] = 2.73e+06  

# "2p64p,J=3/2"
hyperfineA[index("2p64p,J=3/2")] = sc.h * 0.0#### ??
hyperfineB[index("2p64p,J=3/2")] = sc.h * 0.0 #### ??
hyperfineC[index("2p64p,J=3/2")] = sc.h * 0.0  #### ??
stateJ[index("2p64p,J=3/2")]     = 1.5
gfactors[index("2p64p,J=3/2")]   = 0.0 #### ??
state_energies[index("2p64p,J=3/2")] = sc.h*sc.c*30272.58e2

decayrates[index("2p64p,J=3/2"),index("2p63d,J=3/2")] = 1.59e+04 
decayrates[index("2p63d,J=3/2"),index("2p64p,J=3/2")] = 1.59e+04 
decayrates[index("2p64p,J=3/2"),index("2p63d,J=5/2")] = 1.43e+05 
decayrates[index("2p63d,J=5/2"),index("2p64p,J=3/2")] = 1.43e+05 
decayrates[index("2p64p,J=3/2"),index("2p64s,J=1/2")] = 6.64e+06 
decayrates[index("2p64s,J=1/2"),index("2p64p,J=3/2")] = 6.64e+06  
decayrates[index("2p64p,J=1/2"),index("2p63s,J=1/2")] = 2.75e+06 
decayrates[index("2p63s,J=1/2"),index("2p64p,J=1/2")] = 2.75e+06 

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