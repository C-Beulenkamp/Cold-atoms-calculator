#Atomic data of 133 Cs
#
# Sources:
# Steck, Cesium D Line Data
# J. E. Sansonetti, J. Phys. Chem. Ref. Data 38, 761 (2009)
# ???? Hyperfine Interaction, Zeeman and Stark Effects for Excited States in Cesium, G Belin et al 1976 Phys. Scr. 14 39

import scipy.constants as sc
import numpy as np

a0 = sc.hbar/(sc.m_e * sc.c * sc.alpha)

m  = 132.905451931*1.66053906660e-27    # kg, atomic mass
I = 7.0/2.0                             # Nuclear spin
gI = -0.00039885395                     #Nuclear g-factor

states = np.array(["5p6(1S)6s,J=1/2","5p6(1S)6p,J=1/2","5p6(1S)6p,J=3/2","5p6(1S)5d,J=3/2","5p6(1S)5d,J=5/2","5p6(1S)7s,J=1/2","5p6(1S)7p,J=1/2","5p6(1S)7p,J=3/2","5p6(1S)6d,J=3/2","5p6(1S)6d,J=5/2"])
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

#Ground
hyperfineA[index("5p6(1S)6s,J=1/2")] = sc.h *   2.2981579425e9
hyperfineB[index("5p6(1S)6s,J=1/2")] = sc.h *   0.0
stateJ[index("5p6(1S)6s,J=1/2")]     = 0.5
gfactors[index("5p6(1S)6s,J=1/2")]   = 2.00254032
state_energies[index("5p6(1S)6s,J=1/2")] = 0.0

#D1
hyperfineA[index("5p6(1S)6p,J=1/2")] = sc.h *   291.9309e6
hyperfineB[index("5p6(1S)6p,J=1/2")] = sc.h *   0.0
stateJ[index("5p6(1S)6p,J=1/2")]     = 0.5
gfactors[index("5p6(1S)6p,J=1/2")]   = 0.66590
state_energies[index("5p6(1S)6p,J=1/2")] = sc.h*sc.c*11178.2686e2

decayrates[index("5p6(1S)6s,J=1/2"),index("5p6(1S)6p,J=1/2")] = 2.862e7
decayrates[index("5p6(1S)6p,J=1/2"),index("5p6(1S)6s,J=1/2")] = 2.862e7

#D2
hyperfineA[index("5p6(1S)6p,J=3/2")] = sc.h *   50.28825e6
hyperfineB[index("5p6(1S)6p,J=3/2")] = sc.h *   (-0.4940e6)
hyperfineC[index("5p6(1S)6p,J=3/2")] = sc.h *   (0.56e3)
stateJ[index("5p6(1S)6p,J=3/2")]     = 1.5
gfactors[index("5p6(1S)6p,J=3/2")]   = 1.3340
state_energies[index("5p6(1S)6p,J=3/2")] = sc.h*sc.c*11732.3079e2

decayrates[index("5p6(1S)6s,J=1/2"),index("5p6(1S)6p,J=3/2")] = 3.279e7
decayrates[index("5p6(1S)6p,J=3/2"),index("5p6(1S)6s,J=1/2")] = 3.279e7

#5d, 3/2
hyperfineA[index("5p6(1S)5d,J=3/2")] = sc.h *   48.75e6
hyperfineB[index("5p6(1S)5d,J=3/2")] = sc.h *   0.1e6
stateJ[index("5p6(1S)5d,J=3/2")]     = 1.5
state_energies[index("5p6(1S)5d,J=3/2")] = sc.h*sc.c*14499.2584e2

decayrates[index("5p6(1S)6s,J=1/2"),index("5p6(1S)5d,J=3/2")] = 0
decayrates[index("5p6(1S)5d,J=3/2"),index("5p6(1S)6s,J=1/2")] = 0
decayrates[index("5p6(1S)6p,J=1/2"),index("5p6(1S)5d,J=3/2")] = 9.13e5
decayrates[index("5p6(1S)5d,J=3/2"),index("5p6(1S)6p,J=1/2")] = 9.13e5
decayrates[index("5p6(1S)6p,J=3/2"),index("5p6(1S)5d,J=3/2")] = 1.07e5
decayrates[index("5p6(1S)5d,J=3/2"),index("5p6(1S)6p,J=3/2")] = 1.07e5

#5d, 5/2
hyperfineA[index("5p6(1S)5d,J=5/2")] = sc.h *   (-21.24e6)
hyperfineB[index("5p6(1S)5d,J=5/2")] = sc.h *   0.2e6
stateJ[index("5p6(1S)5d,J=5/2")]     = 2.5
state_energies[index("5p6(1S)5d,J=5/2")] = sc.h*sc.c*14596.8423e2

decayrates[index("5p6(1S)6s,J=1/2"),index("5p6(1S)5d,J=5/2")] = 22.2
decayrates[index("5p6(1S)5d,J=5/2"),index("5p6(1S)6s,J=1/2")] = 22.2
decayrates[index("5p6(1S)6p,J=3/2"),index("5p6(1S)5d,J=5/2")] = 7.81e5
decayrates[index("5p6(1S)5d,J=5/2"),index("5p6(1S)6p,J=3/2")] = 7.81e5

#7s
hyperfineA[index("5p6(1S)7s,J=1/2")] = sc.h *   545.90e6
hyperfineB[index("5p6(1S)7s,J=1/2")] = sc.h *   0.0
stateJ[index("5p6(1S)7s,J=1/2")]     = 0.5
state_energies[index("5p6(1S)7s,J=1/2")] = sc.h*sc.c*18535.529e2

#7p 1/2
hyperfineA[index("5p6(1S)7p,J=1/2")] = sc.h *   94.35e6
hyperfineB[index("5p6(1S)7p,J=1/2")] = sc.h *   0.0
stateJ[index("5p6(1S)7p,J=1/2")]     = 0.5
state_energies[index("5p6(1S)7p,J=1/2")] = sc.h*sc.c*21765.35e2

decayrates[index("5p6(1S)6s,J=1/2"),index("5p6(1S)7p,J=1/2")] = 7.93e5
decayrates[index("5p6(1S)7p,J=1/2"),index("5p6(1S)6s,J=1/2")] = 7.93e5

#7p 3/2
hyperfineA[index("5p6(1S)7p,J=3/2")] = sc.h *   16.609e6
hyperfineB[index("5p6(1S)7p,J=3/2")] = sc.h *   0.0
stateJ[index("5p6(1S)7p,J=3/2")]     = 1.5
state_energies[index("5p6(1S)7p,J=3/2")] = sc.h*sc.c*21946.396e2

decayrates[index("5p6(1S)6s,J=1/2"),index("5p6(1S)7p,J=3/2")] = 1.84e6
decayrates[index("5p6(1S)7p,J=3/2"),index("5p6(1S)6s,J=1/2")] = 1.84e6

#6d 3/2
hyperfineA[index("5p6(1S)6d,J=3/2")] = sc.h *   16.34e6
hyperfineB[index("5p6(1S)6d,J=3/2")] = sc.h *   (-0.1e6)
stateJ[index("5p6(1S)6d,J=3/2")]     = 1.5
state_energies[index("5p6(1S)6d,J=3/2")] = sc.h*sc.c*22588.8210e2

#6d 5/2
hyperfineA[index("5p6(1S)6d,J=5/2")] = sc.h *   (-4.66e6)
hyperfineB[index("5p6(1S)6d,J=5/2")] = sc.h *   0.9e6
stateJ[index("5p6(1S)6d,J=5/2")]     = 2.5
state_energies[index("5p6(1S)6d,J=5/2")] = sc.h*sc.c*22631.6863

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

#Cesium D2 line 
J_S = stateJ[0]
J_D1 = stateJ[1]
J_D2 = stateJ[2]
f_D2      = 351.72571850e12 # Hz, transition frequency
lamb_D2   = sc.c/f_D2       # m, transition vacuum wavelength
omega_D2  = 2.0*np.pi* f_D2 # Hz, angular transition frequency
E_D2      = sc.hbar*omega_D2   # J, transition energy
tau_D2    = 30.473e-9       # s, lifetime
gamma_D2  = 1.0/tau_D2      # 1/s  ( = 2 pi 5.2227 MHz), decay rate
vr_D2     = 3.5225e-2       # m/s, recoil velocity
omegar_D2 = 2.0* np.pi * 2.0663e3 # Hzs, recoil energy
Tr_D2     = 198.34e-9       # K, Recoil remperature
TD_D2     = 125.0e-6         # K, Doppler temperature
Isat_D2   = 2.7059 * (0.001/(0.01*0.01)) #mW/cm^2, saturation intensity

A_D2      = 3.279e7   # 1/s , transition probability, from J. Phys. Chem. Ref. Data 38, 761 (2009)
fosc_D2   = A_D2*( (2.0*J_D2+1)*2.0*np.pi*sc.epsilon_0*sc.c*sc.c*sc.c*sc.m_e  )/(  (2.0*J_S+1)* omega_D2*omega_D2 *sc.e*sc.e ) # Oscillator strength
d_D2      = np.sqrt(3.0*sc.hbar*sc.e*sc.e*(2.0*J_S+1.0)*fosc_D2/(2.0*sc.m_e*omega_D2)) #Transition dipole moment

#Cesium D1 line 
f_D1      = 335.116048807e12            # Hz, transition frequency
lamb_D1   = sc.c/f_D1       # m, transition vacuum wavelength
omega_D1  = 2.0*np.pi* f_D1             # Hz, angular transition frequency
E_D1      = sc.hbar*omega_D1            # J, transition energy
tau_D1    = 34.894e-9                   # s, lifetime
gamma_D1  = 1.0/tau_D2                  # 1/s  ( = 2 pi 5.2227 MHz), decay rate
vr_D1     = 3.3561e-2                   # m/s, recoil velocity
omegar_D1 = 2.0* np.pi * 1.8758e3       # Hzs, recoil energy
Tr_D1     = 180.05e-9                   # K, Recoil remperature

A_D1      = 2.862e7                     # 1/s , transition probability, from J. Phys. Chem. Ref. Data 38, 761 (2009)
fosc_D1   = A_D1*( (2.0*J_D1+1)*2.0*np.pi*sc.epsilon_0*sc.c*sc.c*sc.c*sc.m_e  )/(  (2.0*J_S+1)* omega_D1*omega_D1 *sc.e*sc.e )  # Oscillator strength
d_D1      = np.sqrt(3.0*sc.hbar*sc.e*sc.e*(2.0*J_S+1.0)*fosc_D1/(2.0*sc.m_e*omega_D1)) #Transition dipole moment


#Efimov parameters
aminus = -850.0
aplus  = 1060.0
etaminus = 0.06
etaplus  = 0.06
s0 = 1.00624

