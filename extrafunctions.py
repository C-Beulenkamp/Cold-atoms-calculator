import scipy.constants as sc
import importlib
import numpy as np
from sympy.physics.wigner import wigner_3j,wigner_6j

muB  = sc.value(u'Bohr magneton') # Bohr magneton


def QuietWigner6j(j1,j2,j3,j4,j5,j6):
    try:
        return wigner_6j(j1,j2,j3,j4,j5,j6).evalf()
    except:
        return 0.0
    

def BreitRabi(b,m,s,atom):
    atomdata = importlib.import_module("AtomicData."+atom)  
    ehfs = atomdata.hyperfineA[0]*(atomdata.I+0.5)
    x = (atomdata.gfactors[0]-atomdata.gI)*muB*b/ehfs
    if m == atomdata.I+0.5:
        return ehfs*atomdata.I/(2.0*atomdata.I+1.0) + 0.5*(atomdata.gfactors[0] + 2.0*atomdata.I*atomdata.gI)*muB*b
    if 0 == atomdata.I+0.5 + m:
        return ehfs*atomdata.I/(2.0*atomdata.I+1.0) - 0.5*(atomdata.gfactors[0] + 2.0*atomdata.I*atomdata.gI)*muB*b
    return -0.5*ehfs/(2.0*atomdata.I+1.0) + atomdata.gI*muB*m*b + s*0.5*ehfs*np.sqrt(1+4.0*m*x/(2.0*atomdata.I+1.0) + x*x)
    
def scientificnotation(value):
	return str('{:.4e}'.format(value))
    
    
def gPolyLogarithm(z,gamma):
    output = z.astype(np.float64)
    index = 1
    
    while index < 10000.0:
        index += 1
        nextterm = np.power(z,index)/np.power(index,gamma)
        output += nextterm
        #if (nextterm/output < 1.0e-6):
        #    break
        
    return output

def superscript(num):
    tmp = ["\u2070","\u00b9","\u00b2", "\u00b3","\u2074","\u2075","\u2076","\u2077","\u2078","\u2079"]
    return tmp[num]


def subscript(num):
    tmp = ["\u2080","\u2081","\u2082", "\u2083","\u2084","\u2085","\u2086","\u2087","\u2088","\u2089"]
    return tmp[num]
