import tkinter as tk                    
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
import GUIfunctions as GUI
from traps import *
from polarizability import *
from lattices import *
#from examples import *
from magneticfields import *
from feshbach import *
from atomicclouds import *
from other import *
#from scatteringthresholds import *

#Main window 
root = tk.Tk()
root.title("Cold atom calculator")
root.geometry("1190x870")
photo = tk.PhotoImage(file = "icon.png")
root.iconphoto(False, photo)
root.resizable(False,False)

s = ttk.Style()
s.configure('My.TFrame', background='#541a1a')
#s.theme_use('alt')
#print(s.theme_names()) #print a list of available styles.

#Define the tabs
tabControl = ttk.Notebook(root,)

EMtab = ttk.Frame(tabControl)
Potentialstab = ttk.Frame(tabControl)
AtomicScatteringtab = ttk.Frame(tabControl)

tabControl.add(EMtab, text ='E.M fields')
tabControl.add(Potentialstab, text ='Potentials')
tabControl.add(AtomicScatteringtab, text ='Atomic scattering')

AtomicCloudsTab = ttk.Frame(tabControl)
OtherTab = ttk.Frame(tabControl)

tabControl.add(AtomicCloudsTab, text ='Atomic clouds')
tabControl.add(OtherTab, text ='Other')
#tabControl.add(tab1,text ='functions')
tabControl.pack(expand = 1, fill ="both")

emtabcontrol = ttk.Notebook(EMtab,)
HyperfineTab = ttk.Frame(emtabcontrol)
PolarizabilityTab = ttk.Frame(emtabcontrol)
emtabcontrol.add(HyperfineTab, text ='Hyperfine + Zeeman')
emtabcontrol.add(PolarizabilityTab, text ='AC Polarizability')
emtabcontrol.pack(expand = 1, fill ="both")

potentialtabcontrol = ttk.Notebook(Potentialstab,)
TrapsTab = ttk.Frame(Potentialstab)
LatticeTab = ttk.Frame(Potentialstab)
potentialtabcontrol.add(TrapsTab, text ='Magnetic + optical trapping')
potentialtabcontrol.add(LatticeTab, text ='Optical lattices')
potentialtabcontrol.pack(expand = 1, fill ="both")

AtomicScatteringtabcontrol = ttk.Notebook(AtomicScatteringtab,)
ScatteringDataTab = ttk.Frame(AtomicScatteringtab)
ScatteringThresholdsTab = ttk.Frame(AtomicScatteringtab)
AtomicScatteringtabcontrol.add(ScatteringDataTab, text ='Scattering data')
#AtomicScatteringtabcontrol.add(ScatteringThresholdsTab, text ='Scattering thresholds')
AtomicScatteringtabcontrol.pack(expand = 1, fill ="both")


traps(TrapsTab,row=0,column=0)
lattices(LatticeTab,row=0,column=0)
magnetic_fields(HyperfineTab,row=0,column=0)
feshbach_spectrum(ScatteringDataTab,row=0,column=0)
atomicclouds(AtomicCloudsTab,row=0,column=0)
polarizability(PolarizabilityTab,row=0,column=0)
GaussianBeamCalculator(OtherTab,row=0,column=0)
FiberCouplingCalculator(OtherTab,row=1,column=0)
#scatteringthresholds(ScatteringThresholdsTab,row=0,column=0)
#DiffractionGratingCalculator(OtherTab,row=0,column=1)
#CavityCalculator(OtherTab,row=0,column=2)
#GUI_Examples(tab1)


####
root.mainloop()  
