import tkinter as tk                    
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
import GUIfunctions as GUI
import importlib
import scipy
import scipy.constants as sc
from AngularMomentumOperatorClass import Operator,AngMomOperators
from sympy import N
from scipy import optimize
from polarizability import compute_polarizability
from extrafunctions import *
from math import floor
import glob, os


sharedfx = 0.0
sharedfy = 0.0
sharedfz = 0.0
sharedax = 0.0
shareday = 0.0
sharedaz = 0.0
muB  = sc.value(u'Bohr magneton') # Bohr magneton
polarizabilityunit = sc.e*sc.e* sc.physical_constants['Bohr radius'][0]*sc.physical_constants['Bohr radius'][0]/(sc.physical_constants['Hartree energy'][0])


def set_get_frequencies(option,*args):
    global sharedfx
    global sharedfy
    global sharedfz
    global sharedax
    global shareday
    global sharedaz
    if option == 'set':
        sharedfx, sharedfy, sharedfz, sharedax, shareday, sharedaz = args
    
    if option == 'get':
        return sharedfx, sharedfy, sharedfz, sharedax, shareday, sharedaz


def traps(tab,row,column):
    #os.chdir("AtomicData")
    #Constants
    
    #The frame on the left
    settingsframe = GUI.MakeFrame(frame=tab,row=0,column=0,sticky='ns')
    ###########
    #Optical dipole trap parameters
    beams_frame = GUI.MakeLabelFrame(frame=settingsframe,row=1,column=0,width=100,height=100,text='Optical traps')
    #Magnetic trap parameters
    magnetictraps_frame = GUI.MakeLabelFrame(frame=settingsframe,row=2,column=0,width=50,height=50,text='Magnetic traps')
    #Atomic parameters
    atoms_frame = GUI.MakeLabelFrame(frame=settingsframe,row=0,column=0,width=50,height=50,text='Atomic details')
    ###########
    
    
    ########################################################################################################################################################################################
    ########################################################################################################################################################################################
    #############################################################################################The frame on the right
    ########################################################################################################################################################################################
    ########################################################################################################################################################################################
    outputframe = GUI.MakeLabelFrame(frame=tab,row=0,column=1,width=150,height=150,text='Output')
    ###########
    ###############Canvas and figure for plotting
    figure = plt.Figure(figsize=(7.9,5.225), dpi=100)
    figure.set_tight_layout(True)
    ax = figure.add_subplot(111)
    ax.set_title('Cross-sections')
    trapscanvas= FigureCanvasTkAgg(figure, outputframe)
    trapscanvas.get_tk_widget().grid(column = 0,row =0,padx = 10,pady = 10)
    #
    ax.plot([-250,250],[0,0])
    ax.set_xlim([-250,250])
    ax.set_ylim([-1,1])
    ax.set_ylabel("Trap depth (\u03bcK)")
    ax.set_xlabel("Position (\u03bcm)")
    figure.tight_layout()
    
    #subframe under the canvas
    outputsubframe = GUI.MakeFrame(frame=outputframe,row=2,column=0,sticky='ns')
    
    toolbar = NavigationToolbar2Tk(trapscanvas,outputsubframe,pack_toolbar=False)
    toolbar.grid(column = 1,row =0,padx = 0,pady = 0, sticky=tk.EW)
    
    ######################
    #subframe for the plotting settings
    outputsettingsframe = GUI.MakeLabelFrame(frame=outputsubframe,row=1,column=0,width=200,height=220,text='Grid details')
    outputsettingsframe.grid(padx=10)
    
    plotsettingssizes = [22,15,4]
    
    # Set grid properties
    #gridpoints
    gridpointsvariable, entry_gridpoints = GUI.FullEntry(frame=outputsettingsframe,row=0,column=0,sizes=plotsettingssizes,description="Gridpoints = ",unit='')
    entry_gridpoints.insert(tk.END,1000.0)
    #fitpoints
    fitpointsvariable, entry_fitpoints = GUI.FullEntry(frame=outputsettingsframe,row=1,column=0,sizes=plotsettingssizes,description="Harmonic fit points = ",unit='')
    entry_fitpoints.insert(tk.END,10.0)
    #linear fitpoints
    linearfitpointsvariable, entry_linearfitpoints = GUI.FullEntry(frame=outputsettingsframe,row=2,column=0,sizes=plotsettingssizes,description="Linear fit points = ",unit='')
    entry_linearfitpoints.insert(tk.END,100.0)
    #gridsize
    gridsizevariable, entry_gridsize = GUI.FullEntry(frame=outputsettingsframe,row=3,column=0,sizes=plotsettingssizes,description="Computation range = ",unit='\u03bcm')
    entry_gridsize.insert(tk.END,500)

    computebuttonframe = tk.Frame(outputsettingsframe)
    computebuttonframe.grid(column = 0,row =4)
    # toggle to show the linear fit    
    showlinearfiton = tk.IntVar()
    entry_showlinearfit = GUI.MakeCheckButton(frame=computebuttonframe,row=0,column=0,text='Show linear fit',command=lambda: null,state='active',variable=showlinearfiton)
    holdplot = tk.IntVar()
    entry_holdplot = GUI.MakeCheckButton(frame=computebuttonframe,row=0,column=1,text='Hold plot',command=lambda: null,state='active',variable=holdplot)
       
    
    #entry_showlinearfit.grid(row=4,column = 0,sticky=tk.E)
    
    
    ########
    #output textbox 1
    outputtext1frame = GUI.MakeText(outputsubframe,row=0,column=2,width=24,height=15)
    outputtext1frame.grid(rowspan=2)
    outputtext1frame.tag_config('justified', justify='right')
    #output textbox 2
    outputtext2frame = GUI.MakeText(outputsubframe,row=1,column=1,width=45,height=10)
    #outputtext2frame.grid(rowspan=2)
    outputtext2frame.tag_config('justified', justify='right')
    #output textbox 3
    #outputtext3frame = GUI.MakeText(outputsubframe,row=0,column=2,width=24,height=15)
    #outputtext3frame.grid(rowspan=2)
    #outputtext3frame.tag_config('justified', justify='right')
    ###########
    
    outputtext1 = 'calculated parameters: \n'
    outputtext1 += 'fx = '+str(0)+' Hz\n'
    outputtext1 += 'fy = '+str(0)+' Hz\n'
    outputtext1 += 'fz = '+str(0)+' Hz\n'
    outputtext1 += 'f_mean = '+str(0)+' Hz\n'
    outputtext1 += 'vert sag = '+str(0)+' \u03bcm\n'
    outputtext1 += '\n'
    outputtext1 += 'ax = '+str(0)+' m/s^2\n'
    outputtext1 += 'ay = '+str(0)+' m/s^2\n'
    outputtext1 += 'az = '+str(0)+' m/s^2\n'
    outputtext1 += '\n'
    outputtext1 += 'X depth = '+str(0)+' \u03bcK\n'
    outputtext1 += 'Y depth = '+str(0)+' \u03bcK\n'
    outputtext1 += 'Z depth = '+str(0)+' \u03bcK\n'

    outputtext2 = 'Peak photon scattering rates: \n'
    outputtext2 += 'Rayleigh scattering = '+str(0)+' Hz\n'
    #outpt2 += 'Raman scattering = '+str(np.around(0.0,decimals=2))+' Hz\n'
    outputtext2 += 'Heating rate = '+str(0)+' nK/s \n'
    outputtext2 += '\n'
    
    outputtext2 += 'Effective B-field beams = \n '
    outputtext2 += 'Bx  =  '+str(0)+' G \n'
    outputtext2 += 'By  =  '+str(0)+' G \n'
    outputtext2 += 'Bz  =  '+str(0)+' G \n'
    
    outputtext1frame.delete(1.0,tk.END)
    outputtext1frame.insert(tk.END,outputtext1, 'justified')
    outputtext2frame.delete(1.0,tk.END)
    outputtext2frame.insert(tk.END,outputtext2, 'justified')
    
    
    ########################################################################################################################################################################################
    #Optical trap entries 
    opticaltrapssizes = [13,15,5]
    opticaltrapssizes2 = [7,7,5]

    beamparametersubframe = tk.Frame(beams_frame)
    beamparametersubframe.grid(row=1,column=0)
    
    verticalseparator = ttk.Separator(beamparametersubframe, orient='vertical').grid(row=0,column=1,rowspan = 9,sticky='ns')

	
    # Checkbutton for elliptical beams
    ######Beam along the x axis
    Xpowervariable, entry_Xpower = GUI.FullEntry(frame=beamparametersubframe,row=0,column=0,sizes=opticaltrapssizes2,description="X  P= ",unit="mW ")
    entry_Xpower.insert(tk.END,0)
    #
    XwaistYvariable, entry_XwaistY = GUI.FullEntry(frame=beamparametersubframe,row=1,column=0,sizes=opticaltrapssizes2,description="X wy=",unit="\u03bcm ")
    entry_XwaistY.insert(tk.END,100.0)
    #
    XwaistZvariable, entry_XwaistZ = GUI.FullEntry(frame=beamparametersubframe,row=2,column=0,sizes=opticaltrapssizes2,description="X wz=",unit="\u03bcm ")
    entry_XwaistZ.insert(tk.END,100.0)
    #
    XposYvariable, entry_XposY = GUI.FullEntry(frame=beamparametersubframe,row=1,column=2,sizes=opticaltrapssizes2,description="X y0=",unit="\u03bcm ")
    entry_XposY.insert(tk.END,0.0)
    #entry_XposY.config(state= "disabled") 
    #
    XposZvariable, entry_XposZ = GUI.FullEntry(frame=beamparametersubframe,row=2,column=2,sizes=opticaltrapssizes2,description="X z0=",unit="\u03bcm ")
    entry_XposZ.insert(tk.END,0.0)
    #entry_XposZ.config(state= "disabled") 
    ######SecondBeam along the x axis
    Xpowervariable2, entry_Xpower2 = GUI.FullEntry(frame=beamparametersubframe,row=3,column=0,sizes=opticaltrapssizes2,description="X2 P= ",unit="mW ")
    entry_Xpower2.insert(tk.END,0)
    #
    XwaistYvariable2, entry_XwaistY2 = GUI.FullEntry(frame=beamparametersubframe,row=4,column=0,sizes=opticaltrapssizes2,description="X2wy=",unit="\u03bcm ")
    entry_XwaistY2.insert(tk.END,100.0)
    #
    XwaistZvariable2, entry_XwaistZ2 = GUI.FullEntry(frame=beamparametersubframe,row=5,column=0,sizes=opticaltrapssizes2,description="X2wz=",unit="\u03bcm ")
    entry_XwaistZ2.insert(tk.END,100.0)
    #
    XposYvariable2, entry_XposY2 = GUI.FullEntry(frame=beamparametersubframe,row=4,column=2,sizes=opticaltrapssizes2,description="X2y0=",unit="\u03bcm ")
    entry_XposY2.insert(tk.END,0.0)
    #entry_XposY.config(state= "disabled") 
    #
    XposZvariable2, entry_XposZ2 = GUI.FullEntry(frame=beamparametersubframe,row=5,column=2,sizes=opticaltrapssizes2,description="X2z0=",unit="\u03bcm ")
    entry_XposZ2.insert(tk.END,0.0)
    #entry_XposZ.config(state= "disabled") 
    ######Beam along the y axis
    Ypowervariable, entry_Ypower = GUI.FullEntry(frame=beamparametersubframe,row=6,column=0,sizes=opticaltrapssizes2,description="Y  P =",unit="mW ")
    entry_Ypower.insert(tk.END,0)
    #
    YwaistXvariable, entry_YwaistX = GUI.FullEntry(frame=beamparametersubframe,row=7,column=0,sizes=opticaltrapssizes2,description="Y wx=",unit="\u03bcm ")
    entry_YwaistX.insert(tk.END,100.0)
    #
    YwaistZvariable, entry_YwaistZ = GUI.FullEntry(frame=beamparametersubframe,row=8,column=0,sizes=opticaltrapssizes2,description="Y wz=",unit="\u03bcm ")
    entry_YwaistZ.insert(tk.END,100.0)
    #
    YposXvariable, entry_YposX = GUI.FullEntry(frame=beamparametersubframe,row=7,column=2,sizes=opticaltrapssizes2,description="Y x0=",unit="\u03bcm ")
    entry_YposX.insert(tk.END,0.0)
    #entry_YposX.config(state= "disabled") 
    #
    YposZvariable, entry_YposZ = GUI.FullEntry(frame=beamparametersubframe,row=8,column=2,sizes=opticaltrapssizes2,description="Y z0=",unit="\u03bcm ")
    entry_YposZ.insert(tk.END,0.0)
    #entry_YposZ.config(state= "disabled") 
    ######Beam along the z axis
    Zpowervariable, entry_Zpower = GUI.FullEntry(frame=beamparametersubframe,row=9,column=0,sizes=opticaltrapssizes2,description="Z  P =",unit="mW ")
    entry_Zpower.insert(tk.END,0)
    #
    ZwaistYvariable, entry_ZwaistY = GUI.FullEntry(frame=beamparametersubframe,row=10,column=0,sizes=opticaltrapssizes2,description="Z wy=",unit="\u03bcm ")
    entry_ZwaistY.insert(tk.END,100.0)
    #
    ZwaistXvariable, entry_ZwaistX = GUI.FullEntry(frame=beamparametersubframe,row=11,column=0,sizes=opticaltrapssizes2,description="Z wx=",unit="\u03bcm ")
    entry_ZwaistX.insert(tk.END,100.0)
    #
    ZposXvariable, entry_ZposX = GUI.FullEntry(frame=beamparametersubframe,row=10,column=2,sizes=opticaltrapssizes2,description="Z x0=",unit="\u03bcm ")
    entry_ZposX.insert(tk.END,0.0)
    #entry_ZposX.config(state= "disabled") 
    #
    ZposYvariable, entry_ZposY = GUI.FullEntry(frame=beamparametersubframe,row=11,column=2,sizes=opticaltrapssizes2,description="Z y0=",unit="\u03bcm ")
    entry_ZposY.insert(tk.END,0.0)
    #entry_ZposY.config(state= "disabled") 
    
    
    #############setting the polarization
    polarization = ''    
    def set_polarization(pol):
        global polarization
        polarization = pol
        
        if pol ==0:
            Pibutton['relief'] = 'sunken'
            SigmaPlusbutton['relief'] = 'raised'
            SigmaMinusbutton['relief'] = 'raised'
            xpolbutton['relief'] = 'raised'
        if pol ==1:
            Pibutton['relief'] = 'raised'
            SigmaPlusbutton['relief'] = 'sunken'
            SigmaMinusbutton['relief'] = 'raised'
            xpolbutton['relief'] = 'raised'
        if pol ==-1:
            Pibutton['relief'] = 'raised'
            SigmaPlusbutton['relief'] = 'raised'
            SigmaMinusbutton['relief'] = 'sunken'
            xpolbutton['relief'] = 'raised'
        if pol ==[-1,1]:
            Pibutton['relief'] = 'raised'
            SigmaPlusbutton['relief'] = 'raised'
            SigmaMinusbutton['relief'] = 'raised'
            xpolbutton['relief'] = 'sunken'
            
        return
    polarizationbuttonframe= GUI.MakeFrame(frame=beams_frame,row=12,column=0,sticky='ns')
    polarizationbuttonframe.grid(columnspan =3,sticky=tk.N)
    
    #Pibutton = GUI.MakeButton(frame=polarizationbuttonframe,row=0,column=0,width=8,text='(pi) z',command=lambda: set_polarization(0),state='normal',relief='raised')
    #SigmaPlusbutton = GUI.MakeButton(frame=polarizationbuttonframe,row=0,column=1,width=8,text='(+) x+iy',command=lambda: set_polarization(1),state='normal',relief='raised')
    #SigmaMinusbutton = GUI.MakeButton(frame=polarizationbuttonframe,row=0,column=2,width=8,text='(-) x-iy',command=lambda: set_polarization(-1),state='normal',relief='raised')
    #xpolbutton = GUI.MakeButton(frame=polarizationbuttonframe,row=0,column=4,width=8,text='(lin) x or y',command=lambda: set_polarization([-1,1]),state='normal',relief='raised')
    #set_polarization(0)
    
    
    
    #Switching elliptical beams on and off
    def switchwaistentries():
        if entry_ZwaistX["state"] == 'disabled':
            entry_XwaistZ.config(state= "normal") 
            entry_XwaistZ2.config(state= "normal") 
            entry_YwaistZ.config(state= "normal") 
            entry_ZwaistX.config(state= "normal") 
            return
            
        if entry_ZwaistX["state"] == 'normal':
            entry_XwaistZ.config(state= "disabled") 
            entry_XwaistZ2.config(state= "disabled") 
            entry_YwaistZ.config(state= "disabled") 
            entry_ZwaistX.config(state= "disabled") 
            return
    switchwaistentries()
    ellipticalbeams_variable = 0
    ellipticalbutton = GUI.MakeCheckButton(frame=beams_frame,row=0,column=0,text='Elliptical',command=switchwaistentries,state='active',variable=ellipticalbeams_variable)
    # checkbutton.select()
    # checkbutton.deselect()
    
    ########################################################################################################################################################################################
    ###########
    #Magnetic trap entries 
    magnetictrapssizes = [15,15,7]
    #
    magnetictraps_label = GUI.MakeLabel(frame=magnetictraps_frame,row=0,column=0,width=40,text="Circular coil pairs centered on Z-axis",anchor='w')
    magnetictraps_label.grid(columnspan=3)
    #
    Biasvariable, entry_Bias = GUI.FullEntry(frame=magnetictraps_frame,row=1,column=0,sizes=magnetictrapssizes,description="Bias Z = ",unit=" G ")
    entry_Bias.insert(tk.END,0)
    #
    Quadrupolevariable, entry_Quadrupole = GUI.FullEntry(frame=magnetictraps_frame,row=2,column=0,sizes=magnetictrapssizes,description="Grad Z = ",unit=" G/cm ")
    entry_Quadrupole.insert(tk.END,0)
    #
    Curvaturevariable, entry_Curvature = GUI.FullEntry(frame=magnetictraps_frame,row=3,column=0,sizes=magnetictrapssizes,description="curv Z = ",unit=" G/cm\u00b2 ")
    entry_Curvature.insert(tk.END,0)
    #
    # entry_Bias       = name_entry_unit(magnetictrapsframe,1,'Bias Z      = ','   G          ',10)
    # entry_Quadrupole = name_entry_unit(magnetictrapsframe,2,'Grad Z      = ','   G/cm   ',10)
    # entry_curvature  = name_entry_unit(magnetictrapsframe,3,'Curvature Z = ','   G/cm\u00b2  ',10)
    # add_label(magnetictrapsframe, 'Transverse fields', 5, 0,1,1)
    # entry_BiasX       = name_entry_unit(magnetictrapsframe,6,'Bias X      = ','   G          ',10)
    # entry_BiasY       = name_entry_unit(magnetictrapsframe,7,'Bias Y      = ','   G          ',10)
    
    # entry_Quadrupole.insert(tk.END, 0.0)
    # entry_Bias.insert(tk.END, 0.0)
    # entry_curvature.insert(tk.END, 0.0)
    # entry_BiasX.insert(tk.END, 0.0)
    # entry_BiasY.insert(tk.END, 0.0)
    # entry_BiasX.config(state= "disabled") 
    # entry_BiasY.config(state= "disabled") 
    
    ########################################################################################################################################################################################
    ##################
    #Atomic details
    atombuttonwidth= 8
    atompropssizes = [22,12,9]
    
    buttonframe= GUI.MakeFrame(frame=atoms_frame,row=0,column=0,sticky='ns')
   
    statsubframe = tk.Frame(atoms_frame)
    statsubframe.grid(row=10,column=0)
    label_state = GUI.MakeLabel(frame=statsubframe,row=10,column=0,width=25,text="state = ",anchor='e')
    state, statedropdown= GUI.Make_DropDownMenu(statsubframe,row=10,column=1,columnspan=1,rowspan=1,choices=[[0,0]])
    statedropdown.grid(sticky='nw')
    
    species = ''    
    def set_atom(atom,button):
        global species
        nonlocal statedropdown
        nonlocal state
        species = atom
        
        atomdata = importlib.import_module("AtomicData."+atom) 
        entry_mass.delete(0, tk.END)  
        entry_mass.insert(tk.END, GUI.DecimalNotation(atomdata.m/scipy.constants.physical_constants["atomic mass constant"][0],2))
       
        #Enter the hyperfine states into the dropdown menu
        tmpoperator = Operator(atomdata.I,0.5)
        statedropdown['menu'].delete(0, 'end')
        
        lbl = str(tmpoperator.FmFindex[0][0])+','+str(tmpoperator.FmFindex[0][1])
        
        state.set(lbl)
        for choice in tmpoperator.FmFindex:
            lbl = str(choice[0])+','+str(choice[1])
            statedropdown['menu'].add_command(label=lbl, command=tk._setit(state, lbl))
        
        for widget in atomchoicesframe.winfo_children():
            widget.config(relief='raised')                   
        button.config(relief='sunken')    
        return
    
    
    atomchoicesframe = GUI.MakeLabelFrame(frame=buttonframe,row=1,column=0,width=50,height=50,text='')
    atomchoicesframe.grid(columnspan=3)
    #GUI.MakeFrame(frame=buttonframe,row=1,column=0,sticky='ns')
    
    # Li6Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='\u2076Li',command=lambda: set_atom('6Li'),state='normal',relief='raised')
    # Na23Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='\u00b2\u00b3Na',command=lambda: set_atom('23Na'),state='normal',relief='raised')
    # K39Button = GUI.MakeButton(frame=atomchoicesframe,row=2,column=0,width=atombuttonwidth,text='\u00b3\u2079K',command=lambda: set_atom('39K'),state='normal',relief='raised')
    # K40Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='\u2074\u2070K',command=lambda: set_atom('40K'),state='normal',relief='raised')
    # Rb85Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='\u2078\u2075Rb',command=lambda: set_atom('85Rb'),state='normal',relief='raised')
    # Rb87Button = GUI.MakeButton(frame=atomchoicesframe,row=2,column=1,width=atombuttonwidth,text='\u2078\u2077Rb',command=lambda: set_atom('87Rb'),state='normal',relief='raised')
    # CsButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text='\u00B9\u00B3\u00B3Cs',command=lambda: set_atom('133Cs'),state='normal',relief='raised')   
    
    # Ca40Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='40Ca',command=lambda: set_atom('40Ca'),state='normal',relief='raised')
    # SrBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Bos Sr',command=lambda: set_atom('B-Sr'),state='normal',relief='raised')
    # Sr87Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='87Sr',command=lambda: set_atom('87Sr'),state='normal',relief='raised')
    
    # EuBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='Bos Eu',command=lambda: set_atom('B-Eu'),state='normal',relief='raised')
    # DyBosButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='Bos Dy',command=lambda: set_atom('B-Dy'),state='normal',relief='raised')
    # DyFerButton = GUI.MakeButton(frame=atomchoicesframe,row=2,column=0,width=atombuttonwidth,text='Fer Dy',command=lambda: set_atom('F-Dy'),state='normal',relief='raised')
    # ErBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Bos Er',command=lambda: set_atom('B-Er'),state='normal',relief='raised')
    # Er167Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='F-167Er',command=lambda: set_atom('F-167Er'),state='normal',relief='raised')
    # YtBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text='Bos Yt',command=lambda: set_atom('B-Yt'),state='normal',relief='raised')
    # Yt171Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=2,width=atombuttonwidth,text='171Yt',command=lambda: set_atom('171Yt'),state='normal',relief='raised')
    # Yt173Button = GUI.MakeButton(frame=atomchoicesframe,row=2,column=2,width=atombuttonwidth,text='173Yt',command=lambda: set_atom('173Yt'),state='normal',relief='raised')
               
    
    def set_atom_choices(atomset):
        # nonlocal Li6Button
        # nonlocal Na23Button
        # nonlocal K39Button
        # nonlocal K40Button
        # nonlocal Rb85Button
        # nonlocal Rb87Button
        # nonlocal CsButton
        
        # nonlocal Ca40Button
        # nonlocal SrBosButton
        # nonlocal Sr87Button
        
        # nonlocal EuBosButton
        # nonlocal DyBosButton
        # nonlocal DyFerButton
        # nonlocal ErBosButton
        # nonlocal Er167Button
        # nonlocal YtBosButton
        # nonlocal Yt171Button
        # nonlocal Yt173Button
        
        # Li6Button.destroy()
        # Na23Button.destroy()
        # K39Button.destroy()
        # K40Button.destroy()
        # Rb85Button.destroy()
        # Rb87Button.destroy()
        # CsButton.destroy()
        
        # Ca40Button.destroy()
        # SrBosButton.destroy()
        # Sr87Button.destroy()
        
        # EuBosButton.destroy()
        # DyBosButton.destroy()
        # DyFerButton.destroy()
        # ErBosButton.destroy()
        # Er167Button.destroy()
        # YtBosButton.destroy()
        # Yt171Button.destroy()
        # Yt173Button.destroy()
               
        
        for widget in atomchoicesframe.winfo_children():
            widget.destroy()
        for i in np.arange(9):
            fillerbutton = GUI.MakeButton(frame=atomchoicesframe,row=floor(i/3),column=i%3,width=atombuttonwidth,text=' ',command=lambda: null,state='normal',relief='raised')       
               
        if atomset == 'alkali':
            Li6Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='\u2076Li',command=lambda: set_atom('6Li',Li6Button),state='normal',relief='raised')
            Na23Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='\u00b2\u00b3Na',command=lambda: set_atom('23Na',Na23Button),state='normal',relief='raised')
            K39Button = GUI.MakeButton(frame=atomchoicesframe,row=2,column=0,width=atombuttonwidth,text='\u00b3\u2079K',command=lambda: set_atom('39K',K39Button),state='normal',relief='raised')
            K40Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='\u2074\u2070K',command=lambda: set_atom('40K',K40Button),state='normal',relief='raised')
            Rb85Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='\u2078\u2075Rb',command=lambda: set_atom('85Rb',Rb85Button),state='normal',relief='raised')
            Rb87Button = GUI.MakeButton(frame=atomchoicesframe,row=2,column=1,width=atombuttonwidth,text='\u2078\u2077Rb',command=lambda: set_atom('87Rb',Rb87Button),state='normal',relief='raised')
            CsButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text='\u00B9\u00B3\u00B3Cs',command=lambda: set_atom('133Cs',CsButton),state='normal',relief='raised')
            
        if atomset == 'Alkaline earth':
            Ca40Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='40Ca',command=lambda: set_atom('40Ca',Ca40Button),state='normal',relief='raised')
            SrBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Bos Sr',command=lambda: set_atom('B-Sr',SrBosButton),state='normal',relief='raised')
            Sr87Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='87Sr',command=lambda: set_atom('87Sr',Sr87Button),state='normal',relief='raised')
            
        if atomset == 'Lanthanide':
            EuBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='Bos Eu',command=lambda: set_atom('B-Eu',EuBosButton),state='normal',relief='raised')
            DyBosButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='Bos Dy',command=lambda: set_atom('B-Dy',DyBosButton),state='normal',relief='raised')
            DyFerButton = GUI.MakeButton(frame=atomchoicesframe,row=2,column=0,width=atombuttonwidth,text='Fer Dy',command=lambda: set_atom('F-Dy',DyFerButton),state='normal',relief='raised')
            ErBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Bos Er',command=lambda: set_atom('B-Er',ErBosButton),state='normal',relief='raised')
            Er167Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='F-167Er',command=lambda: set_atom('F-167Er',Er167Button),state='normal',relief='raised')
            YtBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text='Bos Yt',command=lambda: set_atom('B-Yt',YtBosButton),state='normal',relief='raised')
            Yt171Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=2,width=atombuttonwidth,text='171Yt',command=lambda: set_atom('171Yt',Yt171Button),state='normal',relief='raised')
            Yt173Button = GUI.MakeButton(frame=atomchoicesframe,row=2,column=2,width=atombuttonwidth,text='173Yt',command=lambda: set_atom('173Yt',Yt173Button),state='normal',relief='raised')
        
        return
    
    alkalibuttonButton = GUI.MakeButton(frame=buttonframe,row=0,column=0,width=10,text='alkali',command=lambda: set_atom_choices('alkali'),state='normal',relief='raised')
    earthalkalibuttonButton = GUI.MakeButton(frame=buttonframe,row=0,column=1,width=10,text='Alkaline earth',command=lambda: set_atom_choices('Alkaline earth'),state='normal',relief='raised')
    LanthanidesbuttonButton = GUI.MakeButton(frame=buttonframe,row=0,column=2,width=10,text='Lanthanide',command=lambda: set_atom_choices('Lanthanide'),state='normal',relief='raised')
    set_atom_choices('alkali')
    
    #Atomic mass
    massvariable, entry_mass = GUI.FullEntry(frame=atoms_frame,row=1,column=0,sizes=atompropssizes,description="Atomic mass = ",unit="  u ")
    entry_mass.insert(tk.END,0)
    
    
    gravon = tk.IntVar()
    entry_gravity = GUI.MakeCheckButton(frame=atoms_frame,row=2,column=0,text='Gravity',command=lambda : None ,state='active',variable=gravon).grid(columnspan=3)
    
    #Trapping wavelength
    
    wavelengthvariable, entry_wavelength = GUI.FullEntry(frame=atoms_frame,row=3,column=0,sizes=opticaltrapssizes,description="Wavelength = ",unit="nm ")
    entry_wavelength.insert(tk.END,1064.0)
    
    #Magnetic moment
    magmomentvariable, entry_magmoment = GUI.FullEntry(frame=atoms_frame,row=13,column=0,sizes=atompropssizes,description="Magnetic moment = ",unit=" muB ")
    entry_magmoment.insert(tk.END,0)
    #Scalar polarizability
    polvariable, entry_pol = GUI.FullEntry(frame=atoms_frame,row=5,column=0,sizes=atompropssizes,description="Scalar polarizability (re)=",unit=' e\u00b2 a\u2080\u00b2 / E\u2095')
    entry_pol.insert(tk.END,0)
    #Vector polarizability
    pol_vecvariable, entry_pol_vec = GUI.FullEntry(frame=atoms_frame,row=6,column=0,sizes=atompropssizes,description="Vector polarizability (re)=",unit=' e\u00b2 a\u2080\u00b2 / E\u2095')
    entry_pol_vec.insert(tk.END,0)
    
    
    
    
    ########################################################################################################################################################################################
    
    def insert_polarizability():
        global species
        atom = importlib.import_module("AtomicData."+species)
              
        wavelength = float(entry_wavelength.get())
        omega = 2.0*np.pi * sc.c/(wavelength)
        alpha = compute_polarizability(species,wavelength,0)
        alpha_s,alpha_v,alpha_t = alpha[0],alpha[1],alpha[2]
        

        entry_pol.delete(0, tk.END)  
        entry_pol.insert(tk.END,np.format_float_scientific( np.real(alpha_s),precision = 4))
        entry_pol_vec.delete(0, tk.END)  
        entry_pol_vec.insert(tk.END, np.format_float_scientific(np.real(alpha_v),precision = 4))
        return 
        
        
    separator = ttk.Separator(atoms_frame, orient='horizontal').grid(row=4,column = 0, sticky="ew" )
    computepolbutton = GUI.MakeButton(frame=atoms_frame,row=4,column=0,width=20,text='Compute polarizability',command=lambda: insert_polarizability(),state='normal',relief='raised')
    computepolbutton.grid(sticky=tk.N)
    ################################################################################################################################################################################
        
    def compute_magneticmoment():
        global species
        atom = importlib.import_module("AtomicData."+species)
        
        chosenstate = state.get()
                
        Bias = float(entry_Bias.get())*1.0e-4
        
        f, mf = float(chosenstate[0:3]), float(chosenstate[4:])
                
        if f == atom.I+atom.stateJ[0]:
            s = +1
        else:
            s = -1
        
        B1 = BreitRabi(Bias,mf,s,species)
        B2 = BreitRabi(Bias+1.0e-6,mf,s,species)
        
        moment = (B1-B2)/(1.0e-6)
        
        
        
        entry_magmoment.delete(0, tk.END)  
        entry_magmoment.insert(tk.END,np.around(moment/muB,decimals=4))
        return 
        
    
    
    separator = ttk.Separator(atoms_frame, orient='horizontal').grid(row=12,column = 0,columnspan=3, sticky="ew" )
    computemagneticmomentbutton = GUI.MakeButton(frame=atoms_frame,row=12,column=0,width=25,text='Compute moment at bias field',command=lambda: compute_magneticmoment(),state='normal',relief='raised')
    computemagneticmomentbutton.grid(sticky=tk.N)
    
    
    ################################################################################################################################################################################
        
    def compute_traps():    
        nonlocal ax 
        nonlocal trapscanvas
        #%% Reading entries
        gridpoints = float(entry_gridpoints.get())   
        fitpoints  = float(entry_fitpoints.get())   
        linearfitpoints  = float(entry_linearfitpoints.get())   
        gridsize   = float(entry_gridsize.get())*1.0e-6
        wavelength = float(entry_wavelength.get())*1.0e-9
        
        # entry_XwaistY.insert(tk.END, 100.0)
        # entry_XwaistZ.insert(tk.END, 100.0)
        # entry_YwaistX.insert(tk.END, 100.0)
        # entry_YwaistZ.insert(tk.END, 100.0)
        # entry_ZwaistX.insert(tk.END, 100.0)
        # entry_ZwaistY.insert(tk.END, 100.0)
        
        XwaistY     = float(entry_XwaistY.get())*1.0e-6
        XwaistY2     = float(entry_XwaistY2.get())*1.0e-6
        YwaistX     = float(entry_YwaistX.get())*1.0e-6
        ZwaistX     = float(entry_ZwaistX.get())*1.0e-6
        
        if entry_ZwaistX["state"] == 'normal':
            XwaistZ     = float(entry_XwaistZ.get())*1.0e-6
            XwaistZ2    = float(entry_XwaistZ2.get())*1.0e-6
            YwaistZ     = float(entry_YwaistZ.get())*1.0e-6
            ZwaistY     = float(entry_ZwaistY.get())*1.0e-6 
        else:
            XwaistZ = 1.0*XwaistY
            XwaistZ2 = 1.0*XwaistY2
            YwaistZ = 1.0*YwaistX
            ZwaistY = 1.0*ZwaistX
            
        
        powerX     = float(entry_Xpower.get())*1.0e-3 
        powerX2    = float(entry_Xpower2.get())*1.0e-3 
        powerY     = float(entry_Ypower.get())*1.0e-3 
        powerZ     = float(entry_Zpower.get())*1.0e-3 
        polarizability = float(entry_pol.get())*polarizabilityunit
        #polarizabilityimag = float(scat_entry_pol_im.get())*polarizabilityunit
        vectorpolarizability = float(entry_pol_vec.get())*polarizabilityunit
        #vectorpolarizabilityimag = float(scat_entry_pol_vec_im.get())*polarizabilityunit
        Quadrupole = float(entry_Quadrupole.get())*0.01
        curvature = float(entry_Curvature.get())
        Bias = float(entry_Bias.get())*1.0e-4
        BiasX = 0.0 #float(entry_BiasX.get())*1.0e-4
        BiasY = 0.0 #float(entry_BiasY.get())*1.0e-4
        m = float(entry_mass.get())*scipy.constants.physical_constants["atomic mass constant"][0]
        magmoment = float(entry_magmoment.get())*sc.value(u'Bohr magneton')

        Erec = sc.hbar*sc.hbar*(2.0*np.pi/wavelength)*(2.0*np.pi/wavelength)/(2.0*m)
        Trec = 2.0*Erec/sc.k
        
        # Make the number of points even for convenience later on.
        if gridpoints % 2 == 1:
            gridpoints += 1
        if fitpoints % 2 == 1:
            fitpoints += 1
        stepsize = gridsize/gridpoints
        
        #Constructing coordinates
        x = np.arange(-int(gridpoints/2),int(gridpoints/2))*stepsize
        y = np.arange(-int(gridpoints/2),int(gridpoints/2))*stepsize
        z = np.arange(-int(gridpoints/2),int(gridpoints/2))*stepsize
        
        #%%Building optical potential
        kdipoletrap   = 2.0*np.pi/wavelength
        omegaTrap     = kdipoletrap*sc.c
        
        zRX_Y = 0.5*kdipoletrap*XwaistY*XwaistY
        zRX_Z = 0.5*kdipoletrap*XwaistZ*XwaistZ
        
        zRX2_Y = 0.5*kdipoletrap*XwaistY2*XwaistY2
        zRX2_Z = 0.5*kdipoletrap*XwaistZ2*XwaistZ2
        
        zRY_X = 0.5*kdipoletrap*YwaistX*YwaistX
        zRY_Z = 0.5*kdipoletrap*YwaistZ*YwaistZ
        
        zRZ_X = 0.5*kdipoletrap*ZwaistX*ZwaistX
        zRZ_Y = 0.5*kdipoletrap*ZwaistY*ZwaistY
        
        IX = 2.0 * powerX/(np.pi * XwaistY*XwaistZ)
        IX2 = 2.0 * powerX2/(np.pi * XwaistY2*XwaistZ2)
        IY = 2.0 * powerY/(np.pi * YwaistX*YwaistZ)
        IZ = 2.0 * powerZ/(np.pi * ZwaistX*ZwaistY)
        
        EsqrdX   = 2.0*IX/(sc.epsilon_0*sc.c)
        EsqrdX2  = 2.0*IX2/(sc.epsilon_0*sc.c)
        EsqrdY   = 2.0*IY/(sc.epsilon_0*sc.c)
        EsqrdZ   = 2.0*IZ/(sc.epsilon_0*sc.c)
        
        trapdepthX  = 0.25*EsqrdX*polarizability
        trapdepthX2 = 0.25*EsqrdX2*polarizability
        trapdepthY  = 0.25*EsqrdY*polarizability
        trapdepthZ  = 0.25*EsqrdZ*polarizability

        #decayrateX = 0.25*EsqrdX*polarizabilityimag
        #decayrateY = 0.25*EsqrdY*polarizabilityimag
        #decayrateZ = 0.25*EsqrdZ*polarizabilityimag
        
        BeffX = vectorpolarizability*EsqrdX/(8.0*muB*0.5*2.0)
        BeffY = vectorpolarizability*EsqrdY/(8.0*muB*0.5*2.0)
        BeffZ = vectorpolarizability*EsqrdZ/(8.0*muB*0.5*2.0)

        Xy0 = float(entry_XposY.get())*1.0e-6
        Xz0 = float(entry_XposZ.get())*1.0e-6
        Xy02 = float(entry_XposY2.get())*1.0e-6
        Xz02 = float(entry_XposZ2.get())*1.0e-6
        Yx0 = float(entry_YposX.get())*1.0e-6
        Yz0 = float(entry_YposZ.get())*1.0e-6
        Zx0 = float(entry_ZposX.get())*1.0e-6
        Zy0 = float(entry_ZposY.get())*1.0e-6
        
        def fullpotential(xx,yy,zz):
            ## Xbeam
            XwaistfactorY = np.sqrt(1+ (xx/zRX_Y)**2)
            XwaistfactorZ = np.sqrt(1+ (xx/zRX_Z)**2)
            
            XwaistZZEntries = (zz-Xz0)*(zz-Xz0)/(XwaistZ*XwaistZ*XwaistfactorZ*XwaistfactorZ)
            XwaistYYEntries = (yy-Xy0)*(yy-Xy0)/(XwaistY*XwaistY*XwaistfactorY*XwaistfactorY)
            
            pot  = -trapdepthX*np.exp(- 2.0*(XwaistZZEntries+XwaistYYEntries))/(XwaistfactorY*XwaistfactorZ)
                        
            ## Xbeam2
            XwaistfactorY2 = np.sqrt(1+ (xx/zRX2_Y)**2)
            XwaistfactorZ2 = np.sqrt(1+ (xx/zRX2_Z)**2)
            
            XwaistZZEntries2 = (zz-Xz02)*(zz-Xz02)/(XwaistZ2*XwaistZ2*XwaistfactorZ2*XwaistfactorZ2)
            XwaistYYEntries2 = (yy-Xy02)*(yy-Xy02)/(XwaistY2*XwaistY2*XwaistfactorY2*XwaistfactorY2)
            
            pot -= trapdepthX2*np.exp(- 2.0*(XwaistZZEntries2+XwaistYYEntries2))/(XwaistfactorY2*XwaistfactorZ2)
                        
            ## Ybeam
            YwaistfactorX = np.sqrt(1+ (yy/zRY_X)**2)
            YwaistfactorZ = np.sqrt(1+ (yy/zRY_Z)**2)
            
            YwaistXXEntries = (xx-Yx0)*(xx-Yx0)/(YwaistX*YwaistX*YwaistfactorX*YwaistfactorX)
            YwaistZZEntries = (zz-Yz0)*(zz-Yz0)/(YwaistZ*YwaistZ*YwaistfactorZ*YwaistfactorZ)
            
            pot -= trapdepthY*np.exp(- 2.0*(YwaistXXEntries+YwaistZZEntries))/(YwaistfactorX*YwaistfactorZ)
                       
            
            ## Zbeam
            
            pot -= trapdepthZ*np.exp(- 2.0* ((yy-Zy0)*(yy-Zy0) + (xx-Zx0)*(xx-Zx0))/(ZwaistX*ZwaistX *(1+ (zz/zRZ_X)**2) ))/(1+ (zz/zRZ_X)**2)
            
            
            #%%Building magnetic potential
            # Reference, Electromagnet Design Basics for Cold Atom Experiments by Todd Meyrath
            Bz = Bias + zz*Quadrupole + curvature*(zz*zz - 0.5* (xx*xx+yy*yy) )
            Br = -0.5*Quadrupole*np.sqrt(xx*xx+yy*yy) - curvature*zz*np.sqrt(xx*xx+yy*yy)
            Bmag = np.sqrt(Bz*Bz+Br*Br)
            pot -= magmoment*Bmag
            if gravon.get() == True:
                pot += m*9.81*zz 
            return pot
        
        #%%start with the potential along the z-axis, and find the sag due to gradients.
        UZaxis = fullpotential(0.0,0.0,z)
        
        #Check the slope to set the direction of the walk
        if (UZaxis[int(gridpoints/2)+1] > UZaxis[int(gridpoints/2)]):
            stepsign = -1
        else:
            stepsign =  1
        
        trapbroken = False
        
        #walk towards the minimum
        keepwalking = True
        walkingindex = int(gridpoints/2)
        while keepwalking == True:
            if (walkingindex == gridpoints-1 or walkingindex == 1):
                trapbroken = True
                break 
            if (UZaxis[walkingindex+stepsign] > UZaxis[walkingindex]):
                trapcentre = walkingindex
                keepwalking = False
            walkingindex += stepsign
            
        if (trapbroken == False):
            # continue the walk until the 
            keepwalking = True
            walkingindex += 5*stepsign
            while keepwalking == True:
                if (UZaxis[walkingindex+stepsign] < UZaxis[walkingindex] or walkingindex+stepsign < 1 or walkingindex+stepsign> gridpoints -2):
                    zcrest = walkingindex
                    keepwalking = False
                walkingindex += stepsign
        else: 
            trapcentre = int(gridpoints/2)
            zcrest = int(gridpoints/2)
                    
        #Potential energy at the bottom of the trap.
        Ucenter = UZaxis[trapcentre]
        #trap depth along the z-axis
        ztrapdepth = UZaxis[zcrest]-Ucenter
        
        #Fit the curvature to determine trap frequencies and a more precise minimum.
        zfitcoords = z[trapcentre-int(fitpoints/2):trapcentre+int(fitpoints/2)]
        zfitarray  = UZaxis[trapcentre-int(fitpoints/2):trapcentre+int(fitpoints/2)]
        
        linzfitcoords = z[trapcentre-int(linearfitpoints/2):trapcentre+int(linearfitpoints/2)]
        linzfitarray  = UZaxis[trapcentre-int(linearfitpoints/2):trapcentre+int(linearfitpoints/2)]
                
        def harmonicpotential(x,a,b,u0):
            return u0 + a*x*x + b*x
        
        def linearpotential(x,a,u0):
            return u0 + a*np.absolute(x)
                
        zfitopt, zfitcov = optimize.curve_fit(harmonicpotential,zfitcoords,zfitarray, p0=[0,0.0,Ucenter])
        
        
        if (trapbroken == False):
            Z0 = -0.5*zfitopt[1]/zfitopt[0]
        else:
            Z0 = z[int(gridpoints/2)]
        
        #%%Now the X and Y axes, crossing through the trap centre at z0
        
        UYaxis = fullpotential(0.0,y,Z0)
        UXaxis = fullpotential(x,0.0,Z0)
                
        startrangex = int((gridpoints-fitpoints)/2)
        endrangex   = int((gridpoints+fitpoints)/2)
        startrangez = trapcentre - int(fitpoints/2)  
        endrangez   = trapcentre + int(fitpoints/2) 
        
        if endrangez > gridpoints:
            endrangez    -= int(fitpoints)
            startrangez  -= int(fitpoints)
        if startrangez < 0:
            endrangez    += int(fitpoints)
            startrangez  += int(fitpoints)
            
        xfitcoords = x[startrangex:endrangex]
        xfitarray  = UXaxis[startrangex:endrangex]
        
        yfitcoords = y[startrangex:endrangex]
        yfitarray  = UYaxis[startrangex:endrangex]
                
        #%%%linear fit
        
        startrangex = int((gridpoints-linearfitpoints)/2)
        endrangex   = int((gridpoints+linearfitpoints)/2)
        startrangez = trapcentre - int(linearfitpoints/2)  
        endrangez   = trapcentre + int(linearfitpoints/2) 
        
        if endrangez > gridpoints:
            endrangez    -= int(linearfitpoints)
            startrangez  -= int(linearfitpoints)
        if startrangez < 0:
            endrangez    += int(linearfitpoints)
            startrangez  += int(linearfitpoints)
        
        linxfitcoords = x[startrangex:endrangex]
        linxfitarray  = UXaxis[startrangex:endrangex]
        
        linyfitcoords = y[startrangex:endrangex]
        linyfitarray  = UYaxis[startrangex:endrangex]
               
        xfitopt, xfitcov = optimize.curve_fit(harmonicpotential,xfitcoords,xfitarray, p0=[0,0,0.0])
        yfitopt, yfitcov = optimize.curve_fit(harmonicpotential,yfitcoords,yfitarray, p0=[0,0,0.0])
        
        #%% 
        
        linzfitopt, linzfitcov = optimize.curve_fit(linearpotential,linzfitcoords,linzfitarray, p0=[0.0,Ucenter])
        linxfitopt, linxfitcov = optimize.curve_fit(linearpotential,linxfitcoords,linxfitarray, p0=[0.0,Ucenter])
        linyfitopt, linyfitcov = optimize.curve_fit(linearpotential,linyfitcoords,linyfitarray, p0=[0.0,Ucenter])
        
        
        #%%
        keepwalking = True
        walkingindex = int(gridpoints/2)+5
        while keepwalking == True:
            if (UYaxis[walkingindex+1] < UYaxis[walkingindex] or walkingindex+3> gridpoints):
                ycrest = walkingindex
                keepwalking = False
            walkingindex += 1
            
        keepwalking = True
        walkingindex = int(gridpoints/2)+5
        while keepwalking == True:
            if (UXaxis[walkingindex+1] < UXaxis[walkingindex] or walkingindex+3> gridpoints):
                xcrest = walkingindex
                keepwalking = False
            walkingindex += 1
        
        if (trapbroken == False):
            xtrapdepth = UXaxis[xcrest]-Ucenter
            ytrapdepth = UYaxis[ycrest]-Ucenter
        else:
            xtrapdepth = 0.0
            ytrapdepth = 0.0
            ztrapdepth = 0.0
        
        
        ############################
        #%%making the plot
        if holdplot.get() == False:
            ax.cla()
        
        
        ax.plot(x*1.0e6,(UXaxis-Ucenter)*1.0e6/sc.k,label='x-axis',color = 'r') 
        ax.plot(x*1.0e6,np.ones(int(gridpoints))*xtrapdepth*1.0e6/sc.k,linestyle = 'dashed',color = 'r')
        ax.plot(xfitcoords*1.0e6,(harmonicpotential(xfitcoords,*xfitopt)-Ucenter)*1.0e6/sc.k,'o',color = 'r')
        
        ax.plot(y*1.0e6,(UYaxis-Ucenter)*1.0e6/sc.k,label='y-axis',color = 'b')
        ax.plot(y*1.0e6,np.ones(int(gridpoints))*ytrapdepth*1.0e6/sc.k,linestyle = 'dashed',color = 'b')
        ax.plot(yfitcoords*1.0e6,(harmonicpotential(yfitcoords,*yfitopt)-Ucenter)*1.0e6/sc.k,'o',color = 'b')
        
        ax.plot(z*1.0e6,(UZaxis-Ucenter)*1.0e6/sc.k,label='z-axis',color = 'g')
        ax.plot(z*1.0e6,np.ones(int(gridpoints))*ztrapdepth*1.0e6/sc.k,linestyle = 'dashed',color = 'g')
        ax.plot(zfitcoords*1.0e6,(harmonicpotential(zfitcoords,*zfitopt)-Ucenter)*1.0e6/sc.k,'o',color = 'g')
        ax.grid(True)
        
        if showlinearfiton.get() == True:
            ax.plot(linxfitcoords*1.0e6,(linearpotential(linxfitcoords,*linxfitopt)-Ucenter)*1.0e6/sc.k,'^',color = 'maroon')
            ax.plot(linyfitcoords*1.0e6,(linearpotential(linyfitcoords,*linyfitopt)-Ucenter)*1.0e6/sc.k,'^',color = 'navy')
            ax.plot(linzfitcoords*1.0e6,(linearpotential(linzfitcoords,*linzfitopt)-Ucenter)*1.0e6/sc.k,'^',color = 'darkgreen')
            
        if (trapbroken == False):
            maxdepth= np.amax([ztrapdepth,xtrapdepth,ytrapdepth])*1.0e6/sc.k
            ax.set_ylim([-0.1*maxdepth,1.1*maxdepth])
        
        
        ax.set_ylabel("Trap depth (\u03bcK)")
        ax.set_xlabel("Position (\u03bcm)")
        
        ax.legend()
        trapscanvas.draw()
        ############################      
               
               
        #%% printing parameters
        if (trapbroken == False):
            fx = np.sqrt(xfitopt[0]/(2.0*np.pi*np.pi*m))
            fy = np.sqrt(yfitopt[0]/(2.0*np.pi*np.pi*m))
            fz = np.sqrt(zfitopt[0]/(2.0*np.pi*np.pi*m))
        else: 
            fx = 1
            fy = 1
            fz = 1
        
        set_get_frequencies('set',fx,fy,fz,linxfitopt[0]/m,linyfitopt[0]/m,linzfitopt[0]/m)
        
        if (trapbroken == False):
            outpt = 'calculated parameters: \n'
            outpt += 'fx = '+str(np.around(fx,decimals=2))+' Hz\n'
            outpt += 'fy = '+str(np.around(fy,decimals=2))+' Hz\n'
            outpt += 'fz = '+str(np.around(fz,decimals=2))+' Hz\n'
            outpt += 'f_mean = '+str(np.around(np.power(fz*fy*fx,1.0/3.0),decimals=2))+' Hz\n'
            outpt += 'vert sag = '+str(np.around(Z0*1.0e6,decimals=2))+' \u03bcm\n'
            outpt += '\n'
            outpt += 'ax = '+str(np.around(linxfitopt[0]/m,decimals=2))+' m/s^2\n'
            outpt += 'ay = '+str(np.around(linyfitopt[0]/m,decimals=2))+' m/s^2\n'
            outpt += 'az = '+str(np.around(linzfitopt[0]/m,decimals=2))+' m/s^2\n'
            outpt += '\n'
            outpt += 'X depth = '+str(np.around(xtrapdepth*1.0e6/sc.k,decimals=2))+' \u03bcK\n'
            outpt += 'Y depth = '+str(np.around(ytrapdepth*1.0e6/sc.k,decimals=2))+' \u03bcK\n'
            outpt += 'Z depth = '+str(np.around(ztrapdepth*1.0e6/sc.k,decimals=2))+' \u03bcK\n'

           # print('')
           # print(trapdepthX/sc.h)
           # print(trapdepthX/sc.k)
            #print(xtrapdepth/sc.h)
            #print(xtrapdepth/sc.k)
           # print(decayrateX/sc.h)
            outpt2 = 'Peak photon scattering rates: \n'
            outpt2 += 'Rayleigh scattering = '+str(np.around(0.0/sc.hbar,decimals=2))+' Hz\n'
            #outpt2 += 'Raman scattering = '+str(np.around(0.0,decimals=2))+' Hz\n'
            outpt2 += 'Heating rate = '+str(np.around((0.0/sc.h)*Trec*1.0e9,decimals=2))+' nK/s \n'
            outpt2 += '\n'
            
            outpt2 += 'Effective B-field beams for R+ pol = \n '
            outpt2 += 'Bx = '+str(np.around(BeffX*1.0e7,decimals=4))+' mG \n'
            outpt2 += 'By = '+str(np.around(BeffY*1.0e7,decimals=4))+' mG \n'
            outpt2 += 'Bz = '+str(np.around(BeffZ*1.0e7,decimals=4))+' mG \n'
			
            
        else:
            outpt = 'No trap !\n'
            outpt2 = ' '
        
        outputtext1frame.delete(1.0,tk.END)
        outputtext1frame.insert(tk.END,outpt, 'justified')
        
        outputtext2frame.delete(1.0,tk.END)
        outputtext2frame.insert(tk.END,outpt2, 'justified')
        return
    
    
    
    computebutton = GUI.MakeButton(frame=computebuttonframe,row=1,column=1,width=12,text='Plot potential',command=lambda: compute_traps(),state='normal',relief='raised')
    
    #GUI.MakeLabel(frame=atoms_frame,row=5,column=0,width=labelwidth_atom,text=" ",anchor='e')
    
    
