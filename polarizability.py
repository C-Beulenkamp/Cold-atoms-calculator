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
import math
from math import floor
import scipy
import scipy.constants as sc
from AngularMomentumOperatorClass import Operator,AngMomOperators
from sympy.physics.wigner import wigner_3j,wigner_6j
from sympy import N
from scipy import optimize
from extrafunctions import *
import glob, os




sharedfx = 0.0
sharedfy = 0.0
sharedfz = 0.0
sharedax = 0.0
shareday = 0.0
sharedaz = 0.0
muB  = sc.value(u'Bohr magneton') # Bohr magneton
polarizabilityunit = sc.e*sc.e* sc.physical_constants['Bohr radius'][0]*sc.physical_constants['Bohr radius'][0]/(sc.physical_constants['Hartree energy'][0])


def compute_polarizability(species,wavelength,stateindex):
    atom = importlib.import_module("AtomicData."+species)
    omega = 2.0*np.pi * sc.c/(wavelength*1.0e-9)
    
    #Scalar polarizability:        
    
    alpha0 = 0.0j
    alpha1 = 0.0j
    alpha2 = 0.0j

    
    for nj in np.arange(atom.statenumber):
        print(nj)
        print(atom.states[nj])
        print(QuietWigner6j(1.0,0.0,1.0,atom.stateJ[stateindex],atom.stateJ[nj],atom.stateJ[stateindex]))

        transitionEnergy = atom.state_energies[nj]-atom.state_energies[stateindex]
        photonEnergy = sc.hbar*omega
        prefactor = np.power(-1.0,atom.stateJ[stateindex]+atom.stateJ[nj]+1.0)
        dipolesqrd = atom.dipolemoments[stateindex,nj]*atom.dipolemoments[stateindex,nj]
        detuningfactoreven = (1.0/(transitionEnergy-photonEnergy-0.5j*sc.hbar*atom.decayrates[nj,stateindex])+ 1.0/(transitionEnergy+photonEnergy+0.5j*sc.hbar*atom.decayrates[nj,stateindex]))
        detuningfactorodd = (1.0/(transitionEnergy-photonEnergy-0.5j*sc.hbar*atom.decayrates[nj,stateindex])- 1.0/(transitionEnergy+photonEnergy+0.5j*sc.hbar*atom.decayrates[nj,stateindex]))
        
        alpha0 += prefactor*complex(QuietWigner6j(1.0,0.0,1.0,atom.stateJ[stateindex],atom.stateJ[nj],atom.stateJ[stateindex]))*dipolesqrd*detuningfactoreven
        alpha1 -= prefactor*np.sqrt(3.0)*complex(QuietWigner6j(1.0,1.0,1.0,atom.stateJ[stateindex],atom.stateJ[nj],atom.stateJ[stateindex]))*dipolesqrd*detuningfactorodd
        alpha2 += prefactor*np.sqrt(5.0)*complex(QuietWigner6j(1.0,2.0,1.0,atom.stateJ[stateindex],atom.stateJ[nj],atom.stateJ[stateindex]))*dipolesqrd*detuningfactoreven
      
    alpha_ss = alpha0/np.sqrt(6.0*atom.stateJ[stateindex]+3.0)
    alpha_vv = -1.0*alpha1*np.sqrt(2.0*atom.stateJ[stateindex])/np.sqrt((atom.stateJ[stateindex]+1.0)*(2.0*atom.stateJ[stateindex]+1.0))
    alpha_tt = -1.0*alpha1*np.sqrt(2.0*atom.stateJ[stateindex]*(2.0*atom.stateJ[stateindex]-1.0))/np.sqrt(3.0*(atom.stateJ[stateindex]+1.0)*(2.0*atom.stateJ[stateindex]+1.0)*(2.0*atom.stateJ[stateindex]+3.0))
      
    return np.array([(alpha_ss)/polarizabilityunit,(alpha_vv)/polarizabilityunit,np.real(alpha_tt)/polarizabilityunit])
        
def polarizability(tab,row,column):
    #os.chdir("AtomicData")
    #Constants
    
    outputframe = GUI.MakeLabelFrame(frame=tab,row=0,column=1,width=200,height=150,text='')
    ###############Canvas and figure for plotting
    figure = plt.Figure(figsize=(7.9,7.4), dpi=100)
    figure.set_tight_layout(True)
    ax = figure.add_subplot(211)
    ax.set_title('Polarizability')
    ax3 = figure.add_subplot(212)
    ax3.set_title('Scattering rate')
    trapscanvas= FigureCanvasTkAgg(figure, outputframe)
    trapscanvas.get_tk_widget().grid(column = 0,row =0,padx = 5,pady = 5)
    #
    ax.plot([300,1600],[0,0])
    ax.set_xlim([300,1600])
    ax3.set_xlim([300,1600])
    ax.set_ylim([-1,1])
    ax.set_ylabel(r"AC polarizability ($e^2 a_0^2 /E_H$)")
    #ax2 = ax.twiny()
    wavelengthticks = ax.get_xticks()
    #ax.set_xticks([])
	#ax2.set_xlim([300,1600])
    #ax2.set_xticks(wavelengthticks,labels=np.around(sc.c*1.0e-12/(wavelengthticks*1.0e-9),decimals=3))
    #ax2.set_xlabel(r"Frequency (THz)")
    ax3.set_xlabel("Wavelength (nm)")
    figure.tight_layout()    
    
    toolbarframe = GUI.MakeFrame(frame=outputframe,row=1,column=0,sticky='ns')
    #toolbarframe.grid(height=50)
    toolbar = NavigationToolbar2Tk(trapscanvas,toolbarframe,pack_toolbar=False)
    toolbar.grid(column = 0,row =0,padx = 0,pady = 0, sticky=tk.EW)
    
    #The frame on the left
    settingsframe = GUI.MakeFrame(frame=tab,row=0,column=0,sticky='ns')
    
    
    
    ###########
    #Optical dipole trap parameters
    #beams_frame = GUI.MakeLabelFrame(frame=settingsframe,row=0,column=0,width=100,height=100,text='Optical traps')
    #Magnetic trap parameters
    #magnetictraps_frame = GUI.MakeLabelFrame(frame=settingsframe,row=1,column=0,width=50,height=50,text='Magnetic traps')
    #Atomic parameters
    ###########
    
    
    ########################################################################################################################################################################################
    ########################################################################################################################################################################################
    #############################################################################################The frame on the bottom 
    ########################################################################################################################################################################################
    ########################################################################################################################################################################################
    ###########
    ######################
    #subframe for the plotting settings
    atoms_frame = GUI.MakeLabelFrame(frame=settingsframe,row=1,column=0,width=50,height=50,text='Atomic details')
    atoms_frame.grid(padx=10)
    state_frame = GUI.MakeLabelFrame(frame=settingsframe,row=2,column=0,width=50,height=50,text='State')
    outputsettingsframe = GUI.MakeLabelFrame(frame=settingsframe,row=3,column=0,width=100,height=220,text='Grid details')
    outputsettingsframe.grid(padx=0)
    
    ##
    #Atomic details
    atombuttonwidth= 8
    statebuttonwidth= 14
    
    for i in np.arange(6):
        statefllerbutton1 = GUI.MakeButton(frame=state_frame,row=i%3,column=floor(i/3),width=statebuttonwidth,text=' ',command=lambda: null,state='normal',relief='raised')
    statevariable = 0
          
    buttonframe= GUI.MakeFrame(frame=atoms_frame,row=0,column=0,sticky='ns')
          
    species = ''    
    statebuttons = []
    
    def set_state(x,y):
        nonlocal statevariable
        nonlocal statebuttons
        statevariable = x
        for widget in state_frame.winfo_children():
            widget.config(relief='raised')     
        statebuttons[y].config(relief='sunken')     
        return
    
    def set_atom(atom,button):
        global species
        nonlocal atomfamily
        nonlocal statebuttons
        statebuttons = []
        species = atom
        atomdata = importlib.import_module("AtomicData."+species)
        
        for widget in state_frame.winfo_children():
            widget.destroy()
       
        choosablestates = atomdata.polarizability_plottable_states
        buttonindex = 0
        for i in np.arange(len(choosablestates)):
            statebuttons.append(GUI.MakeButton(frame=state_frame,row=i%3,column=floor(i/3),width=statebuttonwidth,text=atomdata.states[i],command=lambda i=i: set_state(choosablestates[i],i),state='normal',relief='raised'))
        
        for j in np.arange(len(choosablestates),6):
            statefllerbutton1 = GUI.MakeButton(frame=state_frame,row=j%3,column=floor(j/3),width=statebuttonwidth,text=' ',command=lambda: null,state='normal',relief='raised')
            
        for widget in atomchoicesframe.winfo_children():
            widget.config(relief='raised')                   
        button.config(relief='sunken')
        return
    
    
    atomchoicesframe = GUI.MakeLabelFrame(frame=buttonframe,row=1,column=0,width=50,height=50,text='')
    atomchoicesframe.grid(columnspan=3)
    #GUI.MakeFrame(frame=buttonframe,row=1,column=0,sticky='ns')


    LiButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='Li',command=lambda: set_atom('6Li',LiButton),state='normal',relief='raised')
    NaButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Na',command=lambda: set_atom('23Na',NaButton),state='normal',relief='raised')
    KButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text='K',command=lambda: set_atom('39K',KButton),state='normal',relief='raised')
    RbButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='Rb',command=lambda: set_atom('87Rb',RbButton),state='normal',relief='raised')
    CsButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='Cs',command=lambda: set_atom('133Cs',CsButton),state='normal',relief='raised')   


    CaButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='Ca',command=lambda: set_atom('Ca',CaButton),state='normal',relief='raised')
    SrButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='Sr',command=lambda: set_atom('88Sr',SrButton),state='normal',relief='raised')
   
    EuButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='Eu',command=lambda: set_atom('B-Eu',EuButton),state='normal',relief='raised')
    DyButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='Dy',command=lambda: set_atom('B-Dy',DyButton),state='normal',relief='raised')
    ErButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Er',command=lambda: set_atom('B-Er',ErButton),state='normal',relief='raised')
    YtButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text='Yt',command=lambda: set_atom('B-Yt',YtButton),state='normal',relief='raised')
       
    buttonlist = [LiButton,NaButton,KButton,RbButton,CsButton,CaButton,SrButton,EuButton,DyButton,ErButton,YtButton]
    
    
    atomfamily = ''
    def set_atom_choices(atomset):
        nonlocal atomfamily
        nonlocal LiButton
        nonlocal NaButton
        nonlocal KButton
        nonlocal RbButton
        nonlocal CsButton
        
        nonlocal CaButton
        nonlocal SrButton
        
        nonlocal EuButton
        nonlocal DyButton
        nonlocal ErButton
        nonlocal YtButton
        
        atomfamily = atomset
        
        for widget in atomchoicesframe.winfo_children():
            widget.destroy()
        
        for i in np.arange(6):
            fillerbutton = GUI.MakeButton(frame=atomchoicesframe,row=floor(i/3),column=i%3,width=atombuttonwidth,text=' ',command=lambda: null,state='normal',relief='raised')
                   
        if atomset == 'alkali':
            LiButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='Li',command=lambda: set_atom('6Li',LiButton),state='normal',relief='raised')
            NaButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Na',command=lambda: set_atom('23Na',NaButton),state='normal',relief='raised')
            KButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text='K',command=lambda: set_atom('39K',KButton),state='normal',relief='raised')
            RbButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='Rb',command=lambda: set_atom('87Rb',RbButton),state='normal',relief='raised')
            CsButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='Cs',command=lambda: set_atom('133Cs',CsButton),state='normal',relief='raised') 
                    
        if atomset == 'Alkaline earth':
            CaButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='Ca',command=lambda: set_atom('Ca',CaButton),state='normal',relief='raised')
            SrButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='Sr',command=lambda: set_atom('88Sr',SrButton),state='normal',relief='raised')
                    
        if atomset == 'Lanthanide':
            EuButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='Eu',command=lambda: set_atom('B-Eu',EuButton),state='normal',relief='raised')
            DyButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='Dy',command=lambda: set_atom('B-Dy',DyButton),state='normal',relief='raised')
            ErButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Er',command=lambda: set_atom('B-Er',ErButton),state='normal',relief='raised')
            YtButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text='Yt',command=lambda: set_atom('B-Yt',YtButton),state='normal',relief='raised')


        return
    
    alkalibuttonButton = GUI.MakeButton(frame=buttonframe,row=0,column=0,width=10,text='alkali',command=lambda: set_atom_choices('alkali'),state='normal',relief='raised')
    earthalkalibuttonButton = GUI.MakeButton(frame=buttonframe,row=0,column=1,width=10,text='Alkaline earth',command=lambda: set_atom_choices('Alkaline earth'),state='normal',relief='raised')
    LanthanidesbuttonButton = GUI.MakeButton(frame=buttonframe,row=0,column=2,width=10,text='Lanthanide',command=lambda: set_atom_choices('Lanthanide'),state='normal',relief='raised')
    set_atom_choices('alkali')
    ##############
    
    
    
    
    
    plotsettingssizes = [13,10,8]
    
    # Set grid properties
    #gridpoints
    variable_gridpoints, entry_gridpoints = GUI.FullEntry(frame=outputsettingsframe,row=0,column=0,sizes=plotsettingssizes,description="points = ",unit='')
    entry_gridpoints.insert(tk.END,1300.0)
    #min wavelength
    variable_lambdamin, entry_lambdamin = GUI.FullEntry(frame=outputsettingsframe,row=1,column=0,sizes=plotsettingssizes,description="lambda min = ",unit='nm')
    entry_lambdamin.insert(tk.END,300.0)
    #max wavelength
    variable_lambdamax, entry_lambdamax = GUI.FullEntry(frame=outputsettingsframe,row=2,column=0,sizes=plotsettingssizes,description="lambda max = ",unit='nm')
    entry_lambdamax.insert(tk.END,1600.0)
    #y-range
    variable_yrange, entry_yrange = GUI.FullEntry(frame=outputsettingsframe,row=3,column=0,sizes=plotsettingssizes,description="y range = ",unit=' e\u00b2 a\u2080\u00b2 / E\u2095')
    entry_yrange.insert(tk.END,40000.0)
    #
    variable_fcenter, entry_fcenter = GUI.FullEntry(frame=outputsettingsframe,row=4,column=0,sizes=plotsettingssizes,description="f center",unit='THz')
    entry_fcenter.insert(tk.END,300.0)
    #
    variable_fspan, entry_fspan = GUI.FullEntry(frame=outputsettingsframe,row=5,column=0,sizes=plotsettingssizes,description="f span",unit='GHz')
    entry_fspan.insert(tk.END,100.0)




    computebuttonframe = tk.Frame(outputsettingsframe)
    computebuttonframe.grid(column = 0,row =6)
    
    #
    
    
    
    # toggle to show the linear fit    
    holdplot = tk.IntVar()
    entry_holdplot = GUI.MakeCheckButton(frame=computebuttonframe,row=0,column=1,text='Hold plot',command=lambda: null,state='active',variable=holdplot)
    # toggle for which polarizabilities to plot
    
    scalarplot = tk.IntVar()
    entry_scalarplot = GUI.MakeCheckButton(frame=computebuttonframe,row=0,column=0,text='Scalar',command=lambda: null,state='active',variable=scalarplot)
    entry_scalarplot.select()
    vectorplot = tk.IntVar()
    entry_vectorplot = GUI.MakeCheckButton(frame=computebuttonframe,row=1,column=0,text='Vector',command=lambda: null,state='active',variable=vectorplot)
    tensorplot = tk.IntVar()
    entry_tensorplot = GUI.MakeCheckButton(frame=computebuttonframe,row=2,column=0,text='Tensor',command=lambda: null,state='active',variable=tensorplot)
    
    
    
    colorsubframe = tk.Frame(computebuttonframe)
    colorsubframe.grid(row=1,column=1)
    label_color = GUI.MakeLabel(frame=colorsubframe,row=0,column=0,width=7,text="color = ",anchor='e')
    plotcolor, plotcolordropdown= GUI.Make_DropDownMenu(colorsubframe,row=0,column=1,columnspan=1,rowspan=1,choices=['black','red','blue','green','purple','orange'])
    plotcolordropdown.grid(sticky='nw')
    
    ########################################################################################################################################################################################
            
          
    ################################################################################################################################################################################
        
    def plot_polarizability(xaxis):    
        nonlocal statevariable
        global species
        nonlocal ax 
        nonlocal trapscanvas
        nonlocal plotcolor
        #%% Reading entries
        gridpoints = float(entry_gridpoints.get())  
        minlambda = float(entry_lambdamin.get())  
        maxlambda = float(entry_lambdamax.get())  
        yrange = float(entry_yrange.get())
        fcenter = float(entry_fcenter.get())  
        fspan = float(entry_fspan.get())  
        fmin = fcenter*1.0e12 - 0.5*fspan*1.0e9
        fmax = fmin+fspan*1.0e9
        
        atomdata = importlib.import_module("AtomicData."+species)
        
        if xaxis=="lambda":
            lambdas = np.linspace(minlambda,maxlambda,int(gridpoints))
            freqs = sc.c/(lambdas*1.0e-9)
            alphaarray = compute_polarizability(species,lambdas,statevariable)
                    
            #%%making the plot
            if not holdplot.get():
                ax.cla()
                ax3.cla()
            if scalarplot.get():
                ax.plot(lambdas,np.real(alphaarray[0]),label=species+' '+atomdata.states[statevariable]+' scalar',color = plotcolor.get())
                ax3.plot(lambdas,np.imag(alphaarray[0]),label=species+' '+atomdata.states[statevariable]+' scalar',color = plotcolor.get())
                np.savetxt('alphaS.txt',[lambdas,np.real(alphaarray[0])])
            if vectorplot.get():
                ax.plot(lambdas,np.real(alphaarray[1]),label=species+' '+atomdata.states[statevariable]+' vector',color = plotcolor.get(),linestyle='dashed')
                ax3.plot(lambdas,np.imag(alphaarray[1]),label=species+' '+atomdata.states[statevariable]+' vector',color = plotcolor.get(),linestyle='dashed')
                np.savetxt('alphaV.txt',[lambdas,np.real(alphaarray[1])])
            if tensorplot.get():
                ax.plot(lambdas,np.real(alphaarray[2]),label=species+' '+atomdata.states[statevariable]+' tensor',color = plotcolor.get(),linestyle='dotted') 
            ax.grid(True)
            ax.set_ylim([-yrange,yrange])
            #ax2 = ax.twiny()
            wavelengthticks = ax.get_xticks()
            
            #ax2.set_xlim([minlambda,maxlambda])
            ax.set_xlim([minlambda,maxlambda])
            
            #ax2.set_xlim([freqs[0]*1.0e-12,freqs[len(freqs)-1]*1.0e-12])
            #ax2.set_xlabel(r"Frequency (THz)")
            #ax.set_xlabel("Wavelength (nm)")
            #ax.set_ylabel(r"AC polarizability ($e^2 a_0^2 /E_H$)")
                 
            
            #ax.set_ylabel("Trap depth (\u03bcK)")
            #ax.set_xlabel("Position (\u03bcm)")
            
            ax.legend()
            trapscanvas.draw()
               
        if xaxis=="f":
            freqs =  np.linspace(fmin,fmax,int(gridpoints))
            plotf = np.linspace(-fspan*0.5,fspan*0.5,int(gridpoints))
            
            lambdas = sc.c*1.0e9/(freqs) ## in nanometer
            alphaarray = compute_polarizability(species,lambdas,statevariable)
                    
            #%%making the plot
            if not holdplot.get():
                ax.cla()
                ax3.cla()
            if scalarplot.get():
                ax.plot(plotf,alphaarray[0],label=species+' '+atomdata.states[statevariable]+' scalar',color = plotcolor.get()) 
            if vectorplot.get():
                ax.plot(plotf,alphaarray[1],label=species+' '+atomdata.states[statevariable]+' vector',color = plotcolor.get(),linestyle='dashed') 
            if tensorplot.get():
                ax.plot(plotf,alphaarray[2],label=species+' '+atomdata.states[statevariable]+' tensor',color = plotcolor.get(),linestyle='dotted') 
            ax.grid(True)
            ax.set_ylim([-yrange,yrange])
            #ax2 = ax.twiny()
            #wavelengthticks = ax.get_xticks()
            
            #ax2.set_xlim([minlambda,maxlambda])
            ax.set_xlim([-fspan*0.5,fspan*0.5])
            
            #ax2.set_xlim([freqs[0]*1.0e-12,freqs[len(freqs)-1]*1.0e-12])
            ax.set_xlabel(r"Frequency (GHz)")
            #ax.set_xlabel("Wavelength (nm)")
            #ax.set_ylabel(r"AC polarizability ($e^2 a_0^2 /E_H$)")
                 
            
            #ax.set_ylabel("Trap depth (\u03bcK)")
            #ax.set_xlabel("Position (\u03bcm)")
            
            ax.legend()
            trapscanvas.draw()
        return
    
    
    
    computebutton = GUI.MakeButton(frame=computebuttonframe,row=2,column=1,width=12,text='Plot vs lambda',command=lambda: plot_polarizability("lambda"),state='normal',relief='raised')
    
    computebutton2 = GUI.MakeButton(frame=computebuttonframe,row=3,column=1,width=12,text='Plot vs f',command=lambda: plot_polarizability("f"),state='normal',relief='raised')
    
    #GUI.MakeLabel(frame=atoms_frame,row=5,column=0,width=labelwidth_atom,text=" ",anchor='e')
    
    
