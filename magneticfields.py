import tkinter as tk
from tkinter import ttk
import GUIfunctions as GUI
import numpy as np
from scipy import interpolate
from sympy import S
from numpy import linalg as LA
import scipy.constants as sc
import importlib
from AngularMomentumOperatorClass import Operator,AngMomOperators
from extrafunctions import subscript,superscript
import matplotlib.pyplot as plt
import time
from math import floor
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

def compute_states(I,J,ahf,bhf,chf,gI,gJ,m,Bmax,Bstep):
    #Setting Hilbert space dimension
    dimension = int(2.0*I+1)*int(2.0*J+1)
    
    #Setting up the magnetic field steps
    BTeslamax  = Bmax*1.0e-4
    BTeslastep = Bstep*1.0e-4
    points = int((BTeslamax)/BTeslastep)+1
    
    #Arrays for storing calculated energies, eigenvectors and magnetic moments
    energieslist = np.zeros((points,dimension))
    Eigenvectors = np.zeros((points+1,dimension,dimension),dtype=np.complex128)
    Eigenvectors[0] = (1.0+0.0j)*np.identity(dimension)
    momentslist = np.zeros((points-1,dimension))
        
    Ops = AngMomOperators(I,J)
    
    #Hyperfine hamiltonian
    H_hf = ahf * Ops.IdotJ.FmF
    if J > 0.5:
        QuadrupolePart = 3.0 * np.matmul(Ops.IdotJ.FmF,Ops.IdotJ.FmF) + 1.5 * Ops.IdotJ.FmF - np.matmul(Ops.Isqrd.FmF,Ops.Jsqrd.FmF)
        QuadrupolePart /= (2.0* I * (2.0*I - 1) * J * (2*J-1))
        
        OctupoleConstant = 0.0*Ops.Id.FmF 
        #OctupoleConstant /= (I * (I - 1)*(2.0*I - 1)  * J * (J-1)* (2*J-1))
    else:
        QuadrupolePart   = 0.0*Ops.Id.FmF
        OctupoleConstant = 0.0*Ops.Id.FmF 
    
    if bhf > 0.0:
        tempnumerator   = 3.0 * np.matmul(Ops.IdotJ.FmF,Ops.IdotJ.FmF) + 1.5 * Ops.IdotJ.FmF - np.matmul(Ops.Isqrd.FmF,Ops.Jsqrd.FmF)
        tempdenominator = 2.0 * I * (2.0*I - 1) * J * (2.0*J - 1)
        H_hf += chf * tempnumerator/tempdenominator


    #Eigenvectors = np.zeros((points+1,dimension),dtype=np.complex128)
    #newEigenvectors = (1.0+0.0j)*np.identity(dimension)

    tmpoperator = Operator(I,J)
                
    for i in np.arange(points):
    
        progressbar['value'] = 50.0+50.0*float(i)/points
        progressbar.update_idletasks()
        
        Bz = i*BTeslastep
        
        H = H_hf + muB * (gJ * Ops.Jz.FmF + gI * Ops.Iz.FmF) * Bz
        eigenvals, eigenvecscol = LA.eigh(H)
        eigenvecs = np.transpose(eigenvecscol)
                    
        for num in np.arange(dimension):
            overlaplist = np.absolute(np.dot(np.conjugate(eigenvecs),Eigenvectors[i,num]))
            bestmatchindex = np.argwhere(overlaplist==np.amax(overlaplist)) 
                                
            energieslist[i,num] = np.real(eigenvals[bestmatchindex])/sc.h
            Eigenvectors[i+1,num] = 1.0*eigenvecs[bestmatchindex]
                                       
    
    for j in np.arange(points-1):
        momentslist[j] =  - (sc.h/muB)*(energieslist[j+1]-energieslist[j])/BTeslastep
            
    Blist = np.arange(points)*Bstep + 1.0e-6
    
    levitationgradients = m*9.81/(momentslist*muB)
        
    #%% Set the hyperfine states in the transition frequency calculator
    tmpoperator = Operator(I,J)
    state1dropdown['menu'].delete(0, 'end')
    state2dropdown['menu'].delete(0, 'end')
    state1.set(tmpoperator.FmFindex[0])
    state2.set(tmpoperator.FmFindex[0])
    for choice in tmpoperator.FmFindex:
    
        state1dropdown['menu'].add_command(label=choice, command=tk._setit(state1, choice))
        state2dropdown['menu'].add_command(label=choice, command=tk._setit(state2, choice))

    PlotEnergybutton['state'] = 'normal'
    PlotMomentsbutton['state'] = 'normal'
    PlotLevGradbutton['state'] = 'normal'
    computetransitionenergiesbutton['state'] = 'normal'
    
    
    return energieslist, Eigenvectors, tmpoperator.FmFindex


       
def magnetic_fields(tab,row,column):
    
    magneticfields_frame = GUI.MakeFrame(frame=tab,column=column,row=row)
    
    ##########
    ##########
    muB  = sc.value(u'Bohr magneton') # Bohr magneton
    ##########
    graphframe = GUI.MakeLabelFrame(frame=magneticfields_frame,row=0,column=0,width=150,height=150,text='Output')
    graphframe.grid(rowspan=2)
    #
    figure = plt.Figure(figsize=(8.325,5.0), dpi=100)
    figure.set_tight_layout(True)
    ax = figure.add_subplot(111)
    ax.set_title('')
    trapscanvas= FigureCanvasTkAgg(figure, graphframe)
    trapscanvas.get_tk_widget().grid(column = 0,row =0,padx = 10,pady = 10)
    #
    ax.plot([0,100],[0,0])
    #ax.set_xlim([0,100])
    #ax.set_ylim([-1,1])
    ax.set_ylabel("Trap depth (\u03bcK)")
    ax.set_xlabel("B (G)")
    #figure.tight_layout()
    
    ##########
    settingsframe = GUI.MakeFrame(frame=magneticfields_frame,row=2,column=0)
    buttonframe= GUI.MakeFrame(frame=settingsframe,row=0,column=0,sticky='ns')
    buttonframe.grid(rowspan=2)
    atomchoicesframe = GUI.MakeLabelFrame(frame=buttonframe,row=1,column=0,width=50,height=50,text='')
    atomchoicesframe.grid(columnspan=3)
        
    state_frame = GUI.MakeLabelFrame(frame=buttonframe,row=2,column=0,width=50,height=50,text='State')
    state_frame.grid(columnspan=3)
        # atomsframe = GUI.MakeLabelFrame(frame=settingsframe,row=0,column=0,width=300,height=300,text='Insert atom details')
    
    #Atomic details
    atombuttonwidth= 8
    statebuttonwidth= 14
    
    for i in np.arange(6):
        statefllerbutton1 = GUI.MakeButton(frame=state_frame,row=i%3,column=floor(i/3),width=statebuttonwidth,text=' ',command=lambda: null,state='normal',relief='raised')
    statevariable = 0
    
    species = ''    
    statebuttons = []
    
    def set_state(x,y):
        global species
        nonlocal statevariable
        nonlocal statebuttons
        statevariable = x
        atomdata = importlib.import_module("AtomicData."+species)
        
        entry_I.delete(0, tk.END)  
        entry_J.delete(0, tk.END)  
        entry_ahf.delete(0, tk.END)  
        entry_bhf.delete(0, tk.END)  
        entry_chf.delete(0, tk.END)  
        entry_gI.delete(0, tk.END)  
        entry_gJ.delete(0, tk.END)  
        entry_mass.delete(0, tk.END)  
        
        entry_I.insert(tk.END, atomdata.I)
        entry_mass.insert(tk.END, atomdata.m/sc.physical_constants["atomic mass constant"][0])
        entry_gI.insert(tk.END, atomdata.gI)
        entry_ahf.insert(tk.END, np.around(atomdata.hyperfineA[statevariable]*1.0e-6/sc.h,decimals=6))
        entry_bhf.insert(tk.END, np.around(atomdata.hyperfineB[statevariable]*1.0e-6/sc.h,decimals=8))
        entry_chf.insert(tk.END, atomdata.hyperfineC[statevariable]/sc.h)
        entry_gJ.insert(tk.END, atomdata.gfactors[statevariable])
        entry_J.insert(tk.END, atomdata.stateJ[statevariable])
        
        
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
        
    
    #Atomic details
    atombuttonwidth= 8
    atompropssizes = [22,12,9]
    ##
    
    atomfamily = ''
    def set_atom_choices(atomset):
                
        for widget in atomchoicesframe.winfo_children():
            widget.destroy()
        for i in np.arange(9):
            fillerbutton = GUI.MakeButton(frame=atomchoicesframe,row=floor(i/3),column=i%3,width=atombuttonwidth,text=' ',command=lambda: null,state='normal',relief='raised')       
               
        if atomset == 'alkali':
            Li6Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text=superscript(6)+'Li',command=lambda: set_atom('6Li',Li6Button),state='normal',relief='raised')
            Li7Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text=superscript(7)+'Li',command=lambda: set_atom('7Li',Li6Button),state='normal',relief='raised')
            Na23Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text=superscript(2)+superscript(3)+'Na',command=lambda: set_atom('23Na',Na23Button),state='normal',relief='raised')
            K39Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text=superscript(3)+superscript(9)+'K',command=lambda: set_atom('39K',K39Button),state='normal',relief='raised')
            K40Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text=superscript(4)+superscript(0)+'K',command=lambda: set_atom('40K',K40Button),state='normal',relief='raised')
            K41Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=2,width=atombuttonwidth,text=superscript(4)+superscript(1)+'K',command=lambda: set_atom('41K',K40Button),state='normal',relief='raised')
            Rb85Button = GUI.MakeButton(frame=atomchoicesframe,row=2,column=0,width=atombuttonwidth,text=superscript(8)+superscript(5)+'Rb',command=lambda: set_atom('85Rb',Rb85Button),state='normal',relief='raised')
            Rb87Button = GUI.MakeButton(frame=atomchoicesframe,row=2,column=1,width=atombuttonwidth,text=superscript(8)+superscript(7)+'Rb',command=lambda: set_atom('87Rb',Rb87Button),state='normal',relief='raised')
            CsButton = GUI.MakeButton(frame=atomchoicesframe,row=2,column=2,width=atombuttonwidth,text=superscript(1)+superscript(3)+superscript(3)+'Cs',command=lambda: set_atom('133Cs',CsButton),state='normal',relief='raised')
            
        if atomset == ' ':
            return
            #Ca40Button = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='40Ca',command=lambda: set_atom('40Ca',Ca40Button),state='normal',relief='raised')
            #SrBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Bos Sr',command=lambda: set_atom('B-Sr',SrBosButton),state='normal',relief='raised')
            #Sr87Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='87Sr',command=lambda: set_atom('87Sr',Sr87Button),state='normal',relief='raised')
            
        if atomset == '  ':
            return
            #EuBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=0,width=atombuttonwidth,text='Bos Eu',command=lambda: set_atom('B-Eu',EuBosButton),state='normal',relief='raised')
            #DyBosButton = GUI.MakeButton(frame=atomchoicesframe,row=1,column=0,width=atombuttonwidth,text='Bos Dy',command=lambda: set_atom('B-Dy',DyBosButton),state='normal',relief='raised')
            #DyFerButton = GUI.MakeButton(frame=atomchoicesframe,row=2,column=0,width=atombuttonwidth,text='Fer Dy',command=lambda: set_atom('F-Dy',DyFerButton),state='normal',relief='raised')
            #ErBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=1,width=atombuttonwidth,text='Bos Er',command=lambda: set_atom('B-Er',ErBosButton),state='normal',relief='raised')
            #Er167Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=1,width=atombuttonwidth,text='F-167Er',command=lambda: set_atom('F-167Er',Er167Button),state='normal',relief='raised')
            #YtBosButton = GUI.MakeButton(frame=atomchoicesframe,row=0,column=2,width=atombuttonwidth,text='Bos Yt',command=lambda: set_atom('B-Yt',YtBosButton),state='normal',relief='raised')
            #Yt171Button = GUI.MakeButton(frame=atomchoicesframe,row=1,column=2,width=atombuttonwidth,text='171Yt',command=lambda: set_atom('171Yt',Yt171Button),state='normal',relief='raised')
            #Yt173Button = GUI.MakeButton(frame=atomchoicesframe,row=2,column=2,width=atombuttonwidth,text='173Yt',command=lambda: set_atom('173Yt',Yt173Button),state='normal',relief='raised')
        
        return
    
    alkalibuttonButton = GUI.MakeButton(frame=buttonframe,row=0,column=0,width=10,text='alkali',command=lambda: set_atom_choices('alkali'),state='normal',relief='raised')
    earthalkalibuttonButton = GUI.MakeButton(frame=buttonframe,row=0,column=1,width=10,text=' ',command=lambda: set_atom_choices(' '),state='normal',relief='raised')
    LanthanidesbuttonButton = GUI.MakeButton(frame=buttonframe,row=0,column=2,width=10,text='  ',command=lambda: set_atom_choices('  '),state='normal',relief='raised')
    set_atom_choices('alkali')
    
    ##########
    toolbar = NavigationToolbar2Tk(trapscanvas,settingsframe,pack_toolbar=False)
    toolbar.grid(column = 1,row =0,padx = 10,pady = 10, sticky=tk.EW)
        
    variables_frame = GUI.MakeLabelFrame(frame=settingsframe,row=1,column=1,width=150,height=150,text='Parameters')
    variablessizes = [40,15,6] 
    
    Ivariable, entry_I = GUI.FullEntry(frame=variables_frame,row=0,column=0,sizes=variablessizes,description="Nuclear spin I = ",unit="")
    Jvariable, entry_J = GUI.FullEntry(frame=variables_frame,row=1,column=0,sizes=variablessizes,description="Electronic angular momentum J = ",unit="")
    ahfvariable, entry_ahf = GUI.FullEntry(frame=variables_frame,row=2,column=0,sizes=variablessizes,description="Magnetic dipole constant = ",unit="MHz ")
    bhfvariable, entry_bhf = GUI.FullEntry(frame=variables_frame,row=3,column=0,sizes=variablessizes,description="Electric quadrupole constant = ",unit="MHz ")
    chfvariable, entry_chf = GUI.FullEntry(frame=variables_frame,row=4,column=0,sizes=variablessizes,description="octupole (not used) = ",unit="Hz ")
    gIvariable, entry_gI = GUI.FullEntry(frame=variables_frame,row=5,column=0,sizes=variablessizes,description="Nuclear g-factor = ",unit="")
    gJvariable, entry_gJ = GUI.FullEntry(frame=variables_frame,row=6,column=0,sizes=variablessizes,description="Fine structure Lande g-factor = ",unit="")
    massvariable, entry_mass = GUI.FullEntry(frame=variables_frame,row=6,column=0,sizes=variablessizes,description="Mass   m = ",unit="u  ")
    
    ##########
    plotparamsframe = GUI.MakeLabelFrame(frame=settingsframe,row=0,column=2,width=150,height=150,text='Plot settings')
    plotparamsframe.grid(rowspan=2)
    
    paramssizes = [15,10,3] 
    Bmaxvariable, entry_Bmax = GUI.FullEntry(frame=plotparamsframe,row=0,column=0,sizes=paramssizes,description="Bmax = ",unit=" G")
    entry_Bmax.insert(tk.END,200.0)
    Bstepvariable, entry_Bstep = GUI.FullEntry(frame=plotparamsframe,row=1,column=0,sizes=paramssizes,description="Bstep = ",unit=" G")
    entry_Bstep.insert(tk.END,1.0)
    
    GUI.MakeLabel(plotparamsframe,row=2,column=0,width=10,text='                                 ',anchor='center')
    
    plotbuttonwidth = 20
    
    progressbar = ttk.Progressbar(plotparamsframe, orient='horizontal', length=150, mode='determinate')
    progressbar.grid(row=7)
    
    
    ##########
    transitionsframe = GUI.MakeLabelFrame(frame=magneticfields_frame,row=0,column=1,width=150,height=150,text='Transition frequencies')
    transitionsframe.grid(rowspan=2)
    
    Bvariable, entry_B = GUI.FullEntry(frame=transitionsframe,row=1,column=0,sizes=[3,10,3],description=" B = ",unit=" G")
    entry_B.insert(tk.END,0.0)
    transitionssubframe = tk.Frame(transitionsframe)
    transitionssubframe.grid(row=0,column=0,columnspan=2)
    
    tmp = GUI.MakeLabel(transitionssubframe,row=0,column=1,width=10,text='<->',anchor='center')
    state1, state1dropdown = GUI.Make_DropDownMenu(frame=transitionssubframe,row=0,column=0,rowspan=1,columnspan=1,choices=[[0,0]])
    state2, state2dropdown = GUI.Make_DropDownMenu(frame=transitionssubframe,row=0,column=2,rowspan=1,columnspan=1,choices=[[0,0]])
    
    transitionEnergyoutput = GUI.MakeText(transitionsframe,row=2,column=0,width=35,height=25)
    transitionEnergyoutput.insert('end','')
    transitionEnergyoutput.grid(columnspan=2)
    transitionEnergyoutput.tag_config('justified', justify='left')
    
    #state1, state1dropdown= create_dropdown(transitionssubframe,0,0,1,1,[[0,0]])
    #state2, state2dropdown= create_dropdown(transitionssubframe,0,2,1,1,[[0,0]])
    def compute():
            global I
            global J
            global Blist
            global Eigenvectors
            global energieslist
            global momentslist
            global levitationgradients
            global Ops
            global H_hf
            global I
            global J
            global ahf
            global bhf
            global chf
            global gI
            global gJ
            nonlocal state1
            nonlocal state2
            nonlocal state1dropdown
            nonlocal state2dropdown
            nonlocal progressbar
            
            progressbar['value'] = 0.0
            progressbar.update_idletasks()
            
            I   = float(entry_I.get())
            J   = float(entry_J.get())
            ahf = float(entry_ahf.get())*sc.h*1.0e6
            bhf = float(entry_bhf.get())*sc.h*1.0e6
            chf = float(entry_chf.get())*sc.h
            gI  = float(entry_gI.get())
            gJ  = float(entry_gJ.get())
            m  = float(entry_mass.get())*sc.physical_constants["atomic mass constant"][0]
            
            dimension = int(2.0*I+1)*int(2.0*J+1)
            
            Bmax = float(entry_Bmax.get())
            Bstep = float(entry_Bstep.get())
            BTeslamax  = Bmax*1.0e-4
            BTeslastep = Bstep*1.0e-4
            points = int((BTeslamax)/BTeslastep)
            
            energieslist = np.zeros((points,dimension))
            Eigenvectors = np.zeros((points+1,dimension,dimension),dtype=np.complex128)
            Eigenvectors[0] = (1.0+0.0j)*np.identity(dimension)
            momentslist = np.zeros((points-1,dimension))
            
            progressbar['value'] = 50.0
            progressbar.update_idletasks()
            
            Ops = AngMomOperators(I,J)
            
            tmpoperator = Operator(I,J)
            
            
            #Use Breit-Rabi formula for the ground state
            if False:
                def BreitRabi(b,m,s):
                    ehfs = ahf*(I+0.5)
                    x = (gJ-gI)*muB*b/ehfs
                    if m == I+0.5:
                        return ehfs*I/(2.0*I+1.0) + 0.5*(gJ + 2.0*I*gI)*muB*b
                    if -m == I+0.5:
                        return ehfs*I/(2.0*I+1.0) - 0.5*(gJ + 2.0*I*gI)*muB*b
                    return -0.5*ehfs/(2.0*I+1.0) + gI*muB*m*b + s*0.5*ehfs*np.sqrt(1+4.0*m*x/(2.0*I+1.0) + x*x)
                
				
                for i in np.arange(points):
                    progressbar['value'] = 50.0+50.0*float(i)/points
                    progressbar.update_idletasks()
                    
                    Bz = i*BTeslastep
                    
                    statenum =0
                    for f in np.arange(I-J,I+J+1,1):
                        for mf in np.arange(-f,f+1,1):
                            s = -1
                            if f == I+J:
                                s=1
                            
                            energieslist[i,statenum] = BreitRabi(Bz,mf,s)/sc.h
                            statenum += 1
                    
            
            else:
                #Hyperfine hamiltonian
                if J > 0.5:
                    QuadrupolePart = 3.0 * bhf * np.matmul(Ops.IdotJ.FmF,Ops.IdotJ.FmF) + 1.5 * Ops.IdotJ.FmF - Ops.Id.FmF * I*(I+1.0)*J*(J+1.0)
                    QuadrupolePart /= (2.0* I * (2.0*I - 1) * J * (2*J-1))
                    OctupoleConstant = 0.0*Ops.Id.FmF 
                    #OctupoleConstant /= (I * (I - 1)*(2.0*I - 1)  * J * (J-1)* (2*J-1))
                else:
                    QuadrupolePart   = 0.0*Ops.Id.FmF
                    OctupoleConstant = 0.0*Ops.Id.FmF
                    
                H_hf = ahf * Ops.IdotJ.FmF 
                if bhf > 0.0:
                    tempnumerator   = 3.0 * np.matmul(Ops.IdotJ.FmF,Ops.IdotJ.FmF) + 1.5 * Ops.IdotJ.FmF - Ops.Id.FmF * I*(I+1.0)*J*(J+1.0)
                    tempdenominator = 2.0 * I * (2.0*I - 1) * J * (2.0*J - 1)
                    H_hf += bhf * tempnumerator/tempdenominator

                for i in np.arange(points):
                
                    progressbar['value'] = 50.0+50.0*float(i)/points
                    progressbar.update_idletasks()
                    
                    Bz = i*BTeslastep  + 1.0e-6 #A small initial fields splits the states so the eigenvectors are not mixed up.
                    
                    H = H_hf + muB * (gJ * Ops.Jz.FmF + gI * Ops.Iz.FmF) * Bz
                    eigenvals, eigenvecscol = LA.eigh(H)
                                
                    for num in np.arange(dimension):
                        if i==0:
                            overlaplist = np.absolute(np.dot((1.0+0.0j)*np.identity(dimension)[num],eigenvecscol))
                        else:
                            overlaplist = np.absolute(np.dot(Eigenvectors[i-1,num],eigenvecscol))
                            
                        if np.amax(overlaplist) < 0.5:
                            print("bad overlaps!")
                        bestmatchindex = np.argwhere(overlaplist==np.amax(overlaplist))
                                            
                        energieslist[i,num] = np.real(eigenvals[bestmatchindex[0,0]])/sc.h
                        #print(eigenvecs[bestmatchindex[0,0]])
                        Eigenvectors[i,num] = 1.0*eigenvecscol[:,bestmatchindex[0,0]]
                                       
                
                
                htmp = H_hf + muB * (gJ * Ops.Jz.FmF + gI * Ops.Iz.FmF) * 10.0*BTeslastep
                
                print(np.amax(np.matmul(Ops.Fz.FmF,htmp)-np.matmul(htmp,Ops.Fz.FmF)))
                
                #print(Eigenvectors[1,14,14])
                #plt.clf()
                #plt.plot(Eigenvectors[:,14,14],'o',linestyle='none',label='pop')
                #plt.legend()
                #plt.show()
                
                #print()
                #print(Eigenvectors[1,12])
                #print(tmpoperator.FmFindex[12])

                #propagate(Eigenvectors[0,13],energieslist[0,13], H_hf + muB * (gJ * Ops.Jz.FmF + gI * Ops.Iz.FmF) * BTeslastep,dimension,100.0)

            
            for j in np.arange(points-1):
                momentslist[j] =  - (sc.h/muB)*(energieslist[j+1]-energieslist[j])/BTeslastep
            
            
            #print(H_hf)
            #for i in np.arange(dimension):
            #    print(tmpoperator.FmFindex[i,1])
            #    vect = np.zeros(dimension)
            #    vect[i] = 1
            #    print(np.dot(np.conjugate(vect),np.dot(Ops.Fz.FmF,vect)))
                
            
            
            Blist = np.arange(points)*Bstep + 1.0e-6
            
            levitationgradients = m*9.81/(momentslist*muB)
            
            
            
            #%% Set the hyperfine states in the transition frequency calculator
            #tmpoperator = Operator(I,J)
            state1dropdown['menu'].delete(0, 'end')
            state2dropdown['menu'].delete(0, 'end')
            state1.set(tmpoperator.FmFindex[0])
            state2.set(tmpoperator.FmFindex[0])
            for choice in tmpoperator.FmFindex:
            
                state1dropdown['menu'].add_command(label=choice, command=tk._setit(state1, choice))
                state2dropdown['menu'].add_command(label=choice, command=tk._setit(state2, choice))

            PlotEnergybutton['state'] = 'normal'
            PlotMomentsbutton['state'] = 'normal'
            PlotLevGradbutton['state'] = 'normal'
            computetransitionenergiesbutton['state'] = 'normal'


            return    


    #%% plotting styles
    plotcolors = ['black','maroon','red','darkorange','yellowgreen','darkgreen','deepskyblue','navy','purple']
    plotlinestyles = ['solid','dotted','dashed','dashdot']
    plotstylescolor,plotstyleslines  = np.meshgrid(plotcolors,plotlinestyles)
    
    #%%
    def plot_energies():
        global I
        global J
        dimension = int(2.0*I+1)*int(2.0*J+1) 
                
        ax.cla()
        ax.set_title('Hyperfine + Zeeman Energies')
        ax.set_ylabel('Frequency (MHz)')
        ax.set_xlabel('B (Gauss)')
        ax.grid(True)
        
        tmpoperator = Operator(I,J)
        for d in np.arange(dimension):
            zerofieldstate = tmpoperator.FmFindex[d]
            
            ax.plot(Blist,(1.0e-6)*energieslist[:,d],linestyle = plotstyleslines.flatten()[d],color = plotstylescolor.flatten()[d], label= 'F,mF='+str(zerofieldstate[0])+','+str(zerofieldstate[1]))#,'.',markersize=0.6) 
        
        ax.legend(loc='upper right',fontsize='x-small',ncol=int(dimension/5))
        figure.tight_layout()
        trapscanvas.draw()
        return 

    
    def export():
        global I
        global J
        dimension = int(2.0*I+1)*int(2.0*J+1) 
        
        tmpoperator = Operator(I,J)
        statelabels = tmpoperator.FmFindex[:]
        
        head = ""
        for state in statelabels:
            head += str(state)+'\n'


        exportdata = np.zeros((energieslist.shape[0],dimension+1))
        exportdata[:,0] = Blist
        exportdata[:,1:] = (1.0e-6)*energieslist
		
        np.savetxt('ZeemanHyperfineData.txt',exportdata,header=head,comments='')
        print("Exported")
        return 	

    def plot_moments():
        global I
        global J
        dimension = int(2.0*I+1)*int(2.0*J+1) 
        
        ax.cla()
        ax.set_title('Magnetic moments')
        ax.set_ylabel('Moment (\u03bcB)')
        ax.set_xlabel('B (Gauss)')
        ax.grid(True)
        
        tmpoperator = Operator(I,J)
        for d in np.arange(dimension):
            zerofieldstate = tmpoperator.FmFindex[d]
            
            ax.plot(Blist[:-1],momentslist[:,d],linestyle = plotstyleslines.flatten()[d],color = plotstylescolor.flatten()[d], label= 'F,mF='+str(zerofieldstate[0])+','+str(zerofieldstate[1]))#,'.',markersize=0.6) 
        
        ax.legend(loc='upper right',fontsize='x-small',ncol=int(dimension/5))
        figure.tight_layout()
        trapscanvas.draw()
        return 
        
    
    def plot_levitationgradient():
        global I
        global J
        dimension = int(2.0*I+1)*int(2.0*J+1) 
        
        ax.cla()
        ax.set_title('Levitation gradient')
        ax.set_ylabel('G/cm')
        ax.set_xlabel('B (Gauss)')
        ax.grid(True)
        ax.set_ylim([-100.0,100.0])
        
        tmpoperator = Operator(I,J)
        for d in np.arange(dimension):
            zerofieldstate = tmpoperator.FmFindex[d]
            
            ax.plot(Blist[:-1],levitationgradients[:,d]*100.0,linestyle = plotstyleslines.flatten()[d],color = plotstylescolor.flatten()[d], label= 'F,mF='+str(zerofieldstate[0])+','+str(zerofieldstate[1]))#,'.',markersize=0.6) 
        
        ax.legend(loc='upper right',fontsize='x-small',ncol=int(dimension/5))
        figure.tight_layout()
        trapscanvas.draw()
        return 
            
    
    def give_transition_energy():
        firststate  = state1.get()
        secondstate = state2.get()
        
        tmpoperator = Operator(I,J)
        dimension = len(tmpoperator.FmFindex)
        indices = tmpoperator.FmFindex
                
        BB   = float(entry_B.get())  +1.0e-6
        
        closestBindex = np.argmin(np.absolute(BB-Blist))
        
        for d in np.arange(dimension):
            
            if str(indices[d]) == firststate:
                firstindex = d
            
            if str(indices[d]) == secondstate:
                secondindex = d
        
        interpolatingfunctionfirststate = interpolate.interp1d(Blist, energieslist[:,firstindex])
        interpolatingfunctionsecondstate = interpolate.interp1d(Blist, energieslist[:,secondindex])
                
        interpolatingmomentfunctionfirststate = interpolate.interp1d(Blist[:-1], momentslist[:,firstindex])
        interpolatingmomentfunctionsecondstate = interpolate.interp1d(Blist[:-1], momentslist[:,secondindex])
        
        firstenergy = interpolatingfunctionfirststate(BB)
        secondenergy = interpolatingfunctionsecondstate(BB)
        
        firstmoment = interpolatingmomentfunctionfirststate(BB)
        secondmoment = interpolatingmomentfunctionsecondstate(BB)        
        
		####### 
        state1eigenvector = Eigenvectors[closestBindex+1,firstindex]
        state2eigenvector = Eigenvectors[closestBindex+1,secondindex]
        
        H_hf = ahf * Ops.IdotJ.FmF
        if J > 0.5:
            QuadrupolePart = 3.0 * np.matmul(Ops.IdotJ.FmF,Ops.IdotJ.FmF) + 1.5 * Ops.IdotJ.FmF - np.matmul(Ops.Isqrd.FmF,Ops.Jsqrd.FmF)
            QuadrupolePart /= (2.0* I * (2.0*I - 1) * J * (2*J-1))
            
            OctupoleConstant = 0.0*Ops.Id.FmF 
            #OctupoleConstant /= (I * (I - 1)*(2.0*I - 1)  * J * (J-1)* (2*J-1))
        else:
            QuadrupolePart   = 0.0*Ops.Id.FmF
            OctupoleConstant = 0.0*Ops.Id.FmF
        if bhf > 0.0:
            tempnumerator   = 3.0 * np.matmul(Ops.IdotJ.FmF,Ops.IdotJ.FmF) + 1.5 * Ops.IdotJ.FmF - np.matmul(Ops.Isqrd.FmF,Ops.Jsqrd.FmF)
            tempdenominator = 2.0 * I * (2.0*I - 1) * J * (2.0*J - 1)
            H_hf += chf * tempnumerator/tempdenominator
        H = H_hf + muB * (gJ * Ops.Jz.FmF + gI * Ops.Iz.FmF) * BB*1.0e-4
        
        # firststateeigsindex = np.argmin(np.absolute(eigenvals/sc.h-firstenergy))
        # secondstateeigsindex = np.argmin(np.absolute(eigenvals/sc.h-secondenergy))
        # print(firststateeigsindex,secondstateeigsindex)
        # print('\n')
        # print(eigenvecs[firststateeigsindex])
        Jz1to1 = np.dot(np.conjugate(state1eigenvector),np.dot(Ops.Jz.FmF,state1eigenvector))
        Jz2to2 = np.dot(np.conjugate(state2eigenvector),np.dot(Ops.Jz.FmF,state2eigenvector))
        Jx1to2 = np.dot(np.conjugate(state1eigenvector),np.dot(Ops.Jx.FmF,state2eigenvector))
        Jy1to2 = np.dot(np.conjugate(state1eigenvector),np.dot(Ops.Jy.FmF,state2eigenvector))
        Jz1to2 = np.dot(np.conjugate(state1eigenvector),np.dot(Ops.Jz.FmF,state2eigenvector))
        Ix1to2 = np.dot(np.conjugate(state1eigenvector),np.dot(Ops.Ix.FmF,state2eigenvector))
        Iy1to2 = np.dot(np.conjugate(state1eigenvector),np.dot(Ops.Iy.FmF,state2eigenvector))
        Iz1to2 = np.dot(np.conjugate(state1eigenvector),np.dot(Ops.Iz.FmF,state2eigenvector))
        
        #print(state1eigenvector)
        #print(state2eigenvector)
        
        E1 = np.absolute(np.dot(np.conjugate(state1eigenvector),np.dot(H,state1eigenvector)))
        E2 = np.absolute(np.dot(np.conjugate(state2eigenvector),np.dot(H,state2eigenvector)))
        
        #print('test')
        #print(firstindex,secondindex)
        #print(firstenergy,secondenergy)
        #print(E1/sc.h,E2/sc.h)
        # #print(Jx1to2)
        # print('\n')
		
        outpt = str(firststate)+': '+str(np.around(firstenergy*1.0e-6,decimals=3))+' MHz, '
        outpt += str(np.around(firstmoment,decimals=2))+' \u03bcB \n'
        outpt += 'Zeeman shift = '+str(np.around((firstenergy-interpolatingfunctionfirststate(1.0e-6))*1.0e-6,decimals=3))+' MHz \n'
        outpt += '<1|Jz|1> = '+str(np.around(Jz1to1,decimals=3))+' \n'
        outpt += str(secondstate)+': '+str(np.around(secondenergy*1.0e-6,decimals=3))+' MHz, '
        outpt += str(np.around(secondmoment,decimals=2))+' \u03bcB \n'
        outpt += 'Zeeman shift = '+str(np.around((secondenergy-interpolatingfunctionsecondstate(1.0e-6))*1.0e-6,decimals=3))+' MHz \n'
        outpt += '<2|Jz|2> = '+str(np.around(Jz2to2,decimals=3))+' \n \n'
        #outpt += str(secondstate)+' : '+str(np.around(secondenergy*1.0e-6,decimals=3))+' MHz \n'
        outpt += 'diff = '+str(np.around((secondenergy-firstenergy)*1.0e-6,decimals=3))+' MHz, \n'
        outpt += str(np.around(firstmoment-secondmoment,decimals=2))+' \u03bcB = '+str(np.around((firstmoment-secondmoment)*1.399,decimals=2))+str(' MHz/G \n')
        outpt += '<1|Jx|2> = '+str(np.around(Jx1to2,decimals=3))+' \n'
        outpt += '<1|Jy|2> = '+str(np.around(Jy1to2,decimals=3))+' \n'
        outpt += '<1|Jz|2> = '+str(np.around(Jz1to2,decimals=3))+' \n'
        outpt += '<1|Ix|2> = '+str(np.around(Ix1to2,decimals=3))+' \n'
        outpt += '<1|Iy|2> = '+str(np.around(Iy1to2,decimals=3))+' \n'
        outpt += '<1|Iz|2> = '+str(np.around(Iz1to2,decimals=3))+' \n'
        
        transitionEnergyoutput.delete(1.0,tk.END)
        transitionEnergyoutput.insert(tk.END,outpt, 'justified')
        
        return

    
    computebutton = GUI.MakeButton(frame=plotparamsframe,row=3,column=0,width=plotbuttonwidth,text='Compute',command=lambda: compute(),state='normal',relief='raised')
    PlotEnergybutton = GUI.MakeButton(frame=plotparamsframe,row=4,column=0,width=plotbuttonwidth,text='Plot energies',command=lambda: plot_energies(),state='disable',relief='raised')
    PlotMomentsbutton = GUI.MakeButton(frame=plotparamsframe,row=5,column=0,width=plotbuttonwidth,text='Plot magnetic moments',command=lambda: plot_moments(),state='disable',relief='raised')
    PlotLevGradbutton = GUI.MakeButton(frame=plotparamsframe,row=6,column=0,width=plotbuttonwidth,text='Plot levitation gradient',command=lambda: plot_levitationgradient(),state='disable',relief='raised')
    exportbutton = GUI.MakeButton(frame=plotparamsframe,row=8,column=0,width=plotbuttonwidth,text='Export',command=lambda: export(),state='normal',relief='raised')
    
    computetransitionenergiesbutton = GUI.MakeButton(frame=transitionsframe,row=1,column=1,width=10,text='Print',command=lambda: give_transition_energy(),state='disable',relief='raised')
    
    #add_button(transitionfrequencies,'Print' , give_transition_energy ,3,0,1,1,20)
    
    return
