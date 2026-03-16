import tkinter as tk
from tkinter import ttk
import GUIfunctions as GUI
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
from matplotlib.pyplot import cm
import importlib
import scipy
import scipy.constants as sc
import numpy as np
import importlib
from sympy import N
from scipy import optimize
from numpy import linalg as LA

def lattices(tab,row,column):

    main_lattices_frame = GUI.MakeFrame(frame=tab,column=column,row=row)

    leftframe = GUI.MakeFrame(frame=main_lattices_frame,row=0,column=0,sticky='ns')
   
    ##### Settings for lattice calculation
    settingsframe = GUI.MakeLabelFrame(frame=leftframe,row=0,column=0,width=150,height=150,text='Input')
    
    # entry_depth =           name_entry_unit(latticetrapsframe,1,'Depth           = ','  \u03bcK',10)
    # entry_wavelength =      name_entry_unit(latticetrapsframe,2,'Wavelength      = ',' nm',10)
    # entry_wavelength.insert(tk.END, 1064.0)
    # entry_mass =            name_entry_unit(latticetrapsframe,3,'mass            = ','    u',10)
    
    # separator = ttk.Separator(latticetrapsframe, orient='horizontal').grid(row=4,column = 0, sticky="ew" )
    # entry_momentumcutoff =  name_entry_unit(latticetrapsframe,5,'momentum cutoff =   ',' 2 k_L',5)
    # entry_momentumcutoff.insert(tk.END, 30)
    # entry_unitcells =       name_entry_unit(latticetrapsframe,6,'unit cells      = 2x','+1         ',5)
    # entry_unitcells.insert(tk.END, 100)
    labelwidth_plotsettings = 20
    entrywidth_plotsettings = 8
    unitwidth_plotsettings = 5
    #depth
    activerow = 0
    label_depth = GUI.MakeLabel(frame=settingsframe,row=activerow,column=0,width=labelwidth_plotsettings,text="depth = ",anchor='e')
    depthvariable = ''
    entry_depth = GUI.MakeEntry(frame=settingsframe,row=activerow,column=1,width=entrywidth_plotsettings,variable=depthvariable,state='normal')
    unit_depth = GUI.MakeLabel(frame=settingsframe,row=activerow,column=2,width=unitwidth_plotsettings,text="\u03bcK",anchor='e')
    entry_depth.insert(tk.END,5.0)
    activerow += 1
    #wavelength 
    label_wavelength = GUI.MakeLabel(frame=settingsframe,row=activerow,column=0,width=labelwidth_plotsettings,text="wavelength = ",anchor='e')
    wavelengthvariable = ''
    entry_wavelength = GUI.MakeEntry(frame=settingsframe,row=activerow,column=1,width=entrywidth_plotsettings,variable=wavelengthvariable,state='normal')
    unit_wavelength = GUI.MakeLabel(frame=settingsframe,row=activerow,column=2,width=unitwidth_plotsettings,text="nm",anchor='e')
    entry_wavelength.insert(tk.END,1064.0)
    activerow += 1
    #mass 
    label_mass = GUI.MakeLabel(frame=settingsframe,row=activerow,column=0,width=labelwidth_plotsettings,text="mass = ",anchor='e')
    massvariable = ''
    entry_mass = GUI.MakeEntry(frame=settingsframe,row=activerow,column=1,width=entrywidth_plotsettings,variable=massvariable,state='normal')
    unit_mass = GUI.MakeLabel(frame=settingsframe,row=activerow,column=2,width=unitwidth_plotsettings,text=" u ",anchor='e')
    entry_mass.insert(tk.END,39.0)
    activerow += 1	
    #scattering length 
    label_scatteringlength = GUI.MakeLabel(frame=settingsframe,row=activerow,column=0,width=labelwidth_plotsettings,text="a    = ",anchor='e')
    scatteringlengthvariable = ''
    entry_scatteringlength = GUI.MakeEntry(frame=settingsframe,row=activerow,column=1,width=entrywidth_plotsettings,variable=scatteringlengthvariable,state='normal')
    unit_scatteringlength = GUI.MakeLabel(frame=settingsframe,row=activerow,column=2,width=unitwidth_plotsettings,text=" a0 ",anchor='e')
    entry_scatteringlength.insert(tk.END,0.0)
    activerow += 1	
    #
    ttk.Separator(settingsframe, orient='horizontal').grid(row=activerow,column = 0,columnspan=3, sticky="ew" )
    activerow += 1	
    #momentumcutoff
    label_momentumcutoff = GUI.MakeLabel(frame=settingsframe,row=activerow,column=0,width=labelwidth_plotsettings,text="momentumcutoff = ",anchor='e')
    momentumcutoffvariable = ''
    entry_momentumcutoff = GUI.MakeEntry(frame=settingsframe,row=activerow,column=1,width=entrywidth_plotsettings,variable=momentumcutoffvariable,state='normal')
    unit_momentumcutoff = GUI.MakeLabel(frame=settingsframe,row=activerow,column=2,width=unitwidth_plotsettings,text='2 k\u2097',anchor='e')
    entry_momentumcutoff.insert(tk.END,30)
    activerow += 1	
    #unitcells 
    label_unitcells = GUI.MakeLabel(frame=settingsframe,row=activerow,column=0,width=labelwidth_plotsettings,text="unitcells = ",anchor='e')
    unitcellsvariable = ''
    entry_unitcells = GUI.MakeEntry(frame=settingsframe,row=activerow,column=1,width=entrywidth_plotsettings,variable=unitcellsvariable,state='normal')
    unit_unitcells = GUI.MakeLabel(frame=settingsframe,row=activerow,column=2,width=unitwidth_plotsettings,text="+1    ",anchor='w')
    entry_unitcells.insert(tk.END,100)
    
    
    #####
    textoutputframe = GUI.MakeLabelFrame(frame=leftframe,row=1,column=0,width=200,height=150,text='Output')
    
    scrollbar = GUI.MakeScrollBar(frame=textoutputframe,column=1,row=0,width=20)
    mylist = tk.Listbox(textoutputframe, yscrollcommand = scrollbar.set,width = 35,height=37,justify='right')
    #for line in range(30):
    #   mylist.insert('end', f"Variable {line} = {line*line}" )
       

    mylist.grid(column = 0,row =0)
    scrollbar.config( command = mylist.yview ) 
    
    
    #######################Output graph
    graphframe = GUI.MakeLabelFrame(frame=main_lattices_frame,row=0,column=1,width=150,height=150,text='')

    ##########
    figure = plt.figure()
    ax2 = plt.subplot2grid((1,3), (0, 2))
    ax = plt.subplot2grid((1,3), (0, 0),colspan=2,sharey = ax2)
    # ax,ax2 = plt.subplots(1, 2, sharey=True)
    figure.tight_layout()
    figure.set_size_inches(8.65,7.3)
    figure.subplots_adjust(wspace=0,left=0.1)
    #ax.plot([],[])
    trapscanvas= FigureCanvasTkAgg(figure, graphframe)  # A tk.DrawingArea.
    trapscanvas.get_tk_widget().grid(column = 0,row =0,padx = 10,pady = 10)

    ax.set_ylim([-100,100])


    toolbar = NavigationToolbar2Tk(trapscanvas,graphframe,pack_toolbar=False)
    toolbar.grid(column = 0,row =1,padx = 10,pady = 10, sticky=tk.EW)


    def compute_latticeparameters():
        nonlocal ax
        nonlocal ax2
        nonlocal trapscanvas
    
        depth = float(entry_depth.get())*1.0e-6
        wavelength = float(entry_wavelength.get())*1.0e-9
        m = float(entry_mass.get())*scipy.constants.physical_constants["atomic mass constant"][0]
        a = float(entry_scatteringlength.get())*scipy.constants.physical_constants["Bohr radius"][0]
        momentumstates = int(entry_momentumcutoff.get())
        unitcells= 2*int(entry_unitcells.get())+1
                
        ftrap = np.sqrt(2.0*sc.k*depth/m)/wavelength
        klattice = 2.0*np.pi/wavelength
        
        maxbands = int(3.0*sc.k*depth/(sc.h*ftrap))
    
        plaser = sc.hbar*2.0*np.pi/wavelength
        Erec  = plaser*plaser/(2.0*m)
        
        u0diver = sc.k*depth/Erec
        
        J = 4.0*Erec*pow(u0diver,0.75)*np.exp(-2.0*np.sqrt(u0diver))/np.sqrt(np.pi)
        Jfreq = J/sc.h

        g = (4.0*np.pi*sc.hbar*sc.hbar * a/m)
		
        LambDicke = np.sqrt(Erec/(sc.hbar*2.0*np.pi*ftrap))
        
       
        #mylist.insert('end','')
        
        ######
        mylist.delete(0,'end')
        mylist.insert('end','depth   =   '+str(u0diver)[:5]+' Erec\n')
        mylist.insert('end','   = '+str(np.around(sc.k*depth/(sc.h*1000.0),decimals = 2))+' kHz \n' )     
        mylist.insert('end','Lamb-Dicke \u03B7        = '+str(np.format_float_scientific(LambDicke,precision = 3)+' \n'))
        mylist.insert('end','\n')
        mylist.insert('end','Tight binding limit\n')
        mylist.insert('end','computed f      = '+str(ftrap/1000.0)[:6]+' kHz \n')
        mylist.insert('end','J = h '+str(np.around(Jfreq,decimals = 6))+' Hz \n')
        mylist.insert('end','tunneling time = '+str(np.format_float_scientific( 1.0/(0.001*Jfreq),precision = 3))+' ms \n')
        mylist.insert('end','m g \u03BB/2 = h '+str(np.around(0.5*m*9.81*wavelength/sc.h,decimals = 3))+' Hz \n')
        
        
     
        #%%%calculating vibrational states
        stepsout = int(entry_unitcells.get())
        dy = 3.0/(4.0*stepsout)
        reducedcoords = np.arange(-stepsout,stepsout+1)*dy
        totalpoints = 2*stepsout+1
        coords = reducedcoords*wavelength
        
        V = -sc.k*depth*np.cos(2.0*np.pi*reducedcoords)*np.cos(2.0*np.pi*reducedcoords)
        
        #%%% Daley thesis
        
        momentumindices = np.arange(-momentumstates,momentumstates+1)
        momentumstatesnum = 2*momentumstates+1
        
        BlochCoefficients = np.zeros((unitcells,maxbands,momentumstatesnum),dtype=np.complex128)
        BlochEnergies     = np.zeros((unitcells,maxbands))
        
        def Hq(q):
            diagonal = Erec*(q + 2.0*momentumindices)**2 +  0.5*depth*sc.k
                        
            H = np.diag(diagonal,k=0)
            H -= np.diag(np.ones(momentumstatesnum-1)*depth*0.25*sc.k,k=-1)
            H -= np.diag(np.ones(momentumstatesnum-1)*depth*0.25*sc.k,k= 1)
            
            return H
        
        quasimomenta = np.linspace(-1.0,1.0,num=unitcells)
        qindex = 0
        for q in quasimomenta:    
            eigenvalues, eigenvectors = LA.eigh(Hq(q))
            
            BlochCoefficients[qindex] = np.transpose(eigenvectors[:,:maxbands])
            BlochEnergies[qindex] = np.absolute(eigenvalues[:maxbands])
                        
            qindex += 1
           
        for band in np.arange(maxbands):
            realcoefficients = np.real(BlochCoefficients[:,band,:])
            if band%2 == 0:
                BlochCoefficients[:,band,:]= 1.0*realcoefficients
            
            if band%2 == 1:
                BlochCoefficients[:,band,:]= 1.0j*realcoefficients
                
        coords2 = coords
        
        kx = klattice*quasimomenta
        planewaves = np.exp(2.0j*klattice*np.outer(momentumindices,coords2))
        qplanewaves = np.exp(1.0j*np.outer(coords2,kx))
        
        
        WannierFunctions = np.zeros((maxbands,unitcells),dtype=np.complex128)
        NormalizedWannierFunctions = np.zeros((maxbands,unitcells),dtype=np.complex128)
        bandcenters = np.mean(BlochEnergies,axis=0)
        bandtop = np.amax(BlochEnergies,axis=0)
        bandbottom = np.amin(BlochEnergies,axis=0)
        bandwidths  = bandtop-bandbottom
        
        for i in np.arange(maxbands):
            WannierFunctions[i] = np.real(np.diag(LA.multi_dot([qplanewaves,BlochCoefficients[:,i,:],planewaves])))

            NormalizedWannierFunctions[i] = WannierFunctions[i]/np.linalg.norm(WannierFunctions[i])
            
            WannierFunctions[i] *= 1.0/np.amax(np.absolute(WannierFunctions[i]))
        
        #########
        ####Computing tunneling rates
        qindex = 0
        tmp=0
        for q in quasimomenta:    
            tmp += BlochEnergies[qindex,0]*np.exp(-1.0j*q*np.pi)/(unitcells)
            qindex += 1
        print(tmp/(sc.h))
        #############
        		
	##Computing interaction energy
        U = g*np.sum(NormalizedWannierFunctions[0]**4)*1.5*wavelength/unitcells
        mylist.insert('end','U = h '+str(np.format_float_scientific( U/sc.h,precision = 3))+' Hz \n')
		        
        #Plotting the potential
        ax.cla()
        ax2.cla()
        
        ax.plot(reducedcoords,V/(sc.h*1000.0), color='black')
        trapscanvas.draw()
               
        mylist.insert('end', '\n' )     
        mylist.insert('end', 'Band centers and widths \n')               
        
        ##%% plotting styles
        plotcolors = color = cm.rainbow(np.linspace(1, 0, maxbands))
        
        ['black','maroon','red','darkorange','yellowgreen','darkgreen','deepskyblue','navy','purple','fuchsia','orchid','deeppink','pink']
        
        
        for i in np.arange(maxbands):
            centreofband = (bandcenters[i]-sc.k*depth)/(sc.h*1000.0)
            top = (bandtop[i]-sc.k*depth)/(sc.h*1000.0)
            bottom = (bandbottom[i]-sc.k*depth)/(sc.h*1000.0)
                
            # ax.plot([-0.75,0.75],[centreofband,centreofband],linestyle='dashed',color=plotcolors[i])
            # ax.plot(reducedcoords,-0.003*ftrap*(sortedeigenvectors[:,bandindex])+i/(sc.h*1000.0),linestyle='solid',color=plotcolors[bandindex])
            ax.plot(reducedcoords,0.0003*ftrap*WannierFunctions[i] + centreofband,linestyle='solid',color=plotcolors[i])
            ax.fill_between([-0.75,0.75],[top,top],[bottom,bottom],color=plotcolors[i], alpha=0.3)
            
			#calculate hopping rate J
            ddx = np.diag(np.ones(unitcells-2),k=-2) #- 8.0*np.diag(np.ones(unitcells-1),k=-1)+8.0*np.diag(np.ones(unitcells-1),k=1)-np.diag(np.ones(unitcells-2),k=2)
            #ddx /= (12.0*)
			
			
            mylist.insert('end',str(i)+':'+str(np.around(bandcenters[i]/(sc.h*1000.0),decimals=2)).rjust(7)+' +- '+str(np.around(bandwidths[i]/(sc.h*1000.0),decimals=2)).rjust(5)+' kHz, J = '+'\n')
                
        
        ax.set_ylim([-1.1*sc.k*depth/(sc.h*1000.0),1.0*sc.k*depth/(sc.h*1000.0)])

        
        #ax.set_xlim([-1.0/3.0,1.0/3.0])
        ax.set_xticks([-0.75,-0.5,-0.25,0,0.25,0.5,0.75])
        ax.set_xticklabels(['-3/4','-2/4','-1/4','0','1/4','2/4','3/4'], minor=False)
        
        for i in np.arange(maxbands):
            ax2.plot(quasimomenta,(BlochEnergies[:,i]-depth*sc.k)/(sc.h*1000.0),color=plotcolors[i],label='n = '+str(i))
        
        ax2.set_xlim([-1.1,1.1])
        ax2.set_ylim([-1.1*sc.k*depth/(sc.h*1000.0),1.0*sc.k*depth/(sc.h*1000.0)])
        ax2.set_xticks([-1,-0.5,0,0.5,1])
        #ax2.set_yticks([])
        ax2.set_xticklabels(['-1','1/2','0','1/2','1'], minor=False)
        ax2.legend(loc='upper right')
               
        ax2.set_title('Band structure',fontsize=15)
        ax2.set_xlabel(r'q ($k_{L}$)',fontsize=15)
        ax2.grid(True)
        ax.set_title('Wannier functions',fontsize=15)
        ax.set_ylabel('Energy (kHz)',fontsize=15)
        ax.set_xlabel('x (\u03bb)',fontsize=15)
        ax.grid(True)
               
        #output.delete(1.0,tk.END)
        #output.insert(tk.END,outpt)
        #output.tag_config('justified', justify='left')
        trapscanvas.draw()
        
        return  
    
    
    
    latt_button = GUI.MakeButton(frame=settingsframe,row=7,column=0,width=20,text='Compute',command= lambda: compute_latticeparameters(),state='normal',relief='raised')
    latt_button.grid(columnspan=3)

    compute_latticeparameters()
