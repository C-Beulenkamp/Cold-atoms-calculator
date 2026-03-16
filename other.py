import tkinter as tk                    
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
import GUIfunctions as GUI
import scipy.constants as sc
from sympy import S
from extrafunctions import QuietWigner6j


def GaussianBeamCalculator(tab,row,column):
        
        
        
        GaussianbeamFrame = GUI.MakeLabelFrame(tab,row=row,column=column,width=100,height=100,text='Gaussian beams')
        
        sizes = [12,15,6]
        
        wavelength, wavelengthEntry = GUI.FullEntry(GaussianbeamFrame,row=2,column=0,sizes=sizes,description='\u03bb         = ',unit=' nm')
        power, powerEntry = GUI.FullEntry(GaussianbeamFrame,row=3,column=0,sizes=sizes,description='P         = ',unit=' mW')
        wx, wxEntry = GUI.FullEntry(GaussianbeamFrame,row=4,column=0,sizes=sizes,description='wx         = ',unit=' \u03bcm')
        wy, wyEntry = GUI.FullEntry(GaussianbeamFrame,row=5,column=0,sizes=sizes,description='wy         = ',unit=' \u03bcm')
        
        
        Gaussian_output = GUI.MakeText(GaussianbeamFrame,row=6,column=0,width=30,height=10)
        Gaussian_output.tag_config('justified', justify='left')
        
        def compute_beam():
            lamb = float(wavelengthEntry.get()) * 1.0e-9
            P = float(powerEntry.get()) *1.0e-3
            wx = float(wxEntry.get()) * 1.0e-6
            wy = float(wyEntry.get()) * 1.0e-6
            
            rayleighlengthx = np.pi*wx*wx/lamb
            rayleighlengthy = np.pi*wy*wy/lamb 
            intensity = 2.0*P/(np.pi*wx*wy)
            E0 = np.sqrt(2.0*intensity/(sc.c*sc.epsilon_0))
            
            divergenceX = lamb/(np.pi*wx)
            divergenceY = lamb/(np.pi*wy)
            NAx = wx/rayleighlengthx
            NAy = wy/rayleighlengthy
            
            Gaussian_output.delete(1.0,tk.END)
            gaussianinitoutput  = 'zR_x       = '+str('{:.2e}'.format(rayleighlengthx*1.0e3))[:8]+' mm \n'
            gaussianinitoutput += 'zR_y       = '+str('{:.2e}'.format(rayleighlengthy*1.0e3))[:8]+' mm \n'
            gaussianinitoutput += 'full div X = '+str(np.around(2.0*divergenceX*180.0/np.pi,decimals=2))+' deg \n'
            gaussianinitoutput += 'full div Y = '+str(np.around(2.0*divergenceY*180.0/np.pi,decimals=2))+' deg \n'
            gaussianinitoutput += 'NA X       = '+str(np.around(NAx,decimals=5))+'  \n'     
            gaussianinitoutput += 'NA Y       = '+str(np.around(NAy,decimals=5))+'  \n'            
            gaussianinitoutput += 'I0         = '+str('{:.2e}'.format(intensity*0.1))[:8]+' mW/cm\u00b2\n'
            gaussianinitoutput += 'E0         = '+str('{:.2e}'.format(E0))[:8]+' V/m'
            Gaussian_output.insert(tk.END,gaussianinitoutput, 'justified')
            return
        
        Gaussian_button = GUI.MakeButton(frame=GaussianbeamFrame,row=7,column=0,width=12,text='Compute',command=lambda: compute_beam(),state='normal',relief='raised')
        
        
        # add_label(GaussianbeamFrame, '\u03bb         = ', 2, 0,1,1)
        # Gaussian_entry_lambda = add_entry(GaussianbeamFrame,2,1,1,1,entrywidth,tk.E)
        # add_label(GaussianbeamFrame, ' nm', 2, 2,1,1)
        
        # add_label(GaussianbeamFrame, 'P         = ', 3, 0,1,1)
        # Gaussian_entry_P = add_entry(GaussianbeamFrame,3,1,1,1,entrywidth,tk.E)
        # add_label(GaussianbeamFrame, ' mW', 3, 2,1,1)

        # add_label(GaussianbeamFrame, 'wx         = ', 4, 0,1,1)
        # Gaussian_entry_wx = add_entry(GaussianbeamFrame,4,1,1,1,entrywidth,tk.E)
        # add_label(GaussianbeamFrame, ' \u03bcm', 4, 2,1,1)
        
        # add_label(GaussianbeamFrame, 'wy         = ', 5, 0,1,1)
        # Gaussian_entry_wy = add_entry(GaussianbeamFrame,5,1,1,1,entrywidth,tk.E)
        # add_label(GaussianbeamFrame, ' \u03bcm', 5, 2,1,1)
        
        # gaussianinitoutput = 'zR_x =  mm \n   zR_y =  mm \n I0 = mW/cm\u00b2'       
        # Gaussian_output = create_outputwindow(GaussianbeamFrame,gaussianinitoutput,7, 0,1,3,30,10)
    
    

        

        
        # Gaussian_button = add_button(GaussianbeamFrame,'Compute',compute_beam,6,0,1,3,10)
        
        


def FiberCouplingCalculator(tab,row,column):
        
             
        FiberCouplingFrame = GUI.MakeLabelFrame(tab,row=row,column=column,width=100,height=100,text='Fiber coupling')
        
        sizes = [12,15,6]
        
        wavelength, wavelengthEntry = GUI.FullEntry(FiberCouplingFrame,row=0,column=0,sizes=sizes,description='\u03bb         = ',unit=' nm')
        MFD, MFDEntry = GUI.FullEntry(FiberCouplingFrame,row=1,column=0,sizes=sizes,description='MDF        = ',unit=' \u03bcm')
        beamdiameter, beamdiameterEntry = GUI.FullEntry(FiberCouplingFrame,row=2,column=0,sizes=sizes,description='D         = ',unit=' \u03bcm')
        
        FiberCoupling_output = GUI.MakeText(FiberCouplingFrame,row=6,column=0,width=30,height=3)
        FiberCoupling_output.tag_config('justified', justify='left')
       
        def compute_focallength():
            lamb = float(wavelengthEntry.get()) * 1.0e-9
            MDF = float(MFDEntry.get()) * 1.0e-6
            D = float(beamdiameterEntry.get()) * 1.0e-6
            
            f = D * np.pi * MDF /(4.0*lamb)
            
            FiberCoupling_output.delete(1.0,tk.END)
            gaussianinitoutput  = 'f         = '+str('{:.2e}'.format(f*1.0e3))[:8]+' mm \n'
            gaussianinitoutput += 'approx NA = '+str(np.sin(lamb/(np.pi*0.5*MDF)))[:8]+' \n'
            gaussianinitoutput += ' \n'
            FiberCoupling_output.insert(tk.END,gaussianinitoutput, 'justified')
            return
        
        FiberCoupling_button = GUI.MakeButton(frame=FiberCouplingFrame,row=7,column=0,width=12,text='Compute',command=lambda: compute_focallength(),state='normal',relief='raised')

def DiffractionGratingCalculator(tab,row,column):
        
             
        DiffractionGratingFrame = GUI.MakeLabelFrame(tab,row=row,column=column,width=100,height=100,text='Diffraction Grating')
        
        sizes = [15,15,4]
        
        wavelength, wavelengthEntry = GUI.FullEntry(DiffractionGratingFrame,row=0,column=0,sizes=sizes,description='\u03bb   = ',unit=' nm')
        lines, linesEntry           = GUI.FullEntry(DiffractionGratingFrame,row=1,column=0,sizes=sizes,description='lines/mm = ',unit=' ')
        AngleIn, AngleInEntry       = GUI.FullEntry(DiffractionGratingFrame,row=2,column=0,sizes=sizes,description='angle in = ',unit=' deg')
        
        AngleBlaze, AngleBlazeEntry       = GUI.FullEntry(DiffractionGratingFrame,row=4,column=0,sizes=sizes,description='angle blaze = ',unit=' deg')
        AngleBlazeEntry.config(state= "disabled") 
        ReflectiveGrating = tk.IntVar()
        def SwitchReflective():
            if ReflectiveGrating.get() == True:
                AngleBlazeEntry.config(state= "normal") 
            else:                
                AngleBlazeEntry.config(state= "disabled") 
            return
        entry_ReflectiveGrating = GUI.MakeCheckButton(frame=DiffractionGratingFrame,row=3,column=0,text='Reflective',command=lambda: SwitchReflective(),state='active',variable=ReflectiveGrating)
        
        
        #FiberCoupling_output = GUI.MakeText(FiberCouplingFrame,row=6,column=0,width=30,height=3)
        #FiberCoupling_output.tag_config('justified', justify='left')
       
        def compute_focallength():
            lamb = float(wavelengthEntry.get()) * 1.0e-9
            MDF = float(MFDEntry.get()) * 1.0e-6
            D = float(beamdiameterEntry.get()) * 1.0e-6
            
            f = D * np.pi * MDF /(4.0*lamb)
            
            FiberCoupling_output.delete(1.0,tk.END)
            gaussianinitoutput  = 'f       = '+str('{:.2e}'.format(f*1.0e3))[:8]+' mm \n'
            gaussianinitoutput += ' \n'
            FiberCoupling_output.insert(tk.END,gaussianinitoutput, 'justified')
            return
        
        DiffractionGrating_button = GUI.MakeButton(frame=DiffractionGratingFrame,row=7,column=0,width=12,text='Compute',command=lambda: compute_focallength(),state='normal',relief='raised')
        
        

def CavityCalculator(tab,row,column):
        
             
        CavityFrame = GUI.MakeLabelFrame(tab,row=row,column=column,width=100,height=100,text='Cavities')
        
        sizes = [15,15,4]
        wavelength, wavelengthEntry = GUI.FullEntry(CavityFrame,row=0,column=0,sizes=sizes,description='\u03bb   = ',unit=' nm')
        Length, LengthEntry           = GUI.FullEntry(CavityFrame,row=1,column=0,sizes=sizes,description='Cavity Length = ',unit=' mm')
        Radius1, Radius1Entry       = GUI.FullEntry(CavityFrame,row=2,column=0,sizes=sizes,description='Radius curv 1 = ',unit=' mm')
        Reflectivity1, Reflectivity1Entry       = GUI.FullEntry(CavityFrame,row=3,column=0,sizes=sizes,description='reflectivity 1 ',unit=' ')
        Radius2, Radius2Entry       = GUI.FullEntry(CavityFrame,row=4,column=0,sizes=sizes,description='Radius curv 2 = ',unit=' mm')
        Reflectivity2, Reflectivity2Entry       = GUI.FullEntry(CavityFrame,row=5,column=0,sizes=sizes,description='reflectivity 2 ',unit=' ')
        # AngleIn, AngleInEntry       = GUI.FullEntry(DiffractionGratingFrame,row=2,column=0,sizes=sizes,description='angle in = ',unit=' deg')
        
        # AngleBlaze, AngleBlazeEntry       = GUI.FullEntry(DiffractionGratingFrame,row=4,column=0,sizes=sizes,description='angle blaze = ',unit=' deg')
        # AngleBlazeEntry.config(state= "disabled") 
        # ReflectiveGrating = tk.IntVar()
        # def SwitchReflective():
            # if ReflectiveGrating.get() == True:
                # AngleBlazeEntry.config(state= "normal") 
            # else:                
                # AngleBlazeEntry.config(state= "disabled") 
            # return
        # entry_ReflectiveGrating = GUI.MakeCheckButton(frame=DiffractionGratingFrame,row=3,column=0,text='Reflective',command=lambda: SwitchReflective(),state='active',variable=ReflectiveGrating)
        
        
        
        
        #StabilityParameter
        #Finesse
        #FSR
        
        # FiberCoupling_output = GUI.MakeText(FiberCouplingFrame,row=6,column=0,width=30,height=3)
        # FiberCoupling_output.tag_config('justified', justify='left')
       
        def compute_cavity():
            # lamb = float(wavelengthEntry.get()) * 1.0e-9
            # MDF = float(MFDEntry.get()) * 1.0e-6
            # D = float(beamdiameterEntry.get()) * 1.0e-6
            
            # f = D * np.pi * MDF /(4.0*lamb)
            
            # FiberCoupling_output.delete(1.0,tk.END)
            # gaussianinitoutput  = 'f       = '+str('{:.2e}'.format(f*1.0e3))[:8]+' mm \n'
            # gaussianinitoutput += ' \n'
            # FiberCoupling_output.insert(tk.END,gaussianinitoutput, 'justified')
            return
        
        Cavity_button = GUI.MakeButton(frame=CavityFrame,row=7,column=0,width=12,text='Compute',command=lambda: compute_cavity(),state='normal',relief='raised')
