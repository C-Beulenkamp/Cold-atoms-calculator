import tkinter as tk
from tkinter import ttk
import GUIfunctions as GUI
import scipy
import scipy.constants as sc
import numpy as np
import importlib
from traps import set_get_frequencies
from extrafunctions import *
from scipy.optimize import fsolve

abohr = scipy.constants.physical_constants["Bohr radius"][0]

def atomicclouds(tab,row,column):
    #%% Defining frames
    
    atomiccloudsframe = GUI.MakeFrame(frame=tab,column=column,row=row)
    
    # top and bottom frames
    scattering_top_frame = tk.Frame(atomiccloudsframe)
    scattering_top_frame.grid(row=0,column = 0)
    scattering_bottom_frame = tk.Frame(tab)
    scattering_bottom_frame.grid(row=1,column = 0)
    
    ###########top frame parts
    atomchoiceframe = GUI.MakeLabelFrame(frame=scattering_top_frame,row=0,column=0,width=50,height=50,text='Atom choice')
    mass_variable, mass_entry = GUI.FullEntry(atomchoiceframe,row=1,column=0,sizes=[15,10,2],description='m = ',unit=' u')
       
    #
    thermodynamicsframe = GUI.MakeLabelFrame(frame=scattering_top_frame,row=1,column=0,width=50,height=50,text='Thermodynamic quantities')
    N_variable, N_entry = GUI.FullEntry(thermodynamicsframe,row=0,column=0,sizes=[15,10,5],description='N = ',unit=' e6')
    T_variable, T_entry = GUI.FullEntry(thermodynamicsframe,row=1,column=0,sizes=[15,10,5],description='T  = ',unit=' \u03bcK')
    N_entry.insert(tk.END,1.0)
    T_entry.insert(tk.END,1.0)
    
    #
    collisionssizes = [13,14,6]
    collisionsframe = GUI.MakeLabelFrame(frame=scattering_top_frame,row=0,column=1,width=50,height=50,text='Collisional parameters')
    collisionsframe.grid(rowspan=2)
    Rea_variable, Rea_entry = GUI.FullEntry(collisionsframe,row=0,column=0,sizes=collisionssizes,description='Re(a) = ',unit=' a0')
    k2_variable, k2_entry = GUI.FullEntry(collisionsframe,row=1,column=0,sizes=collisionssizes,description='k\u2082(T=0) = ',unit=' cm\u00B3/s')
    K3_variable, K3_entry = GUI.FullEntry(collisionsframe,row=2,column=0,sizes=collisionssizes,description='k\u2083(T=0) = ',unit=' cm\u2076/s')
    ttk.Separator(collisionsframe, orient='horizontal').grid(row=3,column = 0, sticky="ew" )
    amin_variable, amin_entry = GUI.FullEntry(collisionsframe,row=4,column=0,sizes=collisionssizes,description='a- = ',unit=' a0')
    aplus_variable, aplus_entry = GUI.FullEntry(collisionsframe,row=5,column=0,sizes=collisionssizes,description='a+ = ',unit=' a0')
    etamin_variable, etamin_entry = GUI.FullEntry(collisionsframe,row=6,column=0,sizes=collisionssizes,description='eta- = ',unit=' ')
    etaplus_variable, etaplus_entry = GUI.FullEntry(collisionsframe,row=7,column=0,sizes=collisionssizes,description='eta+ = ',unit=' ')
    s0_variable, s0_entry = GUI.FullEntry(collisionsframe,row=8,column=0,sizes=collisionssizes,description='s0 = ',unit=' ')
    Rea_entry.insert(tk.END,0.0)
    k2_entry.insert(tk.END,0.0)
    K3_entry.insert(tk.END,0.0)
    
    
    #
    trapsframe = GUI.MakeLabelFrame(frame=scattering_top_frame,row=0,column=2,width=50,height=50,text='Trap parameters')
    trapsframe.grid(rowspan=2)
    
    trapssizes=[4,8,3]
    GUI.MakeLabel(frame=trapsframe,row=1,column=0,width=15,text="U = \u00bdm(2\u03c0f)\u00b2 x\u00b2",anchor='n')
    fx_variable, fx_entry = GUI.FullEntry(trapsframe,row=2,column=0,sizes=trapssizes,description=' fx = ',unit=' Hz   ')
    fy_variable, fy_entry = GUI.FullEntry(trapsframe,row=3,column=0,sizes=trapssizes,description=' fy = ',unit=' Hz   ')
    fz_variable, fz_entry = GUI.FullEntry(trapsframe,row=4,column=0,sizes=trapssizes,description=' fz = ',unit=' Hz   ')
    fx_entry.insert(tk.END,1.0)
    fy_entry.insert(tk.END,1.0)
    fz_entry.insert(tk.END,1.0)
  
    
  
    trapssizes=[4,8,5]
    GUI.MakeLabel(frame=trapsframe,row=1,column=2,width=15,text="U = m a |x|",anchor='n')
    ax_variable, ax_entry = GUI.FullEntry(trapsframe,row=2,column=2,sizes=trapssizes,description=' ax = ',unit=' m/s\u00b2')
    ay_variable, ay_entry = GUI.FullEntry(trapsframe,row=3,column=2,sizes=trapssizes,description=' ay = ',unit=' m/s\u00b2')
    az_variable, az_entry = GUI.FullEntry(trapsframe,row=4,column=2,sizes=trapssizes,description=' az = ',unit=' m/s\u00b2')
    ax_entry.insert(tk.END,1.0)
    ay_entry.insert(tk.END,1.0)
    az_entry.insert(tk.END,1.0)
    
    def switchlinharm():
        if xlinear.get() == True:
            fx_entry.config(state= "disabled") 
            ax_entry.config(state= "normal") 
        else:
            fx_entry.config(state= "normal") 
            ax_entry.config(state= "disabled") 
        if ylinear.get() == True:
            fy_entry.config(state= "disabled") 
            ay_entry.config(state= "normal") 
        else:
            fy_entry.config(state= "normal") 
            ay_entry.config(state= "disabled") 
        if zlinear.get() == True:
            fz_entry.config(state= "disabled") 
            az_entry.config(state= "normal") 
        else:
            fz_entry.config(state= "normal") 
            az_entry.config(state= "disabled") 
        return
        
    def get_frequencies(): 
        ffx, ffy, ffz, aax, aay, aaz = set_get_frequencies('get')
        fx_entry.delete(0, tk.END)  
        fx_entry.insert(tk.END, np.around(ffx,decimals=2))
        fy_entry.delete(0, tk.END)  
        fy_entry.insert(tk.END, np.around(ffy,decimals=2))
        fz_entry.delete(0, tk.END)  
        fz_entry.insert(tk.END, np.around(ffz,decimals=2))
            
        ax_entry.delete(0, tk.END)  
        ax_entry.insert(tk.END, np.around(aax,decimals=2))
        ay_entry.delete(0, tk.END)  
        ay_entry.insert(tk.END, np.around(aay,decimals=2))
        az_entry.delete(0, tk.END)  
        az_entry.insert(tk.END, np.around(aaz,decimals=2))
        return
        
    getf_button = GUI.MakeButton(frame=trapsframe,row=5,column=0,width=30,text='Get parameters from "Dipole traps" tab',command=lambda: get_frequencies(),state='normal',relief='raised')
    getf_button.grid(columnspan=2)
    
    GUI.MakeLabel(frame=trapsframe,row=6,column=0,width=15,text=" ",anchor='n') #padding on top and bottom
    GUI.MakeLabel(frame=trapsframe,row=0,column=0,width=15,text=" ",anchor='n')
    
    xlinear = tk.IntVar()
    entry_xlinear = GUI.MakeCheckButton(trapsframe,row=2,column=1,text='x linear',command=lambda: switchlinharm(),state='normal',variable=xlinear)
    ylinear = tk.IntVar()
    entry_ylinear = GUI.MakeCheckButton(trapsframe,row=3,column=1,text='y linear',command=lambda: switchlinharm(),state='normal',variable=ylinear)
    zlinear = tk.IntVar()
    entry_zlinear = GUI.MakeCheckButton(trapsframe,row=4,column=1,text='z linear',command=lambda: switchlinharm(),state='normal',variable=zlinear)
    switchlinharm()
    
    #Setting the atomic properties 
    def set_atom(atom):
        atomdata = importlib.import_module("AtomicData."+atom)        
        mass_entry.delete(0, tk.END)  
        amin_entry.delete(0, tk.END)  
        aplus_entry.delete(0, tk.END)
        etamin_entry.delete(0, tk.END)
        etaplus_entry.delete(0, tk.END)
        s0_entry.delete(0, tk.END)
        
        mass_entry.insert(tk.END,np.around(atomdata.m/scipy.constants.physical_constants["atomic mass constant"][0],decimals = 2))
        amin_entry.insert(tk.END, np.around(atomdata.aminus,decimals = 4))
        aplus_entry.insert(tk.END, np.around(atomdata.aplus,decimals = 4))
        etamin_entry.insert(tk.END, np.around(atomdata.etaminus,decimals = 4))
        etaplus_entry.insert(tk.END, np.around(atomdata.etaplus,decimals = 4))
        s0_entry.insert(tk.END, np.around(atomdata.s0,decimals = 4))
        return
        
    atombuttonwidth= 5
    atombuttonframe = tk.Frame(atomchoiceframe)
    atombuttonframe.grid(row=0,column=0)
    Li6Button = GUI.MakeButton(frame=atombuttonframe,row=0,column=0,width=atombuttonwidth,text='\u2076Li',command=lambda: set_atom('6Li'),state='normal',relief='raised')
    Na23Button = GUI.MakeButton(frame=atombuttonframe,row=0,column=1,width=atombuttonwidth,text='\u00b2\u00b3Na',command=lambda: set_atom('23Na'),state='normal',relief='raised')
    K39Button = GUI.MakeButton(frame=atombuttonframe,row=0,column=2,width=atombuttonwidth,text='\u00b3\u2079K',command=lambda: set_atom('39K'),state='normal',relief='raised')
    #K40Button = GUI.MakeButton(frame=buttonframe,row=1,column=1,width=atombuttonwidth,text='40K',command=lambda: set_atom('40K'),state='normal',relief='raised')
    #Rb85Button = GUI.MakeButton(frame=buttonframe,row=0,column=2,width=atombuttonwidth,text='85Rb',command=lambda: set_atom('85Rb'),state='normal',relief='raised')
    Rb87Button = GUI.MakeButton(frame=atombuttonframe,row=1,column=0,width=atombuttonwidth,text='\u2078\u2077Rb',command=lambda: set_atom('87Rb'),state='normal',relief='raised')
    CsButton = GUI.MakeButton(frame=atombuttonframe,row=1,column=1,width=atombuttonwidth,text='\u00B9\u00B3\u00B3Cs',command=lambda: set_atom('133Cs'),state='normal',relief='raised')
    
    ###########bottom frame parts
    datasetinfo0 = GUI.MakeText(scattering_top_frame,row=0,column=3,width=40,height=7)
    datasetinfo0.grid(rowspan=2)
    datasetinfo1 = GUI.MakeText(scattering_bottom_frame,row=0,column=0,width=45,height=27)
    datasetinfo2 = GUI.MakeText(scattering_bottom_frame,row=0,column=1,width=45,height=27)
    datasetinfo3 = GUI.MakeText(scattering_bottom_frame,row=0,column=2,width=45,height=27)
    
    
    #Desired_font = tk.font.Font( size = 15) #,weight = "bold"
    #datasetinfo0.configure(font = Desired_font)
    #datasetinfo1.configure(font = Desired_font)
    #datasetinfo2.configure(font = Desired_font)
    #datasetinfo3.configure(font = Desired_font)
   
    
    def compute_K3(): 
        K3_entry.delete(0, tk.END)
        
        #Losses  
        m = float(mass_entry.get())*scipy.constants.physical_constants["atomic mass constant"][0]
        a = float(Rea_entry.get())
        aplus    = float(aplus_entry.get())
        amin     = float(amin_entry.get())
        etaplus  = float(etaplus_entry.get())    
        etaminus = float(etamin_entry.get())    
        s0       = float(s0_entry.get())    
        
        #Formulas from Beringer thesis
        if a > 0:     
            Cthreebody = 67.1*np.exp(-2.0*etaplus)*(np.cos(s0 * np.log(a/aplus))**2 + np.sinh(etaplus)**2) +16.8*(1.0-np.exp(-4.0*etaplus))
        if a < 0:
            Cthreebody = 4590.0*np.sinh(2.0*etaminus)/( np.sin(s0*np.log(a/amin))**2 + np.sinh(etaminus)**2)
        if a == 0:
            Cthreebody = 0.0
        LL3 = 3.0 * Cthreebody*sc.hbar*np.power(a*abohr,4)/m
        
        K3_entry.insert(tk.END, np.format_float_scientific(LL3*1.0e12,precision = 3))
        return
    
    L3_button =  GUI.MakeButton(frame=collisionsframe,row=9,column=0,width=20,text='Compute universal L3',command=lambda: compute_K3(),state='normal',relief='raised')
    
       
    # #%% computing stuff
    def compute_atomicclouds():       
        a = float(Rea_entry.get())
        k2 = float(k2_entry.get())*1.0e-6
        L3 = float(K3_entry.get())*1.0e-12
        N = float(N_entry.get())*1.0e6
        T = float(T_entry.get())*1.0e-6
        m = float(mass_entry.get())*scipy.constants.physical_constants["atomic mass constant"][0]
        fx = float(fx_entry.get())
        fy = float(fy_entry.get())
        fz = float(fz_entry.get())         
        ax = float(ax_entry.get())
        ay = float(ay_entry.get())
        az = float(az_entry.get())         
       
        #Reduced mass for use in scattering calculations
        redmass = 0.5*m
        
        #Imaginary part of the scattering length a -i b
        b = redmass*k2/(4.0*np.pi*sc.hbar)
        
        
        omegabar = 2.0*np.pi*np.power(fx*fy*fz,1.0/3.0)  #Mean trap frequencies
        abar = np.sqrt(sc.hbar/(m*omegabar)) 
        
        
        
        #%%BEC stuff
        Nc = np.power(scipy.constants.k*T/(0.94*scipy.constants.h),3.0)/(fx*fy*fz)
        Tc = 0.94*scipy.constants.h*np.power(N*fx*fy*fz,1.0/3.0)/scipy.constants.k    
        
        finiteN_Tc_correction = -0.73*(2.0*np.pi/3.0)*(fx + fy + fz)/( omegabar * np.power(N,1.0/3.0) ) ##Pethick, Smith
        
        Ncondensate = N*(1.0- np.power(T/Tc,3))
        if T > Tc:
            Ncondensate = 0.0
        Nth = N - Ncondensate
        
        interactiong = 4.0*np.pi*sc.hbar*sc.hbar*a*abohr/m
        chemicalpotential = 0.5*np.power(15*Ncondensate * a*abohr/ abar,2.0/5.0) * sc.hbar* omegabar
        N0condensate = chemicalpotential/interactiong
        BECfraction = Ncondensate/N
        
        #Harmonic oscillator lengths 
        lhox = np.sqrt(sc.hbar/(m*2.0*np.pi*fx))
        lhoy = np.sqrt(sc.hbar/(m*2.0*np.pi*fy))
        lhoz = np.sqrt(sc.hbar/(m*2.0*np.pi*fz))
                
        # 
        U0 = sc.h*sc.h*a*abohr/(m*np.pi)
        
        #Thomas-Fermi lengths
        Rx = np.sqrt(2.0*chemicalpotential/(m*4*np.pi*np.pi*fx*fx))
        Ry = np.sqrt(2.0*chemicalpotential/(m*4*np.pi*np.pi*fy*fy))
        Rz = np.sqrt(2.0*chemicalpotential/(m*4*np.pi*np.pi*fz*fz))
        Rbar = np.power(Rx*Ry*Rz,1.0/3.0)

        EkinoverEint = sc.hbar*sc.hbar/(2.0*m*chemicalpotential*Rbar*Rbar)

        sigmafactor = np.sqrt(2.0*scipy.constants.k*T/m)/(2.0*np.pi)
        sigmax = sigmafactor/fx
        sigmay = sigmafactor/fy
        sigmaz = sigmafactor/fz
        
        if xlinear.get() ==True: 
            sigmax = sc.k * T/(m* ax)
        if ylinear.get() ==True: 
            sigmay = sc.k * T/(m* ay)
        if zlinear.get() ==True: 
            sigmaz = sc.k * T/(m* az)
        
        lambT = np.sqrt(2.0*np.pi*scipy.constants.hbar*scipy.constants.hbar/(m*scipy.constants.Boltzmann*T))
        
        n0 = Nth
        if xlinear.get() ==True: 
            n0 /= (2.0*sc.k*T)/(m*ax)
        else:
            n0 *= np.sqrt(m / (2.0*np.pi*sc.k * T)) * 2.0* np.pi*fx
        if ylinear.get() ==True: 
            n0 /= (2.0*sc.k*T)/(m*ay)
        else:
            n0 *= np.sqrt(m / (2.0*np.pi*sc.k * T)) * 2.0* np.pi*fy
        if zlinear.get() ==True: 
            n0 /= (2.0*sc.k*T)/(m*az)
        else:
            n0 *= np.sqrt(m / (2.0*np.pi*sc.k * T)) * 2.0* np.pi*fz
               
        
        if (xlinear.get() ==False) and (ylinear.get() ==False) and (zlinear.get() ==False):
            fugacitytarget = Nth*np.power(sc.hbar * omegabar/(sc.k*T),3.0)
            def tempfunc(z):
                return gPolyLogarithm(z,3.0) - fugacitytarget
            
            fugacity = fsolve(tempfunc, 1.0)
            bosePSD = gPolyLogarithm(fugacity[0],1.5)
            bosen0 = bosePSD/(lambT*lambT*lambT)
        else:
            bosen0 = 0.0
            bosePSD = 0.0
        
        PSD = n0*lambT*lambT*lambT
        
        #Kraemer thesis
        L3unitaritylimit = 844.0*(sc.hbar/m)*(sc.hbar/m)*(sc.hbar/m)*(sc.hbar/(sc.k*T))*(sc.hbar/(sc.k*T))
            
        condensatethreebodylossrate = N0condensate*N0condensate*L3/6.0
        if condensatethreebodylossrate > 0.0:
            condensatethreebodyliftime  = 1.0/condensatethreebodylossrate
        else:
            condensatethreebodyliftime = float('inf')
        
        threebodylossrate = n0*n0*L3
        if threebodylossrate > 0.0:
            threebodyliftime  = 1.0/threebodylossrate
        else:
            threebodyliftime = float('inf')
        
        vthermal = np.sqrt(3.0*sc.k*T/m)
        
        kthermal = m*vthermal/sc.hbar
        
        vbar = np.sqrt(8.0*sc.k*T/(np.pi*m))
        
        k_collisionwavenumber = m*vthermal/sc.hbar#np.sqrt(np.pi*sc.hbar*sc.hbar/(4.0*sc.k*T*m))
        cross_sectionzerotemp = 8.0*np.pi*a*a*abohr*abohr
        cross_section = 8.0*np.pi*a*a*abohr*abohr/(1.0+k_collisionwavenumber*k_collisionwavenumber*a*abohr*abohr*a)
        
        #%%% Averaging elastic and inelastic collision rates over momentum distribution.     
		#Energy dependent scattering rate.        
        beta = 1.0/(T*sc.k)
        kcutoff = np.sqrt(m*20.0 * sc.k*T)/sc.hbar
        dk = 0.001*kcutoff
        integralks = np.arange(1000)*dk
        scatteringrate_integrand = 8.0*a*abohr*a*abohr*np.pi/(1.0 + a*abohr*a*abohr*integralks*integralks+2.0*integralks*b)  # Elastic cross section 
        scatteringrate_integrand *= sc.hbar*integralks/redmass #Multiply by velocity
        scatteringrate_integrand *= 4.0*np.pi*integralks*integralks #Multiply by Jacobian
        scatteringrate_integrand *= np.exp(-0.5*beta*(sc.hbar*integralks)*(sc.hbar*integralks)/(redmass)) #Multiply by Boltzmann factor
        scatteringrate_integrand *= sc.hbar*sc.hbar*sc.hbar*np.power(beta/(np.pi*redmass),1.5) #Normalization factor
        scatteringrate_integrand *= 0.25 # fudge factor to make sure it gets to the right zero temperature result. should look up where this factor belongs.
        
        peakscatteringrate_energydependent = n0*np.sum(scatteringrate_integrand*dk)
        
        inelastic_integrand = 1.0/(1.0 + a*abohr*a*abohr*integralks*integralks+2.0*integralks*b)  # sigma_inel(k)*k = k2(k) = k2(k=0)/(1+a^2 k^2)
        inelastic_integrand *= 4.0*np.pi*integralks*integralks #Multiply by Jacobian
        inelastic_integrand *= np.exp(-0.5*beta*(sc.hbar*integralks)*(sc.hbar*integralks)/(redmass)) #Multiply by Boltzmann factor
        inelastic_integrand *= sc.hbar*sc.hbar*sc.hbar*np.power(beta/(np.pi*redmass),1.5) #Normalization factor
        inelastic_integrand *= 0.25 # fudge factor to make sure it gets to the right zero temperature result. should look up where this factor belongs.
        
        
        k2thermalaveragereductionfactor = np.sum(inelastic_integrand*dk)
        inelasticscatteringrate_energydependent = n0*k2thermalaveragereductionfactor*k2
        
        		
        meanfreepath = 1.0/(n0*cross_section)
        
        #%%% Printing output
        datasetinfo0.delete(1.0,tk.END)
        datasetinfo1.delete(1.0,tk.END)
        datasetinfo2.delete(1.0,tk.END)
        datasetinfo3.delete(1.0,tk.END)
        
        #trap 
        
        outpt = '-----------Trap------------\n' 
        outpt += 'mean f  = '+str(np.around(np.power(fx*fy*fz,1.0/3.0),decimals=1))+' Hz \n'
        outpt += 'l_ho_x = '+str(np.around(lhox*1.0e6,decimals=1))+' \u03bcm\n'
        outpt += 'l_ho_y = '+str(np.around(lhoy*1.0e6,decimals=1))+' \u03bcm\n'
        outpt += 'l_ho_z = '+str(np.around(lhoz*1.0e6,decimals=1))+' \u03bcm\n'
        
        datasetinfo0.insert(tk.END,outpt, 'justified')
        
        #Collisional part
        
        outpt = '-----------Collisions------\n'     
        outpt += 'v_th     = '+str(np.around(vthermal*100.0,decimals=2))+' cm/s \n'
        outpt += 'k_th     = '+str('{:.2e}'.format(kthermal*0.01))[:8]+' 1/cm = '+str('{:.2e}'.format(kthermal*abohr))[:8]+' 1/a0 \n \n'
        outpt += '--Two-body---\n'  
        outpt += '       a  = ('+str(np.around(a,decimals=1))+' - i '+str('{:.2e}'.format(b/abohr))[:8]+') a0 \n'
        outpt += 'k_th * a  = '+str(np.around(kthermal*abohr*a,decimals=1))+'  \n '
        outpt += '\n'
        outpt += 'cross section  \u03c3(v=0) = '+str('{:.2e}'.format(cross_sectionzerotemp*1.0e4))[:8]+' cm\u00b2 \n'
        outpt += 'unitarity limit: \n           \u03c3(v_th,a->inf) = '+str('{:.2e}'.format(8.0*np.pi*1.0e4/(kthermal*kthermal)))[:8]+' cm\u00b2 \n'
        outpt += '\n'
        outpt += 'Finite T collision factor = '+str('{:.2e}'.format(k2thermalaveragereductionfactor))[:8]+'\n'
        outpt += '                 <k2(k)>  = '+str('{:.2e}'.format(k2thermalaveragereductionfactor*k2*1.0e6))[:8]+' cm\u00B3/s \n \n'
         
        outpt += '--Three-body--\n'  
        outpt += ' L3 (T=0)            = '+str('{:.2e}'.format(L3*1.0e12))[:8]+' cm\u2076/s \n'
        outpt += ' L3 unitarity limit  = '+str('{:.2e}'.format(L3unitaritylimit*1.0e12))[:8]+' cm\u2076/s \n'
        
        
        
        
        
        
        datasetinfo1.insert(tk.END,outpt, 'justified')
        
        #Thermal cloud part
        
        outpt = '-----------Thermal---------\n'       
        outpt += 'Nth      = '+str('{:.2e}'.format(Nth))[:8]+' \n'
        outpt += '\u03BB_T     = '+str(np.around(lambT*1.0e9,decimals=2))+' nm\n'
        outpt += 'n0 (class)= '+str('{:.2e}'.format(n0*1.0e-6))[:8]+' 1/cm\u00B3\n'
        outpt += 'n0  (bose)= '+str('{:.2e}'.format(bosen0*1.0e-6))[:8]+' 1/cm\u00B3\n'
        outpt += 'PSD(class)= '+str('{:.2e}'.format(PSD))[:8]+' \n'
        outpt += 'PSD (bose)= '+str('{:.2e}'.format(bosePSD))[:8]+' \n'
        outpt += '\n'
        
        outpt += 'cross section (v_th) = '+str('{:.2e}'.format(cross_section*1.0e4))[:8]+' cm\u00b2 \n'
        outpt += 'peak scattering rate: \n           n0 <\u03c3(v) v> = '+str(np.around(peakscatteringrate_energydependent,decimals=1))+' Hz \n'
        outpt += '\n'
        outpt += 'peak 2-body loss rate = '+str(np.around(inelasticscatteringrate_energydependent,decimals=1))+' Hz  \n'
        outpt += 'peak 2-body lifetime   = '+str(np.around(1.0/inelasticscatteringrate_energydependent,decimals=1))+' s  \n'

        outpt += '\n'
        #outpt += 'C  = '+str('{:.2e}'.format(Cthreebody))[:8]+'  \n'
        outpt += 'peak 3-body loss rate = '+str(np.around(threebodylossrate,decimals=1))+' Hz  \n'
        outpt += '  avg 3-body loss rate = '+str(np.around(threebodylossrate/np.sqrt(27.0),decimals=1))+' Hz  \n'
        outpt += '        3-body lifetime    = '+str(np.around(threebodyliftime*np.sqrt(27.0),decimals=3))+' s \n'
        outpt += '\n'
        outpt += '1/e cloud radii \n'
        
        if xlinear.get() ==True: 
            outpt += 'r_x = '+str(np.around(sigmax*1.0e6,decimals=1))+' \u03bcm :     exp(-x/rx)\n'
        else:
            outpt += 'r_x = '+str(np.around(sigmax*1.0e6,decimals=1))+' \u03bcm : exp(-x\u00b2/rx\u00b2)\n'
        
        if ylinear.get() ==True: 
            outpt += 'r_y = '+str(np.around(sigmay*1.0e6,decimals=1))+' \u03bcm :     exp(-y/ry)\n'
        else:
            outpt += 'r_y = '+str(np.around(sigmay*1.0e6,decimals=1))+' \u03bcm : exp(-y\u00b2/ry\u00b2)\n'
        if zlinear.get() ==True: 
            outpt += 'r_z = '+str(np.around(sigmaz*1.0e6,decimals=1))+' \u03bcm :     exp(-z/rz)\n'
        else:
            outpt += 'r_z = '+str(np.around(sigmaz*1.0e6,decimals=1))+' \u03bcm : exp(-z\u00b2/rz\u00b2)\n'
        
        outpt += 'mean free path (peak) = '+str(np.around(meanfreepath*1.0e6,decimals=1))+' \u03bcm \n'
        
        datasetinfo2.insert(tk.END,outpt, 'justified')
        #BEC cloud part
        
        outpt = '------------BEC------------\n'        
        outpt += 'Tc = '+str('{:.2e}'.format(Tc*1.0e6))[:8]+' \u03bcK\n'
        outpt += 'Nc = '+str('{:.2e}'.format(Nc))[:8]+' \n' 
        
        outpt += 'Finite N corrected Tc = '+str('{:.2e}'.format(Tc*(1.0+finiteN_Tc_correction)*1.0e6))[:8]+' \u03bcK\n'
        
        outpt += '\n'
        outpt += 'BEC frac = '+str(int(np.trunc(BECfraction*100.0))/100)+' \n'
        outpt += 'N_BEC    = '+str('{:.2e}'.format(Ncondensate))[:8]+' \n'
        outpt += '\n'
        outpt += 'N a/ l_ho = '+str(np.around(Ncondensate*a*abohr/np.power(lhox*lhoy*lhoz,1.0/3.0),decimals=1))+' \n'
        outpt += 'TF approximation:  \n'
        outpt += '\u03bc     = '+str('{:.2e}'.format(chemicalpotential/(2.0*np.pi*sc.hbar)))[:8]+' Hz ='+str('{:.2e}'.format(chemicalpotential*1.0e6/sc.k))[:8]+' \u03bcK\n'
        outpt += 'R_x = '+str(np.around(Rx*1.0e6,decimals=1))+' \u03bcm\n'
        outpt += 'R_y = '+str(np.around(Ry*1.0e6,decimals=1))+' \u03bcm\n'
        outpt += 'R_z = '+str(np.around(Rz*1.0e6,decimals=1))+' \u03bcm\n'
        outpt += 'Ekin/Eint = '+str('{:.2e}'.format(EkinoverEint))[:8]+' \n'
        outpt += '\n'
        outpt += 'n_0_BEC = '+str('{:.2e}'.format(N0condensate*1.0e-6))[:8]+' 1/cm\u00B3\n'
        outpt += 'BEC 2-body loss rate = '+str('{:.2e}'.format(0.0))[:8]+' 1/s  \n'
        outpt += 'BEC 2-body lifetime   = '+str('{:.2e}'.format(0.0))[:8]+' s \n'
        outpt += 'BEC 3-body loss rate = '+str('{:.2e}'.format(condensatethreebodylossrate))[:8]+' 1/s  \n'
        outpt += 'BEC 3-body lifetime   = '+str('{:.2e}'.format(condensatethreebodyliftime))[:8]+' s \n'
        
        datasetinfo3.insert(tk.END,outpt, 'justified')
        
        
        
        
    compute_button = GUI.MakeButton(frame=trapsframe,row=5,column=2,width=10,text='Compute',command=lambda: compute_atomicclouds(),state='normal',relief='raised')