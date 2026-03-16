import tkinter as tk
from tkinter import ttk
import glob, os
import numpy as np
import GUIfunctions as GUI
from sympy import S
from numpy import linalg as LA
from functools import partial
import scipy.constants as sc
import importlib
from AngularMomentumOperatorClass import Operator,AngMomOperators
import matplotlib.pyplot as plt
import time
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
    



def feshbach_spectrum(tab,row,column):
    os.chdir("AtomicData")
    #The frame for the whole tab
    main_collisions_frame = GUI.MakeFrame(frame=tab,column=column,row=row)
    
    #The frame on the left
    choiceframe = GUI.MakeFrame(frame=main_collisions_frame,row=0,column=0,sticky='ns')
    atoms_frame = GUI.MakeLabelFrame(frame=choiceframe,row=0,column=0,width=50,height=50,text='Select atom')
    dataset_frame = GUI.MakeLabelFrame(frame=choiceframe,row=1,column=0,width=50,height=50,text='Select data')
    dataset_frame.destroy()
    
    #The frame on the right
    graphframe = GUI.MakeFrame(frame=main_collisions_frame,row=0,column=1,sticky='ns')
    figure = plt.Figure(figsize=(9.25,5.975), dpi=100)
    figure.set_tight_layout(True)
    ax = figure.add_subplot(111)
    ax.set_title('Cross-sections')
    canvas= FigureCanvasTkAgg(figure, graphframe)
    canvas.get_tk_widget().grid(column = 0,row =0,padx = 10,pady = 10)
    #
    ax.plot([-250,250],[0,0])
    ax.set_xlim([-250,250])
    ax.set_ylim([-1,1])
    ax.set_ylabel(" ")
    ax.set_xlabel("B (G)")
    figure.tight_layout()
    # 
    settingsframe = GUI.MakeFrame(frame=graphframe,row=1,column=0,sticky='ns')
    leftsettings = GUI.MakeLabelFrame(frame=settingsframe,row=0,column=0,width=50,height=50,text='Plot settings')
    leftsettings.grid(rowspan=2)
    
    entrysizes = [10,15,3]
    Bmin_variable, Bmin_entry = GUI.FullEntry(leftsettings,row=0,column=0,sizes=entrysizes,description='Bmin',unit=' G')
    Bmax_variable, Bmax_entry = GUI.FullEntry(leftsettings,row=1,column=0,sizes=entrysizes,description='Bmax',unit=' G')
    ymin_variable, ymin_entry = GUI.FullEntry(leftsettings,row=2,column=0,sizes=entrysizes,description='ymin',unit=' a0')
    ymax_variable, ymax_entry = GUI.FullEntry(leftsettings,row=3,column=0,sizes=entrysizes,description='ymax',unit=' a0')
    logymin_variable, logymin_entry = GUI.FullEntry(leftsettings,row=4,column=0,sizes=entrysizes,description='logymin',unit=' a0')
    logymax_variable, logymax_entry = GUI.FullEntry(leftsettings,row=5,column=0,sizes=entrysizes,description='logymax',unit=' a0')
    Bmin_entry.insert(tk.END,0.0)
    Bmax_entry.insert(tk.END,200.0)
    ymin_entry.insert(tk.END,-500.0)
    ymax_entry.insert(tk.END,2000.0)
    logymin_entry.insert(tk.END,-15.0)
    logymax_entry.insert(tk.END,-8.0)
    
    
    
    #
    toolbar = NavigationToolbar2Tk(canvas,settingsframe,pack_toolbar=False)
    toolbar.grid(column = 1,row =0,padx = 0,pady = 0, sticky=tk.EW)
    #
    datasetinfoframe = GUI.MakeLabelFrame(frame=settingsframe,row=1,column=1,width=50,height=50,text='File info')
    datasetinfo = GUI.MakeText(datasetinfoframe,row=0,column=0,width=50,height=7)
    # feshbach_frame = create_labelframe(tab,'',row,column,1,1,200)
    # setfileframe = create_labelframe(feshbach_frame,'',0,0,1,1,200)
    # plotframe = create_labelframe(feshbach_frame,'',0,1,1,1,200)
    
    
    holdplot = tk.IntVar()
    entry_holdplot = tk.Checkbutton(leftsettings, text='Hold plot', variable=holdplot).grid(row=6,column = 0,sticky=tk.E)
    
    
    def set_file(file):
        f = open(file,'r')
        
        fullfile = np.loadtxt(file,max_rows = 2, delimiter = '&', dtype='unicode')
        
        datasetinfo.delete(1.0,tk.END)
        datasetinfo.tag_configure('entry', justify='left')
        datasetinfo.insert(tk.END,fullfile[0,0],'entry')
        
        data = np.loadtxt(file,skiprows=2)
        
        Bmin = float(Bmin_entry.get())
        Bmax = float(Bmax_entry.get())
        amin = float(ymin_entry.get())
        amax = float(ymax_entry.get())
        logamin = float(logymin_entry.get())
        logamax = float(logymax_entry.get())
        
        if holdplot.get() == False:
            ax.cla()
        Meanvalue = np.mean(np.absolute(data[:,1]))
		
        ax.plot(data[:,0],data[:,1],label=file)
        ax.set_xlim([Bmin,Bmax])
        ax.grid(True)
        ax.set_xlabel(fullfile[1,0])
        ax.set_ylabel(fullfile[1,1])
        ax.legend(loc='upper right')
        
        if fullfile[0,1] == 'log':
            ax.set_yscale('log')
            ax.set_ylim([np.exp(logamin*np.log(10)),np.exp(logamax*np.log(10))])
            
        else:
            ax.set_ylim([amin,amax])
            
        figure.tight_layout()
        canvas.draw()
        #feshfigure.tight_layout()
        
        return
    
    def set_filehelper(file):
        set_file(file)
        return
    
    def set_filelist(atom):
        nonlocal dataset_frame
        dataset_frame.destroy()
        dataset_frame = GUI.MakeLabelFrame(frame=choiceframe,row=1,column=0,width=190,height=50,text='Select data')
        
        filelist = glob.glob(atom+"_*.txt")
        for i in np.arange(len(filelist)):
            GUI.MakeButton(dataset_frame,row=i+1,column=0,width=25,text=filelist[i][len(atom)+1:-4],command =partial(set_filehelper,filelist[i]),state='normal',relief='raised')
            
        return 

  
    # atombuttonframe = create_labelframe(setfileframe,'Select atoms',0,0,1,1,200)
    
    buttonwidth = 5
    
    K39button0 = GUI.MakeButton(atoms_frame,row=0,column=0,width=buttonwidth,text='39K',command =lambda: set_filelist('39K'),state='normal',relief='raised')
    K40button0 = GUI.MakeButton(atoms_frame,row=0,column=1,width=buttonwidth,text='40K',command =lambda: set_filelist('40K'),state='normal',relief='raised')
    K41button0 = GUI.MakeButton(atoms_frame,row=0,column=2,width=buttonwidth,text='41K',command =lambda: set_filelist('41K'),state='normal',relief='raised')
    R87bbutton0 = GUI.MakeButton(atoms_frame,row=0,column=3,width=buttonwidth,text='87Rb',command =lambda: set_filelist('87Rb'),state='normal',relief='raised')
    Csbutton0 = GUI.MakeButton(atoms_frame,row=1,column=0,width=buttonwidth,text='133Cs',command =lambda: set_filelist('133Cs'),state='normal',relief='raised')
    KCsbutton0 = GUI.MakeButton(atoms_frame,row=1,column=1,width=buttonwidth,text='39KCs',command =lambda: set_filelist('39KCs'),state='normal',relief='raised')
    
    return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
