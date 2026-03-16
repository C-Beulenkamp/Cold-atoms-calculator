import tkinter as tk                    
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)

def scientificnotation(value):
	return str('{:.4e}'.format(value))
    
    
def DecimalNotation(value,decimals):
	return str(np.around(value,decimals=decimals))

def MakeButton(frame,row,column,width,text,command,state,relief='raised'):
	button = tk.Button(frame,text =text,command = command,bd=3,width=width,state=state,relief=relief)
	button.grid(column = column,row =row,padx = 0,pady = 0)
	return button

def MakeCanvas(frame,row,column,width,height):
	canvas = tk.Canvas(frame,
					   width = width,
					   height= height)
	canvas.grid(column = column,row =row,padx = 0,pady = 0)

	return canvas

def MakeCheckButton(frame,row,column,text,command,state,variable):
	checkbutton = tk.Checkbutton(frame,text = text,command= command,state = state,variable = variable)
	checkbutton.grid(column = column,row =row,padx = 0,pady = 0)
	return checkbutton
	
def MakeEntry(frame,row,column,width,variable,state='normal'):
	entry = tk.Entry(frame,
					 textvariable = variable,
					 width = width,
                     justify = 'right',
					 state=state)
	entry.grid(column = column,row =row,padx = 0,pady = 0)
	return entry

def MakeLabel(frame,row,column,width,text,anchor='center'):
	label = ttk.Label(frame,
					  text =text,
					  anchor = anchor,
					  width = width)
	label.grid(column = column,row =row,padx = 0,pady = 0)
	return label
	

def Make_DropDownMenu(frame,row,column,rowspan,columnspan,choices):
    variable = tk.StringVar(frame,choices[0])
    
    dropdown = tk.OptionMenu(frame, variable, *choices)
    dropdown.grid(row=row,column=column,rowspan=rowspan,columnspan=columnspan, sticky=tk.E)

    return variable,dropdown
    
def MakeDropDownCheckMenu(frame,row,column,width,text):
	menubutton = tk.Menubutton(frame,
                           text = text,
                           width = width)
	menubutton.grid(column = column,row =row,padx = 30,pady = 30)
	menubutton.menu =  tk.Menu ( menubutton, tearoff = 0 )
	menubutton["menu"] =  menubutton.menu
	return menubutton
	
def MakeFrame(frame,row,column,sticky='ns'):
	frame = tk.Frame(frame)
	frame.grid(column = column,row =row,padx = 5,pady = 5,sticky=sticky)
	return frame
	
def MakeScrollBar(frame,row,column,width):
	scrollbar = tk.Scrollbar(frame,width = width)
	scrollbar.grid(column = column,row =row,padx = 0,pady = 0,sticky='ns')
	return scrollbar
		
def MakeText(frame,row,column,width,height):
    text = tk.Text(frame,width = width,height =height,)
    text.grid(column = column,row =row,padx = 0,pady = 0,sticky='ns')
    return text

def MakeLabelFrame(frame,row,column,width,height,text=''):
    labelframe = tk.LabelFrame(frame,text=text,width = width,height =height)
    labelframe.grid(column = column,row =row,padx = 5,pady = 5)
    return labelframe

def FullEntry(frame,row,column,sizes,description='',unit=''):
    subframe = tk.Frame(frame)
    subframe.grid(column = column,row =row)
    labelwidth, entrywidth, unitwidth = sizes[0], sizes[1], sizes[2]
    MakeLabel(frame=subframe,row=0,column=0,width=labelwidth,text=description,anchor='e')
    variable = ''
    entry = MakeEntry(frame=subframe,row=0,column=1,width=entrywidth,variable=variable,state='normal')
    MakeLabel(frame=subframe,row=0,column=2,width=unitwidth,text=unit,anchor='e') 
    return variable, entry