import tkinter as tk                    
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
import GUIfunctions as GUI

def GUI_Examples(frame):


    #Setting a style
    #s = ttk.Style()
    #s.theme_use('alt')
    #print(s.theme_names()) #print a list of available styles.

    #Define the tabs
    # tabControl = ttk.Notebook(root)
    # frame = ttk.Frame(tabControl)
    # tab2 = ttk.Frame(tabControl)
    # tab3 = ttk.Frame(tabControl)
      
    # tabControl.add(frame, text ='Tab 1')
    # tabControl.add(tab2, text ='Tab 2')
    # tabControl.add(tab3, text ='Tab 3')
    # tabControl.pack(expand = 1, fill ="both")

        
    #Tab 1
    #Button
    def dummyfunction():
        return print('pressed button')
    button = GUI.MakeButton(frame=frame,row=0,column=0,width=50,text='example button',command=dummyfunction,relief='raised',state='active')

    #Canvas
    canvas = GUI.MakeCanvas(frame=frame,width=100,height=20,column=0,row=1)

    #checkbutton
    checkbutton_variable = 0
    checkbutton = GUI.MakeCheckButton(frame=frame,row=2,column=0,text='checkbutton',command=dummyfunction,state='active',variable=checkbutton_variable)
    checkbutton.select()
    checkbutton.deselect()

    #Entry
    entryvariable = 'test'
    entry = GUI.MakeEntry(frame=frame,row=3,column=0,width=40,variable=entryvariable)
    entry.insert(0,'example text entry')
    entry.get()
    print(entry.get())


    #Text label
    label = GUI.MakeLabel(frame=frame,row=4,column=0,width=100,text="example label \n second line",anchor='e')

    #menubutton
    menubutton = GUI.MakeDropDownCheckMenu(frame=frame,row=5,column=0,width=100,text='options')
    CsVar = tk.IntVar()
    KVar = tk.IntVar()
    menubutton.menu.add_checkbutton ( label="Cs",variable=CsVar )
    menubutton.menu.add_checkbutton ( label="K",variable=KVar )

    #Seperate subframe
    frame2 = GUI.MakeFrame(frame=frame,column=0,row=6)

    #Scrollbar in this frame
    scrollbar = GUI.MakeScrollBar(frame=frame2,column=1,row=0,width=20)
    mylist = tk.Listbox(frame2, yscrollcommand = scrollbar.set,width = 50,height=5)
    for line in range(100):
       mylist.insert('end', f"Variable {line} = {line*line}" )
       

    mylist.grid(column = 0,row =0)
    scrollbar.config( command = mylist.yview ) 

    #Text
    text = tk.Text(frame,
                   width = 20 ,
                   height =3)
    text.grid(column = 0,row =7,padx = 0,pady = 0,sticky='ns')

    text.insert('end','1 \n')
    text.insert('end','2 \n')
    text.delete(1.0,'end')
    text.insert('end','3 \n')

    #Labelframe
    labelframe = tk.LabelFrame(frame,text='labelframe',
                               width=100,
                               height=50)
    labelframe.grid(column = 0,row =8,padx = 10,pady = 10)
    tk.Label(labelframe,text='labelframe').grid(column=0,row=0)

    #Matplotlib figure
    figureframe = tk.Frame(frame)
    figureframe.grid(column = 0,row =10,padx = 0,pady = 0,sticky='ns')

    #figure and canvas
    figure = plt.Figure(figsize=(3,3), dpi=100)
    figure.set_tight_layout(True)
    ax = figure.add_subplot(111)
    plotcanvas= FigureCanvasTkAgg(figure, figureframe)
    plotcanvas.get_tk_widget().grid(column = 0,row =0,padx = 10,pady = 10)

    #axes
    ax.plot([0,1],[0,1])

    #toolbar
    toolbar = NavigationToolbar2Tk(plotcanvas,figureframe,pack_toolbar=False)
    toolbar.grid(column = 0,row =1,padx = 10,pady = 10)

    ax.cla()
    ax.plot([1,0],[0,1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    tmp = 0.0
    #Update plot button
    def updateplot():
        nonlocal tmp
        tmp += 1
        ax.cla()
        ax.plot([tmp+1,0],[0,tmp+1])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plotcanvas.draw()

        return print('pressed button')

    plotbutton = tk.Button(frame,
                       text ="Update plot",
                       command = updateplot,
                       bd=3,
                       width=50,
                       state='active',
                       relief='raised')
    plotbutton.grid(column = 0,row =11,padx = 0,pady = 0)

    ####
    #root.mainloop()  
