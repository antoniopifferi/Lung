# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 11:56:40 2020

@author: anton
"""
from numpy import *
from scipy import *
from matplotlib.pyplot import *
from pandas import *
close('all')

# OPTIONS
PLOT_TYPE1=True
Opt=['Mua','Mus']
YLABEL={'Mua':'absorption (cm-1)','Mus':'reduced scattering (cm-1)'}
XLABEL='time (s)'

# PARAMETERS
PathAnalysis='E:\\Beta\\Analysis\\Polmone\\FIT\\'
FileData='POLm0080new.txt'
FileKey='keyPOLm0080.txt'

# CONVERSION FUNCTION
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

# LOAD DICT
dataKey=read_csv(PathAnalysis+FileKey,sep='\t')
dcKey=dict(zip(dataKey.Key, dataKey.Value))

# LOAD DATA
Data=read_csv(PathAnalysis+FileData,sep='\t')
Data.rename(columns=dcKey,inplace=True)

# FILT DATA

# PLOT TYPE1
Color=['red','blue']
Linestyle=['-','--']
if PLOT_TYPE1:
    for od in Data.Detector.unique():
        for os in Data[Data.Detector==od].Subject.unique():
            fig=figure(figsize=cm2inch(40, 15))
            protocol=Data.Protocol.unique()
            for ip,op in enumerate(protocol):
                ax=fig.add_subplot(1,len(protocol),1+ip)
                for io,oo in enumerate(Opt):
                    table=Data[(Data.Detector==od)&(Data.Subject==os)&(Data.Protocol==op)].pivot_table(Opt,index='Time',columns='Repetition',aggfunc='mean')
                    table[oo].plot(ax=ax,secondary_y=(oo=='Mus'),style=Linestyle,color=Color[io])
                    ylabel(YLABEL[oo],color=Color[io])
                xlabel(XLABEL)
                title('det='+od+' - subj='+os)
                grid(True)
            fig.tight_layout()
            show()    