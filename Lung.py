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
if PLOT_TYPE1:
    for id in Data.Detector.unique():
        for iss in Data[Data.Detector==id].Subject.unique():
            fig=figure(figsize=cm2inch(40, 15))
            table=Data[(Data.Detector==id)&(Data.Subject==iss)].pivot_table(Opt,index='Time',columns=['Protocol'],aggfunc='mean')
            for io in Opt:
                fig.add_subplot(1,len(Opt),1+Opt.index(io))
                plot(table[io])
                legend(table[io]), xlabel(XLABEL), ylabel(YLABEL[io])
                title('det='+id+' - subj='+iss+' - opt='+io)
                grid(True)
    show()