# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 11:56:40 2020

@author: anton
"""
from numpy import *
#from scipy import *
from matplotlib.pyplot import *
from pandas import *

# OPTIONS
PRODUCTION=True
Prot=2 # 1=PROT10, 2=PROT5
REFOLDING=True
SAVEFIG=True
PLOT_TYPE1=False # only Mua
PLOT_TYPE2=False # only Mua folding average
PLOT_TYPE3=True # all param
PLOT_SPECTRA=True # plot the broadband spectra
PLOT_DTOF=True # plot the DTOF
Opt=['Mua','Mus']
Gates=['DeltaGateNorm02','DeltaGateNorm04','DeltaGateNorm06','DeltaGateNorm08','DeltaGateNorm10']
sMeanGateIn=['MeanGateIn00','MeanGateIn01','MeanGateIn02','MeanGateIn03','MeanGateIn04','MeanGateIn05','MeanGateIn06','MeanGateIn07','MeanGateIn08','MeanGateIn09','MeanGateIn10','MeanGateIn11','MeanGateIn12','MeanGateIn13','MeanGateIn14','MeanGateIn15']
sMeanGateOut=['MeanGateOut00','MeanGateOut01','MeanGateOut02','MeanGateOut03','MeanGateOut04','MeanGateOut05','MeanGateOut06','MeanGateOut07','MeanGateOut08','MeanGateOut09','MeanGateOut10','MeanGateOut11','MeanGateOut12','MeanGateOut13','MeanGateOut14','MeanGateOut15']
sMeanGateDiff=['MeanGateDiff00','MeanGateDiff01','MeanGateDiff02','MeanGateDiff03','MeanGateDiff04','MeanGateDiff05','MeanGateDiff06','MeanGateDiff07','MeanGateDiff08','MeanGateDiff09','MeanGateDiff10','MeanGateDiff11','MeanGateDiff12','MeanGateDiff13','MeanGateDiff14','MeanGateDiff15']
YLABEL={'Mua':'absorption (cm-1)','Mus':'reduced scattering (cm-1)'}
XLABEL='time (s)'
COMP_LIST=['HHb','O2Hb','tHb','SO2','Lipid','H2O','Coll','Tot']
FIRST_LAMBDA=610
LAMBDA0=600
LAMBDA1=700
LAMBDA2=900


# PARAMETERS
PATHBETA='C:\\OneDrivePolimi\\OneDrive - Politecnico di Milano\\Beta\\'
PATHDATA='Data\\Polmone\\'
PATHANALYSIS='Analysis\\Polmone\\FIT\\'
EXTDATA='.dat'
FILEFIT='POLm0080new2.txt'
FILEKEY='keyPOLm0080.txt'
FILEKEYSPECTRA='keySpectra.txt'
FILESPECTRA='FileSpectraAll.txt'
FileComponents='Components0.txt'
FILECOMPOUT='CompOut.txt'
PROT10=1
PROT5=2
sProt=['Prot10','Prot05']
PROT10_PERIOD=10
PROT10_BASE=10
PROT5_PERIOD=5
NUMGATE=16
WIDTHGATE=500 # ps
MINREF=100 # used to avoid divide by zero (could be set higher)
InClock=[]
OutClock=[]
InClock.insert(PROT10,range(2,5))
InClock.insert(PROT5,range(3,10))
OutClock.insert(PROT10,range(7,10))
OutClock.insert(PROT5,range(13,20))

# CONSTANTS
NUMCHAN=4096
GAIN=4
HEADLEN=764
SUBLEN=204
ch2ps=50000/GAIN/NUMCHAN
# if Prot==1, InClock=3:5; else InClock=4:10; end
# if Prot==1, OutClock=8:10; else OutClock=14:20; end
# if strcmp(System,'HYBD'), BkgCh1=330; else BkgCh1=500; end
# if strcmp(System,'HYBD'), BkgCh2=700; else BkgCh2=800; end

# FONT SIZE FOR PLOTS
if PRODUCTION:
    SMALL_SIZE = 10
    MEDIUM_SIZE = 12
    BIGGER_SIZE = 14
    matplotlib.rc('font', size=SMALL_SIZE)          # controls default text sizes
    matplotlib.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    matplotlib.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    matplotlib.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    matplotlib.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    matplotlib.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    matplotlib.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# DIMENSION OF FIGURES
if PRODUCTION:
    FIGWIDTH = 40
else: 
    FIGWIDTH = 50

# CONVERSION FUNCTION
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

# LOAD FIT DATA
dataKey=read_csv(PATHBETA+PATHANALYSIS+FILEKEY,sep='\t')
dcKey=dict(zip(dataKey.Key, dataKey.Value))
Data=read_csv(PATHBETA+PATHANALYSIS+FILEFIT,sep='\t')
Data.rename(columns=dcKey,inplace=True)

# CREATES COL-FIELDS INITIALISED TO 0
for ig in arange(NUMGATE):
    Data['DeltaGateNorm'+str(ig).zfill(2)]=0
# DELETE???
    
# # CREATES COL-FIELDS FOR MEAN GATES
# sMeanGateIn=[]
# sMeanGateOut=[]
# sMeanGateDiff=[]

# for ig in arange(NUMGATE):
#     sMeanGateIn.append(['MeanGateIn'+str(ig)])
#     sMeanGateOut.append(['MeanGateOut'+str(ig)])
#     sMeanGateDiff.append(['MeanGateDiff'+str(ig)])
# sMeanGateIn=tuple(sMeanGateIn)

# LOAD RAW DATA AND PROCESS THEM TO POPULATE DataFrame
sumDtof=zeros([len(Data.Subject.unique()),NUMCHAN])
for od in Data.Data.unique():
    item=Data[Data.Data==od].iloc[0,:]
    
    # Allocate memory
    Dtof=zeros([item.NumRep,item.NumClock,NUMCHAN]) # DTOF for this file
    Gate=zeros([item.NumRep,item.NumClock,NUMGATE]) # GATES for this file
    DeltaGateNorm=zeros([item.NumRep,item.NumClock,NUMGATE]) # GATES normalised to REF(=mean over all instants) for this file
    MeanGateIn=zeros([NUMGATE]) # GATES averaged for all IN instants for this file
    
    # load system
    fid = open(PATHBETA+PATHDATA+item.Syst+EXTDATA,'rb')
    fid.seek(HEADLEN,0) # 0=beginning
    fid.seek(SUBLEN,1) # 1=current
    Syst=fromfile(fid,dtype=uint16,count=NUMCHAN)
    fid.close()
    
    # load data
    fid = open(PATHBETA+PATHDATA+item.Data+EXTDATA,'rb')
    fid.seek(HEADLEN,0) # 0=beginning
    for ir in arange(item.NumRep):
        for ik in arange(item.NumClock):
            fid.seek(SUBLEN,1) # 1=current
            Dtof[ir,ik,:]=fromfile(fid,dtype=uint16,count=NUMCHAN)
            
    # Calc peak
    Chan0=Syst.tolist().index(max(Syst))
    mtime=(arange(NUMCHAN)-Chan0)*ch2ps
    
    # Collect DTOF for Figure
    if (item.Protocol==PROT10)&(item.Position=='UR')&(item.Detector=='HYBD'):
        iSubject=Data.Subject.unique().tolist().index(item.Subject)
        sumDtof[iSubject,:]+=sum(Dtof[1,:,:],axis=0)
    

# % calculate
# DataIn=squeeze(mean(mean(Data(FirstRep:end,InClock,:),2),1));
# DataOut=squeeze(mean(mean(Data(FirstRep:end,OutClock,:),2),1));
# DataClock=squeeze(mean(Data(FirstRep:end,:,:),1));

    # Calc gate
    StopGate=(Chan0+arange(NUMGATE)*(WIDTHGATE/ch2ps)).astype(int)
    StartGate=StopGate[arange(NUMGATE)-1]+1 # stop ends exacly before the start
    StartGate=insert(StartGate,0,Chan0) # add first gate;
    TimeGate=arange(NUMGATE)*WIDTHGATE
    for ig in arange(NUMGATE):
        Gate[:,:,ig]=sum(Dtof[:,:,StartGate[ig]:StopGate[ig]],axis=2)
        
    # Calc NormGate/Ref
    for ir in arange(item.NumRep):
        Ref=mean(Gate[ir,:,:],axis=0) # nore: Gate here is 2D since ir reduce dimensions
        Ref[Ref<MINREF]=0
        for ik in arange(item.NumClock):
            DeltaGateNorm[ir,ik,:]=nan_to_num(-log(Gate[ir,ik,:]/Ref))
    for ig in arange(NUMGATE):    
        Data.loc[Data.Data==od,'DeltaGateNorm'+str(ig).zfill(2)]=DeltaGateNorm[:,:,ig].flatten()    
        
    # Calc MeanGate
    MeanGateIn=mean(mean(DeltaGateNorm[item.FirstRep-1:,InClock[item.Protocol-1],:],axis=1),axis=0)
    MeanGateOut=mean(mean(DeltaGateNorm[item.FirstRep-1:,OutClock[item.Protocol-1],:],axis=1),axis=0)
    MeanGateDiff=MeanGateOut-MeanGateIn
    #DeltaMuaClock=squeeze(mean(DeltaMua(FirstRep:end,:,:),1));
    for ig in arange(NUMGATE):    
        Data.loc[Data.Data==od,'MeanGateIn'+str(ig).zfill(2)]=MeanGateIn[ig]    
        Data.loc[Data.Data==od,'MeanGateOut'+str(ig).zfill(2)]=MeanGateOut[ig]    
        Data.loc[Data.Data==od,'MeanGateDiff'+str(ig).zfill(2)]=MeanGateDiff[ig]    

# FILT DATA
Data.replace([inf, -inf, nan], 0)


# FOLDING AVERAGE
Data['RefTime']=0
if REFOLDING:
    Data.RefTime=(Data.Time)%Data.NumClock
else:
    Data.RefTime=Data.Time

for i in range(20):
    print('\n')

# PLOT TYPE1
Color=['red','blue']
Linestyle=['-','--']
if PLOT_TYPE1:
    for od in Data.Detector.unique():
        for os in Data[Data.Detector==od].Subject.unique():
            fig=figure(figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
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
            
# PLOT TYPE2
Color=['purple','blue']
Linestyle=['-','--']
if PLOT_TYPE2:
    for od in Data.Detector.unique():
        for os in Data[Data.Detector==od].Subject.unique():
            fig=figure(figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
            protocol=Data.Protocol.unique()
            for ip,op in enumerate(protocol):
                ax=fig.add_subplot(1,len(protocol),1+ip)
                for io,oo in enumerate(Opt):
                    table=Data[(Data.Detector==od)&(Data.Subject==os)&(Data.Protocol==op)].pivot_table(Opt,index='RefTime',aggfunc='mean')
                    table[oo].plot(ax=ax,secondary_y=(oo=='Mus'),style=Linestyle,color=Color[io])
                    ylabel(YLABEL[oo],color=Color[io])
                xlabel(XLABEL)
                title('det='+od+' - subj='+os)
                grid(True)
            fig.tight_layout()
            
# PLOT TYPE3
Color=['red','blue']
Linestyle=['-','--']
if PLOT_TYPE3:
    pData=Data[(Data.Detector=='HYBD')&(Data.Protocol==Prot)]
    position=pData.Position.unique()
    subject=pData.Subject.unique()
    np=len(position)
    ns=len(subject)
    figOpt=figure(figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
    figGate=figure(figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
    figMean=figure(figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
    gsOpt=figOpt.add_gridspec(np,ns, hspace=0, wspace=0)
    gsGate=figGate.add_gridspec(np,ns, hspace=0, wspace=0)
    gsMean=figMean.add_gridspec(np,ns, hspace=0, wspace=0)
    axsOpt=gsOpt.subplots(sharex=True)
    axsGate=gsGate.subplots(sharex=True,sharey=True)
    axsMean=gsMean.subplots(sharex=True,sharey=True)
    
    for ip,op in enumerate(position):
        for iss,os in enumerate(subject):
            
            # plot Opt
            sca(axsOpt[ip,iss])
            for io,oo in enumerate(Opt):
                table=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(Opt,index='RefTime',aggfunc='mean')
                #table[oo].plot(ax=axsOpt[ip,iss],secondary_y=(oo=='Mus'),style=Linestyle,color=Color[io])
                table[oo].plot(secondary_y=(oo=='Mus'),style=Linestyle,color=Color[io])
                if ((oo=='Mua')&(iss==0)): ylabel(YLABEL[oo],color=Color[io])
                if ((oo=='Mus')&(iss==(ns-1))): ylabel(YLABEL[oo],color=Color[io])
                tick_params(axis="y",direction="in", pad=-35, labelcolor=Color[io])
                gca().yaxis.set_major_locator(MaxNLocator(3))
                #gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                #gca().yaxis.set_major_locator(MultipleLocator(0.1))
            if ip==(np-1): xlabel(XLABEL)
            if ip==0: title('pos='+op+' - subj='+os)
            grid(True)
            
            # plot Gate
            sca(axsGate[ip,iss])
            table=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(Gates,index='RefTime',aggfunc='mean')
            table.plot(ax=axsGate[ip,iss])
            xlabel(XLABEL)
            ylabel('log ratio to REF')
            title('pos='+op+' - subj='+os)
            grid(True)

            # plot Mean
            sca(axsMean[ip,iss])
            temp=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(sMeanGateIn,index='Detector',aggfunc='mean')
            mgIn=temp.to_numpy().transpose()
            temp=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(sMeanGateOut,index='Detector',aggfunc='mean')
            mgOut=temp.to_numpy().transpose()
            temp=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(sMeanGateDiff,index='Detector',aggfunc='mean')
            # mgDiff=temp.to_numpy().transpose()
            mgDiff=mgOut-mgIn
            # A=column_stack(mgIn,mgOut,mgDiff)
            plot(TimeGate,mgIn,label='IN')
            plot(TimeGate,mgOut,label='OUT')
            plot(TimeGate,mgDiff,label='OUT-IN')
            legend()
            xlabel('time gate (ps)')
            ylabel('log ratio to REF')
            title('pos='+op+' - subj='+os)
            grid(True)
                
    figOpt.suptitle("OPTICAL PROPERTIES - Prot = "+sProt[Prot-1]+", Folding = "+str(REFOLDING), fontsize=16)
    figGate.suptitle("NORMALISED GATES - Prot = "+sProt[Prot-1]+", Folding = "+str(REFOLDING), fontsize=16)
    figMean.suptitle("MEAN GATES OVER PHASES - Prot = "+sProt[Prot-1]+", Folding = "+str(REFOLDING), fontsize=16)
    figOpt.tight_layout()
    figGate.tight_layout()
    figMean.tight_layout()
    if SAVEFIG:
        figOpt.savefig(PATHBETA+PATHANALYSIS+'FigOpt.eps',format='eps')
        figGate.savefig(PATHBETA+PATHANALYSIS+'FigGate.eps',format='eps')
        figMean.savefig(PATHBETA+PATHANALYSIS+'FigMean.eps',format='eps')


# plot DTOF
if PLOT_DTOF:
    figDtof=figure()
    sumDtof=sumDtof.transpose()
    semilogy(mtime,sumDtof)
    legend(Data.Subject.unique())
    xlabel('time (ps)')
    ylabel('counts')
    xlim([0,8000])
    grid(True)
    if SAVEFIG:
        figDtof.savefig(PATHBETA+PATHANALYSIS+'FigDTOF.eps',format='eps')

# plot SPECTRA
if PLOT_SPECTRA:
    
    # LOAD FIT DATA
    dataKey=read_csv(PATHBETA+PATHANALYSIS+FILEKEYSPECTRA,sep='\t')
    dcKey=dict(zip(dataKey.Key, dataKey.Value))
    Spectra=read_csv(PATHBETA+PATHANALYSIS+FILESPECTRA,sep='\t')
    Spectra.rename(columns=dcKey,inplace=True)
    Spectra['Position']=0
    Spectra.loc[(Spectra.Region=='ANTERIORE')&(Spectra.Side=='DX'),'Position']='UR'
    Spectra.loc[(Spectra.Region=='ANTERIORE')&(Spectra.Side=='SX'),'Position']='UL'
    Spectra.loc[(Spectra.Region=='POSTERIORE')&(Spectra.Side=='DX'),'Position']='DR'
    Spectra.loc[(Spectra.Region=='POSTERIORE')&(Spectra.Side=='SX'),'Position']='DL'

    # plot    
    #pSpectra=Spectra[(Spectra.Fit=='-0.8-0.01')&(Spectra.Organ=='POLMONE')]
    pSpectra=Spectra[(Spectra.Fit!='0.3-0.001')&(Spectra.Organ=='POLMONE')]
    pSpectra = pSpectra[pSpectra.Lambda!=920]
    pSpectra = pSpectra[pSpectra.Lambda!=600]
    figSpectra=figure(figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
    # plot Opt
    for io,oo in enumerate(Opt):
        ax=figSpectra.add_subplot(1,2,io+1)
        sca(ax)
        table=pSpectra.pivot_table(values=oo,index='Lambda',columns='Subject',aggfunc='mean')
        table.plot(ax=ax,marker='D')
        xlabel('Wavelength (nm)')
        ylabel(YLABEL[oo])
        if oo=='Mua':
            ylim([0,0.5])
        else:
            ylim([0,15])
        grid(True)       
    figSpectra.tight_layout()
    if SAVEFIG:
        figSpectra.savefig(PATHBETA+PATHANALYSIS+'FigSpectra.eps',format='eps')

#%% CALC COMPONENTS
Components=read_table(PATHBETA+PATHANALYSIS+FileComponents)
figure()
Components.plot(x='Lambda')
yscale('log')
ylim([0,0.5]), title('Components'), xlabel('wavelength (nm)'), ylabel('specific absorption (cm-1)')
show()

table=pSpectra.pivot_table(values=Opt,index='Lambda',columns='Subject',aggfunc='mean')
comp=Components[Components['Lambda'].isin(pSpectra.Lambda.unique())].values
comp=delete(comp,0,1)
aComp=linalg.lstsq(comp,table['Mua'],rcond=None)[0] #[0] to extract m-coeff in lstsq
y=log(table.loc[LAMBDA1:LAMBDA2,'Mus'])
x=log(y.index/LAMBDA0)
model=polyfit(x,y,1)
b=-model[0]
a=exp(model[1])
#A = vstack([x, np.ones(len(x))]).T
dfComp=DataFrame(data=aComp.transpose(),columns=Components.columns[1:],index=table.Mua.columns)
dfComp['tHb']=dfComp['HHb']+dfComp['O2Hb']
dfComp['SO2']=dfComp['O2Hb']/dfComp['tHb']
dfComp['Tot']=dfComp['Lipid']+dfComp['H2O']+dfComp['Coll']
dfComp['FitComp']='LambdaFit'
dfComp['a']=a
dfComp['b']=b
dfComp.plot()
dfComp.to_csv(path_or_buf=PATHBETA+PATHANALYSIS+FILECOMPOUT,sep='\t')
figure()
plot(x,y)
#filtData=merge(filtData,dfComp,on=['Subject','Meas','Rho'])

show()



