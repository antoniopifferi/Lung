# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 11:56:40 2020

@author: anton
"""
from numpy import *
#from scipy import *
from matplotlib.pyplot import *
from pandas import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker

# OPTIONS
PRODUCTION=False
FIG_PW2021=False
FIG_PAPER=True
BKG=False
Prot=2 # 1=PROT5, 2=PROT10 (PROT5=10 rep with 5+5s, PROT10=5 rep with 10+10s)
if Prot==1: VERT_LINE=5
if Prot==2: VERT_LINE=10
MIN_COUNT_RATE=10000 # minimum number of photons for late gate (1s acquisition)
MAX_TIME=8000 # ps time scale
REFOLDING=False
SAVEFIG=False
PLOT_TYPE1=False # only Mua
PLOT_TYPE2=False # only Mua folding average
PLOT_TYPE3=True # all param
PLOT_SPECTRA=True # plot the broadband spectra
PLOT_DTOF=True # plot the DTOF
PLOT_DEPTH=True # plot the DTOF
Opt=['Mua','Mus']
Gates=['DeltaGateNorm02','DeltaGateNorm04','DeltaGateNorm06','DeltaGateNorm08','DeltaGateNorm10']
sMeanGateIn=['MeanGateIn00','MeanGateIn01','MeanGateIn02','MeanGateIn03','MeanGateIn04','MeanGateIn05','MeanGateIn06','MeanGateIn07','MeanGateIn08','MeanGateIn09','MeanGateIn10']
sMeanGateOut=['MeanGateOut00','MeanGateOut01','MeanGateOut02','MeanGateOut03','MeanGateOut04','MeanGateOut05','MeanGateOut06','MeanGateOut07','MeanGateOut08','MeanGateOut09','MeanGateOut10']
sMeanGateDiff=['MeanGateDiff00','MeanGateDiff01','MeanGateDiff02','MeanGateDiff03','MeanGateDiff04','MeanGateDiff05','MeanGateDiff06','MeanGateDiff07','MeanGateDiff08','MeanGateDiff09','MeanGateDiff10']
YLABEL={'Mua':'$\mu_a$ (cm$^{-1}$)','Mus':'$\mu_s\prime$ (cm$^{-1}$)'}
XLABEL='time (s)'
COMP_LIST=['HHb','O2Hb','tHb','SO2','Lipid','H2O','Coll','Tot']
FIRST_LAMBDA=610
LAMBDA0=600
LAMBDA1=700
LAMBDA2=900
LAMBDA_DEPTH=820 # wavelength where the depth is calculated
MUS0_DEPTH=10 # (cm-1) pivot scattering for the empirical formula on Zmax (depth)
BKG_F=2580
BKG_L=2615

# PARAMETERS
PATHBETA='C:\\OneDrivePolimi\\OneDrive - Politecnico di Milano\\Beta\\'
PATHDATA='Data\\Polmone\\'
PATHANALYSIS='Analysis\\Polmone\\FIT\\'
EXTDATA='.dat'
FILEFIT='POLm0080new2.txt'
FILEKEY='keyPOLm0080.txt'
FILEKEYSPECTRA='keySpectra.txt'
FILEKEYSUBJECT='keySubject.txt'
FILESPECTRA='FileSpectraAll.txt'
FileComponents='Components0.txt'
FILECOMPOUT='CompOut.txt'
PROT05=1
PROT10=2
NUM_REPLICA_DTOF=2 #number of repeated measurements for each position
NUMGATE=11
WIDTHGATE=500 # ps
MINREF=100 # used to avoid divide by zero (could be set higher)
InClock=[]
OutClock=[]
InClock.insert(PROT05,range(2,5))
InClock.insert(PROT10,range(3,10))
OutClock.insert(PROT05,range(7,10))
OutClock.insert(PROT10,range(13,20))

# EMPIRICAL FORMULA FOR DEPTH: Zmax=2*a*t^b (n_int=1.4, n_ext=1) NOTE:2* because formula is for Zmax
aDepth=0.1511
bDepth=0.5433

# CONSTANTS
NUMCHAN=4096
GAIN=4
HEADLEN=764
SUBLEN=204
ch2ps=50000/GAIN/NUMCHAN

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

if FIG_PW2021 or FIG_PAPER:
    style.use(['science','no-latex'])
    FIGWIDTH = 35
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 16
    matplotlib.rc('font', size=SMALL_SIZE)          # controls default text sizes
    matplotlib.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    matplotlib.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    matplotlib.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    matplotlib.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    matplotlib.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    matplotlib.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# NAMING FILE
if Prot==1:
    sProt='Prot05s'
else:
    sProt='Prot10s'
if REFOLDING==True:
    sRefold='ReTrue'
else:
    sRefold='ReFals'
sProtRef=sProt+'_'+sRefold
    

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
sbjKey=read_csv(PATHBETA+PATHANALYSIS+FILEKEYSUBJECT,sep='\t')
dcKeyS=dict(zip(sbjKey.Key, sbjKey.Value))
Data['Subject'].replace(dcKeyS, inplace=True)



# CREATES COL-FIELDS INITIALISED TO 0
for ig in arange(NUMGATE):
    Data['DeltaGateNorm'+str(ig).zfill(2)]=0
# DELETE???
    
# LOAD RAW DATA AND PROCESS THEM TO POPULATE DataFrame
meanDtof=zeros([len(Data.Subject.unique()),NUMCHAN])
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
            if BKG: Dtof[ir,ik,:]=Dtof[ir,ik,:]-mean(Dtof[ir,ik,BKG_F:BKG_L])
            
    # Calc peak
    Chan0=Syst.tolist().index(max(Syst))
    mtime=(arange(NUMCHAN)-Chan0)*ch2ps
    
    # Collect DTOF for Figure
    if (item.Protocol==PROT05)&(item.Position=='UR')&(item.Repetition=='a')&(item.Detector=='HYBD'):
        iSubject=Data.Subject.unique().tolist().index(item.Subject)
        meanDtof[iSubject,:]=+sum(Dtof[0,:,:],axis=0)/(item.NumRep)
        
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
            DeltaGateNorm[ir,ik,:]=nan_to_num((Gate[ir,ik,:]-Ref)/Ref,nan=0.0, posinf=0.0, neginf=0.0)
            # DeltaGateNorm[ir,ik,:]=nan_to_num(-log(Gate[ir,ik,:]/Ref))
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
            
#%% PLOT TYPE3
Color=['red','blue']
Linestyle=['-','--']
if PLOT_TYPE3:
    pData=Data[(Data.Detector=='HYBD')&(Data.Protocol==Prot)]
    position=pData.Position.unique()
    subject=pData.Subject.unique()
    np=len(position)
    ns=len(subject)
    figOpt=figure(num='FigOpt',figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
    figGate=figure(num='FigGate',figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
    figMean=figure(num='FigMean',figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
    gsOpt=figOpt.add_gridspec(np,ns, hspace=0, wspace=0.5)
    gsGate=figGate.add_gridspec(np,ns, hspace=0, wspace=0.3)
    gsMean=figMean.add_gridspec(np,ns, hspace=0, wspace=0.3)
    axsOpt=gsOpt.subplots(sharex=True)
    axsGate=gsGate.subplots(sharex=True)
    axsMean=gsMean.subplots(sharex=True)
    
    for ip,op in enumerate(position):
        for iss,os in enumerate(subject):
            
            # plot Opt
            sca(axsOpt[ip,iss])
            for io,oo in enumerate(Opt):
                table=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(Opt,index='RefTime',aggfunc='mean')
                #table[oo].plot(ax=axsOpt[ip,iss],secondary_y=(oo=='Mus'),style=Linestyle,color=Color[io])
                table[oo].plot(secondary_y=(oo=='Mus'),style=Linestyle,color=Color[io])
                tick_params(axis="y",labelcolor=Color[io])
                if oo=='Mus':
                    gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    if ((iss==(ns-1))&(ip==0)): ylabel('UR - '+YLABEL[oo],color=Color[io])
                    if ((iss==(ns-1))&(ip==(np-1))): ylabel('DR - '+YLABEL[oo],color=Color[io])
                    axvline(VERT_LINE,color='gray',linewidth=2)
                else:
                    gca().yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                    if ((iss==0)&(ip==0)): ylabel('UR - '+YLABEL[oo],color=Color[io])
                    if ((iss==0)&(ip==(np-1))): ylabel('DR - '+YLABEL[oo],color=Color[io])
            grid(True)
            if (FIG_PW2021 or FIG_PAPER):
                if ip==(np-1): axsOpt[ip,iss].set_xlabel('time (s)')
                axsOpt[ip,iss].set_xlabel('time (s)')
                if ip==0: title('subject #'+str(iss+1))
            else:
                figOpt.suptitle("OPTICAL PROPERTIES - Prot = "+sProtRef, fontsize=16)
                if iss==(ns-1): xlabel(XLABEL)
                if iss==0: title('pos='+op+' - subj='+os)
                if ((oo=='Mua')&(iss==0)): ylabel(YLABEL[oo],color=Color[io])
                if ((oo=='Mus')&(iss==(ns-1))): ylabel(YLABEL[oo],color=Color[io])
            
            # plot Gate
            sca(axsGate[ip,iss])
            table=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(Gates,index='RefTime',aggfunc='mean')
            table.plot(ax=axsGate[ip,iss],legend=False)
            grid(True)
            if (FIG_PW2021 or FIG_PAPER):
                if ip==(np-1): xlabel('time (s)')
                if(((Prot==1)and(ip==1)and(iss==1))or((Prot==2)and(ip==1)and(iss==4))):
                    legend(['0.5 ns','1.5 ns','2.5 ns','3.5 ns','4.5 ns'])
                else:
                    axvline(VERT_LINE,color='gray',linewidth=2)
                if ip==0: title('subject #'+str(iss+1))
                if iss==0 and ip==0: ylabel('UR - contrast')
                if iss==0 and ip==(np-1): ylabel('DR - contrast')
            else:
                xlabel(XLABEL)
                ylabel('log ratio to REF')
                title('pos='+op+' - subj='+os)
                figGate.suptitle("NORMALISED GATES - Prot = "+sProtRef, fontsize=16)

            # plot Mean
            sca(axsMean[ip,iss])
            temp=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(sMeanGateIn,index='Detector',aggfunc='mean')
            mgIn=temp.to_numpy().transpose()
            temp=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(sMeanGateOut,index='Detector',aggfunc='mean')
            mgOut=temp.to_numpy().transpose()
            temp=pData[(pData.Position==op)&(pData.Subject==os)].pivot_table(sMeanGateDiff,index='Detector',aggfunc='mean')
            mgDiff=mgIn-mgOut
            plot(TimeGate,mgIn,label='IN')
            plot(TimeGate,mgOut,label='OUT')
            plot(TimeGate,mgDiff,label='IN-OUT')
            #legend()
            grid(True)
            if (FIG_PW2021 or FIG_PAPER):
                if((ip==1)and(iss==2)): legend()
                if ip==0: title('subject #'+str(iss+1))
                if ip==(np-1): xlabel('time gate (ps)')
                if iss==0 and ip==0: ylabel('UR - contrast')
                if iss==0 and ip==(np-1): ylabel('DR - contrast')
                # axsMean[ip,iss].yaxis.set_major_formatter(ticker.PercentFormatter())
            else:
                xlabel('time gate (ps)')
                ylabel('log ratio to REF')
                title('pos='+op+' - subj='+os)
                figMean.suptitle("MEAN GATES OVER PHASES - Prot = "+sProtRef, fontsize=16)
                
    figOpt.tight_layout()
    figGate.tight_layout()
    figMean.tight_layout()
    if SAVEFIG:
        figOpt.savefig(PATHBETA+PATHANALYSIS+'FigOpt'+sProtRef+'.eps',format='eps')
        figGate.savefig(PATHBETA+PATHANALYSIS+'FigGate'+sProtRef+'.eps',format='eps')
        figMean.savefig(PATHBETA+PATHANALYSIS+'FigMean'+sProtRef+'.eps',format='eps')
    if (FIG_PW2021 or FIG_PAPER):
        figOpt.savefig(PATHBETA+PATHANALYSIS+'FigOpt'+sProtRef+'.jpg',format='jpg')
        figGate.savefig(PATHBETA+PATHANALYSIS+'FigGate'+sProtRef+'.jpg',format='jpg')
        figMean.savefig(PATHBETA+PATHANALYSIS+'FigMean'+sProtRef+'.jpg',format='jpg')

# plot DTOF
if PLOT_DTOF:
    figDtof=figure(num='FigDtof',figsize=cm2inch(0.4*FIGWIDTH,0.4*FIGWIDTH))
    meanDtof=meanDtof.transpose()
    semilogy(mtime,meanDtof)
    legend(Data.Subject.unique())
    xlabel('time (ps)')
    ylabel('mean counts/s')
    xlim([0,MAX_TIME])
    grid(True)
    legend([r'#1 $\rho=7$ cm',r'#2 $\rho=7$ cm',r'#3 $\rho=7$ cm',r'#4 $\rho=7$ cm',r'#5 $\rho=9$ cm'])
    if SAVEFIG:
        figDtof.savefig(PATHBETA+PATHANALYSIS+'FigDTOF.eps',format='eps')
    if (FIG_PW2021 or FIG_PAPER):
        figDtof.savefig(PATHBETA+PATHANALYSIS+'FigDTOF.jpg',format='jpg')

# plot SPECTRA
if PLOT_SPECTRA:
    
    # LOAD FIT DATA
    dataKey=read_csv(PATHBETA+PATHANALYSIS+FILEKEYSPECTRA,sep='\t')
    dcKey=dict(zip(dataKey.Key, dataKey.Value))
    Spectra=read_csv(PATHBETA+PATHANALYSIS+FILESPECTRA,sep='\t')
    Spectra.rename(columns=dcKey,inplace=True)
    Spectra['Subject'].replace(dcKeyS, inplace=True)
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
    figSpectra=figure(num='FigSpectra',figsize=cm2inch(FIGWIDTH,0.4*FIGWIDTH))
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
    if (FIG_PW2021 or FIG_PAPER):
        figSpectra.savefig(PATHBETA+PATHANALYSIS+'FigSpectra.jpg',format='jpg')

# plot DEPTH (see on top parameters)
if PLOT_DEPTH:
    intDtof=zeros(shape(meanDtof))
    figIntDtof=figure(num='FigIntDtof',figsize=cm2inch(0.4*FIGWIDTH,0.4*FIGWIDTH))
    for ic in arange(NUMCHAN):
        intDtof[ic,:]=sum(meanDtof[ic:,:],axis=0)
    semilogy(mtime,intDtof)
    legend(Data.Subject.unique())
    xlabel('time (ps)')
    ylabel('counts')
    xlim([0,MAX_TIME])
    grid()
    
    maxTimeChan=argmin((intDtof>MIN_COUNT_RATE),axis=0)
    maxTimeTime=mtime[maxTimeChan]

    table=pSpectra.pivot_table(values=Opt,index='Lambda',columns='Subject',aggfunc='mean')
    MusDepth=table['Mus'].loc[LAMBDA_DEPTH]
    MuaDepth=table['Mua'].loc[LAMBDA_DEPTH]
    Zmax=zeros(shape(meanDtof))
    for iss,oss in enumerate(Data.Subject.unique()):
        Zmax[:,iss]=2*aDepth*(mtime*MUS0_DEPTH/MusDepth[oss])**bDepth
    # Zmax=zeros(shape(meanDtof))
    # for iss,oss in enumerate(Data.Subject.unique()):
    #     Zmax[:,iss]=2*aDepth*(mtime)**bDepth
    figZmax=figure(num='FigZmax',figsize=cm2inch(0.4*FIGWIDTH,0.4*FIGWIDTH))
    gca().set_prop_cycle(None)
    plot(mtime,Zmax,linewidth=2)
    # plot(mtime,Zmax,linewidth=4,color='gray')
    
    maxTimeZ=zeros(shape(maxTimeChan))
    # for iss in arange(ns):
    for iss,oss in enumerate(Data.Subject.unique()):
        maxTimeZ[iss]=Zmax[maxTimeChan[iss],iss]
    xx=stack((maxTimeTime,maxTimeTime))
    yy=stack((zeros(shape(maxTimeZ)),maxTimeZ))
    gca().set_prop_cycle(None)
    plot(xx,yy)
    xx=stack((zeros(shape(maxTimeTime)),maxTimeTime))
    yy=stack((maxTimeZ,maxTimeZ))
    gca().set_prop_cycle(None)
    plot(xx,yy)

    legend(Data.Subject.unique())
    xlabel('time (ps)')
    ylabel('Zmax (mm)')
    xlim([0,MAX_TIME])
    ylim([0,60])
    grid()

    if SAVEFIG:
        figZmax.savefig(PATHBETA+PATHANALYSIS+'FigZmax.eps',format='eps')
    if (FIG_PW2021 or FIG_PAPER):
        figZmax.savefig(PATHBETA+PATHANALYSIS+'FigZmax.jpg',format='jpg')
           

#%% CALC COMPONENTS
Components=read_table(PATHBETA+PATHANALYSIS+FileComponents)
figure(num='FigComp')
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
figure(num='FigConc')
plot(x,y)
#filtData=merge(filtData,dfComp,on=['Subject','Meas','Rho'])

show()