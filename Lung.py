# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 11:56:40 2020

@author: anton
"""
from numpy import *
from scipy import *
from matplotlib.pyplot import *
from pandas import *

# OPTIONS
PLOT_TYPE1=False
PLOT_TYPE2=False
Opt=['Mua','Mus']
YLABEL={'Mua':'absorption (cm-1)','Mus':'reduced scattering (cm-1)'}
XLABEL='time (s)'

# PARAMETERS
PATHBETA='C:\\OneDrivePolimi\\OneDrive - Politecnico di Milano\\Beta\\'
PATHDATA='Data\\Polmone\\'
PATHANALYSIS='Analysis\\Polmone\\FIT\\'
EXTDATA='.dat'
FILEFIT='POLm0080new.txt'
FILEKEY='keyPOLm0080.txt'
PROT10=1
PROT5=2
PROT10_PERIOD=10
PROT10_BASE=10
PROT5_PERIOD=5
NUMGATE=16
WIDTHGATE=500 # ps


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


# CONVERSION FUNCTION
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

# LOAD DICT
dataKey=read_csv(PATHBETA+PATHANALYSIS+FILEKEY,sep='\t')
dcKey=dict(zip(dataKey.Key, dataKey.Value))

# LOAD FIT DATA
Data=read_csv(PATHBETA+PATHANALYSIS+FILEFIT,sep='\t')
Data.rename(columns=dcKey,inplace=True)

# LOAD RAW DATA FOR GATE
for ig in arange(NUMGATE):
    Data['Gate'+str(ig)]=0
for od in Data.Data.unique():
    item=Data[Data.Data==od].iloc[0,:]
    
    # load system
    fid = open(PATHBETA+PATHDATA+item.Syst+EXTDATA,'rb')
    fid.seek(HEADLEN,0) # 0=beginning
    fid.seek(SUBLEN,1) # 1=current
    Syst=fromfile(fid,dtype=uint16,count=NUMCHAN)
    fid.close()
    
    # load data
    Dtof=zeros([item.NumRep,item.NumClock,NUMCHAN])
    fid = open(PATHBETA+PATHDATA+item.Data+EXTDATA,'rb')
    fid.seek(HEADLEN,0) # 0=beginning
    for ir in arange(item.NumRep):
        for ik in arange(item.NumClock):
            fid.seek(SUBLEN,1) # 1=current
            Dtof[ir,ik,:]=fromfile(fid,dtype=uint16,count=NUMCHAN)
            
    # Calc peak
    Chan0=Syst.tolist().index(max(Syst))
    mtime=(arange(NUMCHAN)-Chan0)*ch2ps
    

# % calculate
# DataIn=squeeze(mean(mean(Data(FirstRep:end,InClock,:),2),1));
# DataOut=squeeze(mean(mean(Data(FirstRep:end,OutClock,:),2),1));
# DataClock=squeeze(mean(Data(FirstRep:end,:,:),1));

    # Calc gate
    Gate=zeros([item.NumRep,item.NumClock,NUMGATE])
    StopGate=(Chan0+arange(NUMGATE)*(WIDTHGATE/ch2ps)).astype(int)
    StartGate=StopGate[arange(NUMGATE)-1]+1 # stop ends exacly before the start
    StartGate=insert(StartGate,0,Chan0) # add first gate;
    for ig in arange(NUMGATE):
        Gate[:,:,ig]=sum(Dtof[:,:,StartGate[ig]:StopGate[ig]],axis=2)
        
    Data[Data.Data==od].Gate0=Gate[:,:,0].flatten()
        
# TimeGate=(0:NumGate-1)*WidthGate;
# for ir=1:NumRep,
#     Ref=mean(Gate(ir,:,:),2);
#     Ref(Ref<MINREF)=0;
#     for ik=1:NumClock,
#         DeltaMua(ir,ik,:)=-log(Gate(ir,ik,:)./Ref);
#         %if(Ref<1) DeltaMua(ir,ik,:)=0; end
#     end
# end

# %bad=abs(DeltaMua)>1;
# %DeltaMua=DeltaMua.*(1-bad);

# for ig=1:NumGate,
#     temp=(squeeze(DeltaMua(:,:,ig)))';
#     MuaTimeLab(:,ig)=temp(:);
# end
# TimeLab=(1:NumRep*NumClock)*1;

# % plot time course
# for ig=1:NumGate,
#     subplot(4,4,ig);
#     plot(TimeLab,squeeze(MuaTimeLab(:,ig)));
#     xlabel('laboratory time (s)'), ylabel('Delta Mua * l');
#     title([num2str(ig*WidthGate) ' ps']);
#     ax = gca; ax.XTick=[0:NumClock:NumClock*NumRep]; grid on;
# end

# suptitle(LabelMeas);

 
# % PLOT MUA vs GATE

# figure('Name','PLOT MUA vs GATE','Position',[0 0 scrsz(3) 0.83*scrsz(4)]);


# % calculate
# DeltaMuaIn=squeeze(mean(mean(DeltaMua(FirstRep:end,InClock,:),2),1));
# DeltaMuaOut=squeeze(mean(mean(DeltaMua(FirstRep:end,OutClock,:),2),1));
# DeltaMuaOutIn=DeltaMuaOut-DeltaMuaIn;
# DeltaMuaClock=squeeze(mean(DeltaMua(FirstRep:end,:,:),1));

# % plot mean phase
# subplot(1,2,1);
# plot(TimeGate,[DeltaMuaIn DeltaMuaOut DeltaMuaOutIn]);
# legend('IN','OUT', 'OUT/IN'), xlabel('time (ps)'), ylabel('Delta Mua * l'), xlim([TNEG,TPOS]);
# title('mean Phase');

# % plot mean clock
# subplot(1,2,2);
# plot(TimeGate,DeltaMuaClock);
# legend(num2str((1:NumClock)')), xlabel('time (ps)'), ylabel('Delta Mua * l'), xlim([TNEG,TPOS]);
# title('mean Clock');

# suptitle(LabelMeas);


# %% END CYCLE
# %close all;
# end

    

# FILT DATA

# FOLDING AVERAGE
Data['RefTime']=0
Data.RefTime=(Data.Time)%Data.NumClock

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
            
# PLOT TYPE2
Color=['purple','blue']
Linestyle=['-','--']
if PLOT_TYPE2:
    for od in Data.Detector.unique():
        for os in Data[Data.Detector==od].Subject.unique():
            fig=figure(figsize=cm2inch(40, 15))
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
            show()