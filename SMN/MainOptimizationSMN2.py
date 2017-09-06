#%%
def NTC_Processing(Repeates, N_Block, threshold_int, NTC_List):

    # NTC Analysis
    bandwidth = 0.5
    Channels = 2
    
    N_NTC = len(NTC_List)
    EstimatedRobustMode = np.zeros((N_NTC, Channels))
    
    for n_NTC in range(0, N_NTC):
        
        csv = np.genfromtxt (NTC_List[n_NTC], delimiter=",")
       
        Ch1 = csv[1:len(csv),0]
        Ch2 = csv[1:len(csv),1]
        
        Ch1_Ch2 = [Ch1, Ch2]
        for channel in range(0, Channels):
            x = Ch1_Ch2[channel]
         ## Estimate the robust mode          
            M = modeest.mlv(x, method = "hsm", bw = bandwidth)
            M = M[0]
            EstimatedRobustMode[n_NTC, channel] = M[0]
            x = x - EstimatedRobustMode[n_NTC, channel]   
            
            if n_NTC == 0:
                if channel == 0:
                    Channel1 = x
                elif channel == 1:
                    Channel2 = x
            else:
                if channel == 0:
                    Channel1 = np.concatenate((Channel1,x))
                elif channel == 1:
                    Channel2 = np.concatenate((Channel2,x))
                

    Ch1_Ch2 = [Channel1, Channel2]    
    N_Block = N_Block*N_NTC
    Threshold = np.zeros((Channels, Repeates))
    
    for channel in range(0, Channels):     
        x = Ch1_Ch2[channel]
        L_Block = math.floor(len(x)/float(N_Block))
        
        r = -1
        while r < Repeates - 1:
            r = r + 1
            ## Find maxima
            np.random.shuffle(x)
            Maxima = []
            for i in range(0, N_Block):
                Maxima.append(max(x[int(i*L_Block):int((i+1)*L_Block-1)]))
                
            Maxima = np.asarray(Maxima)
            
            ## Fit GEV
            try:                
                droplet_fit = evd.fgev(Maxima, method="L-BFGS-B")
                droplet_fit_str = str(droplet_fit)
                droplet_fit_str =  str.split(droplet_fit_str)
                Index = droplet_fit_str.index('Estimates')
                loc = float(droplet_fit_str[Index+4])
                scale = float(droplet_fit_str[Index+5])
                shape = float(droplet_fit_str[Index+6])
                
                cutoff_quantile = threshold_int
                
                quantgev = evd.qgev(cutoff_quantile,loc,scale,shape)
                Threshold[channel,r] = float(quantgev[0])
            except:
                r = r - 1
                
                
    Mean_Threshold = np.nanmean(Threshold, axis = 1)
    Mean_EstimatedRobustMode = np.mean(EstimatedRobustMode, axis = 0)
    
    return Mean_Threshold, Mean_EstimatedRobustMode

#%%    
def Sample_Processing(Mean_Threshold, Mean_EstimatedRobustMode, Ratio, Sample_Name):      
    # Sample Analysis
   
    Raw_Threshold = Mean_Threshold + Mean_EstimatedRobustMode
    csv = np.genfromtxt (Sample_Name,delimiter=",")    
    
    bandwidth = 0.5    
    Channels = 2
    
    Final_Threshold = np.zeros((Channels))
    Negative = np.zeros((Channels))
    Positive = np.zeros((Channels))

    Ch1 = csv[0:len(csv),0]
    Ch2 = csv[0:len(csv),1]    
    
    Ch1_Ch2 = [Ch1, Ch2]  
    for channel in range(0, Channels):
        x = Ch1_Ch2[channel]
        Raw_Negative = np.asarray([value for index,value in enumerate(x)\
            if value <= Raw_Threshold[channel]])
        M = modeest.mlv(Raw_Negative, method = "hsm", bw = bandwidth)
        M = M[0]
        Raw_Negative_EstimatedRobustMode = M[0]

        Final_Threshold[channel] = Mean_Threshold[channel] + \
            Raw_Negative_EstimatedRobustMode
        
        Negative[channel] = len(np.asarray([value for index,value in \
            enumerate(x) if value <= Final_Threshold[channel]]))
            
        Positive[channel] = len(np.asarray([value for index,value in \
            enumerate(x) if value > Final_Threshold[channel]]))            
 
         
    Negative1 = Negative[0]
    Positive1 = Positive[0]

    Negative2 = Negative[1]
    Positive2 = Positive[1]
    
    C1 = -math.log(Negative1/float(Negative1+Positive1))    
    C2 = -math.log(Negative2/float(Negative2+Positive2))
    
    if C1 == 0 or C2 == 0:
        R = np.NaN
    else:
        R = Ratio*C1/float(C2)            
    
    return R


#%%    
def RunOnAllValidations(x):
    
    Repeates = x[0] 
    N_Block = x[1]
    threshold_int = x[2] 
    Ratio = x[3]
  
    Repeates = int(np.around(Repeates))
    N_Block = int(np.around(N_Block))
    if threshold_int >= 1:        
        threshold_int = 1 - np.finfo(float).eps
    
      
    OperatingSystem = platform.system()
    if OperatingSystem == 'Darwin':    
        FilesMainPath = '"/Mohammad/"'
    elif OperatingSystem == 'Windows':    
        FilesMainPath = '"C:\Mohammad\\"'
    FilesMainPath =  FilesMainPath[1:len(FilesMainPath)-1]
    
    DataPath = \
        'Work/Phosphorus/Simulations/Data/ddPCR_Threshold/SMN2/Validation/Complete/'
    
    AddedPaths = FilesMainPath + DataPath
    
    if OperatingSystem == 'Darwin':    
        AddedPaths = AddedPaths.replace('\\','/')
    elif OperatingSystem == 'Windows':    
        AddedPaths = AddedPaths.replace('/','\\')
    
    if os.path.isfile(AddedPaths + '.DS_Store'):
        os.remove(AddedPaths + '.DS_Store')
    
    All_Sample_Results = []
    for folder in os.listdir(AddedPaths):
        Path = AddedPaths + folder + '/'
        
        if OperatingSystem == 'Darwin':    
            Path = Path.replace('\\','/')
        elif OperatingSystem == 'Windows':    
            Path = Path.replace('/','\\')
        
        NTC_List = []
        Sample_List = []
        for name in os.listdir(Path):
            if 'A03' in name or 'D01' in name and '.csv' in name:
    #        if 'D01' in name and '.csv' in name:    
                NTC_List.append(Path + name)
            elif '.csv' in name: # To exclude NTC files from the analyzed wells
                Sample_List.append(name)            
    #        if '.csv' in name:        
    #            Sample_List.append(name)
            
        Mean_Threshold, Mean_EstimatedRobustMode = NTC_Processing(Repeates, N_Block, threshold_int, NTC_List)    
        for sample in Sample_List:
            Sample_Results = Sample_Processing(Mean_Threshold, Mean_EstimatedRobustMode, Ratio, Path + sample)
            All_Sample_Results.append(Sample_Results)
     
    Results_Method_Round = np.around(np.asarray(All_Sample_Results)) 
    
    ReferenceFilePath = FilesMainPath + \
    'Work/Phosphorus/Simulations/Data/ddPCR_Threshold/SMN2/Validation/Reference/'
    if OperatingSystem == 'Darwin':    
            ReferenceFilePath = ReferenceFilePath.replace('\\','/')
    elif OperatingSystem == 'Windows':    
            ReferenceFilePath = ReferenceFilePath.replace('/','\\')
            
    csv = np.genfromtxt (ReferenceFilePath + 'AllValidationsCleanedSMN2.csv',delimiter=",")
    
    Results_True = csv[1:len(csv),1] 
        
    Results_Error_Method = np.count_nonzero(Results_Method_Round - Results_True)/float(len(Results_True))
    
    return Results_Error_Method
        
#%% Use R
import os
import platform
OperatingSystem = platform.system()
if OperatingSystem == 'Darwin':
    os.environ['R_HOME']='/Library/Frameworks/R.framework/Resources'
elif OperatingSystem == 'Windows':
    RVersion = os.listdir('C:\Program Files\R')
    RVersion = RVersion[len(RVersion) - 1]
    os.environ['R_HOME']='C:\Program Files\R' + '\\' + RVersion  
    os.environ['R_USER']='User'
 
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

#%% Package Import
modeest = importr('modeest')
evd = importr('evd')

import numpy as np
import math
from scipy.optimize import minimize
#%%

Repeates = 50
N_Block = 252
threshold_int = 0.99851459
Ratio = 1.7911230  


x0 = np.zeros(4)
x0[0] = Repeates
x0[1] = N_Block
x0[2] = threshold_int
x0[3] = Ratio


#Results_Error_Method = RunOnAllValidations(x0)     

Methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG','L-BFGS-B', 'TNC',\
    'COBYLA', 'SLSQP', 'dogleg', 'trust-ncg']

OptimizationMethod = Methods[6] 

print(OptimizationMethod)

res = minimize(RunOnAllValidations, x0, method=OptimizationMethod, options={'disp': True})    
     
print(OptimizationMethod)
print(res)








