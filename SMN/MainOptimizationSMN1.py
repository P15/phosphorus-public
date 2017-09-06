#%%
def BatchAnalysis(Repeates, N_Block, threshold_int, Ratio, AddedPaths):

    
    # NTC Analysis
    if os.path.isfile(AddedPaths + 'SMN Template_D01_Amplitude.csv'):
        csv = np.genfromtxt (AddedPaths + 'SMN Template_D01_Amplitude.csv',delimiter=",")
    elif os.path.isfile(AddedPaths + 'SMN Validation2_D01_Amplitude.csv'):
        csv = np.genfromtxt (AddedPaths + 'SMN Validation2_D01_Amplitude.csv',delimiter=",")
    elif os.path.isfile(AddedPaths + 'SMN Validation3_D01_Amplitude.csv'):
        csv = np.genfromtxt (AddedPaths + 'SMN Validation3_D01_Amplitude.csv',delimiter=",")
    elif os.path.isfile(AddedPaths + 'SMN Validation4_D01_Amplitude.csv'):
        csv = np.genfromtxt (AddedPaths + 'SMN Validation4_D01_Amplitude.csv',delimiter=",")
    
    Ch1 = csv[1:len(csv),0]
    Ch2 = csv[1:len(csv),1]
    
    Ch1_Ch2 = [Ch1, Ch2]
    
    bandwidth = 0.5
    Channels = 2
    Threshold = np.zeros((Channels, Repeates))
    EstimatedRobustMode = np.zeros(Channels)
    for channel in range(0, Channels):
        x = Ch1_Ch2[channel]
     ## Estimate the robust mode          
        M = modeest.mlv(x, method = "hsm", bw = bandwidth)
        M = M[0]
        EstimatedRobustMode[channel] = M[0]
        x = x - EstimatedRobustMode[channel]        
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
                
                
    Mean_Threshold = np.nanmean(Threshold, axis =1)
    
    # Sample Analysis
    AllFiles = glob.glob(AddedPaths + "*.csv")
    AllData = []
    for file_ in AllFiles:
        df = pd.read_csv(file_,index_col=None, header=0)
        AllData.append(df)
    
    N_AllData = len(AllData)
    
    Raw_Threshold = Mean_Threshold + EstimatedRobustMode
    Final_Threshold = np.zeros((N_AllData,Channels))
    Negative = np.zeros((N_AllData,Channels))
    Positive = np.zeros((N_AllData,Channels))
    FileName = []
    FileIndex = -1
    for CurrentData in AllData:
        FileIndex = FileIndex + 1
        csv = CurrentData.as_matrix()
        Ch1 = csv[0:len(csv),0]
        Ch2 = csv[0:len(csv),1]
        
        File = AllFiles[FileIndex]
        NameIndex = int(File.find('Amplitude'))
        FileName.append(File[NameIndex-4:NameIndex-1])
    
        Ch1_Ch2 = [Ch1, Ch2]  
        for channel in range(0, Channels):
            x = Ch1_Ch2[channel]
            Raw_Negative = np.asarray([value for index,value in enumerate(x)\
            if value <= Raw_Threshold[channel]])
            M = modeest.mlv(Raw_Negative, method = "hsm", bw = bandwidth)
            M = M[0]
            Raw_Negative_EstimatedRobustMode = M[0]
            Final_Threshold[FileIndex,channel] = Mean_Threshold[channel] + \
            Raw_Negative_EstimatedRobustMode
            
            Negative[FileIndex,channel] = len(np.asarray([value for index,value in \
            enumerate(x) if value <= Final_Threshold[FileIndex,channel]]))
                
            Positive[FileIndex,channel] = len(np.asarray([value for index,value in \
            enumerate(x) if value > Final_Threshold[FileIndex,channel]]))            
     
         
    R = np.zeros((N_AllData))
    for i in range(0, N_AllData):
        Negative1 = Negative[i,0]
        Positive1 = Positive[i,0]
    
        Negative2 = Negative[i,1]
        Positive2 = Positive[i,1]
        
        C1 = -math.log(Negative1/float(Negative1+Positive1))    
        C2 = -math.log(Negative2/float(Negative2+Positive2))
        
        if C1 == 0 or C2 == 0:
            R[i] = np.NaN
        else:
            R[i] = Ratio*C1/float(C2)
                
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
    
    DataPath = ['Work/Phosphorus/Simulations/Data/ddPCR_Threshold/SMN1/Validation/Complete/Validation_1/', \
            'Work/Phosphorus/Simulations/Data/ddPCR_Threshold/SMN1/Validation/Complete/Validation_2/',\
            'Work/Phosphorus/Simulations/Data/ddPCR_Threshold/SMN1/Validation/Complete/Validation_3/',\
            'Work/Phosphorus/Simulations/Data/ddPCR_Threshold/SMN1/Validation/Complete/Validation_4/']
    
    Validation = range(0, 4)
    
    Results_Method = []
    for v in Validation:
        AddedPaths = FilesMainPath + DataPath[v]
        
        if OperatingSystem == 'Darwin':    
            AddedPaths = AddedPaths.replace('\\','/')
        elif OperatingSystem == 'Windows':    
            AddedPaths = AddedPaths.replace('/','\\')
            
        R = BatchAnalysis(Repeates, N_Block, threshold_int, Ratio, AddedPaths)
        R = np.concatenate((R[1:36],R[37:96]))
        if v == 0:
            Results_Method = R
        else:
            Results_Method = np.concatenate((Results_Method,R))

    ReferenceFilePath = FilesMainPath + \
    'Work/Phosphorus/Simulations/Data/ddPCR_Threshold/SMN1/Validation/Reference/'
    if OperatingSystem == 'Darwin':    
            ReferenceFilePath = ReferenceFilePath.replace('\\','/')
    elif OperatingSystem == 'Windows':    
            ReferenceFilePath = ReferenceFilePath.replace('/','\\')
            
    csv = np.genfromtxt (ReferenceFilePath + 'SMN Validation Total Cleaned.csv',delimiter=",")
    
    Results_True = csv[1:len(csv),3]
    
    Results_Method_Round = np.around(Results_Method)
    
    Results_Error_Method = np.nansum(np.absolute(Results_True-Results_Method_Round)) \
        /float(np.count_nonzero(~np.isnan(Results_True-Results_Method_Round)))

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
import glob
import pandas as pd
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

#cons = ({'type': 'ineq', 'fun' : lambda x: np.array([1 - x[2]])})
#          
#res = minimize(RunOnAllValidations, x0, constraints=cons, \
#     method='COBYLA', options={'disp': True})      
 
Methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG','L-BFGS-B', 'TNC',\
    'COBYLA', 'SLSQP', 'dogleg', 'trust-ncg']

OptimizationMethod = Methods[10] 

print(OptimizationMethod)

res = minimize(RunOnAllValidations, x0, method=OptimizationMethod, options={'disp': True})    
     
print(OptimizationMethod)
print(res)












    
    
     










