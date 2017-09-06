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
        
        if C1 ==0 or C2 == 0:
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
    Results_CNV = csv[1:len(csv),2]
    
    Results_Method_Round = np.around(Results_Method)
    Results_CNV_Round = np.around(Results_CNV)
    
    Results_Difference_Method = np.nansum(np.absolute(Results_True-Results_Method_Round))
    Results_Difference_CNV = np.nansum(np.absolute(Results_True-Results_CNV_Round)) 
   
    
    Results_Error_Method = np.nansum(np.absolute(Results_True-Results_Method_Round)) \
        /float(np.count_nonzero(~np.isnan(Results_True-Results_Method_Round)))
    Results_Error_CNV = np.nansum(np.absolute(Results_True-Results_CNV_Round)) \
        /float(np.count_nonzero(~np.isnan(Results_True-Results_CNV_Round)))

    return Results_Error_Method, Results_Error_CNV, Results_Difference_Method, \
        Results_Difference_CNV, Results_Method, Results_Method_Round,\
        Results_CNV, Results_CNV_Round, Results_True, csv


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

Results_Error_Method, Results_Error_CNV, Results_Difference_Method, \
        Results_Difference_CNV, Results_Method, Results_Method_Round,\
        Results_CNV, Results_CNV_Round, Results_True, csv =  RunOnAllValidations(x0)

Error_Reduction_Percent = 100 - Results_Error_Method/float(Results_Error_CNV)*100
print Error_Reduction_Percent    

#Repeates = 50
#N_Block = 250
#threshold_int = 0.99851459
#Ratio = 1.7911230 
## Error_Reduction_Percent = 11.4516129032   
    
#Repeates = 50
#N_Block = 250 
#threshold_int = 0.9991
#Ratio = 1.8   
## Error_Reduction_Percent = 8.17204301075

#%% Accuracy Analysis
Confusion_Matrix_CNV = np.zeros((6,9))
Confusion_Matrix_Method = np.zeros((6,9))
N_Samples = len(Results_True)
for i in range(0, N_Samples):
    for j in range(1, 5):
        
        if Results_True[i] == j:
            
            if Results_CNV_Round[i] == 1:            
                Confusion_Matrix_CNV[j,1] = Confusion_Matrix_CNV[j,1] + 1
            elif Results_CNV_Round[i] == 2:            
                Confusion_Matrix_CNV[j,2] = Confusion_Matrix_CNV[j,2] + 1
            elif Results_CNV_Round[i] == 3:            
                Confusion_Matrix_CNV[j,3] = Confusion_Matrix_CNV[j,3] + 1
            elif Results_CNV_Round[i] == 4:            
                Confusion_Matrix_CNV[j,4] = Confusion_Matrix_CNV[j,4] + 1
            elif Results_CNV_Round[i] == 5:            
                Confusion_Matrix_CNV[j,5] = Confusion_Matrix_CNV[j,5] + 1
            elif Results_CNV_Round[i] == 0:            
                Confusion_Matrix_CNV[j,6] = Confusion_Matrix_CNV[j,6] + 1
            elif np.isnan(Results_CNV_Round[i]):            
                Confusion_Matrix_CNV[j,7] = Confusion_Matrix_CNV[j,7] + 1
             
            if Results_Method_Round[i] == 1:            
                Confusion_Matrix_Method[j,1] = Confusion_Matrix_Method[j,1] + 1
            elif Results_Method_Round[i] == 2:            
                Confusion_Matrix_Method[j,2] = Confusion_Matrix_Method[j,2] + 1
            elif Results_Method_Round[i] == 3:            
                Confusion_Matrix_Method[j,3] = Confusion_Matrix_Method[j,3] + 1
            elif Results_Method_Round[i] == 4:            
                Confusion_Matrix_Method[j,4] = Confusion_Matrix_Method[j,4] + 1
            elif Results_Method_Round[i] == 5:            
                Confusion_Matrix_Method[j,5] = Confusion_Matrix_Method[j,5] + 1
            elif Results_Method_Round[i] == 0:            
                Confusion_Matrix_Method[j,6] = Confusion_Matrix_Method[j,6] + 1
            elif np.isnan(Results_Method_Round[i]):            
                Confusion_Matrix_Method[j,7] = Confusion_Matrix_Method[j,7] + 1          
   

Confusion_Matrix_CNV[5,:] = np.sum(Confusion_Matrix_CNV,axis = 0)
Confusion_Matrix_CNV[:,8] = np.sum(Confusion_Matrix_CNV,axis = 1)

Confusion_Matrix_Method[5,:] = np.sum(Confusion_Matrix_Method,axis = 0)
Confusion_Matrix_Method[:,8] = np.sum(Confusion_Matrix_Method,axis = 1)

Accuracy_Table_CNV = np.zeros((5,6))
Accuracy_Table_Method = np.zeros((5,6))

TP_CNV = np. zeros(5)
TN_CNV = np. zeros(5)
FP_CNV = np. zeros(5)
FN_CNV = np. zeros(5)

TP_Method = np. zeros(5)
TN_Method = np. zeros(5)
FP_Method = np. zeros(5)
FN_Method = np. zeros(5)

for i in range(1, 5):
    TP_CNV[i] = Confusion_Matrix_CNV[i, i]
    TN_CNV[i] = (Confusion_Matrix_CNV[5, 8] - Confusion_Matrix_CNV[5, i]) \
        - (Confusion_Matrix_CNV[i, 8] - Confusion_Matrix_CNV[i, i])
    FP_CNV[i] = Confusion_Matrix_CNV[5, i]  -  Confusion_Matrix_CNV[i, i]
    FN_CNV[i] = (Confusion_Matrix_CNV[i, 8]) - Confusion_Matrix_CNV[i, i]
    
    TP_Method[i] = Confusion_Matrix_Method[i, i]
    TN_Method[i] = (Confusion_Matrix_Method[5, 8] - Confusion_Matrix_Method[5, i]) \
        - (Confusion_Matrix_Method[i, 8] - Confusion_Matrix_Method[i, i])
    FP_Method[i] = Confusion_Matrix_Method[5, i]  -  Confusion_Matrix_Method[i, i]
    FN_Method[i] = (Confusion_Matrix_Method[i, 8]) - Confusion_Matrix_Method[i, i]
    
    # Sensitivity    
    Accuracy_Table_CNV[i, 1] = TP_CNV[i]/float(TP_CNV[i]+FN_CNV[i])
    # Specificity
    Accuracy_Table_CNV[i, 2] = TN_CNV[i]/float(TN_CNV[i]+FP_CNV[i])
    # PPV
    Accuracy_Table_CNV[i, 3] = TP_CNV[i]/float(TP_CNV[i]+FP_CNV[i])
    # NPV
    Accuracy_Table_CNV[i, 4] = TN_CNV[i]/float(TN_CNV[i]+FN_CNV[i])
    # Accuracy
    Accuracy_Table_CNV[i, 5] = (TP_CNV[i]+TN_CNV[i])\
        /float(TP_CNV[i]+TN_CNV[i]+FP_CNV[i]+FN_CNV[i])

    # Sensitivity    
    Accuracy_Table_Method[i, 1] = TP_Method[i]/float(TP_Method[i]+FN_Method[i])
    # Specificity
    Accuracy_Table_Method[i, 2] = TN_Method[i]/float(TN_Method[i]+FP_Method[i])
    # PPV
    Accuracy_Table_Method[i, 3] = TP_Method[i]/float(TP_Method[i]+FP_Method[i])
    # NPV
    Accuracy_Table_Method[i, 4] = TN_Method[i]/float(TN_Method[i]+FN_Method[i])
    # Accuracy
    Accuracy_Table_Method[i, 5] = (TP_Method[i]+TN_Method[i])\
        /float(TP_Method[i]+TN_Method[i]+FP_Method[i]+FN_Method[i])

Overall_Accuracy_CNV = (Confusion_Matrix_CNV[1,1] + Confusion_Matrix_CNV[2,2] \
    + Confusion_Matrix_CNV[3,3] + Confusion_Matrix_CNV[4,4])\
    /float(Confusion_Matrix_CNV[5,8])
Overall_Accuracy_Method = (Confusion_Matrix_Method[1,1] + Confusion_Matrix_Method[2,2] \
    + Confusion_Matrix_Method[3,3] + Confusion_Matrix_Method[4,4])\
    /float(Confusion_Matrix_Method[5,8])


Accuracy_Table_CNV = 100 * Accuracy_Table_CNV
Accuracy_Table_Method = 100 * Accuracy_Table_Method

Overall_Accuracy_CNV = 100 * Overall_Accuracy_CNV
Overall_Accuracy_Method = 100 * Overall_Accuracy_Method
 


     













    
    
     










