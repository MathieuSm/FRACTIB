# 00 Initialization
import os
import numpy as np
import pandas as pd

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)


# 01 Set variables
WorkingDirectory = os.getcwd()
DataPath = os.path.join(WorkingDirectory,'04_Results/02_MorphoAnalysis/Binarized_Scans/')
ResultsPath = os.path.join(WorkingDirectory,'04_Results/02_MorphoAnalysis/')

# 02 Collect data
Files = [File for File in os.listdir(DataPath) if File.endswith('.csv')]
MorphoData = pd.DataFrame()

for File in Files:

    MorphoFile = pd.read_csv(DataPath + File,sep=';')

    BVTV = MorphoFile['$BVTV_voxel'].values[0]
    SMI = MorphoFile['$SMI_voxel'].values[0]
    Tb_N_Mean = MorphoFile['$Tb_N_mean'].values[0]

    Tb_Sp_Mean = MorphoFile['$Tb_Sp_mean'].values[0]
    Tb_Sp_Max = MorphoFile['$Tb_Sp_max'].values[0]
    Tb_Sp_Min = MorphoFile['$Tb_Sp_min'].values[0]
    Tb_Sp_Std = MorphoFile['$Tb_Sp_stddev'].values[0]

    Tb_Th_Mean = MorphoFile['$Tb_Th_mean'].values[0]
    Tb_Th_Max = MorphoFile['$Tb_Th_max'].values[0]
    Tb_Th_Min = MorphoFile['$Tb_Th_min'].values[0]
    Tb_Th_Std = MorphoFile['$Tb_Th_stddev'].values[0]

    m1 = MorphoFile['$DA_lam_1'].values[0]
    m2 = MorphoFile['$DA_lam_2'].values[0]
    m3 = MorphoFile['$DA_lam_3'].values[0]
    DA = MorphoFile['$DA_value'].values[0]

    m1x = MorphoFile['$DA_vec_1x'].values[0]
    m1y = MorphoFile['$DA_vec_1y'].values[0]
    m1z = MorphoFile['$DA_vec_1z'].values[0]

    m2x = MorphoFile['$DA_vec_2x'].values[0]
    m2y = MorphoFile['$DA_vec_2y'].values[0]
    m2z = MorphoFile['$DA_vec_2z'].values[0]

    m3x = MorphoFile['$DA_vec_3x'].values[0]
    m3y = MorphoFile['$DA_vec_3y'].values[0]
    m3z = MorphoFile['$DA_vec_3z'].values[0]

    Ct_Th_Mean = MorphoFile['$Ct_Th_mean'].values[0]
    Ct_Th_Max = MorphoFile['$Ct_Th_max'].values[0]
    Ct_Th_Min = MorphoFile['$Ct_Th_min'].values[0]
    Ct_Th_Std = MorphoFile['$Ct_Th_stddev'].values[0]

    PoresN = MorphoFile['$Po_N_value'].values[0]
    OpenPoresN = MorphoFile['$Po_N_value(op)'].values[0]
    ClosedPoresN = MorphoFile['$Po_N_value(cl)'].values[0]

    PoresVolume = MorphoFile['$Po_V_voxel'].values[0]
    TotalVolume = MorphoFile['$TV_voxel'].values[0]


    Parameters = {'Scan':File[:-14],
                  'BV/TV':BVTV,
                  'SMI':SMI,
                  'Mean Tb N': Tb_N_Mean,
                  'Mean Tb Th':Tb_Th_Mean,
                  'Std Tb Th':Tb_Th_Std,
                  'Min Tb Th':Tb_Th_Min,
                  'Max Tb Th':Tb_Th_Max,
                  'Mean Tb Sp':Tb_Sp_Mean,
                  'Std Tb Sp':Tb_Sp_Std,
                  'Min Tb Sp':Tb_Sp_Min,
                  'Max Tb Sp':Tb_Sp_Max,
                  'm1':m1,
                  'm2':m2,
                  'm3':m3,
                  'm11': m1x,
                  'm12': m1y,
                  'm13': m1z,
                  'm21': m2x,
                  'm22': m2y,
                  'm23': m2z,
                  'm31': m3x,
                  'm32': m3y,
                  'm33': m3z,
                  'Degree of Anisotropy':DA,
                  'Mean Ct Th':Ct_Th_Mean,
                  'Std Ct Th':Ct_Th_Std,
                  'Min Ct Th':Ct_Th_Min,
                  'Max Ct Th':Ct_Th_Max,
                  'PN':PoresN,
                  'oPN':OpenPoresN,
                  'cPN':ClosedPoresN,
                  'PV':PoresVolume,
                  'TV':TotalVolume}

    MorphoData = MorphoData.append(Parameters,ignore_index=True)


# 03 Sort and save to csv
MorphoData = MorphoData.sort_values(by=['Scan'],ignore_index=True)
MorphoData.to_csv(ResultsPath + '00_Data.csv',index=False)
