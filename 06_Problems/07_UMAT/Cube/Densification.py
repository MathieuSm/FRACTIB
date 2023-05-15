#%% Imports
# Modules import

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%% Main
# Main code

def Main():

    # Read data
    Nodes = pd.read_csv('Nodes.csv')
    Elements = pd.read_csv('Elements.csv')

    # Input file values
    BVTV = 0.1
    MM1 = 1.0
    MM2 = 1.0
    MM3 = 1.0
    PBVT = 1.0

    # UMAT values
    SCA1 = 0.9419
    SCA2 = 0.7801

    # Elasticity parameters (FABRGWTB PMUBC, PHZ, 2013) new version 2017 trabecular bone
    E0  = 9759.0*(12000.0/9759.0)*SCA1

    # Strength parameters (FABRGWTB, PHZ, 2013) trabecular bone
    SIGD0N = 73.1*(12000.0/9759.0)*SCA2
    SIGD0P = 57.69*(12000.0/9759.0)*SCA2
    ZETA0  = 0.28
    TAUD0  = 29.61*(12000.0/9759.0)*SCA2
    PP     = 1.82
    QQ     = 0.98

    # Densification
    ECRIT = -0.3
    GAMMAL0 = 1100.0
    RL = 2.928
    GAMMAP0 = 1300.0
    RP = 2.77
    POWERDNS = 6.0

    # Compute densification
    TRACEEN = Elements['E1']+Elements['E2']+Elements['E3']
    CRITERION = TRACEEN-ECRIT
    SDNST = np.zeros((len(Elements), 3,3))
    XIDEN = np.eye(3)

    GAMMAL = GAMMAL0*(BVTV**RL)
    GAMMAP = GAMMAP0*(BVTV**RP)
    for K1 in range(3):
        for K2 in range(3):
            SDNST[:,K1,K2] = (GAMMAL*CRITERION+GAMMAP* CRITERION**(POWERDNS-1.0))*XIDEN[K1,K2]

    # Elastic trial stress
    STR = np.zeros((6,len(Elements)))
    for i in range(6):
        STR[i] = Elements['S' + str(i+1)]

    # Quadric fourth order tensor FFFF
    S0 = (SIGD0P+SIGD0N)/2.0/SIGD0P/SIGD0N
    FFFFT = np.zeros((6,6))
    FFFFT[0,0] = S0**2/((BVTV**PP)*MM1**(2.0*QQ))**2
    FFFFT[1,1] = S0**2/((BVTV**PP)*MM2**(2.0*QQ))**2
    FFFFT[2,2] = S0**2/((BVTV**PP)*MM3**(2.0*QQ))**2
    FFFFT[0,1] = -(ZETA0*(MM1/MM2)**(2.0*QQ))*S0**2/((BVTV**PP)*MM1**(2.0*QQ))**2
    FFFFT[0,2] = -(ZETA0*(MM1/MM3)**(2.0*QQ))*S0**2/((BVTV**PP)*MM1**(2.0*QQ))**2
    FFFFT[1,2] = -(ZETA0*(MM2/MM3)**(2.0*QQ))*S0**2/((BVTV**PP)*MM2**(2.0*QQ))**2
    FFFFT[1,0] = FFFFT[0,1]
    FFFFT[2,0] = FFFFT[0,2]
    FFFFT[2,1] = FFFFT[1,2]
    FFFFT[3,3] = 0.5/((TAUD0*(BVTV**PP)*(MM1*MM2)**QQ)**2)
    FFFFT[4,4] = 0.5/((TAUD0*(BVTV**PP)*(MM1*MM3)**QQ)**2)
    FFFFT[5,5] = 0.5/((TAUD0*(BVTV**PP)*(MM2*MM3)**QQ)**2)
    FFFF = FFFFT*PBVT

    # Quadratic second order tensor FF
    FFT = np.zeros((3,3))
    FFT[0,0] = -(SIGD0P-SIGD0N)/2.0/SIGD0P/SIGD0N/((BVTV**PP)*MM1**(2.0*QQ))
    FFT[1,1] = -(SIGD0P-SIGD0N)/2.0/SIGD0P/SIGD0N/((BVTV**PP)*MM2**(2.0*QQ))
    FFT[2,2] = -(SIGD0P-SIGD0N)/2.0/SIGD0P/SIGD0N/((BVTV**PP)*MM3**(2.0*QQ))
    FFT  = FFT*PBVT

    FF = np.zeros(6)
    FF[0] = FFT[0,0]
    FF[1] = FFT[1,1]
    FF[2] = FFT[2,2]
    FF[3] = 0.0
    FF[4] = 0.0
    FF[5] = 0.0

    # Yield criterion
    RAD  = 1.0 # Perfect plasticity
    FFS  = np.matmul(FFFF,STR)
    SFFS = np.sum(STR * FFS, axis=0)
    YSTR = np.sqrt(SFFS) + np.sum(FF * STR.T, axis=1) - RAD
    TOL = 1.0E-12
    SLim = Elements['S3'][Elements['S3'][YSTR > TOL].index[0]-1]
    ELim = Elements['E3'][Elements['E3'][YSTR > TOL].index[0]-1]

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(Elements['E3'], Elements['S3'], color=(1,0,0), linewidth=2, label='Abaqus Results')
    Axis.plot([ECRIT, ECRIT], [0, -1.5], color=(0,0,0), linewidth=1, linestyle='--', label='E$_c$')
    Axis.plot(ELim, SLim, color=(0,0,0), marker='x', linestyle='none', label='S$y$')
    Axis.set_xlabel('Strain (-)') # Height of 1 mm -> direct conversion
    Axis.set_ylabel('Stress (MPa)') # Surface of 1 mm2 -> direct conversion
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.2))
    plt.show()

    return

#%% Execution part
# Execution as main
if __name__ == '__main__':
    Main()