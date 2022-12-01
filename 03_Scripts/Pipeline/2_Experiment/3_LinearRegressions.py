# 00 Initialization
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats.distributions import t

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

## Define functions
def PlotRegressionResults(Model,Alpha=0.95):

    print(Model.summary())

    ## Plot results
    Y_Obs = Model.model.endog
    Y_Fit = Model.fittedvalues
    N = int(Model.nobs)
    C = np.matrix(Model.cov_params())
    X = np.matrix(Model.model.exog)
    X_Obs = np.sort(np.array(X[:,1]).reshape(len(X)))


    ## Compute R2 and standard error of the estimate
    E = Y_Obs - Y_Fit
    RSS = np.sum(E ** 2)
    SE = np.sqrt(RSS / Model.df_resid)
    TSS = np.sum((Model.model.endog - Model.model.endog.mean()) ** 2)
    RegSS = TSS - RSS
    R2 = RegSS / TSS
    R2adj = 1 - RSS/TSS * (N-1)/(N-X.shape[1]+1-1)

    ## Compute CI lines
    B_0 = np.sqrt(np.diag(np.abs(X * C * X.T)))
    t_Alpha = t.interval(Alpha, N - X.shape[1] - 1)
    CI_Line_u = Y_Fit + t_Alpha[0] * SE * B_0
    CI_Line_o = Y_Fit + t_Alpha[1] * SE * B_0

    t_Alpha2 = t.interval(0.9, N - X.shape[1] - 1)
    CI_Line_u2 = Y_Fit + t_Alpha2[0] * SE * B_0
    CI_Line_o2 = Y_Fit + t_Alpha2[1] * SE * B_0


    ## Plots
    DPI = 100
    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI, sharey=True, sharex=True)
    Axes.plot(X[:,1], Y_Fit, color=(1,0,0))

    if Model.model.endog_names == 'Stiffness':
        Axes.fill_between(X_Obs, np.sort(CI_Line_o2), np.sort(CI_Line_u2), color=(0, 0, 0), alpha=0.1)
        Axes.plot(X_Obs, np.sort(CI_Line_u), color=(0, 0, 1), linestyle='--')
        Axes.plot(X_Obs, np.sort(CI_Line_o), color=(0, 0, 1), linestyle='--')
        Axes.annotate(r'$N$  : ' + str(N), xy=(0.05, 0.875), xycoords='axes fraction')
        Axes.annotate(r'$R^2$ : ' + format(round(R2, 2), '.2f'), xy=(0.05, 0.8), xycoords='axes fraction')
        Axes.annotate(r'$SE$ : ' + format(round(SE, 2), '.2f'), xy=(0.05, 0.725), xycoords='axes fraction')
        Axes.set_ylabel('Loading Max Stiffness (kN/mm)')

    elif Model.model.endog_names == 'Load':
        Axes.fill_between(X_Obs, np.sort(CI_Line_o2)[::-1], np.sort(CI_Line_u2)[::-1], color=(0, 0, 0), alpha=0.1)
        Axes.plot(X_Obs, np.sort(CI_Line_u)[::-1], color=(0, 0, 1), linestyle='--')
        Axes.plot(X_Obs, np.sort(CI_Line_o)[::-1], color=(0, 0, 1), linestyle='--')
        Axes.annotate(r'$N$  : ' + str(N), xy=(0.05, 0.175), xycoords='axes fraction')
        Axes.annotate(r'$R^2$ : ' + format(round(R2, 2), '.2f'), xy=(0.05, 0.1), xycoords='axes fraction')
        Axes.annotate(r'$SE$ : ' + format(round(SE, 2), '.2f'), xy=(0.05, 0.025), xycoords='axes fraction')
        Axes.set_ylabel('Ultimate Load (kN)')

    Axes.plot(X[:,1], Y_Obs, linestyle='none', marker='o', color=(0,0,0), fillstyle='none')
    Axes.set_xlabel('BMC (HA mg)')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.show()
    plt.close(Figure)


# 01 Set variables
WorkingDirectory = os.getcwd()
DataPath = os.path.join(WorkingDirectory,'04_Results/03_Experiment/')

## Load files
Data = pd.read_csv(DataPath + '00_Data.csv')

## Build data to perform fit
Data2Fit = pd.DataFrame()
Data2Fit['BMC'] = Data['BMC (HA mg)'].values
Data2Fit['Stiffness'] = Data['Loading Max Stiffness (N/mm)'].values/1e3
Data2Fit['Load'] = Data['Ultimate Load (N)'].values/1e3




# 02 Perform linear regression
Stiffness_Fit = smf.ols("Stiffness ~ BMC", data=Data2Fit).fit()
Load_Fit = smf.ols("Load ~ BMC", data=Data2Fit).fit()

PlotRegressionResults(Stiffness_Fit)
PlotRegressionResults(Load_Fit)