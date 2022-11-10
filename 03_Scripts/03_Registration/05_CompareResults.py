import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

DataPath = Path.cwd() / '../../04_Results/03_Registration'
B = pd.read_csv(str(DataPath / 'Benjamin.csv'))
M = pd.read_csv(str(DataPath / 'Results.csv'))

Figure, Axis = plt.subplots(1,1)
Axis.plot(B['Sample'],B['Dice'], marker='o', fillstyle='none', linestyle='none', color=(0,0,1), label='Benjamin')
Axis.plot(M['Sample'],M['Dice'], marker='o', fillstyle='none', linestyle='none', color=(1,0,0), label='Mathieu')
Axis.set_xticks(np.arange(len(M)))
Axis.set_xticklabels(M['Sample'], rotation=90)
plt.legend()
plt.show()