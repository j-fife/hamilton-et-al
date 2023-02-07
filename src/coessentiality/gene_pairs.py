import numpy as np, statsmodels.api as sm
from load_screens import load_screens
import pandas as pd 
# Load batch-corrected screens

screens = load_screens()
genes = screens.index

# Compute p-values and effect directions with GLS
# loading lipid data instead, to go back to normal remove the two lines below: 
screens = pd.read_csv("pos_lfc.csv", index_col = 0)
genes = screens.index

cholsigmainv = np.linalg.cholesky(np.linalg.pinv(screens.cov())).T
GLS_beta = np.empty((len(screens), len(screens)))
GLS_p = np.empty((len(screens), len(screens)))
for A_index, A_row in enumerate(screens.values):
    input = cholsigmainv.dot(sm.add_constant(A_row))
    for B_index, B_row in enumerate(screens.values):
        output = cholsigmainv.dot(B_row)
        model = sm.OLS(output, input).fit()
        GLS_beta[A_index, B_index] = model.params[1]
        GLS_p[A_index, B_index] = model.pvalues[1]
np.fill_diagonal(GLS_p, 1)
GLS_sign = np.sign(GLS_beta)

# Save everything

np.save('GLS_p_lipids.npy', GLS_p)
np.save('GLS_sign_lipids.npy', GLS_sign)
genes.to_series().to_csv('genes_lipids.txt', index=False, header=False)
