import numpy as np
import pandas as pd
from uncertainties import unumpy

# The two calibration files for each filter

B = pd.read_csv("Bcalib.csv")
V = pd.read_csv("Vcalib.csv")

# the headings are: 'ID', 'RA', 'Dec', 'm_inst', 'm_inst_err', 'm_sim', 'm_sim_err'

# Creating uncertianty arrays for the instrament magntiude

B_i = unumpy.uarray([B['m_inst']],[B['m_inst_err']])
V_i = unumpy.uarray([V['m_inst']],[V['m_inst_err']])

# Doing the same as above for simbad magntiude

B_s = unumpy.uarray([B['m_sim']],[B['m_sim_err']])
V_s = unumpy.uarray([V['m_sim']],[V['m_sim_err']])

# Subtracting the simbad magnitude form the instriment value to obtain the zero point magnitude

B_z = B_s - B_i
V_z = V_s - V_i

# Taking the mean value for the zeropoint magnitude

print(f"B_z = {np.mean(B_z)} \nV_z = {np.mean(V_z)} \n")

# Creating a table of results of the calibration measurments

print(B.to_latex(index=False, float_format="{:.2f}".format))
print(V.to_latex(index=False, float_format="{:.2f}".format))
