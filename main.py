import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from uncertainties import unumpy

# increase the amount of viewable coloums on dataframe
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# 1 is B filter and 2 is V filter
df = pd.read_csv("matched.csv")
print(df.shape[0])

# chaning names of columns to a more readable format
df = df.rename({'MAG_ISO_1' : 'B', 'MAG_ISO_2' : 'V'}, axis = 'columns')

# removing flagged stars
df = df[df['FLAGS_1'] == 0]
df = df[df['FLAGS_2'] == 0]

# NGC 3201 Ra, Dec and radius in degrees
Ra = 154.403416
Dec = -46.412472
d = 5000

# Radius of: core, half-light and radius of cluster
coreR = 0.022
halfR = 0.052
Radius =  0.152

# calculating the mean Ra and Dec of the matched stars
df['mRa'] = (df['X_WORLD_1'] + df['X_WORLD_2'])/2
df['mDec'] = (df['Y_WORLD_1'] + df['Y_WORLD_2'])/2

# calculating the radius from the center of the cluster 
df['R'] = np.sqrt( (df['mRa'] - Ra)**2 + (df['mDec'] - Dec)**2)

# removing all stars that lie beyound the radius of the cluster
df = df[df['R'] < Radius]

df['V'] = df['V'] - 5*np.log10(d/10)
df['B'] = df['B'] - 5*np.log10(d/10)
df['B_V'] = df['B'] - df['V']

# plotting a colour magntiude diagram
fig1, ax1 = plt.subplots(figsize=(7, 4.5))
sc1 = ax1.scatter(df['B_V'], df['V'], s = 2, c = df['R']*60, cmap = 'copper')
ax1.invert_yaxis()
ax1.set_xlabel('B-V')
ax1.set_ylabel('V')
ax1.set_xlim(-0.5,2.5)
ax1.set_ylim(8,-2)
fig1.colorbar(sc1, label = 'Radius/arcmins')
plt.tight_layout()
plt.savefig("Figures/colourmag.pdf")

# error bars
fig2,ax2 = plt.subplots(figsize=(7, 4.5))
Berrors = unumpy.uarray(df['B'],df['MAGERR_ISO_1'])
Verrors = unumpy.uarray(df['V'],df['MAGERR_ISO_2'])
B_Verr = Berrors - Verrors
ax2.errorbar(unumpy.nominal_values(B_Verr),df['V'], yerr = df['MAGERR_ISO_2'], xerr = unumpy.std_devs(B_Verr), fmt = 'k.',ecolor = 'r', elinewidth = 0.5, capsize = 1.0, markersize=2)
ax2.invert_yaxis()
ax2.set_xlabel('B-V')
ax2.set_ylabel('V')
ax2.set_xlim(-0.5,2.5)
ax2.set_ylim(8,-2)
plt.tight_layout()
plt.savefig("Figures/errobar.pdf")


# loading isochrone data
iso = pd.read_csv("iso2.csv")

# all unique ages in the isochrone data
uniqage = iso['logAge'].unique()

# all unique metalicities in the isochrone data
uniqmh = iso['MH'].unique()

# plotting all isochrones and ages onto the colour magnitude diagram
fig4, ax4 = plt.subplots(figsize=(5, 3.5))
ax4.errorbar(unumpy.nominal_values(B_Verr),df['V'], yerr = df['MAGERR_ISO_2'], xerr = unumpy.std_devs(B_Verr), fmt = 'k.',ecolor = 'r', elinewidth = 0.5, capsize = 1.0, markersize=2, zorder=0)
ax4.invert_yaxis()
a = 0
for j in range(0,len(uniqmh),4):
    hold = iso[iso['MH'] == uniqmh[j]]
    a += 1
    for i in range(0,len(uniqage),4):
        cid = ['tab:blue','tab:orange','tab:green','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        pr = hold[hold['logAge'] == uniqage[i]]
        pr = pr[pr['Vmag']<8]
        pr = pr[pr['Vmag']>1.5]
        pr = pr[(pr['Bmag']-pr['Vmag'])<1]
        ax4.plot(pr['Bmag']-pr['Vmag'],pr['Vmag'], color = cid[a-1])
        print(f"{uniqmh[j]} : {cid[a-1]}")
ax4.set_xlabel('B-V')
ax4.set_ylabel('V')
ax4.set_xlim(-0.5,2.5)
ax4.set_ylim(8,-2)
plt.tight_layout()
plt.savefig("Figures/manyiso.pdf")

# plotting the isochrone selected 
fig5, ax5 = plt.subplots(figsize=(5, 3.5))
ax5.errorbar(unumpy.nominal_values(B_Verr),df['V'], yerr = df['MAGERR_ISO_2'], xerr = unumpy.std_devs(B_Verr), fmt = 'k.',ecolor = 'r', elinewidth = 0.5, capsize = 1.0, markersize=2, zorder = 0)
ax5.invert_yaxis()
pr = iso[iso['MH'] == -0.4]
for i in range(len(uniqage)):
    print(uniqage[i])
hold = pr[pr['logAge'] == 10.07918]
hold = hold[hold['Vmag']<8]
hold = hold[hold['Vmag']>1.5]
hold = hold[(hold['Bmag']-pr['Vmag'])<1.0]
ax5.plot(hold['Bmag']-hold['Vmag'],hold['Vmag'], color = 'tab:purple', marker='.', ms = 1, lw = 3.0, label = f'M/H = -0.4 \nAge = {10**10.07918 / 1E9 : .2f} Gyr')
ax5.set_xlabel('B-V')
ax5.set_ylabel('V')
ax5.set_xlim(-0.5,2.5)
ax5.set_ylim(8,-2)
ax5.legend()
plt.tight_layout()
plt.savefig("Figures/singiso.pdf")

# reading in manual recodings
Bman = pd.read_csv("bmanual2.csv")
Vman = pd.read_csv("vmanual2.csv")

# converting to absoulute magnitude
Vman['v'] = Vman['v'] - 5*np.log10(d/10)
Bman['b'] = Bman['b'] - 5*np.log10(d/10)

# obtaining the errors
Bm = unumpy.uarray(Bman['b'],Bman['berr'])
Vm = unumpy.uarray(Vman['v'],Vman['verr'])
Bm_Vm = Bm - Vm

# plotting the manual magntiudes
fig6, ax6 = plt.subplots(figsize=(7, 4.5))
ax6.errorbar(unumpy.nominal_values(Bm_Vm),Vman['v'], yerr = Vman['verr'], xerr = unumpy.std_devs(Bm_Vm), fmt = 'k.',ecolor = 'r', elinewidth = 0.5, capsize = 1.0, markersize=2, zorder = 0)
ax6.plot(hold['Bmag']-hold['Vmag'],hold['Vmag'], color = 'tab:purple', marker='.', ms = 1, lw = 3.0, label = f'M/H = -0.4 \nAge = {10**10.07918 / 1E9 : .2f} Gyr')
ax6.set_xlabel('B-V')
ax6.set_ylabel('V')
ax6.set_xlim(-0.5,2.5)
ax6.set_ylim(8,-2)
plt.tight_layout()
plt.savefig("Figures/manual.pdf")

print(df.shape[0])
