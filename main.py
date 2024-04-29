import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from uncertainties import unumpy
from astropy.io import fits
from astropy.wcs import WCS

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
df = df[df['R'] > 1.5/60]

df['V'] = df['V'] - 5*np.log10(d/10)
df['B'] = df['B'] - 5*np.log10(d/10)
df['B_V'] = df['B'] - df['V']

# plotting a colour magntiude diagram
fig1, ax1 = plt.subplots(figsize=(7, 4.5))
sc1 = ax1.scatter(df['mRa'], df['mDec'], s = 2, c = df['R']*60, cmap = 'copper')
ax1.invert_yaxis()
ax1.set_xlabel('RA')
ax1.set_ylabel('Dec')
#ax1.set_xlim(-0.5,2.5)
#ax1.set_ylim(8,-2)
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
fig4, ax4 = plt.subplots(figsize=(4, 3))
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
        ax4.plot(pr['Bmag']-pr['Vmag'],pr['Vmag'], color = cid[a-1], alpha=0.8)
        print(f"{uniqmh[j]} : {cid[a-1]}")
ax4.set_xlabel('B-V')
ax4.set_ylabel('V')
ax4.set_xlim(-0.5,2.5)
ax4.set_ylim(8,-2)
plt.tight_layout()
plt.savefig("Figures/manyiso.pdf")

# plotting the isochrone selected 
fig5, ax5 = plt.subplots(figsize=(4, 3))
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
ax5.annotate('[M/H] = −0.4\n12.0 Gyr',(1.5,0))
ax5.set_xlabel('B-V')
ax5.set_ylabel('V')
ax5.set_xlim(-0.5,2.5)
ax5.set_ylim(8,-2)
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
ax6.annotate('[M/H] = −0.4\n12.0 Gyr',(1.5,0))
ax6.set_xlabel('B-V')
ax6.set_ylabel('V')
ax6.set_xlim(-0.5,2.5)
ax6.set_ylim(8,-2)
plt.tight_layout()
plt.savefig("Figures/manual.pdf")

pop1 = df[df['R'] < 0.09]
pop2 = df[df['R'] > 0.09]

Berrors1 = unumpy.uarray(pop1['B'],pop1['MAGERR_ISO_1'])
Verrors1 = unumpy.uarray(pop1['V'],pop1['MAGERR_ISO_2'])
B_Verr1 = Berrors1 - Verrors1

Berrors2 = unumpy.uarray(pop2['B'],pop2['MAGERR_ISO_1'])
Verrors2 = unumpy.uarray(pop2['V'],pop2['MAGERR_ISO_2'])
B_Verr2 = Berrors2 - Verrors2

fig7,ax7 = plt.subplots(figsize=(4, 3))
ax7.errorbar(pop1['B_V'],pop1['V'], yerr = pop1['MAGERR_ISO_2'], xerr = unumpy.std_devs(B_Verr1), fmt = 'k.',ecolor = 'r', elinewidth = 0.5, capsize = 1.0, markersize=2, zorder = 0)
pr = iso[iso['MH'] == -0.4]
hold = pr[pr['logAge'] == 10.11394]
hold = hold[hold['Vmag']<8]
hold = hold[hold['Vmag']>1.5]
hold = hold[(hold['Bmag']-pr['Vmag'])<1.0]
ax7.plot(hold['Bmag']-hold['Vmag'],hold['Vmag'], color = 'tab:purple', marker='.', ms = 1, lw = 3.0, label = f'M/H = -0.4 \nAge = {10**10.07918 / 1E9 : .2f} Gyr')
ax7.annotate('[M/H] = −0.4\n12.0 Gyr',(1.5,0))
ax7.set_xlabel('B-V')
ax7.set_ylabel('V')
ax7.set_xlim(-0.5,2.5)
ax7.set_ylim(8,-2)
plt.tight_layout()
plt.savefig("Figures/pop1.pdf")

fig8,ax8 = plt.subplots(figsize=(4, 3))
ax8.errorbar(pop2['B_V'],pop2['V'], yerr = pop2['MAGERR_ISO_2'], xerr = unumpy.std_devs(B_Verr2), fmt = 'k.',ecolor = 'r', elinewidth = 0.5, capsize = 1.0, markersize=2, zorder = 0)
pr = iso[iso['MH'] == -0.4]
hold = pr[pr['logAge'] == 10.07918]
hold = hold[hold['Vmag']<8]
hold = hold[hold['Vmag']>1.5]
hold = hold[(hold['Bmag']-pr['Vmag'])<1.0]
ax8.plot(hold['Bmag']-hold['Vmag'],hold['Vmag'], color = 'tab:purple', marker='.', ms = 1, lw = 3.0, label = f'M/H = -0.4 \nAge = {10**10.07918 / 1E9 : .2f} Gyr')
ax8.annotate('[M/H] = −0.4\n13.0 Gyr',(1.5,0))
ax8.set_xlabel('B-V')
ax8.set_ylabel('V')
ax8.set_xlim(-0.5,2.5)
ax8.set_ylim(8,-2)
plt.tight_layout()
plt.savefig("Figures/pop2.pdf")

fig9, ax9 = plt.subplots(figsize=(4, 3))
ax9.errorbar(unumpy.nominal_values(Bm_Vm),Vman['v'], yerr = Vman['verr'], xerr = unumpy.std_devs(Bm_Vm), fmt = 'k.',ecolor = 'r', elinewidth = 0.5, capsize = 1.0, markersize=2, zorder = 0)
pr = iso[iso['MH'] == -0.4]
hold = pr[pr['logAge'] == 9.77815]
hold = hold[hold['Vmag']<8]
hold = hold[hold['Vmag']>1.5]
hold = hold[(hold['Bmag']-pr['Vmag'])<1.0]
ax9.plot(hold['Bmag']-hold['Vmag'],hold['Vmag'], color = 'tab:purple', marker='.', ms = 1, lw = 3.0, label = f'M/H = -0.4 \nAge = {10**10.07918 / 1E9 : .2f} Gyr')
ax9.annotate('[M/H] = −0.4\n6.0 Gyr',(1.5,0))
ax9.set_xlabel('B-V')
ax9.set_ylabel('V')
ax9.set_xlim(-0.5,2.5)
ax9.set_ylim(8,-2)
plt.tight_layout()
plt.savefig("Figures/popmanual.pdf")

x, y = Bman['x'].to_numpy(), Bman['y'].to_numpy()
ra = np.zeros(len(Bman['x'].to_numpy()))
dec = np.zeros(len(Bman['x'].to_numpy()))

for i in range(len(Bman['x'].to_numpy())):
    f = fits.open('LCO/B/B_stack.fits')
    mywcs = WCS(f[0].header)
    ra[i], dec[i] = mywcs.all_pix2world([[x[i],y[i]]], 0)[0]
    manR = np.sqrt(ra**2 + dec**2)

fig10,ax10 = plt.subplots(figsize=(4, 3))
ax10.plot(ra,dec,'tab:red',marker='.',ms=2,lw=0,label='Core stars')
ax10.plot(pop1['mRa'], pop1['mDec'], 'tab:orange',marker='.',ms=2,lw=0, label='Population\na stars')
ax10.plot(pop2['mRa'], pop2['mDec'], 'tab:blue',marker='.',ms=2,lw=0, label='Population\nb stars')
plt.axis('scaled')
ax10.set_xlabel('Ra')
ax10.set_ylabel('Dec')
ax10.set_xlim(154.2,154.6)
ax10.set_ylim(-46.6,-46.2)
ax10.legend(loc='center left', bbox_to_anchor=(0.9,0.5), fontsize = 'small', markerscale = 3, frameon = True, shadow=True, edgecolor='0.4')
plt.tight_layout()
plt.savefig("Figures/radec.pdf")

dataman = {'RA$^{\circ}$' : ra,
           'Dec$^{\circ}$' : dec,
           'B' : Bman['b'],
           'B$_{err}$' : Bman['berr'],
           'V' : Vman['v'],
           'V$_{err}$' : Vman['verr']}

df20 = pd.DataFrame(dataman)
#print(df20.to_latex(index = False,float_format="%.3f"))

