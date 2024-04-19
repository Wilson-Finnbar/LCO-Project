import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# increase the amount of viewable coloums on dataframe
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# 1 is B filter and 2 is V filter
df = pd.read_csv("matched.csv")

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
fig1, ax1 = plt.subplots()
sc1 = ax1.scatter(df['B_V'], df['V'], s = 2, c = df['R']*60, cmap = 'copper')
ax1.invert_yaxis()
ax1.set_xlabel('B-V')
ax1.set_ylabel('V')
ax1.set_xlim(-0.5,2.5)
ax1.set_ylim(22,12)
fig1.colorbar(sc1, label = 'Radius/arcmins')
plt.tight_layout()
#plt.show()

# Mean elongation and ellipticity
df['Elong'] = (df['ELONGATION_1'] + df['ELONGATION_2'])/2
df['Elip'] = (df['ELLIPTICITY_1'] + df['ELLIPTICITY_2'])/2

# plotting RA/Dec with elongtion
fig2, ax2 = plt.subplots()
sc2 = ax2.scatter(df['mRa'], df['mDec'], s = 2, c = df['Elong'], cmap = 'copper_r')
ax2.set_xlabel('RA (deg)')
ax2.set_ylabel('Dec (deg)')
plt.axis('scaled')
ax2.set_ylim(-46.6,-46.25)
ax2.set_xlim(154.2,154.6)
fig2.colorbar(sc2)
plt.tight_layout()
#plt.show()

# plotting RA/Dec with ellipticity
fig3, ax3 = plt.subplots()
sc3 = ax3.scatter(df['mRa'], df['mDec'], s = 2, c = df['Elip'], cmap = 'copper_r')
ax3.set_xlabel('RA (deg)')
ax3.set_ylabel('Dec (deg)')
plt.axis('scaled')
ax3.set_ylim(-46.6,-46.25)
ax3.set_xlim(154.2,154.6)
fig3.colorbar(sc3)
plt.tight_layout()
#plt.show()

# loading isochrone data
iso = pd.read_csv("iso.csv")
print(iso.head())
uniq = iso['logAge'].unique()

fig4, ax4 = plt.subplots()
ax4.plot(df['B_V'], df['V'],'r.', ms = 2)
ax4.invert_yaxis()
for i in range(0,len(uniq),2):
    df2 = iso[iso['logAge'] == uniq[i]]
    ax4.plot(df2['Bmag']-df2['Vmag'],df2['Vmag'],'k.')
ax4.set_xlabel('B-V')
ax4.set_ylabel('V')
#ax4.set_xlim(-0.5,2.5)
#ax4.set_ylim(22,12)
plt.tight_layout()
plt.show()


