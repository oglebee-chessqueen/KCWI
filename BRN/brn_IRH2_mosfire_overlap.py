#
#
# Plot Eu (K) of H2 IR lines versus emitting wavelength in IR + shaded regions
# covered by Keck-MOSFIRE grating/filter combos
#
#

import matplotlib.pyplot as plt
#~ from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import numpy
plt.rcParams["font.family"] = "serif"
plt.rcParams['axes.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'
ref_dir = 'C:\\Users\\Keri Hoadley\\Documents\\BRN\\'
fname = ref_dir+'H2_transitions_nearIR.dat'



transition = []
wavelength = []
energy = []
aul = []
LyAcascade = []

with open(fname) as fobj:
	for line in fobj:
		if line.startswith('#'):
			continue
		row = line.split()
		transition.append( str(row[0:1]) )
		wavelength.append( float(row[2]) )
		energy.append( float(row[3]) )
		aul.append( float(row[4]) )
		LyAcascade.append( float(row[5]) )

wavelength = numpy.array(wavelength)*1e-4
energy = numpy.array(energy)
aul = numpy.array(aul)#*1e-7
LyAcascade = numpy.array(LyAcascade)


# Plot energy vs. wavelength, plus mosfire filter ranges
mosfire_y = [0.9716, 1.1250]
mosfire_j = [1.1530, 1.3520]
mosfire_h = [1.4680, 1.8040]
mosfire_k = [1.9540, 2.3970]

# Energy vs. wavelength
fig = plt.figure(figsize=(8,8),facecolor='white')
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'H$_2$ Energy Level (K)',weight='bold',size='x-large')
ax1.set_ylabel(r'Wavelength ($\mu$m)',weight='bold',size='x-large')
ax1.set_ylim(0.95, 2.75)
ax1.set_xlim(0, 45000)
ax1.set_title(r"H$_2$ Emission in the Keck-MOSFIRE Filters",weight='bold',size='x-large')
plt.scatter(energy, wavelength, s=8, color='black')

ax1.fill_betweenx(mosfire_y, 0, 45000, alpha=0.3, color='mediumblue', edgecolor='darkblue', label='MOSFIRE Y')
ax1.fill_betweenx(mosfire_j, 0, 45000, alpha=0.3, color='limegreen', edgecolor='darkolivegreen', label='MOSFIRE J')
ax1.fill_betweenx(mosfire_h, 0, 45000, alpha=0.3, color='deeppink', edgecolor='m', label='MOSFIRE H')
ax1.fill_betweenx(mosfire_k, 0, 45000, alpha=0.3, color='red', edgecolor='firebrick', label='MOSFIRE K')

plt.legend(loc='upper right', fontsize='large')


# Aul vs. wavelength
fig = plt.figure(figsize=(8,8),facecolor='white')
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'A$_{ul}$ (10$^{-7}$ s$^{-1}$)',weight='bold',size='x-large')
ax1.set_ylabel(r'Wavelength ($\mu$m)',weight='bold',size='x-large')
ax1.set_ylim(0.95, 2.75)
ax1.set_xlim(0, 15)
ax1.set_title(r"H$_2$ Transition Probabilities in the Keck-MOSFIRE Filters",weight='bold',size='x-large')
plt.scatter(aul, wavelength, s=8, color='black')

ax1.fill_betweenx(mosfire_y, 0, 15, alpha=0.3, color='mediumblue', edgecolor='darkblue', label='MOSFIRE Y')
ax1.fill_betweenx(mosfire_j, 0, 15, alpha=0.3, color='limegreen', edgecolor='darkolivegreen', label='MOSFIRE J')
ax1.fill_betweenx(mosfire_h, 0, 15, alpha=0.3, color='deeppink', edgecolor='m', label='MOSFIRE H')
ax1.fill_betweenx(mosfire_k, 0, 15, alpha=0.3, color='red', edgecolor='firebrick', label='MOSFIRE K')

plt.legend(loc='upper right', fontsize='large')


# Aul*LyAcascade*0.01 vs. wavelength --> expected relative strength of lines vs wavelength
fig = plt.figure(figsize=(8,8),facecolor='white')
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'Relative Strength (F(Ly$\alpha$) [erg s$^{-1}$])',weight='bold',size='x-large')
ax1.set_ylabel(r'Wavelength ($\mu$m)',weight='bold',size='x-large')
ax1.set_ylim(0.95, 2.75)
ax1.set_xlim(-0.2, 65)
ax1.set_title(r"H$_2$ Relative Emission Strengths in the Keck-MOSFIRE Filters",weight='bold',size='x-large')
plt.scatter(aul*LyAcascade, wavelength, s=8, color='black')

ax1.fill_betweenx(mosfire_y, -0.2, 65, alpha=0.3, color='mediumblue', edgecolor='darkblue', label='MOSFIRE Y')
ax1.fill_betweenx(mosfire_j, -0.2, 65, alpha=0.3, color='limegreen', edgecolor='darkolivegreen', label='MOSFIRE J')
ax1.fill_betweenx(mosfire_h, -0.2, 65, alpha=0.3, color='deeppink', edgecolor='m', label='MOSFIRE H')
ax1.fill_betweenx(mosfire_k, -0.2, 65, alpha=0.3, color='red', edgecolor='firebrick', label='MOSFIRE K')

plt.legend(loc='upper right', fontsize='large')



plt.show()
