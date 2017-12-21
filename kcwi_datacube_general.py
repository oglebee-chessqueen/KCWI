### GEN. SCRIPT FOR KCWI DATA MANIPULATION ###
###
###
##
#
#
#
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy

# Define file
'''
Object						Date			File		Slicer	Grating/Wave	Time	N/S?
------						------		-----		------	------------	----	----
IC 1191						170514		143			Large		BH3,4950			0			Y
														144			Large		BH3						205		Y
														145			Large		BH3						0			Y
														146			Large		BH3						1015	Y
														147			Large		BH3						120		Y
														148			Large		BH3						120		Y
														149			Large		BH3						0			Y
														150			Large		BH3						187		Y
														151			Large		BH3						0			Y
														152			Large		BH3						0			Y
														153			Large		BH3						0			Y
														154			Large		BH3						261		Y
									170517		146			Large		BL,4500				120		Y
														147			Large		BL,4500				120		Y
														148			Large		BL,4500				120		Y
														149			Large		BL,4500				120		Y
														150			Large		BL,4500				120		Y
														151			Large		BL,4500				120		Y
														152			Large		BL,4500				120		Y
														153			Large		BL,4500				120		Y
														154			Large		BL,4500				120		Y
														155			Large		BL,4500				120		Y
														156			Large		BL,4500				120		Y
														157			Large		BL,4500				120		Y

BRN								170620		77			Small		BH2,4850			1			N
														78			Small		BH2,4850			300		N
														79			Small		BH2,4850			300		N
														80			Small		BH2,4850			300		N

HH 32							170614		176			Small		BM,4020				600		N
														177			Small		BM,4020				600		N
														178			Small		BM,4020				600		SKY
														179			Small		BM,4700				600		SKY
														180			Small		BM,4700				600		N
														181			Small		BM,4700				600		N
														186			Small		BM,5950				600		N
														187			Small		BM,5950				600		N
														188			Small		BM,5950				600		SKY

Ring Nebula				170412		210			Small		BH2, Hbeta		60		N
(M57, NGC6720)							211			Small		BH2, Hbeta		60		N
														212			Small		BL,4500				20		N
														213			Small		BL,4500				300		N
									170619		259			Large		BM,4550				1050	Y? (DON'T USE)
														260			Large		BM,4550				1050	Y?
														261			Large		BM,4550				950		Y?
									170620		64			Medium	BM,5200				157		N?
'''
path='C:\Users\Keri Hoadley\Documents\KCWI'
dir = '\\'
redux='redux'
date='170619'	#'170412'
index=260	#211
filename='kb'+date+'_00%03i_icube.fits' % index
print filename
file=path+dir+date+dir+redux+dir+filename
print file
print

# Constants:
c = 3.0e5		# Speed of light: km/s



# List of important line locations for NGC 6720
HI = [4101.73, 4340.47, 4861.35]  	# delta, gamma, beta
HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
HeII = [4540, 4686, 5411]
NI = [5198, 5200]
NII = [5755]
FeIII = [4658]
OI = [5577]
OII = [3726, 3729, 4267, 4661]
OIII = [4363.21, 4958.91, 5006.84, 5592.37]
NeIII = [3869, 3968]
NeIV = [4724, 4725]
SII = [4068, 4076]
ClII = [3679]
ClIII = [5517, 5537]
ClIV = [5323]
ArIII = [5192]
ArIV = [4711, 4739]
ArV = [4626]
KIV = [4511]
KV = [4123, 4163]
CaV = [3996, 5309]
UnKnown = []


hdulist = fits.open(file)		# Open the file

print hdulist.info()				# Print info on the fits file, gill get spatial
														# pix size and wavelength array

data = hdulist[0].data		# Define data cube (3D: wave, y, x)
#wv = hdulist[0].awav

print numpy.shape(data)		# Check size of data cube

wv1 = hdulist[0].header['WAVGOOD0']		# First "good" wavelength defined
wv2 = hdulist[0].header['WAVGOOD1']		# Last "good" wavelength defined
wmid = hdulist[0].header['WAVMID']					# Mid-way wavelength
wv0 = hdulist[0].header['CRVAL3']					# wavelength zeropoint
wv_interval = hdulist[0].header['CD3_3']		# Defines delta(wave)
nwv = numpy.size(data[:,0,0])
#print wv1, wv2, wmid, wv_interval
wlo = wmid - (nwv/2.)*wv_interval
whi = wmid + (nwv/2.)*wv_interval
wvlast = wv0 + nwv*wv_interval

#allwv = numpy.arange(wlo, whi, wv_interval)	# Defined using known delta(wave)
																							# and central wavelength
# allwv = numpy.arange(w1,w2, (wv2-wv1)/nwv)	# Defined using "good" wavelengths
																							# and total number of wavelength
																							# elements
all_wv = numpy.arange(wv0, wvlast, wv_interval)		# Defined using known delta(wave) and central wavelength
print wv0, wvlast, numpy.size(all_wv)

ra0 = hdulist[0].header['RA']#['CD1_1']		# RA (at center?)
dec0 = hdulist[0].header['DEC']#['CD2_2']		# DEC (at center?)

ra0_deg = hdulist[0].header['CRVAL1']	# zeropoint RA (at element 0) in degrees
dec0_deg = hdulist[0].header['CRVAL2']	# zeropoint DEC (at element 0)


delta_ra = hdulist[0].header['CD1_1']*3600.		# change in RA per pix (IN DEGREES)
nra = nwv = numpy.size(data[0,0,:])
ra_lim = (nra/2.)*delta_ra					# Extrema coverage to one end of the image
all_ra = numpy.arange(-1*ra_lim, ra_lim, delta_ra)	# array containing change in RA from center of image

delta_dec = hdulist[0].header['CD2_2']*3600.	# change in DEC per pix (IN DEGREES)
ndec = nwv = numpy.size(data[0,:,0])
dec_lim = (ndec/2.)*delta_dec					# Extrema coverage to one end of the image
all_dec = numpy.arange(-1*dec_lim, dec_lim, delta_dec)	# array containing change in DEC from center of image

print ra0, dec0
print delta_ra, delta_dec
print ra_lim,
print

hdulist.close()


i1 = 12
i2 = 14
# Find wavelength range for an arbitrary line of interest
l1 = HeI[6]	#HI[2]		# OIII[1]
# 1. Print a "narrow-band" region of the cube
# Define a range to extrapolate
dl1 = 2.0		# +/- 5 Ang on either side of l1
ind_l1 = numpy.where((all_wv >= l1-dl1) & (all_wv <= l1+dl1))[0]		# good!
# wave --> velocity array
vobs_l1 = ((all_wv[ind_l1] - l1) / l1)*c
print l1
print all_wv[ind_l1]
print vobs_l1
print vobs_l1[i1], vobs_l1[i2]
print

data_l1 = data[ind_l1,:,:]		# okay
max1= numpy.max(data_l1)
# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax1.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax1.set_title(r'He I: '+ra0+' '+dec0)
plt.imshow(numpy.sum(data_l1[i1:i2,7:135,:],axis=0), vmin=0, vmax=30, origin='lower', interpolation="spline36", cmap='Blues', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good
#plt.imshow(numpy.sum(data_l1[:,7:135,:],axis=0), vmin=1, vmax=max1/2, origin='lower', interpolation="spline36", cmap='Blues', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')
#plt.show()


# Find wavelength range for an arbitrary line of interest
l1 = HI[2]		# OIII[1]
# 1. Print a "narrow-band" region of the cube
# Define a range to extrapolate
dl1 = 2.0		# +/- 5 Ang on either side of l1
ind_l1 = numpy.where((all_wv >= l1-dl1) & (all_wv <= l1+dl1))[0]		# good!
# wave --> velocity array
vobs_l1 = ((all_wv[ind_l1] - l1) / l1)*c
print l1
print all_wv[ind_l1]
print vobs_l1
print vobs_l1[i1], vobs_l1[i2]
print

data_l1b = data[ind_l1,:,:]		# okay
max1= numpy.max(data_l1b)
# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax1.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax1.set_title(r'H$\beta$: '+ra0+' '+dec0)
plt.imshow(numpy.sum(data_l1b[i1:i2,7:135,:],axis=0), origin='lower', interpolation="spline36", cmap='Reds', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good
#plt.imshow(numpy.sum(data_l1[:,7:135,:],axis=0), vmin=1, vmax=max1/2, origin='lower', interpolation="spline36", cmap='Blues', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')


# Find wavelength range for an arbitrary line of interest
l2 = HeII[1]	#OIII[1]
# 1. Print a "narrow-band" region of the cube
# Define a range to extrapolate
dl2 = 2.0		# +/- 5 Ang on either side of l1
ind_l2 = numpy.where((all_wv >= l2-dl2) & (all_wv <= l2+dl2))[0]		# good!
vobs_l2 = ((all_wv[ind_l2] - l2) / l2)*c
print l2
print all_wv[ind_l2]
print vobs_l2
print vobs_l2[i1], vobs_l2[i2]

data_l2 = data[ind_l2,:,:]		# okay
max2 = numpy.max(data_l2)
# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
ax2.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax2.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax2.set_title('He II: '+ra0+' '+dec0)
plt.imshow(numpy.sum(data_l2[i1:i2,7:135,:],axis=0), origin='lower', interpolation="spline36", cmap='Blues', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good
#plt.imshow(numpy.sum(data_l2[:,7:135,:],axis=0), vmin=1, vmax=max2, origin='lower', interpolation="spline36", cmap='Blues', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')


# Ratio map: Show relative enhancement/dehancement of certain elements over
# the covered image
# 		[OIII]/Hbeta
ratio = numpy.sum(data_l1,axis=0) / numpy.sum(data_l2,axis=0)
minr = numpy.min(ratio)
maxr = numpy.max(ratio)
print
print minr, maxr

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title('Ionization Map: '+ra0+' '+dec0)
#~ if maxr > minr+0.1:
	#~ plt.imshow(ratio[7:135,:], vmin=-1*minr/2, vmax=-3*minr/4, origin='lower', interpolation="none", cmap='Blues', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good
#~ else:
plt.imshow(ratio[7:135,:], origin='lower', vmin=0, vmax=0.5, interpolation="none", cmap='viridis', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel(r'He I/He II',weight='bold',size='large')


# Ratio map: Show relative enhancement/dehancement of certain elements over
# the covered image
# 		[OIII]/Hbeta
ratio = numpy.sum(data_l2,axis=0) / numpy.sum(data_l1b,axis=0)
minr = numpy.min(ratio)
maxr = numpy.max(ratio)
print
print minr, maxr

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title('Ionization Map: '+ra0+' '+dec0)
#~ if maxr > minr+0.1:
	#~ plt.imshow(ratio[7:135,:], vmin=-1*minr/2, vmax=-3*minr/4, origin='lower', interpolation="none", cmap='Blues', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good
#~ else:
plt.imshow(ratio[7:135,:], origin='lower', vmin=0, vmax=0.5, interpolation="none", cmap='viridis', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel(r'He II/H$\beta$',weight='bold',size='large')



# Ratio map: Show relative enhancement/dehancement of certain elements over
# the covered image
# 		[OIII]/Hbeta
ratio = numpy.sum(data_l1,axis=0) / numpy.sum(data_l1b,axis=0)
minr = numpy.min(ratio)
maxr = numpy.max(ratio)
print
print minr, maxr

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title('Ionization Map: '+ra0+' '+dec0)
#~ if maxr > minr+0.1:
	#~ plt.imshow(ratio[7:135,:], vmin=-1*minr/2, vmax=-3*minr/4, origin='lower', interpolation="none", cmap='Blues', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good
#~ else:
plt.imshow(ratio[7:135,:], origin='lower', vmin=0, vmax=0.1, interpolation="none", cmap='viridis', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel(r'He I/H$\beta$',weight='bold',size='large')

plt.show()






