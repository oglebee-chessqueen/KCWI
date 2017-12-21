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
import scipy.ndimage

def open_fits(file):
	'''
	Open the fits file.
	Grab the data.
	Get info on wavelength.
	Get info on observing range (RA,DEC).
	'''
	hdulist = fits.open(file)		# Open the file

	print hdulist.info()				# Print info on the fits file, gill get spatial
															# pix size and wavelength array

	data = hdulist[0].data		# Define data cube (3D: wave, y, x)

	#~ wv1 = hdulist[0].header['WAVGOOD0']		# First "good" wavelength defined
	#~ wv2 = hdulist[0].header['WAVGOOD1']		# Last "good" wavelength defined
	#~ wmid = hdulist[0].header['WAVMID']					# Mid-way wavelength
	wv0 = hdulist[0].header['CRVAL3']					# wavelength zeropoint
	wv_interval = hdulist[0].header['CD3_3']		# Defines delta(wave)
	nwv = numpy.size(data[:,0,0])

	#~ wlo = wmid - (nwv/2.)*wv_interval
	#~ whi = wmid + (nwv/2.)*wv_interval
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
	print ra_lim,dec_lim
	print

	hdulist.close()

	return data, all_wv, ra0, dec0, all_ra, all_dec, ra_lim, dec_lim



def narrowband(line, data, wv, dwv_blue, dwv_red, continuum_yn):
	'''
	Use to define a data array with limited wavelength coverage around line.
	continuum_yn is a string of length 1, either "y" or "n".
	If "y", then a continuum band of the same size as the narrowband image is
	subtracted (for example, to get rid of bright sources, if interested in
	diffuse media).
	If "n", then won't.
	Continuum is defined here as being 10lambda away from line
	'''
	ind_line = numpy.where((wv >= line-dwv_blue) & (wv <= line+dwv_red))[0]		# good!
	if continuum_yn == "y":
		cline = line-10
		cont_line = numpy.where((wv >= cline-dwv_blue) & (wv <= cline+dwv_red))[0]
		data_line = data[ind_line,:,:] - data[cont_line,:,:]  		# okay
	else:
		data_line = data[ind_line,:,:]

	return data_line, ind_line



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
									170415		188			Medium	BM,5900				10		N
									170619		259			Large		BM,4550				1050	Y? (DON'T USE)
														260			Large		BM,4550				1050	Y?
														261			Large		BM,4550				950		Y?
									170620		64			Medium	BM,5200				157		N?
'''
path='C:\Users\Keri Hoadley\Documents\KCWI'
dir = '\\'
redux='redux'
date1='170620'	#'170619'
index1=64	#260
filename1='kb'+date1+'_00%03i_icube.fits' % index1
file1=path+dir+date1+dir+redux+dir+filename1

date2='170415'	#
index2=188
filename2='kb'+date2+'_00%03i_icube.fits' % index2
file2=path+dir+date2+dir+redux+dir+filename2

date3='170619'	#'170619'
index3=260	#260
filename3='kb'+date3+'_00%03i_icube.fits' % index3
file3=path+dir+date3+dir+redux+dir+filename3

date4='170412'	#'170619'
index4=213	#260
filename4='kb'+date4+'_00%03i_icube.fits' % index4
file4=path+dir+date4+dir+redux+dir+filename4

date5 = date4
index5=211	#260
filename5='kb'+date5+'_00%03i_icube.fits' % index5
file5=path+dir+date5+dir+redux+dir+filename5
print

# Constants:
c = 3.0e5		# Speed of light: km/s

# List of important line locations for NGC 6720
HI = [4101.73, 4340.47, 4861.35]  	# delta, gamma, beta
HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
HeII = [4540, 4686, 5411]
NI = [5198.5, 5199]
NII = [5755]
FeIII = [4658]
OI = [5577,6300]
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


# Read in file, get important data from files
data, waves, ra, dec, all_ra, all_dec, dra, ddec = open_fits(file1)
data2, waves2, ra2, dec2, all_ra2, all_dec2, dra2, ddec2 = open_fits(file2)
data3, waves3, ra3, dec3, all_ra3, all_dec3, dra3, ddec3 = open_fits(file3)
data4, waves4, ra4, dec4, all_ra4, all_dec4, dra4, ddec4 = open_fits(file4)
data5, waves5, ra5, dec5, all_ra5, all_dec5, dra5, ddec5 = open_fits(file5)

#~ spectra1 = numpy.sum(data, axis=(1,2))
#~ plt.plot(waves,spectra1,drawstyle='steps-mid',color='black')
#~ plt.title(filename1)
#~ plt.xlabel(r'Wavelength ($\AA$)')
#~ plt.ylabel('Counts')

#~ spectra2 = numpy.sum(data2, axis=(1,2))
#~ plt.plot(waves2,spectra2,drawstyle='steps-mid',color='green')
#~ plt.title(filename2)
#~ plt.xlabel(r'Wavelength ($\AA$)')
#~ plt.ylabel('Counts')

#~ spectra3 = numpy.sum(data3, axis=(1,2))
#~ plt.plot(waves3,spectra3,drawstyle='steps-mid',color='blue')
#~ plt.title(filename3)
#~ plt.xlabel(r'Wavelength ($\AA$)')
#~ plt.ylabel('Counts')


#~ spectra5 = numpy.sum(data5, axis=(1,2))
#~ plt.plot(waves5,spectra5,drawstyle='steps-mid',color='purple')
#~ plt.title(filename5)
#~ plt.xlabel(r'Wavelength ($\AA$)')
#~ plt.ylabel('Counts')

#~ spectra4 = numpy.sum(data4, axis=(1,2))
#~ plt.plot(waves4,spectra4,drawstyle='steps-mid',color='orange')
#~ plt.title(filename4)
#~ plt.xlabel(r'Wavelength ($\AA$)')
#~ plt.ylabel('Counts')
#~ plt.show()

# Define small wavelength bandpass to image object
# Provide central wavelength, range around central wavelength, and
# whether to subtract out a continuum band (defined as "y" or "n")
dwv1_red = 0		# In case we only want to look at the blue/red region of a line
dwv1_blue = 3
line = NI[0]
data_wv1, interval_wv1 = narrowband(line, data, waves, dwv1_blue, dwv1_red, "y")
vobs_wv1 = ((waves[interval_wv1] - line) / line)*c

dwv2_red = 3		# In case we only want to look at the blue/red region of a line
dwv2_blue = 0
line = NI[1]
data_wv2, interval_wv2 = narrowband(line, data, waves, dwv2_blue, dwv2_red, "y")
vobs_wv2 = ((waves[interval_wv2] - line) / line)*c


dwv3_red = 0.5		# In case we only want to look at the blue/red region of a line
dwv3_blue = 3
line = NII[0]
data_wv3, interval_wv3 = narrowband(line, data2, waves2, dwv3_blue, dwv3_red, "y")
data_wv3[data_wv3<0] = 0.
vobs_wv3 = ((waves2[interval_wv3] - line) / line)*c
data_wv3[data_wv3<0] = 0.0

dwv4_red = 1.5		# In case we only want to look at the blue/red region of a line
dwv4_blue = 1.5
line = OI[1]
data_wv4, interval_wv4 = narrowband(line, data2, waves2, dwv4_blue, dwv4_red, "y")
data_wv4[data_wv4<0] = 0.
vobs_wv4 = ((waves2[interval_wv4] - line) / line)*c

print numpy.shape(data_wv4)

dwv5_red = 1.5		# In case we only want to look at the blue/red region of a line
dwv5_blue = 1.5
line = OI[0]
data_wv5, interval_wv5 = narrowband(line, data2, waves2, dwv5_blue, dwv5_red, "y")
vobs_wv5 = ((waves2[interval_wv5] - line) / line)*c


dwv6_red = 1.5		# In case we only want to look at the blue/red region of a line
dwv6_blue = 1.5
line = HeII[0]		#HeI[3]
data_wv6, interval_wv6 = narrowband(line, data3, waves3, dwv6_blue, dwv6_red, "y")
vobs_wv6 = ((waves3[interval_wv6] - line) / line)*c

dwv7_red = 1.5		# In case we only want to look at the blue/red region of a line
dwv7_blue = 1.5
line = HeI[3]
data_wv7, interval_wv7 = narrowband(line, data3, waves3, dwv7_blue, dwv7_red, "y")
vobs_wv7 = ((waves3[interval_wv7] - line) / line)*c

#~ print numpy.shape(data_wv6), numpy.shape(data_wv5)
#~ print waves3[interval_wv6]
#~ print numpy.sum(data_wv6[:,:,:],axis=0)
#~ print numpy.sum(data[interval_wv6,:,:],axis=0)


# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax1.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax1.set_title('M57 Central Star: [NI] 5198 - continuum')#r'NI 1598: '+ra0+' '+dec0)
levels=[40]
lw=[2.5]
CS = plt.contour(numpy.sum(data_wv1,axis=0), levels,
									linewidths=lw, colors='k', corner_mask=True,
									extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
plt.imshow(numpy.sum(data_wv1,axis=0), vmin=0, vmax=100, origin='lower',
						interpolation="none", cmap='Blues',
						extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')
#plt.show()


# Try to plot velocity v. declination
# Not quite what I want...
data_v1 = scipy.ndimage.zoom(data_wv1,[2,1,1])
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$\Delta$velocity [km s$^{-1}$]',weight='bold',size='large')
ax1.set_ylabel(r'$\Delta \delta$ [$^{\prime \prime}$]',weight='bold',size='large')
ax1.set_title('M57 Central Star: [NI] 5198 - continuum')#r'NI 1598: '+ra0+' '+dec0)
levels=[40]
lw=[2.5]
#~ CS = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									#~ linewidths=lw, colors='k', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(numpy.sum(data_wv1[:,25:60,10:-1],axis=2), vmin=0, vmax=200,
						origin='lower', interpolation="none", cmap='plasma',
						extent=[vobs_wv1[0],vobs_wv1[-1],all_dec[25],all_dec[60]])
plt.axes().set_aspect('auto')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')


# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
ax2.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax2.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax2.set_title('M57 Central Star: [NI] 5200 - continuum')#+ra0+' '+dec0)
levels2=[40]
lw2=[2.5]
CS = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									linewidths=lw2, colors='k', corner_mask=True,
									extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(numpy.sum(data_wv2[:,:,:],axis=0), vmin=0, vmax=100, origin='lower',
						interpolation="none", cmap='Reds',
						extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')


ratio = numpy.sum(data_wv1,axis=0) / numpy.sum(data_wv2,axis=0)
ratio[ratio<0] = 0.
ratio[ratio>3.5] = 0.

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title('M57 Central Star: Medium Slicer')
CS1 = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									linewidths=lw, colors='b', corner_mask=True,
									extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
CS2 = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									linewidths=lw2, colors='r', corner_mask=True,
									extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(ratio[:,:], origin='lower', vmin=0, vmax=2,
						interpolation="none", cmap='Greys',
						extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel(r'[NI] 5198/5200',weight='bold',size='large')
#~ #cbar.ax.set_ylabel(r'Photon Counts',weight='bold',size='large')


# Try to contourf the ratio
levels=[0.0,0.5,1.0,1.5,2.0,2.5]
fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title('M57 Central Star: Medium Slicer')
cs = plt.contourf(ratio, 25, cmap='gist_heat', corner_mask=True,
						extent=[all_ra[-1],all_ra[0],all_dec[0],all_dec[-1]])#, origin='lower') # inferno - good
#~ CS1 = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									#~ linewidths=lw, colors='b', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
#~ CS2 = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									#~ linewidths=lw2, colors='y', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'[NI] 5198/5200',weight='bold',size='large')


plt.show()


# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax1.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax1.set_title('M57 Central Star: [NII] 5755 - continuum')#r'NI 1598: '+ra0+' '+dec0)

plt.imshow(numpy.sum(data_wv3[:,:,:],axis=0), vmin=0, vmax=50, origin='lower',
						interpolation="none", cmap='Greens',
						extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')



#~ # Plot narrow-band image of data_l1
#~ #norm=LogNorm(vmin=500, vmax=15000),
#~ fig = plt.figure()
#~ ax1 = fig.add_subplot(1,1,1)
#~ ax1.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
#~ ax1.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
#~ ax1.set_title('M57 Central Star: [OI] 6300 - continuum')#r'NI 1598: '+ra0+' '+dec0)

#~ plt.imshow(numpy.sum(data_wv4[:,:,:],axis=0), vmin=0, vmax=100, origin='lower',
						#~ interpolation="none", cmap='RdPu',
						#~ extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
#~ cbar = plt.colorbar()
#~ cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')


# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
#~ #fig = plt.figure()
#~ tnum = numpy.size(data_wv4[:,0,0])+2
#~ col = [1, 1, 1, 2, 3, 3, 3]
#~ row = [1, 2, 3, 2, 1, 2, 3]
#~ for img in range(0,numpy.size(data_wv4[:,0,0])):
	#~ plt.subplot(3,3, img+1)
	#~ plt.xlabel(r'$\Delta \alpha$ [$^{\prime \prime}$]',weight='bold',size='small')
	#~ plt.ylabel(r'$\Delta \delta$ [$^{\prime \prime}$]',weight='bold',size='small')
	#~ plt.title('[OI] 6300: v = %0.2d' % vobs_wv4[img] )#r'NI 1598: '+ra0+' '+dec0)

	#~ plt.imshow(data_wv4[img,:,:], vmin=0, vmax=50, origin='lower',
						#~ interpolation="none", cmap='RdPu',
						#~ extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
	#~ #count += 1

	#~ cbar = plt.colorbar()
	#~ cbar.ax.set_ylabel('Photon Counts',weight='bold',size='small')


#~ tnum = numpy.size(data_wv5[:,0,0])+2
#~ col = [1, 1, 1, 2, 3, 3, 3]
#~ row = [1, 2, 3, 2, 1, 2, 3]
#~ for img in range(0,numpy.size(data_wv5[:,0,0])):
	#~ plt.subplot(3,3, img+1)
	#~ plt.xlabel(r'$\Delta \alpha$ [$^{\prime \prime}$]',weight='bold',size='small')
	#~ plt.ylabel(r'$\Delta \delta$ [$^{\prime \prime}$]',weight='bold',size='small')
	#~ plt.title('[OI] 5577: v = %0.2d' % vobs_wv4[img] )#r'NI 1598: '+ra0+' '+dec0)

	#~ plt.imshow(data_wv5[img,:,:], vmin=0, origin='lower',
						#~ interpolation="none", cmap='RdPu',
						#~ extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
	#~ #count += 1

	#~ cbar = plt.colorbar()
	#~ cbar.ax.set_ylabel('Photon Counts',weight='bold',size='small')
data_wv4[data_wv4>200] = 0.
#~ data_wv4 = scipy.ndimage.zoom(data_wv4,2)
#~ print numpy.size(data_wv3[:,0,0])
#~ data_wv3 = scipy.ndimage.zoom(data_wv3,[1,2,2])
#~ print numpy.size(data_wv3[:,0,0])
#~ levels=[160]
#~ lw=[2.5]
levels=[70]
lw=[2.5]
tnum = numpy.size(data_wv3[:,0,0])+2
col = [1, 1, 1, 2, 3, 3, 3]
row = [1, 2, 3, 2, 1, 2, 3]
for img in range(0,numpy.size(data_wv3[:,0,0])):
	plt.subplot(3,3, img+1)
	if img > 5:
		plt.xlabel(r'$\Delta \alpha$ [$^{\prime \prime}$]',weight='bold',size='small')
	if (img%3) == 0:
		plt.ylabel(r'$\Delta \delta$ [$^{\prime \prime}$]',weight='bold',size='small')

	plt.title(r'[NII] 5755: $\Delta$v = %0.2d km/s' % vobs_wv3[img] )#r'NI 1598: '+ra0+' '+dec0)
	CS = plt.contour(numpy.sum(data_wv4[:,:,:],axis=0), levels,
									linewidths=lw, colors='k', corner_mask=True,
									extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
	plt.imshow(data_wv3[img,:,:], vmin=0, vmax=10, origin='lower',
						interpolation="none", cmap='Blues',
						extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
	#count += 1

	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Photon Counts',weight='bold',size='small')
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)


#data_wv3 = scipy.ndimage.zoom(data_wv3,[1,5,5])
# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax1.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax1.set_title('M57 Central Star: [OI] 6300 - continuum')#r'NI 1598: '+ra0+' '+dec0)
CS = plt.contour(numpy.sum(data_wv4[:,:,:],axis=0), levels,
									linewidths=lw, colors='k', corner_mask=True,
									extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
plt.imshow(numpy.sum(data_wv3[1:3,:,:],axis=0), vmin=0, vmax=15, origin='lower',
						interpolation="none", cmap='Blues',
						extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')

# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax1.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax1.set_title('M57 Central Star: [OI] 6300 - continuum')#r'NI 1598: '+ra0+' '+dec0)
CS = plt.contour(numpy.sum(data_wv4[:,:,:],axis=0), levels,
									linewidths=lw, colors='k', corner_mask=True,
									extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
plt.imshow(numpy.sum(data_wv3[4:6,:,:],axis=0), vmin=0, vmax=20, origin='lower',
						interpolation="none", cmap='Blues',
						extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')

# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax1.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax1.set_title('M57 Central Star: [OI] 6300 - continuum')#r'NI 1598: '+ra0+' '+dec0)
CS = plt.contour(numpy.sum(data_wv4[:,:,:],axis=0), levels,
									linewidths=lw, colors='k', corner_mask=True,
									extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
plt.imshow(numpy.sum(data_wv4[:,:,:],axis=0), vmin=0, vmax=150, origin='lower',
						interpolation="none", cmap='RdPu',
						extent=[all_ra2[0],all_ra2[-1],all_dec2[5],all_dec2[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')




# Plot narrow-band image of data_l1
#norm=LogNorm(vmin=500, vmax=15000),
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax1.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
ax1.set_title(r'M57 Central Star: HeII 4540 - continuum')#r'HeI 4471: '+ra0+' '+dec0)

plt.imshow(numpy.sum(data_wv7[:,:,:],axis=0), vmin=0, vmax=500, origin='lower',
						interpolation="none", cmap='Purples',
						extent=[all_ra3[0],all_ra3[-1],all_dec3[5],all_dec3[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')

# Ratio map: Show relative enhancement/dehancement of certain elements over
# the covered image
# 		[OIII]/Hbeta
ratio = numpy.sum(data_wv6,axis=0) / numpy.sum(data_wv7,axis=0)
#~ ratio = data_wv7/data_wv6


fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.set_xlabel(r'$\Delta$RA (arcsec)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta$DEC (arcsec)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title('M57 Central Star: Large Slicer')
#~ #if maxr > minr+0.1:
	#~ #plt.imshow(ratio[7:135,:], vmin=-1*minr/2, vmax=-3*minr/4, origin='lower', interpolation="none", cmap='Blues', extent=[all_ra[0],all_ra[23],all_dec[5],all_dec[135]])#'hot')#, origin='lower') # inferno - good
#~ #else:
plt.imshow(ratio, origin='lower', vmin=0, vmax=1,
						interpolation="gaussian", cmap='viridis',
						extent=[all_ra3[0],all_ra3[-1],all_dec3[0],all_dec3[-1]])#, origin='lower') # inferno - good
#~ #plt.imshow(numpy.sum(data[interval_wv1,:,:]+data[interval_wv2,:,:],axis=0),
						#~ #vmin=0, vmax=300, origin='lower', interpolation="none",
						#~ #cmap='viridis', extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel(r'HeII 4540/HeI 4471',weight='bold',size='large')
#~ #cbar.ax.set_ylabel(r'Photon Counts',weight='bold',size='large')





plt.show()








