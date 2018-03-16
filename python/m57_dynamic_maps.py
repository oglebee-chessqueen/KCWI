#
#		Dynamic maps of M57 emission lines using predefined velocity bins
#

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy
plt.rcParams["font.family"] = "serif"



def open_fits_err(file,*args):
	'''
	Open the vcube fits file.
	Grab the variance cube + wavelength.
	*args - if set, then read in different RA, DEC header value
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


	hdulist.close()

	return data, all_wv

def open_fits(file,*args):
	'''
	Open the fits file.
	Grab the data.
	Get info on wavelength.
	Get info on observing range (RA,DEC).

	*args - if set, then read in different RA, DEC header value
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
	print len(args)

	if len(args) > 0:
		ra0 = hdulist[0].header['TARGRA']#['CD1_1']		# RA (at center?)
		dec0 = hdulist[0].header['TARGDEC']#['CD2_2']		# DEC (at center?)
	else:
		ra0 = hdulist[0].header['RA']		# RA (at center?)
		dec0 = hdulist[0].header['DEC']		# DEC (at center?)
	print type(ra0)
	print ra0[0:2], ra0[3:5], ra0[6:]
	hr = numpy.int(ra0[0:2])
	min = numpy.int(ra0[3:5])
	sec = numpy.float(ra0[6:])
	print hr, min, sec

	ra0_deg = 360.*( (hr + (min/60.) + (sec/3600.)) / 24.)

	print ra0_deg

	#~ print type(ra0)
	print dec0
	print dec0[1:3], dec0[4:6], dec0[7:]
	hr = numpy.int(dec0[1:3])
	min = numpy.int(dec0[4:6])
	sec = numpy.float(dec0[7:])
	print hr, min, sec

	dec0_deg = ( (hr + (min/60.) + (sec/3600.)) )

	print dec0_deg

	#ra0_deg = hdulist[0].header['CRVAL1']	# zeropoint RA (at element 0) in degrees
	#dec0_deg = hdulist[0].header['CRVAL2']	# zeropoint DEC (at element 0)

	delta_ra = hdulist[0].header['CD1_1']#*3600.		# change in RA per pix (IN DEGREES)
	nra = numpy.size(data[0,0,:])
	ra_lim = (nra/2.)*delta_ra					# Extrema coverage to one end of the image
	all_ra = numpy.arange(ra0_deg-ra_lim, ra0_deg+ra_lim, delta_ra)	# array containing change in RA from center of image
	#all_ra = numpy.arange(-1*ra_lim, ra_lim, delta_ra)	# array containing change in RA from center of image

	delta_dec = hdulist[0].header['CD2_2']#*3600.	# change in DEC per pix (IN DEGREES)
	ndec = numpy.size(data[0,:,0])
	dec_lim = (ndec/2.)*delta_dec
	print dec_lim	# Extrema coverage to one end of the image
	all_dec = numpy.arange(dec0_deg-dec_lim, dec0_deg+dec_lim, delta_dec)	# array containing change in DEC from center of image
	#all_dec = numpy.arange(-1*dec_lim, dec_lim, delta_dec)	# array containing change in DEC from center of image

	print ra0, dec0
	print ra0_deg, dec0_deg
	print ra_lim,dec_lim
	print ra0_deg-ra_lim, ra0_deg+ra_lim
	print dec0_deg-dec_lim, dec0_deg+dec_lim
	print

	hdulist.close()

	return data, all_wv, ra0_deg, dec0_deg, all_ra, all_dec


def narrowband(line, data, wv, dw, continuum_yn):
	'''
	Use to define a data array with limited wavelength coverage around line.
	continuum_yn is a string of length 1, either "y" or "n".
	If "y", then a continuum band of the same size as the narrowband image is
	subtracted (for example, to get rid of bright sources, if interested in
	diffuse media).
	If "n", then won't.
	Continuum is defined here as being 10lambda away from line
	'''
	ind_line = numpy.where((wv >= line-dw) & (wv <= line+dw))[0]		# good!
	if continuum_yn == "y":
		cline = line-10
		cont_line = numpy.where((wv >= cline-dw) & (wv <= cline+dw))[0]
		data_line = data[ind_line,:,:] - data[cont_line,:,:]  		# okay
	else:
		data_line = data[ind_line,:,:]

	return data_line, ind_line


def fit_data_to_vbin(data_in, vdata, vbin, ra, dec, *args):
	'''
	Use to shift narrow-band data bins into defined vbin intervals,
	to better compare flux/ratios in same dynamical view of the images

	*args = line ID to include in plot title
	'''
	# define a new data array to store slices that match with vbin
	data_vbin = numpy.zeros( [numpy.size(vbin), numpy.size(data_in[0,:,0]), numpy.size(data_in[0,0,:])] )
	dv = 25.	# size on either size of vbin[i] to allow data to be within velocity range
	# Loop through each vbin element, and find indices in vdata that are within
	# the defined range to be includes in the vbin data structure
	for i in range(0,numpy.size(vbin)):
		index = numpy.where((vdata >= vbin[i]-dv) & (vdata <= vbin[i]+dv))[0]
		print vbin[i], index
		if numpy.size(index) > 1:
			data_vbin[i] = numpy.sum(data_in[index,:,:],axis=0)/numpy.size(index)
		elif numpy.size(index) == 1:
			data_vbin[i] = data_in[index,:,:]

	# Plot each velocity bin to see how the dynamics look!
	ax = plt.figure()
	for img in range(0,numpy.size(vbin)):
		plt.subplot(2,5, img+1)
		# Only plot labels and title depending on some conditions:
		if img >= 5:
			plt.xlabel(r'$\Delta \alpha$ [$^{\prime \prime}$]',weight='bold',size='small')
		if (img%5) == 0:
			plt.ylabel(r'$\Delta \delta$ [$^{\prime \prime}$]',weight='bold',size='small')
		if len(args) > 0:
			plt.title(r'%s: $\Delta$v = %0.2d km/s' % (args[0],vbin[img]) )
		else:
			plt.title(r'$\Delta$v = %0.2d km/s' % vbin[img] )#r'NI 1598: '+ra0+' '+dec0)
		# Plot the dynamic view of the image!!
		plt.imshow(data_vbin[img,33:-28,12:-12]*1.0e17, origin='lower', vmin=0,
						interpolation="none", cmap='Greens',
						extent=[ra[12],ra[-12],dec[33],dec[-28]])
		cbar = plt.colorbar()		# Each image should have its own colorbar
		ax = plt.gca()
		ax.get_xaxis().get_major_formatter().set_useOffset(False)
		ax.get_yaxis().get_major_formatter().set_useOffset(False)
		# ...But only include label of colorbars at the end of each row
		if (img%5) == 4:
			cbar.ax.set_ylabel(r'Flux (10$^{-17}$ erg cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$)',
												 weight='bold')#,size='small')
	plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)

	#plt.show()

	return data_vbin



def extinction(vbin, data1, data2, I12, f12, ra, dec):
	'''
			Use to determine extinction (c_Hb)
	'''

	ax = plt.figure()
	for img in range(0,numpy.size(data1[:,0,0])):
		plt.subplot(2,5, img+1)
		# Only plot labels and title depending on some conditions:
		if img >= 5:
			plt.xlabel(r'$\Delta \alpha$ [$^{\prime \prime}$]',weight='bold',size='small')
		if (img%5) == 0:
			plt.ylabel(r'$\Delta \delta$ [$^{\prime \prime}$]',weight='bold',size='small')
		else:
			plt.title(r'$\Delta$v = %0.2d km/s' % vbin[img] )#r'NI 1598: '+ra0+' '+dec0)
		# Plot the dynamic view of the image!!
		ratio = data1[img,33:-28,12:-12]/data2[img,33:-28,12:-12]
		c_Hb = (-1./f12)*numpy.log10( (ratio/I12) )
		c_Hb[numpy.isnan(c_Hb)] = 0
		c_Hb[numpy.isneginf(c_Hb)] = 0
		c_Hb[numpy.isposinf(c_Hb)] = 0
		mn = numpy.nanmedian(c_Hb)
		print mn
		plt.imshow(c_Hb, origin='lower', vmin=mn-0.25, vmax=mn+0.25,
						interpolation="none", cmap='viridis',
						extent=[ra[12],ra[-12],dec[33],dec[-28]])
		cbar = plt.colorbar()		# Each image should have its own colorbar
		# ...But only include label of colorbars at the end of each row
		ax = plt.gca()
		ax.get_xaxis().get_major_formatter().set_useOffset(False)
		ax.get_yaxis().get_major_formatter().set_useOffset(False)
		if (img%5) == 4:
			cbar.ax.set_ylabel(r'c$_{H \beta}$',
												 weight='bold')#,size='small')
	plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)


	return






'''
Ring Nebula				170412		210			Small		BH2, Hbeta		60		N
(M57, NGC6720)							211			Small		BH2, Hbeta		60		N
														212			Small		BL,4500				20		N
														213			Small		BL,4500				300		N
									170415		188			Medium	BM,5900				10		N
									170619		259			Large		BM,4550				1050	Y? (DON'T USE)
														260			Large		BM,4550				1050	Y?
														261			Large		BM,4550				950		Y?
									170620		64			Medium	BM,5200				157		N?

Units of flux-cal data: erg/s/cm^2/Angstrom
'''
path='C:\Users\Keri Hoadley\Documents\KCWI'
dir = '\\'
redux='redux'
int = 'icubes'
var = 'vcubes'

date='170412'

index1=210		# not ready
index2=211		# not ready
index3=212		# for OII, OIII lines analysis
index4=213		# all other lines, use this higher S/N data set

# Use flux-calibrated data cubes
#~ filename1='kb'+date+'_00%03i_%s.fits' % (index1,int)
#~ file1=path+dir+date+dir+redux+dir+filename1

#~ filename2='kb'+date+'_00%03i_%s.fits' % (index2,int)
#~ file2=path+dir+date+dir+redux+dir+filename2

intfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,int)
file3 = path+dir+date+dir+redux+dir+intfile3
varfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,var)
vfile3 = path+dir+date+dir+redux+dir+varfile3

print file3

intfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,int)
file4 = path+dir+date+dir+redux+dir+intfile4
varfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,var)
vfile4 = path+dir+date+dir+redux+dir+varfile4

# Constants:
c = 3.0e5		#km/s
I_HgHb = 0.466
I_HdHb = 0.259
f_HgHb = 0.15		#0.13	#
f_HdHb = 0.22		#0.19	#
# Define velocity bin from -250 --> +250 km/s in intervals of 50 km/s (?)
vbin = numpy.arange(-250,250,50)

# List of important line locations for NGC 6720
HI = [4101.73, 4340.47, 4861.35 ]  	# delta, gamma, beta
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
#~ data, waves, ra, dec, all_ra, all_dec = open_fits(file1)
#~ data2, waves2, ra2, dec2, all_ra2, all_dec2 = open_fits(file2)
data3, waves3, ra3, dec3, all_ra3, all_dec3 = open_fits(file3)		# data in units erg cm-2 s-1
data4, waves4, ra4, dec4, all_ra4, all_dec4 = open_fits(file4)		# data in units erg cm-2 s-1
var3, varwv3 = open_fits_err(vfile3)
var4, varwv4 = open_fits_err(vfile4)

# Convert all_ra3 to a delta array away from central RA, and do same for DEC
# Multiply by 360 to convert to arcseconds
dra3 = (all_ra3 - ra3) * 3600.
ddec3 = (all_dec3 - dec3) * 3600
dra4 = (all_ra4 - ra4) * 3600
ddec4 = (all_dec4 - dec4) * 3600

# Define narrow bands around Hb, Hg, Hd
dw = 6

data3 = data4
waves3 = waves4
var3 = var4

line = HeI[5]		# Hgamma
data_wv1, interval_wv1 = narrowband(line, data3, waves3, dw, "n")
data_wv1 = data_wv1[1:,:,:]
interval_wv1 = interval_wv1[1:]
var_wv1 = var3[interval_wv1,:,:]
vobs_wv1 = ((waves3[interval_wv1] - line) / line)*c

data1_vbin = fit_data_to_vbin(data_wv1, vobs_wv1, vbin, all_ra3, all_dec3, r'HeI4922')


line = HI[2]		# Hbeta
data_wv2, interval_wv2 = narrowband(line, data3, waves3, dw, "n")
data_wv2 = data_wv2[1:-1,:,:]
interval_wv2 = interval_wv2[1:-1]
var_wv2 = var3[interval_wv2,:,:]
vobs_wv2 = ((waves3[interval_wv2] - line) / line)*c

data2_vbin = fit_data_to_vbin(data_wv2, vobs_wv2, vbin, all_ra3, all_dec3, r'H$\beta$')


line = OIII[0]		# Hdelta
data_wv3, interval_wv3 = narrowband(line, data3, waves3, dw, "n")
data_wv3 = data_wv3[1:-1,:,:]
interval_wv3 = interval_wv3[1:-1]
var_wv3 = numpy.sqrt(var3[interval_wv3,:,:])
vobs_wv3 = ((waves3[interval_wv3] - line) / line)*c

data3_vbin = fit_data_to_vbin(data_wv3, vobs_wv3, vbin, all_ra3, all_dec3, r'[OIII]4363')

# Calibrate to low flux limit (2.0e0-17)
data1_vbin[data1_vbin<1e-17] = 0
data2_vbin[data2_vbin<1e-17] = 0
data3_vbin[data3_vbin<1e-17] = 0

# What do the ratio maps of the H-lines look like in each velocity bin?
# Plot each velocity bin to see how the dynamics look!
# 1. Hg/Hb
ax = plt.figure()
for img in range(0,numpy.size(vbin)):
	plt.subplot(2,5, img+1)
	# Only plot labels and title depending on some conditions:
	if img >= 5:
		plt.xlabel(r'$\Delta \alpha$ [$^{\prime \prime}$]',weight='bold',size='small')
	if (img%5) == 0:
		plt.ylabel(r'$\Delta \delta$ [$^{\prime \prime}$]',weight='bold',size='small')
	plt.title(r'$\Delta$v = %0.2d km/s' % vbin[img] )#r'NI 1598: '+ra0+' '+dec0)
	# Plot the dynamic view of the image!!
	ratio1 = data1_vbin[img,33:-28,12:-12]/data2_vbin[img,33:-28,12:-12]
	ratio1[numpy.isnan(ratio1)] = 0
	mn = numpy.nanmedian(ratio1[ratio1>0])
	print mn
	#~ if img == 1 or img == 2 or img == 7 or img == 8:
		#~ plt.imshow(ratio1, origin='lower', #vmax=numpy.max(ratio1[img,30:-25,10:-10])/2, vmin=0,
						#~ interpolation="none", cmap='viridis', vmin=mn-0.15, vmax=mn+0.15,
						#~ extent=[all_ra3[12],all_ra3[-12],all_dec3[33],all_dec3[-28]])
	#~ elif img == 3:
		#~ plt.imshow(ratio1, origin='lower', #vmax=numpy.max(ratio1[img,30:-25,10:-10])/2, vmin=0,
						#~ interpolation="none", cmap='viridis', vmin=mn-0.05, vmax=mn+0.05,
						#~ extent=[all_ra3[12],all_ra3[-12],all_dec3[33],all_dec3[-28]])
	#~ else:
		#~ plt.imshow(ratio1, origin='lower', #vmax=numpy.max(ratio1[img,30:-25,10:-10])/2, vmin=0,
						#~ interpolation="none", cmap='viridis', vmin=mn-0.025, vmax=mn+0.025,
						#~ extent=[all_ra3[12],all_ra3[-12],all_dec3[33],all_dec3[-28]])
	if img == 4 or img == 5:
		plt.imshow(ratio1, origin='lower', #vmax=numpy.max(ratio1[img,30:-25,10:-10])/2, vmin=0,
						interpolation="none", cmap='viridis', vmin=mn-0.005, vmax=mn+0.005,
						extent=[all_ra3[12],all_ra3[-12],all_dec3[33],all_dec3[-28]])
	elif img == 6 or img == 3:
		plt.imshow(ratio1, origin='lower', #vmax=numpy.max(ratio1[img,30:-25,10:-10])/2, vmin=0,
						interpolation="none", cmap='viridis', vmin=mn-0.01, vmax=mn+0.01,
						extent=[all_ra3[12],all_ra3[-12],all_dec3[33],all_dec3[-28]])
	else:
		plt.imshow(ratio1, origin='lower', #vmax=numpy.max(ratio1[img,30:-25,10:-10])/2, vmin=0,
						interpolation="none", cmap='viridis', vmin=mn-0.1, vmax=mn+0.1,
						extent=[all_ra3[12],all_ra3[-12],all_dec3[33],all_dec3[-28]])
	cbar = plt.colorbar(format='%.02f')		# Each image should have its own colorbar
	# ...But only include label of colorbars at the end of each row
	if (img%5) == 4:
		cbar.ax.set_ylabel(r'HeI4922/H$\beta$',
											 weight='bold',size='large')
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)

#~ extinction(vbin, data3_vbin, data2_vbin, I_HdHb, f_HdHb, dra3, ddec3)

plt.show()
