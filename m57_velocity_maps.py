from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy
import numpy.polynomial.polynomial as poly
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

	wv0 = hdulist[0].header['CRVAL3']					# wavelength zeropoint
	wv_interval = hdulist[0].header['CD3_3']		# Defines delta(wave)
	nwv = numpy.size(data[:,0,0])
	wvlast = wv0 + nwv*wv_interval
	all_wv = numpy.arange(wv0, wvlast, wv_interval)		# Defined using known delta(wave) and central wavelength

	if len(args) > 0:
		ra0 = hdulist[0].header['TARGRA']#['CD1_1']		# RA (at center?)
		dec0 = hdulist[0].header['TARGDEC']#['CD2_2']		# DEC (at center?)
	else:
		ra0 = hdulist[0].header['RA']		# RA (at center?)
		dec0 = hdulist[0].header['DEC']		# DEC (at center?)

	# Define RA in degrees
	hr = numpy.int(ra0[0:2])
	min = numpy.int(ra0[3:5])
	sec = numpy.float(ra0[6:])
	ra0_deg = 360.*( (hr + (min/60.) + (sec/3600.)) / 24.)
	# Define DEC in degrees
	hr = numpy.int(dec0[1:3])
	min = numpy.int(dec0[4:6])
	sec = numpy.float(dec0[7:])
	dec0_deg = ( (hr + (min/60.) + (sec/3600.)) )

	delta_ra = hdulist[0].header['CD1_1']#*3600.		# change in RA per pix (IN DEGREES)
	nra = numpy.size(data[0,0,:])
	ra_lim = (nra/2.)*delta_ra					# Extrema coverage to one end of the image
	all_ra = numpy.arange(ra0_deg-ra_lim, ra0_deg+ra_lim, delta_ra)	# array containing change in RA from center of image

	delta_dec = hdulist[0].header['CD2_2']#*3600.	# change in DEC per pix (IN DEGREES)
	ndec = numpy.size(data[0,:,0])
	dec_lim = (ndec/2.)*delta_dec
	all_dec = numpy.arange(dec0_deg-dec_lim, dec0_deg+dec_lim, delta_dec)	# array containing change in DEC from center of image

	hdulist.close()

	return data, all_wv, ra0_deg, dec0_deg, all_ra, all_dec



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


def mine_ring():
	'''
	Use to look at 04/12 data sets (co-added for max S/N)
	'''
	############### Main Ring (inner edge) ##################
	# Line list
	HI = [4101.73, 4340.47, 4861.35]  	# delta, gamma, beta
	HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
	HeII = [4540, 4686, 5411]
	NI = [5198.5, 5199]
	OI = [5577,6300]
	OII = [3726, 3729, 4267, 4661]
	OIII = [4363.21, 4958.91, 5006.84, 5592.37]

	date='170412'

	int = 'icubes'
	var = 'vcubes'

	index1=210		# spectral range: 4650 - 5050
	index2=211		#
	index3=212		# spectral range: 3600 - 5500
	index4=213		#

	intfile1 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index1,index2,int)
	file1 = path+dir+date+dir+redux+dir+intfile1
	varfile1 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index1,index2,var)
	vfile1 = path+dir+date+dir+redux+dir+varfile1

	intfile2 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,int)
	file2 = path+dir+date+dir+redux+dir+intfile2
	varfile2 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,var)
	vfile2 = path+dir+date+dir+redux+dir+varfile2


	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1)		# data in units erg cm-2 s-1
	data2, waves2, ra, dec, all_ra, all_dec = open_fits(file2)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)
	var2, varwv2 = open_fits_err(vfile2)

	# Do things here



	return



def central_star():
	'''
	Use to look at 06/2
	0 data set (with [NI] ONLY)
	'''
	############### Main Ring (inner edge) ##################
	# Line list
	NI = [5198.5, 5199]

	############### Central Region 1: 06/20 ##################
	date='170620'

	int = 'icuber'
	var = 'vcuber'

	index1 = 64		# not ready

	intfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,int)
	file1 = path+dir+date+dir+redux+dir+intfile1
	varfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,var)
	vfile1 = path+dir+date+dir+redux+dir+varfile1


	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)

	# Do things here



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
path='C:\\Users\\Keri Hoadley\\Documents\\KCWI'
dir = '\\'
redux='redux'

# Line list
HI = [4101.73, 4340.47, 4861.35]  	# delta, gamma, beta
HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
HeII = [4540, 4686, 5411]
NI = [5198.5, 5199]
OI = [5577,6300]
OII = [3726, 3729, 4267, 4661]
OIII = [4363.21, 4958.91, 5006.84, 5592.37]













############### Central Region 1: 06/20 ##################
date='170620'

int = 'icuber'
var = 'vcuber'

index3 = 64		# not ready

intfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,int)
file3 = path+dir+date+dir+redux+dir+intfile3
varfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,var)
vfile3 = path+dir+date+dir+redux+dir+varfile3



