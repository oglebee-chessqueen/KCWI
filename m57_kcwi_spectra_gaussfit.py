#
#
# Gaussian line shape fitting routine
#
#
# Read in: KCWI fits data cubes (ideally flux-calibrated).
#					 Emission lines DAT file w/ rest wavelengths, IDs, Aul, and statistical weight.
#
# Goal: Calculate the resolving power of CHESS as a function of wavelength
#
#
#

import numpy
import numpy.polynomial.polynomial as poly
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
plt.rcParams["font.family"] = "serif"



def open_fits(file,*args):
	'''
	Open the vcube fits file.
	Grab the variance cube + wavelength.
	*args - if set, then read in different RA, DEC header value (NOT APPLICABLE)
	'''
	hdulist = fits.open(file)		# Open the file
	data = hdulist[0].data		# Define data cube (3D: wave, y, x)
	wv0 = hdulist[0].header['CRVAL3']					# wavelength zeropoint
	wv_interval = hdulist[0].header['CD3_3']		# Defines delta(wave)
	hdulist.close()

	nwv = numpy.size(data[:,0,0])
	wvlast = wv0 + nwv*wv_interval
	all_wv = numpy.arange(wv0, wvlast, wv_interval)		# Defined using known delta(wave) and central wavelength

	# Limit x, y image area to only include "good" image area
	# (i.e., no blank space or edge pixels)
	# Same for wavelength space
	data = data[400:-300,33:-28,12:-12]
	all_wv = all_wv[400:-300]

	if len(args) > 0:
		ra0 = hdulist[0].header['TARGRA']#['CD1_1']		# RA (at center?)
		dec0 = hdulist[0].header['TARGDEC']#['CD2_2']		# DEC (at center?)
	else:
		ra0 = hdulist[0].header['RA']		# RA (at center?)
		dec0 = hdulist[0].header['DEC']		# DEC (at center?)

	hr = numpy.int(ra0[0:2])
	min = numpy.int(ra0[3:5])
	sec = numpy.float(ra0[6:])

	ra0_deg = 360.*( (hr + (min/60.) + (sec/3600.)) / 24.)

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

	#~ print all_dec[33], all_dec[-28], all_ra[12],all_ra[-12]
	#~ print

	return data, all_wv, all_ra, all_dec



def open_lineID(file,*args):
	'''
	Open the DAT file with line IDs of emission per observation.
	Grab the line ID + rest wavelength.
	'''
	hlines = open(file).readlines()

	# Read through each line, and append colymn values to correct arrays
	linewave = []
	lineID = []

	for line in hlines:
		if line.startswith('#'):
			continue
		column = line.split()
		lineID.append( str(column[0]) )
		linewave.append( float(column[2]) )

	return lineID, linewave



def create_spectum(data,*args):
	'''
	Collapse image area to create 1D spectrum.
	No args = collapse entire image area
	args: [RA:RA,DEC:DEC] - define image area to collapse along wavelength
	'''
	if len(args) == 0:
		return numpy.sum(data,axis=(1,2))
	if len(args) == 2:
		return numpy.sum(data4[:,args[0]:args[1],:],axis=(1,2))
	if len(args) == 4:
		return numpy.sum(data4[:,args[0]:args[1],args[2]:args[3]],axis=(1,2))
	else:
		return "INCORRECT NUMBER OF ARGUMENTS - NO SPECTRUM CREATED"



def create_spectum_err(var,*args):
	'''
	SAME AS create_spectrum, but for the error bars
	'''
	if len(args) == 0:
		return numpy.sum(numpy.sqrt(abs(var)), axis=(1,2))
	if len(args) == 2:
		return numpy.sum(numpy.sqrt(abs(var[:,args[0]:args[1],:])),axis=(1,2))
	if len(args) == 4:
		return numpy.sum(numpy.sqrt(abs(var[:,args[0]:args[1],args[2]:args[3]])),axis=(1,2))
	else:
		return "INCORRECT NUMBER OF ARGUMENTS - NO SPECTRUM CREATED"



def fit_continuum_to_spectrum(spectrum, waves, bin_size):
	'''
	Bin the spectra by ~ [bin_size] pixels:
		Use central wavelenth every [bin_size] steps.
		Take median (or avg) of [bin_size] pixels considered.
		These valaues will make up points a long the continuum to fit with
		a polynomial-fitting routine.
			- Then, check polynomial fit with least-squares statistic.
			- When fitting means a certain criteria (fit doesn't improve any more
				with increasing polynomial factor), then note where the fit doesn't
				improve and say that fit is good to nth polynomial
				(THIS SHOULD BE A SEPARATE FUNCTION)

	Input:
	------
	spectrum - 1D array: collapsed data cube as a function of wavelength only
	waves - 1D array: corresponding wavelength array to spectrum
	bin_size - Int (1 number): size of spectral region to find continuum value for

	Output:
	-------
	wave_bin - 1D array: wavelength array where continuum values were found
	flux_bin - 1D array: flux array with continuum values per wave_bin
	'''
	wave_bin = numpy.zeros( numpy.size(waves)/bin_size )
	flux_bin = numpy.zeros( numpy.size(waves)/bin_size )
	count = 0
	for i in range(0, numpy.size(waves)/bin_size):
		if count+bin_size > numpy.size(waves):
			continue
		else:
			tmpwave = waves[count:count+bin_size]
			wave_bin[i] = tmpwave[bin_size/2]
			flux_bin[i] = numpy.median(spectrum[count:count+bin_size])
		count += bin_size
		# check that flux_bin > 0
		if flux_bin[i] <= 0 or numpy.isnan(flux_bin[i]):
			flux_bin[i] = flux_bin[i-1]

	return wave_bin, flux_bin



def continuum_fit_lsq(wave_bin, flux_bin, waves, flux, err):
	'''
	Use continuum binned arrays to find a polynomial fit through the flux continuum.

	Loop through each polynomial fit, and find least-squares statistic between
		total flux array and fitted polynomial spectrum (continuum).
	When statistic does not improve appreciably after the next iteration,
		- note where (index) the statistic did not improve
		- take the previous iteration's polynomial fit as best fit through the continuum.
	Return polynomial fitted array (new "flux" array)

	Inputs:
	-------
	wave_bin - 1D array: wavelength array corresponding to continuum array
	flux_bin - 1D array: continuum array (to fit continuum using polyfit)
	waves - 1D array: entire wavelength array to use to find entire continuum
										and perform least-squares statistics with
	flux - 1D array: entire 1D spectrum to use to find entire continuum
									and perform least-squares statistics with
	err - 1D array: error bars of 1D spectrum

	Outputs:
	--------
	polyfit_vals
	best_continuum_array
	'''
	# Use the numpy polyfit functions
	max_index = 10
	ls = numpy.zeros( max_index )
	#~ ls_min = 0
	#~ continuum_array = numpy.array()
	for i in range(1,max_index):		# increase in i = increase to polynomial factor
		poly_coeffs = poly.polyfit(wave_bin[1:-1], flux_bin[1:-1], i)		# find best coefficients through continuum values
		#cfit = poly.Polynomial(poly_coeffs)			# leave it in functional form (do cfit(waves))
		cfit = poly.polyval(waves,poly_coeffs)	# fit coefficients through wavelength of spectrum
		#~ plt.semilogy(wave_bin[1:-1],flux_bin[1:-1],drawstyle='steps-mid',color='black',lw=3)
		#~ plt.semilogy(waves,cfit,color='red',lw=3)
		#~ plt.semilogy(waves,flux,drawstyle='steps-mid',color='blue',lw=3)
		#~ plt.show()
		# Perform least-squares between cfit and flux:
		# 							[ (flux[i] - cfit[i])**2 ]
		#				LS = SUM|----------------------- |
		#								[       err[i]**2				 ]
		for k in range(0,numpy.size(waves)):
			if flux[k] <= 0 or numpy.isnan(flux[k]):
				continue
			else:
				ls[i] += (flux[k] - cfit[k])**2 / (err[k]**2)
		ls[i] = ls[i]/(numpy.size(waves)) # - numpy.size(poly_coeffs))

		#~ print 'Polynomial: %i - Least-squares: %0.3f' %(i,ls[i])

		# Store lowest least-squares continuum fit to continuum_array
		if i == 1:
			continuum_array = cfit
			ls_min = ls[i]
		else:
			if ls[i] < ls_min:
				continuum_array = cfit
				ls_min = ls[i]

	return continuum_array, ls_min


def subtract_continuum(data,continuum):
	'''
	Subtract the continuum flux from the full spectrum at whatever
	image bin continuum was determined for.
	Will be left with a "zeroed-flux" continuum, which is necessary
	for performing multi-Gaussian profile fitting of all spectral lines.
	'''
	return data - continuum



def gaus(xdata,amp,cen,stdv):
	"""
	This function creats a gaussian curve by taking in values
	for its amplitude (amp), center (cen), and width/standard
	deviation (wid). You also must supply the x value data
	that it will be iterating over.
	input
	-----
	xdata = the independent value array
	amp = the amplitude of the peak emission line (can be an array of amplitudes)
	cen = the central wavelength of the peak (can be an array of centers)
	stdv = the standard deviation of the emission line (can be an array of stdvs)
	"""

	return abs(amp) * numpy.exp(-(xdata-cen)**2 /(2.*stdv**2))



def gaussum(xdata,*params):
	"""
	This function will sum the different gaussian fits together. It takes in the x-data array, as well
	as arrays for the three parameters for a Gaussian curve. Make sure that if you are going to use gassum
	on its own to include in the * when you're inputting the parameter array.
	input
	-----
	xdata = the independent value array
	params = array of all gaussian parameters alternating like [amp1,cen1,stdv1,amp2,cen2,...]
	"""
	amp = numpy.zeros(0)
	cen = numpy.zeros(0)
	stdv = numpy.zeros(0)

	for i in range(0, len(params), 3): #This section is just unpacking the parameter array into amps, cens, and stdvs
		x = params[i]
		amp = numpy.append(amp,x)
		y = params[i+1]
		cen = numpy.append(cen,y)
		z = params[i+2]
		stdv = numpy.append(stdv,z)
	global storage #You may not need storage to be global so think about taking this part out. storage stores the data
	storage = [[0 for x in range(1)] for x in range(len(params)/3)] #from each iteration of the gaussian equation into
	for i in range(len(params)/3):#individual rows. So row one will be the gaussian solutions to the first peak and so on
		storage[i] = gaus(xdata,amp[i],cen[i],stdv[i])
	storage = numpy.asarray(storage)
	return sum(storage)



# Define reference path and fits file, dat file here
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
#~ path='/home/keri/KCWI'
#~ dir = '/'
redux='redux'
int = 'icubes'
var = 'vcubes'

date='170412'

index1=210		# not ready
index2=211		# not ready
index3=212		# for OII, OIII lines analysis
index4=213		# all other lines, use this higher S/N data set

intfile3 = 'kb'+date+'_00%03i_%s_extcorr.fits' % (index3,int)
file3 = path+dir+date+dir+redux+dir+intfile3
varfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,var)
vfile3 = path+dir+date+dir+redux+dir+varfile3


intfile4 = 'kb'+date+'_00%03i_%s_extcorr.fits' % (index4,int)
file4 = path+dir+date+dir+redux+dir+intfile4
varfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,var)
vfile4 = path+dir+date+dir+redux+dir+varfile4


# Read in file, get data + waves from files
data3, waves3, ra, dec = open_fits(file3)		# data in units erg cm-2 s-1
data4, waves4, ra, dec = open_fits(file4)		# data in units erg cm-2 s-1
var3, varwv3, vra, vdec = open_fits(vfile3)
var4, varwv4, vra, vdec = open_fits(vfile4)

# Open emission line template file:
# Important columns in file: Line ID (col 0), Lab Wave (col 2)
linesfile = date+"_212-213_lines.dat"

# Open DAT file with emission line info
lineID, linewave = open_lineID(path+dir+linesfile)
i_Hb = lineID.index('Hbeta')
print lineID
print lineID.index('Hbeta')
print


# Define image "bins" to loop through and find spectra,
# continuum level(s), and emission line characteristics.
# Image size is: (wave, 132L, 22L)
# TEST BINS:
ra_bin = numpy.size(data4[0,0,:])/2			# 2
dec_bin = numpy.size(data4[0,:,0])/2		# 33
#REAL BINS:
#~ ra_bin = 2
#~ dec_bin = 2

bin_size = 20		# Bin size of continuum determination in each spectrum
#~ print 'RA bin: ', ra_bin
#~ print 'DEC bin: ', dec_bin
#~ print

# Make a continuum flux 2D array with size of total # of elements in ra, dec_bin
continuum_img = numpy.zeros( [numpy.size(data4[0,:,0])/dec_bin,
															numpy.size(data4[0,0,:])/ra_bin] )
chi2 = numpy.zeros( [numpy.size(data4[0,:,0])/dec_bin,
										 numpy.size(data4[0,0,:])/ra_bin] )
#~ print numpy.shape(continuum_img)
#~ print numpy.size(data4[0,:,0])-1
#~ print numpy.size(data4[0,0,:])-1
#~ print

ind_i = 0

# Loop through ra, dec bin sizes and do spectral operations
for i in range(0,numpy.size(data4[0,:,0])-1,dec_bin):
	ind_j = 0
	for j in range(0,numpy.size(data4[0,0,:])-1,ra_bin):
		#~ print 'Analyzing for pixels: RA [%i, %i], DEC [%i, %i]' %
		#~ 				(j, j+ra_bin, i, i+dec_bin)

		data_bin = data4[:,i:i+dec_bin,j:j+ra_bin]
		var_bin = var4[:,i:i+dec_bin,j:j+ra_bin]

		# Create spectra (Collapse along x,y axes)
		spectra = create_spectum(data_bin)
		err = create_spectum_err(var_bin)

		# Find continuum fit:
		# 1. Define a bin size to find medium flux through and extrac new wave, flux
		#			arrays for continuum
		wave_bin, flux_bin = fit_continuum_to_spectrum(spectra, waves4, bin_size)
		#
		# 2. Use new continuum (flux) array to fit polynomial through entire
		# 		wavelength range.
		#			Return continuum array and least-squares min chi2 value of fit to data.
		continuum_flux, chi2[ind_i,ind_j] = continuum_fit_lsq(wave_bin, flux_bin, waves4, spectra, err)
		continuum_img[ind_i,ind_j] = numpy.sum(continuum_flux)
		print 'At [%i - %i, %i - %i] -- Total continuum flux: %.3e   Chi2 = %.3f' % (i, i+dec_bin, j, j+ra_bin, numpy.sum(continuum_flux),chi2[ind_i,ind_j])

		# Subtract continuum flux from spectrum
		spectra_0 = subtract_continuum(spectra, continuum_flux)

		#~ # Plot spectra (log)
		#~ fig1 = plt.figure()
		#~ ax1 = fig1.add_subplot(2,1,1)
		#~ plt.semilogy(waves4, numpy.sum(data_bin,axis=(1,2)), drawstyle='steps-mid', lw=3,
						 #~ color='black')
		#~ plt.semilogy(waves4, continuum_flux, drawstyle='steps-mid', lw=3, color='blue')
		#~ plt.semilogy(wave_bin, flux_bin, 'ro', ms = 5)
		#~ plt.ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$ $\AA ^{-1}$)')
		#~ ax2 = fig1.add_subplot(2,1,2)
		#~ plt.plot(waves4, numpy.sum(data_bin,axis=(1,2)), drawstyle='steps-mid', lw=3,
				 #~ color='black')
		#~ plt.plot(waves4, continuum_flux, drawstyle='steps-mid', lw=3, color='blue')
		#~ plt.plot(wave_bin, flux_bin, 'ro', ms = 5)
		#~ plt.xlabel(r'Wavelength ($\AA$)')
		#~ plt.ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$ $\AA ^{-1}$)')

		ind_j += 1

	ind_i += 1

# Plot total continuum flux image
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,2,1)
ax1.set_xlabel(r'$\alpha$ ($^{\circ}$)',weight='bold',size='large')
ax1.set_ylabel(r'$\delta$ ($^{\circ}$)',weight='bold',size='large')
plt.imshow(continuum_img, origin='lower',
						interpolation="none", cmap='CMRmap',#cmap='nipy_spectral',
						extent=[ra[0],ra[-1],dec[0],dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'Total continuum flux (erg cm$^{-2}$ s$^{-1}$)',weight='bold',size='large')
ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)


# Plot chi2 statistics map of continuum fit to data
#fig1 = plt.figure()
ax1 = fig1.add_subplot(1,2,2)
ax1.set_xlabel(r'$\alpha$ ($^{\circ}$)',weight='bold',size='large')
ax1.set_ylabel(r'$\delta$ ($^{\circ}$)',weight='bold',size='large')
plt.imshow(chi2, origin='lower',
						interpolation="none", cmap='nipy_spectral',
						extent=[ra[0],ra[-1],dec[0],dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$\chi ^{2}$ statistic: continuum fit to data',weight='bold',size='large')
ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)


plt.show()
