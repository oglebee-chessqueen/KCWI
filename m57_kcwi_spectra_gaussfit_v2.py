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
#matplotlib.rc('text', usetex=True)



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


def write_fits(file,oldfile,data):
	# First, read-in previous list of events
	try:
		print "Creating and writing fits file..."
		# grab header from original fits file
		hdulist_old = fits.open(oldfile)
		hdr = hdulist_old[0].header
		hdulist_old.close()
		# Open (create) new fits file to store header and continuum-subtracted datacube.
		hdu = fits.PrimaryHDU(data=data, header=hdr)
		hdulist = fits.HDUList([hdu])
		hdulist.writeto(file)
	except:
		print "File already exists - overwriting!"
		hdulist = fits.open(file,mode='update')
		hdulist[0].data = data

		hdulist.flush()		# Should update file with new events array
		hdulist.close()
		print "Updating %s SUCCESS!" % (file)

	return


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



def define_amplitudes(lines, waves, dwave):
	'''
	Finds the amplitube (max) of the emission line from the line list provided
	by finding the line in the spectrum and determining the max value of the line.
	Provides an intial guess of the amplitude for the gaussian fitting routine.
	'''
	i_ww = 0
	lineamp = numpy.zeros( numpy.size(lines) )

	for ww in lines:
		i_line = numpy.where((waves <= ww+0.25) & (waves >= ww-0.25))[0]
		#print ww, i_linewave
		amp = abs(numpy.max(spectra_0[i_line[0]-dwave:i_line[0]+dwave]))
		if numpy.isnan(amp):
			lineamp[i_ww] = lineamp[i_ww-1]
		else:
			lineamp[i_ww] = amp
		i_ww += 1

	return lineamp



def define_parms(cen, amp, stdv):
	'''
	Use lists of wavelength centers, line amplitudes, and stardard deviation (sigma)
	guesses to create the parameter list needed for the multi-gaussian fitting.
	'''
	guess = numpy.zeros(0)
	for k in range(len(amp)): 	 #This is where the guessed values are all packed up into one array looking like
		#~ xi = lineamp[k]         #[amp1,cen1,stdv1,amp2,cen2,...]
		#~ guess = numpy.append(guess,xi)
		#~ yi = linewave[k]
		#~ guess = numpy.append(guess,yi)
		#~ zi = linestdv[k]
		#~ guess = numpy.append(guess,zi)
		guess = numpy.append(guess,[amp[k],cen[k],stdv[k]])
	return guess



def velocity_info(parms, cen, ID, c):
	'''
	Determine the velocity and dispersion (dynamics) of each line per pixel,
	based on the emission line info extracted from the gaussian fits.
	'''
	vr = numpy.zeros( numpy.size(cen) )			# km/s
	fwhm = numpy.zeros( numpy.size(cen) )		# convert to km/s
	for l in range(0, len(parms[:,0])):
		if numpy.any(ID.index('Buffer')) == l:
			continue
		else:
			vr[l] = ((parms[l,1] - cen[l]) / cen[l])*c
			fwhm[l] = abs(2.3548*parms[l,2])*(c/cen[l])		# in wavelength (Ang) units --> km/s
	return vr, fwhm




def total_flux(waves, parms, ID):
	'''
	Calculate total flux in each emission line found in gaussum and
	return the array of total flux values (to be stored in main program).
	'''
	f_tot = numpy.zeros( len(ID) )			# ergs cm-2 s-1
	for l in range(0, len(parms[:,0])):
		if numpy.any(ID.index('Buffer')) == l:
			continue
		else:
			lineflux = gaus(waves,parms[l,0],parms[l,1],parms[l,2])
			f_tot[l] = numpy.sum(lineflux)
	return vr, fwhm




def plot_gaussfit_show(waves, spectra, fit, parms):
	'''
	SHOW plot of Gaussian fits over data,
	this way can manipulate and see how fits went.
	'''
	fig1 = plt.figure()
	ax2 = fig1.add_subplot(1,1,1)
	plt.plot(waves, spectra, drawstyle='steps-mid', lw=6,
					 color='black')
	plt.plot(waves, fit, color='green', drawstyle='steps-mid', lw=3)
	# also plot individual Gaussian fits
	# Use this loop to determine v_r, fwhm from Gaussian fits
	for l in range(0, len(parms[:,0])):
		gfit = gaus(waves,parms[l,0],parms[l,1],parms[l,2])
		plt.plot(waves, gfit, drawstyle='steps-mid', lw=2)
	plt.xlabel(r'Wavelength ($\AA$)')
	plt.ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$ $\AA ^{-1}$)')
	plt.show()
	return




def plot_gaussfit_save(waves, spectra, fit, parms, index, file):
	'''
	SAVE Gaussian fits of emission lines over data, defined at
	each pixel location.
	'''
	# Save a plot with subfigures cutting the spectrum up into chunks
	img_tot = len(index)-1
	fig1 = plt.figure(figsize=(11,10))
	for img in range(0, len(index)-1):
		i1 = index[img]-20
		i2 = index[img+1]+20
		# total number of figs = wave_split_index - 1
		ax2 = fig1.add_subplot(img_tot,1,img+1)
		# plot spectrum and gaussfit over spectrum
		plt.plot(waves[i1:i2], spectra[i1:i2]*1e15, drawstyle='steps-mid',
						 lw=5, color='black')
		plt.plot(waves[i1:i2], fit[i1:i2]*1e15, color='blue',
						 drawstyle='steps-mid', lw=3)
		for k in range(0, len(parms[:,0])):
			gfit = gaus(waves,parms[k,0],parms[k,1],parms[k,2])
			plt.plot(waves[i1:i2], gfit[i1:i2]*1e15, drawstyle='steps-mid', lw=1.5)
		plt.xlim( (waves[i1],waves[i2]) )
		# add axes at appropriate positions
		if img == (len(index)-2)/2:
			ax2.set_ylabel(r'Flux (10$^{-15}$ erg cm$^{-2}$ s$^{-1}$ $\AA ^{-1}$)',
										 weight='bold',size='large')
		if img == len(index)-2:
			ax2.set_xlabel(r'Wavelength ($\AA$)',weight='bold',size='large')
	plt.subplots_adjust(bottom=0.15, right=0.9, left=0.1, top=0.9)
	plt.savefig(file, orientation='landscape', format='ps')
	plt.close(fig1)		# clear figure
	return





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
	#~ global storage #You may not need storage to be global so think about taking this part out. storage stores the data
	storage = numpy.zeros( shape=[len(xdata),len(params)/3]) #[[0 for x in range(1)] for x in range(len(params)/3)] #from each iteration of the gaussian equation into
	for i in range(len(params)/3):#individual rows. So row one will be the gaussian solutions to the first peak and so on
		storage[:,i] = gaus(xdata,amp[i],cen[i],stdv[i])
	#storage = numpy.asarray(storage)

	# see what the plot looks like?
	#~ plt.plot(xdata, numpy.sum(storage, axis=1), drawstyle='steps-mid')
	#~ plt.show()

	return numpy.sum(storage, axis=1)






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
c = 3.0E5
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

intfile3 = 'kb'+date+'_00%03i_%s_extcorr.fits' % (index3,int)		#_extcorr.
#~ intfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,int)		#_extcorr.
file3 = path+dir+date+dir+redux+dir+intfile3
varfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,var)
vfile3 = path+dir+date+dir+redux+dir+varfile3

intfile4 = 'kb'+date+'_00%03i_%s_extcorr.fits' % (index4,int)		#_extcorr
#~ intfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,int)		#_extcorr
file4 = path+dir+date+dir+redux+dir+intfile4
varfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,var)
vfile4 = path+dir+date+dir+redux+dir+varfile4

# Write new datafile (no continuum) out to new file: (at end of program)
newfile = path+dir+date+dir+redux+dir+'kb'+date+'_00%03i_%s_gaussfit_parms.npz' % (index4,int)

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

# For plotting spectrum fits later
wave_split_index = [numpy.where(waves4 == 3600.)[0][0],
										numpy.where(waves4 == 4000.)[0][0],
										numpy.where(waves4 == 4400.)[0][0],
										numpy.where(waves4 == 4800.)[0][0],
										numpy.where(waves4 == 5200.)[0][0],
										numpy.where(waves4 == 5600.)[0][0]]

# Change to data3 file outputs
#~ data4 = data3
#~ var4 = var3
#~ waves4 = waves3
#~ newfile = path+dir+date+dir+redux+dir+'kb'+date+'_00%03i_%s_nocontinuum.fits' % (index3,int)


# Define image "bins" to loop through and find spectra,
# continuum level(s), and emission line characteristics.
# Image size is: (wave, 132L, 22L)
# TEST BINS:
ra_bin = numpy.size(data4[0,0,:])/2			# 2				 #
dec_bin = numpy.size(data4[0,:,0])/2		# 33 			 #
#REAL BINS:
ra_bin = 1
dec_bin = 1

bin_size = 20		# Bin size of continuum determination in each spectrum
#~ print 'RA bin: ', ra_bin
#~ print 'DEC bin: ', dec_bin
#~ print

# Make a continuum flux 2D array with size of total # of elements in ra, dec_bin
continuum_img = numpy.zeros( [numpy.size(data4[0,:,0])/dec_bin,
															numpy.size(data4[0,0,:])/ra_bin] )
chi2 = numpy.zeros( [numpy.size(data4[0,:,0])/dec_bin,
										 numpy.size(data4[0,0,:])/ra_bin] )

# Make arrays to hold Gauss-fit results: v_r, fwhm, and total line flux_bin
vr = numpy.zeros( [len(linewave),
									 numpy.size(data4[0,:,0])/dec_bin,
									 numpy.size(data4[0,0,:])/ra_bin] )
fwhm = numpy.zeros( numpy.shape(vr) )
flux_tot = numpy.zeros( numpy.shape(vr) )

# Make an additional array to store Gaussfit parameters per spatial pixels
# (to re-create Gaussian plots later or get more info as needed, without re-running)
g_parms = numpy.zeros( [len(linewave)*3,
												numpy.size(data4[0,:,0])/dec_bin,
												numpy.size(data4[0,0,:])/ra_bin] )


# Plot data image before and after continuum subtraction
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,2,1)
ax1.set_xlabel(r'$\alpha$ ($^{\circ}$)',weight='bold',size='large')
ax1.set_ylabel(r'$\delta$ ($^{\circ}$)',weight='bold',size='large')
plt.imshow(numpy.sum(data4,axis=0), origin='lower',
			interpolation="none", cmap='nipy_spectral',
			extent=[ra[0],ra[-1],dec[0],dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'Total continuum flux (erg cm$^{-2}$ s$^{-1}$)',weight='bold',size='large')
ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)


ind_i = 0

# Loop through ra, dec bin sizes and do spectral operations
for i in range(0,numpy.size(data4[0,:,0]),dec_bin):
	ind_j = 0
	for j in range(0,numpy.size(data4[0,0,:]),ra_bin):
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

		# 3. Subtract continuum flux from spectrum
		spectra_0 = subtract_continuum(spectra, continuum_flux)
		# Check for nan's and get right of them
		for w in range(0,len(spectra_0)):
			if spectra_0[w] <= 0 or numpy.isnan(spectra_0[w]):
				spectra_0[w] = abs(numpy.nanmedian(spectra_0))
				err[0] = 0

		# 4. Find amplitudes of each emission line in spectra_0
		# DO EACH ITERATION - if not and use "guess" from previous iteration,
		# will think chi2 fit is good enough and won't search parameter space
		# for better fits...
		lineamp = define_amplitudes(linewave, waves4, 4)			# dwave = 4, # of indices to look around linecenter for amplitude
		linestdv = numpy.ones( numpy.size(linewave) )*0.5 	# Std. dev. array w/ guess:
																												# linewidth = 1

		# 5. Now that linewave, lineamp, and linestdv exist, can call gaussum
		# through the curve_fit routine
		gauss_parms = define_parms(linewave, lineamp, linestdv)
		####### Set Boundaries (But doesn't work!)
		#~ bounds_lower = define_parms(numpy.array(linewave)-5, -5*lineamp, numpy.ones( numpy.size(linestdv) )*-numpy.inf)		# defines the lower limit bounds for curve_fit,
																																#~ # assuming the same for all elements
		#~ bounds_lower = tuple( bounds_lower )
		#~ bounds_upper = define_parms(numpy.array(linewave)+5, 5*lineamp, numpy.ones( numpy.size(lineamp) )*numpy.inf)
		#~ bounds_upper = tuple( bounds_upper )

		# This is where we take the gaussum function and use curvefit to find the best-fit
		# multi-gaussian function across the spectrum!
		# OPTIMIZATION!
		#~ #try:
		popt, pcov = curve_fit(gaussum, waves4, spectra_0, p0=gauss_parms)#, bounds=(bounds_lower,bounds_upper))#, sigma = err)

		fit = gaussum(waves4, *popt)	# Use the best-fit parms to create the emission line spectrum!

		#The optimized parameters for each emission peak is stored here in a 3xn array.
		g_parms[:,ind_i,ind_j] = popt
		popt = popt.reshape((len(popt)/3,3))
		print popt
		print pcov
		print

		# Determine velocity (v_r) and dispersion (fwhm) of each line
		vr[:,ind_i,ind_j], fwhm[:,ind_i,ind_j] = velocity_info(popt, linewave, lineID, c)

		# Determine total flux in Gaussian lines
		flux_tot[:,ind_i,ind_j] = total_flux(waves4, popt, lineID)

		#~ # plot the fit!
		#~ plot_gaussfit_show(waves4, spectra_0, fit, popt)

		#~ # plot to save to file!
		#~ figfile = path+dir+date+'_gaussfit_ra%02i%02i_dec%03i%03i.ps' % (j, j+ra_bin, i, i+dec_bin)
		#~ plot_gaussfit_save(waves4, spectra_0, fit, popt, wave_split_index, figfile)


		#~ except:
			#~ print "Gaussfit didn't work - trying scipy.optimize.leastsq routine instead"



		# (save to separate fits file - after end of loop)
		# take the average of continuum flux over bin sizes
		for m in range(i,i+dec_bin):
			for l in range(j,j+ra_bin):
				data4[:,m,l] = data4[:,m,l] - continuum_flux/(dec_bin*ra_bin)

		ind_j += 1

	ind_i += 1


# Save v_r, fwhm, tot_flux, g_parms arrays to numpy SAVE file
# Should  be kept in array format
# Save keywords:
# vr = v_r
# fwhm = fwhm
# lineflux = tot_flux
# parms = g_parms
numpy.savez(newfile, vr=v_r, fwhm=fwhm, lineflux=tot_flux, parms=g_parms)
