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

import os
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
	#~ data = data[:,33:-28,12:-12]		#400:-300
	#~ all_wv = all_wv[400:-300]
	#~ all_wv = all_wv + (all_wv[400] - all_wv[0])

	if len(args) > 0:
		ra0 = hdulist[0].header['TARGRA']#['CD1_1']		# RA (at center?)
		dec0 = hdulist[0].header['TARGDEC']#['CD2_2']		# DEC (at center?)
	else:
		ra0 = hdulist[0].header['RA']		# RA (at center?)
		dec0 = hdulist[0].header['DEC']		# DEC (at center?)

	#~ ra0_deg = hdulist[0].header['CRVAL1']		# RA (at center?)
	#~ dec0_deg = hdulist[0].header['CRVAL2']		# DEC (at center?)

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
	print ra0_deg, delta_ra, ra_lim
	print all_ra

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



def continuum_fit_lsq(wave_bin, flux_bin, waves, flux, err, wv_range):
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
	max_index = 15
	ls = numpy.zeros( max_index )
	#~ ls_min = 0
	#~ continuum_array = numpy.array()
	for i in range(1,max_index):		# increase in i = increase to polynomial factor
		poly_coeffs = poly.polyfit(wave_bin[1:-1], flux_bin[1:-1], i)		# find best coefficients through continuum values
		#cfit = poly.Polynomial(poly_coeffs)			# leave it in functional form (do cfit(waves))
		cfit = poly.polyval(waves,poly_coeffs)	# fit coefficients through wavelength of spectrum

		#~ fig1 = plt.figure()
		#~ ax1 = fig1.add_subplot(1,1,1)
		#~ plt.semilogy(wave_bin[1:-1],flux_bin[1:-1],drawstyle='steps-mid',color='black',lw=3)
		#~ plt.semilogy(waves,flux,drawstyle='steps-mid',color='blue',lw=3)
		#~ plt.semilogy(waves,cfit,color='red',lw=3)
		#~ plt.axvline(waves[300],color='red')
		#~ plt.axvline(waves[-300],color='red')
		#~ plt.show()
		# Perform least-squares between cfit and flux:
		# 							[ (flux[i] - cfit[i])**2 ]
		#				LS = SUM|----------------------- |
		#								[       err[i]**2				 ]
		for k in range(800,numpy.size(waves)-900):		#wv_range[0],wv_range[-1]):
			if flux[k] <= 0 or numpy.isnan(flux[k]):
				continue
			else:
				ls[i] += (flux[k] - cfit[k])**2 / (err[k]**2)
		ls[i] = ls[i]/(numpy.size(waves)-800) # - numpy.size(poly_coeffs))

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



def calc_EW(waves, i_Hb, EW_Hb, continuum, spectrum, w_Hb):
	'''
	Find the total area of continuum needed to match the total flux in Hbeta.
	Should output d_continuum (in wavelength, Angstroms) as one scalar.
	'''
	fcont = numpy.mean(continuum[w_Hb])
	#~ print EW_Hb/(fcont), 1400. - EW_Hb/(fcont)
	EW = 1400. - EW_Hb/(fcont)
	# Loop through +/- 1 dwave through the continuum flux
	#~ for i in range(1,numpy.size(waves),1):
		#~ # Find total flux in continuum between waves[i_Hb]+/-i
		#~ fcont = numpy.sum(continuum[i_Hb-i:i_Hb+i])
		#~ if i == 1:
			#~ print EW_Hb/(fcont/3)
			#~ EW = EW_Hb/(fcont/3)
			#~ break
		#~ # Check if fcont >= EW_Hb. If so, store dwave and break from FOR loop
		#~ if fcont*(waves[i_Hb+i] - waves[i_Hb-i]) >= EW_Hb*waves[i_Hb]:
			#~ EW = waves[i_Hb+i] - waves[i_Hb-i]
			#~ break
		#~ else:
			#~ continue
	#~ print EW

	# Plot the spectrum to show areas of Hb and continuum
	#~ fig1 = plt.figure()
	#~ ax1 = fig1.add_subplot(1,1,1)
	#~ plt.plot(waves, spectrum, drawstyle='steps-mid', color='black', lw=5)
	#~ plt.plot(waves[w_Hb], spectrum[w_Hb], color='red', lw=3)
	#~ plt.plot(waves[i_Hb-i:i_Hb+i], continuum[i_Hb-i:i_Hb+i],color='blue')
	#~ plt.xlim(waves[i_Hb]-200,waves[i_Hb]+200)
	#~ plt.show()

	return EW





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
path='C:\\Users\\Keri Hoadley\\Documents\\KCWI'
dir = '\\'
#~ path='/home/keri/KCWI'
#~ dir = '/'
redux='redux'
int = 'icuber'
var = 'vcuber'

#~ date='170412'
#~ date='170415'
date='170620'

#~ index1=210		# not ready
#~ index2=211		# not ready
#~ index3=212		# for OII, OIII lines analysis
#~ index4=213		# all other lines, use this higher S/N data set
#~ index4=188		# all other lines, use this higher S/N data set
index4=64		# all other lines, use this higher S/N data set

#~ intfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,int)		#_extcorr.
#~ intfile3 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index1,index2,int)		#_extcorr.
#~ file3 = path+dir+date+dir+redux+dir+intfile3
#~ varfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,var)
#~ varfile3 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index1,index2,var)
#~ vfile3 = path+dir+date+dir+redux+dir+varfile3


intfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,int)		#_extcorr
#~ intfile4 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,int)		#_extcorr
file4 = path+dir+date+dir+redux+dir+intfile4
varfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,var)
#~ varfile4 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,var)
vfile4 = path+dir+date+dir+redux+dir+varfile4

# Write new datafile (no continuum) out to new file: (at end of program)
newfile = file4.replace('.fits','.nocontinuum.fits')
newvfile = vfile4.replace('.fits','.nocontinuum.fits')

# Read in file, get data + waves from files
#~ data3, waves3, ra, dec = open_fits(file3)		# data in units erg cm-2 s-1
data4, waves4, ra, dec = open_fits(file4)		# data in units erg cm-2 s-1
#~ var3, varwv3, vra, vdec = open_fits(vfile3)
var4, varwv4, vra, vdec = open_fits(vfile4)

# Open emission line template file:
# Important columns in file: Line ID (col 0), Lab Wave (col 2)
#~ linesfile = date+"_212-213_lines.dat"
#~ linesfile = date+"_210-211_lines.dat"

# Open DAT file with emission line info
#~ lineID, linewave = open_lineID(path+dir+linesfile)

# For plotting spectrum fits later
#~ wave_split_index = [numpy.where(waves4 == 3600.)[0][0],
										#~ numpy.where(waves4 == 4000.)[0][0],
										#~ numpy.where(waves4 == 4400.)[0][0],
										#~ numpy.where(waves4 == 4800.)[0][0],
										#~ numpy.where(waves4 == 5200.)[0][0],
										#~ numpy.where(waves4 == 5600.)[0][0]]

# Change to data3 file outputs
#~ data4 = data3
#~ var4 = var3
#~ waves4 = waves3
#~ file4 = file3
#~ vfile4 = vfile3
#~ newfile = file4.replace('.fits','.nocontinuum.fits')
#~ newvfile = vfile4.replace('.fits','.nocontinuum.fits')

data4[numpy.isnan(data4)] = 0
data4[data4<0] = 0	#numpy.average(data4)
data4[data4>100*numpy.average(data4)] = numpy.average(data4)

# Define regioin between wave_good_start and wave_good_end
# will change per data set, so keep track of this!
# 170415
wave_good_start = 5500	#5060
wave_good_end = 6350	# 5300
# 170620
wave_good_start = 5060
wave_good_end = 5300
wv_range = numpy.where((waves4 > wave_good_start) & (waves4 < wave_good_end))[0]
print wv_range
#~ waves4 = waves4[wv_range]
#~ data4 = data4[wv_range,:,:]
#~ var4 = var4[wv_range,:,:]


# Find where Hbeta is to
#~ Hb = linewave[lineID.index('Hbeta')]
#~ wv_Hb = numpy.where( (waves4 <= Hb+10.) & (waves4 >= Hb-10.) )[0]
#~ cen_Hb = numpy.where( (waves4 <= Hb+0.3) & (waves4 >= Hb-0.3) )[0]
#~ cen_Hb = cen_Hb[0]

# Define image "bins" to loop through and find spectra,
# continuum level(s), and emission line characteristics.
# Image size is: (wave, 132L, 22L)
# TEST BINS:
ra_bin = 1 #numpy.size(data4[0,0,:])/2			# 2
dec_bin = 1	#33 #numpy.size(data4[0,:,0])/2		# 33
#REAL BINS:
#~ ra_bin = 2
#~ dec_bin = 2

bin_size = 20		# Bin size of continuum determination in each spectrum
#~ print 'RA bin: ', ra_bin
#~ print 'DEC bin: ', dec_bin
#~ print

# Make a continuum flux 2D array with size of total # of elements in ra, dec_bin
continuum_img = numpy.zeros( [numpy.size(data4[:,0,0]),
															numpy.size(data4[0,:,0])/dec_bin,
															numpy.size(data4[0,0,:])/ra_bin] )
chi2 = numpy.zeros( [numpy.size(data4[0,:,0])/dec_bin,
										 numpy.size(data4[0,0,:])/ra_bin] )
# Create an EW array, which compares the total flux in Hbeta to the total
# continuum area needed to match the flux in Hbeta.
# Gives a measure of the temperature of (cooler?) gas.
#~ EW_continuum = numpy.zeros( [numpy.size(data4[0,:,0])/dec_bin,
										 #~ numpy.size(data4[0,0,:])/ra_bin] )




# Plot data image before and after continuum subtraction
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,2,1)
ax1.set_xlabel(r'$\alpha$ ($^{\circ}$)',weight='bold',size='large')
ax1.set_ylabel(r'$\delta$ ($^{\circ}$)',weight='bold',size='large')
plt.imshow(numpy.sum(data4[wv_range,:,:],axis=0), origin='lower',
			interpolation="none", cmap='CMRmap',#cmap='nipy_spectral',
			extent=[ra[0],ra[-1],dec[0],dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'Total continuum flux (erg cm$^{-2}$ s$^{-1}$)',weight='bold',size='large')
ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)
#plt.show()

print ra[0], ra[-1]
ind_i = 0

# Loop through ra, dec bin sizes and do spectral operations
for i in range(0,numpy.size(data4[0,:,0]),dec_bin):
	ind_j = 0
	for j in range(0,numpy.size(data4[0,0,:]),ra_bin):
		#~ print 'Analyzing for pixels: RA [%i, %i], DEC [%i, %i]' %
		#~ 				(j, j+ra_bin, i, i+dec_bin)

		#~ data_bin = data4[:,i:i+dec_bin,j:j+ra_bin]
		#~ var_bin = var4[:,i:i+dec_bin,j:j+ra_bin]

		data_bin = data4[:,i,j]
		var_bin = var4[:,i,j]

		# Create spectra (Collapse along x,y axes)
		spectra = data_bin#create_spectum(data_bin)
		err = numpy.sqrt(abs(var_bin)) #create_spectum_err(var_bin)

		# Find continuum fit:
		# 1. Define a bin size to find medium flux through and extrac new wave, flux
		#			arrays for continuum
		wave_bin, flux_bin = fit_continuum_to_spectrum(spectra[wv_range], waves4[wv_range], bin_size)
		#
		# 2. Use new continuum (flux) array to fit polynomial through entire
		# 		wavelength range.
		#			Return continuum array and least-squares min chi2 value of fit to data.
		continuum_flux, chi2[ind_i,ind_j] = continuum_fit_lsq(wave_bin, flux_bin, waves4, spectra, err, wv_range)
		continuum_img[:,ind_i,ind_j] =  continuum_flux
		print 'At [%i - %i, %i - %i] -- Total continuum flux: %.3e   Chi2 = %.3f' % (i, i+dec_bin, j, j+ra_bin, numpy.sum(continuum_flux),chi2[ind_i,ind_j])


		# 3. Subtract continuum flux from spectrum
		spectra_0 = subtract_continuum(spectra, continuum_flux)
		# Check for nan's and get right of them
		for w in range(0,len(spectra_0)):
			if spectra_0[w] <= 0 or numpy.isnan(spectra_0[w]):
				spectra_0[w] = abs(numpy.nanmedian(spectra_0))
				err[0] = 0

		# Use continuum_flux array to find delta_wave where EW_continuum == EW_Hb
		#~ EW_Hb = numpy.sum(spectra_0[wv_Hb])		# Total flux in Hb line
		#~ EW_continuum[ind_i,ind_j] = calc_EW(waves4, cen_Hb, EW_Hb, continuum_flux, spectra, wv_Hb)


		#### PLOTTING OPTIONS (to check what continuum subtraction/create is doing) #####
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


		#~ # Plot the before and after continuum subtracted
		#~ fig1 = plt.figure()
		#~ ax2 = fig1.add_subplot(1,1,1)
		#~ plt.plot(waves4, spectra, drawstyle='steps-mid', lw=3,
				 #~ color='black')
		#~ plt.plot(waves4, spectra_0, drawstyle='steps-mid', lw=3, color='blue')
		#~ plt.axhline(0, linestyle='dashed', color='red', lw=3)
		#~ plt.xlabel(r'Wavelength ($\AA$)')
		#~ plt.ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$ $\AA ^{-1}$)')



		# (save to separate fits file - after end of loop)
		# take the average of continuum flux over bin sizes
		for m in range(i,i+dec_bin):
			for l in range(j,j+ra_bin):
				data4[:,m,l] = data4[:,m,l] - continuum_flux/(dec_bin*ra_bin)

		ind_j += 1

	ind_i += 1

# Again, make sure no parts of data4 are nan or left out
data4[numpy.isnan(data4)] = 0
data4[numpy.isinf(data4)] = 0
data4[data4<0] = 0
var4[numpy.isnan(var4)] = 0
var4[numpy.isinf(var4)] = 0
var4[var4<0] = 0
continuum_img[numpy.isnan(continuum_img)] = 0



# Plot after continuum subtraction image of data4
ax1 = fig1.add_subplot(1,2,2)
ax1.set_xlabel(r'$\alpha$ ($^{\circ}$)',weight='bold',size='large')
ax1.set_ylabel(r'$\delta$ ($^{\circ}$)',weight='bold',size='large')
plt.imshow(numpy.sum(data4[wv_range,:,:],axis=0), origin='lower',
			interpolation="none", cmap='CMRmap',# cmap='nipy_spectral',
			extent=[ra[0],ra[-1],dec[0],dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'Total continuum flux (erg cm$^{-2}$ s$^{-1}$)',weight='bold',size='large')
ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)



# Plot total continuum flux image + chi2 fit
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,2,1)
ax1.set_xlabel(r'$\alpha$ ($^{\circ}$)',weight='bold',size='large')
ax1.set_ylabel(r'$\delta$ ($^{\circ}$)',weight='bold',size='large')
plt.imshow(numpy.sum(continuum_img[300:-400],axis=0), origin='lower',
						interpolation="none", cmap='CMRmap',#cmap='nipy_spectral',
						extent=[ra[0],ra[-1],dec[0],dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'Total continuum flux (erg cm$^{-2}$ s$^{-1}$)',
									 weight='bold',size='large')
ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)

# Plot chi2 statistics map of continuum fit to data
#fig1 = plt.figure()
ax1 = fig1.add_subplot(1,2,2)
ax1.set_xlabel(r'$\alpha$ ($^{\circ}$)',weight='bold',size='large')
ax1.set_ylabel(r'$\delta$ ($^{\circ}$)',weight='bold',size='large')
plt.imshow(chi2, origin='lower', vmin = 0, #vmax = 250,
						interpolation="none", cmap='nipy_spectral',
						extent=[ra[0],ra[-1],dec[0],dec[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$\chi ^{2}$ statistic: continuum fit to data',
									 weight='bold',size='large')
ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)



#~ # Plot the equivalent width measurement through the ring
#~ fig1 = plt.figure()
#~ ax1 = fig1.add_subplot(1,1,1)
#~ ax1.set_xlabel(r'$\alpha$ ($^{\circ}$)',weight='bold',size='large')
#~ ax1.set_ylabel(r'$\delta$ ($^{\circ}$)',weight='bold',size='large')
#~ plt.imshow(EW_continuum, origin='lower',
						#~ interpolation="none", cmap='nipy_spectral', #cmap='CMRmap',#
						#~ extent=[ra[0],ra[-1],dec[0],dec[-1]])
#~ cbar = plt.colorbar()
#~ cbar.ax.set_ylabel(r'H$\beta$ Equivalent Width ($\AA$)',
									 #~ weight='bold',size='large')
#~ ax = plt.gca()
#~ ax.get_xaxis().get_major_formatter().set_useOffset(False)
#~ ax.get_yaxis().get_major_formatter().set_useOffset(False)


plt.show()


# Write new datafile (no continuum) out to new file:
write_fits(newfile,file4,data4)
# Write continuum to a new file
newcfile = file4.replace('.fits','.continuum.fits')
write_fits(newcfile,file4,continuum_img)
write_fits(newvfile,vfile4,var4)		# also save variance cube with correct size

