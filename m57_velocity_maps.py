'''
Ring Nebula				170412		210			Small		BH2, Hbeta		60		N
(M57, NGC6720)						211			Small		BH2, Hbeta		60		N
									212			Small		BL,4500			20		N
									213			Small		BL,4500			300		N
						170415		188			Medium		BM,5900			10		N
						170619		259			Large		BM,4550			1050	Y? (DON'T USE)
									260			Large		BM,4550			1050	Y?
									261			Large		BM,4550			950		Y?
						170620		64			Medium		BM,5200			157		N?
Units of flux-cal data: erg/s/cm^2/Angstrom
'''
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.ndimage.filters import gaussian_filter
import scipy.ndimage
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


def subtract_continuum(waves,flux):
	'''
	Fit a line through the continuum around the region of interest.
	'''
	inds = [0,1,2,3,4,-5,-4,-3,-2,-1]
	wv_fit = numpy.zeros( len(inds) )
	fl_fit = numpy.zeros( len(inds) )
	for i in range(0,len(inds)):
		wv_fit[i] = waves[inds[i]]
		fl_fit[i] = flux[inds[i]]
	#~ wv_fit = waves[0:4,-5:-1] #numpy.array( [waves[0:5],waves[-6:-1]] )
	#~ numpy.reshape(wv_fit,1)
	#~ fl_fit = flux[0,1,2,3,4:-5,-4,-3,-2,-1] #numpy.array( [flux[0:5],flux[-6:-1]] )
	#~ print wv_fit, fl_fit
	poly_coeffs = poly.polyfit(wv_fit, fl_fit, 1)		# find best coefficients through continuum values
	cfit = poly.polyval(waves,poly_coeffs)

	# Subtract cfit from flux
	subtracted_flux = flux - cfit

	# Plot the flux, continuum fit, and sutracted flux
	#~ fig1 = plt.figure()
	#~ ax1 = fig1.add_subplot(1,2,1)
	#~ ax1.set_ylabel('Flux')
	#~ plt.plot(waves,flux,drawstyle='steps-mid',color='black',lw=3)
	#~ plt.plot(waves,cfit,color='red',lw=3)
	#~ ax2 = fig1.add_subplot(1,2,2)
	#~ ax2.set_ylabel('Flux')
	#~ ax2.set_xlabel(r'Wavelength ($\AA$)')
	#~ plt.plot(waves, subtracted_flux, drawstyle='steps-mid', color='black', lw=3)
	#~ plt.show()
	return subtracted_flux



def sn_cut(data_in, var, sigma):
	'''
	Defines the signal-to-noise ratio (SN) cut of the narrowband image!
	'''
	#~ snr = data_in/numpy.sqrt(var)
	data_sn = numpy.zeros( [numpy.size(data_in[:,0,0]), numpy.size(data_in[0,:,0]), numpy.size(data_in[0,0,:])] )
	#~ data_sn[snr>=sigma] = 0.0
	for i in range(0,numpy.size(data_in[:,0,0])):
		for j in range(0,numpy.size(data_in[0,:,0])):
			for k in range(0,numpy.size(data_in[0,0,:])):
				snr = data_in[i,j,k]/numpy.sqrt(var[i,j,k])	#numpy.sqrt(data_in[i,j,k])
				if snr > sigma:
					data_sn[i,j,k] = data_in[i,j,k]
				else:
					continue

	return data_sn



def define_amplitudes(lines, waves, spectra_0, dwave):
	'''
	Finds the amplitube (max) of the emission line from the line list provided
	by finding the line in the spectrum and determining the max value of the line.
	Provides an intial guess of the amplitude for the gaussian fitting routine.
	'''
	i_ww = 0
	lineamp = numpy.zeros( numpy.size(lines) )

	for ww in lines:
		i_line = numpy.where((waves <= ww+1) & (waves >= ww-1))[0]
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



def gaussfit(parms, x0, y0):
	'''
	This function does all the heavy (time consuming) work:
	Does a try/except:
		- try curve_fit optimization first.
		- if that doesn't work (usually gets to max iterations without converging):
			- except to leastsq() optimization, which should at least provided
			  a "best-guess" set of parms from the last iteration
	'''
	#~ plt.plot(x0,y0,drawstyle='steps-mid',color='blue')
	#~ plt.show()

	#~ try:
		#popt, pcov = curve_fit(gaussum, waves4, spectra_0, p0=gauss_parms)#, bounds=(bounds_lower,bounds_upper))#, sigma = err)
		#~ popt, pcov = curve_fit(gaussum, x0, y0, p0=parms, maxfev = 15000)
		#~ print popt.reshape((len(popt)/3,3))
		#~ print pcov
		#~ print
	#~ except RuntimeError:
		#~ print "RuntimeError encountered: trying leastsq"
		#~ #popt = curve_fit(gaussum, x0, y0, p0=parms, maxfev = 100000000)
		#~ #popt, pcov, infodict, errmsg, ier = curve_fit(gaussum, x0, y0, p0=parms, full_output=True)		# nope
	func = lambda p,x,y: gaussum(x,*p) - y
	popt, pcov, infodict, errmsg, ier = leastsq(func,parms,args=(x0,y0),full_output=True, maxfev = 5000)
	#~ print popt.reshape((len(popt)/3,3))
	#~ print pcov
		#~ print infodict
		#~ print errmsg, ier
		#~ print
	return popt


def velocity_info_doublet(parms, cen, c):
	'''
	Determine the velocity and dispersion (dynamics) of each line per pixel,
	based on the emission line info extracted from the gaussian fits.
	'''
	vr = numpy.zeros( numpy.size(cen) )			# km/s
	vdisp = numpy.zeros( numpy.size(cen) )		# convert to km/s
	for l in range(0, len(parms[:,0])):
		vr[l] = ((parms[l,1] - cen[l]) / cen[l])*c
		vdisp[l] = abs(2.3548*parms[l,2])*(c/cen[l])		# in wavelength (Ang) units --> km/s
		if (vr[l]<-500) or (vr[l]>500):
			vr[l] = -numpy.inf
		if vdisp[l]>65:
			vdisp[l] = -numpy.inf	#0
	if ((vr[1] - vr[0]) < -75) or (abs(vr[0]-vr[1]) > 75):
		vr = -numpy.inf	#0
	if abs(vdisp[0] - vdisp[1]) > 30:
		vdisp[0] = vdisp[1] = numpy.min(vdisp)	#0
	return vr, vdisp




def main_ring():
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
	Use to look at 06/20 data set (with [NI] ONLY)
	'''
	############### Main Ring (inner edge) ##################
	# Line list
	NI = [5197.902, 5200.257]
	ID = ['[NI]','[NI]']
	c = 3.0e5		# km/s

	############### Central Region 1: 06/20 ##################
	date='170620'

	int = 'icuber'
	var = 'vcuber'

	index1 = 64		# Spectral Range: ~ 5100 - 5300 (only [NI])

	intfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,int)
	file1 = path+dir+date+dir+redux+dir+intfile1
	varfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,var)
	vfile1 = path+dir+date+dir+redux+dir+varfile1

	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)

	all_ra = (all_ra - ra)*3600.
	all_dec = (all_dec - dec)*3600.

	# RA/DEC of central star(s)
	#~ ra_cs1 =
	#~ dec_cs1 =

	# Do things here
	# 1. Define emission line(s) to define velocity maps
	lines = NI
	dlam = 7.5	#10	#5. 		# +/- extra wavelength coverage from min/max of lines
	wv_range = numpy.where((waves1 > (numpy.min(lines)-dlam)) & (waves1 < (numpy.max(lines)+dlam)))[0]
	#~ wv_range = numpy.where((waves1 > 5060) & (waves1 < 5300))[0]
	waves_lines = waves1[wv_range]
	data_lines = data1[wv_range,:,:]
	var_lines = var1[wv_range,:,:]

	flux_nocont = numpy.zeros( numpy.shape(data_lines) )

	# 2. Loop through each image pixel to do continuu subtraction
	for i in range(0,numpy.size(data_lines[0,:,0])):
		for j in range(0,numpy.size(data_lines[0,0,:])):
			flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

	# 3. Define S/N ration cut
	sigma = 4.5			# 2.75 - good with var
	sigma = 2.5#2.2#33
	data_cut = sn_cut(flux_nocont, var_lines, sigma)
	#~ data_cut_contours = scipy.ndimage.zoom(numpy.sum(data_cut,axis=0), 3)
	data_cut_contours = gaussian_filter(numpy.sum(data_cut,axis=0), 0.9)

	# 4. Fit Gaussian profile(s) to line(s) of interest
	velocity_lines = numpy.zeros( [len(lines),numpy.size(data_cut[0,:,0]),numpy.size(data_cut[0,0,:])] )
	vdisp_lines = numpy.zeros( [len(lines),numpy.size(data_cut[0,:,0]),numpy.size(data_cut[0,0,:])] )
	for i in range(0,numpy.size(data_lines[0,:,0])):
		for j in range(0,numpy.size(data_lines[0,0,:])):
			if numpy.sum(data_cut[:,i,j]) <= 0:
				print "No lines detected at [%i,%i]" % (i,j)
				for ww in range(0,len(lines)):
					velocity_lines[ww,i,j] = -numpy.inf
					vdisp_lines[ww,i,j] = -numpy.inf
				continue
			else:
				# Define amplitude, st. dev, and parameter array for Gaussfit
				dwave = 1
				amp_lines = define_amplitudes(lines, waves_lines, data_cut[:,i,j], dwave)
				linestdv = numpy.ones( numpy.size(waves_lines) )*0.5 	# Std. dev. array w/ guess:
				gauss_parms = define_parms(lines, amp_lines, linestdv)

				# fit gaussian profile(s) to the line(s)
				# get the list of best-fit parameters back
				popt = gaussfit(gauss_parms, waves_lines, data_cut[:,i,j])
				popt = popt.reshape((len(popt)/3,3))

				# Determine velocity (v_r) and dispersion (fwhm) of each line
				vel, disp = velocity_info_doublet(popt, lines, c)
				velocity_lines[:,i,j] = vel
				vdisp_lines[:,i,j] = disp

	#print all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]
	#~ velocity_lines = gaussian_filter(velocity_lines, 0.1)
	#~ vdisp_lines = gaussian_filter(vdisp_lines, 0.1)
	print numpy.shape(velocity_lines)

	fig1 = plt.figure()
	ax1 = fig1.add_subplot(1,2,1)
	#~ plt.plot(waves_lines,numpy.sum(data_lines,axis=(1,2)),drawstyle='steps-mid')
	#~ plt.show()
	#~ levels=[20,30,40,50]
	#~ lw=[1, 2, 3, 4]
	levels=[25]
	lw=[2.5]
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)',weight='bold',size='large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\circ}$)',weight='bold',size='large')
	CS = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(velocity_lines[0,:,:], cmap='bwr', origin='lower',
			   vmin = -100, vmax = 100, #interpolation='bicubic',
			   extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ cbar = plt.colorbar()
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)

	#~ fig2 = plt.figure()
	ax2 = fig1.add_subplot(1,2,2)
	ax2.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)',weight='bold',size='large')
	#~ ax2.set_ylabel(r'$\delta$ ($^{\circ}$)',weight='bold',size='large')
	CS = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(velocity_lines[1,:,:], cmap='bwr', origin='lower',
			   vmin = -100, vmax = 100, #interpolation='bicubic',
			   extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ plt.imshow(numpy.sum(data_ni,axis=0), cmap='Blues', alpha=0.8, origin='lower',
				#~ extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	cbar = plt.colorbar()		# Each image should have its own colorbar
	cbar.ax.set_ylabel(r'Velocity (km/s)',
					   weight='bold',size='large')
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	#plt.show()



	fig1 = plt.figure()
	ax1 = fig1.add_subplot(1,2,1)
	#~ plt.plot(waves_lines,numpy.sum(data_lines,axis=(1,2)),drawstyle='steps-mid')
	#~ plt.show()
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)',weight='bold',size='large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\circ}$)',weight='bold',size='large')
	CS = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(vdisp_lines[0,:,:], cmap='BuGn', origin='lower',
			   vmin = 0, vmax = 50,
			   extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ cbar = plt.colorbar()
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)

	#~ fig2 = plt.figure()
	ax2 = fig1.add_subplot(1,2,2)
	ax2.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)',weight='bold',size='large')
	CS = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(vdisp_lines[1,:,:],cmap='BuGn', origin='lower',
			   vmin = 0, vmax = 50,
			   extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ plt.imshow(numpy.sum(data_ni,axis=0), cmap='Blues', alpha=0.8, origin='lower',
				#~ extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	cbar = plt.colorbar()		# Each image should have its own colorbar
	cbar.ax.set_ylabel(r'Velocity Dispersion (km/s)',
					    weight='bold',size='large')
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	#plt.show()



	fig2 = plt.figure()
	ax2 = fig2.add_subplot(1,1,1)
	ax2.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)',weight='bold',size='large')
	ax2.set_ylabel(r'$\Delta \delta$ ($^{\circ}$)',weight='bold',size='large')
	CS = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(numpy.sum(data_cut,axis=0), cmap='Blues', origin='lower',	#cmap='gnuplot'
				#vmin = 0, vmax = 100,
				extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	cbar = plt.colorbar()		# Each image should have its own colorbar
	cbar.ax.set_ylabel(r'[NI] Intensity (arb. units)',
					    weight='bold',size='large')
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	plt.show()


	return



def central_star_offset_medium():
	'''
	Use to look at 04/15 data set (with [OI], maybe [NII] around stars,
	and [OI], [NII], [OIII], HeII, and maybe other lines? in Inner Ring)
	'''
	############### Main Ring (inner edge) ##################
	# Line list
	HeII = [5411]
	NII = [5755]
	OI = [5577,6300]
	OIII = [5592.37]
	ClIII = [5517, 5537]

	############### Central Region 1: 06/20 ##################
	date='170415'

	int = 'icubes'
	var = 'vcubes'

	index1 = 188		# not ready

	intfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,int)
	file1 = path+dir+date+dir+redux+dir+intfile1
	varfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,var)
	vfile1 = path+dir+date+dir+redux+dir+varfile1


	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)

	# Do things here



	return



def central_star_offset_large():
	'''
	Use to look at 06/19 data set of central stars (no clouds) and Inner Ring
	'''
	############### Main Ring (inner edge) ##################
	# Line list
	NI = [5198.5, 5199]

	############### Central Region 1: 06/20 ##################
	date='170415'

	int = 'icubes'
	var = 'vcubes'

	index1 = 188		# not ready

	intfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,int)
	file1 = path+dir+date+dir+redux+dir+intfile1
	varfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,var)
	vfile1 = path+dir+date+dir+redux+dir+varfile1


	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)

	# Do things here



	return





if __name__ == "__main__":
	global path, dir, redux

	#path='C:\\Users\\Keri Hoadley\\Documents\\KCWI'
	#dir = '\\'
	path='/home/keri/KCWI/'
	dir = '/'
	redux='redux'

	# Line list
	HI = [4101.73, 4340.47, 4861.35]  	# delta, gamma, beta
	HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
	HeII = [4540, 4686, 5411]
	NI = [5198.5, 5199]
	OI = [5577,6300]
	OII = [3726, 3729, 4267, 4661]
	OIII = [4363.21, 4958.91, 5006.84, 5592.37]


	# Run each function per region probed
	central_star() #path,dir,redux)
