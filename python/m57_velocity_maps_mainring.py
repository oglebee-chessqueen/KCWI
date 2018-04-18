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
#~ from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import numpy
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import ks_2samp, ttest_ind, kde
import scipy.ndimage
plt.rcParams["font.family"] = "serif"
plt.rcParams['axes.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'



def save_outputs(ID, velocity, vdisp, region_ID, region):
	'''
	SAVE ONE ION AT A TIME! Can do multiple lines from same ion is desired
	(i.e., [NI]5197+5200, [OII]3726+3729, [OIII]4959,5007

	Use to save line IDs, velocity, and dispersion arrays to use for comparisons
	in either other programs or in the main function of this program.

	Collapse each velocity array to one dimension before saving.(?)
	'''
	file = "%s_kcwi_velocities_regions.npz" % (ID[0])
	numpy.savez(file, ID=ID, velocity=velocity, vdisp=vdisp, region_ID=region_ID, region=region)

	return "Arrays saved as: velocity, vdisp, region_ID, region"


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

	#~ ra0_deg = hdulist[0].header['CRVAL1']
	#~ dec0_deg = hdulist[0].header['CRVAL2']

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


def plot_velocity_double(velocity, contours, lvls, lwds, ra, dec, ID, line, min, max, color, cbar_label):
	'''
	Plot doublet velocity distribution(s), including contour regions.
	'''
	velocity[numpy.isinf(velocity)] = 0
	# Define contour levels
	levels = lvls
	lw = lwds
	fig1 = plt.figure(figsize=(8,8),facecolor='white')
	# Plot 1
	ax1 = fig1.add_subplot(1,2,1)
	ax1.set_facecolor('white')
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax1.set_title(ID[0]+'%i'%(line[0]),weight='bold',size='x-large')
	#~ ax1.set_title(ID[0],weight='bold',size='x-large')
	CS = plt.contour(contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[ra[0],ra[-1],dec[0],dec[-1]])
	im = plt.imshow(velocity[0,:,:], cmap=color, origin='lower',
			   vmin = min, vmax = max,
			   extent=[ra[0],ra[-1],dec[0],dec[-1]])
	#~ plt.scatter(ra_stars,dec_stars, s=1000, facecolor='lime', marker='*', #alpha=0.5,
							#~ edgecolor='limegreen', linewidths=2.5)
	ax = plt.gca()
	plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
	plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
	# Plot 2
	ax2 = fig1.add_subplot(1,2,2)
	ax2.set_facecolor('white')
	ax2.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax2.set_title(ID[1]+'%i'%(line[1]),weight='bold',size='x-large')
	#~ ax2.set_title(ID[1],weight='bold',size='x-large')
	CS = plt.contour(contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[ra[0],ra[-1],dec[0],dec[-1]])
	im = plt.imshow(velocity[1,:,:], cmap=color, origin='lower',
			   vmin = min, vmax = max,
			   extent=[ra[0],ra[-1],dec[0],dec[-1]])
	#~ plt.scatter(ra_stars,dec_stars, s=1000, facecolor='lime', marker='*', #alpha=0.5,
							#~ edgecolor='limegreen', linewidths=2.5)
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = fig1.colorbar(im, cax=cax)
	cbar.ax.set_ylabel(cbar_label, weight='bold', size='x-large')
	#~ plt.tight_layout()
	#~ plt.savefig('m57_%s_%s_wstars.pdf'%(ID[0],cbar_label[0:-7]),orientation='landscape')#,transparent=True)
	plt.show()
	return



def plot_velocity(velocity, contours, lvls, lwds, ra, dec, cont_ra, cont_dec, ID, line, min, max, color, cbar_label, *args):
	'''
	Plot doublet velocity distribution(s), including contour regions.
	Velocity - should only be 2D, of one line or avg. velocity of multiple lines
	'''
	velocity[numpy.isinf(velocity)] = 0

	levels = lvls
	lw = lwds
	fig = plt.figure(figsize=(8,8),facecolor='white')
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)', weight='bold', size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)', weight='bold', size='x-large')
	#~ ax1.set_title(ID[0], weight='bold', size='x-large')
	CS = plt.contour(contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[cont_ra[0],cont_ra[-1],cont_dec[0],cont_dec[-1]])
	im = plt.imshow(velocity, cmap=color, origin='lower',	#cmap='gnuplot'
					vmin = min, vmax = max,
					extent=[ra[0],ra[-1],dec[0],dec[-1]])
	#~ cbar = plt.colorbar()		# Each image should have its own colorbar
	#~ plt.scatter(ra_stars,dec_stars, s=1000, facecolor='lime', marker='*', #alpha=0.5,
							#~ edgecolor='limegreen', linewidths=2.5)
	ax = plt.gca()
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = fig.colorbar(im, cax=cax)
	cbar.ax.set_ylabel(cbar_label, weight='bold', size='x-large')
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)


	#~ if len(args) == 1:
		#~ plt.savefig('m57_%s_wstars.eps'%(args[0]))#,transparent=True)
	#~ else:
		#~ plt.savefig('m57_%s_%s_wstars.eps'%(ID[0],cbar_label[0:-7]))#,transparent=True)
	plt.show()
	return



def plot_intensity_map(data, contours, lvls, lwds, ra, dec, cont_ra, cont_dec, ID, line, min, max, color, cbar_label):
	'''
	Plot doublet velocity distribution(s), including contour regions.
	'''
	levels = lvls
	lw = lwds
	fig = plt.figure(figsize=(8,8))
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)', weight='bold', size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)', weight='bold', size='x-large')
	CS = plt.contour(contours, levels,
					 linewidths=lw, colors='w', corner_mask=True,
					 extent=[cont_ra[0],cont_ra[-1],cont_dec[0],cont_dec[-1]])
	im = plt.imshow(numpy.sum(data,axis=0), cmap=color, origin='lower',	#cmap='gnuplot'
				#vmin = 0, vmax = 100,
				extent=[ra[0],ra[-1],dec[0],dec[-1]])
	#~ cbar = plt.colorbar()		# Each image should have its own colorbar
	plt.scatter(ra_stars,dec_stars, s=1000, facecolor='lime', marker='*', #alpha=0.5,
							edgecolor='limegreen', linewidths=2.5)
	ax = plt.gca()
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = fig.colorbar(im, cax=cax)
	cbar.ax.set_ylabel(cbar_label, weight='bold', size='x-large')
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)

	plt.show()

	#~ plt.savefig('m57_%s_wstars.eps'%(cbar_label[0:-13]))#,transparent=True)
	return


def plot_intensity_map_labeled(data, contours, lvls, lwds, ra, dec, cont_ra, cont_dec, ID, line, min, max, color, cbar_label, regions, region_labels):
	'''
	Plot doublet velocity distribution(s), including contour regions.
	'''
	levels = lvls
	lw = lwds
	fig = plt.figure(figsize=(8,8))
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)', weight='bold', size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)', weight='bold', size='x-large')
	CS = plt.contour(contours, levels,
					 linewidths=lw, colors='w', corner_mask=True,
					 extent=[cont_ra[0],cont_ra[-1],cont_dec[0],cont_dec[-1]])
	im = plt.imshow(numpy.sum(data,axis=0), cmap=color, origin='lower',	#cmap='gnuplot'
				vmin = min, vmax = max,
				extent=[ra[0],ra[-1],dec[0],dec[-1]])
	# Plot region labels to reference in paper
	# For [NI] (central region)
	#~ plt.text(ra[regions[0,2]]+7.0, dec[regions[0,3]], region_labels[0],
					 #~ weight='bold',size='x-large', color='white')	# Cloud 1
	#~ plt.text(ra[regions[1,0]]+1.5, dec[regions[1,3]], region_labels[1],
					 #~ weight='bold',size='x-large', color='white')	# Cloud 2
	#~ plt.text(ra[regions[2,0]]+1.0, dec[regions[2,3]]+0.5, region_labels[2],
					 #~ weight='bold',size='x-large', color='white')	# Cloud 3
	#~ plt.text(ra[regions[3,2]]+0.5, dec[regions[3,3]]-0.5, region_labels[3],
					 #~ weight='bold',size='x-large', color='white')	# Cloud 4
	# For [OI] (central_region_medium)
	plt.text(ra[regions[0,2]]+8.0, dec[regions[0,3]]-1.0, region_labels[0],
					 weight='bold',size='x-large', color='white')	# Cloud 1
	plt.text(ra[regions[1,0]], dec[regions[1,3]]+0.5, region_labels[1],
					 weight='bold',size='x-large', color='white')	# Cloud 2
	plt.text(ra[regions[2,0]]+4.7, dec[regions[2,3]]-2.0, region_labels[2],
					 weight='bold',size='x-large', color='white')	# Inner Ring

	plt.scatter(ra_stars,dec_stars, s=1000, facecolor='lime', marker='*', #alpha=0.5,
							edgecolor='limegreen', linewidths=2.5)
	#~ cbar = plt.colorbar()		# Each image should have its own colorbar
	ax = plt.gca()
	#~ divider = make_axes_locatable(ax)
	#~ cax = divider.append_axes("right", size="5%", pad=0.05)
	#~ cbar = fig.colorbar(im, cax=cax)
	#~ cbar.ax.set_ylabel(cbar_label, weight='bold', size='x-large')
	plt.xlim(ra[0],ra[-1])
	plt.ylim(dec[0],dec[-1])
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)

	plt.show()

	#~ plt.savefig('m57_%s_labeled_wstars_nocbar.pdf'%(cbar_label[0:-13]))#,transparent=True)
	return


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


def define_amplitudes_2comp(lines, waves, spectra_0, dwave):
	'''
	Finds the amplitube (max) of the emission line from the line list provided
	by finding the line in the spectrum and determining the max value of the line.
	Provides an intial guess of the amplitude for the gaussian fitting routine.
	'''
	#~ i_ww = 0
	lineamp = numpy.zeros( numpy.size(lines)*2 )

	for ww in lines:
		i_line = numpy.where((waves <= ww+1) & (waves >= ww-1))[0]
		amp = abs(numpy.max(spectra_0[i_line[0]-dwave:i_line[0]+dwave]))
		#~ if numpy.isnan(amp):
			#~ lineamp[i_ww] = lineamp[i_ww-1]
		#~ else:
		lineamp[0] = amp
		lineamp[1] = amp*0.3		# Make smaller than first component
		#~ i_ww += 1

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
		#~ if (vr[l]<-30) or (vr[l]>50):
			#~ vr[l] = -numpy.inf
		if vdisp[l]>300:
			vdisp[l] = -numpy.inf	#0
	if ((vr[1] - vr[0]) < -75) or (abs(vr[0]-vr[1]) > 75):
		vr = -numpy.inf	#0
	if abs(vdisp[0] - vdisp[1]) > 70:
		vdisp[0] = vdisp[1] = numpy.min(vdisp)	#0
	return vr, vdisp



def velocity_info_single(parms, cen, c):
	'''
	Determine the velocity and dispersion (dynamics) of each line per pixel,
	based on the emission line info extracted from the gaussian fits.
	For single lines, don't have another line to compare behavior to!
	'''
	vr = numpy.zeros( numpy.size(cen) )			# km/s
	vdisp = numpy.zeros( numpy.size(cen) )		# convert to km/s
	for l in range(0, len(parms[:,0])):
		if (parms[l,1] == cen[l]) and (parms[l,2] <= 0):
			vr[l] = -numpy.inf
			vdisp[l] = -numpy.inf
		else:
			vr[l] = ((parms[l,1] - cen[l]) / cen[l])*c
			vdisp[l] = abs(2.3548*parms[l,2])*(c/cen[l])		# in wavelength (Ang) units --> km/s
			if (vr[l]<-500) or (vr[l]>500):
				vr[l] = -numpy.inf
		#~ if vdisp[l]>150:
			#~ vdisp[l] = -numpy.inf	#0

	return vr, vdisp



def velocity_distributions(velocity, regions, bin, legend, color, alp):
	'''
	Use the results from the Gaussian fits over high S/N regions to find
	the distribution of velocities over different regions in the image.
	Find the distribution of velocities per region with their own histograms.
	Region array is 2D, where: [ Region 1:[x1, y1, x2, y2],
								 Region 2:[x1, y1, x2, y2],
								....
								 Region N:[x1, y1, x2, y2] ]
	'''
	# Take out nans and infs
	velocity[numpy.isinf(velocity)] = 5000000
	velocity[numpy.isnan(velocity)] = 5000000

	regs = legend

	for reg in range(0,1): #numpy.size(regs)):	#regions[:,0])):
		if numpy.size(regions) == 4:
			y1 = regions[0]
			y2 = regions[2]
			x1 = regions[1]
			x2 = regions[3]
		else:
			y1 = regions[reg,0]
			y2 = regions[reg,2]
			x1 = regions[reg,1]
			x2 = regions[reg,3]
		# For when multiple lines are read into this function:
		if velocity.ndim > 2:
			vel_hist = numpy.zeros( [(numpy.size(velocity[0,y1:y2,0])*numpy.size(velocity[0,0,x1:x2])),numpy.size(velocity[:,0,0])] )
			for i in range(0,numpy.size(velocity[:,0,0])):
				vel_hist[:,i] = velocity[i,y1:y2,x1:x2].flatten()
			plt.hist(vel_hist, bin, align='mid', histtype='step', fill='True', #density=True, stacked='True',
							 alpha=alp, label=legend, linewidth=5)
			vel = velocity[:,y1:y2,x1:x2].flatten()
		# Just one emission line at once
		else:
			vel = velocity[y1:y2,x1:x2].flatten()
			plt.hist(vel, bin, align='mid', histtype='step', fill='True', #density=True, stacked='True',
							 alpha=alp, color=color, label=legend, linewidth=5)

	if numpy.any(vel < 0):
		plt.axvline(x=v_m57, color='black', linestyle='dashed', lw=3)#, label='M57')
	#~ plt.legend(loc='upper right', fontsize='large')
	#~ ax = plt.gca()
	#~ ax.axes.get_yaxis().set_ticks([])
	#~ plt.savefig('m57_[OI]_Velocity_histogram.pdf',rasterized=True)		# this looks bad as eps; save as PDF
	#~ plt.show()
	return



def velocity_distributions_smoothed(velocity, regions, bin, legend, color, alp):
	'''
	Use the results from the Gaussian fits over high S/N regions to find
	the distribution of velocities over different regions in the image.
	Find the distribution of velocities per region with their own histograms.
	Region array is 2D, where: [ Region 1:[x1, y1, x2, y2],
								 Region 2:[x1, y1, x2, y2],
								....
								 Region N:[x1, y1, x2, y2] ]
	'''
	# Take out nans and infs
	velocity[numpy.isinf(velocity)] = 5000000
	velocity[numpy.isnan(velocity)] = 5000000
	regs = legend	#['Cloud 1','Cloud 2','Cloud 3','Cloud 4']
	xx = numpy.arange(-100,100,1)

	for reg in range(0,1): #numpy.size(regs)):	#regions[:,0])):
		if numpy.size(regions) == 4:
			y1 = regions[0]
			y2 = regions[2]
			x1 = regions[1]
			x2 = regions[3]
		else:
			y1 = regions[reg,0]
			y2 = regions[reg,2]
			x1 = regions[reg,1]
			x2 = regions[reg,3]
		# For when multiple lines are read into this function:
		#~ print velocity.ndim
		if velocity.ndim > 2:
			vel_hist = numpy.zeros( [(numpy.size(velocity[0,y1:y2,0])*numpy.size(velocity[0,0,x1:x2])),numpy.size(velocity[:,0,0])] )
			for i in range(0,numpy.size(velocity[:,0,0])):
				vel_hist[:,i] = velocity[i,y1:y2,x1:x2].flatten()
			density = kde.gaussian_kde(vel_hist)
			plt.plot(xx, density(xx), color=color, label=legend, linewidth=5)
			vel = velocity[:,y1:y2,x1:x2].flatten()
		# Just one emission line at once
		else:
			vel = velocity[y1:y2,x1:x2].flatten()
			density = kde.gaussian_kde(vel)
			plt.plot(xx, density(xx), color=color, label=legend, linewidth=5)

	if numpy.any(vel < 0):
		plt.axvline(x=v_m57, color='black', linestyle='dashed', lw=3)#, label='M57')
	#~ plt.legend(loc='upper right', fontsize='large')
	#~ ax = plt.gca()
	#~ ax.axes.get_yaxis().set_ticks([])
	#~ plt.savefig('m57_[OI]_Velocity_histogram.pdf',rasterized=True)		# this looks bad as eps; save as PDF
	#~ plt.show()
	return



def velocity_distributions_sameregion(velocity, regions, bin, legend, title):
	'''
	*** Same as velocity_distributions, but separating same region by
	    velocity component ***
	Use the results from the Gaussian fits over high S/N regions to find
	the distribution of velocities over different regions in the image.
	Find the distribution of velocities per region with their own histograms.
	Region array is 2D, where: [ Region 1:[x1, y1, x2, y2],
								 Region 2:[x1, y1, x2, y2],
								....
								 Region N:[x1, y1, x2, y2] ]
	'''
	# Take out nans and infs
	velocity[numpy.isinf(velocity)] = -500
	velocity[numpy.isnan(velocity)] = -500

	fig = plt.figure(figsize=(8,8))
	fig.set_rasterized(True)
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_rasterized(True)
	ax1.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
	ax1.set_title(title, weight='bold', size='xx-large')
	#~ ax1.set_yticklabels([])

	#~ clr = ['blue','green','red','yellow']		# [NI]
	#~ clr = ['blue','green','cyan']		# [OI]
	#~ clr = ['purple','gold']		# [OII]
	#~ clr = ['deepskyblue','magenta']		# HeII
	#~ clr = ['deepskyblue','magenta']		# Hbeta
	clr = ['limegreen','orange']		# [OIII]	# [OIII]
	regs = legend	#['Cloud 1','Cloud 2','Cloud 3','Cloud 4']
	alp = [0.6, 0.4, 0.5, 0.6]
	y1 = regions[0]
	y2 = regions[2]
	x1 = regions[1]
	x2 = regions[3]
	# Loop through each box and make a histogram of velocity distributi
	for reg in range(0,numpy.size(velocity[:,0,0])):	#regions[:,0])):
		vel = velocity[reg,x1:x2,y1:y2].flatten()
		plt.hist(vel, bin, align='mid', histtype='step', stacked='True', fill='True', #density=True,
				 color=clr[reg], alpha=alp[reg], label=regs[reg],
				 edgecolor=clr[reg], linewidth=5)

	plt.legend(loc='upper right', fontsize='large')
	ax = plt.gca()
	ax.axes.get_yaxis().set_ticks([])
	#~ plt.savefig('m57_[OI]_Velocity_histogram.pdf',rasterized=True)		# this looks bad as eps; save as PDF
	plt.show()
	return



def ks_test(velocity, region1, region2, bin, legend):
	'''
	Use defined regions to determine the likelihood that two clouds
	share the same velocity distribution, because if they are from the same
	distirbution, they likely
	'''
	# Take out nans and infs
	velocity[numpy.isinf(velocity)] = -100000

	# Define vel1, vel2 regions to make histogram
	# Cloud 1:
	y1 = region1[0]
	y2 = region1[2]
	x1 = region1[1]
	x2 = region1[3]
	#~ vel1 = []
	#~ for i in range(x1,x2):
		#~ for j in range(y1,y2):
			#~ if velocity[:,i,j].any() == numpy.inf:
				#~ print "Infinite value - continue"
				#~ continue
			#~ else:
				#~ vel1.append(velocity[0,i,j],velocity[1,i,j])
	vel1 = velocity[:,x1:x2,y1:y2].flatten()
	indv = numpy.where(vel1 != -100000)[0]
	vel1 = vel1[indv]
	print numpy.size(vel1)
	#~ vel1 = velocity[x1:x2,y1:y2].flatten()
	#~ hist1, egdes1 = numpy.histogram(vel1,bins=50,range=(-100,100))#,density=True)
	# Cloud 2
	y1 = region2[0]
	y2 = region2[2]
	x1 = region2[1]
	x2 = region2[3]
	#~ vel2 = []
	#~ for i in range(x1,x2):
		#~ for j in range(y1,y2):
			#~ if velocity[:,i,j].any() == numpy.inf:
				#~ print "Infinite value - continue"
				#~ continue
			#~ else:
				#~ vel2.append(velocity[:,i,j])
	#~ vel2 = numpy.array(vel2)
	vel2 = velocity[:,x1:x2,y1:y2].flatten()
	indv = numpy.where(vel2 != -100000)[0]
	vel2 = vel2[indv]
	print numpy.size(vel2)
	#~ vel2 = velocity[x1:x2,y1:y2].flatten()
	#~ hist2, edges2 = numpy.histogram(vel2,bins=50,range=(-100,100))#,density=True)


	# Find K-S Test between two samples
	ks_stat, p_ks = ks_2samp(vel1, vel2)
	print "K-S statistic = ", ks_stat,", p-value = ", p_ks," for ", legend[0]," and ", legend[1]
	# Find Student's T-test statistic between the two samples
	tt_stat, p_tt = ttest_ind(vel1, vel2, equal_var = False)
	print "T-Test statistic = ", tt_stat,", p-value = ", p_tt," for ", legend[0]," and ", legend[1]
	print
	vel1[numpy.isinf(vel1)] = -100000
	vel2[numpy.isinf(vel2)] = -100000

	# Plot histograms with density
	fig = plt.figure(figsize=(8,8))
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_xlabel('Velocity (km/s)')
	plt.hist(vel1, bins=20, range=(-100,100), density=True, color='blue', edgecolor='blue',
			 linewidth=5, alpha=0.7, label=legend[0])
	plt.hist(vel2, bins=20, range=(-100,100), density=True, color='green', edgecolor='green',
			 linewidth=5, alpha=0.7, label=legend[1])

	plt.legend(loc='upper right')
	#~ plt.show()


	return





def low_res():
	'''
	Use to look at 04/12 data sets (co-added for max S/N) low-resolution files
	(212, 213)
	'''
	############### Main Ring (inner edge) ##################
	c = 3.0E5
	# Line list
	HI = [4101.73, 4340.47, 4861.35]  	# delta, gamma, beta
	HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
	HeII = [4540, 4685.804, 5411.50]
	NI = [5197.902, 5200.257]
	OII = [3726.032, 3728.815]
	OIII = [4363.209, 4958.911, 5006.843]

	hi_res_ID = ['HeII 4686',r'H$\beta$','[OIII] 4959','[OIII] 5007']
	hi_res_ID = ['HeII','HI','[OIII]','[OIII]']
	hi_res = [4685.804, 4861.350, 4958.911, 5006.843]


	#~ lo_res = [5197.902, 5200.257]
	#~ lo_res_ID= ['[NI]','[NI]']

	date='170412'

	int = 'icubes'
	var = 'vcubes'

	index1=210		# spectral range: 4650 - 5050
	index2=211		#
	index3=212		# spectral range: 3600 - 5500
	index4=213		#

	#~ intfile1 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index1,index2,int)
	#~ file1 = path+dir+date+dir+redux+dir+intfile1
	#~ varfile1 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index1,index2,var)
	#~ vfile1 = path+dir+date+dir+redux+dir+varfile

	intfile2 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,int)
	file2 = path+dir+date+dir+redux+dir+intfile2
	varfile2 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,var)
	vfile2 = path+dir+date+dir+redux+dir+varfile2


	# Read in file, get important data from files
	#~ data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1,1)		# data in units erg cm-2 s-1
	data2, waves2, ra, dec, all_ra, all_dec = open_fits(file2,1)		# data in units erg cm-2 s-1
	waves2 = waves2 - 0.24 #0.18		# offset in lines from BH2
	#~ var1, varwv1 = open_fits_err(vfile1)
	var2, varwv2 = open_fits_err(vfile2)

	# Do things here
	ra_cen = ra - (ra_ref-ra)-0.00394
	dec_cen = dec - (dec_ref-dec)+0.00292
	all_ra = (all_ra - ra_cen)*3600.
	all_dec = (all_dec - dec_cen)*3600.

	#~ data1 = data1[:,32:-31,12:-12]
	#~ var1 = var1[:,32:-31,12:-12]
	data2 = data2[:,32:-30,11:-12]
	var2 = var2[:,32:-30,11:-12]
	all_dec = all_dec[32:-31]
	all_ra = all_ra[12:-12]

	clouds = [0,0,130,22]
	cloud_labels = ['Main Ring']

	# Loop through all lines to make a large velocity/dispersion array for
	# a bunch of lines at once
	# ALL LINES
	#~ lo_res_ID= ['[OII]','[OII]','[NeIII]','[SII]','Hdelta','Hgamma',
							#~ 'HeI','HeII','Hbeta','[OIII]','[OIII]',
							#~ '[NI]','[NI]','[ClIII]','[ClIII]']
	lo_res = [[3726.032, 3728.815], 3868.760, [4068.600, 4076.349], 4101.734, 4340.472,
						4471.48, 4685.804, 4861.350, 4958.911, 5006.843, #[4958.911, 5006.843],	#4713.15,
						[5197.902, 5200.257], [5517.709, 5537.873]]#, 5577.339, [4562.6, 4571.1]]
	lo_res_ID= ['[OII]','[NeIII]','[SII]',r'H$\delta$',r'H$\gamma$',
							'HeI','HeII',r'H$\beta$','[OIII]','[OIII]',
							'[NI]','[ClIII]']#,'[OI]','[MgI']
	# Split into recombination and collision species
	colors = ['red','gold','forestgreen','mediumblue','darkviolet','darkorange','silver']
	# Collisional vs recombination
	cel = [0,1,2,8,9,10]#,11,12]	#  collisonally-excited
	rel = [3,4,5,6,7]	# recombination
	# Same sampling in high-res data set as low-res
	cel = [0,1,2,10,11,3,4]#,11,12]	#  low-res sampling only
	rel = [6,7,5,8,9]	# also sampled in high-res
	# TEST
	#~ lo_res = [3868.760, 4068.600, 4101.734, 4340.472,
						#~ 4713.15, 4685.804, 4861.350, 5006.843, 5517.709]
	#~ lo_res_ID= ['[NeIII]','[SII]','Hdelta','Hgamma',
							#~ 'HeI','HeII','Hbeta','[OIII]','[ClIII]']
	# TEST
	#~ lo_res = [[3726.032, 3728.815], 3868.760, 4068.600, 4101.734, 4861.350]
	#~ lo_res_ID= ['[OII]','[NeIII]','[SII]','Hdelta','Hbeta']
	#~ cel = [0,1,2]
	#~ rel = [3,4]


	lines = lo_res	#hi_res[0]
	velocity_lines = numpy.zeros( [len(lines),numpy.size(data2[0,:,0]),numpy.size(data2[0,0,:])] )
	vdisp_lines = numpy.zeros( [len(lo_res_ID),numpy.size(data2[0,:,0]),numpy.size(data2[0,0,:])] )

	# Loop through the lines
	for line in range(0,len(lines)):

		dlam = 10. 		# +/- extra wavelength coverage from min/max of lines
		wv_range = numpy.where((waves2 > (numpy.min(lines[line])-dlam)) & (waves2 < (numpy.max(lines[line])+dlam)))[0]
		waves_lines = waves2[wv_range]
		data_lines = data2[wv_range,:,:]
		var_lines = var2[wv_range,:,:]

		flux_nocont = numpy.zeros( numpy.shape(data_lines) )

		# 2. Loop through each image pixel to do continuu subtraction
		for i in range(0,numpy.size(data_lines[0,:,0])):
			for j in range(0,numpy.size(data_lines[0,0,:])):
				flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

		# 3. Define S/N ration cut
		sigma = 10		# 2.75 - good with var
		data_cut = sn_cut(flux_nocont, var_lines, sigma)

		#~ # 4. Fit Gaussian profile(s) to line(s) of interest
		for i in range(0,numpy.size(data_lines[0,:,0])):
			for j in range(0,numpy.size(data_lines[0,0,:])):
				if numpy.sum(data_cut[:,i,j]) <= 0:
					print "No lines detected at [%i,%i]" % (i,j)
					#for ww in range(0,len(lines)):
					velocity_lines[line,i,j] = -numpy.inf
					vdisp_lines[line,i,j] = -numpy.inf
					continue
				else:
					# Define amplitude, st. dev, and parameter array for Gaussfit
					dwave = 1
					if type(lines[line]) == list:	#len(lines[line]) == 2:		# for doublet lines
						#~ print "Doublet lines."
						amp_lines = define_amplitudes(lines[line], waves_lines, data_cut[:,i,j], dwave)
						linestdv = numpy.ones( numpy.size(waves_lines) )*0.5 	# Std. dev. array w/ guess:
						gauss_parms = define_parms(lines[line], amp_lines, linestdv)

						# fit gaussian profile(s) to the line(s)
						# get the list of best-fit parameters back
						popt = gaussfit(gauss_parms, waves_lines, data_cut[:,i,j])
						popt = popt.reshape((len(popt)/3,3))

						# Determine velocity (v_r) and dispersion (fwhm) of each line
						vel, disp = velocity_info_doublet(popt, lines[line], c)
						#~ vel, disp = velocity_info_single(popt, [lines[line]], c)
						velocity_lines[line,i,j] = numpy.ma.average(vel)	#,axis=0
						vdisp_lines[line,i,j] = numpy.ma.average(disp)		#,axis=0
					else:													# for single lines
						amp_lines = define_amplitudes([lines[line]], waves_lines, data_cut[:,i,j], dwave)
						linestdv = numpy.ones( numpy.size(waves_lines) )*0.5 	# Std. dev. array w/ guess:
						gauss_parms = define_parms([lines[line]], amp_lines, linestdv)

						# fit gaussian profile(s) to the line(s)
						# get the list of best-fit parameters back
						popt = gaussfit(gauss_parms, waves_lines, data_cut[:,i,j])
						popt = popt.reshape((len(popt)/3,3))

						# Determine velocity (v_r) and dispersion (fwhm) of each line
						#~ vel, disp = velocity_info_doublet(popt, lines, c)
						vel, disp = velocity_info_single(popt, [lines[line]], c)
						velocity_lines[line,i,j] = vel
						vdisp_lines[line,i,j] = disp

	bin = numpy.linspace(-100,100,40)
	bin_disp =  numpy.linspace(0,200,100)
	#~ print velocity_lines[0,:,:]
	#~ print
	#~ print velocity_lines[1,:,:]
	#~ vdisp_lines[vdisp_lines>100] = vdisp_lines[vdisp_lines>100] - 100
	#~ vdisp_lines[vdisp_lines>60] = numpy.median(vdisp_lines)

	#save_outputs(lo_res_ID, velocity_lines, vdisp_lines, cloud_labels, clouds)
	#print "[OII] analysis complete!"

	#~ plot_velocity(velocity_lines[0,:,:], data_cut_contours, levels, lw, all_ra, all_dec, all_ra, all_dec, lo_res_ID, lines, -40, 40, 'bwr', 'Velocity (km/s)')
	#~ plot_velocity(vdisp_lines[0,:,:], data_cut_contours, levels, lw, all_ra, all_dec, all_ra, all_dec, lo_res_ID, lines, 90, 140, 'BuGn', 'Dispersion (km/s)')
	#~ plot_velocity_double(velocity_lines, data_cut_contours, levels, lw, all_ra, all_dec, lo_res_ID, lines, -50, 50, 'bwr', 'Velocity (km/s)')
	#~ plot_velocity_double(vdisp_lines, data_cut_contours, levels, lw, all_ra, all_dec, lo_res_ID, lines, 0, 80, 'BuGn', 'Dispersion (km/s)')

	# Plots and statistics and stuff
	# 1. radial velocity (projected motion)
	# All in one plot (CEL vs RL)
	#~ fig = plt.figure(figsize=(8,8))
	#~ # CEL
	#~ ax1 = fig.add_subplot(1,2,1)
	#~ ax1.set_title('Collisionally-Excited Lines', weight='bold', size='x-large')
	#~ ax1.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
	#~ alp = 1.0
	#~ img = 1
	#~ vel_cel = numpy.zeros( [len(cel),numpy.size(data_cut[0,:,0]),numpy.size(data_cut[0,0,:])] )
	#~ ID_cel = []
	#~ for c in range(0,len(cel)):
		#~ vel_cel[c,:,:] = velocity_lines[cel[c],:,:]
		#~ ID_cel.append( lo_res_ID[cel[c]] )
	#~ velocity_distributions(vel_cel, clouds, bin, ID_cel, 'Collisionally-Excited Lines', alp)
	#~ plt.xlim(-75,29)
	#~ plt.legend(loc='upper right', fontsize='large')
	#~ ax = plt.gca()
	#~ ax.axes.get_yaxis().set_ticks([])

	#~ # REL
	#~ ax2 = fig.add_subplot(1,2,2)
	#~ ax2.set_title('Recombination Lines', weight='bold', size='x-large')
	#~ ax2.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
	#~ vel_rel = numpy.zeros( [len(rel),numpy.size(data_cut[0,:,0]),numpy.size(data_cut[0,0,:])] )
	#~ ID_rel = []
	#~ for r in range(0,len(rel)):
		#~ vel_rel[r,:,:] = velocity_lines[rel[r],:,:]
		#~ ID_rel.append( lo_res_ID[rel[r]] )
	#~ velocity_distributions(vel_rel, clouds, bin, ID_rel, 'Recombination Lines', alp)
	#~ plt.xlim(-75,29)
	#~ plt.legend(loc='upper right', fontsize='large')
	#~ ax = plt.gca()
	#~ ax.axes.get_yaxis().set_ticks([])

	#~ plt.show()


	# 2. radial velocity (projected motion)
	# Split up histrograms into individual histograms, stacked
	fig = plt.figure(figsize=(8,8))
	alp = 0.5
	img = 1
	for c in range(0,len(cel)):
		# CEL
		if c == 0:
			ax1 = fig.add_subplot(len(cel),1,1)
			ax1.set_title('Collisionally-Excited Lines', weight='bold', size='x-large')
			velocity_distributions(velocity_lines[cel[c],:,:], clouds, bin, lo_res_ID[cel[c]], colors[c], alp)
			#~ velocity_distributions_smoothed(velocity_lines[cel[c],:,:], clouds, bin, lo_res_ID[cel[c]], colors[c], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(-75,50)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
		elif c == len(cel)-1:
			ax2 = fig.add_subplot(len(cel),1,img+c, sharex=ax1)
			ax2.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
			velocity_distributions(velocity_lines[cel[c],:,:], clouds, bin, lo_res_ID[cel[c]], colors[c], alp)
			#~ velocity_distributions_smoothed(velocity_lines[cel[c],:,:], clouds, bin, lo_res_ID[cel[c]], colors[c], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(-75,50)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
			fig.subplots_adjust(hspace=0)
			plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
		else:
			ax2 = fig.add_subplot(len(cel),1,img+c, sharex=ax1)
			velocity_distributions(velocity_lines[cel[c],:,:], clouds, bin, lo_res_ID[cel[c]], colors[c], alp)
			#~ velocity_distributions_smoothed(velocity_lines[cel[c],:,:], clouds, bin, lo_res_ID[cel[c]], colors[c], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(-75,50)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])

	#~ plt.show()
	# REL
	bin = numpy.linspace(-100,100,100)
	#~ bin_disp =  numpy.linspace(0,50,50)

	fig = plt.figure(figsize=(8,8))
	alp = 0.5
	img = 1
	for r in range(0,len(rel)):
		# RL
		if r == 0:
			ax1 = fig.add_subplot(len(rel),1,1)
			ax1.set_title('Resonance Lines', weight='bold', size='x-large')
			velocity_distributions(velocity_lines[rel[r],:,:], clouds, bin, lo_res_ID[rel[r]], colors[r], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(-75,50)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
		elif r == len(rel)-1:
			ax2 = fig.add_subplot(len(rel),1,img+r, sharex=ax1)
			ax2.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
			velocity_distributions(velocity_lines[rel[r],:,:], clouds, bin, lo_res_ID[rel[r]], colors[r], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(-75,50)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
			fig.subplots_adjust(hspace=0)
			plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
		else:
			ax2 = fig.add_subplot(len(rel),1,img+r, sharex=ax1)
			velocity_distributions(velocity_lines[rel[r],:,:], clouds, bin, lo_res_ID[rel[r]], colors[r], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(-75,50)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])

	#~ plt.show()

	#~ cel = [0,1,2,8,9,10,11]#,11,12]	#  collisonally-excited
	#~ rel = [3,4,5,6,7]	# recombination

	#~ fig = plt.figure(figsize=(8,8))
	#~ ax1 = fig.add_subplot(2,1,1)
	#~ ax1.set_title('Dispersion', weight='bold', size='x-large')
	#~ ax1.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
	#~ alp = 0.3
	#~ img = 1
	#~ velocity_distributions(vdisp_lines, clouds, bin_disp, lo_res_ID, colors, alp)
	#~ plt.legend(loc='upper left', fontsize='large')
	#~ ax = plt.gca()
	#~ ax.axes.get_yaxis().set_ticks([])



	fig = plt.figure(figsize=(8,8))
	alp = 0.5
	img = 1
	for c in range(0,len(cel)):
		# CEL
		if c == 0:
			ax1 = fig.add_subplot(len(cel),1,1)
			ax1.set_title('Collisionally-Excited Lines', weight='bold', size='x-large')
			velocity_distributions(vdisp_lines[cel[c],:,:], clouds, bin_disp, lo_res_ID[cel[c]], colors[c], alp)
			#~ velocity_distributions_smoothed(velocity_lines[cel[c],:,:], clouds, bin, lo_res_ID[cel[c]], colors[c], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(0,200)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
		elif c == len(cel)-1:
			ax2 = fig.add_subplot(len(cel),1,img+c, sharex=ax1)
			ax2.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
			velocity_distributions(vdisp_lines[cel[c],:,:], clouds, bin_disp, lo_res_ID[cel[c]], colors[c], alp)
			#~ velocity_distributions_smoothed(velocity_lines[cel[c],:,:], clouds, bin, lo_res_ID[cel[c]], colors[c], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(0,200)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
			fig.subplots_adjust(hspace=0)
			plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
		else:
			ax2 = fig.add_subplot(len(cel),1,img+c, sharex=ax1)
			velocity_distributions(vdisp_lines[cel[c],:,:], clouds, bin_disp, lo_res_ID[cel[c]], colors[c], alp)
			#~ velocity_distributions_smoothed(velocity_lines[cel[c],:,:], clouds, bin, lo_res_ID[cel[c]], colors[c], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(0,200)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])


	fig = plt.figure(figsize=(8,8))
	alp = 0.5
	img = 1
	for r in range(0,len(rel)):
		# RL
		if r == 0:
			ax1 = fig.add_subplot(len(rel),1,1)
			ax1.set_title('Resonance Lines', weight='bold', size='x-large')
			velocity_distributions(vdisp_lines[rel[r],:,:], clouds, bin_disp, lo_res_ID[rel[r]], colors[r], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(0,200)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
		elif r == len(rel)-1:
			ax2 = fig.add_subplot(len(rel),1,img+r, sharex=ax1)
			ax2.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
			velocity_distributions(vdisp_lines[rel[r],:,:], clouds, bin_disp, lo_res_ID[rel[r]], colors[r], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(0,200)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
			fig.subplots_adjust(hspace=0)
			plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
		else:
			ax2 = fig.add_subplot(len(rel),1,img+r, sharex=ax1)
			velocity_distributions(vdisp_lines[rel[r],:,:], clouds, bin_disp, lo_res_ID[rel[r]], colors[r], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(0,200)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])

	plt.show()
	#~ ks_test(velocity_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(velocity_lines, clouds[0,:], clouds[2,:], bin, [cloud_labels[0],cloud_labels[2]])
	#~ ks_test(velocity_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])

	#~ ks_test(vdisp_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(vdisp_lines, clouds[0,:], clouds[2,:], bin, [cloud_labels[0],cloud_labels[2]])
	#~ ks_test(vdisp_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])

	return




def hi_res():
	'''
	Use to look at 04/12 data sets (co-added for max S/N) low-resolution files
	(212, 213)
	'''
	############### Main Ring (inner edge) ##################
	c = 3.0E5
	# Line list
	HI = [4101.73, 4340.47, 4861.35]  	# delta, gamma, beta
	HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
	HeII = [4540, 4685.804, 5411.50]
	NI = [5197.902, 5200.257]
	OII = [3726.032, 3728.815]
	OIII = [4363.209, 4958.911, 5006.843]

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


	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1,1)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)

	# Do things here
	ra_cen = ra - (ra_ref-ra)-0.00394
	dec_cen = dec - (dec_ref-dec)+0.00292
	all_ra = (all_ra - ra_cen)*3600.
	all_dec = (all_dec - dec_cen)*3600.

	data1 = data1[:,32:-31,12:-12]
	var1 = var1[:,32:-31,12:-12]
	all_dec = all_dec[32:-31]
	all_ra = all_ra[12:-12]

	clouds = [0,0,130,22]
	cloud_labels = ['Main Ring']

	# Loop through all lines to make a large velocity/dispersion array for
	# a bunch of lines at once
	# ALL LINES
	hi_res_ID = ['HeII 4686',r'H$\beta$','[OIII] 4959','[OIII] 5007', 'HeI 5015']
	hi_res_ID = ['HeII',r'H$\beta$','HeI','[OIII] 4959','[OIII] 5007']
	hi_res = [4685.804, 4861.350, 5015.678, 4958.911, 5006.843]
	# Split into recombination and collision species
	colors = ['red','gold','forestgreen','mediumblue','darkviolet']
	#~ cel = [2]	#  collisonally-excited
	#~ rel = [0,1,3]	# recombination

	lines = hi_res
	velocity_lines = numpy.zeros( [len(lines),numpy.size(data1[0,:,0]),numpy.size(data1[0,0,:])] )
	vdisp_lines = numpy.zeros( [len(hi_res_ID),numpy.size(data1[0,:,0]),numpy.size(data1[0,0,:])] )

	# Loop through the lines
	for line in range(0,len(lines)):
		if line == 2:
			dlam = 1.5
			sigma = 1		# 2.75 - good with var
			#~ wv_range = numpy.where((waves1 > (numpy.min(lines[line])-dlam)) & (waves1 < (numpy.max(lines[line])+3*dlam)))[0]
		else:
			dlam = 5. 		# +/- extra wavelength coverage from min/max of lines
			sigma = 10		# 2.75 - good with var
		wv_range = numpy.where((waves1 > (numpy.min(lines[line])-dlam)) & (waves1 < (numpy.max(lines[line])+dlam)))[0]
		waves_lines = waves1[wv_range]
		data_lines = data1[wv_range,:,:]
		var_lines = var1[wv_range,:,:]

		flux_nocont = numpy.zeros( numpy.shape(data_lines) )

		# 2. Loop through each image pixel to do continuu subtraction
		for i in range(0,numpy.size(data_lines[0,:,0])):
			for j in range(0,numpy.size(data_lines[0,0,:])):
				flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

		# 3. Define S/N ration cut

		data_cut = sn_cut(flux_nocont, var_lines, sigma)

		#~ fig2 = plt.figure(figsize=(8,8))
		#~ ax1 = fig2.add_subplot(1,1,1)
		#~ ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
		#~ ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='x-large')
		#~ ax1.set_title(hi_res_ID[line],weight='bold',size='x-large')
		#~ plt.imshow(numpy.sum(data_cut,axis=0), cmap='nipy_spectral', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
					#~ #vmin = 65, vmax = 150,
					#~ extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
		#~ ax = plt.gca()
		#~ ax.get_xaxis().get_major_formatter().set_useOffset(False)
		#~ ax.get_yaxis().get_major_formatter().set_useOffset(False)
		#~ plt.show()

		#~ # 4. Fit Gaussian profile(s) to line(s) of interest
		for i in range(0,numpy.size(data_lines[0,:,0])):
			for j in range(0,numpy.size(data_lines[0,0,:])):
				if numpy.sum(data_cut[:,i,j]) <= 0:
					print "No lines detected at [%i,%i]" % (i,j)
					#for ww in range(0,len(lines)):
					velocity_lines[line,i,j] = -numpy.inf
					vdisp_lines[line,i,j] = -numpy.inf
					continue
				else:
					# Define amplitude, st. dev, and parameter array for Gaussfit
					dwave = 1
					if type(lines[line]) == list:	#len(lines[line]) == 2:		# for doublet lines
						#~ print "Doublet lines."
						amp_lines = define_amplitudes(lines[line], waves_lines, data_cut[:,i,j], dwave)
						linestdv = numpy.ones( numpy.size(waves_lines) )*0.5 	# Std. dev. array w/ guess:
						gauss_parms = define_parms(lines[line], amp_lines, linestdv)

						# fit gaussian profile(s) to the line(s)
						# get the list of best-fit parameters back
						popt = gaussfit(gauss_parms, waves_lines, data_cut[:,i,j])
						popt = popt.reshape((len(popt)/3,3))

						# Determine velocity (v_r) and dispersion (fwhm) of each line
						vel, disp = velocity_info_doublet(popt, lines[line], c)
						#~ vel, disp = velocity_info_single(popt, [lines[line]], c)
						velocity_lines[line,i,j] = numpy.ma.average(vel)	#,axis=0
						vdisp_lines[line,i,j] = numpy.ma.average(disp)		#,axis=0
					else:													# for single lines
						amp_lines = define_amplitudes([lines[line]], waves_lines, data_cut[:,i,j], dwave)
						linestdv = numpy.ones( numpy.size(waves_lines) )*0.5 	# Std. dev. array w/ guess:
						gauss_parms = define_parms([lines[line]], amp_lines, linestdv)

						# fit gaussian profile(s) to the line(s)
						# get the list of best-fit parameters back
						popt = gaussfit(gauss_parms, waves_lines, data_cut[:,i,j])
						popt = popt.reshape((len(popt)/3,3))

						# Determine velocity (v_r) and dispersion (fwhm) of each line
						#~ vel, disp = velocity_info_doublet(popt, lines, c)
						vel, disp = velocity_info_single(popt, [lines[line]], c)
						velocity_lines[line,i,j] = vel
						vdisp_lines[line,i,j] = disp

	bin = numpy.linspace(-100,100,100)
	bin_disp =  numpy.linspace(0,50,50)

	# 1. radial velocity (projected motion)
	# Split up histrograms into individual histograms, stacked
	fig = plt.figure(figsize=(8,8))
	alp = 0.5
	img = 1
	# ALL LINES
	for c in range(0,len(hi_res_ID)):
		if c == 0:
			ax1 = fig.add_subplot(len(hi_res_ID),1,1)
			ax1.set_title('High-Res', weight='bold', size='x-large')
			velocity_distributions(velocity_lines[c,:,:], clouds, bin, hi_res_ID[c], colors[c], alp)
			#~ velocity_distributions_smoothed(velocity_lines[c,:,:], clouds, bin, hi_res_ID[c], colors[c], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(-75,50)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
		elif c == len(hi_res_ID)-1:
			ax2 = fig.add_subplot(len(hi_res_ID),1,img+c, sharex=ax1)
			ax2.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
			velocity_distributions(velocity_lines[c,:,:], clouds, bin, hi_res_ID[c], colors[c], alp)
			#~ velocity_distributions_smoothed(velocity_lines[c,:,:], clouds, bin, hi_res_ID[c], colors[c], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(-75,50)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])
			fig.subplots_adjust(hspace=0)
			plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
		else:
			ax2 = fig.add_subplot(len(hi_res_ID),1,img+c, sharex=ax1)
			velocity_distributions(velocity_lines[c,:,:], clouds, bin, hi_res_ID[c], colors[c], alp)
			#~ velocity_distributions_smoothed(velocity_lines[c,:,:], clouds, bin, hi_res_ID[c], colors[c], alp)
			plt.legend(loc='upper right', fontsize='large')
			plt.xlim(-75,50)
			ax = plt.gca()
			ax.axes.get_yaxis().set_ticks([])

	#~ plt.show()

	fig = plt.figure(figsize=(8,8))
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_title('Dispersion', weight='bold', size='x-large')
	ax1.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
	alp = 0.3
	img = 1
	velocity_distributions(vdisp_lines, clouds, bin_disp, hi_res_ID, colors, alp)
	#~ velocity_distributions_sameregion(velocity_lines, clouds, bin, lo_res_ID, '[NI] 5197+5200: Cloud Velocities')
	#~ velocity_distributions_sameregion(vdisp_lines, clouds, bin_disp, lo_res_ID, '[NI] 5197+5200: Cloud Dispersion')
	plt.legend(loc='upper right', fontsize='large')
	ax = plt.gca()
	ax.axes.get_yaxis().set_ticks([])
	plt.show()

	return





if __name__ == "__main__":
	global path, dir, redux
	global ra_ref, dec_ref		# Central star coordinates - to subtract all all_ra, all_dec from

	global ra_stars, dec_stars, v_m57
	path='C:\\Users\\Keri Hoadley\\Documents\\KCWI'
	dir = '\\'
	#path='/home/keri/KCWI/'
	#dir = '/'
	redux='redux'

	# Line list
	HI = [4101.73, 4340.47, 4861.35]  	# delta, gamma, beta
	HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
	HeII = [4540, 4686, 5411]
	NI = [5198.5, 5199]
	OI = [5577,6300]
	OII = [3726, 3729, 4267, 4661]
	OIII = [4363.21, 4958.91, 5006.84, 5592.37]


	# RA/DEC of central star
	v_m57 = -19.0		# km/s, moving towards observer
	ra_cen = 283.319
	dec_cen = 33.01775
	ra_stars = [0,-7]		# RA of central stars, in deltaRA
	dec_stars = [0, 4]		# DEC of central stars, in deltaDEC

	ra_ref = 360.*( (18. + (53/60.) + (35.079/3600.)) / 24.)			# 18:53:35.079
	dec_ref = 33. + (01/60.) + (45.03/3600.)		#+33:01:45.03



	# Functions to run high and low resolution data to find v_r and v_disp
	low_res()
	hi_res()
