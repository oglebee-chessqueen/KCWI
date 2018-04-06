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
from scipy.stats import ks_2samp, ttest_ind
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
		vr[l] = ((parms[l,1] - cen[l]) / cen[l])*c
		vdisp[l] = abs(2.3548*parms[l,2])*(c/cen[l])		# in wavelength (Ang) units --> km/s
		if (vr[l]<-500) or (vr[l]>500):
			vr[l] = -numpy.inf
		#~ if vdisp[l]>150:
			#~ vdisp[l] = -numpy.inf	#0

	return vr, vdisp



def velocity_distributions(velocity, regions, bin, legend, title):
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

	fig = plt.figure(figsize=(8,8))
	fig.set_rasterized(True)
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_rasterized(True)
	ax1.set_xlabel('Velocity (km/s)', weight='bold', size='x-large')
	ax1.set_title(title, weight='bold', size='xx-large')
	#~ ax1.set_yticklabels([])

	#~ clr = ['blue','green','red','yellow']		# [NI]
	clr = ['blue','green','cyan']		# [OI]
	#~ clr = ['teal']#'purple']		# [OII]
	#~ clr = ['magenta']		# HeII
	#~ clr = ['deepskyblue','magenta']		# Hbeta
	#~ clr = ['limegreen','orange']		# [OIII]
	regs = legend	#['Cloud 1','Cloud 2','Cloud 3','Cloud 4']
	alp = [0.3, 0.4, 0.5, 0.6]
	#~ print numpy.size(regions)
	# Loop through each box and make a histogram of velocity distributi
	for reg in range(0,numpy.size(regs)):	#regions[:,0])):
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
		vel = velocity[:,x1:x2,y1:y2].flatten()
		plt.hist(vel, bin, align='mid', histtype='step', stacked='True', fill='True', density=True,
				 color=clr[reg], alpha=alp[reg], label=regs[reg],
				 edgecolor=clr[reg], linewidth=5)

	if numpy.any(vel < 0):
		plt.axvline(x=v_m57, color='black', linestyle='dashed', lw=3, label='M57')
	plt.legend(loc='upper right', fontsize='large')
	ax = plt.gca()
	ax.axes.get_yaxis().set_ticks([])
	#~ plt.savefig('m57_[OI]_Velocity_histogram.pdf',rasterized=True)		# this looks bad as eps; save as PDF
	plt.show()
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




def main_ring_contour_plots():
	'''
	Use to look at 04/12 data sets (co-added for max S/N)
	'''
	############### Main Ring (inner edge) ##################
	# Line list
	HI = [4101.73, 4340.47, 4861.35]  	# delta, gamma, beta
	HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
	HeII = [4540, 4685.804, 5411.50]
	NI = [5197.902, 5200.257]
	OII = [3726.032, 3728.815]
	OIII = [4363.209, 4958.911, 5006.843]

	hi_res_ID = ['HeII 4686',r'H$\beta$','[OIII] 4959','[OIII] 5007']
	hi_res = [4685.804, 4861.350, 4958.911, 5006.843]

	lo_res_ID = ['[OII] 3726','[OII] 3729']
	lo_res = [3726.032, 3728.815]

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
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1,1)		# data in units erg cm-2 s-1
	data2, waves2, ra2, dec2, all_ra2, all_dec2 = open_fits(file2,1)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)
	var2, varwv2 = open_fits_err(vfile2)


	# Do things here
	#~ all_ra = (all_ra - ra)*3600.
	#~ all_dec = (all_dec - dec)*3600.
	ra_cen = ra - (ra_ref-ra)-0.00394
	dec_cen = dec - (dec_ref-dec)+0.00292
	all_ra = (all_ra - ra_cen)*3600.
	all_dec = (all_dec - dec_cen)*3600.

	data1 = data1[:,32:-31,12:-12]
	var1 = var1[:,32:-31,12:-12]
	data2 = data2[:,32:-30,11:-12]
	var2 = var2[:,32:-30,11:-12]
	all_dec = all_dec[32:-31]
	all_ra = all_ra[12:-12]

	#~ waves1 = waves2
	#~ data1 = data2
	#~ var1 = var2


	# Do things here
	# 1. Define emission line(s) to define velocity maps
	lines = lo_res	#hi_res[0]
	dlam = 5. 		# +/- extra wavelength coverage from min/max of lines
	wv_range = numpy.where((waves2 > (numpy.min(lines)-dlam)) & (waves2 < (numpy.max(lines)+dlam)))[0]
	#~ wv_range = numpy.where((waves1 > 5060) & (waves1 < 5300))[0]
	waves_lines = waves2[wv_range]
	data_lines = data2[wv_range,:,:]
	var_lines = var2[wv_range,:,:]

	flux_nocont = numpy.zeros( numpy.shape(data_lines) )

	# 2. Loop through each image pixel to do continuu subtraction
	for i in range(0,numpy.size(data_lines[0,:,0])):
		for j in range(0,numpy.size(data_lines[0,0,:])):
			flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

	# 3. Define S/N ration cut
	sigma = 10	#2.5			# 2.75 - good with var
	#~ sigma = 2.5#2.2#33
	data_cut_lo = sn_cut(flux_nocont, var_lines, sigma)
	data_cut_contours_lo = gaussian_filter(numpy.sum(data_cut_lo,axis=0), 1.0)


	levels = [30,80,93,110,130]
	lw = [1,1,2,3,4,5]
	#~ levels = [15,30,68,80,95,112,130,140]
	#~ lw = [1,1,1,1.5,2,2.5,3,4]


	fig2 = plt.figure(figsize=(16,8))
	ax1 = fig2.add_subplot(1,4,1)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax1.set_title(lo_res_ID[0]+'+3729',weight='bold',size='x-large')
	CS2 = plt.contour(data_cut_contours_lo, levels,
					 linewidths=lw, colors='white', #corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(numpy.sum(data_cut_lo,axis=0), cmap='nipy_spectral', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
				#~ vmin = 30, vmax = 190,
				extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)



	# Do things here
	# 1. Define emission line(s) to define velocity maps
	lines = hi_res[1]	#[hi_res[2],hi_res[3]]
	dlam = 10	#5. 		# +/- extra wavelength coverage from min/max of lines
	wv_range = numpy.where((waves1 > (numpy.min(lines)-dlam)) & (waves1 < (numpy.max(lines)+dlam)))[0]
	waves_lines = waves1[wv_range]
	data_lines = data1[wv_range,:,:]
	var_lines = var1[wv_range,:,:]

	flux_nocont = numpy.zeros( numpy.shape(data_lines) )

	# 2. Loop through each image pixel to do continuu subtraction
	for i in range(0,numpy.size(data_lines[0,:,0])):
		for j in range(0,numpy.size(data_lines[0,0,:])):
			flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

	# 3. Define S/N ration cut
	sigma = 3	#5# 2.5			# 2.75 - good with var
	#~ sigma = 2.5#2.2#33
	data_cut = sn_cut(flux_nocont, var_lines, sigma)
	#~ data_cut_contours = gaussian_filter(numpy.sum(data_cut,axis=0), 0.3)
	data_cut_contours = gaussian_filter(numpy.sum(data_cut,axis=0), 1.1)

	levels = [75,100,110,120,140]
	lw = [1,1,2,3,4]

	ax2 = fig2.add_subplot(1,4,2)
	ax2.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	#~ ax2.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax2.set_title(hi_res_ID[1],weight='bold',size='x-large')
	CS2 = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='black', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(numpy.sum(data_cut,axis=0), cmap='nipy_spectral', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
				vmin = 65, vmax = 150,
				extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)


	# Do things here
	# 1. Define emission line(s) to define velocity maps
	lines = [hi_res[2],hi_res[3]]
	dlam = 10	#5. 		# +/- extra wavelength coverage from min/max of lines
	wv_range = numpy.where((waves1 > (numpy.min(lines)-dlam)) & (waves1 < (numpy.max(lines)+dlam)))[0]
	waves_lines = waves1[wv_range]
	data_lines = data1[wv_range,:,:]
	var_lines = var1[wv_range,:,:]

	flux_nocont = numpy.zeros( numpy.shape(data_lines) )

	# 2. Loop through each image pixel to do continuu subtraction
	for i in range(0,numpy.size(data_lines[0,:,0])):
		for j in range(0,numpy.size(data_lines[0,0,:])):
			flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

	# 3. Define S/N ration cut
	sigma = 10	#5			# 2.75 - good with var
	#~ sigma = 2.5#2.2#33
	data_cut = sn_cut(flux_nocont, var_lines, sigma)

	# Define regions around the contours, and plot the rectangles to make sure they
	# encapsulate the clouds of interest
	#~ levels = [20,30,35]
	#~ lw = [2.5,3.5,5]

	ax3 = fig2.add_subplot(1,4,3)
	ax3.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	#~ ax2.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax3.set_title(hi_res_ID[2]+'+5007',weight='bold',size='x-large')
	levels = [75,100,110,120,140]
	lw = [1,1,2,3,4,5]
	CS1 = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	levels = [30,80,93,110,130]
	#~ lw = [1,1,2,2.5,3,4]
	CS2 = plt.contour(data_cut_contours_lo, levels,
					 linewidths=lw, colors='white', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(numpy.sum(data_cut,axis=0), cmap='gnuplot', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
				vmin = 1200,
				extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])

	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)


	# Do things here
	# 1. Define emission line(s) to define velocity maps
	lines = hi_res[0]
	dlam = 5. 		# +/- extra wavelength coverage from min/max of lines
	wv_range = numpy.where((waves1 > (numpy.min(lines)-dlam)) & (waves1 < (numpy.max(lines)+dlam)))[0]
	waves_lines = waves1[wv_range]
	data_lines = data1[wv_range,:,:]
	var_lines = var1[wv_range,:,:]

	flux_nocont = numpy.zeros( numpy.shape(data_lines) )

	# 2. Loop through each image pixel to do continuu subtraction
	for i in range(0,numpy.size(data_lines[0,:,0])):
		for j in range(0,numpy.size(data_lines[0,0,:])):
			flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

	# 3. Define S/N ration cut
	sigma = 5			# 2.75 - good with var
	#~ sigma = 2.5#2.2#33
	data_cut = sn_cut(flux_nocont, var_lines, sigma)

	# Define regions around the contours, and plot the rectangles to make sure they
	# encapsulate the clouds of interest
	#~ levels = [20,30,35]
	#~ lw = [2.5,3.5,5]

	ax3 = fig2.add_subplot(1,4,4)
	ax3.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	#~ ax2.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax3.set_title(hi_res_ID[0],weight='bold',size='x-large')
	levels = [75,100,110,120,140]
	lw = [1,1,2,3,4,5]
	CS1 = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	levels = [30,80,93,110,130]
	#~ lw = [1,1,2,2.5,3,4]
	CS2 = plt.contour(data_cut_contours_lo, levels,
					 linewidths=lw, colors='white', corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(numpy.sum(data_cut,axis=0), cmap='gnuplot', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
				vmin = 3, vmax = 40,
				extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])

	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	plt.tight_layout()
	plt.savefig('m57_Ring_emission_contours.pdf')#,transparent=True)
	plt.show()

	return


def main_ring():
	'''
	Use to look at 04/12 data sets (co-added for max S/N)
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

	lo_res_ID = ['[OII]','[OII]']
	lo_res = [3726.032, 3728.815]

	#~ lo_res = [5197.902, 5200.257]
	#~ lo_res_ID= ['[NI]','[NI]']

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

	#~ intfile2 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,int)
	intfile2 = 'kb'+date+'_00%03i_%s.fits' % (index3,int)
	file2 = path+dir+date+dir+redux+dir+intfile2
	#~ varfile2 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,var)
	varfile2 = 'kb'+date+'_00%03i_%s.fits' % (index3,var)
	vfile2 = path+dir+date+dir+redux+dir+varfile2


	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1,1)		# data in units erg cm-2 s-1
	data2, waves2, ra2, dec2, all_ra2, all_dec2 = open_fits(file2,1)		# data in units erg cm-2 s-1
	waves2 = waves2 - 0.25		# offset in lines from BH2
	var1, varwv1 = open_fits_err(vfile1)
	var2, varwv2 = open_fits_err(vfile2)

	# Do things here
	ra_cen = ra - (ra_ref-ra)-0.00394
	dec_cen = dec - (dec_ref-dec)+0.00292
	all_ra = (all_ra - ra_cen)*3600.
	all_dec = (all_dec - dec_cen)*3600.

	data1 = data1[:,32:-31,12:-12]
	var1 = var1[:,32:-31,12:-12]
	data2 = data2[:,32:-30,11:-12]
	var2 = var2[:,32:-30,11:-12]
	all_dec = all_dec[32:-31]
	all_ra = all_ra[12:-12]

	clouds = [0,0,130,23]
	cloud_labels = ['Main Ring']

	# Do things here
	# A) find velocity behavior of lower-res data sets ([OII])
	# 1. Define emission line(s) to define velocity maps
	lines = lo_res	#hi_res[0]
	dlam = 5. 		# +/- extra wavelength coverage from min/max of lines
	wv_range = numpy.where((waves2 > (numpy.min(lines)-dlam)) & (waves2 < (numpy.max(lines)+dlam)))[0]
	#~ wv_range = numpy.where((waves1 > 5060) & (waves1 < 5300))[0]
	waves_lines = waves2[wv_range]
	data_lines = data2[wv_range,:,:]
	var_lines = var2[wv_range,:,:]

	flux_nocont = numpy.zeros( numpy.shape(data_lines) )

	# 2. Loop through each image pixel to do continuu subtraction
	for i in range(0,numpy.size(data_lines[0,:,0])):
		for j in range(0,numpy.size(data_lines[0,0,:])):
			flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

	# 3. Define S/N ration cut
	sigma = 10			# 2.75 - good with var
	#~ sigma = 2.5#2.2#33
	data_cut = sn_cut(flux_nocont, var_lines, sigma)
	data_cut_contours = gaussian_filter(numpy.sum(data_cut,axis=0), 1.0)

	levels = [30,80,93,110,130]
	lw = [1,2,3,4,5]
	fig2 = plt.figure(figsize=(8,8))
	ax1 = fig2.add_subplot(1,2,1)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax1.set_title(lo_res_ID[0],weight='bold',size='x-large')
	CS2 = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='white', #corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(numpy.sum(data_cut,axis=0), cmap='nipy_spectral', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
				#vmin = 30, vmax = 190,
				extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)


	intfile2 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,int)
	file2 = path+dir+date+dir+redux+dir+intfile2
	varfile2 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index3,index4,var)
	vfile2 = path+dir+date+dir+redux+dir+varfile2


	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1,1)		# data in units erg cm-2 s-1
	data2, waves2, ra2, dec2, all_ra2, all_dec2 = open_fits(file2,1)		# data in units erg cm-2 s-1
	waves2 = waves2 - 0.18		# offset in lines from BH2
	var1, varwv1 = open_fits_err(vfile1)
	var2, varwv2 = open_fits_err(vfile2)

	# Do things here
	ra_cen = ra - (ra_ref-ra)-0.00394
	dec_cen = dec - (dec_ref-dec)+0.00292
	all_ra = (all_ra - ra_cen)*3600.
	all_dec = (all_dec - dec_cen)*3600.

	data1 = data1[:,32:-31,12:-12]
	var1 = var1[:,32:-31,12:-12]
	data2 = data2[:,32:-30,11:-12]
	var2 = var2[:,32:-30,11:-12]
	all_dec = all_dec[32:-31]
	all_ra = all_ra[12:-12]

	lo_res = [5197.902, 5200.257]
	lo_res_ID= ['[NI]','[NI]']
	#~ lo_res = [4068.600]#, 4076.349]
	#~ lo_res_ID= ['[SII]']#,'[SII]']
	#~ lo_res = [5517.709, 5537.873]
	#~ lo_res_ID= ['[ClIII]','[ClIII]']
	#~ lo_res = [3868.760]
	#~ lo_res_ID= ['[NeIII]']
	#~ lo_res = [4340.472]
	#~ lo_res_ID= ['Hgam']
	#~ lo_res = [4861.350]
	#~ lo_res_ID= ['Hbeta']
	#~ lo_res = [4101.734]
	#~ lo_res_ID= ['Hdelta']
	#~ lo_res_ID = ['[OII]','[OII]']
	#~ lo_res = [3726.032, 3728.815]


	lines = lo_res	#hi_res[0]
	dlam = 5. 		# +/- extra wavelength coverage from min/max of lines
	wv_range = numpy.where((waves2 > (numpy.min(lines)-dlam)) & (waves2 < (numpy.max(lines)+dlam)))[0]
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



	ax1 = fig2.add_subplot(1,2,2)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax1.set_title(lo_res_ID[0],weight='bold',size='x-large')
	CS2 = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='white', #corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(numpy.sum(data_cut,axis=0), cmap='nipy_spectral', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
				#vmin = 30, vmax = 190,
				extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	plt.show()

	#~ # 4. Fit Gaussian profile(s) to line(s) of interest
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
				#~ vel, disp = velocity_info_doublet(popt, lines, c)
				vel, disp = velocity_info_single(popt, lines, c)
				velocity_lines[:,i,j] = vel
				vdisp_lines[:,i,j] = disp

	bin = numpy.linspace(-50,50,40)
	bin = numpy.linspace(-100,100,40)
	bin_disp =  numpy.linspace(0,200,80)
	#~ vdisp_lines[vdisp_lines>100] = vdisp_lines[vdisp_lines>100] - 100
	#~ vdisp_lines[vdisp_lines>60] = numpy.median(vdisp_lines)

	#save_outputs(lo_res_ID, velocity_lines, vdisp_lines, cloud_labels, clouds)
	#print "[OII] analysis complete!"

	#~ plot_velocity(velocity_lines[0,:,:], data_cut_contours, levels, lw, all_ra, all_dec, all_ra, all_dec, lo_res_ID, lines, -40, 40, 'bwr', 'Velocity (km/s)')
	#~ plot_velocity(vdisp_lines[0,:,:], data_cut_contours, levels, lw, all_ra, all_dec, all_ra, all_dec, lo_res_ID, lines, 90, 140, 'BuGn', 'Dispersion (km/s)')
	#~ plot_velocity_double(velocity_lines, data_cut_contours, levels, lw, all_ra, all_dec, lo_res_ID, lines, -50, 50, 'bwr', 'Velocity (km/s)')
	#~ plot_velocity_double(vdisp_lines, data_cut_contours, levels, lw, all_ra, all_dec, lo_res_ID, lines, 0, 80, 'BuGn', 'Dispersion (km/s)')

	# Plots and statistics and stuff
	velocity_distributions(velocity_lines, clouds, bin, cloud_labels, '%s: Velocities'%(lo_res_ID[0]))
	velocity_distributions(vdisp_lines, clouds, bin_disp, cloud_labels, '%s: Dispersion'%(lo_res_ID[0]))
	#~ velocity_distributions_sameregion(velocity_lines, clouds, bin, lo_res_ID, '[NI] 5197+5200: Cloud Velocities')
	#~ velocity_distributions_sameregion(vdisp_lines, clouds, bin_disp, lo_res_ID, '[NI] 5197+5200: Cloud Dispersion')

	#~ ks_test(velocity_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(velocity_lines, clouds[0,:], clouds[2,:], bin, [cloud_labels[0],cloud_labels[2]])
	#~ ks_test(velocity_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])

	#~ ks_test(vdisp_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(vdisp_lines, clouds[0,:], clouds[2,:], bin, [cloud_labels[0],cloud_labels[2]])
	#~ ks_test(vdisp_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])


	# B) find velocity behavior of high-res data sets (Hbeta, HeII, [OIII])
	# 1. Define emission line(s) to define velocity maps
	#~ lines = [hi_res[1]]	# Hbeta
	#~ ID = [hi_res_ID[1]]	#[hi_res_ID[1]+' Comp1', hi_res_ID[1]+' Comp2']	#
	lines = [hi_res[0]]	# HeII
	ID = [hi_res_ID[0]]	#[hi_res_ID[0]+' Comp1', hi_res_ID[0]+' Comp2']
	#~ lines = [hi_res[2],hi_res[3]]	# [OIII]
	#~ ID = [hi_res_ID[2], hi_res_ID[3]]	#[hi_res_ID[2]+' Comp1', hi_res_ID[2]+' Comp2', hi_res_ID[3]+' Comp1', hi_res_ID[3]+' Comp2']
	dlam = 10	#5. 		# +/- extra wavelength coverage from min/max of lines
	wv_range = numpy.where((waves1 > (numpy.min(lines)-dlam)) & (waves1 < (numpy.max(lines)+dlam)))[0]
	waves_lines = waves1[wv_range]
	data_lines = data1[wv_range,:,:]
	var_lines = var1[wv_range,:,:]

	flux_nocont = numpy.zeros( numpy.shape(data_lines) )

	# 2. Loop through each image pixel to do continuu subtraction
	for i in range(0,numpy.size(data_lines[0,:,0])):
		for j in range(0,numpy.size(data_lines[0,0,:])):
			flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

	# 3. Define S/N ration cut
	sigma =  10	#5	#5# 2.5			# HeII, OIII
	#~ sigma = 3	#2.2#33				# Hbeta
	data_cut = sn_cut(flux_nocont, var_lines, sigma)
	#~ data_cut_contours = gaussian_filter(numpy.sum(data_cut,axis=0), 0.3)
	data_cut_contours = gaussian_filter(numpy.sum(data_cut,axis=0), 1.1)

	levels = [5,15,25,29,32]	# He II
	#~ levels = [75,100,110,120,135]	# Hbeta
	#~ levels = [1200,1450,1700,2100,2400]	# [OIII]
	lw = [1,1,2,3,4]
	#~ fig2 = plt.figure(figsize=(8,8))
	#~ ax1 = fig2.add_subplot(1,1,1)
	#~ ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	#~ ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	#~ ax1.set_title(hi_res_ID[2],weight='bold',size='x-large')
	#~ CS2 = plt.contour(data_cut_contours, levels,
					 #~ linewidths=lw, colors='white', #corner_mask=True,
					 #~ extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ plt.imshow(numpy.sum(data_cut,axis=0), cmap='nipy_spectral', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
				#~ #vmin = 65, vmax = 150,
				#~ extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ ax = plt.gca()
	#~ ax.get_xaxis().get_major_formatter().set_useOffset(False)
	#~ ax.get_yaxis().get_major_formatter().set_useOffset(False)
	#~ plt.show()

	# 4. Fit Gaussian profile(s) to line(s) of interest

	#~ # *** SET UP FOR 2-component FIT!!! ***
	#~ ID = [hi_res_ID[1]+' Comp1', hi_res_ID[1]+' Comp2']
	#~ velocity_lines = numpy.zeros( [len(lines)*2,numpy.size(data_cut[0,:,0]),numpy.size(data_cut[0,0,:])] )
	#~ vdisp_lines = numpy.zeros( [len(lines)*2,numpy.size(data_cut[0,:,0]),numpy.size(data_cut[0,0,:])] )
	#~ for i in range(0,numpy.size(data_lines[0,:,0])):
		#~ for j in range(0,numpy.size(data_lines[0,0,:])):
			#~ if numpy.sum(data_cut[:,i,j]) <= 0:
				#~ print "No lines detected at [%i,%i]" % (i,j)
				#~ for ww in range(0,len(lines)):
					#~ velocity_lines[ww,i,j] = -numpy.inf
					#~ vdisp_lines[ww,i,j] = -numpy.inf
				#~ continue
			#~ else:
				#~ # Define amplitude, st. dev, and parameter array for Gaussfit
				#~ dwave = 1
				#~ amp_lines = define_amplitudes_2comp(lines, waves_lines, data_cut[:,i,j], dwave)
				#~ linestdv = numpy.ones( numpy.size(velocity_lines[:,0,0]) )*0.5 	# Std. dev. array w/ guess:
				#~ for stv in range(0,numpy.size(linestdv)):
					#~ if (linestdv[stv]%2) == 0:
						#~ linestdv[stv] = 2.5
				#~ gauss_parms = define_parms(numpy.array([lines,lines]), amp_lines, linestdv)

				#~ # fit gaussian profile(s) to the line(s)
				#~ # get the list of best-fit parameters back
				#~ popt = gaussfit(gauss_parms, waves_lines, data_cut[:,i,j])
				#~ popt = popt.reshape((len(popt)/3,3))

				#~ # Determine velocity (v_r) and dispersion (fwhm) of each line
				#~ vel, disp = velocity_info_doublet(popt, numpy.array([lines,lines]), c)
				#~ velocity_lines[:,i,j] = vel
				#~ vdisp_lines[:,i,j] = disp

	#~ ID = [hi_res_ID[1]]	#[hi_res_ID[1]+' Comp1', hi_res_ID[1]+' Comp2']
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
				linestdv = numpy.ones( numpy.size(waves_lines) )*2.5 	# Std. dev. array w/ guess:
				gauss_parms = define_parms(lines, amp_lines, linestdv)

				# fit gaussian profile(s) to the line(s)
				# get the list of best-fit parameters back
				popt = gaussfit(gauss_parms, waves_lines, data_cut[:,i,j])
				popt = popt.reshape((len(popt)/3,3))

				# Determine velocity (v_r) and dispersion (fwhm) of each line
				vel, disp = velocity_info_single(popt, lines, c)
				velocity_lines[:,i,j] = vel
				vdisp_lines[:,i,j] = disp

	bin = numpy.linspace(-50,-20,40)				# Hbeta, HeII
	bin_disp =  numpy.linspace(20,60,40)		# Hbeta
	bin_disp =  numpy.linspace(0,40,40)			# HeII

	save_outputs(ID, velocity_lines, vdisp_lines, cloud_labels, clouds)
	#print "%s analysis complete!" % (ID[0])
	#~ bin_disp =  numpy.linspace(10,50,40)			# OIII
	#~ vdisp_lines[vdisp_lines>100] = vdisp_lines[vdisp_lines>100] - 100
	#~ vdisp_lines[vdisp_lines>60] = numpy.median(vdisp_lines)
	#~ velocity_avg = numpy.average(velocity_lines,axis=0)
	#~ plot_velocity(velocity_lines[0,:,:], data_cut_contours, levels, lw, all_ra, all_dec, all_ra, all_dec, ID, lines, -40, -25, 'Blues_r', 'Velocity (km/s)','HI_Velocity_v2')
	#~ plot_velocity(vdisp_lines[0,:,:], data_cut_contours,  levels, lw, all_ra, all_dec, all_ra, all_dec, ID, lines, 35, 50, 'BuGn', 'Dispersion (km/s)','HI_Dispersion_v2')
	#~ plot_velocity_double(velocity_lines, data_cut_contours, levels, lw, all_ra, all_dec, ID, lines, -35, -20, 'Blues_r', 'Velocity (km/s)')
	#~ plot_velocity_double(vdisp_lines, data_cut_contours, levels, lw, all_ra, all_dec, ID, lines, 20, 40, 'BuGn', 'Dispersion (km/s)')

	# Plots and statistics and stuff
	#~ clouds = [0,0,130,23]
	#~ cloud_labels = ['Main Ring']
	#~ velocity_distributions(velocity_lines, clouds, bin, cloud_labels, 'HeII: Velocities')
	#~ velocity_distributions(vdisp_lines, clouds, bin_disp, cloud_labels, 'HeII: Dispersion')
	#~ velocity_distributions_sameregion(velocity_lines, clouds, bin, ID, '[OIII]: Velocities')
	#~ velocity_distributions_sameregion(vdisp_lines, clouds, bin_disp, ID, '[OIII]: Dispersion')
	#~ ks_test(velocity_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(velocity_lines, clouds[0,:], clouds[2,:], bin, [cloud_labels[0],cloud_labels[2]])
	#~ ks_test(velocity_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])

	#~ ks_test(vdisp_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(vdisp_lines, clouds[0,:], clouds[2,:], bin, [cloud_labels[0],cloud_labels[2]])
	#~ ks_test(vdisp_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])




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
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1,1)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)

	#~ all_ra = (all_ra - ra)*3600.
	#~ all_dec = (all_dec - dec)*3600.
	ra_cen = ra - (ra_ref-ra)
	dec_cen = dec - (dec_ref-dec)+0.00065#*10
	all_ra = (all_ra - ra_cen)*3600.
	all_dec = (all_dec - dec_cen)*3600.


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


	# Define regions around the contours, and plot the rectangles to make sure they
	# encapsulate the clouds of interest
	#~ levels = [25]
	#~ lw = [2.5]

	clouds = [ [6,35,23,58],[11,23,18,35],[20,0,23,10],[0,59,5,69] ]
	clouds = numpy.array(clouds)
	cloud_labels = ['Cloud 1','Cloud 2','Cloud 3','Cloud 4']

	#~ fig2 = plt.figure(figsize=(8,8))
	#~ ax2 = fig2.add_subplot(1,1,1)
	#~ ax2.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='large')
	#~ ax2.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='large')
	#~ CS = plt.contour(data_cut_contours, levels,
					 #~ linewidths=lw, colors='k', corner_mask=True,
					 #~ extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ plt.imshow(numpy.sum(data_cut,axis=0), cmap='Blues', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
				#~ #vmin = 0, vmax = 100,
				#~ extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ # Overplot region locations with cloud labels
	#~ plt.text(all_ra[clouds[0,2]]+7, all_dec[clouds[0,3]], cloud_labels[0],weight='bold',size='x-large')	# Cloud 1
	#~ plt.text(all_ra[clouds[1,0]]+1.5, all_dec[clouds[1,3]], cloud_labels[1],weight='bold',size='x-large')	# Cloud 2
	#~ plt.text(all_ra[clouds[2,0]]+1.0, all_dec[clouds[2,3]]+0.5, cloud_labels[2],weight='bold',size='x-large')	# Cloud 3
	#~ plt.text(all_ra[clouds[3,2]]+0.5, all_dec[clouds[3,3]]-0.5, cloud_labels[3],weight='bold',size='x-large')	# Cloud 4

	#~ for cl in range(0,numpy.size(clouds[:,0])):
		#~ for p in [
		#~ patches.Rectangle(
			#~ (clouds[cl,0], clouds[cl,1]),   # (x,y)
			#~ numpy.abs( clouds[cl,2] - clouds[cl,0] ),          # width
			#~ numpy.abs( clouds[cl,3] - clouds[cl,1] ),          # height
			#~ fill=False, edgecolor="red"      # remove background
		#~ ),
		#~ ]:
		#patches.Rectangle(
			#(all_ra[clouds[cl,0]], all_dec[clouds[cl,1]]),   # (x,y)
			#numpy.abs( all_ra[clouds[cl,2]] - all_ra[clouds[cl,0]] ),          # width
			#numpy.abs( all_dec[clouds[cl,3]] - all_dec[clouds[cl,1]] ),          # height
			#fill=False, edgecolor="red"      # remove background
		#),
		#]:
			#~ ax2.add_patch(p)
	#~ #cbar = plt.colorbar()		# Each image should have its own colorbar
	#~ #cbar.ax.set_ylabel(r'[NI] Intensity (arb. units)',
					    #~ #weight='bold',size='large')
	#~ ax = plt.gca()
	#~ ax.get_xaxis().get_major_formatter().set_useOffset(False)
	#~ ax.get_yaxis().get_major_formatter().set_useOffset(False)
	#~ plt.show()


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
	#~ #velocity_lines = gaussian_filter(velocity_lines, 0.1)
	#~ #vdisp_lines = gaussian_filter(vdisp_lines, 0.1)
	bin = numpy.linspace(-100,100,30)
	bin_disp = numpy.linspace(0,80,20)
	#~ velocity_avg = numpy.average(velocity_lines,axis=0)

	# Plots and statistics and stuff
	#~ velocity_distributions(velocity_lines, clouds, bin, cloud_labels, '[NI] 5197+5200: Cloud Velocities')
	#~ velocity_distributions(vdisp_lines, clouds, bin_disp, cloud_labels, '[NI] 5197+5200: Cloud Dispersion')

	#~ plot_velocity(velocity_avg, data_cut_contours, [25], [2.5], all_ra, all_dec, ID, lines, -100, 100, 'bwr', '[NI] Avg. Velocity (km/s)')
	#~ ks_test(velocity_lines, clouds[0,:], clouds[3,:], bin, [cloud_labels[0],cloud_labels[3]])
	#~ ks_test(velocity_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])
	#~ ks_test(velocity_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(velocity_lines, clouds[2,:], clouds[3,:], bin, [cloud_labels[2],cloud_labels[3]])
	#~ plt.show()

	#~ ks_test(vdisp_lines, clouds[0,:], clouds[3,:], bin, [cloud_labels[0],cloud_labels[3]])
	#~ ks_test(vdisp_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])
	#~ ks_test(vdisp_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(vdisp_lines, clouds[2,:], clouds[3,:], bin, [cloud_labels[2],cloud_labels[3]])
	#~ plt.show()

	# Plot velocity maps:
	#~ plot_velocity_double(velocity_lines, data_cut_contours, [25], [2.5], all_ra, all_dec, ID, lines, -100, 100, 'bwr', 'Velocity (km/s)')
	#~ plot_velocity_double(vdisp_lines, data_cut_contours, [25], [2.5], all_ra, all_dec, ID, lines, 0, 20, 'BuGn', 'Dispersion (km/s)')

	# Plot S/N cut data:
	#~ plot_intensity_map(data_cut, data_cut_contours, [25], [2.5], all_ra, all_dec, all_ra, all_dec, ID, lines, 0, 100, 'gnuplot', '[NI] Intensity (arb. units)')
	#~ plot_intensity_map_labeled(data_cut, data_cut_contours, [25], [2.5], all_ra, all_dec, all_ra, all_dec, ID, lines, 0, 100, 'gnuplot', '[NI] Intensity (arb. units)', clouds, cloud_labels)

	#~ save_outputs(ID, velocity_lines, vdisp_lines, cloud_labels, clouds)

	print "[NI] analysis complete!"



	return data_cut_contours, all_ra, all_dec




def central_star_offset_medium(data_cut_contours, NIra, NIdec):
	'''
	Use to look at 04/15 data set (with [OI], maybe [NII] around stars,
	and [OI], [NII], [OIII], HeII, and maybe other lines? in Inner Ring)
	'''
	############### Main Ring (inner edge) ##################
	# Line list
	HeII = [5411]
	NII = [5755]
	OI = [5577.339, 6300.00]
	OIII = [5592.37]
	ClIII = [5517, 5537]

	c = 3.0e5		# km/s


	############### Central Region 1: 06/20 ##################
	date='170415'

	int = 'icuber'
	var = 'vcuber'

	index1 = 188		# not ready

	intfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,int)
	file1 = path+dir+date+dir+redux+dir+intfile1
	varfile1 = 'kb'+date+'_00%03i_%s.fits' % (index1,var)
	vfile1 = path+dir+date+dir+redux+dir+varfile1


	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)

	# Do things here
	ra_cen = ra - ((ra_ref-ra)/-1.5) + 0.0002
	dec_cen = dec - (dec_ref-dec)-0.00225
	all_ra = (all_ra - ra_cen)*3600.
	all_dec = (all_dec - dec_cen)*3600.
	data1 = data1[:,1:-2,:]
	var1 = var1[:,1:-2,:]
	all_dec = all_dec[1:-2]


	# Do things here
	# 1. Define emission line(s) to define velocity maps
	lines = [OI[1]]
	ID = ['[OI]']#,'[OI]']
	dlam = 10. 		# +/- extra wavelength coverage from min/max of lines
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
	sigma = 3.5		# 2.75 - good with var
	#~ sigma = 2.5#2.2#33
	data_cut = sn_cut(flux_nocont, var_lines, sigma)
	print numpy.average(data_cut[data_cut>10])
	data_cut[data_cut>100] = 0
	#~ data_cut_contours = scipy.ndimage.zoom(numpy.sum(data_cut,axis=0), 3)
	#~ data_cut_contours = gaussian_filter(numpy.sum(data_cut,axis=0), 0.9)


	# Define regions around the contours, and plot the rectangles to make sure they
	# encapsulate the clouds of interest
	levels = [25]
	lw = [2.5]

	clouds = [ [0,10,15,38],[1,0,8,10],[20,45,23,66] ]
	clouds = numpy.array(clouds)
	cloud_labels = ['Cloud 1','Cloud 2','Inner Ring']

	#~ plot_intensity_map_labeled(data_cut, data_cut_contours, [25], [2.5], all_ra, all_dec, NIra, NIdec, ID, lines, 0, 100, 'gnuplot', '[OI] Intensity (arb. units)', clouds, cloud_labels)

	#~ fig2 = plt.figure(figsize=(8,8))
	#~ ax2 = fig2.add_subplot(1,1,1)
	#~ ax2.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='large')
	#~ ax2.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='large')
	#~ CS = plt.contour(data_cut_contours, levels,
					 #~ linewidths=lw, colors='k', corner_mask=True,
					 #~ extent=[NIra[0],NIra[-1],NIdec[0],NIdec[-1]])
					 #~ #extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ plt.imshow(numpy.sum(data_cut,axis=0), cmap='Blues', origin='lower', aspect='auto',	#cmap='gnuplot'
				#~ vmin = 0, vmax = 75,
				#~ extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ # Overplot region locations with cloud labels
	#~ plt.text(all_ra[clouds[0,2]]+8, all_dec[clouds[0,3]]-1.0, cloud_labels[0],weight='bold',size='x-large')	# Cloud 1
	#~ plt.text(all_ra[clouds[1,0]], all_dec[clouds[1,3]]+0.5, cloud_labels[1],weight='bold',size='x-large')	# Cloud 2
	#~ plt.text(all_ra[clouds[2,0]]+2.0, all_dec[clouds[2,3]]-2.0, cloud_labels[2],weight='bold',size='x-large')	# Cloud 3

	#~ for cl in range(0,numpy.size(clouds[:,0])):
		#~ for p in [
		#~ patches.Rectangle(
			#~ (clouds[cl,0], clouds[cl,1]),   # (x,y)
			#~ numpy.abs( clouds[cl,2] - clouds[cl,0] ),          # width
			#~ numpy.abs( clouds[cl,3] - clouds[cl,1] ),          # height
			#~ fill=False, edgecolor="red"      # remove background
		#~ ),
		#~ ]:
		#~ patches.Rectangle(
			#~ (all_ra[clouds[cl,0]], all_dec[clouds[cl,1]]),   # (x,y)
			#~ numpy.abs( all_ra[clouds[cl,2]] - all_ra[clouds[cl,0]] ),          # width
			#~ numpy.abs( all_dec[clouds[cl,3]] - all_dec[clouds[cl,1]] ),          # height
			#~ fill=False, edgecolor="red"      # remove background
		#~ ),
		#~ ]:
			#~ ax2.add_patch(p)
	#~ #cbar = plt.colorbar()		# Each image should have its own colorbar
	#~ #cbar.ax.set_ylabel(r'[NI] Intensity (arb. units)',
					    #~ #weight='bold',size='large')
	#~ ax = plt.gca()
	#~ ax.get_xaxis().get_major_formatter().set_useOffset(False)
	#~ ax.get_yaxis().get_major_formatter().set_useOffset(False)
	#~ plt.tight_layout()
	#~ plt.show()

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
				vel, disp = velocity_info_single(popt, lines, c)
				velocity_lines[:,i,j] = vel
				vdisp_lines[:,i,j] = disp

	bin = numpy.linspace(-100,100,25)
	bin_disp = numpy.linspace(0,80,20)
	#~ velocity_avg = numpy.average(velocity_lines,axis=0)
	plot_velocity(velocity_lines[0,:,:], data_cut_contours, [25], [2.5], all_ra, all_dec, NIra, NIdec, ID, lines, -100, 100, 'bwr', '[OI]6300 Velocity (km/s)')
	plot_velocity(vdisp_lines[0,:,:], data_cut_contours, [25], [2.5], all_ra, all_dec, NIra, NIdec, ID, lines,  0, 20, 'BuGn', '[OI]6300 Dispersion (km/s)')

	# Plots and statistics and stuff
	#~ velocity_distributions(velocity_lines, clouds, bin, cloud_labels, '[OI] 6300: Cloud Velocities')
	#~ velocity_distributions(vdisp_lines, clouds, bin_disp, cloud_labels, '[OI] 6300: Cloud Dispersion')

	#~ ks_test(velocity_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(velocity_lines, clouds[0,:], clouds[2,:], bin, [cloud_labels[0],cloud_labels[2]])
	#~ ks_test(velocity_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])

	#~ ks_test(vdisp_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(vdisp_lines, clouds[0,:], clouds[2,:], bin, [cloud_labels[0],cloud_labels[2]])
	#~ ks_test(vdisp_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])
	#~ plt.show()

	#~ save_outputs(ID, velocity_lines, vdisp_lines, cloud_labels, clouds)
	print "[OI] analysis complete!"



	return



def central_star_offset_large():
	'''
	Use to look at 06/19 data set of central stars (no clouds) and Inner Ring
	'''
	############### Main Ring (inner edge) ##################
	# Line list
	# Line list
	linelist = [4471.480, 4541.6, 4634.14, 4640.64, 4649.135, 4661.633, 4562.6, 4571.1]	#4650.839
	lineIDs = ['HeI', 'HeII', 'NIII', 'NIII', 'OII', 'OII', '[MgI]','[MgI]']	#, 'OII'
	c = 3.0e5

	############### Central Region 1: 06/20 ##################
	date='170619'

	int = 'icube'
	var = 'vcube'

	index1 = 260		# not ready
	index2 = 261

	intfile1 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index1,index2,int)
	file1 = path+dir+date+dir+redux+dir+intfile1
	varfile1 = 'kb'+date+'_00%03i+00%03i_%s_coadd.fits' % (index1,index2,var)
	vfile1 = path+dir+date+dir+redux+dir+varfile1


	# Read in file, get important data from files
	data1, waves1, ra, dec, all_ra, all_dec = open_fits(file1)		# data in units erg cm-2 s-1
	var1, varwv1 = open_fits_err(vfile1)

	# Do things here
	ra_cen = ra - (ra_ref-ra)+0.10325
	dec_cen = dec - (dec_ref-dec)-0.00399
	all_ra = (all_ra - ra_cen)*3600.
	all_dec = (all_dec - dec_cen)*3600.
	data1 = data1[:,2:-2,:]
	var1 = var1[:,2:-2,:]
	all_dec = all_dec[2:-2]

	clouds = [0,0,66,23]
	cloud_labels = ['Central Region + Inner Ring']



	# 1. Define emission line(s) to define velocity maps
	#~ lines = [linelist[0]]	# HeI
	#~ ID = [lineIDs[0]]	#[hi_res_ID[2]+' Comp1', hi_res_ID[2]+' Comp2', hi_res_ID[3]+' Comp1', hi_res_ID[3]+' Comp2']
	#~ lines = [linelist[1]]	# HeII
	#~ ID = [lineIDs[1]]	#[hi_res_ID[2]+' Comp1', hi_res_ID[2]+' Comp2', hi_res_ID[3]+' Comp1', hi_res_ID[3]+' Comp2']
	#~ lines = [linelist[2],linelist[3]]	# NIII
	#~ ID = [lineIDs[2],lineIDs[3]]	#[hi_res_ID[2]+' Comp1', hi_res_ID[2]+' Comp2', hi_res_ID[3]+' Comp1', hi_res_ID[3]+' Comp2']
	#~ lines = [linelist[4],linelist[5]]	# OII
	#~ ID = [lineIDs[4],lineIDs[5]]	#[hi_res_ID[2]+' Comp1', hi_res_ID[2]+' Comp2', hi_res_ID[3]+' Comp1', hi_res_ID[3]+' Comp2']
	lines = [linelist[6],linelist[7]]	# OII
	ID = [lineIDs[6],lineIDs[7]]	#[hi_res_ID[2]+' Comp1', hi_res_ID[2]+' Comp2', hi_res_ID[3]+' Comp1', hi_res_ID[3]+' Comp2']


	dlam = 5	#5. 		# +/- extra wavelength coverage from min/max of lines
	wv_range = numpy.where((waves1 > (numpy.min(lines)-dlam)) & (waves1 < (numpy.max(lines)+dlam)))[0]
	waves_lines = waves1[wv_range]
	data_lines = data1[wv_range,:,:]
	var_lines = var1[wv_range,:,:]

	flux_nocont = numpy.zeros( numpy.shape(data_lines) )

	# 2. Loop through each image pixel to do continuu subtraction
	for i in range(0,numpy.size(data_lines[0,:,0])):
		for j in range(0,numpy.size(data_lines[0,0,:])):
			flux_nocont[:,i,j] = subtract_continuum(waves_lines,data_lines[:,i,j])

	# 3. Define S/N ration cut
	sigma = 2.5	#5	#5# 2.5			# HeII, OIII
	data_cut = sn_cut(flux_nocont, var_lines, sigma)
	data_cut[data_cut>90] = 0
	#~ data_cut = sn_cut(data_lines, var_lines, sigma)
	#~ data_cut_contours = gaussian_filter(numpy.sum(data_cut,axis=0), 0.3)
	data_cut_contours = gaussian_filter(numpy.sum(data_cut,axis=0), 1.5)

	#~ levels = [100, 300,600,1000,1500]#,1750]	# He I
	#~ levels = [310]#,500]#,120,135]	# He II
	#~ levels = [390,500,700]#,500]#,120,135]	# NIII
	#~ levels = [30,100,200]	# OII
	levels = [90,200,370]	# MgI
	#~ lw = [1,1,2,3,4]
	#~ lw = [3,4]
	lw = [1,3,4]
	fig2 = plt.figure(figsize=(8,8))
	ax1 = fig2.add_subplot(1,1,1)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='x-large')
	ax1.set_title(ID[0],weight='bold',size='x-large')
	CS2 = plt.contour(data_cut_contours, levels,
					 linewidths=lw, colors='white', #corner_mask=True,
					 extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.imshow(numpy.sum(data_cut,axis=0), cmap='nipy_spectral', origin='lower', #aspect='auto')#,	#cmap='gnuplot'
				#~ vmin = 1000, vmax = 2000,
				extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	plt.scatter(ra_stars,dec_stars, s=1000, facecolor='black', marker='*', #alpha=0.5,
							edgecolor='white', linewidths=2.5)
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	plt.show()

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
				#~ vel, disp = velocity_info_single(popt, lines, c)
				vel, disp = velocity_info_doublet(popt, lines, c)
				velocity_lines[:,i,j] = vel
				vdisp_lines[:,i,j] = disp

	bin = numpy.linspace(-100,100,20)
	bin_disp =  numpy.linspace(0,120,20)
	# For HeI, HeII lines
	#~ plot_velocity(velocity_lines[0,:,:], data_cut_contours, levels, lw, all_ra, all_dec, all_ra, all_dec, ID, lines, -30, 0, 'Blues_r', 'HeI Velocity (km/s)')
	#~ plot_velocity(vdisp_lines[0,:,:], data_cut_contours, levels, lw, all_ra, all_dec, all_ra, all_dec, ID, lines,  0, 150, 'BuGn', 'HeI Dispersion (km/s)')
	# For NIII, OII lines
	#~ plot_velocity_double(velocity_lines, data_cut_contours, levels, lw, all_ra, all_dec, ID, lines, -40, 40, 'bwr', 'Velocity (km/s)')
	#~ plot_velocity_double(vdisp_lines, data_cut_contours, levels, lw, all_ra, all_dec, ID, lines, 0, 150, 'BuGn', 'Dispersion (km/s)')

	velocity_distributions(velocity_lines, clouds, bin, cloud_labels, '%s: Cloud Velocities'%(ID[0]))
	velocity_distributions(vdisp_lines, clouds, bin_disp, cloud_labels, '%s: Cloud Dispersion'%(ID[0]))

	save_outputs(ID, velocity_lines, vdisp_lines, cloud_labels, clouds)

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
	npzfiles = ['[NI]_kcwi_velocities_regions.npz',
							'[OI]_kcwi_velocities_regions.npz',
							'[OII]_kcwi_velocities_regions.npz',
							'[OIII]_kcwi_velocities_regions.npz',
							'HI_kcwi_velocities_regions.npz',
							'HeII_kcwi_velocities_regions.npz']

	# RA/DEC of central star
	v_m57 = -19.0		# km/s, moving towards observer
	ra_cen = 283.319
	dec_cen = 33.01775
	ra_stars = [0,-7]		# RA of central stars, in deltaRA
	dec_stars = [0, 4]		# DEC of central stars, in deltaDEC

	ra_ref = 360.*( (18. + (53/60.) + (35.079/3600.)) / 24.)			# 18:53:35.079
	dec_ref = 33. + (01/60.) + (45.03/3600.)		#+33:01:45.03



	# Run each function per region probed
	''' THESE SHOULD BE DONE! '''
	#~ NI_contours, NIra, NIdec = central_star()		#path,dir,redux)
	#~ central_star_offset_medium(NI_contours, NIra, NIdec)
	main_ring()
	#~ main_ring_contour_plots()		# Makes 4-panel plot of Main Ring emission,
															# including contours to highlight regions of
															# interest ([OII] and Hbeta)
	#~ central_star_offset_large()
