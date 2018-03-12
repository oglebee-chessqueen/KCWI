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


def plot_velocity_double(velocity, contours, lvls, lwds, ra, dec, ID, line, min, max, color, cbar_label):
	'''
	Plot doublet velocity distribution(s), including contour regions.
	'''
	velocity[numpy.isinf(velocity)] = 0
	# Define contour levels
	levels = lvls
	lw = lwds
	fig1 = plt.figure(figsize=(15,8),facecolor='white')
	# Plot 1
	ax1 = fig1.add_subplot(1,2,1)
	ax1.set_facecolor('white')
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)',weight='bold',size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\circ}$)',weight='bold',size='x-large')
	ax1.set_title(ID[0]+'%i'%(line[0]),weight='bold',size='x-large')
	CS = plt.contour(contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[ra[0],ra[-1],dec[0],dec[-1]])
	im = plt.imshow(velocity[0,:,:], cmap=color, origin='lower',
			   vmin = min, vmax = max,
			   extent=[ra[0],ra[-1],dec[0],dec[-1]])
	ax = plt.gca()
	plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
	plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
	# Plot 2
	ax2 = fig1.add_subplot(1,2,2)
	ax2.set_facecolor('white')
	ax2.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)',weight='bold',size='x-large')
	ax2.set_title(ID[1]+'%i'%(line[1]),weight='bold',size='x-large')
	CS = plt.contour(contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[ra[0],ra[-1],dec[0],dec[-1]])
	im = plt.imshow(velocity[1,:,:], cmap=color, origin='lower',
			   vmin = min, vmax = max,
			   extent=[ra[0],ra[-1],dec[0],dec[-1]])
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = fig1.colorbar(im, cax=cax)
	cbar.ax.set_ylabel(cbar_label, weight='bold', size='x-large')

	plt.savefig('m57_%s_%s.eps'%(ID[0],cbar_label[0:-7]))#,transparent=True)
	return



def plot_velocity(velocity, contours, lvls, lwds, ra, dec, ID, line, min, max, color, cbar_label):
	'''
	Plot doublet velocity distribution(s), including contour regions.
	Velocity - should only be 2D, of one line or avg. velocity of multiple lines
	'''
	velocity[numpy.isinf(velocity)] = 0

	levels = lvls
	lw = lwds
	fig = plt.figure(figsize=(8,8),facecolor='white')
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)', weight='bold', size='large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\circ}$)', weight='bold', size='large')
	CS = plt.contour(contours, levels,
					 linewidths=lw, colors='k', corner_mask=True,
					 extent=[ra[0],ra[-1],dec[0],dec[-1]])
	im = plt.imshow(velocity, cmap=color, origin='lower',	#cmap='gnuplot'
					vmin = min, vmax = max,
					extent=[ra[0],ra[-1],dec[0],dec[-1]])
	#~ cbar = plt.colorbar()		# Each image should have its own colorbar
	ax = plt.gca()
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = fig.colorbar(im, cax=cax)
	cbar.ax.set_ylabel(cbar_label, weight='bold', size='large')
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)

	#~ plt.show()

	plt.savefig('m57_%s.eps'%(cbar_label[0:-7]))#,transparent=True)
	return



def plot_intensity_map(data, contours, lvls, lwds, ra, dec, ID, line, min, max, color, cbar_label):
	'''
	Plot doublet velocity distribution(s), including contour regions.
	'''
	levels = lvls
	lw = lwds
	fig = plt.figure(figsize=(8,8))
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)', weight='bold', size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\circ}$)', weight='bold', size='x-large')
	CS = plt.contour(contours, levels,
					 linewidths=lw, colors='w', corner_mask=True,
					 extent=[ra[0],ra[-1],dec[0],dec[-1]])
	im = plt.imshow(numpy.sum(data,axis=0), cmap=color, origin='lower',	#cmap='gnuplot'
				#vmin = 0, vmax = 100,
				extent=[ra[0],ra[-1],dec[0],dec[-1]])
	#~ cbar = plt.colorbar()		# Each image should have its own colorbar
	ax = plt.gca()
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = fig.colorbar(im, cax=cax)
	cbar.ax.set_ylabel(cbar_label, weight='bold', size='x-large')
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)

	#~ plt.show()

	plt.savefig('m57_%s.eps'%(cbar_label[0:-13]))#,transparent=True)
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





def velocity_distributions(velocity, regions, bin, legend):
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
	velocity[numpy.isinf(velocity)] = -500

	fig = plt.figure(figsize=(8,8))
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_xlabel('Velocity (km/s)', weight='bold', size='large')

	clr = ['blue','green','red','yellow','magenta','cyan']
	regs = legend	#['Cloud 1','Cloud 2','Cloud 3','Cloud 4']
	alp = [0.3, 0.5, 0.7, 0.9]
	# Loop through each box and make a histogram of velocity distributi
	for reg in range(0,numpy.size(regions[:,0])):
		y1 = regions[reg,0]
		y2 = regions[reg,2]
		x1 = regions[reg,1]
		x2 = regions[reg,3]
		vel = velocity[:,x1:x2,y1:y2].flatten()
		plt.hist(vel, bin, align='mid', #density=True,
				 color=clr[reg], alpha=alp[reg], label=regs[reg],
				 edgecolor=clr[reg], linewidth=5)

	plt.legend(loc='upper right')
	ax = plt.gca()
	ax.axes.get_yaxis().set_ticks([])
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
	plt.hist(vel1, bins=50, range=(-100,100), density=True, color='blue', edgecolor='blue',
			 linewidth=5, alpha=0.7, label=legend[0])
	plt.hist(vel2, bins=50, range=(-100,100), density=True, color='green', edgecolor='green',
			 linewidth=5, alpha=0.7, label=legend[1])

	plt.legend(loc='upper right')
	#~ plt.show()


	return




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


	# Define regions around the contours, and plot the rectangles to make sure they
	# encapsulate the clouds of interest
	levels = [25]
	lw = [2.5]

	clouds = [ [6,35,23,58],[11,23,18,35],[20,0,23,10],[0,59,5,69] ]
	clouds = numpy.array(clouds)
	cloud_labels = ['Cloud 1','Cloud 2','Cloud 3','Cloud 4']

	#~ fig2 = plt.figure(figsize=(8,8))
	#~ ax2 = fig2.add_subplot(1,1,1)
	#~ ax2.set_xlabel(r'$\Delta \alpha$ ($^{\circ}$)',weight='bold',size='large')
	#~ ax2.set_ylabel(r'$\Delta \delta$ ($^{\circ}$)',weight='bold',size='large')
	#~ CS = plt.contour(data_cut_contours, levels,
					 #~ linewidths=lw, colors='k', corner_mask=True)#,
					 #~ #extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
	#~ plt.imshow(numpy.sum(data_cut,axis=0), cmap='Blues', origin='lower')#,	#cmap='gnuplot'
				#~ #vmin = 0, vmax = 100,
				#extent=[all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]])
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
			#(all_ra[clouds[cl,0]]-2.5, all_dec[clouds[cl,1]]),   # (x,y)
			#numpy.abs( all_ra[clouds[cl,2]] - all_ra[clouds[cl,0]] ),          # width
			#numpy.abs( all_dec[clouds[cl,3]] - all_dec[clouds[cl,1]] ),          # height
			#fill=False, edgecolor="red"      # remove background
		#),
		#]:
			#~ ax2.add_patch(p)
	#cbar = plt.colorbar()		# Each image should have its own colorbar
	#cbar.ax.set_ylabel(r'[NI] Intensity (arb. units)',
					    #weight='bold',size='large')
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
				vel, disp = velocity_info_doublet(popt, lines, c)
				velocity_lines[:,i,j] = vel
				vdisp_lines[:,i,j] = disp

	#print all_ra[0],all_ra[-1],all_dec[0],all_dec[-1]
	#~ velocity_lines = gaussian_filter(velocity_lines, 0.1)
	#~ vdisp_lines = gaussian_filter(vdisp_lines, 0.1)
	bin = numpy.linspace(-100,100,30)
	velocity_avg = numpy.average(velocity_lines,axis=0)
	#~ velocity_distributions(velocity_lines, clouds, bin, cloud_labels)
	#~ plot_velocity(velocity_avg, data_cut_contours, [25], [2.5], all_ra, all_dec, ID, lines, -100, 100, 'bwr', '[NI] Avg. Velocity (km/s)')
	#~ ks_test(velocity_lines, clouds[0,:], clouds[3,:], bin, [cloud_labels[0],cloud_labels[3]])
	#~ ks_test(velocity_lines, clouds[1,:], clouds[2,:], bin, [cloud_labels[1],cloud_labels[2]])
	#~ ks_test(velocity_lines, clouds[0,:], clouds[1,:], bin, [cloud_labels[0],cloud_labels[1]])
	#~ ks_test(velocity_lines, clouds[2,:], clouds[3,:], bin, [cloud_labels[2],cloud_labels[3]])
	plt.show()

	# Plot velocity maps:
	plot_velocity_double(velocity_lines, data_cut_contours, [25], [2.5], all_ra, all_dec, ID, lines, -100, 100, 'bwr', 'Velocity (km/s)')
	plot_velocity_double(vdisp_lines, data_cut_contours, [25], [2.5], all_ra, all_dec, ID, lines, 0, 20, 'BuGn', 'Dispersion (km/s)')

	# Plot S/N cut data:
	plot_intensity_map(data_cut, data_cut_contours, [25], [2.5], all_ra, all_dec, ID, lines, 0, 100, 'gnuplot', '[NI] Intensity (arb. units)')




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
