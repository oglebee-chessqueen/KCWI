def onedgaussian(wave,flux,err,dwave):
	'''
	READ IN SUBSECTION SURROUNDING SPECTRAL FEATURE!
	Make a 2D Gaussian fit to model the spectral feature.
	Return the peak amplitude (A) and its uncertainty in the covariance matrix.
	Also return a plot of the best fit Gaussian in a 2D image.
	Finally, return the chi2 value and the PTE of the fit.

	Parameters
	----------
	data (wave,flux,err,dwave) - ndarray
	wave - wavelength/velocity arrarry (x-axis)
	flux - flux array
	err - error in flux array
	dwave - offset from center (x0)
	The imag eof the star that we want to fit with our 2D Gaussian.

	Parameters
	----------
	N - int
	Number of degrees of freedom in our array
	N = size(data) - len(theta)
	theta - array
	The array of parameters we want to explore for our Gauss fit
	covar - ndarray
	Covarience matrix determined from the Gaussian best fit model
	chi2 - float
	chi2 value
	pte - float
	Our probability to exceed (PTE), telling us how well our
	Gaussian fit fits with our data.
	'''
	# import modules
	from scipy.optimize import leastsq
	import scipy.stats
	import numpy as num
	import pylab
	import math

	# Initialize theta parameters
	# theta = [A, B, x0, sigma]  # A = 1.0; normalized
	theta_init = [1.0, flux.min(), dwave, 50.0]

	x = wave
	y = flux

	# Determine uncertainties from the data.
	sigk = err

	# Define our Gaussian function, using the parameters v[0] -> v[4]
	# and including the array sigy
	g = lambda theta_init, x: theta_init[0] + theta_init[1]*num.exp( (-1)*( (x-theta_init[2])**2 )/( 2*theta_init[3]**2 ) )

	# Make the residuals for the least squares fit
	residual = lambda theta_init, x, y, sigk: ( (g(theta_init,x) - y)**2 ) / sigk**2


	# Call leastsq function
	theta, covar, infodict, mesg, ier = leastsq(residual,theta_init,args=(x,y,sigk),full_output=1)
	print 'A:', theta[0]
	print 'B:', theta[1]
	print 'x0', theta[2]
	print 'sigma', abs(theta[3])
	print
	print 'Covarience matrix:'
	print covar
#  print 'Marginal uncertainty in B:', num.sqrt(covar[1,1])
	print
	print 'FWHM = ',abs(theta[3])*2.3548
	print



	# Find chi2 of the gaussian with the theta parameters
	gaussian = []
	area_under_gaussian = 0.0

	for i in range(0,len(y)):
		gauss_fnct = theta[0] + theta[1]*math.exp( (-1)*( (x[i]-theta[2])**2 )/( 2*theta[3]**2 ))
		area_under_gaussian += (1.0 - gauss_fnct)
		print x[-1]-x[0], gauss_fnct, area_under_gaussian*(x[-1]-x[0])
		gaussian.append( gauss_fnct )


	chisq = scipy.stats.chisquare(gaussian)

	print 'Chi**2: ', chisq

	# Find the CDF of chi2, and determine the PTE based on this CDF value
	# Define the DOF of the data set, where DOF = # of data points - # of parameters
	dof = len(y) - len(theta)
	chicdf = scipy.stats.chi2.cdf(chisq,dof)
	pte = 1. - chicdf
	print "DOF:", dof
	print 'PTE: ', pte

	# Plot the estimated Gaussian absorption line over
	# the observed absorption line
	gaussian = num.array(gaussian)

	#~ pylab.errorbar(x,y,yerr=sigk,errorevery=5,color='k',ecolor='k',elinewidth=3)
	#~ pylab.plot(x,gaussian,'b--',lw=3)
	#~ pylab.plot((-500,500),(0,0),'r--',linewidth=5)	# plot 1, 0 normalized lines
	#~ pylab.plot((-500,500),(1,1),'r--',linewidth=5)
	#~ pylab.plot((equiv_width+dwave,equiv_width+dwave),(-0.5,1.5),'r--',linewidth=5)	# Plot +ve equiv width
	#~ pylab.plot((equiv_width-dwave,equiv_width-dwave),(-0.5,1.5),'r--',linewidth=5)	# Plot -ve equiv width
	#~ pylab.xlabel('Velocity (km/s)')
	#~ pylab.ylabel('Normalized flux')
	#~ pylab.xlim(-40,80)
	#~ pylab.ylim(-0.1,1.2)
	#~ pylab.title('Best-fit Gaussian absorption line over the data')
	#~ pylab.grid('on')
	#~ pylab.show()


	return theta, gaussian, area_under_gaussian




def open_txt(filename,plot_title):

	# Import modules
	import astropy.io.fits as fits
	import numpy as num
	import pylab
	from matplotlib.colors import LogNorm

	# Open the text file, and read the lines
	# File is set up with columns as: wavelength (Ang), flux(ergs/cm2/s/Ang), err(ers/cm2/s/Ang)
	flines = open(filename).readlines()

	# Read through each line, and append colymn values to correct arrays
	wave = []
	flux = []
	err = []

	for line in flines:
		column = line.split()
		wave.append( float(column[0]) )
		flux.append( float(column[1]) )
		err.append( float(column[2]) )

	# Print img quick to see star pattern
	pylab.errorbar(wave,flux,yerr=err,errorevery=100,color='k',ecolor='k',elinewidth=3)
	pylab.xlabel('Wavelength (Ang)')
	pylab.ylabel('Flux (ergs/cm^2/s/Ang)')
	pylab.xlim(1170,1350)
	#~ pylab.ylim(0,1.4E-10)
	pylab.title(plot_title)
	pylab.grid('on')
	pylab.show()

	return wave, flux, err



if __name__=='__main__':
	# Define variables, including reference path and file name
	import astropy
	import math
	import pylab
	from matplotlib.colors import LogNorm
	import numpy as num
	from scipy.integrate import simps

	# Constants
	c_kms = 3.*10**5		# speed of light, km s-1
	c_cms = 3.0*10**10      # speed of light, cm s-1
	m_e = 9.11*10**-28      # electron mass, g
	ee = 4.80*10**-10       # electrostatic charge, esu

	sig0 = math.pi*ee*ee / (m_e*c_cms)      # classical cross section, cm^2 Hz
	sig0 = sig0*10**-8 / c_cms              # convert cm^2 Hz --> cm^2 / Ang


	# Spectrum file for HD 37061
	#~ ref_path = "C:\Users\keho8439\Documents\ASTR7500 - Space Mission Design\hmwk6"
	#~ filename = "\HD37061_UVdata_012915.txt"   # Make sure this is the fixed image using idl
	#~ plot_title = 'HD 37061 UV Spectrum - HST/STIS E140H'
	ref_path = "C:\\Users\\Keri Hoadley\\Documents\\CHESS\\CHESS2_flight\\"
	filename = "epsPer_1dspectr_v2_256res_Aeffcorr.dat"   # Make sure this is the fixed image using idl
	filename = "iue\\ePer_IUE.dat"   # Make sure this is the fixed image using idl
	plot_title = 'e Per: CHESS-2'

	# Open the fits file, returning the 2d image array in the specified resolution
	wave, flux, err = open_txt(ref_path+filename,plot_title)

	wave = num.array(wave)
	flux = num.array(flux)
	err = num.array(err)

	# Neutral carbon lines, in Angstroms
	CI_waves = [1334.53, 1335.71]#[1193.031, 1277.245]	# Lines within STIS spectrum
	f_osc = [0.129, 0.115]	#[0.0409, 0.0923]      # oscillator strength for each line
	dAng = 5.0		# +/- dAng around the CI lines

	CI_index = 0
	CI_wave = CI_waves[CI_index]
	CI_title = str(CI_wave)
	f_osc = f_osc[CI_index]

	print 'Cross section for C I line at',CI_title,'Ang = ',sig0*f_osc
	print

	# Find indices where wave is between CI_waves +/- dAng
	i_start = num.where( wave > CI_wave-dAng )
	i_start = i_start[0][0]
	i_end = num.where( wave < CI_wave+dAng )
	i_end = i_end[0]
	i_end = i_end[::-1]
	i_end = i_end[0]

	# Store all wave, flux, err in the absorption line
	abs_wave = wave[i_start:i_end]
	abs_flux = flux[i_start:i_end]
	abs_err = err[i_start:i_end]

	# Print around the absorption line
	pylab.errorbar(abs_wave,abs_flux,yerr=abs_err,errorevery=10,color='k',ecolor='k',elinewidth=3)
	pylab.plot((CI_wave-dAng, CI_wave+dAng),(0,0),'r--',linewidth=5)
	pylab.plot( (CI_wave,CI_wave),(2.0E-10,-2.0E-10),'g--',linewidth=5 )
	pylab.xlabel('Wavelength (Ang)')
	pylab.ylabel('Flux (ergs/cm^2/s/Ang)')
	pylab.xlim(CI_wave-dAng, CI_wave+dAng)
	#~ pylab.ylim(0,1.4E-10)
	pylab.title('HD 37061 - CI Absorption line at '+CI_title+' Ang')
	pylab.grid('on')


	# Fit a line/quadratic through the continuum level around the absorption line
	# to normalize C I line.
	cont_wave = num.concatenate( (abs_wave[0:len(abs_wave)/2], abs_wave[3*len(abs_wave)/4:len(abs_wave)]) )
	cont_flux = num.concatenate( (abs_flux[0:len(abs_wave)/2], abs_flux[3*len(abs_wave)/4:len(abs_wave)]) )

	# Use the "continuum" from CI_wave-dAng to CI_wave to normalize the spectrum
	a, b, c = num.polyfit(cont_wave,cont_flux, 2)

	# plot linefit over the absorption profile, to make sure line fit is ok.
	cont_fit = a*abs_wave**2 + b*abs_wave + c
	pylab.plot(abs_wave, cont_fit, 'b--', linewidth=3)
	pylab.show()

	# Divide spectrum by continuum fit to normalize
	norm_flux = abs_flux / cont_fit
	norm_err = abs_err / cont_fit

	pylab.errorbar(abs_wave,norm_flux,yerr=norm_err,errorevery=10,color='k',ecolor='k',elinewidth=3)
	pylab.plot((CI_wave-dAng, CI_wave+dAng),(0,0),'r--',linewidth=5)
	pylab.plot((CI_wave-dAng, CI_wave+dAng),(1,1),'r--',linewidth=5)
	pylab.plot( (CI_wave,CI_wave),(-0.5,1.5),'g--',linewidth=5 )
	pylab.xlabel('Wavelength (Ang)')
	pylab.ylabel('Normalized flux')
	pylab.xlim(CI_wave-dAng, CI_wave+dAng)
	pylab.ylim(-0.1,1.2)
	pylab.title('HD 37061 - Normalized CI Absorption line at '+CI_title+' Ang')
	pylab.grid('on')
	pylab.show()

	# Find shift in abs. lines from the lab wavlength
	# Do by converting abs_wave to velocity space.
	abs_vel = num.zeros(len(abs_wave))
	for i in range(0,len(abs_wave)):
		vel = (abs_wave[i] - CI_wave)/CI_wave
		abs_vel[i] = vel*c_kms

	# Plot again, in velocity space
	pylab.errorbar(abs_vel,norm_flux,yerr=norm_err,errorevery=10,color='k',ecolor='k',elinewidth=3)
	pylab.plot((-500,500),(0,0),'r--',linewidth=5)
	pylab.plot((-500,500),(1,1),'r--',linewidth=5)
	pylab.plot( (0,0),(-0.5,1.5),'g--',linewidth=5 )
	pylab.xlabel('Velocity (km/s)')
	pylab.ylabel('Normalized flux')
	pylab.xlim(-200,200)
	pylab.ylim(-0.1,1.2)
	pylab.title('HD 37061 - Normalized CI Absorption line at '+CI_title+' Ang')
	pylab.grid('on')
	pylab.show()

	# Line center (1277) at ~ 27.5 km s-1
	# Line cneter (1193) at ~ 20.0 km s-1

	v_shift = [0,0]	#[20.0, 27.5]		# km/s, for each line
	v_shift = v_shift[CI_index]

	# Fit a Gaussian absorption line to the spectra for the CI line
	# For CI_index = 0, make abs_vel = 1 for all abs_vel > 40.0
	# Center around line: Do +/- 50-70 km/s around abs line
	#~ if CI_index == 0:
		#~ for vv in range(0,len(abs_vel)):
			#~ if abs_vel[vv] > 39.0:
				#~ norm_flux[vv] = 1.0	# And keep the errors as they are
##        elif CI_index == 1:
##                for vv in range(0,len(abs_vel)):
##                        if abs_vel[vv] > 85.0:
##                                norm_flux[vv] = 1.0

	# Find indices where wave is between v_shift +/- 60 km/s
	dvel = 150.	#60.

	i_start = num.where( abs_vel > v_shift - dvel )
	i_start = i_start[0][0]
	i_end = num.where( abs_vel < v_shift + dvel )
	i_end = i_end[0]
	i_end = i_end[::-1]
	i_end = i_end[0]


	# Overplot fixed flux for abs line fit
	pylab.errorbar(abs_vel,norm_flux,yerr=norm_err,errorevery=6,color='b',ecolor='b',elinewidth=3,ls='-.',lw=3)
	#pylab.show()

	# Define abs_vel, norm_flux, and flux_err within i_start, i_end
	# for Gaussian fitting routine
	# Then plot to check absorption line and enough continuum
	# is in the view
	abs_vel = abs_vel[i_start:i_end]
	abs_wave = abs_wave[i_start:i_end]
	norm_flux = norm_flux[i_start:i_end]
	norm_err = norm_err[i_start:i_end]

	pylab.errorbar(abs_vel,norm_flux,yerr=norm_err,errorevery=5,color='k',ecolor='k',elinewidth=3)
	pylab.plot((-500,500),(0,0),'r--',linewidth=5)
	pylab.plot((-500,500),(1,1),'r--',linewidth=5)
	pylab.plot( (0,0),(-0.5,1.5),'g--',linewidth=5 )
	pylab.xlabel('Velocity (km/s)')
	pylab.ylabel('Normalized flux')
	#~ pylab.xlim(-40,80)
	#~ pylab.ylim(-0.05,1.05)
	pylab.title('HD 37061 - Normalized CI Absorption line at '+CI_title+' Ang')
	pylab.grid('on')
	pylab.show()

	# Now, use Gaussian fitting function to find a best-fit
	# absorption line to sum up all counts
	theta, gaussian, area_under_gaussian = onedgaussian(abs_vel,norm_flux,norm_err,v_shift)

	# Using the area_under_gaussian, find the equivalen with
	# in wavelength space
	# Try finding area using num.trapz and scipy.integrate.simps
	trapz_area = num.trapz(gaussian,dx=200)
	simps_area = simps(gaussian,dx=200)
	print "Area (trap) = ", trapz_area
	print "Area (simp) = ", simps_area
	equiv_width = area_under_gaussian           # W = (dwave) * 1.0 (normalized)
	equiv_width_ang = equiv_width**2*CI_wave/c_kms
	print
	print 'Area under the Gaussian: ',area_under_gaussian
	print 'Equivalent width: ',equiv_width**2.,' km/s --> ',equiv_width_ang,' Angstoms'

	# Solve for the column density: N_CI = equiv_width / (wave**2 * sig0 * f_osc)
	############ FOR tau << 1 (Linear Curve of Growth) ###############
	N_CI = equiv_width**2 / (sig0*f_osc*CI_wave*CI_wave)
	print
	print 'Column density of C I from absorption line at ',CI_title,' = ',N_CI,' cm-2'

	dv = theta[2]

	pylab.errorbar(abs_vel,norm_flux,yerr=norm_err,errorevery=5,color='k',ecolor='k',elinewidth=3)
	pylab.plot(abs_vel,gaussian,'b--',lw=3)
	pylab.plot((-500,500),(0,0),'r--',linewidth=5)	# plot 1, 0 normalized lines
	pylab.plot((-500,500),(1,1),'r--',linewidth=5)
	pylab.plot((dv+equiv_width**2.,dv+equiv_width**2.),(-0.5,1.5),'g--',linewidth=2)	# Plot +ve equiv width
	pylab.plot((dv-equiv_width**2.,dv-equiv_width**2.),(-0.5,1.5),'g--',linewidth=2)	# Plot -ve equiv width
	pylab.xlabel('Velocity (km/s)')
	pylab.ylabel('Normalized flux')
	#~ pylab.xlim(-40,80)
	#~ pylab.ylim(-0.01,1.05)
	pylab.title('Best-fit Gaussian absorption line over CI Absorption line at '+CI_title+' Ang')
	pylab.grid('on')
	pylab.show()







