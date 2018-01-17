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
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
plt.rcParams["font.family"] = "serif"




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
ref_path = 'C:\Users\keho8439\Documents\CHESS\Alignments\chess2_vacuum_alignments'
wavefile = "\chess2_wave_soln_pre-skin.pkl"
#~ ref_path = 'C:\Users\keho8439\Documents\CHESS\Alignments\chess2_postflight_cals'
#~ wavefile = "\chess2_wave_soln_post-flight.pkl"
dats = pkl.load(open(ref_path+wavefile,'r'))

## Order Spectra is a list of all of the individual order spectra.
indv = dats["Order Spectra"]
err = num.sqrt(indv)		# Take Poisson errors only for now.
## Wavelength solution per order.
wave = dats["Order Wavelengths"]
## List of pixel values per order.
pixs = dats["Order Pixels"]

## Wavelength Fit Parms has a list containing the components of an n-polynomial
## fit describing the relationship between pixels and wavelength.
## We can read it out to do some manipulation on the function, if the
## spectrum is shifted. but we expect the same pixel/wavelength relation.
parms = dats["Wavelength Fit Parms"]


# Open H2 emission lamp spectra (simulated)
# File is set up with columns as: wavelength (Ang), flux(ergs/cm2/s/Ang), convolved flux(ers/cm2/s/Ang)
h2file = "C:\Users\keho8439\Documents\IDL code\kf_cos_fit\h2lamp.txt"
h2lines = open(h2file).readlines()

# Read through each line, and append colymn values to correct arrays
waveh2 = []
fluxh2 = []
convh2 = []

for line in h2lines:
	if line.startswith('#'):
		continue
	column = line.split()
	waveh2.append( float(column[0]) )
	fluxh2.append( float(column[1]) )
	convh2.append( float(column[2]) )

flux_mx = max(fluxh2)
fluxh2 = map(lambda x: 1000*(x/flux_mx) - 100, fluxh2)
