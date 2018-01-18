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
	print wv0, wvlast, numpy.size(all_wv)

	# Limit x, y image area to only include "good" image area
	# (i.e., no blank space or edge pixels)
	# Same for wavelength space
	data = data[400:-300,33:-18,12:-12]
	all_wv = aa_wv[400:-300]

	return data, all_wv



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

print file3

intfile4 = 'kb'+date+'_00%03i_%s_extcorr.fits' % (index4,int)
file4 = path+dir+date+dir+redux+dir+intfile4
varfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,var)
vfile4 = path+dir+date+dir+redux+dir+varfile4


# Read in file, get data + waves from files
data3, waves3, ra3, dec3, all_ra3, all_dec3 = open_fits(file3)		# data in units erg cm-2 s-1
data4, waves4, ra4, dec4, all_ra4, all_dec4 = open_fits(file4)		# data in units erg cm-2 s-1
var3, varwv3 = open_fits_err(vfile3)
var4, varwv4 = open_fits_err(vfile4)


# Open emission line template file:
# Important columns in file: Line ID (col 0), Lab Wave (col 2)
linesfile = date+"_212-213_lines.dat"
lines = open(path+dir+linesfile).readlines()

# Read through each line, and append colymn values to correct arrays
waveh2 = []
fluxh2 = []
convh2 = []

for line in lines:
	if line.startswith('#'):
		continue
	column = line.split()
	waveh2.append( float(column[0]) )
	fluxh2.append( float(column[1]) )
	convh2.append( float(column[2]) )

flux_mx = max(fluxh2)
fluxh2 = map(lambda x: 1000*(x/flux_mx) - 100, fluxh2)
