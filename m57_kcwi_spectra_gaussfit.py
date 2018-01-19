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
	print numpy.shape(data)
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
	all_wv = all_wv[400:-300]
	print numpy.shape(data)

	return data, all_wv



def open_lineID(file,*args):
	'''
	Open the DAT file with line IDs of emission per observation.
	Grab the line ID + rest wavelength.
	'''
	hlines = open(file).readlines()

	# Read through each line, and append colymn values to correct arrays
	linewave = []
	lineID = []

	for line in lines:
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
		return numpy.sum(data4[:,args[0]:args[1],args[2]:args[4]],axis=(1,2))
	else:
		return "INCORRECT NUMBER OF ARGUMENTS - NO SPECTRUM CREATED"



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
#~ path='C:\Users\Keri Hoadley\Documents\KCWI'
#~ dir = '\\'
path='/home/keri/KCWI'
dir = '/'
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
data3, waves3 = open_fits(file3)		# data in units erg cm-2 s-1
data4, waves4 = open_fits(file4)		# data in units erg cm-2 s-1
var3, varwv3 = open_fits(vfile3)
var4, varwv4 = open_fits(vfile4)


# Open emission line template file:
# Important columns in file: Line ID (col 0), Lab Wave (col 2)
linesfile = date+"_212-213_lines.dat"
lineID, linewave = open_lineID(path+dir+linesfile)


spectra = numpy.sum(data4,axis=(1,2))

# bin the spectra by ~ 20 pixels?
# Use central wavelenth every 20 steps
# Take median (or avg) of 20 pixels considered
bin_size = 100
wave_bin = numpy.zeros( numpy.size(waves4)/bin_size )
flux_bin = numpy.zeros( numpy.size(waves4)/bin_size )
count = 0
for i in range(0, numpy.size(waves4)/bin_size):
	if count+bin_size > numpy.size(waves4):
		continue
	else:
		tmpwave = waves4[count:count+bin_size]
		wave_bin[i] = tmpwave[9]
		flux_bin[i] = numpy.median(spectra[count:count+bin_size])
	count += bin_size

# Plot spectra (log)
fig1 = plt.figure()
ax1 = fig1.add_subplot(2,1,1)
plt.semilogy(waves4, numpy.sum(data4,axis=(1,2)), drawstyle='steps-mid', lw=3)
plt.semilogy(wave_bin, flux_bin, 'ro', ms = 5)
plt.ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$ $\AA ^{-1}$)')
ax2 = fig1.add_subplot(2,1,2)
plt.plot(waves4, numpy.sum(data4,axis=(1,2)), drawstyle='steps-mid', lw=3)
plt.plot(wave_bin, flux_bin, 'ro', ms = 5)
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$ $\AA ^{-1}$)')
plt.show()
