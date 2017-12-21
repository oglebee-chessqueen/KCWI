### GEN. SCRIPT FOR TAKING OUT NAY CONTINUUM SOURCES IN KCWI N/S SKY ###
###
###
##
#
#
#
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy
import scipy.ndimage


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

	hr = int(ra0[0:2])
	min = int(ra0[3:5])
	sec = float(ra0[6:])

	ra0_deg = 360.*( (hr + (min/60.) + (sec/3600.)) / 24.)

	hr = int(dec0[1:3])
	min = int(dec0[4:6])
	sec = float(dec0[7:])

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



def write_fits(file,data):
	# First, read-in previous list of events
	try:
		hdulist=fits.open(file,mode='update')
		hdulist[0].data = data

		hdulist.flush()		# Should update file with new events array
		hdulist.close()
		print "Updating %s SUCCESS!" % (file)

	except:
		print "Update to %s FAILED" % (file)

	return



def scube_hist(data, wv, ra, dec):
	'''
	Create a histogram of the flattened data array by counts.
	Pull the median value (per wavelength) out, along with the Gaussian width.
	Create a mask around pixels that are many sigma away from
	'''
	fig = plt.figure()
	ax2 = fig.add_subplot(111)
	im1 = plt.imshow(numpy.sum(data,axis=0), vmin=0,vmax=16000, origin='lower',
						interpolation="none", cmap='Reds',
						extent=[ra[0],ra[-1],dec[0],dec[-1]])
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')


	max = numpy.max(data)
	print max
	print numpy.shape( numpy.ravel(data) )
	x_range = numpy.arange(0,max+1)
	y_range = numpy.zeros( int(max+1)+1 )
	print numpy.size(x_range), numpy.size(y_range)
	for i in numpy.ravel(data):
		i = int(i)
		y_range[i] += 1
	y_range2=numpy.log10(y_range)
	ymx = numpy.max(y_range2)
	for i in range(0, numpy.size(y_range2)):
		if  numpy.isnan(y_range2[i])=='True' or numpy.isinf(y_range2[i])=='True' :
			y_range2[i]=0
		elif  y_range2[i]==float('-inf') or y_range2[i]==float('inf'):
			y_range2[i]=0

	#~ plt.plot(x_range,y_range2,drawstyle='steps-mid')
	#~ #plt.plot(x_range_gauss, y_range_gauss, drawstyle='steps-mid',color='r')
	#~ #x0, sig = gaus(x_range_gauss, y_range_gauss)
	#~ #plt.xlim(0,max)
	#~ #plt.xlim(10000,30000)
	#~ plt.ylim(0,numpy.max(y_range2))
	#~ #plt.title(filename)
	#~ plt.show()

	data1 = data
	median = numpy.median(data[data>0])
	i_conx = numpy.where(numpy.sum(data,axis=0) > 15000)[0]
	i_cony = numpy.where(numpy.sum(data,axis=0) > 15000)[1]
	#~ data1[:,i_conx,i_cony] = median
	print i_conx, i_cony
	for w in range(0,len(wv)):
		med = numpy.median(data[w,:,:])
		data1[w,i_conx,i_cony] = med*1.1
		#~ d1 = data1[w,:,:]
		#~ mean = numpy.mean(data[w,:,:])
		#~ #print mean
		#~ d1[data[w,:,:]>mean/1.7] = med
		#~ data1[w,:,:] = d1
	print
	print median



	fig = plt.figure()
	ax2 = fig.add_subplot(111)
	im1 = plt.imshow(numpy.sum(data1,axis=0), vmin=0, vmax=16000, origin='lower',
						interpolation="none", cmap='Reds',
						extent=[ra[0],ra[-1],dec[0],dec[-1]])
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')
	plt.show()

	return data


path = 'C:\Users\Keri Hoadley\Documents\KCWI'
dir = '\\'
redux = 'redux'
date1 = '170619'
index1 = 261
filename_sky = 'kb'+date1+'_00%03i_scube.fits' % index1
filename_obj = 'kb'+date1+'_00%03i_ocube.fits' % index1
filename_int = 'kb'+date1+'_00%03i_icube.fits' % index1
file_sky = path+dir+date1+dir+redux+dir+filename_sky
file_obj = path+dir+date1+dir+redux+dir+filename_obj
file_int = path+dir+date1+dir+redux+dir+filename_int


data, waves, ra, dec, all_ra, all_dec = open_fits(file_sky,1)
data_o, waves_o, ra_o, dec_o, all_ra_o, all_dec_o = open_fits(file_obj,1)

# Creat a fixed sky cube to subtract from the object cube to create
# a new icube
data_sky = scube_hist(data, waves, all_ra, all_dec)

icube = data_o - (data_sky/3.)


fig = plt.figure()
ax2 = fig.add_subplot(111)
im1 = plt.imshow(numpy.sum(data_o[800:1300,:,:],axis=0), vmin=0, vmax=30000, origin='lower',
					interpolation="none", cmap='Greens',
					extent=[all_ra_o[0],all_ra_o[-1],all_dec_o[0],all_dec_o[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')

fig = plt.figure()
ax2 = fig.add_subplot(111)
im1 = plt.imshow(numpy.sum(icube[800:1300,:,:],axis=0), vmin=0, vmax=25000, origin='lower',
					interpolation="none", cmap='Blues',
					extent=[all_ra_o[0],all_ra_o[-1],all_dec_o[0],all_dec_o[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')
plt.show()


# save icube to fits file
#~ write_fits(file_int,icube)
