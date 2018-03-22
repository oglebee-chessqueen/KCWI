### Cross-correlate KCWI images to overlay similar features ###
###
###
##
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick
import numpy
plt.rcParams["font.family"] = "serif"
plt.rcParams['axes.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'

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
	print len(args)

	if len(args) > 0:
		ra0 = hdulist[0].header['TARGRA']#['CD1_1']		# RA (at center?)
		dec0 = hdulist[0].header['TARGDEC']#['CD2_2']		# DEC (at center?)
	else:
		ra0 = hdulist[0].header['RA']		# RA (at center?)
		dec0 = hdulist[0].header['DEC']		# DEC (at center?)
	#~ print type(ra0)
	#~ print ra0[0:2], ra0[3:5], ra0[6:]
	hr = int(ra0[0:2])
	min = int(ra0[3:5])
	sec = float(ra0[6:])
	#~ print hr, min, sec

	ra0_deg = 360.*( (hr + (min/60.) + (sec/3600.)) / 24.)

	print ra0_deg

	#~ print type(ra0)
	print dec0
	print dec0[1:3], dec0[4:6], dec0[7:]
	hr = int(dec0[1:3])
	min = int(dec0[4:6])
	sec = float(dec0[7:])
	print hr, min, sec

	dec0_deg = ( (hr + (min/60.) + (sec/3600.)) )

	print dec0_deg

	#ra0_deg = hdulist[0].header['CRVAL1']	# zeropoint RA (at element 0) in degrees
	#dec0_deg = hdulist[0].header['CRVAL2']	# zeropoint DEC (at element 0)

	delta_ra = hdulist[0].header['CD1_1']#*3600.		# change in RA per pix (IN DEGREES)
	nra = numpy.size(data[0,0,:])
	ra_lim = (nra/2.)*delta_ra					# Extrema coverage to one end of the image
	all_ra = numpy.arange(ra0_deg-ra_lim, ra0_deg+ra_lim, delta_ra)	# array containing change in RA from center of image
	#all_ra = numpy.arange(-1*ra_lim, ra_lim, delta_ra)	# array containing change in RA from center of image

	delta_dec = hdulist[0].header['CD2_2']#*3600.	# change in DEC per pix (IN DEGREES)
	ndec = numpy.size(data[0,:,0])
	dec_lim = (ndec/2.)*delta_dec
	print dec_lim	# Extrema coverage to one end of the image
	all_dec = numpy.arange(dec0_deg-dec_lim, dec0_deg+dec_lim, delta_dec)	# array containing change in DEC from center of image
	#all_dec = numpy.arange(-1*dec_lim, dec_lim, delta_dec)	# array containing change in DEC from center of image

	print ra0, dec0
	print ra0_deg, dec0_deg
	print ra_lim,dec_lim
	print ra0_deg-ra_lim, ra0_deg+ra_lim
	print dec0_deg-dec_lim, dec0_deg+dec_lim
	print

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



def onclick2(event): #This is a function that will store the coordinates of each click
	global ix2, iy2
	ix2,iy2 = event.xdata, event.ydata
	print 'x = %.2f, y = %.2f' %(ix2, iy2)

	global radec
	radec.append((ix2,iy2))
	#~ print radec

	if len(radec) == 2 or len(radec) == 4:
		#~ fig.canvas.mpl_disconnect(cid)
		plt.close()

	return   	# background values at end points only


def kcwi_image_shift(data1,data2,dra1,ddec1,dra2,ddec2):
	'''
	Use to find the offset of one image compared to the other using
	delta ra/dec arrays.
	Inputs
	------
		data1 - 3D data cube of img 1
		data2 - 3D data cube of img 2, to move to overlay on img 1
		dra1, ddec1 - 1D arrays of delta RA/DEC of img1 from central RA/DEC
		dra1, ddec1 - 1D arrays of delta RA/DEC of img2 from central RA/DEC
	Outputs
	-------

	'''
	# We need to set a reference point, so plot each image and click on two
	# points where the two images have shared features.
	# 1. Plot img1, click on two points in the image (preferably stars),
	#		 store clicked points for further use later.
	fig = plt.figure()
	ax = fig.add_subplot(111)
	levels=[100,1000]
	lw=[2.5,3.5]
	plt.imshow(numpy.sum(data1,axis=0), vmin=0, vmax=500, origin='lower',
						interpolation="none", cmap='Blues',
						extent=[dra1[0],dra1[-1],ddec1[0],ddec1[-1]])
	CS = plt.contour(numpy.sum(data1,axis=0), levels,
									linewidths=lw, colors='k', corner_mask=True,
									extent=[dra1[0],dra1[-1],ddec1[0],ddec1[-1]])
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')
	ax.set_title('Click on two stars shared by the two images:')

	cid = fig.canvas.mpl_connect('button_press_event',onclick2)
	print dra1[0],dra1[-1],ddec1[5],ddec1[-1]
	plt.show()

	#~ print 'x1: ', radec[0][0], radec[1][0]
	#~ print 'y1: ', radec[0][1], radec[1][1]

	x1 = [radec[0][0], radec[1][0]]
	y1 = [radec[0][1], radec[1][1]]


	# 2. Plot img2, click on two points in the image (preferably stars),
	#		 store clicked points for further use later.
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.imshow(numpy.sum(data2,axis=0), vmin=0, vmax=2000, origin='lower',
						interpolation="none", cmap='Blues',
						extent=[dra2[0],dra2[-1],ddec2[0],ddec2[-1]])
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')
	ax.set_title('Click on two stars shared by the two images:')

	cid = fig.canvas.mpl_connect('button_press_event',onclick2)
	plt.show()
	#~ if len(radec) == 2:
		#~ fig.canvas.mpl_disconnect(cid)
		#~ plt.close()

	#~ print dra1[0],dra2[0]
	x2 = [radec[2][0], radec[3][0]]
	y2 = [radec[2][1], radec[3][1]]

	dx = x2[0]-x1[0]
	dy = y2[0]-y1[0]

	# 3. Determine the shift between the two stars in img1 and img2
	print 'dx: ', x2[0]-x1[0]
	print 'dy: ', y2[0]-y1[0]

	#~ dra2_0 = (dra2[0]/scale_dx)+dx/1.75
	#~ dra2_1 = (dra2[-1]/scale_dx)+dx/1.75
	#~ ddec2_0 = (ddec2[5]/scale_dy)+dy*2.95
	#~ ddec2_1 = (ddec2[-1]/scale_dy)+dy*2.95

	# Re-image img1 and img2 over one another, but with img2 shifted to
	# fit over img1
	levels=[100,1000]
	lw=[2.5,3.5]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	im1 = plt.imshow(numpy.sum(data2,axis=0), vmin=0, vmax=2500, origin='lower',
						interpolation="none", cmap='Reds',
						extent=[dra2[0]-dx,dra2[-1]-dx,ddec2[0]-dy,ddec2[-1]-dy])
	im2 = plt.imshow(numpy.sum(data1,axis=0), vmin=0, vmax=1000, origin='lower',
						interpolation="none", cmap='Blues', alpha=0.7,
						extent=[dra1[0],dra1[-1],ddec1[0],ddec1[-1]])
	CS = plt.contour(numpy.sum(data1,axis=0), levels,
									linewidths=lw, colors='k', corner_mask=True,
									extent=[dra1[0],dra1[-1],ddec1[0],ddec1[-1]])
	#~ cbar = plt.colorbar()
	#~ cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')
	plt.show()

	return



def plot_image(data1,data2,data3,data4,dra1,ddec1,dra2,ddec2,dra3,ddec3,dra4,ddec4):
	'''
	Use to find the offset of one image compared to the other using
	delta ra/dec arrays.
	Inputs
	------
		data1 - 3D data cube of img 1
		data2 - 3D data cube of img 2, to move to overlay on img 1
		dra1, ddec1 - 1D arrays of delta RA/DEC of img1 from central RA/DEC
		dra1, ddec1 - 1D arrays of delta RA/DEC of img2 from central RA/DEC
	Outputs
	-------

	'''
	# Plot img1 and img2 over one another, but with img2 shifted to
	# fit over img1
	labels = ['[OI] 6300','[NI] 5200','HeI 4471','[OII] 3729']#'HI 4861']
	colors = ['red','grey','blue','green']
	levels=[100,1000]
	levels2 = [1750,2500]
	lw=[2.5,3.5]
	fig = plt.figure(figsize=(10,10))
	#~ fig.set_rasterized(True)
	ax1 = fig.add_subplot(1,1,1)
	#~ ax1.set_rasterized(True)
	ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)', weight='bold', size='x-large')
	ax1.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)', weight='bold', size='x-large')
	plt.imshow(numpy.sum(data2,axis=0), vmin=200, vmax=400, origin='lower',
						interpolation="none", cmap='Reds',
						extent=[dra2[0],dra2[-1],ddec2[0],ddec2[-1]],label=labels[0])

	plt.imshow(numpy.sum(data1,axis=0), vmin=0, vmax=200, origin='lower',
						interpolation="none", cmap='Greys', alpha=0.8,
						extent=[dra1[0],dra1[-1],ddec1[0],ddec1[-1]],label=labels[0])

	plt.imshow(numpy.sum(data3,axis=0), vmin=500, vmax=1500, origin='lower',
						interpolation="none", cmap='Blues', alpha=0.7,
						extent=[dra3[0],dra3[-1],ddec3[0],ddec3[-1]],label=labels[0])

	plt.imshow(numpy.sum(data4[:,32:-30,11:-12],axis=0), vmin=0, vmax=190, origin='lower',
						interpolation="none", cmap='Greens', alpha=1,
						extent=[dra4[11],dra4[-12],ddec4[32],ddec4[-30]],label=labels[0])

	CS1 = plt.contour(numpy.sum(data1,axis=0), levels,
									linewidths=lw, colors='k', corner_mask=True,
									extent=[dra1[0],dra1[-1],ddec1[0],ddec1[-1]])
	patches = [ mpatches.Patch( color=colors[i], label=labels[i] ) for i in range(len(labels)) ]
	# put those patched as legend-handles into the legend
	plt.legend(handles=patches, bbox_to_anchor=(0.02, 0.98), loc=2, borderaxespad=0., fontsize='x-large')
	ax1.get_xaxis().get_major_formatter().set_useOffset(False)
	ax1.get_yaxis().get_major_formatter().set_useOffset(False)
	#~ plt.tight_layout()
	plt.savefig('m57_kcwi_images_combined.pdf',rasterized=True)
	plt.show()

	return




# Set a reference coordinate
# (central star)
ra_ref = 360.*( (18. + (53/60.) + (35.079/3600.)) / 24.)			# 18:53:35.079
dec_ref = 33. + (01/60.) + (45.03/3600.)		#+33:01:45.03

print ra_ref, dec_ref
print


path='C:\Users\Keri Hoadley\Documents\KCWI'
dir = '\\'
redux='redux'
date1='170620'	#'170619'
index1=64	#260
filename1='kb'+date1+'_00%03i_icube.fits' % index1
file1=path+dir+date1+dir+redux+dir+filename1

date2='170415'	#
index2=188
#~ date2='170619'	#
#~ index2=260
#~ date2='170412'	#
#~ index2=213
filename2='kb'+date2+'_00%03i_icube.fits' % index2
file2=path+dir+date2+dir+redux+dir+filename2

date3='170619'	#
index3=260
filename3='kb'+date3+'_00%03i_ocube.fits' % index3
file3=path+dir+date3+dir+redux+dir+filename3

date4='170412'	#
index4=212			# low-res
index5 = 213
filename4='kb'+date4+'_00%03i+00%03i_icubes_coadd.fits' % (index4,index5)
file4=path+dir+date4+dir+redux+dir+filename4
index6=210			# hi-res
index7 = 211
filename5='kb'+date4+'_00%03i+00%03i_icubes_coadd.fits' % (index6,index7)
file5=path+dir+date4+dir+redux+dir+filename5


data1, waves1, ra1, dec1, all_ra1, all_dec1 = open_fits(file1,1)	#N/S mixed up RA,DEC, so add 1 for args to know to switch them
data2, waves2, ra2, dec2, all_ra2, all_dec2 = open_fits(file2)
data3, waves3, ra3, dec3, all_ra3, all_dec3 = open_fits(file3,1)
data4, waves4, ra4, dec4, all_ra4, all_dec4 = open_fits(file4)
data5, waves5, ra5, dec5, all_ra5, all_dec5 = open_fits(file5)


# plot collapsed 1d spectra
spec1 = numpy.sum(data4,axis=(1,2))	# 170412
wv1 = waves4
spec2 = numpy.sum(data5,axis=(1,2))	# 170412
wv2 = waves5
spec3 = numpy.sum(data2,axis=(1,2))	# 170415
wv3 = waves2
spec4 = numpy.sum(data3,axis=(1,2))	# 170619
wv4 = waves3
spec5 = numpy.sum(data1,axis=(1,2))	# 170620
wv5 = waves1

fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(4,1,1)
#~ ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
ax1.set_ylabel(r'Flux (arb.)',weight='bold',size='x-large')
ax1.set_title(date4,weight='bold',size='x-large')
plt.semilogy(wv1,spec1,drawstyle='steps-mid',color='black',label='low-res')
plt.semilogy(wv2,spec2,drawstyle='steps-mid',color='crimson',label='high-res')
plt.legend()
plt.ylim(10**2,1.5*10**6)
plt.xlim(3500,5600)


ax2 = fig.add_subplot(4,1,2)
#~ ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
ax2.set_ylabel(r'Flux (arb.)',weight='bold',size='x-large')
ax2.set_title(date2,weight='bold',size='x-large')
plt.plot(wv3,spec3,drawstyle='steps-mid',color='black')
plt.xlim(5450,6400)
plt.ylim(10**2,2*10**5)
f = mtick.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
plt.gca().yaxis.set_major_formatter(mtick.FuncFormatter(g))
#~ ax2.get_yaxis().set_major_formatter(plt.LogFormatter(10,  labelOnlyBase=True))


ax3 = fig.add_subplot(4,1,3)
#~ ax1.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='x-large')
ax3.set_ylabel(r'Flux (arb.)',weight='bold',size='x-large')
ax3.set_title(date3,weight='bold',size='x-large')
plt.semilogy(wv4,spec4,drawstyle='steps-mid',color='black')
plt.ylim(10**3,3*10**6)
plt.xlim(4300,4800)

ax4 = fig.add_subplot(4,1,4)
ax4.set_xlabel(r'Wavelength ($\AA$)',weight='bold',size='x-large')
ax4.set_ylabel(r'Flux (arb.)',weight='bold',size='x-large')
ax4.set_title(date1,weight='bold',size='x-large')
plt.semilogy(wv5,spec5,drawstyle='steps-mid',color='black')
plt.ylim(4*10**3,2*10**4)
plt.xlim(4950,5450)

plt.tight_layout()
plt.show()


dra1 = ra_ref-ra1
ddec1 = (dec_ref-dec1)*10
dra2 = (ra_ref-ra2)/-3
ddec2 = dec_ref-dec2

dra3 = (ra_ref-ra3)*-560
ddec3 = (dec_ref-dec3)*40

dra4 = (ra_ref-ra4)/2
ddec4 = (dec_ref-dec4)/-1.5

print ra_ref-ra1, dec_ref-dec1
print ra_ref-ra2, dec_ref-dec2
print ra_ref-ra3, dec_ref-dec3
print
print dra1, ddec1
print dra2, ddec2
print dra4, ddec4
print dra3, ddec3


radec = []

c = 3.0e+5
dwv1_red = 3		# In case we only want to look at the blue/red region of a line
dwv1_blue = 3
line = 5199
data_wv1, interval_wv1 = narrowband(line, data1, waves1, dwv1_blue, dwv1_red, "n")
vobs_wv1 = ((waves1[interval_wv1] - line) / line)*c

dwv2_red = 5		# In case we only want to look at the blue/red region of a line
dwv2_blue = 5
line = 6300 	#4861 #4500
data_wv2, interval_wv2 = narrowband(line, data2, waves2, dwv2_blue, dwv2_red, "n")
data_wv2[data_wv2<0] = 0.
vobs_wv2 = ((waves2[interval_wv2] - line) / line)*c

dwv3_red = 3	#15		# In case we only want to look at the blue/red region of a line
dwv3_blue = 3 #15
line = 4471.3	#4640	#4471.5
data_wv3, interval_wv3 = narrowband(line, data3, waves3, dwv3_blue, dwv3_red, "n")
data_wv3[data_wv3<0] = 0.
vobs_wv3 = ((waves3[interval_wv3] - line) / line)*c

dwv4_red = 5	#2.5		# In case we only want to look at the blue/red region of a line
dwv4_blue = 5	#2.5
line = 3728 #4861 #4500
data_wv4, interval_wv4 = narrowband(line, data4, waves4, dwv4_blue, dwv4_red, "n")
data_wv4[data_wv4<0] = 0.
vobs_wv4 = ((waves4[interval_wv4] - line) / line)*c

ra_all1 = ((all_ra1-dra1) - ra_ref)*3600 - 1.129
ra_all2 = ((all_ra2-dra2) - ra_ref)*3600 - 1.129
ra_all3 = ((all_ra3-dra3) - ra_ref)*3600 - 1.129
ra_all4 = ((all_ra4-dra4) - ra_ref)*3600 - 1.129

dec_all1 = ((all_dec1-ddec1) - dec_ref)*3600 - 1.50766
dec_all2 = ((all_dec2-ddec2) - dec_ref)*3600 - 1.50766
dec_all3 = ((all_dec3-ddec3) - dec_ref)*3600 - 1.50766
dec_all4 = ((all_dec4-ddec4) - dec_ref)*3600 - 1.50766

# Set up dictionary with data, ra, dec as keys and arrays as values
data_dict ={'data': [data_wv1, data_wv2, data_wv3, data_wv4],
						'ra': [ra_all1, ra_all2, ra_all3, ra_all4],
						'dec': [dec_all1, dec_all2, dec_all3, dec_all4],}



plot_image(data_wv1,data_wv2,data_wv3,data_wv4,
					 ra_all1, dec_all1,
					 ra_all2, dec_all2,
					 ra_all3, dec_all3,
					 ra_all4, dec_all4)

#kcwi_image_shift(data_wv1,data_wv2,all_ra1,all_dec1,all_ra2,all_dec2)
