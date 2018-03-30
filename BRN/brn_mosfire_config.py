# Create radial flux plot (or cts/sec) of BRN over radius.
# Basically recreate Fig 3a in Sahai+2010
#
#

import numpy
import numpy.polynomial.polynomial as poly
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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
	print hdulist.info()
	data = hdulist[0].data		# Define data cube (3D: wave, y, x)
	try:
		# For Galex fits files
		dra = hdulist[0].header['CDELT1']		# Degrees
		ddec = hdulist[0].header['CDELT2']
		ra0_deg = hdulist[0].header['CRVAL1']	# Degrees
		dec0_deg = hdulist[0].header['CRVAL2']
	except:
		# For Halpha fits file
		dra = hdulist[0].header['CD1_1']		# Degrees
		ddec = hdulist[0].header['CD2_2']
		ra0_deg = hdulist[0].header['CRVAL1']	# Degrees
		dec0_deg = hdulist[0].header['CRVAL2']

	nra = numpy.size(data[0,:])
	ra_lim = (nra/2.)*dra					# Extrema coverage to one end of the image
	all_ra = numpy.arange(ra0_deg-ra_lim, ra0_deg+ra_lim, dra)

	ndec = numpy.size(data[:,0])
	dec_lim = (ndec/2.)*ddec
	all_dec = numpy.arange(dec0_deg-dec_lim, dec0_deg+dec_lim, ddec)

	hdulist.close()

	return data, all_ra, all_dec



path='C:\Users\Keri Hoadley\\Documents\\BRN\\39BRN\\brn_fuv_nuv_ha.bck.dir'
dir = '\\'
frames=['Frame1','Frame2','Frame3']
fitsfile = ['BRN-fd-int.fits','BRN-nd-int.fits','BRN_Ha_final.fits']

fuv_file = path+dir+frames[0]+dir+fitsfile[0]
nuv_file = path+dir+frames[1]+dir+fitsfile[1]
ha_file = path+dir+frames[2]+dir+fitsfile[2]

# Central star ra, dec (Simbad: TYC 2597-735-1)
# 		16:48:37.4122 +35:12:09.452
ra_cs = 252.1558841666
dec_cs = 35.2026255555


# Read in file, get data from files
fuvdata, ra, dec = open_fits(fuv_file)		# data in counts per second
nuvdata, ra, dec = open_fits(nuv_file)		# data in counts per second
hadata, hara, hadec = open_fits(ha_file)			# data in counts per second
#~ hadata[hadata<-210] = 0

ra_rel = (ra - ra_cs)*3600
dec_rel = (dec - dec_cs)*3600
hara_rel = (hara - ra_cs)*3600
hadec_rel = (hadec - dec_cs)*3600

# print size of data images
print 'FUV image size: ',numpy.shape(fuvdata)
print 'NUV image size: ',numpy.shape(nuvdata)
print 'Ha image size: ',numpy.shape(hadata)


fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
plt.xlabel(r'$\Delta \alpha$ [$^{\prime \prime}$]',weight='bold')#,size='normal')
plt.ylabel(r'$\Delta \delta$ [$^{\prime \prime}$]',weight='bold')
plt.imshow(fuvdata, vmin=0, vmax=0.001, #origin='lower',
						interpolation="none", cmap='Blues',
						extent=[ra_rel[0],ra_rel[-1],dec_rel[0],dec_rel[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('GALEX/FUV cts/sec')#,weight='bold',size='large')
plt.imshow(nuvdata, vmin=0, vmax=0.01, #origin='lower',
						interpolation="none", cmap='Greens', alpha=0.25,
						extent=[ra_rel[0],ra_rel[-1],dec_rel[0],dec_rel[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel('GALEX/NUV cts/sec')#,weight='bold',size='large')
#~ plt.scatter(ra_cs, dec_cs, s=100, c='yellow', marker='*')
plt.scatter(0,0, s=100, c='yellow', marker='*')
# Draw a box around region to sum (or find average cts)
plt.imshow(hadata, vmin=-215, vmax=-190, #origin='lower',
						interpolation="none", cmap='Reds', alpha=0.45,
						extent=[hara_rel[0],hara_rel[-1],hadec_rel[0],hadec_rel[-1]])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'H$\alpha$ (cts)')#,weight='bold',size='large')


print 'RA/DEC of middle of image: ', hara[numpy.size(hara)/2], hadec[numpy.size(hadec)/2]

# Determine where in the data arrays the central star is
dec_step = 7.5/3600  #0.01	#0.00277
ra_off = -0.005


ra_neb = 128
dec_neb = -183
# MOSFIRE patch - around CS
rect = patches.Rectangle((0.35,-30.5),-0.35,2*30.5,
												 linewidth=2,edgecolor='darkgreen',facecolor='none')
# Add the patch to the Axes
ax1.add_patch(rect)
# MOSFIRE patch - nebula (6.1' x 3') (dec, ra)
rect = patches.Rectangle((ra_neb+90,dec_neb+50),-180,2*183,
												 linewidth=2,edgecolor='darkblue',facecolor='none')
# Add the patch to the Axes
ax1.add_patch(rect)



plt.show()
