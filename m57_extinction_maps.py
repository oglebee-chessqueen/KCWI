'''

	Create extinction maps using H-recombination lines observed
	along a portion of the Main Ring
	(not able to do this for any other regions -- no Hbeta)

'''


from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy
plt.rcParams["font.family"] = "serif"



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
	print type(ra0)
	print ra0[0:2], ra0[3:5], ra0[6:]
	hr = numpy.int(ra0[0:2])
	min = numpy.int(ra0[3:5])
	sec = numpy.float(ra0[6:])
	print hr, min, sec

	ra0_deg = 360.*( (hr + (min/60.) + (sec/3600.)) / 24.)

	print ra0_deg

	#~ print type(ra0)
	print dec0
	print dec0[1:3], dec0[4:6], dec0[7:]
	hr = numpy.int(dec0[1:3])
	min = numpy.int(dec0[4:6])
	sec = numpy.float(dec0[7:])
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


def narrowband(line, data, wv, dw, continuum_yn):
	'''
	Use to define a data array with limited wavelength coverage around line.
	continuum_yn is a string of length 1, either "y" or "n".
	If "y", then a continuum band of the same size as the narrowband image is
	subtracted (for example, to get rid of bright sources, if interested in
	diffuse media).
	If "n", then won't.
	Continuum is defined here as being 10lambda away from line
	'''
	ind_line = numpy.where((wv >= line-dw) & (wv <= line+dw))[0]		# good!
	if continuum_yn == "y":
		cline = line-10
		cont_line = numpy.where((wv >= cline-dw) & (wv <= cline+dw))[0]
		data_line = data[ind_line,:,:] - data[cont_line,:,:]  		# okay
	else:
		data_line = data[ind_line,:,:]

	return data_line, ind_line



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
index3=212
index4=213

# Use flux-calibrated data cubes
#~ filename1='kb'+date+'_00%03i_%s.fits' % (index1,int)
#~ file1=path+dir+date+dir+redux+dir+filename1

#~ filename2='kb'+date+'_00%03i_%s.fits' % (index2,int)
#~ file2=path+dir+date+dir+redux+dir+filename2

intfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,int)
file3 = path+dir+date+dir+redux+dir+intfile3
varfile3 = 'kb'+date+'_00%03i_%s.fits' % (index3,var)
vfile3 = path+dir+date+dir+redux+dir+varfile3

print file3

intfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,int)
file4 = path+dir+date+dir+redux+dir+intfile4
varfile4 = 'kb'+date+'_00%03i_%s.fits' % (index4,var)
vfile4 = path+dir+date+dir+redux+dir+varfile4

# Constants:
c = 3.0e5		# Speed of light: km/s
I_HgHb = 0.466
I_HdHb = 0.259
f_HgHb = 0.15		#0.13	#
f_HdHb = 0.22		#0.19	#

# List of important line locations for NGC 6720
HI = [4101.73, 4340.47, 4861.3, ]  	# delta, gamma, beta
HeI = [3888.65, 3972.02, 4026.19, 4471.48, 4713.15, 4921.93, 5015.68]
HeII = [4540, 4686, 5411]
NI = [5198.5, 5199]
NII = [5755]
FeIII = [4658]
OI = [5577,6300]
OII = [3726, 3729, 4267, 4661]
OIII = [4363.21, 4958.91, 5006.84, 5592.37]
NeIII = [3869, 3968]
NeIV = [4724, 4725]
SII = [4068, 4076]
ClII = [3679]
ClIII = [5517, 5537]
ClIV = [5323]
ArIII = [5192]
ArIV = [4711, 4739]
ArV = [4626]
KIV = [4511]
KV = [4123, 4163]
CaV = [3996, 5309]
UnKnown = []


# Read in file, get important data from files
#~ data, waves, ra, dec, all_ra, all_dec = open_fits(file1)
#~ data2, waves2, ra2, dec2, all_ra2, all_dec2 = open_fits(file2)
data3, waves3, ra3, dec3, all_ra3, all_dec3 = open_fits(file3)		# data in units erg cm-2 s-1
data4, waves4, ra4, dec4, all_ra4, all_dec4 = open_fits(file4)		# data in units erg cm-2 s-1
var3, varwv3 = open_fits_err(vfile3)
var4, varwv4 = open_fits_err(vfile4)

# Convert all_ra3 to a delta array away from central RA, and do same for DEC
# Multiply by 360 to convert to arcseconds
dra3 = (all_ra3 - ra3) * 3600.
ddec3 = (all_dec3 - dec3) * 3600
dra4 = (all_ra4 - ra4) * 3600
ddec4 = (all_dec4 - dec4) * 3600

# Define narrow bands around Hb, Hg, Hd
dw = 4

data3 = data4
waves3 = waves4
var3 = var4

line = HI[1]		# Hgamma
data_wv1, interval_wv1 = narrowband(line, data3, waves3, dw, "n")
data_wv1 = data_wv1[1:,:,:]
interval_wv1 = interval_wv1[1:]
var_wv1 = var3[interval_wv1,:,:]
vobs_wv1 = ((waves3[interval_wv1] - line) / line)*c
#~ data_wv1 = data_wv1/line
#~ for i in range(0,len(waves3[interval_wv1])):
	#~ data_wv1[i,:,:] = data_wv1[i,:,:]/waves3[interval_wv1[i]]

line = HI[2]		# Hbeta
data_wv2, interval_wv2 = narrowband(line, data3, waves3, dw, "n")
data_wv2 = data_wv2[1:-1,:,:]
interval_wv2 = interval_wv2[1:-1]
var_wv2 = var3[interval_wv2,:,:]
vobs_wv2 = ((waves3[interval_wv2] - line) / line)*c
#~ data_wv2 = data_wv2/line
#~ for i in range(0,len(waves3[interval_wv2])):
	#~ data_wv2[i,:,:] = data_wv2[i,:,:]/waves3[interval_wv2[i]]


line = HI[0]		# Hdelta
data_wv3, interval_wv3 = narrowband(line, data3, waves3, dw, "n")
data_wv3 = data_wv3[1:-1,:,:]
interval_wv3 = interval_wv3[1:-1]
var_wv3 = numpy.sqrt(var3[interval_wv3,:,:])
vobs_wv3 = ((waves3[interval_wv3] - line) / line)*c
#~ data_wv3 = data_wv3/line
#~ for i in range(0,len(waves3[interval_wv3])):
	#~ data_wv3[i,:,:] = data_wv3[i,:,:]/waves3[interval_wv3[i]]




# Compare
spectra3 = numpy.sum(data3[:,33:-28,12:-12], axis=(1,2))
err3 = numpy.sum(numpy.sqrt(var3[:,33:-28,12:-12]), axis=(1,2))
#~ plt.plot(waves3,spectra3,drawstyle='steps-mid',color='blue')
plt.errorbar(waves3,spectra3,yerr=err3,drawstyle='steps-mid',
						 color='blue',ecolor='blue',errorevery=10)

spectra4 = numpy.sum(data4[:,33:-28,12:-12], axis=(1,2))
err4 = numpy.sum(numpy.sqrt(var4[:,33:-28,12:-12]), axis=(1,2))
#~ plt.plot(waves4,spectra4,drawstyle='steps-mid',color='orange')
plt.errorbar(waves4,spectra4,yerr=err4,drawstyle='steps-mid',
						 color='orange',ecolor='orange',errorevery=10)

plt.plot(waves3[interval_wv2],numpy.sum(data_wv2,axis=(1,2)),'r-',lw=5)
plt.plot(waves3[interval_wv1],numpy.sum(data_wv1,axis=(1,2)),'r-',lw=5)
plt.plot(waves3[interval_wv3],numpy.sum(data_wv3,axis=(1,2)),'r-',lw=5)
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Counts')
#plt.show()





fig3 = plt.figure()
#~ fig3.set_size_inches(18.5, 10.5, forward=True)
ax3 = fig3.add_subplot(1,3,1)
ax3.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
#~ ax3.set_title(r'H$\gamma$/H$\beta$')
#~ CS1 = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									#~ linewidths=lw, colors='b', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
#~ CS2 = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									#~ linewidths=lw2, colors='r', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(numpy.sum(data_wv2[:,30:-25,10:-10],axis=0), origin='lower', #vmin=0.4, vmax=0.5,
						interpolation="None", cmap='cubehelix',
						extent=[dra3[10],dra3[-10],ddec3[30],ddec3[-25]])#, origin='lower') # inferno - good
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'H$\beta$',weight='bold',size='large')

ax3 = fig3.add_subplot(1,3,2)
ax3.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
#~ ax3.set_title(r'H$\delta$/H$\beta$')
#~ CS1 = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									#~ linewidths=lw, colors='b', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
#~ CS2 = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									#~ linewidths=lw2, colors='r', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(numpy.sum(data_wv1[:,30:-25,10:-10],axis=0), origin='lower', #vmin=0.2, vmax=0.3,
						interpolation="None", cmap='cubehelix',
						extent=[dra3[10],dra3[-10],ddec3[30],ddec3[-25]])#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel(r'H$\gamma$',weight='bold',size='large')

ax3 = fig3.add_subplot(1,3,3)
ax3.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
#~ ax3.set_title(r'H$\delta$/H$\beta$')
#~ CS1 = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									#~ linewidths=lw, colors='b', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
#~ CS2 = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									#~ linewidths=lw2, colors='r', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(numpy.sum(data_wv3[:,30:-25,10:-10],axis=0), origin='lower', #vmin=0.2, vmax=0.3,
						interpolation="None", cmap='cubehelix',
						extent=[dra3[10],dra3[-10],ddec3[30],ddec3[-25]])#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel(r'H$\delta$',weight='bold',size='large')
ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()






ratio = numpy.sum(data_wv1,axis=0) / numpy.sum(data_wv2,axis=0)
ratio2 = numpy.sum(data_wv3,axis=0) / numpy.sum(data_wv2,axis=0)
#~ ratio2[ratio2>I_HdHb] = I_HdHb
#~ ratio[ratio>I_HgHb] = I_HgHb

fig3 = plt.figure()
#~ fig3.set_size_inches(18.5, 10.5, forward=True)
ax3 = fig3.add_subplot(1,2,1)
ax3.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title(r'H$\gamma$/H$\beta$')
#~ CS1 = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									#~ linewidths=lw, colors='b', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
#~ CS2 = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									#~ linewidths=lw2, colors='r', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(ratio[30:-25,10:-10], origin='lower', vmin=0.4, vmax=0.5,
						interpolation="None", cmap='cubehelix',
						extent=[dra3[10],dra3[-10],ddec3[30],ddec3[-25]])#, origin='lower') # inferno - good
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'H$\gamma$/H$\beta$',weight='bold',size='large')

ax3 = fig3.add_subplot(1,2,2)
ax3.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title(r'H$\delta$/H$\beta$')
#~ CS1 = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									#~ linewidths=lw, colors='b', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
#~ CS2 = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									#~ linewidths=lw2, colors='r', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(ratio2[30:-25,10:-10], origin='lower', vmin=0.2, vmax=0.3,	#[33:-28,12:-12]
						interpolation="None", cmap='cubehelix',
						extent=[dra3[10],dra3[-10],ddec3[30],ddec3[-25]])#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel(r'H$\delta$/H$\beta$',weight='bold',size='large')
ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()



# Do extinction correction on ratios, using expected line flux ratios
# and relative interstellar extinction from Osterbock+2006

c_HgHb = numpy.zeros( shape=[numpy.size(ratio[:,0]),numpy.size(ratio[0,:])] )
c_HdHb = numpy.zeros( shape=[numpy.size(ratio[:,0]),numpy.size(ratio[0,:])] )

for i in range(0,numpy.size(ratio[:,0])):
	for j in range(0,numpy.size(ratio[0,:])):
		c_HgHb[i,j] = (-1./f_HgHb)*numpy.log10( (ratio[i,j]/I_HgHb) )#*0.434
		c_HdHb[i,j] = (-1./f_HdHb)*numpy.log10( (ratio2[i,j]/I_HdHb) )#*0.434
		if numpy.isnan(c_HgHb[i,j]):
			c_HgHb[i,j] = 0
		if numpy.isposinf(c_HgHb[i,j]):
			c_HgHb[i,j] = 0
		if numpy.isneginf(c_HgHb[i,j]):
			c_HgHb[i,j] = 0

		if numpy.isnan(c_HdHb[i,j]):
			c_HdHb[i,j] = 0
		if numpy.isposinf(c_HdHb[i,j]):
			c_HdHb[i,j] = 0
		if numpy.isneginf(c_HdHb[i,j]):
			c_HdHb[i,j] = 0

		#print c_HgHb[i,j], c_HdHb[i,j]

c_HgHb = c_HgHb[33:-28,12:-12]
c_HdHb = c_HdHb[33:-28,12:-12]
# plot c (extinction coefficients) for both Hd and Hg ratios
print numpy.median(c_HgHb[c_HgHb>0.05]), numpy.mean(c_HgHb[c_HgHb>0.05]), numpy.var(c_HgHb[c_HgHb>0.05])
print numpy.median(c_HdHb[c_HdHb>0.05]), numpy.mean(c_HdHb[c_HdHb>0.05]), numpy.var(c_HdHb[c_HdHb>0.05])




fig3 = plt.figure()
#~ fig3.set_size_inches(18.5, 10.5, forward=True)
ax3 = fig3.add_subplot(1,2,1)
ax3.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title(r'c$_{H \gamma /H \beta}$',weight='bold',size='large')
#~ CS1 = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									#~ linewidths=lw, colors='b', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
#~ CS2 = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									#~ linewidths=lw2, colors='r', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(c_HgHb, origin='lower', vmin=0.25, vmax=0.5,
						interpolation="None", cmap='cubehelix',
						extent=[dra3[12],dra3[-12],ddec3[33],ddec3[-28]])#, origin='lower') # inferno - good
#~ cbar = plt.colorbar()
#~ cbar.ax.set_ylabel(r'c$_{H \beta}$',weight='bold')

ax3 = fig3.add_subplot(1,2,2)
ax3.set_xlabel(r'$\Delta \alpha$ ($^{\prime \prime}$)',weight='bold',size='large')
ax3.set_ylabel(r'$\Delta \delta$ ($^{\prime \prime}$)',weight='bold',size='large')
#ax3.set_title(r'H$\beta$/[OIII]: '+ra0+' '+dec0)
ax3.set_title(r'c$_{H \delta /H \beta}$',weight='bold',size='large')
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)
#~ CS1 = plt.contour(numpy.sum(data_wv1[:,:,:],axis=0), levels,
									#~ linewidths=lw, colors='b', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
#~ CS2 = plt.contour(numpy.sum(data_wv2[:,:,:],axis=0), levels2,
									#~ linewidths=lw2, colors='r', corner_mask=True,
									#~ extent=[all_ra[0],all_ra[-1],all_dec[5],all_dec[-1]])
plt.imshow(c_HdHb, origin='lower', vmin=0.25, vmax=0.5,
						interpolation="None", cmap='cubehelix',
						extent=[dra3[12],dra3[-12],ddec3[33],ddec3[-28]])#, origin='lower') # inferno - good

cbar = plt.colorbar()
cbar.ax.set_ylabel(r'c$_{H \beta}$',weight='bold',size='large')
ax = plt.gca()
plt.tight_layout()



'''

	TO ADD IN:
	----------
		- Correct the orientation of the x-axis (or the orientation of the y-axis elements.
		- Change d_arcsec to reflect distance from central star.

	...
	... DONE!!!

'''

ra_star = 283.3962
dec_star = 33.0292


line_ra = dra3[12:-12]
line_dec = ddec3[33:-28]
line_dec = line_dec[::-1]

radec_slope = (line_ra[-1] - line_ra[0])/(line_dec[-1] - line_dec[0])
d_arcsec = numpy.zeros( len(line_dec) )

line_ra2 = numpy.zeros( len(line_dec) )
# Use line_dec, line_ra2 to find extinction values across diagonal of images
c_HgHb_line = numpy.zeros( len(line_dec) )
c_HdHb_line = numpy.zeros( len(line_dec) )
ii = 0
for i in range(0,len(line_ra2)):
	if line_ra[ii] > radec_slope*line_dec[i]:
		line_ra2[i] = line_ra[ii]
	else:
		ii += 1
		line_ra2[i] = line_ra[ii]
	#~ c_HgHb_line[i] = numpy.median( c_HgHb[i,ii:ii+5] )
	#~ c_HdHb_line[i] = numpy.median( c_HdHb[i,ii:ii+5] )
	c_HgHb_line[i] = c_HgHb[i,ii]
	c_HdHb_line[i] = c_HdHb[i,ii]
	d_arcsec[i] = numpy.sqrt( (all_ra3[ii+12] - ra_star)**2 + (all_dec3[i+33] - dec_star)**2 )*3600.



fig = plt.figure()
#~ fig3.set_size_inches(18.5, 10.5, forward=True)
ax3 = fig.add_subplot(1,1,1)
ax3.set_xlabel(r'Distance from Central Star ($^{\prime \prime}$)',weight='bold',size='large')
ax3.set_ylabel(r'c$_{H \beta}$',weight='bold',size='large')
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)
chghb, = plt.plot(d_arcsec[1:-2],c_HgHb_line[1:-2],drawstyle='steps-mid',
									label=r'H$\gamma$/H$\beta$')
chdhb, = plt.plot(d_arcsec[1:-2],c_HdHb_line[1:-2],drawstyle='steps-mid',
									color='black', label = r'H$\delta$/H$\beta$')
plt.legend(handles=[chghb, chdhb])
#~ plt.xlim([10,-10])

plt.show()
