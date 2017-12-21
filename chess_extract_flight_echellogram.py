#
#
# Pull 1D spectra stitches from CHESS-2 flight echellogram
#
#

import math
import astropy.io.fits as fits
import pickle as pkl
import numpy as num
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import rotate
from matplotlib.colors import LogNorm
from astropy.convolution import convolve, Box1DKernel

fits_file = 'C:\Users\Keri Hoadley\\Documents\\CHESS\\CHESS2_flight\\chess2_flight_data.fits'
res = 2
bin = 2048/res

width = 5/res		# Defined for 1D extraction


# Open the fits file,  pull out the data from t > 19400
# hdulist elements: 'X', 'Y', 'T', 'P' (not 'PHD')
#
hdulist = fits.open(fits_file)
photon_list = hdulist[1].data
hdulist.close()


x = photon_list['X']
y = photon_list['Y']
t = photon_list['T']



#~ tstart = 19200.	# 19430.
#~ tend = 19220.		# num.size(t)
#~ tlist = num.where(t >= tstart)
#~ start = tlist[0][0]
#~ tlist = num.where(t >= tend)
#~ end = tlist[0][0]

#~ x = x[start:end]
#~ y = y[start:end]



tstart = 19207.	# Very start of observation
tend = 19608.	# Very end of observation

tstart = 19395.
tend = 19608.
tlist = num.where(t >= tstart)
tlist_end = num.where(t <= tend)
start1 = tlist[0][0]
end1 = tlist_end[0][-1]


x = x[start1:end1]
y = y[start1:end1]



print num.size(x)

# Make an image of the science data, then rotate by theta
theta = -0.908173796791
# 2d histogram of photon list
img, xegdes, yedges = num.histogram2d(x,y,bins=bin,range=[[0,16384],[0,16384]])

#img = rotate(img,theta)
img = num.rot90(img)
img = rotate(img,theta,order=0,reshape=False)


# Trim detector edges; edge effects, extra counts, noise, etc.
def get_edges(arr,res):

	lis = arr.tolist()
	r_lis = list(reversed(lis))

	first_non0 = next((i for i, x in enumerate(lis) if x),res)
	last_non0 = next((i for i, x in enumerate(r_lis) if x),res)

	return len(lis)-first_non0-last_non0

trim=60/res
# Nick K: Hardcoded in for ease and repeatability.
xmid = bin/2 - 1
ymid = bin/2 - 1

xcentercol=img[int(ymid)]
ycenterrow=img[:,int(xmid)]

xarraydif = get_edges(xcentercol,bin)
yarraydif = get_edges(ycenterrow,bin)

avgrad=(xarraydif+yarraydif)/4


for i in range(0,bin):
	for j in range(0,bin):
		if num.sqrt((i-ymid)**2+(j-xmid)**2)>(avgrad-trim):
			img[i,j]=0




# PLot the resulting 2D image to pick out a spectral line to collapse to 1d spectrum
# plot the 2d histogram
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel('x',weight='bold',size='large')
ax1.set_ylabel('y',weight='bold',size='large')
#ax1.set_facecolor('black')
ax1.tick_params(labelsize=12)
#ax1.set_title('CLICK ON SPECTRAL LINE TO ANALYZE')

max = num.max(img)
print "Max counts in image: ",max

#~ plt.imshow(num.rot90(img,k=3), norm=LogNorm(vmin=5, vmax=max/2), cmap='viridis')#, origin='lower')
plt.imshow(num.rot90(img,k=3), vmin=5, vmax=15, cmap='hot')#'hot')#, origin='lower') # inferno - good
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts',weight='bold',size='large')
plt.show()





