#
#
# CHESS-2 Spectral Resolution
#
#
# Read in: Pickle file with wavelength and spectra of orders extracted from
#          CHESS-2 echellogram.
#
# Goal: Calculate the resolving power of CHESS as a function of wavelength
#
#
#

import pickle as pkl
import pylab
import numpy as num
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq




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

	return abs(amp) * num.exp(-(xdata-cen)**2 /(2.*stdv**2))



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
	amp = num.zeros(0)
	cen = num.zeros(0)
	stdv = num.zeros(0)

	for i in range(0, len(params), 3): #This section is just unpacking the parameter array into amps, cens, and stdvs
		x = params[i]
		amp = num.append(amp,x)
		y = params[i+1]
		cen = num.append(cen,y)
		z = params[i+2]
		stdv = num.append(stdv,z)
	global storage #You may not need storage to be global so think about taking this part out. storage stores the data
	storage = [[0 for x in range(1)] for x in range(len(params)/3)] #from each iteration of the gaussian equation into
	for i in range(len(params)/3):#individual rows. So row one will be the gaussian solutions to the first peak and so on
		storage[i] = gaus(xdata,amp[i],cen[i],stdv[i])
	storage = num.asarray(storage)
	return sum(storage)


def onclick(event):
	"""
	Onclick is a function that will record the x and y values of points on the plot that have been clicked on.
	This is used in this script to record the guesses for the central wavelengths and amplitudes. I did not create
	this function, it was taken from stackoverflow.
	"""
	#~ cords = [] #This is an empty list which will store the x and y coordinates of each click on the graph
	#It's fine to keep this as a list because we won't be operating on it
	global ix, iy
	ix,iy = event.xdata, event.ydata
	print 'x = %.5f, y = %.2e' %(ix,iy) #This will print out the x and y values so you can check that no shifting occured

	global cords
	cords.append((ix,iy)) #Stores the x and y click in the array

	return




def onclick2(event): #This is a function that will store the coordinates of each click
	global ix2, iy2
	ix2,iy2 = event.xdata, event.ydata
	print 'x = %d, y = %d' %(ix2, iy2)

	global backg
	backg.append((ix2,iy2))
	#~ print backg

	if len(backg) == 2:
		fig.canvas.mpl_disconnect(cid)
		plt.close()

	return   	# background values at end points only



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

#~ # plot h2 lines to see
#~ pylab.plot(waveh2,fluxh2,'k')
#~ pylab.plot(waveh2,convh2,'b')
#~ pylab.show()


# Plot each order stitch individually (+H2 lamp flux spectra) and click on features
# to fit Gaussian profiles to (using Jacob's NGaussfit code)
# Then, re-plot with best-fit lines and ouput Gaussian parameters, including
# line centers and sigmas (to calculate FWHM).
index = 62
x = wave[index]
y = indv[index]
err = err[index]

'''
Function to calculate the multi-Gaussian fit to the order spectra of CHESS-2.
'''

# Plot the data and have user click on continuum across the order
# to subtract from data, for nicer Gaussian fit
# (though should not be necessary)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,y,'k',drawstyle='steps-mid')
ax.plot(waveh2,fluxh2,'b')
ax.set_xlim([min(x),max(x)])
ax.set_ylim([-20,max(y)+20.])
ax.set_title('Click on background value')

# Calls onclick2 function
backg = []

cid = fig.canvas.mpl_connect('button_press_event',onclick2)
plt.show()

print 'x: ', backg[0][0], backg[1][0]
print 'y: ', backg[0][1], backg[1][1]

bx = [backg[0][0], backg[1][0]]
by = [backg[0][1], backg[1][1]]

bfit = num.polyfit(bx,by,1)

# Calculate and subtract the scatter background from the data
backgr = bfit[0]*x + bfit[1]

# Before subtraction, show the linefit through the background, to make sure
# fit looks good by-eye
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,y,'k',drawstyle='steps-mid')
ax.plot(waveh2,fluxh2,'b')
ax.plot(x,backgr,'g--')
ax.set_xlim([min(x),max(x)])
ax.set_ylim([-20,max(y)+20.])
plt.show()

# Now, subtract the background
y = y - backgr # This subtracts out the average background value so
			   # the gaussian functions flow into the continuum.

# Next, user must click on ~ centers of emission features for
# Gaussian fitting to start
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,y,'k',drawstyle='steps-mid')
ax.plot(waveh2,fluxh2,'b')
ax.axhline(y=0,linestyle='--', color='g')
ax.set_title('Click on Emission Line Peaks and Estimate Widths')
ax.set_xlim([min(x),max(x)])
ax.set_ylim([-20,max(y)+20.])

# Calls onclick function
cords = []

cid = fig.canvas.mpl_connect('button_press_event',onclick)
plt.show()


# Now, set up all the gaussian parameters into arrays to read into gaussum
#~ stdv = raw_input('Guesses for Widths?') #Make sure to input in the format [#,#,#] with the []
#~ stdv = [0.05, 0.1, 0.05, 0.1, 0.05, 0.1, 0.05, 0.15, 0.05, 0.15, 0.05, 0.15, 0.1, 0.05, 0.15, 0.15, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

stdv = num.empty(len(cords))
stdv.fill(0.1)

amp=num.zeros(0)
cen = num.zeros(0)
for i in range(len(cords)):
	#This creates an array that stores the guessed line centers
	c = cords[i][0]
	cen = num.append(cen,c)

	#This creates an array that stores the guessed amplitudes
	a = cords[i][1]
	amp = num.append(amp,a)

#~ stdv2 = eval(stdv)
stdv2 = num.array(stdv)
wid = num.zeros(0)
for i in range(len(stdv2)):
	#This part takes in the guessed widths and converts it to a guessed sigma
	e = stdv2[i]
	wid = num.append(wid,e)
stdv = wid#/(2*num.sqrt(2*num.log(2)))


guess = num.zeros(0)
for i in range(len(amp)): #This is where the guessed values are all packed up into one array looking like
	xi = amp[i]            #[amp1,cen1,stdv1,amp2,cen2,...]
	guess = num.append(guess,xi)
	yi = cen[i]
	guess = num.append(guess,yi)
	zi = stdv[i]
	guess = num.append(guess,zi)

# This is where we take the gaussum function and use curvefit to find the best-fit
# multi-gaussian function across the order
popt, pcov = curve_fit(gaussum, x, y, p0=guess)#,sigma = err)# this optimizes the summed
															   # gaussian function
fit = gaussum(x, *popt)

#The optimized parameters for each emission peak is stored here in a 3xn array.
popt = popt.reshape((len(popt)/3,3))


# Calculates chi^2 of the fit over the order spectrum
chitest = lambda fobs,fexp: sum((fobs-fexp)**2/(fexp**2)) #I had to create my own chisquared function because scipy.stats' wasn't working
chi = num.zeros(0)
chireduced = num.zeros(0)
for i in range(len(popt)):
	opfit = gaus(x,popt[i][0],popt[i][1],popt[i][2])
	fexp = opfit + backgr
	fobs = x
	chil = chitest(fobs,fexp)
	chi = num.append(chi,chil)
	chireduced = num.append(chireduced,chil/(num.size(x)-3))
#print 'Line shape fits for order ',index
print 'chi^2: ', chi
print
print 'Reduced-chi^2: ', chireduced

# Convert reduced-chi^2 to some measure of error
fwhm_err = num.sqrt(chireduced)*2.0*num.sqrt(2.*num.log(2.))

fwhm_p = num.empty(len(popt))
res_p = num.empty(len(popt))
wav_p = num.empty(len(popt))

#~ # Write out to output file with order & readin pickle file name
#~ ref_path2 = ref_path+'\chess2_postflight_res'
#ref_path2 = 'C:\Users\keho8439\Documents\CHESS\Alignments\chess2_res'
#~ resfile = ref_path2+wavefile+'_orderindex=%03d_wave=%.2f-%.2f' %(index,min(x),max(x))
#~ fout = open(resfile+'.dat','w')
#~ fout.write('Order \t Amplitude \t Wavelength \t FWHM (um) \t err (um) \t Resolving Power \n')
#~ fout.write('----- \t --------- \t ---------- \t --------- \t -------- \t --------------- \n')

print 'Order \t Amplitude \t Wavelength \t FWHM (um) \t err (um) \t Resolving Power'
print '----- \t --------- \t ---------- \t --------- \t -------- \t ---------------'
for i in range(len(popt)):
	res_p[i] = popt[i][1]/(abs(popt[i][2])*2.0*num.sqrt(2.*num.log(2.)))
	wav_p[i] = popt[i][1]
	fwhm_p[i] = abs(popt[i][2])*2.0*num.sqrt(2.*num.log(2.))*10**4
	print '%03d \t %.1f \t\t %.4f \t %.4f \t %.4f \t %.1f' %(index,popt[i][0],popt[i][1],abs(popt[i][2]*10**4),fwhm_err[i],res_p[i])
	#~ # Write out to file, too
	#~ fout.write('%03d \t %.1f \t\t %.4f \t %.4f \t %.4f \t %.1f \n' %(index,popt[i][0],popt[i][1],abs(popt[i][2]*10**4),fwhm_err[i],res_p[i]))
#~ fout.close()



#This part plots all of the gaussian solutions as well as the summed solution

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.tick_params(labelsize=16)

for i in range(len(popt)):
	if i%2 == 0:
		plt.plot(x,gaus(x,popt[i][0],popt[i][1],popt[i][2]),'r')
	else:
		plt.plot(x,gaus(x,popt[i][0],popt[i][1],popt[i][2]),'b')
plt.plot(x,fit,'g',lw=3)
plt.plot(x,y,'k',drawstyle='steps-mid')
plt.plot(waveh2,fluxh2,'m')
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
#~ plt.title(r'Full LSF For Order %03d, Wavelength Range: %.2f - %.2f $\AA$' %(index,min(x),max(x)))
plt.xlabel(r'Wavelength ($\AA$)',weight='bold',size='x-large')
plt.ylabel('Counts',weight='bold',size='x-large')
plt.xlim(min(x),max(x))
plt.ylim(-120,max(y)+20.)
#~ ax.set_ylim([-20,max(y)+20.])
#plt.savefig(resfile+'.ps',orientation='landscape')
#print "Saved LSF Fit"

plt.show()



# Calculate resolving power as a function of fwhm.
# Plot results (wave v. fwhm, with R lines)
R1 = 25000
R2 = 33333
R3 = 50000
R4 = 66666
R5 = 100000

wave_R = num.arange(950.0,1700.0,10.0)	# wavlength array, Angstroms
fwhm_R1 = map(lambda w: (w/R1)*10**4, wave_R)		# R = l/dl --> dl = l/R (Ang--> um)
fwhm_R2 = map(lambda w: (w/R2)*10**4, wave_R)
fwhm_R3 = map(lambda w: (w/R3)*10**4, wave_R)
fwhm_R4 = map(lambda w: (w/R4)*10**4, wave_R)
fwhm_R5 = map(lambda w: (w/R5)*10**4, wave_R)

# Define starting x,y positions to put text about each R line
xstart = 1000
ystart = 410


plt.errorbar(wav_p,fwhm_p,yerr=fwhm_err,fmt='o',color='orange',ecolor='orange',ms=5,fillstyle='full')
plt.plot(wave_R,fwhm_R1,'k--',lw=2, label=str(R1))
plt.plot(wave_R,fwhm_R2,'b--',lw=2, label=str(R2))
plt.plot(wave_R,fwhm_R3,'g--',lw=2, label=str(R3))
plt.plot(wave_R,fwhm_R4,'r--',lw=2, label=str(R4))
plt.plot(wave_R,fwhm_R5,'m--',lw=2, label=str(R5))
plt.text(xstart,ystart,'R = 25000',color='black',rotation=17)
plt.text(xstart+10,ystart-100,'R = 33333',color='blue',rotation=14)
plt.text(xstart+20,ystart-210,'R = 50000',color='green',rotation=10.5)
plt.text(xstart+30,ystart-265,'R = 66666',color='red',rotation=7.5)
plt.text(xstart+40,ystart-320,'R = 100000',color='magenta',rotation=5)
plt.xlim(min(wave_R),max(wave_R))
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('FWHM (microns)')
plt.ylim(0,800)
plt.title('CHESS2: 100um pinhole Pre-Skin')
plt.show()
