'''
			 KCRM Gratings: Plots of various properties of the gratings,
			                based on pre-defined parameters
	Inputs:
		- line density (rho)
		- AOI/AOD (from Matt's IDL: kcrm_gratings.pro)
		------------------------
		- maybe tilt?
		- indices of refraction: fused silica, gelatin
		- modulation of gelatin index
	Things to plot/output:
		- central wavelength
		- min/max wavelength (assuming angle subtended by detector +/- 5.2 deg)
		- instantaneous bandpass (dwave)
		- Efficiency envelope vs. wavelength (per angle)
		- Channel envelope (over all AOI) vs. wavelength
		- Resolution v. wavelength
'''

import math
import numpy
import matplotlib.pyplot as plt

d2r = math.pi/180.
gtg = 'rh4'
# lines/micron
r = 1.495		# RH4, KH v1 (05/08/2018)
r = 1.435		# RH4, MM v1
r = 1.220		# RM1
r = 0.921		# RM2
r = 0.514 	# RL
alpha = [56., 53., 50., 45., 42.]						# RH AOI
alpha = [34.5, 32.0, 28.5, 25.0]						# RM AOI
alpha = [19.5, 18.0, 16.0, 14.0, 11.0, 9.0]						# RL AOI
beta = [48.63, 46.17, 43.55, 39.48, 37.50]	# RH4, KH v1
beta = [48.45, 46.00, 43.4, 39.2, 37.9]			# RH4, KH v1
beta = [26.5, 24.23, 21.09, 18.22]					# RM1
beta = [26.44, 24.2, 21.46, 18.69]					# RM2
beta = [6.66, 6.77, 6.88, 7.56, 9.09, 10.0]					# RL
tilt = 2.5			# Tilt of grating fringes from normal
nref = 1.45			# index of refraction of the gelatin
nair = 1.0			# index of refraction of air
thick = 7.07		# thickness of the gelatin layer, in microns(?)
glass = 1.5			# index of refraction of glass (likely not needed)
dn = 0.05			# index of refraction modulation through gelatin

dbeta = 5.2			# +/- angle subtended by detector to intercept spectrum
# For resolution:
F = 13.6
D = 155		# mm
waves = numpy.arange(500.,1100.,1)*0.001		# wavelength coverage from 500 - 1100 nm

#
# 1.
# Output central wavelength, wavelength min/max, and dwave per angle config.
print "For grating %s, rho = %0.3f lines/um:" % (gtg,r)
print "AOI \t wv cen (nm) \t wv max \t wv min \t dwv"
print "--- \t ----------- \t ------ \t ------ \t ---"
for i in range(0,len(alpha)):
	sina = math.sin(alpha[i]*d2r)
	sinb = math.sin(beta[i]*d2r)
	sinbp = math.sin((beta[i]+dbeta)*d2r)
	sinbm = math.sin((beta[i]-dbeta)*d2r)

print "%0.1f \t %0.2f \t %0.2f \t %0.2f \t %0.3f" % (alpha[i], 1000*(sina+sinb)/r, 1000*(sina+sinbp)/r, 1000*(sina+sinbm)/r, (((sina+sinbp)/r) - ((sina+sinbm)/r))*1000)


K = 2*3.14159*r
phi = (90.+tilt) / 180.*math.pi
delta = phi - math.pi/2
bta = 2*3.14159*nref / waves
beta = numpy.zeros( len(alpha) )
#
# 2.
# Deteermine efficiency per angle setting, and over al angles per wavelength
# A: Per angle setting (set alpha; use tilt to determine beta and efficiency)
for i in range(0,len(alpha)):
	# Define important parameters from Kogelnik approximation
	theta = math.asin( math.sin(alpha[i]*d2r)/nref )
	cr = math.cos(theta)
	cs = cr - K / bta *math.cos(phi)
	nus = math.pi*dn*thick / waves / numpy.sqrt(cs*cr)
	nup = nus*(-1*math.cos(2*(theta - phi)))
	xi = thick*( K*math.cos(theta - phi) - K*K*waves / (4*math.pi*nref)) / (2*cs)
	etas = ( numpy.sin( numpy.sqrt( nus**2 + xi**2 ) ) )**2 / (1 + (xi**2 / nus**2))
	etap = ( numpy.sin( numpy.sqrt( nup**2 + xi**2 ) ) )**2 / (1 + (xi**2 / nup**2))
	eta = 0.5*(etas + etap)
	wv0 = waves[eta==numpy.max(eta)]
	#print wv0
	beta[i] = math.asin( wv0*r - math.sin(alpha[i]*d2r + delta) ) / d2r
	#print "alpha = %0.1f \t beta = %0.1f" % (alpha[i], beta[i])

	plt.plot(waves*1000., eta)	#, '-', waves*1000., etas, ':', waves*1000., etap,'--')

#
# 3.
# At each wavelength, find the Braggs angle and determine efficiency at that angles
braggs_eff = numpy.zeros( numpy.size(waves) )
resol = numpy.zeros( numpy.size(waves) )
for i in range(0,numpy.size(waves)):
	#~ braggs = math.asin( nref * math.sin( math.asin(waves[i] * r/2./nref) - delta ) )
	braggs = math.asin(waves[i]*r/2./nref + delta)
	theta = braggs
	cr = math.cos(theta)
	cs = cr - K / bta[i] *math.cos(phi)
	nus = math.pi*dn*thick / waves[i] / numpy.sqrt(cs*cr)
	nup = nus*(-1*math.cos(2*(theta - phi)))
	xi = thick*( K*math.cos(theta - phi) - K*K*waves[i] / (4*math.pi*nref)) / (2*cs)
	etas = ( numpy.sin( numpy.sqrt( nus**2 + xi**2 ) ) )**2 / (1 + (xi**2 / nus**2))
	etap = ( numpy.sin( numpy.sqrt( nup**2 + xi**2 ) ) )**2 / (1 + (xi**2 / nup**2))
	braggs_eff[i] = 0.5*(etas + etap)
	# Also determine resolution here
	resol[i] = (F*D*waves[i]*r)/math.cos(theta)
# Print the max efficiency angle
# Print 70% encapsulated area (where eff >= 70%)
i_mx = numpy.where(braggs_eff == numpy.max(braggs_eff))[0]
i_70 = numpy.where(braggs_eff>0.7)[0]
i_70_min = i_70[0]
i_70_max = i_70[-1]
print "Peak Eff:"
print "wave = %0.2f \t eff = %0.6f" % (waves[i_mx]*1000, braggs_eff[i_mx])
print "70% Eff:"
print "wave = %0.2f \t eff = %0.6f" % (waves[i_70_min]*1000, braggs_eff[i_70_min])
print "wave = %0.2f \t eff = %0.6f" % (waves[i_70_max]*1000, braggs_eff[i_70_max])

plt.plot(waves*1000,braggs_eff,color='black',linewidth=5)
plt.xlim(500.,1100.)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Efficiency')
plt.axhline(0.7, color='black', linestyle='dashed', linewidth=3)
plt.axvline(530., color='purple', linestyle=':')
plt.axvline(1080., color='purple', linestyle=':')
# Add grating parms in text
#~ plt.text()
plt.show()


# resolution plot
plt.plot(waves*1000.,resol,color='black')
plt.xlim(500.,1100.)
plt.xlabel('Wavelength (nm)')
plt.ylabel(r'R ($\Delta \lambda$/$\lambda$)')
plt.axvline(530., color='purple', linestyle=':')
plt.axvline(1080., color='purple', linestyle=':')
plt.axhline(800, color='blue', linestyle='dashed', linewidth=3)
plt.show()
