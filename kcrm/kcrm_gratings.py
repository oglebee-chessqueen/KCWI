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
alpha = [19.5, 16.0, 14.0, 11.0]						# RL AOI
beta = [48.63, 46.17, 43.55, 39.48, 37.50]	# RH4, KH v1
beta = [48.45, 46.00, 43.4, 39.2, 37.9]			# RH4, KH v1
beta = [26.5, 24.23, 21.09, 18.22]					# RM1
beta = [26.44, 24.2, 21.46, 18.69]					# RM2
beta = [6.66, 6.88, 7.56, 9.09]					# RM2
dbeta = 5.2			# +/- angle subtended by detector to intercept spectrum
# For resolution:
F = 13.6
D = 155		# mm
waves = numpy.arange(500.,1100.,1)		# wavelength coverage from 500 - 1100 nm

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
