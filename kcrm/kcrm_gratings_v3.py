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


def efficiency(aoi, nref, bta, dn, thick, waves):
	'''
	Determine theoretical efficiency for the set of grating parms.
	For a given angle.
	aoi = angle of indicence with grating, not gelatin (int)
	'''
	# Define important parameters from Kogelnik approximation
	cr = math.cos(aoi)
	cs = cr - K / bta *math.cos(phi)
	nus = math.pi*dn*thick / waves / numpy.sqrt(cs*cr)
	nup = nus*(-1*math.cos(2*(aoi - phi)))
	xi = thick*( K*math.cos(aoi - phi) - K*K*waves / (4*math.pi*nref)) / (2*cs)
	etas = ( numpy.sin( numpy.sqrt( nus**2 + xi**2 ) ) )**2 / (1 + (xi**2 / nus**2))
	etap = ( numpy.sin( numpy.sqrt( nup**2 + xi**2 ) ) )**2 / (1 + (xi**2 / nup**2))
	eff = 0.5*(etas + etap)
	return eff








import math
import numpy
import matplotlib.pyplot as plt

####### Grating constants ########
nref = 1.4			# index of refraction of the gelatin
nair = 1.0			# index of refraction of air
glass = 1.5			# index of refraction of glass (likely not needed)
dbeta = 5.2			# +/- angle subtended by detector to intercept spectrum
waves = numpy.arange(500.,1100.,1)*0.001		# wavelength coverage from 500 - 1100 nm
# For resolution:
F = 13.6
D = 155		# mm


####### Grating parameters #######
### RL1 ###
#~ gtg = 'KCRM RL1 v1'
#~ res_avg = 800
#~ r = 0.514 	# RL lines/micron
#~ #alpha = [19.5, 18.0, 16.0, 14.0, 11.0, 9.0]						# RL AOI
#~ #alpha = [18.0, 16.0, 13.79]						# RL AOI
#~ alpha = [16.0, 14.0, 11.5]						# RL AOI
#~ tilt = 0.96			# Tilt of grating fringes from normal
#~ dn = 0.05			# index of refraction modulation through gelatin
#~ thick = 7.07		# thickness of the gelatin layer, in microns(?)

### RM1 ###
#~ gtg = 'KCRM RM1 v1'
#~ res_avg = 1900
#~ r = 1.220		# RM1 lines/micron
#~ alpha = [32.0, 29.0, 26.0, 23.0]						# RM AOI
#~ tilt = 1.98	#2.5			# Tilt of grating fringes from normal
#~ dn = 0.04			# index of refraction modulation through gelatin
#~ thick = 8.42		# thickness of the gelatin layer, in microns(?)

### RM2 ###
#~ gtg = 'KCRM RM2 v1'
#~ res_avg = 1900
#~ r = 0.921  # RM2 lines/micron
#~ #alpha = [32.0, 28.5, 25.0, 22.5]						# RM AOI
#~ alpha = [32.0, 29.0, 26.0, 23.0]						# RM AOI
#~ tilt = 1.08	#2.5			# Tilt of grating fringes from normal
#~ dn = 0.05			# index of refraction modulation through gelatin
#~ thick = 9.47		# thickness of the gelatin layer, in microns(?)

### RH1 v1 ###
#~ gtg = 'KCRM RH1 v1'
#~ res_avg = 4600
#~ r = 2.42		# RH1 lines/micron
#~ alpha = [56., 53., 50.]#, 45.]						# RH AOI
#~ tilt = 3.27	#2.5			# Tilt of grating fringes from normal
#~ dn = 0.08			# index of refraction modulation through gelatin
#~ thick = 8.88		# thickness of the gelatin layer, in microns(?)

### RH42v1 ###
#~ gtg = 'KCRM RH2 v1'
#~ res_avg = 4600
#~ r = 2.03	# RH2 lines/micron
#~ alpha = [56., 53., 50.]#, 45.]						# RH AOI
#~ tilt = 4.05	#2.5			# Tilt of grating fringes from normal
#~ dn = 0.09			# index of refraction modulation through gelatin
#~ thick = 8.98		# thickness of the gelatin layer, in microns(?)

### RH3 v1 ###
#~ gtg = 'KCRM RH3 v1'
#~ res_avg = 4600
#~ r = 1.705		# RH3 lines/micron
#~ alpha = [56., 53., 50.]#, 45.]						# RH AOI
#~ tilt = 4.36	#2.5			# Tilt of grating fringes from normal
#~ dn = 0.09			# index of refraction modulation through gelatin
#~ thick = 10.48		# thickness of the gelatin layer, in microns(?)

### RH4 v1 ###
#~ gtg = 'KCRM RH4 v1'
#~ res_avg = 4600
#~ r = 1.435		# RH4, MM v1
#~ #alpha = [56., 53., 50., 45.]						# RH AOI
#~ alpha = [54., 52., 50.]#, 45.]						# RH AOI
#~ tilt = 2.9			# Tilt of grating fringes from normal
#~ dn = 0.09			# index of refraction modulation through gelatin
#~ thick = 13.61		# thickness of the gelatin layer, in microns(?)

### RH4 v2 ###
#~ gtg = 'KCRM RH4 v2'
#~ res_avg = 4600
#~ r = 1.5		# RH4, KH v1 (05/08/2018)
#~ alpha = [56., 53., 50.]#, 45.]						# RH AOI
#~ tilt = 2.2	#2.5			# Tilt of grating fringes from normal
#~ dn = 0.09			# index of refraction modulation through gelatin
#~ thick = 13.61		# thickness of the gelatin layer, in microns(?)
#~ print numpy.size(thick)


########### Grating Sets ############
##
## Medium: RM1 + RM2
gtg = 'KCRM RM v1'
res_avg = 1900
r_suite = numpy.array( [1.220, 0.921] )	# RM lines/micron
alpha = numpy. array( [32.0, 29.0, 26.0, 23.0] )					# RM AOI
tilt_suite = numpy.array( [1.98, 1.08] )	#2.5			# Tilt of grating fringes from normal
dn_suite = numpy.array( [0.04, 0.05] )		# index of refraction modulation through gelatin
thick_suite = numpy.array( [8.42, 9.47] )	# thickness of the gelatin layer, in microns(?)

##
## High: RH1 + RH2 + RH3 + RH4
gtg = 'KCRM RH v1'
res_avg = 4600
r_suite = numpy.array( [2.42, 2.03, 1.705, 1.435] )		# RH1 lines/micron
alpha = numpy.array( [56., 53., 50., 48.] )						# RH AOI
tilt_suite = numpy.array( [3.27, 4.05, 4.36, 2.9] )-1.0	#2.5			# Tilt of grating fringes from normal
dn_suite = numpy.array( [0.08, 0.09, 0.09, 0.09] )			# index of refraction modulation through gelatin
thick_suite = numpy.array( [8.88, 8.98,10.48, 13.61] )



beta = numpy.zeros( [len(r_suite),len(alpha)] )		# Determine output angle
wv0 = numpy.zeros( [len(r_suite),len(alpha)] )
eta = numpy.zeros( [len(r_suite),len(alpha),numpy.size(waves)] )
braggs_eff = numpy.zeros( [len(r_suite),numpy.size(waves)] )
resol = numpy.zeros( [len(r_suite),numpy.size(waves)] )

fig1 = plt.figure(figsize=(16,6))

# Loop through each grating in suite:
for g in range(0,numpy.size(r_suite)):
	r = r_suite[g]
	tilt = tilt_suite[g]
	dn = dn_suite[g]
	thick = thick_suite[g]

	# Define global variables to call in functions
	global d2r, K, phi, delta
	d2r = math.pi/180.
	K = 2*3.14159*r
	phi = (90.+tilt) / 180.*math.pi
	delta = phi - math.pi/2

	bta = 2*3.14159*nref / waves

	#
	# 1.
	# Output central wavelength, wavelength min/max, and dwave per angle config.
	#~ print "For grating %s, rho = %0.3f lines/um:" % (gtg,r)
	#~ print "AOI \t wv cen (nm) \t wv max \t wv min \t dwv"
#	~ print "--- \t ----------- \t ------ \t ------ \t ---"
	#~ for i in range(0,len(alpha)):
		#~ sina = math.sin(alpha[i]*d2r)
		#~ sinb = math.sin(beta[i]*d2r)
		#~ sinbp = math.sin((beta[i]+dbeta)*d2r)
		#~ sinbm = math.sin((beta[i]-dbeta)*d2r)

	#~ print "%0.1f \t %0.2f \t %0.2f \t %0.2f \t %0.3f" % (alpha[i], 1000*(sina+sinb)/r, 1000*(sina+sinbp)/r, 1000*(sina+sinbm)/r, (((sina+sinbp)/r) - ((sina+sinbm)/r))*1000)


#
# 2.
# Deteermine efficiency per angle setting, and over al angles per wavelength
# A: Per angle setting (set alpha; use tilt to determine beta and efficiency)
	for i in range(0,len(alpha)):
		theta = math.asin( math.sin(alpha[i]*d2r)/nref )
		eta[g,i] = efficiency(theta, nref, bta, dn, thick, waves)
		w0 = waves[eta[g,i]==numpy.max(eta[g,i])]
		wv0[g,i] = w0[0]
		#print wv0
		#~ beta[i] = math.asin( nref * math.sin( (math.asin(wv0*r/nref) - math.sin(theta + delta) ) ) ) / d2r
		beta[g,i] = math.asin( (wv0[g,i]*r - math.sin(alpha[i]*d2r)) ) / d2r
		print "wave: %0.1f , alpha = %0.1f \t beta = %0.1f" % (wv0[g,i]*1000,alpha[i], beta[g,i])

		plt.plot(waves*1000., eta[g,i])	#, '-', waves*1000., etas, ':', waves*1000., etap,'--')

	#
	# 3.
	# At each wavelength, find the Braggs angle and determine efficiency at that angles
	for i in range(0,numpy.size(waves)):
		#~ braggs = math.asin(waves[i]*r/2./nref + delta)
		if (nref * math.sin( math.asin(waves[i] * r/2./nref) + delta )) >= 1.0:
			continue
		else:
			braggs = math.asin( nref * math.sin( math.asin(waves[i] * r/2./nref) + delta ) )
			theta = math.asin( math.sin(braggs)/nref )
			#~ braggs = math.asin( nref * math.sin( waves[i] * r/2./nref + delta ) )
			braggs_eff[g,i] = efficiency(theta, nref, bta[i], dn, thick, waves[i])
			#~ print waves[i], braggs/d2r, braggs_eff[i]
			# Also determine resolution here
			resol[g,i] = (F*D*waves[i]*r)/math.cos(braggs)

	# Print the max efficiency angle
	# Print 70% encapsulated area (where eff >= 70%)
	i_mx = numpy.where(braggs_eff[g] == numpy.max(braggs_eff[g]))[0]
	brg_mx = math.asin( nref * math.sin( math.asin(waves[i_mx] * r/2./nref) + delta ) )/d2r	#math.asin(waves[i_mx]*r/2./nref + delta)/d2r
	i_70 = numpy.where(braggs_eff[g]>0.7)[0]
	if len(i_70) == 0:
		i_70 = numpy.where(braggs_eff[g] > numpy.max(braggs_eff[g])*0.7)[0]
	i_70_min = i_70[0]
	i_70_max = i_70[-1]
	wv_mn = waves[i_70_min]
	if wv_mn < 0.53:
		wv_mn = 0.53
	wv_mx = waves[i_70_max]
	if wv_mx > 1.08:
		wv_mx = 1.08
	print "Peak Eff:"
	print "wave = %0.2f \t Braggs = %0.2f \t eff = %0.6f" % (waves[i_mx]*1000,  brg_mx, braggs_eff[g,i_mx])
	print "70% Eff:"
	print "wave = %0.2f \t eff = %0.6f" % (wv_mn*1000, braggs_eff[g,i_70_min])
	print "wave = %0.2f \t eff = %0.6f" % (wv_mx*1000, braggs_eff[g,i_70_max])



	plt.plot(waves*1000,braggs_eff[g],color='black',linewidth=5)
	plt.xlim(500.,1100.)
	plt.ylim(0.0,1.1)
	plt.xlabel('Wavelength (nm)',fontsize=14)
	plt.ylabel('Theoretical Efficiency',fontsize=14)
	plt.axhline(0.7, color='darkblue', linestyle='dashed', linewidth=3)
	plt.axvline(530., color='purple', linestyle=':')
	plt.axvline(1080., color='purple', linestyle=':')
	plt.title(gtg,fontsize=20)
	#~ plt.grid('True')
	# Angles
	for k in range(0,len(alpha)):
		plt.text(wv0[g,k]*1000 - 10, braggs_eff[g,numpy.where(waves==wv0[g,k])[0]]*1.055,
					 r'$\alpha$ = %0.2f$^{\circ}$'%(alpha[k]), color='grey' )
		plt.text(wv0[g,k]*1000 - 10, braggs_eff[g,numpy.where(waves==wv0[g,k])[0]]*1.02,
					 r'$\beta$ = %0.2f$^{\circ}$'%(beta[g,k]), color='grey' )
	# Add grating parms in text
	if len(r_suite) == 1:
		plt.text(550,0.425, r'$\rho$ = %i lines/mm'%(r*1000.),
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(550,0.35, r'$\lambda_c$ = %i nm'%(waves[i_mx]*1000),
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(650,0.35, r'$\lambda_{min}$ = %i nm'%(wv_mn*1000),
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(550,0.275, r'd = %0.2f $\mu$m'%(thick),
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(650,0.275, r'$\lambda_{max}$ = %i nm'%(wv_mx*1000),
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(550,0.20, r'$\phi$ = %0.2f$^{\circ}$'%(tilt),
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(550,0.125, r'n = %0.2f'%(nref),
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(550,0.05, r'$\Delta$n = %0.2f'%(dn),
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
plt.show()

# Make channel curve for grating suite (mode):
channel_eff = numpy.zeros( numpy.size(waves) )	# 1D array, so store highest efficiency of all gratings
for i in range(0,len(waves)):
	if i < 20:
		channel_eff[i] = braggs_eff[0,i]
	else:
		channel_eff[i] = numpy.max(braggs_eff[:,i])

fig3 = plt.figure(figsize=(16,6))
plt.plot(waves*1000,channel_eff,color='black',linewidth=5)
plt.xlim(500.,1100.)
plt.ylim(0.0,1.1)
plt.xlabel('Wavelength (nm)',fontsize=14)
plt.ylabel('Theoretical Efficiency',fontsize=14)
plt.axhline(0.7, color='darkblue', linestyle='dashed', linewidth=3)
plt.axvline(530., color='purple', linestyle=':')
plt.axvline(1080., color='purple', linestyle=':')
plt.title(gtg,fontsize=20)


# resolution plot
fig2 = plt.figure(figsize=(16,6))
for g in range(0,len(r_suite)):
	i_mx = numpy.where(braggs_eff[g] == numpy.max(braggs_eff[g]))[0]
	brg_mx = math.asin( nref * math.sin( math.asin(waves[i_mx] * r/2./nref) + delta ) )/d2r	#math.asin(waves[i_mx]*r/2./nref + delta)/d2r
	i_70 = numpy.where(braggs_eff[g]>0.7)[0]
	if len(i_70) == 0:
		i_70 = numpy.where(braggs_eff[g] > numpy.max(braggs_eff[g])*0.7)[0]
	i_70_min = i_70[0]
	i_70_max = i_70[-1]
	wv_mn = waves[i_70_min]
	if wv_mn < 0.53:
		wv_mn = 0.53
	wv_mx = waves[i_70_max]
	if wv_mx > 1.08:
		wv_mx = 1.08

	plt.plot(waves[i_70_min:i_70_max]*1000.,resol[g,i_70_min:i_70_max],color='black')
	plt.xlim(500.,1100.)
	plt.title(gtg,fontsize=20)
	plt.xlabel('Wavelength (nm)',fontsize=14)
	plt.ylabel(r'R ($\Delta \lambda$/$\lambda$)',fontsize=14)
	plt.axvline(530., color='purple', linestyle=':')
	plt.axvline(1080., color='purple', linestyle=':')
	plt.axhline(res_avg, color='blue', linestyle='dashed', linewidth=3)

plt.show()
