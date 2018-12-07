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
rpath = "C:\\Users\\Keri Hoadley\\Documents\\KCRM\\Gratings\\"


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
## Low: RL1
gtg = 'KCRM RL v1'
gtg_ind = ['KCRM RL1']
res_avg = 800
r_suite = numpy.array( [0.514] )	# RL lines/micron
#alpha = numpy. array( [16.0, 14.0, 11.5] )# RL AOI
alpha = numpy. array( [21.2, 18.0, 16.0, 13.79] )# RL AOI
#tilt_suite = numpy.array( [0.98] )				# Tilt of grating fringes from normal
tilt_suite = numpy.array( [2.4] )				# Tilt of grating fringes from normal
dn_suite = numpy.array( [0.05] )		# index of refraction modulation through gelatin
thick_suite = numpy.array( [7.07] )	# thickness of the gelatin layer, in microns(?)
x0 = [520]
x1 = [1090]

## Low: RL1 (lowest; full  bandpass)
gtg = 'KCRM RL (lowest res)'
gtg_ind = ['KCRM RL-low']
res_avg = 575
r_suite = numpy.array( [0.331] )	# RL lines/micron
#alpha = numpy. array( [16.0, 14.0, 11.5] )# RL AOI
alpha = numpy. array( [11.0] )# RL AOI
#tilt_suite = numpy.array( [0.98] )				# Tilt of grating fringes from normal
tilt_suite = numpy.array( [0.0] )				# Tilt of grating fringes from normal
dn_suite = numpy.array( [0.05] )		# index of refraction modulation through gelatin
thick_suite = numpy.array( [7.77] )	# thickness of the gelatin layer, in microns(?)
x0 = [520]
x1 = [1090]

## Low: RL1 (highest; ~BL)
gtg = r'KCRM RL (R $\sim$ KCWI-BL)'
gtg_ind = ['KCRM RL-BL']
res_avg = 1000
r_suite = numpy.array( [0.604] )	# RL lines/micron
#alpha = numpy. array( [16.0, 14.0, 11.5] )# RL AOI
alpha = numpy. array( [24.2, 20.0, 15.2] )# RL AOI
#tilt_suite = numpy.array( [0.98] )				# Tilt of grating fringes from normal
tilt_suite = numpy.array( [3.8] )				# Tilt of grating fringes from normal
dn_suite = numpy.array( [0.05] )		# index of refraction modulation through gelatin
thick_suite = numpy.array( [7.51] )	# thickness of the gelatin layer, in microns(?)
x0 = [520]
x1 = [1090]

##
## Medium: RM1 + RM2
#~ gtg = 'KCRM RM v1'
#~ gtg_ind = ['KCRM RM1','KCRM RM2']
#~ res_avg = 1900
#~ r_suite = numpy.array( [1.220, 0.921] )	# RM lines/micron
#~ alpha = numpy. array( [32.0, 28.0, 24.0] )				#32.0, 	# RM AOI
#~ #alpha = numpy. array( [29.0, 26.0, 23.0] )				#32.0, 	# RM AOI
#~ tilt_suite = numpy.array( [2.5, 2.7] )	#2.5			# Tilt of grating fringes from normal
#~ #tilt_suite = numpy.array( [1.85, 0.93] )	#2.5			# Tilt of grating fringes from normal
#~ dn_suite = numpy.array( [0.04, 0.05] )		# index of refraction modulation through gelatin
#~ thick_suite = numpy.array( [8.42, 9.47] )	# thickness of the gelatin layer, in microns(?)
#~ x0 = [520, 690]
#~ x1 = [830, 1090]

##
## High: RH1 + RH2 + RH3 + RH4
#~ gtg = 'KCRM RH v1'
#~ gtg_ind = ['KCRM RH1','KCRM RH2', 'KCRM RH3','KCRM RH4']
#~ res_avg = 4600
#~ r_suite = numpy.array( [2.52, 2.10, 1.725, 1.45] )		# RH1 lines/micron
#~ alpha = numpy.array( [56., 52.5, 50., 47.5] )						# RH AOI
#~ #alpha = numpy.array( [56., 53., 50., 48.] )						# RH AOI
#~ tilt_suite = numpy.array( [1.0,1.0,1.0,1.0] )+1
#~ #tilt_suite = numpy.array( [3.39, 4.165, 4.44, 3.01] )-1.0	#2.5			# Tilt of grating fringes from normal
#~ dn_suite = numpy.array( [0.08, 0.09, 0.09, 0.09] )			# index of refraction modulation through gelatin
#~ thick_suite = numpy.array( [8.5, 9.08, 11.08, 13.18] )
#~ x0 = [520, 620, 760, 910]
#~ x1 = [660, 790, 960, 1090]



beta = numpy.zeros( [len(r_suite),len(alpha)] )		# Determine output angle
wv0 = numpy.zeros( [len(r_suite),len(alpha)] )
eta = numpy.zeros( [len(r_suite),len(alpha),numpy.size(waves)] )
braggs_eff = numpy.zeros( [len(r_suite),numpy.size(waves)] )
resol = numpy.zeros( [len(r_suite),numpy.size(waves)] )

#~ fig1 = plt.figure(figsize=(6*numpy.size(r_suite),8))
#~ fig1 = plt.figure(figsize=(16,6))

# Loop through each grating in suite:
for g in range(0,numpy.size(r_suite)):
	fig1 = plt.figure(figsize=(16,6))
	r = r_suite[g]
	tilt = tilt_suite[g]
	dn = dn_suite[g]
	thick = thick_suite[g]

	# make subplot for index in grating suite
	#~ ax = fig1.add_subplot(1,numpy.size(r_suite),g+1)

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
# Determine efficiency per angle setting, and over al angles per wavelength
# A: Per angle setting (set alpha; use tilt to determine beta and efficiency)
	for i in range(0,len(alpha)):
		# Add in file output per AOI here
		#~ print rpath+gtg_ind[g]+' angle = '+str(alpha[i])+'.dat'
		fout = open(rpath+gtg_ind[g]+' angle = '+str(alpha[i])+'.dat','w')

		theta = math.asin( math.sin(alpha[i]*d2r)/nref )
		eta[g,i] = efficiency(theta, nref, bta, dn, thick, waves)
		w0 = waves[eta[g,i]==numpy.max(eta[g,i])]
		wv0[g,i] = w0[0]
		#print wv0
		#~ beta[i] = math.asin( nref * math.sin( (math.asin(wv0*r/nref) - math.sin(theta + delta) ) ) ) / d2r
		beta[g,i] = math.asin( (wv0[g,i]*r - math.sin(alpha[i]*d2r)) ) / d2r
		print "wave: %0.1f , alpha = %0.1f \t beta = %0.1f" % (wv0[g,i]*1000,alpha[i], beta[g,i])
		# +/- wavelength range falling on detector
		sinbp = math.sin((beta[g,i]+dbeta)*d2r)
		sinbm = math.sin((beta[g,i]-dbeta)*d2r)
		# +/- wavelength range covered by +/- angular extent on detector
		wp = ( sinbp + math.sin(alpha[i]*d2r) ) / r
		wm = ( sinbm + math.sin(alpha[i]*d2r) ) / r
		print wm, wp
		print len(eta[g,i])

		i_wm = numpy.where((waves < wm+0.0005) & (waves > wm-0.0005))[0]
		i_wm = i_wm[0]
		i_wp = numpy.where((waves < wp+0.0005) & (waves > wp-0.0005))[0]
		i_wp = i_wp[0]
		print i_wm, i_wp

		for w in range(0,numpy.size(waves)):
			if w == 0:
				fout.write('# Wavelength \t Efficiency \n')
				fout.write('# ---------- \t ---------- \n')
				fout.write('%.1f \t %.4f \n' %(waves[w]*1000,eta[g,i,w]) )
			else:
				fout.write('%.1f \t %.4f \n' %(waves[w]*1000,eta[g,i,w]) )

		plt.plot(waves*1000., eta[g,i])	#, '-', waves*1000., etas, ':', waves*1000., etap,'--')
		det_eta = eta[g,i]
		plt.plot(waves[i_wm:i_wp]*1000., eta[g,i,i_wm:i_wp], color='r', lw=8)
		fout.close()
		#~ plt.show()

	#
	# 3.
	# At each wavelength, find the Braggs angle and determine efficiency at that angles
	# Add file output here
	#~ fout = open(rpath+gtg_ind[g]+'.dat','w')
	#~ fout.write('# Wavelength \t Efficiency \n')
	#~ fout.write('# ---------- \t ---------- \n')
	for i in range(0,numpy.size(waves)):
		#~ braggs = math.asin(waves[i]*r/2./nref + delta)
		if (nref * math.sin( math.asin(waves[i] * r/2./nref) + delta )) >= 1.0:
			continue
		else:
			braggs = math.asin( nref * math.sin( math.asin(waves[i] * r/2./nref) + delta ) )
			# Add condition: If AOI > 56.0 degrees, set to 56 degrees
			if (braggs*(180/math.pi)) > 56.0:
				braggs = 56.0 * (math.pi/180)
			theta = math.asin( math.sin(braggs)/nref )
			#~ braggs = math.asin( nref * math.sin( waves[i] * r/2./nref + delta ) )
			braggs_eff[g,i] = efficiency(theta, nref, bta[i], dn, thick, waves[i])
			#~ print waves[i], braggs/d2r, braggs_eff[i]
			# Also determine resolution here
			resol[g,i] = (F*D*waves[i]*r)/math.cos(braggs)
			# Write out to file, too
			#~ fout.write('%.1f \t %.4f \n' %(waves[i]*1000,braggs_eff[g,i]) )
	# Close file
	#~ fout.close()

	# Print the max efficiency angle
	# Print 70% encapsulated area (where eff >= 70%)
	i_mx = numpy.where(braggs_eff[g] == numpy.max(braggs_eff[g]))[0]
	brg_mx = math.asin( nref * math.sin( math.asin(waves[i_mx] * r/2./nref) + delta ) )/d2r	#math.asin(waves[i_mx]*r/2./nref + delta)/d2r
	if numpy.size(r_suite) == 4:
		cutoff = 0.5
	else:
		cutoff = 0.7
	i_70 = numpy.where(braggs_eff[g]>cutoff)[0]
	if len(i_70) == 0:
		i_70 = numpy.where(braggs_eff[g] > numpy.max(braggs_eff[g])*cutoff)[0]
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
	print "Resolution at Braggs: %0.1f" % ( (F*D*waves[i_mx]*r)/math.cos(brg_mx*d2r) )



	plt.plot(waves*1000,braggs_eff[g],color='black',linewidth=5)
	plt.xlim(500.,1100.)
	plt.xlim(x0[g],x1[g])
	plt.ylim(0.0,1.1)
	plt.xlabel('Wavelength (nm)',fontsize=18)
	#~ if g == 0:
	plt.ylabel('Theoretical Efficiency',fontsize=18)
	plt.axhline(0.7, color='darkblue', linestyle='dashed', linewidth=3)
	plt.axhline(0.5, color='darkblue', linestyle='dashed', linewidth=1)
	plt.axvline(530., color='purple', linestyle=':')
	plt.axvline(1080., color='purple', linestyle=':')
	plt.title(gtg_ind[g],fontsize=18)#(gtg,fontsize=20)
	#~ plt.grid('True')
	# Angles
	for k in range(0,len(alpha)):
		if abs(alpha[k] - beta[g,k]) < 6.5:
			clr = 'crimson'
		else:
			clr = 'darkgreen'#'grey'
		plt.text(wv0[g,k]*1000 - 10, braggs_eff[g,numpy.where(waves==wv0[g,k])[0]]*1.07,
					 r'$\alpha$ = %0.2f$^{\circ}$'%(alpha[k]), color=clr, fontsize = 12 )
		plt.text(wv0[g,k]*1000 - 10, braggs_eff[g,numpy.where(waves==wv0[g,k])[0]]*1.025,
					 r'$\beta$ = %0.2f$^{\circ}$'%(beta[g,k]), color=clr, fontsize = 12 )
	# Add grating parms in text
	dx = (x1[g] - x0[g])
	xx1 = x0[g]+20
	xx2 = xx1+dx/4
	if len(r_suite) > 0:	#== 1:
		plt.text(xx1,0.425, r'$\rho$ = %i lines/mm'%(r*1000.),	# 550
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(xx1,0.35, r'$\lambda_c$ = %i nm'%(waves[i_mx]*1000),	# 550
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(xx2,0.35, r'$\lambda_{min}$ = %i nm'%(wv_mn*1000),	# 750
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(xx1,0.275, r'd = %0.2f $\mu$m'%(thick),	# 550
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(xx2,0.275, r'$\lambda_{max}$ = %i nm'%(wv_mx*1000),	# 750
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(xx1,0.20, r'$\phi$ = %0.2f$^{\circ}$'%(tilt),	# 550
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(xx1,0.125, r'n = %0.2f'%(nref),					# 550
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)
		plt.text(xx1,0.05, r'$\Delta$n = %0.2f'%(dn),				# 550
				 bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 fontsize = 12)

	#~ elif len(r_suite) == 1:
		#~ plt.text(550,0.425, r'$\rho$ = %i lines/mm'%(r*1000.),
				 #~ bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 #~ fontsize = 12)
		#~ plt.text(550,0.35, r'$\lambda_c$ = %i nm'%(waves[i_mx]*1000),
				 #~ bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 #~ fontsize = 12)
		#~ plt.text(650,0.35, r'$\lambda_{min}$ = %i nm'%(wv_mn*1000),
				 #~ bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 #~ fontsize = 12)
		#~ plt.text(550,0.275, r'd = %0.2f $\mu$m'%(thick),
				 #~ bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 #~ fontsize = 12)
		#~ plt.text(650,0.275, r'$\lambda_{max}$ = %i nm'%(wv_mx*1000),
				 #~ bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 #~ fontsize = 12)
		#~ plt.text(550,0.20, r'$\phi$ = %0.2f$^{\circ}$'%(tilt),
				 #~ bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 #~ fontsize = 12)
		#~ plt.text(550,0.125, r'n = %0.2f'%(nref),
				 #~ bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 #~ fontsize = 12)
		#~ plt.text(550,0.05, r'$\Delta$n = %0.2f'%(dn),
				 #~ bbox=dict(facecolor='white', edgecolor='white', alpha=0.9),
				 #~ fontsize = 12)
	plt.show()



# Make channel curve for grating suite (mode):
channel_eff = numpy.zeros( numpy.size(waves) )	# 1D array, so store highest efficiency of all gratings
for i in range(0,len(waves)):
	if i < 20:
		channel_eff[i] = braggs_eff[0,i]
	else:
		channel_eff[i] = numpy.max(braggs_eff[:,i])

fig3 = plt.figure(figsize=(16,6))	#8,8))
#~ ax = fig3.add_subplot(2,1,1)
#~ ax = fig3.add_subplot(1,1,1)
plt.plot(waves*1000,channel_eff,color='black',linewidth=5)
plt.xlim(500.,1100.)
plt.ylim(0.0,1.1)
plt.xlabel('Wavelength (nm)',fontsize=14)
plt.ylabel('Theoretical Efficiency',fontsize=14)
plt.axhline(0.7, color='darkblue', linestyle='dashed', linewidth=3)
plt.axhline(0.5, color='darkblue', linestyle='dashed', linewidth=1)
plt.axvline(530., color='purple', linestyle=':')
plt.axvline(1080., color='purple', linestyle=':')
plt.title(gtg+': Channel Envelope',fontsize=18)


# resolution plot
fig2 = plt.figure(figsize=(16,6))
#~ ax = fig3.add_subplot(2,1,2)
for g in range(0,len(r_suite)):
	i_mx = numpy.where(braggs_eff[g] == numpy.max(braggs_eff[g]))[0]
	brg_mx = math.asin( nref * math.sin( math.asin(waves[i_mx] * r/2./nref) + delta ) )/d2r	#math.asin(waves[i_mx]*r/2./nref + delta)/d2r
	i_70 = numpy.where(braggs_eff[g]>0.68)[0]
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
	if g == 0:
		plt.xlabel('Wavelength (nm)',fontsize=14)
		plt.ylabel(r'Resolution ($\Delta \lambda$/$\lambda$)',fontsize=14)
		plt.title(gtg+': Resolution (70% Bandwidth)',fontsize=18)
	plt.axvline(530., color='purple', linestyle=':')
	plt.axvline(1080., color='purple', linestyle=':')
	plt.axhline(res_avg, color='blue', linestyle='dashed', linewidth=3)

plt.show()
