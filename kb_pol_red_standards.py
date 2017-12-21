#!/usr/bin/env python

import numpy as np
import os, sys, pyfits
from optparse import OptionParser
import scipy as sp
from time import time
from scipy import signal
import glob
import matplotlib.mlab as mlab
from scipy.stats import norm
from scipy.optimize import curve_fit
from glob import glob
from StringIO import StringIO
import math
from astropy.io import fits

plot_flag = True
try:
    import matplotlib.pylab as plt
except ImportError:
    plot_flag = False


data_dir = '/home/gillian/data/otherdata/kcwi/170619/'
log_dir = '/home/gillian/data/otherdata/kcwi/170517-polar/'

def loaddata(log):

	polar = 'Polar'
        imgs = []
	counts = []
	polang = []
	rotpdest = []

	f = open(log,'r')
	for line in f:
	    if line.find(polar) != -1:
                #print line
		idx = line.find('Polar')
		#print (np.array([list(line)])[0][idx+6:]).tostring()
                num = (np.array([list(line)])[0][0:3]).tostring()
                imgs.append(num)
	#print imgs


	for line in imgs:
            imagecube = 'redux/kb'+str(data_dir[-7:-1])+'_'+str(line).zfill(5)+'_icube.fits'
	    #print line
	    #print imagecube
            hdulist = pyfits.open(data_dir + imagecube)
            sci_arr = pyfits.getdata(data_dir + imagecube)[:,26:45,9:15] #box aperature around the star BD323739
	    sky_arr = pyfits.getdata(data_dir + imagecube)[:,26:45,17:23] #box aperature around the sky BD323739
	    sci_arr = pyfits.getdata(data_dir + imagecube)[:,29:47,9:15] #box aperature around the star Hiltner960
	    sky_arr = pyfits.getdata(data_dir + imagecube)[:,29:47,17:23] #box aperature around the sky Hiltner960
	    #sky_arr = pyfits.getdata(data_dir + imagecube)[:,2:21,9:17] #box aperature around the sky
	    counts.append(np.mean(np.mean(sci_arr,axis=2),axis=1)-np.mean(np.mean(sky_arr,axis=2),axis=1))

            sz = sci_arr.shape
 
	    #exptime = hdulist[0].header['EXPTIME']
	    polang.append(round(hdulist[0].header['CALLANG']))
            rotpdest.append((hdulist[0].header['ROTPDEST']))
	    target = (hdulist[0].header['TARGNAME']).replace(" ","_")
            dw = hdulist[0].header['CD3_3']
	    wv0 = hdulist[0].header['CRVAL3']
            wp0 = hdulist[0].header['CRPIX3']
	    wg0 = hdulist[0].header['WAVGOOD0']
    	    wg1 = hdulist[0].header['WAVGOOD1']
 	    hdulist.close()

	mean_rotpdest = np.mean(rotpdest)
	#print mean_rotpdest
	allwaves = np.array([(np.arange(0,(sz[0]-wp0),1))*dw + wv0])
	#print wg0, wg1
	print polang
	idxw = (allwaves >= wg0)*(allwaves <= wg1) 
        allwaves = allwaves[idxw]
	#print np.array(counts)[:,idxw[0]]
	
	indata = log,imgs,num,counts,idxw,allwaves,polang,mean_rotpdest,target

	return (indata)

def stokes(indata):

	logs,imgs,num,counts,idxw,allwaves,polang,rotpdest,target = indata

	counts = np.array(counts)[:,idxw[0]]

	print counts.shape

	tk = 0.849475
	lk = 0.0010602
        Ak = (1./2.)*tk*(1+lk)
        ek = (1-lk)/(1+lk)        
	

	kidx = 5
	kang = 15*kidx
	angle1 = 0+13*kidx
        angle2 = 4+13*kidx
        angle3 = 8+13*kidx

	I = (1./(3.*Ak))*(counts[angle1,:] + counts[angle2,:] + counts[angle3,:])
	Q = (1./(3.*Ak))*(2.*counts[angle1,:] - counts[angle2,:] - counts[angle3,:])
	U = (1./(math.sqrt(3.)*ek*Ak))*(counts[angle2,:] - counts[angle3,:])

	I_mean = (1./(3.*Ak))*(np.mean(counts[angle1,:]) + np.mean(counts[angle2,:]) + np.mean(counts[angle3,:]))
	Q_mean = (1./(3.*Ak))*(2.*np.mean(counts[angle1,:]) - np.mean(counts[angle2,:]) - np.mean(counts[angle3,:]))
	U_mean = (1./(math.sqrt(3.)*ek*Ak))*(np.mean(counts[angle2,:]) - np.mean(counts[angle3,:]))

	#print I,Q,U

	#print I.shape, Q.shape, U.shape

        dolp = 100*np.sqrt(Q**2 + U**2)/I
        #Is U or Q supposed to have a negative sign for whatever reflections might be taking place from sky
        pa = (180./math.pi)*(1./2.)/(np.arctan2(U,Q)) % 180

	dolp_mean = 100*np.sqrt(Q_mean**2 + U_mean**2)/I_mean
        #Is U or Q supposed to have a negative sign for whatever reflections might be taking place from sky
        pa_mean = (180./math.pi)*(1./2.)/(np.arctan2(U_mean,Q_mean)) % 180

        #print allwaves_flatb.shape
        print dolp
        print pa

	print dolp_mean, pa_mean

	indata = logs,counts,dolp,pa,angle1,angle2,angle3,idxw,allwaves,kang

	return (indata)


def singlestokes(indata):

	log,imgs,num,counts,idxw,allwaves,polang,mean_rotpdest,target = indata

	counts = np.array(counts)[:,idxw[0]]

	tk = 0.849475
	lk = 0.0010602
        Ak = (1./2.)*tk*(1+lk)
        ek = (1-lk)/(1+lk)   
        dolp_mean = []
        pa_mean = []    
        idxmw = (allwaves >= 4600)*(allwaves <= 4800)
        #print allwaves[idxmw]
	
	pol0 = np.where(np.array(polang) == 0)
	pol60 = np.where(np.array(polang) == 60)
	pol120 = np.where(np.array(polang) == 120)
	angle1 = pol0[0][0]
	angle2 = pol60[0][0]
	angle3 = pol120[0][0]

	I = (1./(3.*Ak))*(counts[angle1,:] + counts[angle2,:] + counts[angle3,:])
	Q = (1./(3.*Ak))*(2.*counts[angle1,:] - counts[angle2,:] - counts[angle3,:])
	U = (1./(math.sqrt(3.)*ek*Ak))*(counts[angle2,:] - counts[angle3,:])

	I_mean = (1./(3.*Ak))*(np.mean(counts[angle1,idxmw]) + np.mean(counts[angle2,idxmw]) + np.mean(counts[angle3,idxmw]))
	Q_mean = (1./(3.*Ak))*(2.*np.mean(counts[angle1,idxmw]) - np.mean(counts[angle2,idxmw]) - np.mean(counts[angle3,idxmw]))
	U_mean = (1./(math.sqrt(3.)*ek*Ak))*(np.mean(counts[angle2,idxmw]) - np.mean(counts[angle3,idxmw]))

	#print I_mean,Q_mean,U_mean
	#print I.shape, Q.shape, U.shape

	dolp = 100*np.sqrt(Q**2 + U**2)/I
	#Is U or Q supposed to have a negative sign for whatever reflections might be taking place from sky
	pa = (180./math.pi)*(1./2.)/(np.arctan2(U,Q)) % 180

	dolp_mean.append(100*np.sqrt(Q_mean**2 + U_mean**2)/I_mean)
	#Is U or Q supposed to have a negative sign for whatever reflections might be taking place from sky
	pa_mean.append((180./math.pi)*(1./2.)/(np.arctan2(U_mean,Q_mean)) % 180)

        #print allwaves_flatb.shape
        #print dolp
        #print pa

	#print dolp_mean
        #print pa_mean

	indata = log,counts,dolp,pa,angle1,angle2,angle3,idxw,allwaves,dolp_mean,pa_mean,mean_rotpdest,target

	return (indata)

def standards():

	glob_pattern = log_dir + 'kpol-170619-140916-*.dat'
	logs = sorted(set(glob(glob_pattern)))
	#print logs

	sky_dolp = []
	sky_pa = []
	skypa = []
	rotpdest = []

        for log in logs:
	    log,counts,dolp,pa,angle1,angle2,angle3,idxw,allwaves,dolp_mean,pa_mean,mean_rotpdest,target = singlestokes(loaddata(log))
	    log_idx_s = log.find('skypa-')
	    log_idx_e = log.find('.dat')
            skypa.append(int(log[log_idx_s+6:log_idx_e]))
	    sky_dolp.append(dolp_mean)
	    sky_pa.append(pa_mean)
	    rotpdest.append(mean_rotpdest)

	print sky_dolp
	print sky_pa
	print skypa
	print rotpdest

	indata = log,counts,dolp,pa,angle1,angle2,angle3,idxw,allwaves,dolp_mean,pa_mean,rotpdest,sky_dolp,sky_pa,skypa,target
	return (indata)

def plot_pol(indata):

	logs,counts,dolp,pa,angle1,angle2,angle3,idxw,allwaves,kang = indata

	if plot_flag:
		plt.close("all")
		plt.clf()
		fig = plt.figure(figsize=(20,15))
		plt.rc('text', usetex=True)
		plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})
		plt.xlabel(r"\textbf{Wavelength [nm]}",fontsize=25)
		plt.ylabel(r"\textbf{DoLP [\%]}",fontsize=25)
		plt.grid(b=True, which='major', color='0.75', linestyle='--')
		plt.grid(b=True, which='minor', color='0.75', linestyle='--')
		plt.tick_params(axis='x', labelsize=25)
		plt.tick_params(axis='y', labelsize=25)
		#plt.axis([0,np.max(n_log),0,bins[bins.size-1]])
		
		axes = plt.gca()
		axes.set_ylim([0,10])
		#axes.set_xlim([10**0,10**5])
		plt.plot(allwaves, dolp, color='darkgreen',linestyle='-',linewidth=1,label='Degree of Linear Polarisation')
		plt.legend(loc=1,prop={'size':20})
		#figtext(.65, .15,'Pre-amp = %0.4f e-/ADU.\nReadnoise = %0.4f e-' % (preamp_gain, readnoise),fontsize=20)
		fig_name = logs.replace('.log','_dolp_'+str(kang))
		fn = '/home/gillian/data/otherdata/kcwi/figures/'
		fig_dir = os.path.dirname(fn)
		if not os.path.exists(fig_dir):# create data directory if needed
		    os.makedirs(fig_dir)
		plt.savefig(fn + fig_name + '.png')


	if plot_flag:
		plt.close("all")
		plt.clf()
		fig = plt.figure(figsize=(20,15))
		plt.rc('text', usetex=True)
		plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})
		plt.xlabel(r"\textbf{Wavelength [nm]}",fontsize=25)
		plt.ylabel(r"\textbf{PA [$^\circ$]}",fontsize=25)
		plt.grid(b=True, which='major', color='0.75', linestyle='--')
		plt.grid(b=True, which='minor', color='0.75', linestyle='--')
		plt.tick_params(axis='x', labelsize=25)
		plt.tick_params(axis='y', labelsize=25)
		#plt.axis([0,np.max(n_log),0,bins[bins.size-1]])
		
		axes = plt.gca()
		#axes.set_ylim([0,50])
		#axes.set_xlim([10**0,10**5])
		plt.plot(allwaves, pa, color='navy',linestyle='-',linewidth=1,label='Polarisation Angle')
		plt.legend(loc=1,prop={'size':20})
		#figtext(.65, .15,'Pre-amp = %0.4f e-/ADU.\nReadnoise = %0.4f e-' % (preamp_gain, readnoise),fontsize=20)
		fig_name = logs.replace('.log','_pa_'+str(kang))		
		#fig_name = logs.replace('.runlog','_pa_'+str(flata))
		fn = '/home/gillian/data/otherdata/kcwi/figures/'
		fig_dir = os.path.dirname(fn)
		if not os.path.exists(fig_dir):# create data directory if needed
		    os.makedirs(fig_dir)
		plt.savefig(fn + fig_name + '.png')

def plot_meanpol(indata):

	logs,counts,dolp,pa,angle1,angle2,angle3,idxw,allwaves,kang,karr,dolp_mean,pa_mean = indata

        print np.array(dolp_mean).shape
        print np.array(karr).shape
	if plot_flag:
		plt.close("all")
		plt.clf()
		fig = plt.figure(figsize=(20,15))
		plt.rc('text', usetex=True)
		plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})
		plt.xlabel(r"\textbf{K-mirror angle [$^\circ$]}",fontsize=25)
		plt.ylabel(r"\textbf{DoLP [\%]}",fontsize=25)
		plt.grid(b=True, which='major', color='0.75', linestyle='--')
		plt.grid(b=True, which='minor', color='0.75', linestyle='--')
		plt.tick_params(axis='x', labelsize=25)
		plt.tick_params(axis='y', labelsize=25)
		#plt.axis([0,np.max(n_log),0,bins[bins.size-1]])
		
		axes = plt.gca()
		axes.set_ylim([0,10])
		#axes.set_xlim([10**0,10**5])
		plt.plot(karr, dolp_mean, color='darkgreen',marker='.',markersize=8,linestyle='--',linewidth=1,label='Degree of Linear Polarisation: 4600-4800A')
		plt.legend(loc=1,prop={'size':20})
		#figtext(.65, .15,'Pre-amp = %0.4f e-/ADU.\nReadnoise = %0.4f e-' % (preamp_gain, readnoise),fontsize=20)
		fig_name = logs.replace('.log','_meandolp')
		fn = '/home/gillian/data/otherdata/kcwi/figures/'
		fig_dir = os.path.dirname(fn)
		if not os.path.exists(fig_dir):# create data directory if needed
		    os.makedirs(fig_dir)
		plt.savefig(fn + fig_name + '.png')


	if plot_flag:
		plt.close("all")
		plt.clf()
		fig = plt.figure(figsize=(20,15))
		plt.rc('text', usetex=True)
		plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})
		plt.xlabel(r"\textbf{K-mirror angle [$^\circ$]}",fontsize=25)
		plt.ylabel(r"\textbf{PA [$^\circ$]}",fontsize=25)
		plt.grid(b=True, which='major', color='0.75', linestyle='--')
		plt.grid(b=True, which='minor', color='0.75', linestyle='--')
		plt.tick_params(axis='x', labelsize=25)
		plt.tick_params(axis='y', labelsize=25)
		#plt.axis([0,np.max(n_log),0,bins[bins.size-1]])
		
		axes = plt.gca()
		#axes.set_ylim([0,50])
		#axes.set_xlim([10**0,10**5])
		plt.plot(karr, pa_mean, color='navy',marker='.',markersize=8,linestyle='--',linewidth=1,label='Polarisation Angle: 4600-4800A')
		plt.legend(loc=1,prop={'size':20})
		#figtext(.65, .15,'Pre-amp = %0.4f e-/ADU.\nReadnoise = %0.4f e-' % (preamp_gain, readnoise),fontsize=20)
		fig_name = logs.replace('.log','_meanpa')		
		#fig_name = logs.replace('.runlog','_pa_'+str(flata))
		fn = '/home/gillian/data/otherdata/kcwi/figures/'
		fig_dir = os.path.dirname(fn)
		if not os.path.exists(fig_dir):# create data directory if needed
		    os.makedirs(fig_dir)
		plt.savefig(fn + fig_name + '.png')


def plot_standardpol(indata):

	log,counts,dolp,pa,angle1,angle2,angle3,idxw,allwaves,dolp_mean,pa_mean,rotpdest,sky_dolp,sky_pa,skypa,target = indata

	print rotpdest
	print sky_dolp

	if plot_flag:
		plt.close("all")
		plt.clf()
		fig = plt.figure(figsize=(20,15))
		plt.rc('text', usetex=True)
		plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})
		plt.xlabel(r"\textbf{K-mirror angle [$^\circ$]}",fontsize=25)
		plt.ylabel(r"\textbf{DoLP [\%]}",fontsize=25)
		plt.grid(b=True, which='major', color='0.75', linestyle='--')
		plt.grid(b=True, which='minor', color='0.75', linestyle='--')
		plt.tick_params(axis='x', labelsize=25)
		plt.tick_params(axis='y', labelsize=25)
		#plt.axis([0,np.max(n_log),0,bins[bins.size-1]])
		
		axes = plt.gca()
		axes.set_ylim([0,25])
		#axes.set_xlim([10**0,10**5])
		plt.plot(rotpdest, sky_dolp, color='darkgreen',marker='.',markersize=8,linestyle='--',linewidth=1,label='Degree of Linear Polarisation: 4600-4800A')
		plt.legend(loc=1,prop={'size':20})
		#figtext(.65, .15,'Pre-amp = %0.4f e-/ADU.\nReadnoise = %0.4f e-' % (preamp_gain, readnoise),fontsize=20)
		fig_name = data_dir[-7:-1]+str('_')+str(target)+str('_meandolp')
		#fig_name = logs.replace('.log','_meandolp')
		fn = '/home/gillian/data/otherdata/kcwi/figures/'
		fig_dir = os.path.dirname(fn)
		if not os.path.exists(fig_dir):# create data directory if needed
		    os.makedirs(fig_dir)
		plt.savefig(fn + fig_name + '.png')


	if plot_flag:
		plt.close("all")
		plt.clf()
		fig = plt.figure(figsize=(20,15))
		plt.rc('text', usetex=True)
		plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})
		plt.xlabel(r"\textbf{K-mirror angle [$^\circ$]}",fontsize=25)
		plt.ylabel(r"\textbf{PA [$^\circ$]}",fontsize=25)
		plt.grid(b=True, which='major', color='0.75', linestyle='--')
		plt.grid(b=True, which='minor', color='0.75', linestyle='--')
		plt.tick_params(axis='x', labelsize=25)
		plt.tick_params(axis='y', labelsize=25)
		#plt.axis([0,np.max(n_log),0,bins[bins.size-1]])
		
		axes = plt.gca()
		#axes.set_ylim([0,50])
		#axes.set_xlim([10**0,10**5])
		plt.plot(rotpdest, sky_pa, color='navy',marker='.',markersize=8,linestyle='--',linewidth=1,label='Polarisation Angle: 4600-4800A')
		plt.legend(loc=1,prop={'size':20})
		#figtext(.65, .15,'Pre-amp = %0.4f e-/ADU.\nReadnoise = %0.4f e-' % (preamp_gain, readnoise),fontsize=20)
		fig_name = data_dir[-7:-1]+str('_')+str(target)+str('_meanpa')
		#fig_name = logs.replace('.log','_meanpa')		
		#fig_name = logs.replace('.runlog','_pa_'+str(flata))
		fn = '/home/gillian/data/otherdata/kcwi/figures/'
		fig_dir = os.path.dirname(fn)
		if not os.path.exists(fig_dir):# create data directory if needed
		    os.makedirs(fig_dir)
		plt.savefig(fn + fig_name + '.png')

def plot_sine(indata):

	logs,counts,dolp,pa,angle1,angle2,angle3 = indata

	if plot_flag:
		plt.close("all")
		plt.clf()
		fig = plt.figure(figsize=(20,15))
		plt.rc('text', usetex=True)
		plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})
		plt.xlabel(r"\textbf{Polariser [deg]}",fontsize=25)
		plt.ylabel(r"\textbf{Intensity [e]}",fontsize=25)
		plt.grid(b=True, which='major', color='0.75', linestyle='--')
		plt.grid(b=True, which='minor', color='0.75', linestyle='--')
		plt.tick_params(axis='x', labelsize=25)
		plt.tick_params(axis='y', labelsize=25)
		#plt.axis([0,np.max(n_log),0,bins[bins.size-1]])
		
		axes = plt.gca()
		#axes.set_ylim([0,150])
		#axes.set_xlim([10**0,10**5])
		plt.plot(counts[0:13,2000], color='darkgreen',linestyle='-',linewidth=1,label='Sine plot')
		#plt.legend(loc=1,prop={'size':20})
		#figtext(.65, .15,'Pre-amp = %0.4f e-/ADU.\nReadnoise = %0.4f e-' % (preamp_gain, readnoise),fontsize=20)
		fig_name = logs.replace('.log','.sine.2000')
		fn = '/home/gillian/data/otherdata/kcwi/figures/'
		fig_dir = os.path.dirname(fn)
		if not os.path.exists(fig_dir):# create data directory if needed
		    os.makedirs(fig_dir)
		plt.savefig(fn + fig_name + '.png')



## ----- Entry point for script -----

if __name__ == "__main__":

        plot_standardpol(standards())
        #kmirror(loaddata(sys.argv[1]))
        #stokes(loaddata(sys.argv[1]))
	#plot_meanpol(kmirror(loaddata(sys.argv[1])))
	#plot_sine(stokes(loaddata(sys.argv[1])))



