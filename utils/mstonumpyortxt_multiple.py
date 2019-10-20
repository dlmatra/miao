#!/usr/bin/env python
# coding: utf-8

# ## Convert visibilities from CASA MS format to Python save file
# To be executed in CASA.

# In[ ]:


import numpy as np
import pickle
import os
sourcetag,workingdir,vis,nvis,mosaic,phasecenter,weighting,robust,uvtaper,interactive = pickle.load(open('../imaging/imagepars.npy','rb'))
vis=[x.encode('ascii') for x in vis]

tb = casac.table()
ms = casac.ms()
cc=2.9979e10 #cm/s
    
#Use CASA table tools to get columns of UVW, DATA, WEIGHT, etc.
outputfilename=[x[:-3]+'.npy' for x in vis]
tb = casac.table()
ms = casac.ms()
cc=2.9979e10 #cm/s
    
#Use CASA table tools to get columns of UVW, DATA, WEIGHT, etc.
filename=vis

for ii in np.arange(len(filename)):
	tb.open(filename[ii])
	#NB for this to run smooth we need to have the same number of channels for ALL scans. So no spectral windows with different number of channels, otherwise it gets complicated.
	#See https://safe.nrao.edu/wiki/pub/Main/RadioTutorial/BandwidthSmearing.pdf to choose how much to smooth a dataset in frequency.
	#Errors are taking into account when time averaging in split: https://casa.nrao.edu/casadocs/casa-5.1.1/uv-manipulation/time-average
	#And when channel averaging: https://casa.nrao.edu/casadocs/casa-5.1.1/uv-manipulation/channel-average
	data    = tb.getcol("DATA")
	uvw     = tb.getcol("UVW")
	flags   = tb.getcol("FLAG")
	spwid   = tb.getcol("DATA_DESC_ID")
	weight  = tb.getcol("WEIGHT")
	ant1    = tb.getcol("ANTENNA1")
	ant2    = tb.getcol("ANTENNA2")
	tb.close()
	if np.any(flags):
		print "Note: some of the data is FLAGGED"
	print "Found data with "+str(data.shape[-1])+" uv points per channel"


	#Use CASA ms tools to get the channel/spw info
	ms.open(filename[ii])
	spwstuff = ms.getspectralwindowinfo()
	nchan = spwstuff["0"]["NumChan"]
	npol = spwstuff["0"]["NumCorr"]
	ms.close()
	print "with "+str(nchan)+" channels per SPW and "+str(npol)+" polarizations,"

	#Use CASA table tools to get frequencies, which are needed to calculate u-v points from baseline lengths
	tb.open(filename[ii]+"/SPECTRAL_WINDOW")
	freqs = tb.getcol("CHAN_FREQ")
	rfreq = tb.getcol("REF_FREQUENCY")
	tb.close()
	print str(freqs.shape[1])+" SPWs and Channel 0 frequency of 1st SPW of "+str(rfreq[0]/1e9)+" GHz"
	print "corresponding to "+str(2.9979e8/rfreq[0]*1e3)+" mm"

	print "Datasets has baselines between "+str(np.min(np.sqrt(uvw[0,:]**2.0+uvw[1,:]**2.0)))+" and "+str(np.max(np.sqrt(uvw[0,:]**2.0+uvw[1,:]**2.0)))+" m"  

	#Initialize u and v arrays (coordinates in Fourier space)
	uu=np.zeros((freqs.shape[0],uvw[0,:].size))
	vv=np.zeros((freqs.shape[0],uvw[0,:].size))

	#Fill u and v arrays appropriately from data values.
	for i in np.arange(freqs.shape[0]):        
		for j in np.arange(uvw.shape[1]):
			uu[i,j]=uvw[0,j]*freqs[i,spwid[j]]/(cc/100.0)
			vv[i,j]=uvw[1,j]*freqs[i,spwid[j]]/(cc/100.0)

	#Extract real and imaginary part of the visibilities at all u-v coordinates, for both polarization states (XX and YY), extract weights which correspond to 1/(uncertainty)^2
	Re_xx = data[0,:,:].real
	Im_xx = data[0,:,:].imag
	weight_xx = weight[0,:]
	if npol>=2:
		Re_yy = data[1,:,:].real
		Im_yy = data[1,:,:].imag
		weight_yy = weight[1,:]
		#print('Ciao')

		#Since we don't care about polarization, combine polarization states (average them together) and fix the weights accordingly. Also if any of the two polarization states is flagged, flag the outcome of the combination.
		flags = flags[0,:,:]*flags[1,:,:]
		Re = np.where((weight_xx + weight_yy) != 0, (Re_xx*weight_xx + Re_yy*weight_yy) / (weight_xx + weight_yy), 0.)
		Im = np.where((weight_xx + weight_yy) != 0, (Im_xx*weight_xx + Im_yy*weight_yy) / (weight_xx + weight_yy), 0.)
		wgts = (weight_xx + weight_yy)
	else:
		Re=Re_xx
		Im=Im_xx
		wgts=weight_xx
		flags=flags[0,:,:]

	# Find which of the data represents cross-correlation between two antennas as opposed to auto-correlation of a single antenna.
	# We don't care about the latter so we don't want it.
	xc = np.where(ant1 != ant2)[0]
	
	#Select only cross-correlation data
	data_real = Re[:,xc]
	data_imag = Im[:,xc]
	flags = flags[:,xc]
	data_wgts = wgts[xc]
	data_uu = uu[:,xc]
	data_vv = vv[:,xc]
	data_wgts=np.reshape(np.repeat(wgts[xc], uu.shape[0]), data_uu.shape)
	#Delete previously used (and not needed) variables (to free up some memory?)
	del Re
	del Im
	del wgts
	del uu
	del vv

	#Select only data that is NOT flagged
	data_real = data_real[np.logical_not(flags)]
	data_imag = data_imag[np.logical_not(flags)]
	flagss = flags[np.logical_not(flags)]
	data_wgts = data_wgts[np.logical_not(flags)]
	data_uu = data_uu[np.logical_not(flags)]
	data_vv = data_vv[np.logical_not(flags)]

	#Wrap up all the arrays/matrices we need, (u-v coordinates, complex visibilities, and weights for each visibility) and save them all together in a numpy file
	u, v, Re, Im, w = data_uu, data_vv, data_real, data_imag, data_wgts
	#print(filename[ii][:-3]+'.npy')
	np.save(filename[ii][:-3]+'.npy', [u, v, Re, Im, w])

