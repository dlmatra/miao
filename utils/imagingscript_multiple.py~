#!/usr/bin/env python
# coding: utf-8

# ## Carry out imaging in CASA 
# given parameters in imagepars.npy file, which were fed by the main.py routine
# 
# Needs to be run within imaging folder in the presence of imagepars.npy

# In[2]:


#!/usr/bin/env python
# coding: utf-8

# ## Carry out imaging in CASA 
# given parameters in imagepars.npy file, which were fed by the main.py routine
# 
# Needs to be run within imaging folder in the presence of imagepars.npy

# In[ ]:


import numpy as np
import pickle
import os

#Read in imaging parameters
sourcetag,workingdir,vis,nvis,mosaic,phasecenter,weighting,robust,uvtaper,interactive = pickle.load(open('./imagepars.npy','rb'))
sourcetag=sourcetag.encode('ascii')
workingdir=workingdir.encode('ascii')
vis=[x.encode('ascii') for x in vis]
phasecenter=phasecenter.encode('ascii')
weighting=weighting.encode('ascii')
#robust=robust.encode('ascii')
uvtaper=[x.encode('ascii') for x in uvtaper]

#Read in pixel parameters
dxy, nxy = pickle.load(open(workingdir+'/'+sourcetag+'/calibratedms/pixinfo.npy','rb'))

#print(workingdir)
#Run this within imaging folder
concatvis=sourcetag+'_calibratedvis_cont_concat.ms'
if any('_res' in x for x in vis):
    concatvis=sourcetag+'_calibratedvis_cont_concat_res.ms'
if any('_model' in x for x in vis):   
    concatvis=sourcetag+'_calibratedvis_cont_concat_model.ms'
if any('_rwdat' in x for x in vis):   
    concatvis=sourcetag+'_calibratedvis_cont_concat_rwdat.ms'
if nvis>1:
    concaten=True
    weightfacts=[1.0 for x in np.arange(len(vis))]
    imageconcat=True
    imagesingles=True
else:
    concaten=False
    imagesingles=True
    imageconcat=False

if concaten:
    os.system('rm -r '+workingdir+'/'+sourcetag+'/calibratedms/'+concatvis)
    concat(vis=[workingdir+'/'+sourcetag+'/'+'calibratedms/'+x for x in vis], concatvis=workingdir+'/'+sourcetag+'/calibratedms/'+concatvis, visweightscale=weightfacts, copypointing=False)
    #If needed, manually fix coordinates to match coordinates of first observation, otherwise mosaic won't work. Phase centers are all aligned already from fit
    #fixplanets(vis=concatvis, field='0,1', direction='J2000 02h26m16.337489 +06d17m32.38948')
    listobs(workingdir+'/'+sourcetag+'/calibratedms/'+concatvis)

if imageconcat:
    imagename=concatvis[:-3]+'_'+weighting+robust
    #clean parameters
    imsize=[nxy, nxy]
    cell=[str(dxy*180.0/np.pi*3600.0)+'arcsec']
    pblimit=1e-5
    if mosaic:
        gridder='mosaic'
    else:
        gridder='standard'
    deconvolver='multiscale'
    #Scales should be roughly [0, n where n*cell~expected syntesized beam size, 3n, 9n, etc.]
    scales=[0,72,216,648]
    niter=1
    specmode='mfs'
    if robust=='':
        robusttask=0.5
    
    #Remove image if it exists
    os.system('rm -r '+imagename+'.*')
    
    #Run iterative tclean with manual masking
    tclean(vis=workingdir+'/'+sourcetag+'/calibratedms/'+concatvis, interactive=interactive, imsize=imsize, cell=cell, weighting=weighting, niter=niter, specmode=specmode, gridder=gridder, deconvolver=deconvolver, scales=scales, imagename=imagename, uvtaper=uvtaper,  robust=robusttask, pblimit=pblimit)
    
    #Export image to FITS
    exportfits(imagename=imagename+'.image', fitsimage=imagename+'.fits', overwrite=True)
    #Export primary beam to FITS
    exportfits(imagename=imagename+'.pb', fitsimage=imagename+'_pb.fits', overwrite=True)
    
    #View result
    viewer(imagename+'.image')

if imagesingles:
    for i in np.arange(nvis):
        imagename=vis[i][:-3]+'_'+weighting+robust
        #print(imagename)
        imsize=[nxy,nxy]
        cell=[str(dxy*180.0/np.pi*3600.0)+'arcsec']
        pblimit=1e-5
        gridder='standard'
        deconvolver='multiscale'
        #Scales should be roughly [0, n where n*cell~expected syntesized beam size, 3n, 9n, etc.]
        scales=[0,10,30,90]
        niter=1
        specmode='mfs'
        if robust=='':
            robusttask=0.5

        #Remove image if it exists
        os.system('rm -r '+imagename+'.*')

        #Run iterative tclean with manual masking
        #print(workingdir+'/'+sourcetag+'/'+'calibratedms/'+vis[i])
        tclean(vis=workingdir+'/'+sourcetag+'/'+'calibratedms/'+vis[i], interactive=interactive, imsize=imsize, cell=cell,
               weighting=str(weighting), niter=niter, specmode=specmode, gridder=gridder, deconvolver=deconvolver, scales=scales, 
               imagename=str(imagename), uvtaper=uvtaper,  robust=robusttask)

        #Export image to FITS
        exportfits(imagename=imagename+'.image', fitsimage=imagename+'.fits', overwrite=True)
        #Export primary beam to FITS
        exportfits(imagename=imagename+'.pb', fitsimage=imagename+'_pb.fits', overwrite=True)

        #View result
        viewer(imagename+'.image')

        #RMS X uJy for X" taper, beam X" x X" @Xdeg PA


# In[ ]:





# In[ ]:





