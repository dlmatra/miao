#!/usr/bin/env python
# coding: utf-8

# ## Analyzing and fitting interferometric visibility datasets of circumstellar rings

# Bringing your data from archive to paper-ready plots requires several steps, which we will follow and execute in this Notebook. Built in **Python 3 and with CASA 5.4.0**, but may handle other versions.
# 
# **GALARIO package**:
# https://github.com/mtazzari/galario
# **needs to be installed within local Python 3 installation**, e.g. with conda with a command like: 'conda install -c conda-forge galario'. This code was tested with version 1.2
# 
# **RADMC-3D version 0.41**:
# http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/index.html
# comes with the package, in a version that has been slightly modified to avoid print statements if nothing goes wrong, and hence save time. **This however needs to be installed** through: 'cd radmc-gala/radmc-3d/version_0.41/srcnoprint' and then 'make'. It is strongly advised that the user familiarises themselves with the ray-tracing capabilities of radmc-3d, although this knowledge should not be needed to run radmc-gala.
# The **radmc3dPy** add-on: https://www.ast.cam.ac.uk/~juhasz/radmc3dPyDoc/index.html is used to read radmc-3d files into Python, and is included within the radmc-3d package.
# 
# **Astropy**:
# https://www.astropy.org/
# **needs to be installed within local Python 3 installation**, to enable read and write of FITS format files. Installation is as simple as: 'conda install astropy' or 'pip install astropy'.

# ## Step 0: Create list with all requested inputs

# Here you need to decide if you are doing imaging and/or fitting and/or 
# input some parameters that are needed for the imaging and/or fitting and/or postprocessing.
# 

# In[1]:


import numpy as np
import os
import pickle

#Decide which you want to do, imaging and/or fitting and/or postprocessing
imaging=True
fit=True
postproc=True

radmcgalapath='/d1/boudica1/lmatra/radmc-gala'
casapath='/d1/boudica1/casa-release-5.4.0-68.el6/bin'
sourcetag='GJ14'
workingdir='/d1/boudica1/lmatra'
vis=['/sma/SMAusers/lmatra/REASONS_ALMA/cycle5/GJ14/calibratedms/GJ14_calibratedvis_cont.ms']
nvis=len(vis)
if not os.path.exists(workingdir+'/'+sourcetag):
    os.mkdir(workingdir+'/'+sourcetag)
    print('Creating directory for object '+sourcetag+' at '+workingdir)
else:
    print('Directory for object '+sourcetag+' already exists at '+workingdir)

if imaging:
    #Imaging parameters
    mosaic=False
    if mosaic:
        mosaic=True
        phasecenter='J2000 22h57m39.449700 -29d37m22.68681'
    else:
        phasecenter=''
    weighting='natural'
    if weighting=='briggs':
        robust=0.5
    else:
        robust=''
    uvtaper=['']
    interactive=True

if fit:
    #Parameters for fit
    Lstar=0.111  #Solar luminosities
    dist=14.7  #pc
    imsize=24.0     #arcsec, used for radiative transfer. >>2x belt outer radius, but not too large or it will slow down computation.
    #imsize above is also the size of the grid over which the model is setup in RADMC. No disk will be put beyond this.

    #Add disk parameters
    fluxdensity=1.8e-3 #Jy
    fluxdensity_dwn, fluxdensity_up = [0.1e-3, 50e-3]
    rmid=99.0/dist # Radial peak location of Gaussian surface density, arcsec
    rmid_dwn, rmid_up = [0.3, 12.0]
    sigma=32.0/2.35/dist #standard deviation of radial Gaussian surface density
    sigma_dwn, sigma_up = [0.05, 'halfrmid'] #halfrmid is special setup that enables maximum sigma to be half whatever rmid
    #is in current model evaluation.
    useh=True
    if useh:
        h=0.05 #aspect ratio of belt, constant with radius, whose vertical density structure at radius r is a Gaussian with width hr.
        h_dwn, h_up = [0.005, 0.2]
    else:
        h=0.03
    incl=65.0 #inclination, degrees from face-on
    incl_dwn, incl_up = [0.1, 89.99]
    posang=5.5 #position angle, East of North. **Needs to be between 0 and 180 degrees!**
    if posang<45.0 or posang>135.0:
        posang_dwn, posang_up = [-90.0, 90.0]
    else:
        posang_dwn, posang_up = [0.0, 180.0]
       
    
    #Add star if wanted/needed
    star=True
    if star:
        fstar=4e-5 #Jy
        fstar_dwn, fstar_up = [1e-6, 2e-3]
    
    #Add extra parameters independently for each of the visibility datasets.
    if nvis>1:
        dRA=[-0.07,0.1,0.1] #RA offset of star+disk from phase center of observations
        dRA_dwn=[-3.0,-3.0,-3.0]
        dRA_up=[3.0,3.0,3.0]
        dDec=[0.14,-0.1,-0.1] #Dec offset of star+disk from phase center of observations
        dDec_dwn=[-3.0,-3.0,-3.0]
        dDec_up=[3.0,3.0,3.0]
        wtfact=[0.289,0.289,1e-4] #factors by which incorrect weights should be multiplied by\
        wtfact_dwn=[1e-6,1e-6,1e-6]
        wtfact_up=[10.0,10.0,10.0]
    else:
        dRA=[-0.07] #RA offset of star+disk from phase center of observations
        dRA_dwn=[-3.0]
        dRA_up=[3.0]
        dDec=[0.14] #Dec offset of star+disk from phase center of observations
        dDec_dwn=[-3.0]
        dDec_up=[3.0]
        wtfact=[0.289] #factors by which incorrect weights should be multiplied by
        wtfact_dwn=[1e-6]
        wtfact_up=[10.0]

    #Add galaxies if needed
    ngal=0
    if ngal>=1:
        resolved=[True, False, False] # if definitely resolved, use 2D Gaussian as galaxy model (6 free parameters), and set resolved=True for that galaxy. 
    #Otherwise, use point source (3 free parameters) by setting resolved=False.
        fbkg=[285e-6,500e-6,500e-6] #Flux (Jy)
        fbkg_dwn=[1e-6,1e-6,1e-6]
        fbkg_up=[10e-3,10e-3,10e-3]
        dRAbkg=[2.97,2.97,2.97] # RA offset (")
        dRAbkg_dwn=[x-3.0 for x in dRAbkg]
        dRAbkg_up=[x+3.0 for x in dRAbkg]
        dDecbkg=[1.81,1.81,1.81] # Dec offset (")
        dDecbkg_dwn=[x-3.0 for x in dDecbkg]
        dDecbkg_up=[x+3.0 for x in dDecbkg]
        sigmagal=[0.3,0,0] # sigma (")
        sigmagal_dwn=[0.05,0.05,0.05]
        sigmagal_up=[1.0,1.0,1.0]
        PAgal=[28.0,0,0] # PA (deg), East of North. **Needs to be between 0 and 180 degrees!**
        PAgal_dwn=[]
        PAgal_up=[]
        for i in PAgal:
            if PAgal<45.0 or PAgal>135.0:
                PAgal_dwn.append(-90.0)
                PAgal_up.append(90.0)
            else:
                PAgal_dwn.append(0.0)
                PAgal_up.append(180.0)   
        incgal=[49.0,0,0] # inc (deg)
        incgal_dwn=[0.1,0.1,0.1]
        incgal_up=[89.99,89.99,89.99]
    else: 
        resolved=None
        fbkg=None
        dRAbkg=None
        dDecbkg=None
        sigmagal=None
        PAgal=None
        incgal=None
        
    #Now define array containing all initial parameters
    ## NB DO NOT MODIFY ORDER AS SETUP_MCMC CODE WILL NOT RECOGNISE PARAMETERS CORRECTLY
    pars_init=[fluxdensity, rmid, sigma, incl, posang]
    priors_dwn=[fluxdensity_dwn, rmid_dwn, sigma_dwn, incl_dwn, posang_dwn]
    priors_up=[fluxdensity_up, rmid_up, sigma_up, incl_up, posang_up]
    if useh:
        pars_init.append(h)
        priors_dwn.append(h_dwn)
        priors_up.append(h_up)
    if star:
        pars_init.append(fstar)
        priors_dwn.append(fstar_dwn)
        priors_up.append(fstar_up)
    for i in np.arange(nvis):
        pars_init.append(dRA[i])
        priors_dwn.append(dRA_dwn[i])
        priors_up.append(dRA_up[i])
        pars_init.append(dDec[i])
        priors_dwn.append(dDec_dwn[i])
        priors_up.append(dDec_up[i])
        pars_init.append(wtfact[i])
        priors_dwn.append(wtfact_dwn[i])
        priors_up.append(wtfact_up[i])
    if ngal>=1:
        for i in np.arange(ngal):
            pars_init.append(fbkg[i])
            priors_dwn.append(fbkg_dwn[i])
            priors_up.append(fbkg_up[i])
            pars_init.append(dRAbkg[i])
            priors_dwn.append(dRAbkg_dwn[i])
            priors_up.append(dRAbkg_up[i])
            pars_init.append(dDecbkg[i])
            priors_dwn.append(dDecbkg_dwn[i])
            priors_up.append(dDecbkg_up[i])
            if resolved[i]:
                pars_init.append(sigmagal[i])
                priors_dwn.append(sigmagal_dwn[i])
                priors_up.append(sigmagal_up[i])
                pars_init.append(PAgal[i])
                priors_dwn.append(PAgal_dwn[i])
                priors_up.append(PAgal_up[i])
                pars_init.append(incgal[i])
                priors_dwn.append(incgal_dwn[i])
                priors_up.append(incgal_up[i])


# ## Step 1: Create directory structure

# In[2]:


print('Creating directory structure')
os.chdir(workingdir+'/'+sourcetag)

#Imaging
if imaging:
    print('Will be carrying out imaging')
    for i in ['calibratedms', 'imaging']:
        if not os.path.exists(workingdir+'/'+sourcetag+'/'+i):
            os.mkdir(workingdir+'/'+sourcetag+'/'+i)
    get_ipython().system('cp -r {radmcgalapath}/utils/mstonumpyortxt_multiple.py {workingdir}/{sourcetag}/calibratedms/.')
    for i in np.arange(len(vis)):
        if not os.path.exists('calibratedms/'+vis[i].rsplit('/',1)[1]):
            get_ipython().system('cp -r {vis[i]} {workingdir}/{sourcetag}/calibratedms/.')
        vis[i]=vis[i].rsplit('/',1)[1]
        #print(vis[i])
    get_ipython().system('cp -r {radmcgalapath}/utils/imagingscript_multiple.py {workingdir}/{sourcetag}/imaging/.  ')
    #Save imaging parameters
    pickle.dump([sourcetag,workingdir,vis,nvis,mosaic,phasecenter,weighting,robust,uvtaper,interactive], open(workingdir+'/'+sourcetag+'/imaging/imagepars.npy', 'wb'), protocol=2)

    
#Fitting
if fit: 
    print('Will be carrying out visibility fit')
    if not os.path.exists(workingdir+'/'+sourcetag+'/'+'imaging'):
        sys.exit('Carry out imaging first: need .pb primary beam image')
    for i in ['uvfit']:
        if not os.path.exists(workingdir+'/'+sourcetag+'/'+i):
            os.mkdir(workingdir+'/'+sourcetag+'/'+i)
    get_ipython().system('cp -r {radmcgalapath}/utils/setup_mcmc.py {workingdir}/{sourcetag}/uvfit/. ')
    get_ipython().system('cp -r {radmcgalapath}/utils/problem_setup_cont_gauss.py {workingdir}/{sourcetag}/uvfit/.  ')
    get_ipython().system('cp -r {radmcgalapath}/utils/dustkappa_10445.micr.inp {workingdir}/{sourcetag}/uvfit/. ')
    #Predict names of primary beam files according to standard naming convention
    pbfilenames=[[] for x in vis]
    for i in np.arange(nvis):
        if not imaging:
            vis[i]=vis[i].rsplit('/',1)[1]
        pbfilenames[i]=vis[i][:-3]+'_'+weighting+robust+'_pb.fits'
    #Save fit parameters
    pickle.dump([pbfilenames,vis,nvis,radmcgalapath,workingdir,sourcetag,Lstar,dist,imsize,useh,star
                 ,ngal,resolved, pars_init, priors_dwn, priors_up], open(workingdir+'/'+sourcetag+'/uvfit/fitpars.npy', 'wb'), protocol=2)

#Postprocessing
if postproc:
    if not (imaging or fit):
        print ('Will be carrying out postprocessing ONLY')
    else:
        print('Will be carrying out postprocessing')
    if not os.path.exists(workingdir+'/'+sourcetag+'/'+'uvfit'):
        sys.exit('Carry out fit before POSTprocessing!')
    for i in ['analysis', 'plots', 'uvfit/evaluation']:
        if not os.path.exists(workingdir+'/'+sourcetag+'/'+i):
            os.mkdir(workingdir+'/'+sourcetag+'/'+i)
    get_ipython().system('cp -r {radmcgalapath}/utils/evaluatemodel_radmc3d.py {workingdir}/{sourcetag}/uvfit/. ')
    get_ipython().system('cp -r {radmcgalapath}/utils/uvresidualtoms.py {workingdir}/{sourcetag}/uvfit/. ')
    get_ipython().system('cp -r {radmcgalapath}/utils/makeuvdeprojplot_simple_multiple.py {workingdir}/{sourcetag}/uvfit/.')
    get_ipython().system('cp -r {radmcgalapath}/utils/plotimage.py {workingdir}/{sourcetag}/analysis/.  ')
    get_ipython().system('cp -r {radmcgalapath}/utils/imagecombo.py {workingdir}/{sourcetag}/analysis/. ')
   


# ## Step 2: Carry out imaging via imagingscript_multiple CASA script
# First convert visibilities in CASA MS format to a python save file

# In[3]:


os.chdir('calibratedms')
get_ipython().system('{casapath}/casa -c mstonumpyortxt_multiple.py')


# Read in visibilities into Python, and use Galario to figure out ideal cell size and image size.
# Typically half the suggested image size is OK.

# In[4]:


from galario.double import get_image_size
u=[[] for x in vis]
v=[[] for x in vis]
Re=[[] for x in vis]
Im=[[] for x in vis]
w=[[] for x in vis]
nxy=[[] for x in vis]
dxy=[[] for x in vis]
for i in np.arange(nvis):
    u[i], v[i], Re[i], Im[i], w[i] = np.load(vis[i][:-3]+'.npy')
    nxy[i], dxy[i] = get_image_size(u, v)
    nxy[i]/=2
    print('Pixel size (arcsec) and number of pixels required for dataset '+vis[i]+':')
    print(dxy[i]*180.0/np.pi*3600, nxy[i])    

#Figure out pix size and number for concatenated image.
dxyall=np.min(dxy)
nxyall=np.int(np.ceil(np.max(np.asarray(dxy)*np.asarray(nxy))/dxyall/2.0)*2)
print('Pixel size (arcsec) and number of pixels that will be used:')
print(dxyall*180.0/np.pi*3600, nxyall)    


    
#Save pixel sizes and number of pixels, for imaging and fitting
pickle.dump([dxyall, np.int(nxyall)], open(workingdir+'/'+sourcetag+'/calibratedms/pixinfo.npy', 'wb'), protocol=2)


# Then run the CLEANing using CASA's tclean.
# 
# Run this within CASA locally. There is no way around this if you want to use the interactive cleaning mode!

# In[5]:


os.chdir('../imaging')
#THE COMMAND BELOW to be run on local computer for CASA to bring up interactive prompt
execfile('imagingscript_multiple.py')


# ## Step 3: Run visibility fit by forward-modeling with RADMC-3D (ray tracing) and GALARIO (for FFT) through emcee package

# In[6]:


os.chdir('../uvfit')
backendaddress='backend_'+sourcetag+'_todaysdate_computersname.pkl'
#np.save('backendaddress.npy', backendaddress)
pickle.dump(backendaddress, open('backendaddress.npy', 'wb'), protocol=2)


#print(np.load('backendaddress.npy'))
#!ls

#import setup_mcmc
get_ipython().run_line_magic('run', '-i setup_mcmc')
nsteps=10

#print(priors_up)
# This step actually starts the MCMC, at which point a progress bar should come up (may need installation of a python package).
# The MCMC runs for nsteps and starts the walkers for each parameter at position 'pos' defined above
sampler.run_mcmc(pos,nsteps,progress=True)
# If number of steps is insufficient, i.e. if the chains have not converged for each parameter, run for more steps using:
#sampler.run_mcmc(None,nsteps,progress=True)
# where nsteps is the number of steps you want to go further by.


# In[ ]:





# In[ ]:





# In[ ]:




