#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def get_sampler(backendaddress):
    #Import general python packages
    import numpy as np
    import matplotlib.pyplot as pl
    import os
    os.environ["OMP_NUM_THREADS"] = "1"
    import astropy.io.fits as pf
    import uuid
    from galario.double import chi2Image, sampleImage, threads, get_image_size, sweep
    threads(1)

    c = 299792458.0

    #Import fit parameters printed within main.py
    pbfilenames,vis,nvis,radmcgalapath,workingdir,sourcetag,Lstar,dist,imsize,useh,star,
                 ngal(",resolved,", "pars_init,", "priors_dwn,", "priors_up", "=", "pickle.load(open('./fitpars.npy','rb'))")
    rmid_init=pars_init[1]
    sigma_init=pars_init[2]
        

    #Import pixel and image sizes needed by galario
    dxy, nxy = pickle.load(open('../calibratedms/pixinfo.npy','rb')) #NB dxy in radians
    dxyarcsec=dxy*180.0/np.pi*3600.0

    #Import necessary scripts and libraries within radmc-gala
    os.sys.path.append(radmcgalapath+'/utils')
    import problem_setup_cont_gauss
    os.sys.path.append(radmcgalapath+'/src')
    import emcee3
    os.sys.path.append(radmcgalapath+'/radmc-3d/version_0.41/python')
    import radmc3dPy

    #### SYSTEM-SPECIFIC INPUTS ####
    tag=sourcetag
    npix=np.int(np.ceil(imsize/(dxyarcsec)/2.0)*2) #Needs to be EVEN. Pixels you actually want to compute radiative transfer over, depends on size of disk.
    imszarcsec=imsize
    sizeau=imszarcsec*dist

    #### GRID INPUT TO BE THE SAME AS IN PROBLEM_SETUP CODE (used for temperature grid) ####
    nr=8.0*sigma_init/(np.min(dxyarcsec)) #Make sure to cover entire width of disk at radial step equal to smallest pixel size of any dataset.
    AU  = 1.49598e13     # Astronomical Unit       [cm]
    rmin=np.max([1.0,rmid_init-4.*sigma_init])*dist*AU
    rmax=np.min([np.min(dxyarcsec*nxy)/2.0,rmid_init+4.*sigma_init])*dist*AU
    tcsize=20
    ri=(np.arange(nr+1, dtype=float))/nr*(rmax-rmin)+rmin
    dr=(np.subtract(ri[1:],ri[:-1]))
    rc=np.asarray(ri[:-1]+dr/2e0)
    rr=np.reshape(np.repeat(rc, tcsize), (rc.size, tcsize))

    ### THIS GOES OUTSIDE MODEL RUN
    # Write dummy wavelength file
    lambda1 = 0.1e0
    lambda2 = 7.0e0
    lambda3 = 25.e0
    lambda4 = 1.0e4
    n12     = 20
    n23     = 100
    n34     = 30
    lam12   = lambda1 * (lambda2/lambda1)**(np.arange(n12, dtype=float)/(1.e0*n12))
    lam23   = lambda2 * (lambda3/lambda2)**(np.arange(n23, dtype=float)/(1.e0*n23))
    lam34   = lambda3 * (lambda4/lambda3)**(np.arange(n34, dtype=float)/(1.e0*(n34-1.e0)))
    lambd = np.concatenate((lam12,lam23,lam34))
    nlam    = len(lambd)
    wfile = open('wavelength_micron.inp', 'w')
    wfile.write('%d\n'%nlam)
    for ilam in range(nlam): wfile.write('%.9e\n'%lambd[ilam])
    wfile.close()

    # Write the radmc3d.inp control file
    wfile = open('./radmc3d.inp', 'w')
    #wfile.write('%s %d\n'%('nphot'+' =',nphot))
    wfile.write('%s %s\n'%('scattering_mode_max'+' =','0'))
    wfile.write('%s %s\n'%('iranfreqmode '+' =','1'))
    wfile.write('%s %s\n'%('istar_sphere '+' = ','0' ))
    wfile.close()

    # Dust opacity control file, check name of dustkappa_*** file is the appropriate one.
    wfile = open('./dustopac.inp', 'w')
    wfile.write('%-15s %s\n'%('2', 'Format number of this file'))
    wfile.write('%-15s %s\n'%('1', 'Nr of dust species'))
    wfile.write('%s\n'%'============================================================================')
    wfile.write('%-15s %s\n'%('1', 'Way in which this dust species is read'))
    wfile.write('%-15s %s\n'%('0', '0=Thermal grain, 1=Quantum heated'))
    wfile.write('%s %s %s\n'%('10445.micr', '    ', 'Extension of name of dustkappa_***.inp file'))
    wfile.write('%s\n'%'----------------------------------------------------------------------------')
    wfile.close()

    #Make and print temperature grid, assumed blackbody, **This GRID needs to be THE SAME AS IN PROBLEM_SETUP FILE**
    t=278.3*(Lstar**0.25)/np.sqrt(rr/AU)
    #print t[:,0,0]
    #print rr[:,0,0]/AU
    wfile=open('./dust_temperature.dat', 'w')
    wfile.write('%d\n'%1)    # Format number
    wfile.write('%d\n'%(rc.size*tcsize))    # Nr of cells
    wfile.write('%d\n'%1)   # Nr of dust species
    for ith in range(tcsize):
        for ir in range(rc.size):
            wfile.write('%.9e\n'%(t[ir,ith]))
    wfile.close()
    

    #Read file containing visibility data, extracted in CASA through 'mstonumpyortxt.py' script, and read primary
    #beam images
    u=[[] for x in vis]
    v=[[] for x in vis]
    Re=[[] for x in vis]
    Im=[[] for x in vis]
    w=[[] for x in vis]
    pbpad=[[] for x in vis]
    for i in np.arange(nvis):
        u[i], v[i], Re[i], Im[i], w[i] = np.load(vis[i][:-3]+'.npy')
        #Import primary beam (pb) to take its effect into account in the model fitting. 
        pbcov, header_pbcov = pf.getdata(pbfilenames[i],0, header=True)
        #CASA images come as 4D cubes (polarization, frequency, spatial y (north-south) direction, spatial x (east-west) 
        #direction).
        #In the case of dust continuum images we have no polarization or frequency, so removing those dimensions below
        pbcov=pbcov[0,0,:,:]
        #Randomly a few pixels from the PB computed by CASA have nans. Remove them
        pbcov[np.where(np.isnan(pbcov))]=0.0
        #Get wavelength which will be a dummy for radmc3d run
        wav=1e6*c/(header_pbcov['CRVAL3'])
        #Pad PB coverage if needed
        pbpad[i]=pbcov #Should be of shape [nxy[i],nxy[i]] as pre-defined internally

    if ngal>=1:
        if np.any(resolved):
    
            # Radial grid parameters for galaxy, grid over which the dust density will be computed. Make sure the step is smaller or equal to the pixel size, and that we roughly cover the extent of the image in arcsec
            Rmin = 0.00001  # arcsec
            dR = 0.01    # arcsec
            nR = 256
            print('Grid goes out to '+str(dR*nR)+' arcsec radius')
            print('Image goes out to '+str(nxy*dxy*180.0/np.pi*3600.0/2.0)+' arcsec radius')
            
            Rminrad=np.radians(Rmin)
            dRrad=np.radians(dR)
            
            #Set up Gaussian radial profile function for bkg galaxy 
            def GaussianProfileGalaxy(f0, sigma, Rmin, dR, nR):
            sigma *= arcsec
            Rmin *= arcsec
            dR *= arcsec
            R = np.linspace(Rmin, Rmin + dR*nR, nR, endpoint=False)
            return f0 * np.exp(-0.5*(R/sigma)**2.0)

    
    #Function that defines likelihood function and multiplies it by the function computing the prior (defined further below).
    #Takes as input the radial grid defined above, and the parameters p that define the model and that will be varied in the fit.
    #Prints ln(L)+ln(prior), where L is the likelihood function, in this case just the chi squared of the data.
    #If, when calling it, locfiles is specified as an address on the machine, the function prints a numpy save file with the reweighting 
    #factor, one with the uv table containing model visibilities, and one containing geometry of the belt (inclination, position angle) 
    #plus RA and Dec offset from the image center. It also returns a model image and lnL.
    def lnpostfn(p, locfiles=None):

        #Calculate prior.
        lnprior = lnpriorfn(p)
        if not np.isfinite(lnprior):
            return -np.inf
    
        countpars=0
        if useh:
            h=p[5]
            countpars+=1
        else:
            h=0.03
        
        if star:
            fstar=p[5+countpars]
            countpars+=1

        ciao = str(uuid.uuid4())
        os.system('mkdir '+'./'+ciao)
        os.chdir('./'+ciao)
        os.system('cp '+'../dustopac.inp .')
        os.system('cp '+'../wavelength_micron.inp .')
        os.system('cp '+'../radmc3d.inp .')
        os.system('cp '+'../dust_temperature.dat .')
        os.system('cp '+'../dustkappa_10445.micr.inp .')
    
        mdisk=1e-7 #Mearth, dummy just to make sure it's optically thin, which it is anyway.
        problem_setup_cont_gauss.problem_setup([mdisk, p[1]*dist, p[2]*dist, h, rmin, rmax, nr, nth])

        ##############
        ## RUN RADMC AND READ OUTPUT
        ##############
        
        # Run
        os.system(radmcgalapath+'/radmc-3d/version_0.41/srcnoprint/radmc3d image lambda '+str(wav)+' incl '+str(p[3])+
                  ' posang '+str(90.0+p[4])+' sizeau '+str(sizeau)+' npix '+str(npix)+' imageunform nostar')
        # Read
        imag     =     radmc3dPy.image.readImage(binary=True)

        #Normalize model using total flux density which is free parameter
        imagjypixdistscaled = imag.imageJyppix[:,:,0]*p[0]/np.sum(imag.imageJyppix[:,:,0])

        #Pad model so that image is large enough to sample shortest visibility spacings
        modpad=np.zeros((nxy,nxy))
        modpad[(nxy/2-npix/2):(nxy/2+npix/2),(nxy/2-npix/2):(nxy/2+npix/2)]+=imagjypixdistscaled


        # Here compute model complex visibilities stored in complex array vis by Fourier Transforming the model image multiplied 
        # by the primary beam (simulating the response of the antennas)	and evaluating this FT at the u-v points sampled 
        # by the interferometer. Then, rotate the visibilities (equivalent to rotating the model image) according to the position angle PA.
        # Then,  phase-shift the visibilities (equivalent to shifting the star from the image center), according to dRA and dDec.
        # Finally, return the chisquared of the model visibilities as:
        # [(Real_data-Real_model)^2+(Imaginary_data-Imaginary_model)^2]*weight summed over all u-v points.
        vismodel=[[] for x in vis]
        warr=[[] for x in vis]
        dRArad=[[] for x in vis]
        dDecrad=[[] for x in vis]
        for i in np.arange(nvis):
            dRArad[i], dDecrad[i] = p[5+countpars]/3600.0*np.pi/180.0, p[6+countpars]/3600.0*np.pi/180.0
            vismodel[i] = sampleImage(np.ascontiguousarray(np.flip(modpad*pbpad[i], axis=0)), dxy, u[i], v[i]
                                      , dRA=dRArad[i], dDec=dDecrad[i])
            warr[i] = p[7+countpars]
            countpars+=3
        
            if star:
                #Then add the star
                vismodelstar=np.zeros(u[i].size, dtype=np.complex_)
                vismodelstar.real+=fstar
                #Phase shift star in the visibilities in the same way as model is being shifted by Galario. So star is in geometric center of ring.
                theta = u[i]*2.0*np.pi*(dRArad[i]) + v[i]*2.0*np.pi*(dDecrad[i])
                vismodelstar = (np.real(vismodelstar) + 1j*np.imag(vismodelstar)) * (np.cos(theta) + 1j*np.sin(theta)) #NB This is a multiplication in a complex number sense! Which is the same as a rotation in [Re, Im] space
                vismodel[i]+=vismodelstar
                
        if ngal>=1:
            for i in np.arange(ngal):
                
                fbkg, dRAbkgrad, dDecbkgrad = p[5+countpars], p[6+countpars]/3600.0*np.pi/180.0, p[7+countpars]/3600.0*np.pi/180.0
                
                if resolved[i]:
                    sigmagal, PAgal, incgal = p[8+countpars], p[9+countpars]*deg, p[10+countpars]*deg
                    #Then add the bkg source as a 2D Gaussian inclined with Galario functions
                    fgal = GaussianProfileGalaxy(1.0, sigmagal, Rmin, dR, nR)
                    #Make a separate image with offset source
                    imbkgonly = sweep(fgal, Rminrad, dRrad, nxy, dxy, inc=incgal)
                    for j in np.arange(nvis):
                        visbkg = sampleImage(imbkgonly*fbkg*
                                             pbpad[j][np.int(pbpad.shape[0]/2.0+p[7+countpars]/pxsz), 
                                                   np.int(pbpad.shape[1]/2.0-p[6+countpars]/pxsz)]/np.sum(imbkgonly), 
                                             dxy[j], u[j], v[j], PA=PAgal, dRA=dRAbkgrad+dRArad[j], dDec=dDecbkgrad+dDecrad[j])
                        vismodel[j]+=visbkg
                    countpars+=6
                else:
                    for j in np.arange(nvis):
                        #Then add bkg point source
                        visbkg=np.zeros(u[j].size, dtype=np.complex_)
                        visbkg.real+=fbkg*pbpad[np.int(pbpad.shape[0]/2.0+p[7+countpars]/pxsz), np.int(pbpad.shape[1]/2.0-p[6+countpars]/pxsz)]
                        #Phase shift bkg in the visibilities in the same way as model is being shifted by Galario.
                        theta = u[j]*2.0*np.pi*(dRAbkgrad+dRArad[j]) + v[j]*2.0*np.pi*(dDecbkgrad+dDecrad[j])
                        visbkg = (np.real(visbkg) + 1j*np.imag(visbkg)) * (np.cos(theta) + 1j*np.sin(theta))
                        vismodel[j]+=visbkg
                    countpars+=3  
        
        chi2=[[] for x in vismodel]
        lnL=lnprior
        for i in np.arange(vismodel):
            chi2[i]=np.sum(((np.real(vismodel[i])-Re[i])**2.0+(np.imag(vismodel[i])-Im[i])**2.0)*w[i])
            lnL+=-0.5*(chi2[i]*warr[i] + np.sum(2.0*np.log(2.0*np.pi/warr[i]/w[i])))

        # Here print useful things if needed for postprocessing, as described above
        if locfiles:
            warr=[wtfact]
            np.save('../'+locfiles+'/warr.npy', warr)
            Re_resid=Re-vismodel.real
            Im_resid=Im-vismodel.imag
            np.save('../'+locfiles+'/'+tag+'_uvtable_model.npy', [u, v, vismodel.real, vismodel.imag, w])
            np.save('../'+locfiles+'/'+tag+'_uvtable_resid.npy', [u, v, Re_resid, Im_resid, w])
            np.save('../'+locfiles+'/bestfitshiftPAinc_rad.npy', [incl*np.pi/180.0, posang*np.pi/180.0, dRArad, dDecrad])
            os.chdir('..')
            os.system('rm -r '+ciao)
            return modpad, pxsz, -0.5*(chi2*wtfact) #-0.5*(chi2*wtfact + np.sum(2.0*np.log(2.0*np.pi/wtfact/w)))
    
    
        os.chdir('..')
    
        os.system('rm -r '+ciao)

        # To conclude, return -chisquare/2 + ln(prior), where chisquare is rescaled by a reweighting factor wtfact which we leave as free parameter.
        # The latter is needed as it has been found that the weights delivered by ALMA and extracted from CASA are incorrect by varying factors
        # which depend on the dataset, and this influences the uncertainty on the model parameters that we obtain from fitting the model to the data.
        # The chisquare is also added to np.sum(2.0*np.log(2.0*np.pi/wtfact/w)) which acts to normalise the probability so that its integral is one -
        # -  in other words exp(-chisquare/2)=model probability is 1 when integrated over the entire parameter space. See emcee documentation, but 
        # note that the latter does not affect the outcome of the fit, except (perhaps?) the value of the reweighting factor.
        return lnL




    # Here define prior probability distributions for each parameter. 
    # We choose uninformative priors, so the probability is 1 for a model where all parameters are within their allowed ranges (which we define here)
    # And 0 if we fall outside of these ranges. In other words, ln(prior)=0 if parameters fall in the allowed ranges (regardless of what they are)
    # and ln(prior)=-infinity if any of them falls outside the allowed ranges.
    def lnpriorfn(p):
        for i in range(len(p)):
            if p[i] > priors_up[i] or p[i] < priors_dwn[i]:
                return -np.inf
        return 0.0

    # Starting point in our parameter space search, for each parameter. Order is
    # fstar, f0, R0, sigma, inc, PA, dRA, dDec, wtfact
    p0 = pars_init

    print('Setting up sampler...')
    ndim = len(pars_init)       	    # number of dimensions of the parameter space = number of model parameters
    nwalkers = ndim*10              # number of walkers that will undergo the MCMC
    #nthreads = 1               # CPU threads that emcee should use. Commented out as 'Pool' computes that automatically
    nsteps = 10    	    # total number of steps to run the MCMC for.

    # Starting point in our parameter space search, for each parameter, here expanded FOR EACH OF THE nwalkers walkers.
    # For each parameter, it's good to start walkers in a 'Gaussian distribution' around the first guess p0, where the first 
    # guess should be 'well guessed' (e.g. by looking at the imaging of the data). 
    pos = [p0 + np.asarray([0.1]*ndim)*np.asarray(p0)*np.random.randn(ndim) for i in range(nwalkers)]


    # Define 'backend' which is practically a file which holds the result of the MCMC computation AS IT IS RUNNING.
    # This ensures we can recover the result not only in the future, but also if something goes wrong and the MCMC crashes
    backend=emcee3.backends.HDFBackend(backendaddress)
    # This command wipes what's currently stored in the backend opened above - so make sure you don't use this command on a backend 
    # containing something important!
    #backend.reset(nwalkers,ndim)

    #import the Pool function which allows the computation to be distributed amongst different cores of the machine we are running on,
    #saving a considerable amount of time.
    #Usage of different cores can be seen through the terminal command 'mpstat -P ALL 1' while the MCMC is running.
    from multiprocessing import Pool

    # Define object of EnsembleSampler class from emcee package. Basically this is the object the computation gets run onto.
    # See https://emcee.readthedocs.io/en/latest for a full description of the object, what function/parameters it has, and
    # in general how the MCMC fit works.
    sampler = emcee3.EnsembleSampler(nwalkers, ndim, lnpostfn, pool=Pool(), backend=backend)#, pool=Pool())

    return pos, sampler

