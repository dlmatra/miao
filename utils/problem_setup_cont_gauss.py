#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def problem_setup(theta):

    import numpy as np
    #IMPORTANT NOTE. CODE ONLY DEALS WITH HEIGHTS UP TO 30 DEG ABOVE MIDPLANE AT THE MOMENT. THAT SHOULD BE OK AS
    #h<~0.15
    Mdisk=theta[0]
    rmid=theta[1] #AU
    sigma=theta[2] #AU
    h=theta[3]
    rmin=theta[4] #cm
    rmax=theta[5] #cm
    nr=theta[6]
    nth=theta[7]-1

    #
    # Some natural constants
    #
    AU  = 1.49598e13     # Astronomical Unit       [cm]
    ME  = 5.97219e27 # Earth mass[g]

    ###################################################################################
    ############################# MAKE POSITIONAL GRID #######################
    ####################################################################################

    #Radial
    if rmax<(rmid*AU+2*sigma*AU):
        print('WARNING - Grid may not be big enough. Max radius set to '+str(rmax/AU)+' but outer 2sigma of Gaussian ring is '+str(rmid+2*sigma))
    ri  = (np.arange(nr+1, dtype=float))/nr*(rmax-rmin)+rmin
    nr  = len(ri)-1
    dr  = (np.subtract(ri[1:],ri[:-1]))
    rc  = np.asarray(ri[:-1]+dr/2e0)
    
    #Vertical
    vertcutoff= 5e0
    be=1.0 #Means aspect ratio constant with radius
    H1AU=h
    
    totalthetacutoff=abs(np.arcsin((np.min([0.99,vertcutoff*H1AU*np.asarray([(rmin/AU)**(be-1.0), (rmax/AU)**(be-1.0)]).max()]))))
    thmax=np.pi/2e0+totalthetacutoff 
    closeto0=1e-2
    logscalconstth=((totalthetacutoff/np.radians(closeto0))**(1e0/(1e0*nth)))-1

    thi  = np.radians(closeto0)*(logscalconstth+1)**np.arange(nth+1, dtype=float)
    thi  = np.concatenate([-thi[::-1],[0.0],thi])+np.pi/2e0
    
    nth = len(thi)-1
    dt  = (np.subtract(thi[1:],thi[:-1]))
    tc  = np.asarray(thi[:-1]+dt/2e0)
    
    wfile=open('./amr_grid.inp', 'w')
    wfile.write('%d\n'%1)                      # Format number
    wfile.write('%d\n'%0)                    # AMR self.style (0=regular self. NO AMR)
    wfile.write('%d\n'%150)                  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical
    wfile.write('%d\n'%0)                    # Gridinfo
    wfile.write('%d %d %d \n'%(1,1,0))       # Which dimension is active
    wfile.write('%d %d %d \n'%(rc.size, tc.size))    # Grid size (x,y,z or r,theta, phi, or r,z, phi)
    for i in range(rc.size+1): wfile.write('%.9e\n'%ri[i])			# X coordinates (cell walls)
    for i in range(tc.size+1): wfile.write('%.9e\n'%thi[i])		# Y coordinates (cell walls)
    wfile.close()							


    ###################################################################################
    ############################# MAKE NUMBER DENSITY GRID #######################
    ####################################################################################
    
    M=Mdisk*ME #IN GRAMS!!
    #gamma	=	gamma
    #rmid=rmid*AU # for Gaussian
    #std=sigma*AU # for Gaussian
    
    rhod= np.zeros([rc.size, tc.size])

    #Create surface density array
    sd = np.exp(-((rc/AU-rmid)**2.0/(2.0*(sigma)**2.0))) # for Gaussian
    #sd = (rc/(rminpowlaw*AU))**gamma #for pow law
    #sd[rc>rmaxpowlaw*AU]=0.0 #for pow law
    #sd[rc<rminpowlaw*AU]=0.0 #for pow law
    #sd = ((rc/(rt*AU))**(-gamma))*((1.0+(rc/(rt*AU))**alpha)**((gamma-betan)/alpha))#for nuker
    
    #Normalize given total mass M
    surfmassdens= sd*M/np.sum(sd*2.0*np.pi*rc*dr)
    
    #pl.figure()
    #pl.plot(rc/AU, surfmassdens)
    
    #Now calculate array of scale heights for each of given radius rc given aspect ratio
    H=(H1AU*((rc/AU)**(be)))*AU
    
    for i in range(rc.size):
        for j in range(tc.size):
            
            #HERE USE VERTICAL GAUSSIAN:
            rhod[i,j]    = np.exp(-(rc[i]*np.sin(((np.pi/2e0)-tc[j]))/(np.sqrt(2.0)*H[i]))**2e0)/(np.sqrt(2.0*np.pi)*H[i])        
        

    #
    # Write the density file
    #
    wfile=open('./dust_density.inp', 'w')
    wfile.write('%d\n'%1)   # Format number
    wfile.write('%d\n'%(rc.size*tc.size))    # Nr of cells
    wfile.write('%d\n'%1)    # Nr of dust species
    for ith in range(tc.size):
        for ir in range(rc.size):
            wfile.write('%.9e\n'%(rhod[ir,ith]))
    wfile.close()


    

    return 0.0

