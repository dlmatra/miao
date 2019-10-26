from matplotlib.patches import Ellipse, Arrow, Rectangle
import astropy.io.fits as pf
import numpy as np
import matplotlib.pyplot as pl
import sys
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import AnchoredText, AnchoredOffsetbox, AuxTransformBox, VPacker, TextArea
from matplotlib import ticker
c = 299792458.0	   #m/s
AU= 1.495978707e13 #cm


#Image setup
fig=pl.figure(figsize=(30,8.0))
grid=pl.GridSpec(4,4,wspace=0.1, hspace=0.1, height_ratios=[0.9,8,4,0.21])
ax1 = fig.add_subplot(grid[:,0])
ax2 = fig.add_subplot(grid[:,1])
ax3 = fig.add_subplot(grid[:,2])
ax41 = fig.add_subplot(grid[1,3])
ax42 = fig.add_subplot(grid[2,3])


# Define class to draw ellipse corresponding to size of the resolution element in the right place in the image
class AnchoredEllipse(AnchoredOffsetbox):
    def __init__(self, transform, width, height, angle, loc,
                 pad=0.1, borderpad=0.1, prop=None, frameon=True):
        """
        Draw an ellipse the size in data coordinate of the give axes.
        pad, borderpad in fraction of the legend font size (or prop)
        """
        self._box = AuxTransformBox(transform)
        self.ellipse = Ellipse((0, 0), width, height, angle, linewidth=2, fill=False, color='white', alpha=0.9)
        self._box.add_artist(self.ellipse)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=self._box,
                                   prop=prop,
                                   frameon=frameon)
# Define class to draw 'scale bar' in AU in the right place in the image
class AnchoredSizeBar(AnchoredOffsetbox):
    def __init__(self, transform, size, label, loc,
                 pad=0.1, borderpad=0.1, sep=2, prop=None, frameon=True):
        """
        Draw a horizontal bar with the size in data coordinate of the given
        axes. A label will be drawn underneath (center-aligned).

        pad, borderpad in fraction of the legend font size (or prop)
        sep in points.
        """
        self.size_bar = AuxTransformBox(transform)
        self.size_bar.add_artist(Rectangle((0, 0), size, 0, ec="white", linewidth=2.0, alpha=0.9))

        self.txt_label = TextArea(label, minimumdescent=False, textprops={'color':'white', 'alpha':0.9})

        self._box = VPacker(children=[self.txt_label, self.size_bar],
                            align="center",
                            pad=0, sep=sep)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=self._box,
                                   prop=prop,
                                   frameon=frameon)


#'Beauty' plot parameters
rcParams['font.family'] = 'Times New Roman'
font={'size':25}	
matplotlib.rc('font', **font)
matplotlib.rc('axes', linewidth=3)
matplotlib.rc('xtick.major', width=2, size=8)
matplotlib.rc('xtick.minor', width=2, size=4)
matplotlib.rc('ytick.major', width=2, size=8)
matplotlib.rc('ytick.minor', width=2, size=4)
matplotlib.rc('xtick', labelsize=22)
matplotlib.rc('ytick', labelsize=22)
matplotlib.rcParams.update({'figure.autolayout': True})


########################################
################## PART 1: CLEANED IMAGE
########################################
print('***Processing PART 1: CLEANED IMAGE')

### Read in dust image from fits file, and also read its header which contains information about resolution element, pixel size in arcseconds, etc.
acont, header_cont =  np.asarray(pf.getdata(data,0, header=True))
#CASA images come as 4D cubes (polarization, frequency, spatial y (north-south) direction, spatial x (east-west) direction).
#In the case of dust continuum images we have no polarization or frequency, so removing those dimensions below
imcont=acont[0,0,:,:]
#Number of pixels along x axis = number of pixels along y axis:
npixcont=imcont.shape[0]
#Extract frequency of the observations from header
restfreqcont=header_cont['CRVAL3']
#Calculate pixels per beam, In case later I want to plot in Jy/pix rather than Jy/beam (I don't actually, though).
#This is basically calculated as the area of the ellipse of the resolution element in the CASA viewer in arcseconds,
#divided by the area of a pixel in arcseconds (length of its side in arcseconds, squared)
pixperbeamcont=(np.pi*3.6e3*abs(header_cont['BMAJ'])*3.6e3*abs(header_cont['BMIN'])/4e0/np.log(2e0))/(3.6e3*abs(header_cont['CDELT1']))**2e0
#1D array containing x coordinates of the image, in arcseconds, for each pixel along x direction
xvaucont=(np.arange(npixcont)-(npixcont-1.0)/2.0)*np.abs(header_cont['CDELT1'])*3600.0	
#1D array containing y coordinates of the image, in arcseconds, for each pixel along y direction
yvaucont=(np.arange(npixcont)-(npixcont-1.0)/2.0)*np.abs(header_cont['CDELT1'])*3600.0	
rmscontim=rms #mJy /beam (from input given above)


# set image size that we want (half of the diameter)
maxdistau=halfimsizearcsec #arcsec not au (ignore name) 
# Redefine x coordinates in arcsec for this cropped version of the image
xvaucontcrop=xvaucont[np.where(np.abs(xvaucont)<maxdistau)[0]]
# Redefine y coordinates in arcsec for this cropped version of the image
yvaucontcrop=yvaucont[np.where(np.abs(yvaucont)<maxdistau)[0]]
# Crop image
imcontcrop=imcont[np.where(np.abs(yvaucont)<maxdistau)[0][0]:np.where(np.abs(yvaucont)<maxdistau)[0][-1], np.where(np.abs(xvaucont)<maxdistau)[0][0]:np.where(np.abs(xvaucont)<maxdistau)[0][-1]]

#Stuff below is to draw North and East arrow, ignoring here
#northangle=0.0 #(180.0-pa)*np.pi/180.0
#length=5.0
#arrN=Arrow(17,-19,length*np.cos(northangle), length*np.sin(northangle), width=2, color='yellow', alpha=0.7)
#arrE=Arrow(17,-19,1.1*length*np.cos(northangle+np.pi/2.0), 1.1*length*np.sin(northangle+np.pi/2.0), width=2, color='yellow', alpha=0.7)


# Plot dust continuum image
# Crop image. multiply by 1e3 so that unit of pixel values is mJy/beam rather than Jy/beam.	
m1=ax1.imshow(1e3*imcontcrop, origin='lower', cmap='inferno', extent=[xvaucontcrop[-1]+np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, xvaucontcrop[0]-np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, yvaucontcrop[0]-np.abs(yvaucontcrop[1]-
yvaucontcrop[0])/2.0, yvaucontcrop[-1]+np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0], interpolation='bicubic', vmax=vmax)

# Below is a bunch of customizable commands to draw a nice colorbar.
divider=make_axes_locatable(ax1)
cax1=divider.append_axes('top', size='5%', pad=0.05)
cbar2=fig.colorbar(m1, cax1, label=r'Flux (mJy beam$^{-1}$)', orientation='horizontal', ticklocation='top')
tick_locator=ticker.MaxNLocator(nbins=5)
cbar2.locator=tick_locator
cbar2.update_ticks()	

#Get contours going. First define contour levels you want (here as a function of number of sigmas):
#Positive contour levels:
poscontlevs=np.array(poslevs)*rmscontim
#Overplot them on top of image as solid lines
ax1.contour(1e3*imcontcrop, poscontlevs, linestyles='solid', colors='black', alpha=0.3, extent=[xvaucontcrop[-1]+np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, xvaucontcrop[0]-np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, yvaucontcrop[0]-np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0, yvaucontcrop[-1]+np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0])
#Negative contour levels:
negcontlevs=np.array(neglevs)*rmscontim
#Overplot them on top of image as dashed lines
ax1.contour(1e3*imcontcrop, negcontlevs, linestyles='dashed', colors='white', alpha=0.3, extent=[xvaucontcrop[-1]+np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, xvaucontcrop[0]-np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, yvaucontcrop[0]-np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0, yvaucontcrop[-1]+np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0])
	
#Draw and plot beamsize using class defined at the top of this script
ae = AnchoredEllipse(ax1.transData, width=3.6e3*abs(header_cont['BMAJ']), height=3.6e3*abs(header_cont['BMIN']), angle=90.0-abs(header_cont['BPA']), loc=3, pad=0.5, borderpad=0.4, frameon=False)
ax1.add_artist(ae)

#Set up how you want the x and y ticks and labels to look like
ax1.set_ylabel('North offset (")')
ax1.tick_params(axis='y', direction='in', labelleft='on', right='on', color='white')
ax1.set_xlabel('East offset (")')
ax1.tick_params(axis='x', direction='in', labelbottom='on', top='on', color='white')
#Set which ticks you want on the x axis if you need to
#ax1.set_xticks([20,10,0,-10,-20])

#Draw a bar to give an idea of the physical scale in AU of the belt, as opposed to the size on the sky in arcseconds.
#Here I decided to make it 50 AU in size
asb=AnchoredSizeBar(ax1.transData, sizebarau/dist, str(np.int(sizebarau))+" au", loc=4, pad=0.5,borderpad=0.4, sep=6, frameon=False)
ax1.add_artist(asb)


#For North and East Arrows, ignoring here:
#ax1.add_patch(arrN)
#ax1.text(14.5,-14, 'N', rotation=0.0, color='yellow', size='small', alpha=0.9)
#ax1.text(21,-22, 'E', rotation=0.0, color='yellow', size='small', alpha=0.9)
#ax1.text(23,22, 'HD127821', color='white', alpha=0.9)

#This adds the name of the belt to the top left using the class defined at top of this script
anchored_text=AnchoredText(sourcetag, loc=2, frameon=False, prop={'color':'white', 'alpha':0.9})
ax1.add_artist(anchored_text)

#Add wavelength of observation if you want
#anchored_text=AnchoredText('%03.2f'%(c/restfreqcont*1e3)+' mm', loc=1, frameon=False, prop={'color':'white', 'alpha':0.9})
#ax1.add_artist(anchored_text)




########################################
################## PART 2: MODEL IMAGE
########################################
print('***Processing PART 2: MODEL IMAGE')


#Number of pixels along x axis = number of pixels along y axis:
npixcontmod=modelim.shape[0]
#1D array containing x coordinates of the image, in arcseconds, for each pixel along x direction
xvaucontmod=(np.arange(npixcontmod)-(npixcontmod-1.0)/2.0)*pxsz	
#1D array containing y coordinates of the image, in arcseconds, for each pixel along y direction
yvaucontmod=(np.arange(npixcontmod)-(npixcontmod-1.0)/2.0)*pxsz

# set image size that we want (half of the diameter)
maxdistau=halfimsizearcsec 
# Redefine x coordinates in arcsec for this cropped version of the image
xvaucontmodcrop=xvaucontmod[np.where(np.abs(xvaucontmod)<maxdistau)[0]]
# Redefine y coordinates in arcsec for this cropped version of the image
yvaucontmodcrop=yvaucontmod[np.where(np.abs(yvaucontmod)<maxdistau)[0]]
# Crop image
imcontmodcrop=modelim[np.where(np.abs(yvaucontmod)<maxdistau)[0][0]:np.where(np.abs(yvaucontmod)<maxdistau)[0][-1], np.where(np.abs(xvaucontmod)<maxdistau)[0][0]:np.where(np.abs(xvaucontmod)<maxdistau)[0][-1]]

#Stuff below is to draw North and East arrow, ignoring here
#northangle=0.0 #(180.0-pa)*np.pi/180.0
#length=5.0
#arrN=Arrow(17,-19,length*np.cos(northangle), length*np.sin(northangle), width=2, color='yellow', alpha=0.7)
#arrE=Arrow(17,-19,1.1*length*np.cos(northangle+np.pi/2.0), 1.1*length*np.sin(northangle+np.pi/2.0), width=2, color='yellow', alpha=0.7)

# Plot dust continuum image
# Crop image. multiply by 1e6 so that unit of pixel values is uJy/pix rather than Jy/pix.	
m2=ax2.imshow(1e6*imcontmodcrop, origin='lower', cmap='inferno', extent=[xvaucontmodcrop[-1]+np.abs(xvaucontmodcrop[1]-xvaucontmodcrop[0])/2.0, xvaucontmodcrop[0]-np.abs(xvaucontmodcrop[1]-xvaucontmodcrop[0])/2.0, yvaucontmodcrop[0]-np.abs(yvaucontmodcrop[1]-yvaucontmodcrop[0])/2.0, yvaucontmodcrop[-1]+np.abs(yvaucontmodcrop[1]-yvaucontmodcrop[0])/2.0], interpolation='bicubic')

# Below is a bunch of customizable commands to draw a nice colorbar.
divider=make_axes_locatable(ax2)
cax2=divider.append_axes('top', size='5%', pad=0.05)
cbar1=fig.colorbar(m2, cax2, label=r'Flux ($\mu$Jy pix$^{-1}$)', orientation='horizontal', ticklocation='top')
tick_locator=ticker.MaxNLocator(nbins=5)
cbar1.locator=tick_locator
cbar1.update_ticks()	
	
#Set up how you want the x and y ticks and labels to look like
#Turn off y axis
ax2.set_yticklabels([])
#ax2.set_ylabel('North offset (")')
ax2.tick_params(axis='y', direction='in', labelleft='on', right='on', color='white')
ax2.set_xlabel('East offset (")')
ax2.tick_params(axis='x', direction='in', labelbottom='on', top='on', color='white')
#Set which ticks you want on the x axis if you need to
#ax2.set_xticks([20,10,0,-10,-20])





########################################
################## PART 3: RESIDUAL IMAGE
########################################
print('***Processing PART 3: RESIDUAL IMAGE')

### Read in dust image from fits file, and also read its header which contains information about resolution element, pixel size in arcseconds, etc.
acontres, header_contres =  np.asarray(pf.getdata(residdata,0, header=True))
#CASA images come as 4D cubes (polarization, frequency, spatial y (north-south) direction, spatial x (east-west) direction).
#In the case of dust continuum images we have no polarization or frequency, so removing those dimensions below
imcontres=acontres[0,0,:,:]
# Crop image
imcontrescrop=imcontres[np.where(np.abs(yvaucont)<maxdistau)[0][0]:np.where(np.abs(yvaucont)<maxdistau)[0][-1], np.where(np.abs(xvaucont)<maxdistau)[0][0]:np.where(np.abs(xvaucont)<maxdistau)[0][-1]]

# Plot dust continuum residual image
# multiply by 1e6 so that unit of pixel values is uJy/beam rather than Jy/beam.	
m3=ax3.imshow(1e6*imcontrescrop, origin='lower', cmap='inferno', extent=[xvaucontcrop[-1]+np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, xvaucontcrop[0]-np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, yvaucontcrop[0]-np.abs(yvaucontcrop[1]-
yvaucontcrop[0])/2.0, yvaucontcrop[-1]+np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0], interpolation='bicubic')

# Below is a bunch of customizable commands to draw a nice colorbar.
divider=make_axes_locatable(ax3)
cax3=divider.append_axes('top', size='5%', pad=0.05)
cbar3=fig.colorbar(m3, cax3, label=r'Flux ($\mu$Jy beam$^{-1}$)', orientation='horizontal', ticklocation='top')
tick_locator=ticker.MaxNLocator(nbins=5)
cbar3.locator=tick_locator
cbar3.update_ticks()	

#Get contours going. First define contour levels you want (here as a function of number of sigmas):
#Positive contour levels:
poscontlevs=np.array(poslevs)*rmscontim
#Overplot them on top of image as solid lines
ax3.contour(1e3*imcontrescrop, poscontlevs, linestyles='solid', colors='black', alpha=0.3, extent=[xvaucontcrop[-1]+np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, xvaucontcrop[0]-np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, yvaucontcrop[0]-np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0, yvaucontcrop[-1]+np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0])
#Negative contour levels:
negcontlevs=np.array(neglevs)*rmscontim
#Overplot them on top of image as dashed lines
ax3.contour(1e3*imcontrescrop, negcontlevs, linestyles='dashed', colors='white', alpha=0.3, extent=[xvaucontcrop[-1]+np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, xvaucontcrop[0]-np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, yvaucontcrop[0]-np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0, yvaucontcrop[-1]+np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0])

#Set up how you want the x and y ticks and labels to look like
#Turn off y axis
ax3.set_yticklabels([])
#ax3.set_ylabel('North offset (")')
ax3.tick_params(axis='y', direction='in', labelleft='on', right='on', color='white')
ax3.set_xlabel('East offset (")')
ax3.tick_params(axis='x', direction='in', labelbottom='on', top='on', color='white')
#Set which ticks you want on the x axis if you need to
#ax3.set_xticks([20,10,0,-10,-20])


########################################
################## PART 4: VISIBILITY PROFILE
########################################
print('***Processing PART 4: VISIBILITY PROFILE')


#'Beauty' plot parameters
rcParams['font.family'] = 'Times New Roman'
font={'size':25}	
matplotlib.rc('font', **font)
matplotlib.rc('axes', linewidth=3)
matplotlib.rc('xtick.major', width=2, size=8)
matplotlib.rc('xtick.minor', width=2, size=4)
matplotlib.rc('ytick.major', width=2, size=8)
matplotlib.rc('ytick.minor', width=2, size=4)
matplotlib.rc('xtick', labelsize=22)
matplotlib.rc('ytick', labelsize=22)
matplotlib.rcParams.update({'figure.autolayout': True})



ax41.errorbar(uvbins, 1e3*bin_re, yerr=1e3*bin_re_err, fmt='.', linewidth=2.5, capsize=2.5, markersize=13, elinewidth=0., zorder=-32, label='Data')
ax41.plot(uvbins_mod, 1e3*bin_re_modfine, color='orange', alpha=0.7, label='Model', linewidth=6)
ax41.set_xlim(0,maxuvdistklam/1e3)
ax41.set_ylim(ylimreal[0],ylimreal[1])
#ax41.axvline(uvdist.min()/1.e3, alpha=0.3, linestyle='dashed')
ax41.axhline(0.0, alpha=0.3, linestyle='solid')
ax41.tick_params(axis='both', direction='in', top='on', right='on', labelbottom='on', labelleft='off', labelright='on')
ax41.set_ylabel('Real (mJy)')
ax41.yaxis.set_label_position("right")
#ax41.set_xlabel(r'R$_{\mathrm{uv}}$ (k$\lambda$)')
ax41.set_xticklabels([])
ax41.legend(loc=1,frameon=False)
#Here create a model ring using bessel function cf MacGregor 2015 eps Eri paper
#radring=2.03*38.48 #AU
#widthring=0.23*2.0*np.sqrt(2.0*np.log(2.0))*38.48 #AU
#stardist=38.48
#totfluxring=1.4e-3 #Jy
#rho=2*pi*bin_uvdist
#radring=np.radians(radring/stardist/3600.0)
#widthring=np.radians(widthring/stardist/3600.0)
#sigring=widthring/2.0/np.sqrt(2.0*np.log(2))
#A=totfluxring/2.0/np.pi/sigring**2.0
#remod=2.0*np.pi*A*sigring**2.0*np.exp(-(rho*np.sqrt(2.0)*sigring)**2.0/4)*j0(rho*radring)
#pl.plot(bin_uvdist/1e3, remod, linewidth=2.5, color='red', alpha=0.7)

ax42.errorbar(uvbins, 1e3*bin_im, yerr=1e3*bin_im_err, fmt='.', linewidth=2.5, capsize=2.5, markersize=13, elinewidth=0., zorder=-32)
ax42.plot(uvbins_mod, 1e3*bin_im_modfine, color='orange', alpha=0.7, linewidth=6)
ax42.set_xlim(0,maxuvdistklam/1e3)
ax42.set_ylim(ylimimag[0],ylimimag[1])
#ax42.axvline(uvdist.min()/1.e3, alpha=0.3, linestyle='dashed')
ax42.axhline(0.0, alpha=0.3, linestyle='solid')
ax42.tick_params(axis='both', direction='in', top='on', right='on', labelbottom='on', labelleft='off', labelright='on')
ax42.set_xlabel(r'R$_{\mathrm{uv}}$ (k$\lambda$)')
ax42.set_ylabel('Imaginary (mJy)')
ax42.yaxis.set_label_position("right")



#Save as PDF
pl.savefig(imagecomboname)
#Save as PNG
#pl.savefig(name+'_0p5astap_cont.png')
#Close Python image
#pl.close()
# Open image externally from python using terminal command (which can be called from within python using os.system function)
import os
#os.system('evince '+name+'_0p5astap_cont.pdf &')
os.system('cp -r '+imagecomboname+' '+workingdir+'/'+sourcetag+'/plots/.')

