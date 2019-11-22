from matplotlib.patches import Ellipse, Arrow, Rectangle
import astropy.io.fits as pf
import numpy as np
import matplotlib
import matplotlib.pyplot as pl
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import AnchoredText, AnchoredOffsetbox, AuxTransformBox, VPacker, TextArea
from matplotlib import ticker
c = 299792458.0	   #m/s
AU= 1.495978707e13 #cm


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
from matplotlib import rcParams
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
dist=distobj #pc (from input given above)


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
#Initialise figure within Python
fig=pl.figure(figsize=(6.5,7.0))
# Initialise axes
ax2=fig.add_subplot(111)
# Crop image. multiply by 1e3 so that unit of pixel values is mJy/beam rather than Jy/beam.	
m2=ax2.imshow(1e3*imcontcrop, origin='lower', cmap='inferno', extent=[xvaucontcrop[-1]+np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, xvaucontcrop[0]-np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, yvaucontcrop[0]-np.abs(yvaucontcrop[1]-
yvaucontcrop[0])/2.0, yvaucontcrop[-1]+np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0], interpolation='bicubic', vmax=vmax)

# Below is a bunch of customizable commands to draw a nice colorbar.
divider=make_axes_locatable(ax2)
cax2=divider.append_axes('top', size='5%', pad=0.05)
cbar2=fig.colorbar(m2, cax2, label=r'Flux (mJy beam$^{-1}$)', orientation='horizontal', ticklocation='top')
tick_locator=ticker.MaxNLocator(nbins=5)
cbar2.locator=tick_locator
cbar2.update_ticks()	

#Get contours going. First define contour levels you want (here as a function of number of sigmas):
#Positive contour levels:
poscontlevs=np.array(poslevs)*rmscontim
#Overplot them on top of image as solid lines
ax2.contour(1e3*imcontcrop, poscontlevs, linestyles='solid', colors='black', alpha=0.3, extent=[xvaucontcrop[-1]+np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, xvaucontcrop[0]-np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, yvaucontcrop[0]-np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0, yvaucontcrop[-1]+np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0])
#Negative contour levels:
negcontlevs=np.array(neglevs)*rmscontim
#Overplot them on top of image as dashed lines
ax2.contour(1e3*imcontcrop, negcontlevs, linestyles='dashed', colors='white', alpha=0.3, extent=[xvaucontcrop[-1]+np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, xvaucontcrop[0]-np.abs(xvaucontcrop[1]-xvaucontcrop[0])/2.0, yvaucontcrop[0]-np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0, yvaucontcrop[-1]+np.abs(yvaucontcrop[1]-yvaucontcrop[0])/2.0])
	
#Draw and plot beamsize using class defined at the top of this script
ae = AnchoredEllipse(ax2.transData, width=3.6e3*abs(header_cont['BMAJ']), height=3.6e3*abs(header_cont['BMIN']), angle=90.0+abs(header_cont['BPA']), loc=3, pad=0.5, borderpad=0.4, frameon=False)
ax2.add_artist(ae)

#Set up how you want the x and y ticks and labels to look like
ax2.set_ylabel('North offset (")')
ax2.tick_params(axis='y', direction='in', labelleft='on', right='on', color='white')
ax2.set_xlabel('East offset (")')
ax2.tick_params(axis='x', direction='in', labelbottom='on', top='on', color='white')
#Set which ticks you want on the x axis if you need to
#ax2.set_xticks([20,10,0,-10,-20])

#Draw a bar to give an idea of the physical scale in AU of the belt, as opposed to the size on the sky in arcseconds.
#Here I decided to make it 50 AU in size
asb=AnchoredSizeBar(ax2.transData, sizebarau/dist, str(np.int(sizebarau))+" au", loc=4, pad=0.5,borderpad=0.4, sep=6, frameon=False)
ax2.add_artist(asb)


#For North and East Arrows, ignoring here:
#ax2.add_patch(arrN)
#ax2.text(14.5,-14, 'N', rotation=0.0, color='yellow', size='small', alpha=0.9)
#ax2.text(21,-22, 'E', rotation=0.0, color='yellow', size='small', alpha=0.9)
#ax2.text(23,22, 'HD127821', color='white', alpha=0.9)

#This adds the name of the belt to the top left using the class defined at top of this script
anchored_text=AnchoredText(sourcetag, loc=2, frameon=False, prop={'color':'white', 'alpha':0.9})
ax2.add_artist(anchored_text)

#Add wavelength of observation if you want
#anchored_text=AnchoredText('%03.2f'%(c/restfreqcont*1e3)+' mm', loc=1, frameon=False, prop={'color':'white', 'alpha':0.9})
#ax2.add_artist(anchored_text)


#Save as PDF
pl.savefig(imageplotname)
#Save as PNG
#pl.savefig(name+'_0p5astap_cont.png')
#Close Python image
#pl.close()
# Open image externally from python using terminal command (which can be called from within python using os.system function)
import os
#os.system('evince '+name+'_0p5astap_cont.pdf &')
os.system('cp -r '+imageplotname+' '+workingdir+'/'+sourcetag+'/plots/.')
#os.system('mkdir /data/wdocs/lmatra/www-docs/reasons/HD14055')
#os.system('cp -r ./*.png /data/wdocs/lmatra/www-docs/reasons/HD14055/.')
