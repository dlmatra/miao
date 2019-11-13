#Nice fonts
font={'family':'Times New Roman', 'size':24}	
matplotlib.rc('font', **font)
matplotlib.rc('axes', linewidth=2)
matplotlib.rc('xtick.major', width=2, size=8)
matplotlib.rc('xtick.minor', width=2, size=4)
matplotlib.rc('ytick.major', width=2, size=8)
matplotlib.rc('ytick.minor', width=2, size=4)
matplotlib.rc('xtick', labelsize=22)
matplotlib.rc('ytick', labelsize=22)
matplotlib.rcParams.update({'figure.autolayout': True})


burn=burnin


nwalkers=chain.shape[1]
ndim=chain.shape[2]
#Check that the two names below roughly match for consistency
locfiles='./evaluation'
#Bring star and disk flux to mJy
print('Changing units for plotting')
#chain[:,:,0]*=1e6
chain[:,:,[x for x in range(len(labels)) if labels[x]=='fnu'][0]]*=1e3
labelparams[[x for x in range(len(labels)) if labels[x]=='fnu'][0]]=r"$F_{\nu_{\rm belt}}$ (mJy)"
if star:
	chain[:,:,[x for x in range(len(labels)) if labels[x]=='fnustar'][0]]*=1e6
	labelparams[[x for x in range(len(labels)) if labels[x]=='fnustar'][0]]=r'$F_{\nu_{\ast}}$ (uJy)'
chain[:,:,[x for x in range(len(labels)) if labels[x]=='r'][0]]*=dist
labelparams[[x for x in range(len(labels)) if labels[x]=='r'][0]]=r"$R$ (au)"
chain[:,:,[x for x in range(len(labels)) if labels[x]=='sigr'][0]]*=dist*2.0*np.sqrt(2.0*np.log(2.0))
labelparams[[x for x in range(len(labels)) if labels[x]=='sigr'][0]]=r"$\Delta R$ (au)"
pachange=0
if np.median(chain[burn:,:,5])<0:
	chain[burn:,:,[x for x in range(len(labels)) if labels[x]=='pa'][0]]+=180.0
	pachange=1
if ngal>=1:
	print('Leaving bkg fluxes in Jy, for now')
	#chain[:,:,[x for x in range(len(labels)) if 'fbkg' in labels[x]]*=1e3
	

print('Setting labels for ndim parameters')
if useh:
    labelssmall=[labelparams[[x for x in range(len(labels)) if labels[x]=='fnu'][0]], labelparams[[x for x in range(len(labels)) if labels[x]=='r'][0]], labelparams[[x for x in range(len(labels)) if labels[x]=='sigr'][0]], labelparams[[x for x in range(len(labels)) if labels[x]=='h'][0]], labelparams[[x for x in range(len(labels)) if labels[x]=='i'][0]], labelparams[[x for x in range(len(labels)) if labels[x]=='pa'][0]]]
    correspondingparameterindicesforchain=[[x for x in range(len(labels)) if labels[x]=='fnu'][0],[x for x in range(len(labels)) if labels[x]=='r'][0],[x for x in range(len(labels)) if labels[x]=='sigr'][0],[x for x in range(len(labels)) if labels[x]=='h'][0],[x for x in range(len(labels)) if labels[x]=='i'][0],[x for x in range(len(labels)) if labels[x]=='pa'][0]]
else:
    labelssmall=[labelparams[[x for x in range(len(labels)) if labels[x]=='fnu'][0]], labelparams[[x for x in range(len(labels)) if labels[x]=='r'][0]], labelparams[[x for x in range(len(labels)) if labels[x]=='sigr'][0]], labelparams[[x for x in range(len(labels)) if labels[x]=='i'][0]], labelparams[[x for x in range(len(labels)) if labels[x]=='pa'][0]]]
    correspondingparameterindicesforchain=[[x for x in range(len(labels)) if labels[x]=='fnu'][0],[x for x in range(len(labels)) if labels[x]=='r'][0],[x for x in range(len(labels)) if labels[x]=='sigr'][0],[x for x in range(len(labels)) if labels[x]=='i'][0],[x for x in range(len(labels)) if labels[x]=='pa'][0]]
    

###########
#make trace plots, select samples
###########


# Trace plots
fig,ax = pl.subplots(ndim,2,figsize=(18,25),sharex='col',sharey=False)
for j in range(chain.shape[1]):
    #ax[-1,0].plot(sampler.lnprobability[j,:burn])
    for i in range(chain.shape[2]):
        ax[i,0].plot(chain[:burn,j,i], alpha=0.1)
        ax[i,0].set_ylabel(labelparams[i])

for j in range(chain.shape[1]):
    #ax[-1,1].plot(sampler.lnprobability[j,burn:])
    for i in range(chain.shape[2]):
        ax[i,1].plot(chain[burn:,j,i], alpha=0.1)
        ax[i,1].set_ylabel(labelparams[i])

ax[-1,0].set_xlabel('burn in')
ax[-1,1].set_xlabel('sampling')
fig.savefig(traceplotname+'.pdf')
#pl.close()
#os.system('evince '+traceplotname+'.pdf &')
os.system('cp -r '+traceplotname+'.pdf ../../plots/.')

#
print('NB: Code is assuming that you chose the end of burn-in phase correctly.')






###########
# Make corner plot, select best-fit parameters
###########
# Remove burn-in
samples=chain[burn:, :, :].reshape((-1, ndim))

# Make fonts smaller for corner plot
font={'family':'Times New Roman', 'size':16}	
matplotlib.rc('font', **font)
matplotlib.rc('axes', linewidth=2)
matplotlib.rc('xtick.major', width=2, size=4)
matplotlib.rc('xtick.minor', width=2, size=4)
matplotlib.rc('ytick.major', width=2, size=4)
matplotlib.rc('ytick.minor', width=2, size=4)
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams.update({'figure.autolayout': True})


#Produce full corner with shifts and all
fig, ax=pl.subplots(ndim,ndim,figsize=(15,15), gridspec_kw={'wspace': 0.03})
corner.corner(samples, labels=labelparams, fig=fig, quantiles=[0.16, 0.50, 0.84], label_kwargs={'labelpad':20, 'fontsize':12}, show_titles=False, title_kwargs={'fontsize':12}, levels=1.0 - np.exp(-0.5 * (np.arange(3)+1.0) ** 2))
for i in np.arange(ndim):
	for j in np.arange(ndim):
		ax[i,j].tick_params(axis='both', direction='in', labelsize=12, length=6)
#pl.tight_layout()
fig.savefig(cornername+'.pdf')
#pl.close()
#os.system('evince '+cornername+'.pdf &')
os.system('cp -r '+cornername+'.pdf ../../plots/.')

#Produce smaller corner
fig, ax=pl.subplots(len(correspondingparameterindicesforchain),len(correspondingparameterindicesforchain),figsize=(20*len(correspondingparameterindicesforchain)/(ndim-1.0),20*len(correspondingparameterindicesforchain)/(ndim-1.0)), gridspec_kw={'wspace': 0.03})
corner.corner(samples[:,correspondingparameterindicesforchain], labels=labelssmall, fig=fig, quantiles=[0.16, 0.50, 0.84], label_kwargs={'labelpad':20, 'fontsize':12}, show_titles=False, title_kwargs={'fontsize':12}, levels=1.0 - np.exp(-0.5 * (np.arange(3)+1.0) ** 2))
for i in np.arange(len(correspondingparameterindicesforchain)):
	for j in np.arange(len(correspondingparameterindicesforchain)):
		ax[i,j].tick_params(axis='both', direction='in', labelsize=12, length=6)
#pl.tight_layout()
pl.savefig(cornername+'_small.pdf')
#pl.close()
#os.system('evince '+cornername+'_small.pdf &')
os.system('cp -r '+cornername+'_small.pdf ../../plots/.')

# Reset font to standard
font={'family':'Times New Roman', 'size':24}	
matplotlib.rc('font', **font)
matplotlib.rc('axes', linewidth=2)
matplotlib.rc('xtick.major', width=2, size=8)
matplotlib.rc('xtick.minor', width=2, size=4)
matplotlib.rc('ytick.major', width=2, size=8)
matplotlib.rc('ytick.minor', width=2, size=4)
matplotlib.rc('xtick', labelsize=22)
matplotlib.rc('ytick', labelsize=22)
matplotlib.rcParams.update({'figure.autolayout': True})


