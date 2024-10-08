def makeuvdeprojplot(datavisloc, modelvisloc, paradoffsetloc, uvmax, uvbin_size, warr, ylimreal, ylimimag, outfile, sourcetag, phaseshift):
	import matplotlib.pyplot as pl 
	import sys
	import numpy as np
	import os
	tag=sourcetag

	u=[[] for x in datavisloc]
	v=[[] for x in datavisloc]
	Re=[[] for x in datavisloc]
	Im=[[] for x in datavisloc]
	w=[[] for x in datavisloc]
	Re_mod=[[] for x in datavisloc]
	Im_mod=[[] for x in datavisloc]

	for i in np.arange(len(datavisloc)):

		u[i], v[i], Re[i], Im[i], w[i] = np.load(datavisloc[i])
		u[i], v[i], Re_mod[i], Im_mod[i], w_mod = np.load(modelvisloc[i])
		inc_mod, PA_mod, dRA_mod, dDec_mod = np.load(paradoffsetloc[i])
	
		inc=inc_mod
		PA=PA_mod
	
		if phaseshift:
			#########
			# Here shift both model and data using best-fit shifts. This ensures that visibilities are well-centered when making uvplots (and also, importantly, that no relative shifts are present among datasets for multiple 	datasets).
			#########
			theta = np.asarray(u[i])*2.0*np.pi*dRA_mod + np.asarray(v[i])*2.0*np.pi*dDec_mod

			## Minus sign below as we are shifting by the negative of the shift that is usually applied to the model
			vmodelsh=(np.asarray(Re_mod[i])+1j*np.asarray(Im_mod[i]))*(np.cos(-theta)+1j*np.sin(-theta))
			vdatash=(np.asarray(Re[i])+1j*np.asarray(Im[i])) * (np.cos(-theta)+1j*np.sin(-theta))
			Re_mod[i]=np.real(vmodelsh)
			Im_mod[i]=np.imag(vmodelsh)
			Re[i]=np.real(vdatash)
			Im[i]=np.imag(vdatash)
		else:
			print('Not shifting dataset '+str(i))
		#Here fixing weights by factor derived from fitting
		w[i]*=warr[i]
		#Save SHIFTED visibilities
		np.save('./evaluation/'+tag+'_uvtable_modelsh'+str(i)+'.npy', [u[i], v[i], Re_mod[i], Im_mod[i], w[i]])
		np.save('./evaluation/'+tag+'_uvtable_datash'+str(i)+'.npy', [u[i], v[i], Re[i], Im[i], w[i]])

	#Concatenate all datasets for plotting
	u=np.asarray([y for x in u for y in x])#np.concatenate((u1,u2,u3,u4,u5,u6,u7))
	v=np.asarray([y for x in v for y in x])#np.concatenate((v1,v2,v3,v4,v5,v6,v7))
	Re=np.asarray([y for x in Re for y in x])#np.concatenate((Re1,Re2,Re3,Re4,Re5,Re6,Re7))
	Re_mod=np.asarray([y for x in Re_mod for y in x])#np.concatenate((Re_mod1,Re_mod2,Re_mod3,Re_mod4,Re_mod5,Re_mod6,Re_mod7))
	Im=np.asarray([y for x in Im for y in x])#np.concatenate((Im1,Im2,Im3,Im4,Im5,Im6,Im7))
	Im_mod=np.asarray([y for x in Im_mod for y in x])#np.concatenate((Im_mod1,Im_mod2,Im_mod3,Im_mod4,Im_mod5,Im_mod6,Im_mod7))
	w=np.asarray([y for x in w for y in x])#np.concatenate((w1,w2,w3,w4,w5,w6,w7))



	#Now recalculate SHIFTED residuals for plotting. REMEMBER: printed residuals for residual imaging are NOT SHIFTED, i.e., the phase center is not at the best-fit stellar location.
	Re_resid=Re-Re_mod
	Im_resid=Im-Im_mod

	#Carry out deprojection
	cos_inc = np.cos(inc)

	if PA == 0:
		u_deproj=u
		v_deproj=v
	else:
		cospa=np.cos(PA)
		sinpa=np.sin(PA)
		#Rotation by PA here, not of the data, but of u and v!!
		u_deproj=u*cospa-v*sinpa
		v_deproj=u*sinpa+v*cospa

	#NB! Rotation by PA useless if we are not deprojecting by inclination below, as we are just measuring uvdist which doesn't change if we rotate the disk in uv plane.
	u_deproj*=cos_inc
	#print('Deprojection turned ON.')
	#Calculate deprojected u-v distance
	uvdist=np.sqrt(u_deproj**2.0+v_deproj**2.0)

	#Figure out min and max boundaries in uv distance for plotting
	uvdistmin=uvdist.min()
	uvdistmax=np.min([uvdist.max(),uvmax])

	#Bin uvdistances
	nbins=np.ceil((uvdistmax-uvdistmin)/uvbin_size).astype('int')
	nbins_mod=50
	bin_uvdist = np.zeros(nbins)
	bin_weights = np.zeros(nbins)
	bin_uvdist_mod = np.zeros(nbins_mod)
	bin_weights_mod = np.zeros(nbins_mod)
	bin_count = np.zeros(nbins, dtype='int')
	bin_count_mod = np.zeros(nbins_mod, dtype='int')
	uv_intervals=[]
	uv_intervals_mod=[]
	uv_bin_edges=np.arange(nbins+1, dtype='float64')*uvbin_size+uvdistmin
	uv_bin_edges_mod=np.arange(nbins_mod+1, dtype='float64')*((uvdistmax-uvdistmin)/np.float(nbins_mod))+uvdistmin

	#Calculate which data points go in which bin, and figure out average uv distance for that bin for more appropriate plotting
	for i in range(nbins):                                     
		uv_interval = np.where((uvdist >= uv_bin_edges[i]) & (uvdist < uv_bin_edges[i+1]))
		bin_count[i] = len(uv_interval[0])
		if bin_count[i] != 0:                         
			bin_uvdist[i] = uvdist[uv_interval].sum()/bin_count[i]
			bin_weights[i] = np.sum(w[uv_interval])
		else:
			bin_uvdist[i] = uv_bin_edges[i]+0.5*uvbin_size
		uv_intervals.append(uv_interval)
	#Do it for model on a finer grid
	for i in range(nbins_mod):                                     
		uv_interval_mod = np.where((uvdist >= uv_bin_edges_mod[i]) & (uvdist < uv_bin_edges_mod[i+1]))
		bin_count_mod[i] = len(uv_interval_mod[0])
		if bin_count_mod[i] != 0:                         
			bin_uvdist_mod[i] = uvdist[uv_interval_mod].sum()/bin_count_mod[i]
			bin_weights_mod[i] = np.sum(w[uv_interval_mod])
		else:
			bin_uvdist_mod[i] = uv_bin_edges_mod[i]+0.5*((uvdistmax-uvdistmin)/np.float(nbins_mod))
		uv_intervals_mod.append(uv_interval_mod)

	#Calculate real and imaginary of data, model and errors for each bin 
	bin_re, bin_re_err, bin_re_mod, bin_re_modfine, bin_re_resid = np.zeros(nbins), np.zeros(nbins), np.zeros(nbins), np.zeros(nbins_mod), np.zeros(nbins)
	bin_im, bin_im_err, bin_im_mod, bin_im_modfine, bin_im_resid = np.zeros(nbins), np.zeros(nbins), np.zeros(nbins), np.zeros(nbins_mod), np.zeros(nbins) 

	for i in range(nbins):                                     
		if bin_count[i] != 0:   
			bin_re[i] = np.sum(Re[uv_intervals[i]]*w[uv_intervals[i]])/bin_weights[i]
			bin_re_err[i] = 1./np.sqrt(bin_weights[i])                                                          
			bin_im[i] = np.sum(Im[uv_intervals[i]]*w[uv_intervals[i]])/bin_weights[i]
			bin_im_err[i] = 1./np.sqrt(bin_weights[i])
			bin_re_mod[i] = np.sum(Re_mod[uv_intervals[i]]*w[uv_intervals[i]])/bin_weights[i]
			bin_im_mod[i] = np.sum(Im_mod[uv_intervals[i]]*w[uv_intervals[i]])/bin_weights[i]
			bin_re_resid[i] = np.sum(Re_resid[uv_intervals[i]]*w[uv_intervals[i]])/bin_weights[i]
			bin_im_resid[i] = np.sum(Im_resid[uv_intervals[i]]*w[uv_intervals[i]])/bin_weights[i]


	for i in range(nbins_mod):                                     
		if bin_count_mod[i] != 0:   
			bin_re_modfine[i] = np.sum(Re_mod[uv_intervals_mod[i]]*w[uv_intervals_mod[i]])/bin_weights_mod[i]
			bin_im_modfine[i] = np.sum(Im_mod[uv_intervals_mod[i]]*w[uv_intervals_mod[i]])/bin_weights_mod[i]

	#Turn bin x coords to klambda	      
	uvbins = bin_uvdist / 1.e3
	uvbins_mod = bin_uvdist_mod / 1.e3
	maxuvdistklam=uvdistmax/1e3


	#Plot figure
	fig, [ax1, ax2] = pl.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]},figsize=(12.0,12.0))

	ax1.errorbar(uvbins, 1e3*bin_re, yerr=1e3*bin_re_err, fmt='.', linewidth=2.5, capsize=2.5, markersize=13, elinewidth=0., zorder=-32, label='Data')
	#ax1.plot(uvbins, 1e3*bin_re_mod, color='orange', alpha=0.5, label='Model', linewidth=4)
	ax1.plot(uvbins_mod, 1e3*bin_re_modfine, color='orange', alpha=0.5, label='Model', linewidth=4)
	ax1.set_xlim(0,maxuvdistklam)
	ax1.set_ylim(ylimreal[0],ylimreal[1])
	ax1.axvline(uvdist.min()/1.e3, alpha=0.3, linestyle='dashed')
	ax1.axhline(0.0, alpha=0.3, linestyle='solid')
	ax1.tick_params(axis='both', direction='in', top='on', right='on', labelbottom='on')
	ax1.set_ylabel('Real (mJy)')
	#ax1.set_xlabel(r'R$_{\mathrm{uv}}$ (k$\lambda$)')
	ax1.set_xticklabels([])
	ax1.legend(loc=1,frameon=False)
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


	ax2.errorbar(uvbins, 1e3*bin_im, yerr=1e3*bin_im_err, fmt='.', linewidth=2.5, capsize=2.5, markersize=13, elinewidth=0., zorder=-32)
	#ax2.plot(uvbins, 1e3*bin_im_mod, color='orange', alpha=0.5, linewidth=4)
	ax2.plot(uvbins_mod, 1e3*bin_im_modfine, color='orange', alpha=0.5, linewidth=4)
	ax2.set_xlim(0,maxuvdistklam)
	ax2.set_ylim(ylimimag[0],ylimimag[1])
	ax2.axvline(uvdist.min()/1.e3, alpha=0.3, linestyle='dashed')
	ax2.axhline(0.0, alpha=0.3, linestyle='solid')
	ax2.tick_params(axis='both', direction='in', top='on', right='on', labelbottom='on')
	ax2.set_xlabel(r'R$_{\mathrm{uv}}$ (k$\lambda$)')
	ax2.set_ylabel('Imaginary (mJy)')

	pl.savefig(outfile)
	#pl.close()
	#os.system('evince '+outfile+' &')


	return uvbins, uvbins_mod, bin_re, bin_re_err, bin_re_modfine, bin_im, bin_im_err, bin_im_modfine
