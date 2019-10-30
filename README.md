## MIAO - Modelling Interferometric Array Observations
CASA imaging and MCMC fitting of interferometric visibilities using - for now - optically thin 3D ring models

**Instructions**: 
Bringing your data from archive to paper-ready plots requires several steps, which should be followed and executed within the 4 Notebooks in the tutorial folder, which will guide us, respectively, through imaging, visibility modelling, postprocessing and final plotting (in that order). It is probably a good idea to follow these in order, even if, for example, you have already carried out your own imaging elsewhere.
As well as **installing the required packages** below, it is also a good idea (to keep things tidy) to create a new folder for the target you are modelling, and copy the tutorial folder of this package within them, so you can execute things from there.

**NOTE**:
This has been tested on single-dataset analysis, and not yet on multiple-dataset analysis. **Use at own risk!**



## Requirements

Before getting started, here is a list of required packages for your Python installation, including versions that this code was tested with. The code was built in **Python 3.7.3**, but may handle other versions with some tweaks here and there. Required packages:

**CASA**:
https://casa.nrao.edu/casa_obtaining.shtml
Needed for handling visibility data from interferometric observatories such as ALMA and the SMA.
**Needs to be installed on your local machine**, and its /bin folder path specified within the first tutorial (imaging, Step 1) This code was tested with version 5.4.0.

**galario**:
https://github.com/mtazzari/galario
**needs to be installed within local Python 3 installation**, e.g. with conda with a command like: 'conda install -c conda-forge galario'. This code was tested with version 1.2

**RADMC-3D version 0.41**:
http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/index.html
comes with the package, in a version that has been slightly modified to avoid print statements if nothing goes wrong, and hence save time. **This however needs to be installed** through: 'cd miao/radmc-3d/version_0.41/srcnoprint' and then 'make'. It is strongly advised that the user familiarises themselves with the ray-tracing capabilities of radmc-3d, although this knowledge is not strictly needed to run miao.
The **radmc3dPy** add-on: https://www.ast.cam.ac.uk/~juhasz/radmc3dPyDoc/index.html is used to read radmc-3d files into Python, and is included within the radmc-3d package.

**emcee**:
https://emcee.readthedocs.io/en/stable/
**needs to be installed within local Python 3 installation**, and is used for MCMC fitting of the data.
Installation is as simple as: 'conda install -c astropy emcee'. This code was tested with version 3.0.0

**Astropy**:
https://www.astropy.org/
**needs to be installed within local Python 3 installation**, to enable read and write of FITS format files. Installation is as simple as: 'conda install astropy' or 'pip install astropy'. This code was tested with version 3.0.4

**Tqdm**:
https://pypi.org/project/tqdm/
**needs to be installed within local Python 3 installation**, to enable a progress bar to appear during MCMC run.
Installation is as simple as: 'conda install -c conda-forge tqdm'. This code was tested with version 4.36.1

**Corner**:
https://corner.readthedocs.io/en/latest/
**needs to be installed within local Python 3 installation**, to enable production of corner plots.
Installation is as simple as: 'conda install -c astropy corner'. This code was tested with version 2.0.1
