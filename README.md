# radmc-gala
CASA imaging and MCMC fitting of interferometric visibilities using optically thin 3D disk models

Bringing your data from archive to paper-ready plots requires several steps, which we will follow and execute in a series of Notebooks. Built in **Python 3.7.3 and with CASA 5.4.0**, but is likely to handle other versions.

**GALARIO package**:
https://github.com/mtazzari/galario
**needs to be installed within local Python 3 installation**, e.g. with conda with a command like: 'conda install -c conda-forge galario'. This code was tested with version 1.2

**RADMC-3D version 0.41**:
http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/index.html
comes with the package, in a version that has been slightly modified to avoid print statements if nothing goes wrong, and hence save time. **This however needs to be installed** through: 'cd radmc-gala/radmc-3d/version_0.41/srcnoprint' and then 'make'. It is strongly advised that the user familiarises themselves with the ray-tracing capabilities of radmc-3d, although this knowledge is not strictly needed to run radmc-gala.
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
