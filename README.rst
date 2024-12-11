******************************************************************
Joint multi-band deconvolution for Euclid and Vera C. Rubin images
******************************************************************

We introduce a novel multi-band deconvolution technique aimed at improving the resolution of ground-based astronomical images by leveraging higher-resolution space-based observations. The method capitalises on the fortunate fact that the Vera C. Rubin r-, i-, and z-bands lie within the Euclid VIS band. The algorithm jointly deconvolves all the data to turn the the r-, i-, and z-band Vera C. Rubin images to the resolution of Euclid.

Installation
============

1) `Download and install Miniconda <http://conda.pydata.org/miniconda.html>`_. Choose the Python 3.x version for your platform.

2) Open a Terminal (Linux/macOS) or Command Prompt (Windows) and run the following commands:

    .. code-block:: bash

        conda update conda
        conda install git
        git clone https://github.com/utsav-akhaury/Multiband-Deconvolution
        cd Multiband-Deconvolution

3) Create a conda environment and install all the required dependencies by running the following commands:

    .. code-block:: bash

        conda env create -f conda_env.yml

Code Overview
=============

.. code-block:: bash

    Multiband-Deconvolution/
        Data/
            euclid.npy
            noisemap_LSST.npy
            noisy_LSST.npy
            psf_euclid_vis.npy
            psf_LSST.npy
            sed.npy
            target_HST.npy
        README.rst
        conda_env.yml
        MBDeconv_FISTA.py
        run_MCDeconv.ipynb

* `Data <https://github.com/utsav-akhaury/Multiband-Deconvolution/tree/main/Data>`_ is the directory containing the test images used in the tutorial notebook.
    * ``euclid.npy`` is the Euclid VIS-band image.
    * ``noisemap_LSST.npy`` is the noise map of the LSST r-, i-, and z-band images.
    * ``noisy_LSST.npy`` is the low-resolution LSST image in r-, i-, and z-bands.
    * ``psf_euclid_vis.npy`` is the Euclid VIS-band PSF.
    * ``psf_LSST.npy`` is the LSST PSF in each LSST band at Euclid resolution.
    * ``sed.npy`` is the fractional contribution of each LSST band to the Euclid VIS band.
    * ``target_HST.npy`` is the target high-resolution HST image.
* `README.rst <https://github.com/utsav-akhaury/Multiband-Deconvolution/blob/main/README.rst>`_ contains getting started information on installation and usage.
* `conda_env.yml <https://github.com/utsav-akhaury/Multiband-Deconvolution/blob/main/conda_env.yml>`_ is a configuration file for Anaconda (Miniconda) that sets up a Python environment with all the required Python packages for using the Multi-band Deconvolution code.
* `MBDeconv_FISTA.py <https://github.com/utsav-akhaury/Multiband-Deconvolution/blob/main/MBDeconv_FISTA.py>`_ contains the implementation of the Multi-band Deconvolution algorithm.
* `run_MBDeconv.ipynb <https://github.com/utsav-akhaury/Multiband-Deconvolution/blob/main/run_MCDeconv.ipynb>`_ is a Jupyter notebook that demonstrates an example of how to deconvolve the simulated LSST images using the Euclid VIS-band image as a high-resolution prior.

Usage
=====

1) Activate the ``mbdeconv`` conda environment:

    .. code-block:: bash

        conda activate mbdeconv

2) Run the `run_MBDeconv.ipynb <https://github.com/utsav-akhaury/Multiband-Deconvolution/blob/main/run_MBDeconv.ipynb>`_ notebook, which will guide you through the deconvolution process.