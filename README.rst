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

Usage
=====

1) Activate the ``mbdeconv`` conda environment:

    .. code-block:: bash

        conda activate mbdeconv