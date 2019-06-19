# PyFluxPro

# NOTE: This version has been deprecated and is no longer supported.  See https://github.com/OzFlux/PyFluxPro.git for the new, integrated version of PyFluxPro.

Welcome to the repository for the "classic" version of PyFluxPro code.

PyFluxPro is a suite of Python scripts integrated into a single GUI that is designed to simplify and standardise the quality control, post-processing, gap filling and partitioning of data from flux towers.  PyFluxPro has been developed by the OzFlux (http://ozflux.org.au) community in Australia and is used by the community as the operational tool for processing data from the OzFlux network of flux towers.  PyFluxPro is not limited to Australia and can be used for flux tower data collected anywhere in the world.  Using PyFluxPro does not require any knowledge of Python (though we would always recommend people learn Python anyway!), all aspects of the processing can be controlled via the simple GUI and by editing text files.  PyFluxPro can read data from Excel workbooks and CSV files and uses netCDF files (http://www.unidata.ucar.edu/software/netcdf/) for storing intermediate and final output data.

The following documentation gives basic information on how to install and use PyFluxPro.  Details on the algorithms used in PyFluxPro and results for several sites in the Australian OzFlux network can be found in a recent publication in the OzFlux Special Issue of Biogeosciences (https://www.biogeosciences.net/14/2903/2017/).

# Installation and Updating
There are 3 steps to installing PyFluxPro:
* 1. Install Python.
* 2. Install the "git" version control software.
* 3. Clone the PyFluxPro repository using the "git" version control software.

## Installing Python
PyFluxPro is written for Python V2.7 and uses a number of standard and 3rd party Python modules.

OzFlux uses and recommends the Anaconda (https://www.continuum.io/) Python V2.7 distribution.  This Python distribution comes with all of the modules used by PyFluxPro and all except 1 are installed by default.  Adding the 1 required module that is not installed by default is very easy, thanks to the conda package manager, and is explained below.  Using the Anaconda distribution is not essential, just very convenient, and it is possible to use any Python V2.7 environment provided the required modules are installed.  There is a list of the required modules in the /docs folder of this repository.

To install the Anaconda Python V2.7 distribution, follow these steps:
* 1. Download the Anaconda Python V2.7 installer for your operating system from https://www.continuum.io/downloads.
* 2. Follow the instructions on the Anaconda web page to install the Anaconda Python V2.7 distribution.
* 3. Accept all the defaults during the installation, including having Anaconda append the path to this Python installation to your system PATH environment variable.
* 4. The default installation provides everything needed to run PyFluxPro with 1 exception, the module required to read and write netCDF files.  This can be installed as follows:
  1. Open a command line window or terminal session.
  2. At the command prompt, type "conda install netcdf4" and follow the instructions.  Accept all of the defaults and the netCDF module will be installed.
* 5. As of October 2017, the default Anaconda installation no longer installs the ConfigObj module and this will need to be installed using the same method as that used for the netCDF module above:
  1. Open a command line window or terminal session.
  2. At the command prompt, type "conda install configobj" and follow the instructions.  Accept all of the defaults and the ConfigObj module will be installed.

At the end of this process, you should have a functioning installation of the Python language interpeter.

## Installing "git"
The version control program "git" provides a convenient way to install OzFluxQC and to update OzFluxQC once it has been installed.

To install "git", follow these steps:
* 1. Download the "git" installer for your operating system from https://git-scm.com/downloads.
* 2. Follow the instructions on the "git" web page to install the "git" version control software.
* 3. Accept all the defaults during the installation.

## Installing PyFluxPro
PyFluxPro is easily installed using the "git" version control software.  This process is refered to as "cloning" the PyFluxPro repository (this web page).  When PyFluxPro is installed using "git" then "git" can also be used to easily update PyFluxPro to make sure you are always using the most recent version.  This is a good idea because PyFluxPro is frequently updated to fix bugs and add new features.

To install PyFluxPro, follow these steps:
* 1. Open a command line window or terminal session and use the "cd" (shorthand for "change directory") command to navigate to the directory into which you want to install PyFluxPro.  Note that the installation process will create a subdirectory called PyFluxPro in the directory from which the install is run.
* 2. Clone the PyFluxPro repository by typing "git clone https://github.com/OzFlux/PyFluxPro.git" at the command prompt.
* 3. PyFluxPro is now ready to use.

## Updating PyFluxPro
PyFluxPro is still being actively developed and there are frequent changes to fix bugs and add new features.  It is always a good idea to update your installation to the latest version every few days.  Updating PyFluxPro is easy when the installation was done using the "git" version control software.

To update a PyFluxPro installation done by "git", follow these steps:
* 1. Open a command line window or terminal session and use the "cd" command to navigate to the PyFluxPro directory created during the installation step above.  Note that while the install is done from the directory one level above the PyFluxPro directory, the update is done from the PyFluxPro directory.
* 2. Type "git pull origin master" at the command prompt in the PyFluxPro directory.  This will update the PyFluxPro installation.

# Running PyFluxPro
The simplest way to run PyFluxPro is from the command line.

To run PyFluxPro, follow these steps:
* 1. Open a command line window or terminal session and use the "cd" command to navigate to the PyFluxPro directory.
* 2. Type "python PyFluxPro.py" at the command prompt in the PyFluxPro directory.
* 3. After a short time, the PyFluxPro GUI will appear.  This can take a couple of minutes when the program is run for the first time.

# Using PyFluxPro
Coming to the Wiki soon ...
