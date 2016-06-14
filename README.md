==============
ElecSus v2.2.0
==============

A program to calculate the electric susceptibility of an atomic ensemble.
The program is designed to model weak probe laser spectra on the D-lines
of thermal alkali metal vapour cells. The program also includes fitting 
routines which allow experimental parameters to be extracted from 
experimental spectroscopic data.

An alternative GUI is provided by this fork. 

2D-Plot:
<img src=https://github.com/matwid/ElecSus/tree/master/elecsus/images/ubuntu_2d_cesium.png>

3D: Poincare Sphere:
[![Screenshot](https://github.com/matwid/ElecSus/tree/master/elecsus/images/ubuntu_3d_sodium.png)]

1D: Single Spectrum:
[![Screenshot](https://github.com/matwid/ElecSus/tree/master/elecsus/images/windows7_cesium_ix.png)]

--------------------
New in version 2.0
--------------------

	-	Significantly improved user-friendliness in the form of a 
		GUI to aid in calculating theory spectra and fitting 
		experimental data. 
		
		Works on Windows and Linux, tested on Windows 7, 8.1, 10
		and Ubuntu 14.04. Currently not tested on Mac.

	-	Rewritten fitting routines MLFittingRouine.py,
		RRFittingRoutine.py and SAFittingRoutine.py to support
		keyword arguments, passed to scipy.curve_fit / leastsq
		methods

	-	Added new support modules:
		- elecsus_methods.py
		- libs/data_proc.py
		
		elecsus_new.py contains two simplified methods for either
		calculating spectra or fitting data, and should be easier to 
		interface with external code for, e.g., batch processing / 
		fitting of data or generating 2D plots.

		data_proc.py contains methods for binning (reducing the 
		number of data points by local averaging) and moving-average 
		smoothing data traces

		both of these new modules are used by the GUI program
		
	-	Renamed the old elecsus.py module for added clarity
	
		elecsus.py --> elecsus_runcard.py
		
		This is the old method of calling elecsus with <runcard>.py files as system arguments.
		This way is now obsolete, being replaced by either the GUI or the methods contained in
		elecsus_methods.py. For backwards compatibility, the elecsus_runcard.py module allows
		the runcards to be used in the same way as before.
		
		The example runcards and example data have been moved to sub-directories, 
		/runcard and /sample_data, respectively.
		
	-	Added a new module, spectra_Efield.py
		
		This module allows calcualtion of spectra with the output of electric field vectors, rather
		than spectroscopic quantities. This should allow calculation of spectra in cells with non-uniform
		magnetic fields, by splitting the cell into sufficiently small parts that the field variation
		across any one part is negligible. Spectroscopic quantities can be calculated from the electric
		field vector by using Jones matrices.
		

-------------
Prerequisites
-------------

Must have the python programming language installed with the following 
packages:

	- Scipy version 0.12.0 or later
	- Numpy
	- Matplotlib
	- wxPython >= 2.8 (for GUI)


------------
Installation
------------

Python and required packages must be installed prior to installing ElecSus.

	- Download the zip file and extract the ElecSus directory.

	- For windows, there are pre-built binaries which will install ElecSus for you. Simply double click on the installer exe/msi file and follow the instructions.
	
	- For linux-based systems, download or clone this repository and navigate to the download location in a terminal window. Install using the setup.py file by typing
		
		python setup.py install
	
	- Note the GUI part of ElecSus is currently untested on Mac OSX!

-----
Usage
-----

	1. For GUI operation:

	- After package installation, from the python interpreter type:
	
		from elecsus import elecsus_gui
		
		elecsus_gui.start()

	- In windows, double-click on the run_gui.bat file in the elecsus directory

	- Alternately, open a terminal or command-line window and navigate to the ElecSus directory. Type:

		python elecsus_gui.py
			
	2. For Runcard operation:

	- Open a terminal window and move to the directory where the files have been extracted to.

	- To run the program taking parameters from runcard.py type:

		python elecsus_runcard.py

	- To run using parameters from a particular runcard type

		python elecsus_runcard.py <run card file name>

	- So to run a the first D1 example, type

		python elecsus.py runcard_D1sample.py

	- To run the second example, type

		python elecsus.py runcard_D2sample.py
		

	3. For integration into external code:
	
	- The elecsus_methods.py module contains two methods, calculate() and fit_data(), 
	  which allow for easy integration into external codes. See the elecsus_methods.py source
	  for more details

------
Manual
------

For GUI documentation, see docs/ElecSus_GUI_UserGuide.pdf

For the ElecSus paper, go to http://dx.doi.org/10.1016/j.cpc.2014.11.023
and download the paper. It is published open-access and therefore freely available.

-------
License
-------

All the files distributed with this program are provided subject to the
Apache License, Version 2.0. A Copy of the license is provided.

-----------
Change Log
-----------

V 2.2.0

	- GUI version number and ElecSus version number are now the same
	- Since Enthought Canopy now ships with wxPython version 3.0, GUI has been
		updated to work with this version of wxPython. All previous functionality should 
		remain the same, but a few things have changed:
			- Theory/Fit plot selections are no longer Transient Popups - they are now Dialog Boxes
			- Default plot scaling may look a bit odd when viewing on small resolution monitors -
				not sure what the real fix for this is, but as a workaround, resizing the GUI window
				should fix this
	- Added ability to use experimental detuning axis on theory curve, 
		for direct comparison when outputting data (Menu > Edit > Use Experimental Detuning Axis)
	- Added ability to turn off automatic scaling of axes (Menu > View > Autoscale)
	- Fixed an issue where save files would not save correctly after performing a fit
	- Minor fix for an issue where starting from the python interpreter would not exit properly on clicking the 'X'
	- Corrected some incorrect tooltips
	- Added show_versions() method to elecsus_gui.py, which shows the currently used version numbers of 
		modules important to running this program
		
V 2.1.0

	- Cleaned up a lot of code in the spectra.py module. Added new methods (calc_chi, get_spectra) to spectra.py that allow users to easily create wrapper methods to return custom data types that aren't returned by default. Separated the method to calculate electric field propagation (get_Efield), for use with non-uniform B fields. Backward compatibility is preserved with the spectrum() method, which is now just a wrapper for get_spectra().
	
V 2.0.3

	- Bug fix for the runcard method of using ElecSus. Now supports runcards in directories other than the local directory.
	
V 2.0.2

	- Minor bug fix for the GUI - fixed an issue where the Phi plots would not be plotted
	
V 2.0.1

	- Updated Eigensystem.py with correct fine structure constant for Rb.