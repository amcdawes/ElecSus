#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
ElecSus-GUI - a TraitsUI interface for ElecSus

Copyleft (C) 2015 Ilja Gerhardt <ilja@quantumlah.org>
                  Matthias Widmann <m.widmann@physik.uni-stuttgart.de>
                  Matthias Niethammer <m.niethammer@physik.uni-stuttgart.de>

This source code is free software; you can redistribute it and/or
modify it under the terms of the GNU Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

This source code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
Please refer to the GNU Public License for more details.

You should have received a copy of the GNU Public License along with
this source code; if not, write to:
Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

----- installation notes ----

You will need ElecSus:

https://github.com/jameskeaveney/ElecSus

Ubuntu 14.04:
sudo apt-get install python-numpy python-wxgtk2.8 python-chaco mayavi2 python-matplotlib python-scipy

---------------------------------

ARCH LINUX:
sudo pacman -S wxpython2.8 python-numpy python-matplotlib python-scipy vtk python2-traitsui python2-pyface python2-apptools python2-envisage python2-configobj python2-setuptools
yaourt -S python2-chaco
pip install mayavi2

There is no package for kiwisolver neither in official repo nor in aur, so download it here:
https://pypi.python.org/pypi/kiwisolver/0.1.2
unpack and install it via command (open shell in unpacked folder in dolphin (Shift+F4)):
sudo python2 setup.py install

simple run in Linux: ./gui.py
(Ubuntu) maybe, you will need a 'chmod 755 gui.py' before 
(Arch)   ETS_TOOLKIT=wx ipython2 then run gui.py


-------------------------------------
Supported Data Formats:
fits : convenient data storing and viewing, for Mathematica (drag and drop) 
       and image viewers. 2x smaller than ascii in this program. http://fits.gsfc.nasa.gov/
       (note: fits is not compatible with Origin, use txt instead)
txt, ascii: the data standard
png, pdf 

"""


#------OS specific imports -------------
import sys,os, platform, codecs
if os.name == 'posix':
   from time import time as timing # Timing for linux or apple
   os.chdir(os.path.dirname(__file__))
else:
   from time import clock as timing # Timing for windoze

if 'Windows' in platform.platform():
    from traits.etsconfig.api import ETSConfig
    ETSConfig.toolkit = 'wx'

sys.path.insert(0, "libs")
if 'ARCH' in platform.platform():
    import wxversion
    wxversion.select("2.8")




#------numpy for sure -------------
import numpy as np
#------chaco plot imports -------------
from   chaco.shell  import *
from   chaco_addons import SavePlot as Plot, SaveTool
from   traits.api   import SingletonHasTraits, Instance, Range, Bool, Array, Str, Enum, Button, on_trait_change, Trait, Float, Dict
from   traitsui.api import View, Item, Group, HGroup, VGroup, VSplit, Tabbed, EnumEditor
from   enable.api   import ComponentEditor, Component, ColorTrait
from   chaco.api    import PlotAxis, CMapImagePlot, ColorBar, LinearMapper, ArrayPlotData, Spectral,prism, PRGn, gray, PlotLabel, summer, hot, RdBu,reverse, DataRange1D, marker_trait
from   chaco.tools  import simple_zoom
import chaco.api
#------3D Mayavi imports -------------
from tvtk.pyface.scene_editor      import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene   import MayaviScene
#------saving and timing imports -----
from utility  import GetSetItemsMixin # for saving 1D plot
from datetime import datetime, timedelta
import time
#------matplotlib imports -------------
from mpltraits import MPL2DImageTrait, MPLInitHandler


import spectra
class Gui (GetSetItemsMixin ):
#class Gui():
    """
    A standard gui.
    Suggested by: m.widmann@physik.uni-stuttgart.de: 2015-04-23
    last modified: 2015-08-26  
    """
    #unicode shorts UTF-16 (hex), because of windows 16bit standard
    sigma  = unichr(0x03C3) # to display sigma sign
    degree = unichr(0x00B0) # to display degree sign
    # Buttons:
    calculate_button          =  Button(label='Calculate', desc='Calculates spectrum')
    calculate_matrix_button   =  Button(label='Calculate matrix', desc='Calculates 2D spectrum')
    save_matrix_image_button  =  Button(label='Save image', desc='Saves matrix image')
    save_matrix_data_button   =  Button(label='Save data', desc='Saves matrix data')


    # parameters here
    frequency_begin           = Float(default_value=-20.,    desc='Begin frequency [GHz]',       label='Begin frequency [GHz]',        mode='text', auto_set=False, enter_set=True)
    frequency_end             = Float(default_value=20.,     desc='End frequency [GHz]',         label='End frequency [GHz]',          mode='text', auto_set=False, enter_set=True)
    frequency_delta           = Float(default_value=10,      desc='Frequency resolution [MHz]',       label='Delta frequency [MHz]',        mode='text', auto_set=False, enter_set=True)
    frequency_range           = Float(default_value=40,      desc='Frequency range[GHz]',       label='Frequency range [GHz]',        mode='text', auto_set=False, enter_set=True)
    frequency_shift           = Float(default_value=0,       desc='Frequency shift [MHz]',       label='Frequency shift [MHz]',        mode='text', auto_set=False, enter_set=True)
    frequency_broadening      = Float(default_value=0,       desc='Frequency broadening [MHz]',  label='Frequency broadening [MHz]',   mode='text', auto_set=False, enter_set=True)
    frequency_width           = Float(default_value=100.,       desc='Frequency broadening [MHz]',  label='Frequency width [MHz]',        mode='text', auto_set=False, enter_set=True)
    length                    = Range(low=0,      high=1e10, value=100.0, desc='Length of the cell [mm]',   label='Length of the cell [mm]',)
    bfield                    = Range(low=-1e10,  high=1e10, value=350.0, desc='B-field [Gauss]',            label='B-field [Gauss]',)
    temp1                     = Range(low=-273.15,high=1e10, value=35.0,  desc='Temperature 1 ['+degree+'C]\nthis determines the number of atoms in the vapor\n(coldest point in the cell)',      label='Temperature 1 ['+degree+'C]',)
    temp2                     = Range(low=-273.15,high=1e10, value=35.0,  desc='Temperature 2 ['+degree+'C]\nthis determined the spectral width of the vapor\n(temperature & velocity of the measured atoms)',      label='Temperature 2 ['+degree+'C]',)
    inputpol                  = Range(low=0.     ,high=360., value=90.0,  desc='Input Polarization ['+degree+']',  label='Input polarization ['+degree+']',)
    detpol                    = Range(low=0.     ,high=360., value=50.0,  desc='Detection Polarization [%]\n0% means that only '+sigma+'+ is detected\n50% means fully linear polarized light\n100% means that only '+sigma+'- light is detected',label='Detection Polarization [%]',)
    element                   = Enum(("Na","K","Rb","Cs"),desc='Atomic vapor under study')
    line                      = Enum(("D1","D2"),value="D2", desc='Line selection (D1 or D2)')
    spectype                  = Enum(("S0","S1","S2","S3","Ix","Iy","RI+","RI-","GI+","GI-","T-","T+","Rot"))
    rb85frac                  = Float(default_value=72.17,     desc='Amount of the isotope\n Rb 85 in the atomic vapor [%]',  label='Rb 85 [%]',  mode='text', auto_set=False, enter_set=True, )
    k40frac                   = Float(default_value=6.7302,    desc='Amount of the isotope\n K 40 in the atomic vapor [%]',   label='K 40 [%]',   mode='text', auto_set=False, enter_set=True, )
    k41frac                   = Float(default_value=0.0117,    desc='Amount of the isotope\n K 41 in the atomic vapor [%]',   label='K 41 [%]',   mode='text', auto_set=False, enter_set=True, )
    state                     = Str('idle', desc='State of the calculations',     label='State',)
    last_calc_time            = Float(default_value=0, desc='Calculation time [ms]',  label='Calculation time [ms]',)
    est_calc_time             = Str('00:00:00:00', desc='Estimated time for matrix calculation [hh:mm:ss]',  label='Estimated calculation time [hh:mm:ss]')
    connect_temp              = Bool(True, desc='Connects temperature 1 and 2', label='connect',)

#------- checkboxes for spectra  ---------------------------------------
    show_S0_spec              = Bool(False,label='S0',desc='Stokes Vector: S0 = I = P_0 + P_90')
    show_S1_spec              = Bool(False,label='S1',desc='Stokes Vector: +1 horizontal, -1 vertical \n S1 =  P_45 - P_90')
    show_S2_spec              = Bool(False,label='S2',desc='Stokes Vector: S2 = P_45 - P_135')
    show_S3_spec              = Bool(False,label='S3',desc='Stokes Vector: Circular component: S3 =  P_'+sigma+'+ -P_'+sigma+'- ')
    show_Ix_spec              = Bool(False,label='Ix',desc='Ix=(S0+S1)/2.0, light intesity transmitted by a PBS')
    show_Iy_spec              = Bool(False,label='Iy',desc='Iy=(S0-S1)/2.0, light intesity reflected by a PBS')
    show_RIp_spec             = Bool(False,label='RI+',desc='Refractive index '+sigma+'+ (real part)')
    show_RIm_spec             = Bool(False,label='RI-',desc='Refractive index '+sigma+'- (real part)')
    show_GIm_spec             = Bool(False,label='GI-',desc='Refractive index '+sigma+'- (imaginary part)')
    show_GIp_spec             = Bool(False,label='GI+',desc='Refractive index '+sigma+'+ (imaginary part)')
    show_Tm_spec              = Bool(False,label='T-',desc='Transmission of the '+sigma+'- component only')
    show_Tp_spec              = Bool(False,label='T+',desc='Transmission of the '+sigma+'+ component only')
    show_Rot_spec             = Bool(False,label='Rot',desc='Rotation of the light against a linear polarizer')
    normalize                 = Bool(False,label='Normalize',desc='Normalize the spectrum')
    normalize_begin           = Float(default_value=0.,     desc='Starting value of normalize',  label='norm. begin',  mode='text', auto_set=False, enter_set=True, )
    normalize_end             = Float(default_value=1.,     desc='End value of normalize',  label='norm. stop',  mode='text', auto_set=False, enter_set=True, )
#------- 3D Sphere parameters  ---------------------------------------
    axes_alpha                = Float(default_value=0.8,  desc='axes_alpha',  label='axes_alpha',  mode='text', auto_set=False, enter_set=True, )
    frame_alpha               = Float(default_value=0.1,  desc='frame_alpha', label='frame_alpha', mode='text', auto_set=False, enter_set=True, )
    frame_radius              = Float(default_value=0.005,desc='frame_radius',label='frame_radius',mode='text', auto_set=False, enter_set=True, )
    axes_radius               = Float(default_value=0.01, desc='axis_radius', label='axis_radius', mode='text', auto_set=False, enter_set=True, )
#------- Boolean parameters for matrix plot ---------------------------------------
    bfield_bool               = Bool(True, label='B-field', desc='B-field [Gauss]')
    temp1_bool                = Bool(False,label='Temperature 1',desc='Temperature 1 ['+degree+'C]\nthis determines the number of atoms in the vapor\n(coldest point in the cell)')
    temp2_bool                = Bool(False,label='Temperature 2',desc='Temperature 2 ['+degree+'C]\nthis determined the spectral width of the vapor\n(temperature & velocity of the measured atoms)')
    length_bool               = Bool(False,label='Length of cell',desc='Length of the cell [mm]')
    inputpol_bool             = Bool(False,label='Input polarization',desc='Input polarization ['+degree+']')
    detpol_bool               = Bool(False,label='Detection polarization',desc='Detection polarization [%]\n0% means that only σ+ is detected\n50% means fully linear polarized light\n100% means that only σ- light is detected')
#------- some necessary or unessary stuff, only god knows -------------------------
    frequency = Array()
#------- Look and feel, line widths etc...-----------------------------------------
    plot_line_width           = Float(default_value=3.,  desc='The width of all lines in the plot',  label='Plot line width',  mode='text', auto_set=False, enter_set=True, )
#---------- Line Plot Data --------------------------------------------------------
    begin           = Float(default_value=0.,   desc='Begin',  label='Begin',   mode='text', auto_set=False, enter_set=True)
    end             = Float(default_value=500., desc='End',    label='End',     mode='text', auto_set=False, enter_set=True)
    delta           = Float(default_value=50.,  desc='Delta',  label='Delta',   mode='text', auto_set=False, enter_set=True)
#---------- automatic change listener --------------------------------------------------------
    @on_trait_change('calculate_button','length','bfield','frequency_range','frequency_delta','frequency_shift','frequency_broadening','length','bfield','temp1','temp2','inputpol','detpol','element','line','rb85frac', 'k40frac','k41frac')
    def calculate(self):
        self.frequency_begin = -self.frequency_range*.5
        self.frequency_end   = self.frequency_range*.5
        self._calculate_all()

    @on_trait_change('calculate_matrix_button')
    def calculate_matrix(self):
        self._calculate_matrix()

    @on_trait_change('save_matrix_image_button')
    def save_matrix_image(self):
        self.matrix_2_plot.save_matrix_image(self.matrix_image_filename+'.png')

    
    @on_trait_change('save_matrix_data_button')
    def save_matrix_data(self):
        self.matrix_2_plot.save_matrix_data(self.matrix_dictionary2,self.matrix_data_filename+'.fits')

#---------- Line Plot Data --------------------------------------------------------
    # first: provide empty plot data
    y_data            = Array()
    x_data            = Array()
    # second: call the plots
    plot_data         = Instance( ArrayPlotData )
    plot              = Instance( Plot )
    # provide data for all plots
    S0_spec           = Array()
    S1_spec           = Array()
    S2_spec           = Array()
    S3_spec           = Array()
    Ix_spec           = Array()
    Iy_spec           = Array()
    RIp_spec          = Array()
    RIm_spec          = Array()
    GIm_spec          = Array()
    GIp_spec          = Array()
    Tm_spec           = Array()
    Tp_spec           = Array()
    Rot_spec          = Array()
    S3_diff_spec      = Array()
    # Instance the labels
    S0_label          = Instance( PlotLabel )
    S1_label          = Instance( PlotLabel )
    S2_label          = Instance( PlotLabel )
#---------- Matrix 1 --------------------------------------------------------------
    matrix_2_values   = Array()
    matrix_2_data     = Array()
    matrix_2_plot     = MPL2DImageTrait(handler=MPLInitHandler())
    matrix_2_label    = Str()
    matrix_dimension  = Float()
    matrix_dictionary = Dict() # all 2D spectra will be saved here
    matrix_dictionary2= Dict() # all 2D spectra will be saved here but in the right order for saving
    matrix_image_filename   = Str()
    matrix_data_filename = Str()

#---------- mayavi Sphere --------------------------------------------------------------
    scene             = Instance(MlabSceneModel, ())
#---------- Save as pys -----------------------------------------------------------
    # here we need to put in everything which need to be saved
    get_set_items=['x_data', 'S0_spec', 'S1_spec', 'S2_spec','S3_spec','Ix_spec','Iy_spec','RIp_spec','RIm_spec' ,'GIm_spec','GIp_spec','Tm_spec','Tp_spec','Rot_spec','S3_diff_spec','matrix_2_values','matrix_dictionary']

#-------- look and feel------------------------------------------------------------
    colormap    = Enum('PRGn','summer','Spectral','gray', 'prism','hot', 'seismic', 'afmhot', 'gnuplot', 'ocean') # colormaps for a 2D plot

#------- GRAPHICAL USER INTERFACE  ------------------------------------------------
    # the actual gui
    traits_view = View(VGroup(HGroup(
                                     Item('state',style='readonly',width=50),
                                     Item('last_calc_time',style='readonly'),
                                     ),
                              HGroup(Item('element'),
                                     Item('line'),

                                     Item('inputpol'),
                                     Item('detpol')),
                              HGroup(Item('rb85frac',enabled_when='element == "Rb"'),
                                     Item('k40frac', enabled_when='element == "K"'),
                                     Item('k41frac', enabled_when='element == "K"')),
                              HGroup(Item('length'),
                                     Item('bfield'),
                                     ),
                              HGroup(Item('temp1'),
                                     Item('connect_temp'),
                                     Item('temp2',enabled_when='connect_temp != True'),
                                     ),
                              Group(VGroup(
                                     Item('calculate_button',show_label=False),
                                     HGroup(
                                             Item('show_S0_spec'),
                                             Item('show_S1_spec'),
                                             Item('show_S2_spec'),
                                             Item('show_S3_spec'),
                                             Item('show_Ix_spec'),
                                             Item('show_Iy_spec'),
                                             Item('show_RIp_spec'),
                                             Item('show_RIm_spec'),
                                             Item('show_GIp_spec'),
                                             Item('show_GIm_spec'),
                                             Item('show_Tm_spec'),
                                             Item('show_Tp_spec'),
                                             Item('show_Rot_spec'),
                                             Item('normalize'),
                                            # Item('normalize_begin'),
                                            # Item('normalize_end'),
                                             Item('plot_line_width'),),

                                     Item('plot', editor=ComponentEditor(), show_label=False, resizable=True),
                                     HGroup(Item('filename',    springy=True),
                                            Item('save_button', show_label=False),
                                            Item('load_button', show_label=False)
                                     ),
                                     label='Single Spectrum',id='1D Plot',),
                              VGroup(
                                     VGroup(
                                        HGroup(Item('begin'),
                                               Item('end'),
                                               Item('delta'),
                                               Item('bfield_bool',  enabled_when='temp1_bool != True and temp2_bool != True and length_bool != True and inputpol_bool != True  and detpol_bool != True' ),
                                               Item('temp1_bool',   enabled_when='bfield_bool!= True and temp2_bool != True and length_bool != True and inputpol_bool != True  and detpol_bool != True' ),
                                               Item('temp2_bool',   enabled_when='bfield_bool!= True and temp1_bool != True and length_bool != True and inputpol_bool != True  and detpol_bool != True' ),
                                               Item('length_bool',  enabled_when='bfield_bool!= True and temp1_bool != True and temp2_bool  != True and inputpol_bool != True  and detpol_bool != True' ),
                                               Item('inputpol_bool',enabled_when='bfield_bool!= True and temp1_bool != True and temp2_bool  != True and length_bool   != True  and detpol_bool != True' ),
                                               Item('detpol_bool',  enabled_when='bfield_bool!= True and temp1_bool != True and temp2_bool  != True and length_bool   != True  and inputpol_bool != True' ),
                                               ),
                                        HGroup(Item('calculate_matrix_button',show_label=False),
                                               Item('est_calc_time',show_label=True,style='readonly'), ),
                                        HGroup(Item('spectype'),
                                               Item('colormap',),
                                               Item('save_matrix_image_button',show_label=False),
                                               Item('save_matrix_data_button',show_label=False),),
                                      ),

                                     HGroup(
                                             Item('matrix_2_plot',show_label=False, enabled_when='normalize != True',),
                                            ), label='2D Plot', id='2D Plot',

                                   ),
                                   Item('scene', editor=SceneEditor(scene_class=MayaviScene),show_label=False),                                  layout='tabbed',
                                   ),

                              HGroup(Item('frequency_range'),
                                     Item('frequency_delta',),
                                     Item('frequency_shift',),
                                     Item('frequency_broadening',),

                                     ),

                              ),
                       title='ElecSus GUI', buttons=[], resizable=True
                       )

#-------  INIT  ------------------------------------------------
    def __init__(self,  **kwargs): #**kwargs means you can overload the class with as much as you like or need
        super(Gui, self).__init__(**kwargs)


        # we need create the plots here, methods are found below
        self._create_plot() # Create line plot
        self.on_trait_change(self._update_index,       'x_data',             dispatch='ui')
        self.on_trait_change(self._update_value,       'y_data',             dispatch='ui')
        self.on_trait_change(self._update_S0_plot,     'S0_spec,normalize,normalize_begin,normalize_end,plot_line_width',  dispatch='ui')
        self.on_trait_change(self._update_S1_plot,     'S1_spec,normalize,normalize_begin,normalize_end,plot_line_width',  dispatch='ui')
        self.on_trait_change(self._update_S2_plot,     'S2_spec,normalize,normalize_begin,normalize_end,plot_line_width',  dispatch='ui')
        self.on_trait_change(self._update_S3_plot,     'S3_spec,normalize,normalize_begin,normalize_end,plot_line_width',  dispatch='ui')
        self.on_trait_change(self._update_Ix_plot,     'Ix_spec,normalize,normalize_begin,normalize_end,plot_line_width',  dispatch='ui')
        self.on_trait_change(self._update_Iy_plot,     'Iy_spec,normalize,normalize_begin,normalize_end,plot_line_width' , dispatch='ui')
        self.on_trait_change(self._update_RIp_plot,    'RIp_spec,normalize,normalize_begin,normalize_end,plot_line_width', dispatch='ui')
        self.on_trait_change(self._update_RIm_plot,    'RIm_spec,normalize,normalize_begin,normalize_end,plot_line_width', dispatch='ui')
        self.on_trait_change(self._update_GIm_plot,    'GIm_spec,normalize,normalize_begin,normalize_end,plot_line_width', dispatch='ui')
        self.on_trait_change(self._update_GIp_plot,    'GIp_spec,normalize,normalize_begin,normalize_end,plot_line_width', dispatch='ui')
        self.on_trait_change(self._update_Tm_plot,     'Tm_spec,normalize,normalize_begin,normalize_end,plot_line_width',  dispatch='ui')
        self.on_trait_change(self._update_Tp_plot,     'Tp_spec,normalize,normalize_begin,normalize_end,plot_line_width',  dispatch='ui')
        self.on_trait_change(self._update_Rot_plot,    'Rot_spec,normalize,normalize_begin,normalize_end,plot_line_width', dispatch='ui')
        self.on_trait_change(self.estimate_calc_time,  'begin,end, delta,',  dispatch='ui')
        self._create_sphere() # Create sphere plot
        self.on_trait_change(self._update_plotsphere,  'Rot_spec',  dispatch='ui')
        self._create_matrix_2_plot() # create matrix plot        
        self.on_trait_change(self._update_matrix_2_data_value,  'y_data,matrix_2_values',     dispatch='ui')
        self.on_trait_change(self._update_matrix_2_data_index,  'matrix_dimension',    dispatch='ui')
        self.on_trait_change(self._spectype_changed,'spectype',    dispatch='ui')


        # default values
        self.show_S0_spec = True
        self.element = 'Cs'
        self.colormap = 'hot'


    def _calculate(self):


        try:
            startTime           = timing()
            self.state          = 'calculating'
            self.x_data         = self.generate_frequency()
            self.y_data         = self.return_spectrum(self.x_data)
            self.state          = 'done'
            self.last_calc_time = 1000*(timing()-startTime)
            file_prefix         = self.return_file_prefix()            
        except:
            self.state = 'error'
        finally:
            self.state = 'done'

    def _calculate_all(self):
        """calculates all spetra and sets their values"""

        try:
            startTime = timing()
            self.state='calculating'
            self.x_data                        = self.generate_frequency()
            S0,S1,S2,S3,Ix,Iy,RIp,RIm,GIm,GIp,Tm,Tp,Rot  = self.return_spectrum_all(self.x_data)
            self.S0_spec,self.S1_spec,self.S2_spec,self.S3_spec,self.Ix_spec,self.Iy_spec,self.RIp_spec,self.RIm_spec,self.GIm_spec,self.GIp_spec,self.Tm_spec,self.Tp_spec,self.Rot_spec  = S0,S1,S2,S3,Ix,Iy,RIp,RIm,GIm,GIp,Tm,Tp,Rot
            self.set_y_data() # thats for the matrix plot data setting
            self.State                         = 'done'
            self.last_calc_time                = 1000*(timing()-startTime)
            self.estimate_calc_time()
            file_prefix                        = self.return_file_prefix()
            self.filename                      = file_prefix+str(self.element).lower()+'_'+str(self.line).lower()+'_'+str(self.bfield)+'_'+str(self.temp1)+'.fits'
        except:
            self.state = 'error'
            print 'some error occured in calculation, check: _calculate_all()'
        finally:

            self.state = 'done'
    def _calculate_matrix(self):

        self.matrix_2_values = self._matrix_2_values_default() # clear matrix data
        my_reset_values=self.return_reset_values() # store the values to reset them at the end
        try:
            startTime             = timing()
            self.state            = 'calculating matrix data'
            self.x_data           = self.generate_frequency()
            myfield               = np.arange(self.begin, self.end, self.delta)
            self.matrix_dimension = len(myfield) # adjust matrix lines here
    
    
    
    
            for i,value in enumerate(myfield):
                self.state = 'calculating'
                if i==0:
                    self.old_data = [self.return_hstacked(self.return_spec_data())]#save all old calculated data here
                if self.bfield_bool  :
                  self.bfield         = value
                  self.matrix_2_label = 'B-Field[Gauss]'
                if self.temp1_bool   :
                  self.temp1          = value
                  self.matrix_2_label = 'Temperature 1 ['+self.degree+'C]'
                if self.temp2_bool   :
                  self.temp2          = value
                  self.matrix_2_label = 'Temperature 2 ['+self.degree+'C]'
                if self.length_bool  :
                  self.length         = value
                  self.matrix_2_label = 'Cell Length [mm]'
                if self.inputpol_bool:
                  self.inputpol       = value
                  self.matrix_2_label = 'Input polarization ['+self.degree+']'
                if self.detpol_bool  :
                  self.detpol         = value
                  self.matrix_2_label = 'Output polarization [%]'
                self.new_data = [self.return_hstacked(self.return_spec_data())]
                self.old_data = np.vstack((self.new_data,self.old_data))
                self.matrix_2_values  = np.vstack( (self.y_data, self.matrix_2_values[:-1,:]) )
    
            cutted_data = self.old_data[:-1,:,:]#one entry has to be removed
            self.matrix_dictionary = self.create_2D_spectra_dict(cutted_data)#for plotting
            self.matrix_dictionary2= self.create_2D_spectra_dict(cutted_data[::-1,:,:])#reverse it for saving
    
            self.state                = 'done'
            self.last_calc_time       = 1000*(timing()-startTime)
            file_prefix               = self.return_file_prefix()
            self.filename             = file_prefix+str(self.element).lower()+'_'+str(self.line).lower()+'_'+str(self.spectype).lower()+'_'+str(self.bfield)+'_'+str(self.temp1)+'_matrix.txt'
            self.matrix_image_filename= file_prefix+str(self.element).lower()+'_'+str(self.line).lower()+'_'+str(self.spectype).lower()+'_'+self.matrix_2_label+'_'+str(self.begin)+'_'+str(self.end)
            self.matrix_data_filename = file_prefix+str(self.element).lower()+'_'+str(self.line).lower()+'_'+self.matrix_2_label+'_'+str(self.begin)+'_'+str(self.end)
        except:
            self.state = 'error'
            print 'uups....something went wrong in _calculate_matrix()'
            self.reset_values(my_reset_values) # reset iterated values to start value
        finally:
            self.state = 'done'
            self.reset_values(my_reset_values)

    def set_y_data(self):
        """method to get the desired y-data for matrix plot"""
        self.y_data = self.get_spec_type(self.spectype)


    def create_dict_matrix(self,key_list,value_list):
        """
        creates a dictionary from matrix calculation
        key_list = spectypes (1D list)
        values_list = 3D array (x,y,z)
        """
        mydict = {}
        i=0
        for key in key_list:
          mydict.update({key:value_list[::,i,::]})
          i+=1
        return mydict

    def return_spec_data(self):
        """returns all specs as list"""
        return [self.S0_spec,self.S1_spec,self.S2_spec,self.S3_spec,self.Ix_spec,self.Iy_spec,self.RIp_spec,self.RIm_spec,self.GIm_spec,self.GIp_spec,self.Tm_spec,self.Tp_spec,self.Rot_spec]

    def return_hstacked(self,mylist):
        stacked_list=[]
        for entry in mylist:
            stacked_list.append(entry)
        return stacked_list



    def create_2D_spectra_dict(self,all_matrix_data):
        value_list = all_matrix_data
        key_list = ['S0','S1','S2','S3','Ix','Iy','RI+','RI-','GI-','GI+','T-','T+','Rot']
        k=self.create_dict_matrix(key_list,value_list)
        return k

    def _spectype_changed(self):
        s = self.spectype
        data=self.matrix_dictionary[s]
        self.matrix_2_values = data




    def get_spec_type(self,spectype):
        """method to get the selected spectypes"""
        if spectype == 'S0':
            s = self.S0_spec
        if spectype == 'S1':
            s = self.S1_spec
        if spectype == 'S2':
            s = self.S2_spec
        if spectype == 'S3':
            s = self.S3_spec
        if spectype == 'Ix':
            s = self.Ix_spec
        if spectype == 'Iy':
            s = self.Iy_spec
        if spectype == 'RI+':
            s = self.RIp_spec
        if spectype == 'RI-':
            s = self.RIm_spec
        if spectype == 'GI-':
            s = self.GIm_spec
        if spectype == 'GI+':
            s = self.GIp_spec
        if spectype == 'T-':
            s = self.Tm_spec
        if spectype == 'T+':
            s = self.Tp_spec
        if spectype == 'Rot':
            s = self.Rot_spec
        return s

    def estimate_calc_time(self):
        """estimate calculation time for matrix plot\n
        sets (global): est_calc_time
        """
        myfield   = np.arange(self.begin, self.end, self.delta) # it depends on length of the field
        seconds = len(myfield)*self.last_calc_time/1000 # our seconds in total,based on previous calculation
        m, s = divmod(seconds, 60) # convert to minutes and second
        h, m = divmod(m, 60) # convert to hours and minutes
        k= "%d:%02d:%02d" % (h, m, s) # convert to string in a proper format
        self.est_calc_time = k # sets the estimated time

    def CalculateStringTime(self, seconds):
        """Calculates current time and returns time as string as\n
        return value: [d.day-1, d.hour, d.minute, d.second]
        """
        try:
            d = datetime(1,1,1) + seconds
            k = "%d:%d:%d:%d" % (d.day-1, d.hour, d.minute, d.second)
        except:
            k='0:0:0:0'
        return k

    #################################################################
    # Helper Methods
    #################################################################

    def generate_frequency(self):
        """simple frequency array generation"""
        mesh = np.arange(self.frequency_begin, self.frequency_end, self.frequency_delta/1000.0)
        return mesh

    def return_spectrum(self, frequency, spectype):
        """here you will get only one single spectrum back, \n
        we prefer to use return_sectrum_all, because it doesnt cost much more computing power"""
        if self.connect_temp:
          self.temp2 = self.temp1
        return spectra.spectrum(frequency,self.element,spectype,self.bfield,self.temp1,self.length,self.rb85frac,self.temp2,self.inputpol,self.detpol,self.frequency_shift,self.frequency_broadening,self.connect_temp,self.line,self.frequency_delta,self.k40frac,self.k41frac)

    def return_spectrum_all(self, frequency):
        """here you will get all spectra back"""        
        if self.connect_temp:
          self.temp2 = self.temp1
        my_dict = self.create_dictionary()
        #return spectra.get_spectra(frequency,self.element,self.bfield,self.temp1,self.length,self.rb85frac,self.temp2,self.inputpol,self.detpol,self.frequency_shift,self.frequency_broadening,self.connect_temp,self.line,self.frequency_delta,self.k40frac,self.k41frac)
        return spectra.get_spectra(frequency*1000.,my_dict,outputs='All')
    def create_dictionary(self):
        my_dict = {}
        key_list = ['Elem','Dline','Bfield','T','GammaBuf','Shift','DoppTemp','Constrain','rb85frac','K40frac','K41frac','lcell','theta0','Pol']
        entries  = [self.element, self.line,self.bfield,self.temp1,self.frequency_broadening, self.frequency_shift,self.temp2,self.connect_temp,self.rb85frac, self.k40frac,self.k41frac, self.length/1000.,self.inputpol,self.detpol]
        for key, entry in zip(key_list,entries):
            my_dict[key] = entry
        return my_dict

    def reset_values(self,reset_values):
        """resets these values after matrix run:\n
        [bfield,temp1,temp2,length,inputpol,detpol]\n
        they need to be stored beforehand
        """
        self.bfield,self.temp1,self.temp2,self.length,self.inputpol,self.detpol = reset_values

    def return_reset_values(self):
        """returns the current parameter set: \n
        [bfield,temp1,temp2,length,inputpol,detpol]
        """
        reset_value = self.bfield,self.temp1,self.temp2,self.length,self.inputpol,self.detpol
        return reset_value
    def return_file_prefix(self):
        """returns the the file prefix for different OS"""
        if 'Windows' in platform.platform():
          file_prefix = 'c:\\'
        else:
          file_prefix = '/tmp/'
        return file_prefix

    def return_normalized_data(self,spec):
        """returns one given spec as normalized spectra from user defined range:\n
        uses (global) normalize_end AND (global normalize_begin)
        """
        if self.normalize:
          if self.normalize_end < self.normalize_begin: # correct for switched normalized factors
            ne,nb = self.normalize_end,self.normalize_begin
            self.normalize_begin, self.normalize_end = ne,nb
          spec_zero = spec-min(spec) # set minimum to zero
          to_one = (spec_zero)/max(spec_zero) # normalize max to 1
          normalized = (to_one*(abs(self.normalize_end)-self.normalize_begin))+self.normalize_begin # scale to user defined range
          return normalized
        else:
         return spec



    #################################################################
    # LINE PLOT DEFINITIONS
    #################################################################
    def _create_plot(self):
        plot_data = ArrayPlotData(x_data=np.array(()), y_data=np.array(()), S0_spec=np.array(()), S1_spec=np.array(()),Ix_spec=np.array(()))
        plot = Plot(plot_data, padding=8, padding_left=64, padding_bottom=64)
        plot.index_axis.title = 'Frequency Detuning [GHz]'
        plot.value_axis.title = 'Spectrum'
        plot.tools.append(SaveTool(plot))
        zoom=simple_zoom.SimpleZoom(plot) # now we can zoom in with mouse wheel
        plot.overlays.append(zoom)
        #------labels-------------
        S0_label  = PlotLabel(text='', hjustify='left', vjustify='bottom', position=[64,128])
        plot.overlays.append(S0_label)
        #------make all global -------------
        self.S0_label         = S0_label
        self.plot_data        = plot_data
        # self.plot_data.set_data('S0_spec', self.S0_spec)
        self.plot             = plot

    def _update_index(self, new):
        self.plot_data.set_data('x_data', new)

    def _update_value(self):
      pass

#-------  S0  ------------
    def _show_S0_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','S0_spec'), style='line', line_width=self.plot_line_width, color='red', name='S0_spec')

        elif not self.show_S0_spec:
            plot.delplot('S0_spec')
        plot.request_redraw()

    def _update_S0_plot(self):
        data = self.return_normalized_data(self.S0_spec)
        self.plot_data.set_data('S0_spec', data)
#-------  S1  ------------
    def _show_S1_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','S1_spec'), style='line', line_width=self.plot_line_width, color='blue', name='S1_spec')
        elif not self.show_S1_spec:
            plot.delplot('S1_spec')
        plot.request_redraw()

    def _update_S1_plot(self):
        data = self.return_normalized_data(self.S1_spec)
        self.plot_data.set_data('S1_spec', data)

#-------  S2  ------------
    def _show_S2_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','S2_spec'), style='line', line_width=self.plot_line_width, color='black', name='S2_spec')
        elif not self.show_S2_spec:
            plot.delplot('S2_spec')
        plot.request_redraw()

    def _update_S2_plot(self):
        data = self.return_normalized_data(self.S2_spec)
        self.plot_data.set_data('S2_spec', data)

#-------  S3  ------------
    def _show_S3_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','S3_spec'), style='line', line_width=self.plot_line_width, color='grey', name='S3_spec')
        elif not self.show_S3_spec:
            plot.delplot('S3_spec')
        plot.request_redraw()

    def _update_S3_plot(self):
        data = self.return_normalized_data(self.S3_spec)
        self.plot_data.set_data('S3_spec', data)

#-------  Ix  ------------
    def _show_Ix_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','Ix_spec'), style='line', line_width=self.plot_line_width, color='green', name='Ix_spec')
        elif not self.show_Ix_spec:
            plot.delplot('Ix_spec')
        plot.request_redraw()

    def _update_Ix_plot(self):
        data = self.return_normalized_data(self.Ix_spec)
        self.plot_data.set_data('Ix_spec', data)

#-------  Iy  ------------
    def _show_Iy_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','Iy_spec'), style='line', line_width=self.plot_line_width, color='orange', name='Iy_spec')
        elif not self.show_Iy_spec:
            plot.delplot('Iy_spec')
        plot.request_redraw()

    def _update_Iy_plot(self):
        data = self.return_normalized_data(self.Iy_spec)
        self.plot_data.set_data('Iy_spec', data)


#-------  RIp  ------------
    def _show_RIp_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','RIp_spec'), style='line', line_width=self.plot_line_width, color='yellowgreen', name='RIp_spec')
        elif not self.show_RIp_spec:
            plot.delplot('RIp_spec')
        plot.request_redraw()

    def _update_RIp_plot(self):
        data = self.return_normalized_data(self.RIp_spec)
        self.plot_data.set_data('RIp_spec', data)
#-------  RIm  ------------
    def _show_RIm_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','RIm_spec'), style='line', line_width=self.plot_line_width, color='black', name='RIm_spec')
        elif not self.show_RIm_spec:
            plot.delplot('RIm_spec')
        plot.request_redraw()

    def _update_RIm_plot(self):
        data = self.return_normalized_data(self.RIm_spec)
        self.plot_data.set_data('RIm_spec', data)
#-------  GIm  ------------
    def _show_GIm_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','GIm_spec'), style='line', line_width=self.plot_line_width, color='blue', name='GIm_spec')
        elif not self.show_GIm_spec:
            plot.delplot('GIm_spec')
        plot.request_redraw()

    def _update_GIm_plot(self):
        data = self.return_normalized_data(self.GIm_spec)
        self.plot_data.set_data('GIm_spec', data)
#-------  GIp  ------------
    def _show_GIp_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','GIp_spec'), style='line', line_width=self.plot_line_width, color='red', name='GIp_spec')
        elif not self.show_GIp_spec:
            plot.delplot('GIp_spec')
        plot.request_redraw()

    def _update_GIp_plot(self):
        data = self.return_normalized_data(self.GIp_spec)
        self.plot_data.set_data('GIp_spec', data)
#-------  Tm  ------------
    def _show_Tm_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','Tm_spec'), style='line', line_width=self.plot_line_width, color='cyan', name='Tm_spec')
        elif not self.show_Tm_spec:
            plot.delplot('Tm_spec')
        plot.request_redraw()

    def _update_Tm_plot(self):
        data = self.return_normalized_data(self.Tm_spec)
        self.plot_data.set_data('Tm_spec', data)
#-------  Tp  ------------
    def _show_Tp_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','Tp_spec'), style='line', line_width=self.plot_line_width, color='darkgreen', name='Tp_spec')
        elif not self.show_Tp_spec:
            plot.delplot('Tp_spec')
        plot.request_redraw()

    def _update_Tp_plot(self):
        data = self.return_normalized_data(self.Tp_spec)
        self.plot_data.set_data('Tp_spec', data)
#-------  Rot  ------------
    def _show_Rot_spec_changed(self,new):
        plot = self.plot
        if new:
            plot.plot(('x_data','Rot_spec'), style='line', line_width=self.plot_line_width, color='darkblue', name='Rot_spec')
        elif not self.show_Rot_spec:
            plot.delplot('Rot_spec')
        plot.request_redraw()

    def _update_Rot_plot(self):
        data = self.return_normalized_data(self.Rot_spec)
        self.plot_data.set_data('Rot_spec', data)

    #################################################################
    # Save definitions, save_all not used right now
    #################################################################

    def save_plot(self, filename):
        save_figure(self.plot, filename)
    def save_plot(self, filename):
        save_figure(self.plot, filename)
    def save_all(self, filename):
        self.plot.save(filename+'_plot01.png')
        self.save(filename+'.pys')
        self.save(filename+'-ACSII.pys')

    #################################################################
    # Matrix PLOT DEFINITIONS
    #################################################################

    def _create_matrix_2_plot(self):
        frequency = self.generate_frequency()
        matrix_2_data = np.zeros((2,2))
        self.matrix_2_plot._create_image()
        self.matrix_2_plot.enable_colorbar = False
        self.matrix_2_plot.cmap = self.colormap
        self.matrix_2_data = matrix_2_data
        
    def setup_mpl_events(self):
        self.matrix_2_plot.figure.canvas.mpl_connect('key_press_event', self.key_event)        

    def _matrix_2_values_default(self):
        self.matrix_dimension = len(np.arange(self.begin,self.end+self.delta, self.delta))
        frequency_dimension = len(self.generate_frequency())
        return np.zeros( (self.matrix_dimension, frequency_dimension) )

    def _update_matrix_2_data_value(self):
        self.matrix_2_plot.image = self.matrix_2_values
        self.matrix_2_plot.axes.set_ylabel(self.matrix_2_label)
        self.matrix_2_plot.axes.set_xlabel('Frequency [GHz]')
        self.matrix_2_plot.update_image()

    def _update_matrix_2_data_index(self):
        self.matrix_2_plot.set_extent(self.frequency_begin,self.frequency_end,self.begin,self.end)
        self.matrix_2_values = self.matrix_2_values[:self.matrix_dimension]

    def _spectra_changed(self):
      pass

    def _colormap_changed(self, new):
      """you want more colormaps? visit: http://matplotlib.org/examples/color/colormaps_reference.html"""
      self.matrix_2_plot.cmap = self.colormap
      self.matrix_2_plot.update_figurecmap()

    #################################################################
    # SPHERE PLOT DEFINITIONS
    #################################################################

    def _create_sphere(self):
        sphere = self.scene.mlab.points3d(0, 0, 0, scale_mode='none',
                                        scale_factor=2,
                                        color=(0.67, 0.77, 0.93),
                                        resolution=50,
                                        opacity=0.2,
                                        name='poin')

        sphere.actor.property.specular = 0.45
        sphere.actor.property.specular_power = 5
        sphere.actor.property.backface_culling = True

        theta = np.linspace(0, 2 * np.pi, 100)
        for angle in np.linspace(-np.pi, np.pi, 9):
            xlat = np.cos(theta) * np.cos(angle)
            ylat = np.sin(theta) * np.cos(angle)
            zlat = np.ones_like(theta) * np.sin(angle)
            xlon = np.sin(angle) * np.sin(theta)
            ylon = np.cos(angle) * np.sin(theta)
            zlon = np.cos(theta)
            self.scene.mlab.plot3d(
                xlat, ylat, zlat,
                color=(0.67, 0.77, 1),
                opacity=self.frame_alpha, tube_radius=self.frame_radius)
            self.scene.mlab.plot3d(
                xlon, ylon, zlon,
                color=(0.67, 0.77, 1),
                opacity=self.frame_alpha, tube_radius=self.frame_radius)

        axis = np.linspace(-1.2, 1.2, 10)
        other = np.zeros_like(axis)
        self.scene.mlab.plot3d(
            axis, other, other,
            color=(1, 0.77, 1),
            tube_radius=self.axes_radius, opacity=self.axes_alpha)
        self.scene.mlab.plot3d(
            other, axis, other,
            color=(1, 0.77, 1),
            tube_radius=self.axes_radius, opacity=self.axes_alpha)
        self.scene.mlab.plot3d(
            other, other, axis,
            color=(1, 0.77, 1),
            tube_radius=self.axes_radius, opacity=self.axes_alpha)

        mesh = np.arange(self.frequency_begin, self.frequency_end, self.frequency_delta/1000.)
        #ss=spectra.all_spec(mesh,self.element,self.bfield,self.temp1,self.length,self.rb85frac,self.temp2,self.inputpol,self.detpol,self.frequency_shift,self.frequency_broadening,True,self.line,self.frequency_delta,self.k40frac,self.k41frac)
        my_dict=self.create_dictionary()
        ss=spectra.get_spectra(mesh*1000.,my_dict,outputs='All')
        t, x, z, y = ss[0:4]#self.S0_spec,self.S1_spec,self.S2_spec,self.S3_spec
        self.plotsphere = self.scene.mlab.plot3d(x, y, z, t,tube_radius=0.01)

    def _update_plotsphere(self):
        t, x, z, y = self.S0_spec,self.S1_spec,self.S2_spec,self.S3_spec
        self.plotsphere.mlab_source.reset(x=x, y=y, z=z, scalars=t)
        self.plotsphere.mlab_source.set(x=x, y=y, z=z, scalars=t) # we need to do set it once again to have the colorbar



if __name__=='__main__':
    elecsus_gui = Gui()
    elecsus_gui.configure_traits()