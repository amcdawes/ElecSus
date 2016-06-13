"""
Copyleft (C) 2015 Ilja Gerhardt <ilja@quantumlah.org>
                  Matthias Widmann <m.widmann@physik.uni-stuttgart.de>
                  Matthias Niethammer <m.niethammer@physik.uni-stuttgart.de>

Much credit to Gael Varoquaux for his 
[tutorial on using matplotlib within traitsui]
(http://docs.enthought.com/traitsui/tutorials/traits_ui_scientific_app.html) 
from which the following Editor classes are derived.

"""
import wx 
import matplotlib, sys
matplotlib.use('WXAgg')
matplotlib.rcParams.update({'figure.autolayout': True})
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from matplotlib.axes import Axes
from matplotlib.widgets import AxesWidget
import matplotlib.pyplot as plt
from matplotlib import cm

 
from traits.api import Instance
from traitsui.wx.editor import Editor
from traitsui.wx.basic_editor_factory import BasicEditorFactory
import numpy as np
 
from traits.api import HasTraits, Array, CInt, Int, CFloat, Float, Str, on_trait_change, Range, Enum, List, Dict, Bool, Trait
from traitsui.api import View, Item, Handler, HGroup, VGroup, StatusItem, TextEditor
from traitsui.ui_info import UIInfo
  
from utility import save_file # to show up a gui for saving
import data_toolbox # to save matrix data as pys,fits, txt etc...

# TODO:  someday come back and play more with resizing, including:
#        - maintaining proportion of image panel
#        - minimum size for each element enforced (i.e. no clipping of UI and no resizing of plot figure) so window can't be made tiny
 
# TODO:  someday come back and play with how gui closes itself out. (if run from python environment the window does not close itself, although control is returned to the python prompt and it can be relaunched just fine again.)  I suspect this is something left dangling by wx or matplotlib, so probably answer is to clean up and delete a bunch of those object/event connections.
 
def _clear_ticks_and_frame_from_axes(ax=None):
    if ax is None:
        ax = plt.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for spine in ax.spines.itervalues():
        spine.set_visible(False)
 
def _set_ticks_and_frame_from_axes(ax=None):
    if ax is None:
        ax = plt.gca()
    ax.xaxis.set_visible(True)
    ax.yaxis.set_visible(True)
    for spine in ax.spines.itervalues():
        spine.set_visible(True)
 

 
 
class _ResizingMPLFigureEditor(Editor):
 
    scrollable = True
    panel = Instance(wx._windows.Panel)
    canvas = Instance(FigureCanvas)
 
    def init(self, parent): 
        self.control = self._create_canvas(parent)

 
    def update_editor(self):
        pass
 
    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        self.panel = panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        self.canvas = FigureCanvas(panel, -1, self.value.figure)
        sizer.Add(self.canvas, 1.0, wx.LEFT | wx.TOP | wx.GROW)
        return panel
 
 
class ResizingMPLFigureEditor(BasicEditorFactory):
 
    klass = _ResizingMPLFigureEditor
 
 
class _NonResizingMPLFigureEditor(Editor):
 
    scrollable = False
    panel = Instance(wx._windows.Panel)
    canvas = Instance(FigureCanvas)
 
    def init(self, parent):
        self.control = self._create_canvas(parent)
 
    def update_editor(self):
        pass
 
    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        self.panel = panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        self.canvas = FigureCanvas(panel, -1, self.value.figure)
        sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        return panel
 
 
class NonResizingMPLFigureEditor(BasicEditorFactory):
 
    klass = _NonResizingMPLFigureEditor
 
 


class MPLInitHandler(Handler):
 
    ui_info = Instance(UIInfo)
 
    def init(self, info):
        """
        This method gets called after the controls have all been
        created but before they are displayed.
        """
        self.ui_info = info
        self.ui_info.object.figure.canvas.mpl_connect('key_press_event', self.ui_info.object.key_event)
        return True
 

class MPLArrayData(HasTraits):

    data = Dict()


    ##ATTACH TO DATA CHANGED EVENTS OF THE STUFF?

    def set_data(self,**kwargs):
        for key,value in kwargs.iteritems():
            try:
                value = np.asarray(value)
                self.data[key] = value
            except:
                raise TypeError("Cannot convert given data to numpy array.")

    def clear_data(self):
        self.data = {}

    def delete_data(self, name):
        try:
            del self.data[name]
        except KeyError:
            print "Data with", name, " not stored in this MPLArrayData instance."
    

class MPL1DPlot(HasTraits):
 

    figure = Instance(Figure, (),editor=ResizingMPLFigureEditor(), show_label=False, width=800, height=200)

 
    def __init__(self):
        super(MPL1DPlot, self).__init__()
        self.figure = Figure(facecolor='white')
        self.axes = self.figure.add_subplot(111)


    def plot(self,x,y,*args,**kwargs):
        """Convenience wrapper for self.axes.plot"""
        self.axes.plot(x,y,*args,**kwargs)

    def configure(self):
        """Overwrite this function for genearl plot configuration.
        Figure is available as self.figure, axes as self.axes.
        This function is called on every redraw event."""
        pass

    def autoscale(self):
        self.axes.relim()
        self.axes.autoscale_view()        

    def draw(self):
        self.figure.canvas.draw()

    def clear(self):
        for axis in self.figure.axes:
            axis.clear()
 
    def default_traits_view(self):
        return View(Item('figure',  editor=ResizingMPLFigureEditor(), show_label=False,
                                              width=400, height=200))
 
    def _fresh_plot_data(self):
        x = np.linspace(0,1,50)
        return x, np.zeros_like(x)
 
    @on_trait_change('plot_num_points, plot_coeff+, plot_gain')
    def update_plot_data(self):
        self.axes.lines = []
        x, perfect_y, observed_y = self._fresh_plot_data()
        self.axes.plot(x, perfect_y, 'bo', markersize=3.0)
        self.axes.plot(x, observed_y, 'm+', markersize=10.0)
        self.axes.relim()
        self.axes.autoscale_view()
        self.figure.canvas.draw()
 
        
 
class MPL2DImage(HasTraits):
 
    image = Array()
    figure = Instance(Figure, (), editor=ResizingMPLFigureEditor(), show_label=False)#, width=800, height=200)
    axes = Instance(Axes)
    axesimage = Instance(AxesImage)
    xsize = Int(2)
    ysize = Int(2)
    x_begin = Float(0.)
    x_end =Float(1.)
    y_begin = Float(0.)
    y_end=Float(1.)
    cmap = Str('hot')
    clim_lower = Trait( 'auto', Str('auto'), Float(0.), desc='Low Limit of image plot', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float))
    clim_upper =  Trait( 'auto', Str('auto'), Float(10000.), desc='High Limit of image plot', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float))
    enable_colorbar = Bool(True)
    cbar_orientation = Str('vertical')

    def __init__(self):
        super(MPL2DImage, self).__init__()
        self.figure = Figure(facecolor='white')
        self.image = self._fresh_image()
        self._create_image()


    def _create_image(self):
        self.axes = self.figure.add_subplot(111,frameon=False)
        
        self.axesimage = self.axes.imshow(self.image, extent=[self.x_begin,self.x_end,self.y_begin,self.y_end], cmap=self.cmap, interpolation='nearest',aspect='auto')   

    def default_traits_view(self):
        return View(HGroup(VGroup(Item('figure', editor=ResizingMPLFigureEditor(), show_label=False),
                                  show_border=False)),
                    resizable=True,
                    handler=MPLInitHandler())
 
    def _fresh_image(self):    
        return (np.outer(np.ones(self.ysize), np.arange(self.xsize) /
                (self.xsize - 1.0)))
 
    def update_image(self):
        self.axesimage.set_data(self.image)
        self.set_color_limits()
        if self.enable_colorbar:
            if not hasattr(self, 'colorbar'):
                self.colorbar = self.figure.colorbar(self.axesimage)
        else:
            try:
                self.figure.delaxes(self.figure.axes[1])
                self.figure.subplots_adjust(right=0.9)           
                del self.colorbar
            except:
                pass
        self.redraw()

    def set_color_limits(self):
        if str(self.clim_lower) == "auto":
            cmin = np.min(self.image)
        else:
            cmin = self.clim_lower
        if str(self.clim_upper) == "auto":
            cmax = np.max(self.image)
        else:
            cmax = self.clim_upper
        self.axesimage.set_clim(cmin, cmax)

    def set_extent(self,x_begin,x_end,y_begin,y_end):
        #print "Setting extent to", x_begin, ",",x_end, ",",y_begin, ",",y_end
        self.x_begin = x_begin
        self.x_end   = y_end
        self.y_begin = y_begin
        self.y_end   = y_end
        self.axesimage.set_extent([x_begin,x_end,y_begin,y_end])
 
    def update_figurecmap(self):
        self.axesimage.set_cmap(self.cmap)
        self.redraw()

    def update_figurecolorbar(self):
        self.figure.colorbar(self.figureaxesimage)
        #self.figure.canvas.draw()#does not work!?

    def key_event(self,event):
        """the event """
        sys.stdout.flush()
        #print 'key: ', event.key
        if event.key != None and event.key == 'ctrl+s':
            filename = save_file(title='Save plot as png or pdf.') #uses utility.save_file for pop up gui
            if filename:
                self.figure.savefig(filename) # now save the plot
    
    def save_matrix_image(self,afilename):
        """similar to key event, just gives the opportunity to use eg. a button"""
        filename = save_file(title='Save plot as png or pdf.', afilename=afilename) #uses utility.save_file for pop up gui
        #filename.dialog_box.filename = 'c:\\'
        if filename:
            self.figure.savefig(filename) # now save the plot 

    def save_matrix_data(self,dict,afilename):
        """pass over a dictionary, which will be saved"""
        filename = save_file(title='Save ALL matrix data as pys, fits, ascii or txt.',afilename=afilename) #uses utility.save_file for pop up gui
        if filename:
            data_toolbox.writeDictToFile(dict,filename)

    def setup_mpl_events(self):
        """setup event listeners, want more? Go to: http://matplotlib.org/users/event_handling.html"""
        #print 'setting up mpl events'
        self.figure.canvas.mpl_connect('key_press_event', self.key_event)


    def redraw(self):
        """redraws the plot, purpose is to avoid error message and to reduce code"""
        try:
            self.figure.canvas.draw()#avoid error message
        except:
            pass#do nothing


   

MPL1DPlotTrait = Trait(MPL1DPlot(),editor=ResizingMPLFigureEditor())
MPL2DImageTrait = Trait(MPL2DImage(),editor=ResizingMPLFigureEditor(),handler=MPLInitHandler())
 
if __name__ == "__main__":
    test=Test()
    test.configure_traits()



