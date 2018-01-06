import numpy as _np
import matplotlib.pyplot as _plt
#import _processData as _process
#from matplotlib.colors import LinearSegmentedColormap as _lsc


# color sequence for plotting.  add more colors if you need more than 7.  
_cSequence=['b', 'r', 'g', 'k', 'm', 'c', 'y']


class plot:
    """
    Class structure that contains a single plot window.  It can support 
    multiple "plots" (data arrays) on it.  
    
    Attributes
    ----------
    title : str
        Title to be displayed at above the figure.  an empty str prints 
        nothing
    subtitle : str
        Subtitle to be displayed at the top (but within) of the figure.  an
        empty str prints nothing
    xLabel : str
        the x-label to be displayed.  an empty str prints nothing
    yLabel : str
        the y-label to be displayed.  an empty str prints nothing
    zLabel : str
        the y-label to be displayed.  an empty str prints nothing.  this is 
        only displayed if there is a z-axis of the data
    xData : list (of numpy.ndarray)
        the x-data of the plot.  must be stored as a list of arrays.  if 
        multiple plots are to be in this same figure, an array for each 
        must be included
        e.g. xData = [np.arange(0,10)] 
             or 
             xData = [np.arange(0,10),np.arange(0,10),np.arange(0,10) ]
    yData : list (of numpy.ndarray)
        the y-data of the plot.  must be stored as a list of arrays.  if 
        multiple plots are to be in this same figure, an array for each 
        must be included
        e.g. yData = [np.arange(0,10)] 
             or 
             yData = [np.arange(0,10),np.arange(0,10),np.arange(0,10) ]
    zData : list (of numpy.ndarray)
        the z-data of the plot.  applicable for contour and scatter plots. 
        must be stored as a list of arrays.  multiple entries isn't supported 
        because having multiple 3D plots on a single plot doesn't make sense.  
        e.g. zData = [np.arange(0,10)] 
        Leave empty if not in use.
    yLegendLabel : list (of numpy.ndarray)
        The legend label for each data array.  
        e.g. yLegendLabel = ['smoothed data','raw data'] 
    linestyle : list (or str)
        List of strings, one for each data array.  Leave empty is you want all
        data to be plotted as standard lines.  Alternatively, specify each 
        for something different.
        E.g. linestyles = ['-',':','']
    marker : list (of str)
        marker to be plotted for its associated data array.  default is none.
    xLim : ilst (of two floats)
        x-limit of the plot
    yLim : ilst (of two floats)
        y-limit of the plot
    legendLoc : str
        Location on on the plot for the legend.  'upper right' is default.  
        Other examples include 'lower left', 'center right', etc.
    showGrid : bool
        True - creates grid on the plot.  True is default
    yErData : list (of numpy.ndarray)
        y-error data to be used on errorbar and errorribbon plots
    axvspan : list (of floats)
        TODO(john) this needs an overhaul
        contains x-boundaries within which to highlight with axvspanColor
    axvspanColor : list (of floats)
        list of colors associated with each pair of x-boundaries in axvspan
    color : list (of char)
        marker and line color for associated data
    plotType : str
        'standard' - default. standard plot with lines.
        'errorBar' - same as standard but with error bars on y-data
        'errorRibbon' - same as standard but with colored error ribbons
                        on y-data
        'scatter' - scatter plot.  similar to standard above with markers only
                    but it also allows for the zData to dictate the color of
                    the points.  The attribute, cmap, is the color map.
        'contour' - contour plot.  zData governs the color
                    TODO needs an update to include colormap
    alpha : list (of floats between 0.0 and 1.0)
        Transparency of plots.  1 = opaque.  0 = invisible.
    fileName : str
        If not an empty string (which is default), plot will automatically 
        save the plot to the active directory using this filename
    aspect : str
        '' - Default.  Does nothing
        'equal' - Sets the x and y spacing on the graph equal.
    cmap : str
        TODO needs an overhaul.  should apply to both scatter and contour plots
        TODO also include link to list of other pre-defined colormaps
    shotno : list (of floats or ints)
        TODO needs an overhaul
        if not empty, these numbers will be printed in small font to the 
        lower-right of the plot.  if a subplot, in the top-most subplot
        
    Example #1
    ----------
    # import library correctly.  I'm using hbt.plot
    fig=hbt.plot.plot()
    x = np.arange(0,10,.1)
    fig.xData=[x, x]
    fig.yData=[np.cos(x),np.sin(x)]
    fig.yLegendLabel=['cos(x)','sin(x)']
    fig.xLabel = 'x'
    fig.yLabel = 'y'
    fig.plot()
    
    Example #2
    ----------
    import hbtepLib as hbt
    import numpy as np
    x=np.linspace(0,4*np.pi,1000)
    y=np.linspace(0,4*np.pi,1000)
    z=np.zeros((len(x),len(y)))
    for i in range(0,len(x)):
        for j in range(0,len(y)):
            z[i,j]=np.cos(x[i])*np.sin(y[j])
    p1=hbt.plot.plot()
    p1.xData=[x]
    p1.yData=[y]
    p1.zData=[z]
    p1.plotType='contour'
    p1.plot()
        
    """
    
    def __init__(self):
        self.title = ''
        self.subtitle = ''
        self.xLabel = ''
        self.yLabel = ''
        self.zLabel = ''
        self.xData = [] # possibility of multiple data arrays.  stored as lists.
        self.yData = [] # possibility of multiple data arrays.  stored as lists.
        self.zData = [] # for contour and scatter plot.  
        self.yLegendLabel = [] # possibility of multiple data arrays.  stored as lists.
        self.linestyle=[] # '-' is default
        self.marker=[] # Line2D.markers for list of markers.  '.' should be default ??
        self.xLim = []
        self.yLim = []
        self.legendLoc='upper right'; #'bottom left' 'center right' etc...
        self.showGrid = True
        self.yErData=[]
        self.axvspan=[]  # http://stackoverflow.com/questions/8270981/in-a-matplotlib-plot-can-i-highlight-specific-x-value-ranges
        self.axvspanColor=[]
        self.color=[]
        self.plotType='' # 'standard', 'errorbar', 'scatter', 'contour'
        self.alpha=[]#[1.0]
        self.fileName=''
        self.aspect=None  # "equal"
        self.cmap='nipy_spectral'
        self.shotno=[]
            
        
    def plot(self):
        subPlot([self])
        
    # TODO(John) add the other subfunctions from the previous prePlot class
        
        
        
        


class subPlot:
    """
    A subplot function.  Constructs a subplot from several _plotTools.plot 
    functions
    
    Parameters
    ----------
    subPlots : list (of _plotTools.plot)
        the list of plot functions that are composed into a subplot function
    fileName : str
        if not empty, the subplot is automatically saved
    plot : bool
        instructs the subplot to automatically plot itself
        
    Attributes
    ----------
    shareX : bool
        True - x-axes are shared.  
    shareY : bool
        True - y-axes are shared.  
    fileName : str
        file name that the subfig is saved as
    subPlots : list (of _plotTools.plot)
        the list of plot functions that are composed into a subplot function
        
    Subfunctions
    ------------
    plot : 
        plots itself
        
    """
    
    def __init__(self,subPlots,fileName='',plot=True):
        
        self.shareX=True;
        self.shareY=False;
        
        # show x label on only bottom-most subplot
        self._showOnlyBottomXLabel=True
        
        self.fileName=fileName;
        self.subPlots=subPlots;
        
#        # show title on only top-most title
#        self.showOnlyTopTitle==True
        
        if plot==True:
            self.plot(plotMe=plot)
        
    def plot(self,plotMe=True):
        
        # make sure data is a list
        if type(self.subPlots) is not list:
            self.subPlots=[self.subPlots];
            
        # determine the number of plots (rows and columns)
        m=len(self.subPlots);
        if type(self.subPlots[0]) is list:
            n=len(self.subPlots[0]);
        else:
            n=1;
                
        # initialize subplot
        fig, axarr = _plt.subplots(nrows=m,
                                   ncols=n, 
                                   sharex=self.shareX,#True,
                                   sharey=self.shareY,
                                   facecolor='w', 
                                   edgecolor='k', 
                                   # TODO either 
                                   # 1) find a way to maximize the figure window such that the tools in the bottom left remain 
                                   # or
                                   # 2) settle into a standard window size.  prev code:   figsize=(15*n, 2.5*m)
                                   figsize=(16, 8), # units in inches
                                    dpi=80);
                                    
        # vertical space between sub figures
        fig.subplots_adjust(hspace=0);
        
        # horizontal space between sub figures
        fig.subplots_adjust(wspace=0.1);
        
        ## adjust margins
        marginWidth=0.075;
        fig.subplots_adjust(top=1-marginWidth,bottom=marginWidth,
                            left=marginWidth, right=1-marginWidth)
                                
                              
        # iterate through sub figures (rows and columns)
        for j in range(0,m): # iterate through rows
            for i in range(0,n): # iterate through columns
                
                # plot instance
                if n==1:
                    data = self.subPlots[j];
                else:
                    data = self.subPlots[j][i]; # TODO issue here. some plots only work with [j][i] instead of [i][j]
                
                # axis handle.  plt.subplots has trouble indexing when the 
                # subplot changes from 0D to 1D to 2D.  these next lines take
                # care of this
                if n==1 and m==1:
                    ax=axarr;
                elif n==1 and m!=1:
                    ax=axarr[j];
                else:
                    ax=axarr[j][i]; # TODO issue here. some plots only work with [j][i] instead of [i][j]
                    
                # check xData formatting
                if type(data.xData) is _np.ndarray:
                    data.xData=[data.xData]
                    data.yData=[data.yData]
                    
                # iterate through data within the same plot
                for k in range(0,len(data.xData)):
                    
                    # set marker
                    if type(data.marker) is str:
                        data.marker=[data.marker]
                    try:
                        marker=data.marker[k]
                    except IndexError:
                        marker=''
                        
                    # set linestyle
                    if type(data.linestyle) is str:
                        data.linestyle=[data.linestyle]
                    try:
                        linestyle=data.linestyle[k]
                    except IndexError:
                        linestyle='-'
                        
                    # set label
                    if type(data.yLegendLabel) is str:
                        data.yLegendLabel=[data.yLegendLabel]
                    try:
                        label=data.yLegendLabel[k]
                    except IndexError:
                        label=''
                      
                    # set alpha
                    try:
                        alpha=data.alpha[k]
                    except IndexError:
                        alpha=1.0;
                    except TypeError:
                        alpha=data.alpha
                        
                    # set color
                    if type(data.color) is str:
                        data.color=[data.color]
                    try:
                        color=data.color[k];
                    except IndexError:
                        color=_cSequence[k]
                        
                    # standard plot
                    if (data.plotType == '') or (data.plotType == 'standard'):
                        ax.plot(data.xData[k], data.yData[k], marker=marker, 
                                linestyle=linestyle,label=label,alpha=alpha,
                                color=color) # , 

                    # error bar plot
                    elif (data.plotType == 'errorbar') or (data.plotType == 'errorBar'):
                        ax.errorbar(data.xData[k], data.yData[k], 
                                    yerr=data.yErData[k], marker=marker, 
                                    linestyle=linestyle,label=label) # 
                    
                    # shaded error bar plot
                    elif (data.plotType == 'errorribbon') or (data.plotType == 'errorRibbon'):
                        ax.fill_between(data.xData[k], 
                                        data.yData[k]-data.yErData[k],
                                        data.yData[k]+data.yErData[k],
                                        alpha=0.3,facecolor=data.color[k]) # 
                        ax.plot(data.xData[k], data.yData[k], marker=marker, 
                                linestyle=linestyle,label=label,color=color) 
                        
                    # scatter plot
                    elif (data.plotType == 'scatter'):
                        cmap='BuPu';
                        markerSize=35;
                        lineWidth=0.1; # lineWidth=0.5;
                        
                        cm = _plt.cm.get_cmap(cmap)
                        sc = _plt.scatter(data.xData[k], data.yData[k], 
                                          c=color, s=markerSize, cmap=cm,
                                          lw=lineWidth,alpha=data.alpha) 
                                          #vmin=0, vmax=20, 
    #                    f.colorbar(sc, axarr[i])
                        
                        #cax = f.add_axes([0.6, 0.05, 0.3, 0.02])  # [left, bottom, width, height]
                        #f.colorbar(sc,cax,orientation='horizontal')
                        a=ax.get_position()
                        cax = fig.add_axes([a.x0+a.width+0.005, a.y0, 0.01, 
                                            a.height*1.13])  
                                            # [left, bottom, width, height]
                        b=fig.colorbar(sc,cax)
                        if data.zLabel!=None:
                            b.set_label(data.zLabel)
                       
                    # contour plot
                    elif (data.plotType == 'contour'):       
                        # TODO needs an update to allow for color map
                        # TODO update so that zData can be a 1D array
                        # TODO need to check if zData is a list or an array
                        X, Y = _np.meshgrid(data.xData[k], data.yData[k])
                        X=_np.transpose(X) # because the plot function is backwards from how i think it should be
                        Y=_np.transpose(Y) # because the plot function is backwards from how i think it should be
                        Z=data.zData[k]
                        # Z[Z==0]=_np.nan  ## any value with z=0 is turned white.  
#                        cp = _plt.contourf(X, Y, Z,100)
                        _plt.contourf(X, Y, Z,100)
    #                    _plt.colorbar(cp)
    #                    f.colorbar(sc, axarr[i])
                        _plt.show()
                        
                ## annotate shot number in bottom right corner of top most plot
                if  data.shotno != [] and i==0:
                    name="";
                    for k in range (0,len(data.shotno)):
                        if name == "":
                            name = str(data.shotno[k]);
                        else:
                            name += ", " + str(data.shotno[k]);
                    axarr[0].annotate(name, xy=(0.999,0.01), 
                        xycoords='axes fraction', fontsize=6,
                        horizontalalignment='right', 
                        verticalalignment='bottom')
                        
                ## create shaded regions
                # TODO(john) update this to be a list of 2 element np.arrays
                if len(data.axvspan)!= 0:
                    for j in range(0,len(data.axvspanColor[0])):
                        ax.axvspan(data.axvspan[j*2],data.axvspan[j*2+1], color=data.axvspanColor[0][j], alpha=0.25)
            
                ## create grid    
                if data.showGrid==True:
                    ax.grid()
                    
                ## create yaxis label
                if data.yLabel!=None:
                    ax.set_ylabel(data.yLabel,fontsize=16)
                    
                # implement equal aspect ratio on plot if requested.  
                if data.aspect == "equal" or data.aspect == "Equal": 
                    # TODO(John) Not working at 100%.  Seems to have issues with xlim.
                    try:
                        ax.set_aspect('equal',adjustable='box') 
                    except ValueError:
                        print "error: share x axes must be turned off to adjust aspect ratio"
    
                ## create y limit    
                if data.yLim!=[]:
                    ax.set_ylim(data.yLim)
                    
                ## create x limit    
                if data.xLim!=[]:
                    ax.set_xlim(data.xLim)
                        
                ## create legend
                if data.yLegendLabel!=[]:
                    ax.legend(loc=data.legendLoc) #, fontsize=10 # not compatible with HBTEP server???
                    
                ## hide the top y tick label on each subplot
                if n!=1 and m!=1:
                    ax.get_yticklabels()[-1].set_visible(False)
                
                ## create x label 
                if data.xLabel != []:
                    # on bottom most plot only
                    if self._showOnlyBottomXLabel==True:
                        if j==n-1:
                            ax.set_xlabel(data.xLabel,fontsize=16);
                    # on every subplot        
                    else:
                        ax.set_xlabel(data.xLabel,fontsize=16);
                
                ## create subtitle
                if not (data.subtitle==None or data.subtitle==''):

                    ax.annotate(data.subtitle, xy=(0.5, 0.98), xycoords='axes fraction', fontsize=14,
                        horizontalalignment='center', verticalalignment='top', bbox=dict(pad=5.0, facecolor="w"))
                    
            
                ## create title on top most plot only
                if i==0 and j==0:
                    if not (data.title==None or data.title==''):
                        fig.suptitle(data.title, #verticalalignment='bottom',
                                     fontsize=18)

        ## save figure
        if self.fileName != '':
            fig.savefig(self.fileName+'.png')   
                
        # plot
        if plotMe==True:
            _plt.show()
        
#
#
##!/usr/bin/env python
#
#import pygtk
#pygtk.require('2.0')
#import gtk
#
## Global variables
#b_entry_checkbox = True
#
#class Checkbox:
#
#    def entry_checkbox(self, widget, checkbox):
#        global b_entry_checkbox
#        b_entry_checkbox = checkbox.get_active()
#        if b_entry_checkbox:
#            print "Box checked"
#        else:
#            print "Not checked"
#        return
#
#    # Main program to draw GUI
#    def __init__(self):
#        global b_entry_checkbox
#
#        # create a new window
#        app_window = gtk.Window(gtk.WINDOW_TOPLEVEL)
#        app_window.set_size_request(500, 100)
#        app_window.set_border_width(10)
#        app_window.set_title("My program title")
#        app_window.connect("delete_event", lambda w,e: gtk.main_quit())
#
#        frame_checkbox = gtk.Frame("Check for true:")
##        frame_checkbox.set_shadow_type(gtk.SHADOW_IN)
#
#        app_window.add(frame_checkbox)
#
#        check_box = gtk.CheckButton("Checkbox text string")
#        check_box.connect("toggled", self.entry_checkbox, check_box)
#        check_box.set_active(True)  # Set the defaut
#        check_box.show()
#
#        frame_checkbox.add(check_box)
#        frame_checkbox.show()
#        
#        
#        frame_checkbox2 = gtk.Frame("asdf")
#        app_window.add(frame_checkbox2)
#        frame_checkbox2.show()
#        
#        app_window.show()
#        return
#
#def main():
#    gtk.main()
#    return 0
#
#if __name__ == "__main__":
#    Checkbox()
#    main()

