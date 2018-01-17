import numpy as _np
import matplotlib.pyplot as _plt
from matplotlib.colors import LinearSegmentedColormap as _lsc # for creating your own colormaps
#import sys as _sys
#import _processData as _process
#from matplotlib.colors import LinearSegmentedColormap as _lsc


# color sequence for plotting.  add more colors if you need more than 7.  
#_cSequence=['b', 'r', 'g', 'k', 'm', 'c', 'y', 'brown']
_cSequence=['red', 'black',"#1f77b4", "#ff7f0e", "#2ca02c",  "#9467bd", "#8c564b", "#d62728","#e377c2", "#7f7f7f", "#bcbd22", "#17becf"] #https://github.com/vega/vega/wiki/Scales#scale-range-literals


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
        the z-data of the plot.  applicable for contour, scatter plots and both error plots. 
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
#    yErData : list (of numpy.ndarray)
#        y-error data to be used on errorbar and errorribbon plots
    axvspan : list (of floats)
        TODO (john) this needs an overhaul
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
    p1=hbt.plot.plot(xLabel='x',yLabel='y',zLabel='z',title='title',
                     subtitle='subtitle')
    x = np.arange(0,10,.1)
    p1.addTrace(xData=x,yData=np.cos(x),yLegendLabel='cos(x)')
    p1.addTrace(xData=x,yData=np.sin(x),yLegendLabel='sin(x)')
    p1.plot()
    
    Example #2
    ----------
    # import libraries
    import hbtepLib as hbt; reload(hbt)
    import numpy as np
    import matplotlib.pyplot as plt   
    # generate data
    delta = 0.025
    x = np.arange(-3.0, 3.01, delta)
    y = np.arange(-1.0, 4.01, delta)
    X, Y = np.meshgrid(x, y)
    Z1 = plt.mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    Z2 = plt.mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    Z = 10 * (Z1 - Z2)    
    # plot contour
    p1=hbt.plot.plot(xLabel='x',yLabel='y',zLabel='z',title='title',
                     subtitle='subtitle',plotType='contour')
    p1.addTrace(x,y,Z)
    # alternatively, use: p1.addTrace(X,Y,Z)
    p1.plot()
        
    """
    
    def __init__(self,title='',titleFontSize=18,
                 subtitle = '',subtitleFontSize=14,
                 xLabel = '',yLabel = '',zLabel = '',axisLabelFontSize=16,
                 legendFontSize=10,legendLoc='upper right',
                 xLim=[],yLim=[],
                 showGrid=True,
                 axvspan=[],axvspanColor=[],
                 plotType='standard',
                 fileName='',
                 aspect=None,
                 colorMap='nipy_spectral',centerColorMapAroundZero=False,
                 shotno=[],
                 shotnoFontSize=8):
        self.title = title
        self.titleFontSize=titleFontSize
        self.subtitle = subtitle
        self.subtitleFontSize=subtitleFontSize
        self.xLabel = xLabel
        self.yLabel = yLabel
        self.zLabel = zLabel
        self.axisLabelFontSize=axisLabelFontSize
        self.xData = [] # possibility of multiple data arrays.  stored as lists.
        self.yData = [] # possibility of multiple data arrays.  stored as lists.
        self.zData = [] # for contour and scatter plot.  
        self.yLegendLabel = [] # possibility of multiple data arrays.  stored as lists.
        self.legendFontSize=legendFontSize
        self.linestyle=[] # '-' is default
        self.marker=[] # Line2D.markers for list of markers.  '.' should be default ??
        self.xLim = xLim
        self.yLim = yLim
        self.legendLoc=legendLoc; #'bottom left' 'center right' etc...
        self.showGrid = showGrid
        self.axvspan=axvspan  # http://stackoverflow.com/questions/8270981/in-a-matplotlib-plot-can-i-highlight-specific-x-value-ranges
        self.axvspanColor=axvspanColor
        self.color=[]
        self.plotType=plotType # 'standard', 'errorbar', 'scatter', 'contour'
        self.alpha=[]#[1.0]
        self.fileName=fileName
        self.aspect=aspect  # "equal"
        self.colorMap=colorMap #https://matplotlib.org/examples/color/colormaps_reference.html
        self.centerColorMapAroundZero=centerColorMapAroundZero
        self.shotno=shotno
        self.shotnoFontSize=shotnoFontSize
            
        
    def plot(self):
        subPlot([self])
        
    # TODO (John) add the other subfunctions from the previous prePlot class
        
    def addTrace(self,xData,yData,zData=[],yLegendLabel='',alpha=1.0,
                 linestyle='-',marker='',color=''):
        """
        Adds a new trace to the plot
        """
        m=len(self.xData)
        
        self.xData.append(xData);
        self.yData.append(yData);
        self.zData.append(zData);
        self.yLegendLabel.append(yLegendLabel)
        self.alpha.append(alpha)
        self.linestyle.append(linestyle)
        self.marker.append(marker)
        if color=='':
            self.color.append(_cSequence[m])
        else:
            self.color.append(color)
            
    def removeTrace(self,index):
        """
        Remove one or more traces from the plot
        
        Parameters
        ----------
        index : list or numpy.ndarray
            list or array of one of more indices within the plot to be removed
        """
        # make sure index is a numpy.ndarray
        if type(index) is not _np.ndarray and type(index) is not list:
            index=_np.array([index]);
        elif type(index) is list:
            index=_np.array(index);
        
        m=len(index)
        
        for i in range(0,m):
            del self.xData[index[i]]
            del self.yData[index[i]]
            del self.zData[index[i]]
            del self.yLegendLabel[index[i]]
            del self.alpha[index[i]]
            del self.linestyle[index[i]]
            del self.marker[index[i]]
            del self.color[index[i]]
            index-=1;
        
    def mergePlots(self,newPlot):
        """
        Merges a new plot into the existing plot
        
        Parameters
        ----------
        newPlot : _plotTools.plot
            the plot to be merged into this one
        """
        
        if type(newPlot.xData)!=list:
            newPlot.xData=[newPlot.xData]
            
        m=len(newPlot.xData); 
        n=len(self.xData)
        
        for i in range(0,m):
            if self.alpha!=[]:
                if newPlot.alpha!=[]:
                    self.alpha.append(newPlot.alpha[i])
                else:
                    self.alpha=[];
                
            if self.color!=[]:
                if newPlot.color!=[]:
                    if newPlot.color[i] in self.color:
                        self.color.append(_cSequence[n+i])
                    else:
                        self.color.append(newPlot.color[i])
                else:
                    self.color=[];
                    
            if self.marker!=[]:
                if newPlot.marker!=[]:
                    self.marker.append(newPlot.marker[i])
                else:
                    self.marker=[];
                
            if self.linestyle!=[]:
                if newPlot.linestyle!=[]:
                    self.linestyle.append(newPlot.linestyle[i])
                else:
                    self.linestyle=[];
                
            if self.yLegendLabel!=[]:
                if newPlot.yLegendLabel!=[]:
                    self.yLegendLabel.append(newPlot.yLegendLabel[i])
                else:
                    self.yLegendLabel=[];
                    
            if self.zData!=[]:
                if newPlot.zData!=[]:
                    self.zData.append(newPlot.zData[i])
                else:
                    self.zData=[];
                
            self.xData.append(newPlot.xData[i])
            self.yData.append(newPlot.yData[i])
        


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
        file name that the subfig is saved as.  if not an empty string, the 
        image is saved.  
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
        
        self.fileName=fileName;
        self.subPlots=subPlots;
        
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
                                   figsize=(16, 8), # units in inches.  not sure how this actually maps to a screen though since it doesn't actually measure 16 inches ...
                                    dpi=80);
                                    
        # vertical space between sub figures
        if self.shareX==True:
            fig.subplots_adjust(hspace=0);
        else:
            fig.subplots_adjust(hspace=0.25);
        
        # horizontal space between sub figures
        fig.subplots_adjust(wspace=0.1);
        
        ## adjust margins
        marginWidth=0.075; # 0 to 1.0
        fig.subplots_adjust(top=1-marginWidth,bottom=marginWidth,
                            left=marginWidth, right=1-marginWidth)
                              
        # iterate through sub figures (rows and columns)
        for j in range(0,m): # iterate through rows
            for i in range(0,n): # iterate through columns
                
                # plot instance
                if n==1:
                    data = self.subPlots[j];
                else:
                    data = self.subPlots[j][i];
                
                # axis handle.  plt.subplots has trouble indexing when the 
                # subplot changes from 0D to 1D to 2D.  these next lines take
                # care of this
                if n==1 and m==1:
                    ax=axarr;
                elif n==1 and m!=1:
                    ax=axarr[j];
                else:
                    ax=axarr[j][i]; 
                    
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
                                    yerr=data.zData[k], marker=marker, 
                                    linestyle=linestyle,label=label) # 
                    
                    # shaded error bar plot
                    elif (data.plotType == 'errorribbon') or (data.plotType == 'errorRibbon'):
                        ax.fill_between(data.xData[k], 
                                        data.yData[k]-data.zData[k],
                                        data.yData[k]+data.zData[k],
                                        alpha=0.3,facecolor=data.color[k]) # 
                        ax.plot(data.xData[k], data.yData[k], marker=marker, 
                                linestyle=linestyle,label=label,color=color) 
                        
                    # scatter plot
                    elif (data.plotType == 'scatter'):
                        # NOTE this will not work well if there are multiple columns of subplots
                        
                        # the color can be an array (from zData) or a single value (from color)
                        if data.zData!=[]:
                            c=data.zData[k]
                        else:
                            c=color
                            
                        # if there is more than 1 column, increase horizontal spacing between subfigures
                        if m>1:
                            if i==0 and j==0:
                                fig.subplots_adjust(wspace=0.35);
                        
                        # colormap
                        cmap=data.colorMap;
                        cm = _plt.cm.get_cmap(cmap)
                        
                        # marker parameters
                        markerSize=35;
                        lineWidth=0.1; # the circles look silly without at least a very small outline
                        
                        # make scatter plot
                        print c
                        p1 = ax.scatter(data.xData[k], data.yData[k], 
                                          c=c, s=markerSize, cmap=cm,
                                          lw=lineWidth,alpha=alpha) 
                                          
                        # place color bar                  
                        a=ax.get_position() # position of subfigure
                        cax = fig.add_axes([a.x0+a.width+0.010, # left
                                            a.y0+0.01, # bottom
                                            0.01, # width
                                            a.height-0.02])  # height
                                            # [left, bottom, width, height]
                        b=fig.colorbar(p1,cax)
                        
                        # zLabel
                        if data.zLabel!=None:
                            b.set_label(data.zLabel,fontsize=data.axisLabelFontSize)
                       
                    # contour plot
                    elif (data.plotType == 'contour'):       
                        # TODO update so that zData can be a 1D array
                    
                        # check to see if meshgrid has already been applied to xData and yData
                        if _np.shape(_np.shape(data.xData[k]))[0]==1:
                            # it's 1D data so apply meshgrid
                            X, Y = _np.meshgrid(data.xData[k], data.yData[k])
                        elif _np.shape(_np.shape(data.xData[k]))[0]==2:
                            # it's 2D data and likely already been meshed
                            X=data.xData[k]
                            Y=data.yData[k]
                        Z=data.zData[k]
                        
                        # create colormap      
                        if type(data.colorMap) is str:
                            cm = _plt.cm.get_cmap(data.colorMap)
                        else:
                            cm = data.colorMap
                        
                        # find min and max of zData.  Center colormap around 
                        # 0 if requested
                        if data.centerColorMapAroundZero==True:
                            vmax=_np.abs(Z).max()
                            vmin=-vmax
                        else:
                            vmax=Z.max()
                            vmin=Z.min()
                        
                        # create contour plot
                        p1=_plt.contourf(X, Y, Z,100,cmap=cm,# 100 is equal to the number of color deviations
                                         vmin=vmin, vmax=vmax) 
                        
                        # place color bar                  
                        a=ax.get_position() # position of subfigure
                        cax = fig.add_axes([a.x0+a.width+0.010, # left
                                            a.y0+0.01, # bottom
                                            0.01, # width
                                            a.height-0.02])  # height
                        b=fig.colorbar(p1,cax)
                        
                        # zLabel
                        if data.zLabel!=None:
                            b.set_label(data.zLabel,fontsize=data.axisLabelFontSize)
                        
                ## annotate shot number(s) in bottom right corner of top-left most plot
                if  data.shotno != [] and i==0 and j==0:
                    
                    # make sure shotno is a list
                    if type(data.shotno) is not list:
                        data.shotno=[data.shotno]
                        
                    # construct string
                    name="";
                    for k in range (0,len(data.shotno)):
                        if name == "":
                            name = str(data.shotno[k]);
                        else:
                            name += ", " + str(data.shotno[k]);
                            
                    # place string
                    ax.annotate(name, xy=(0.999,0.01), 
                        xycoords='axes fraction', 
                        fontsize=data.shotnoFontSize, # 6
                        horizontalalignment='right', 
                        verticalalignment='bottom')
                        
                ## create shaded regions
                # TODO (john) update this to be a list of 2 element np.arrays
                if len(data.axvspan)!= 0:
                    for j in range(0,len(data.axvspanColor[0])):
                        ax.axvspan(data.axvspan[j*2],data.axvspan[j*2+1], color=data.axvspanColor[0][j], alpha=0.25)
            
                ## create grid    
                if data.showGrid==True:
                    ax.grid()
                    
                ## create yaxis label
                if data.yLabel!=None:
                    if self.shareY==False:
                        ax.set_ylabel(data.yLabel,fontsize=data.axisLabelFontSize)
                    elif self.shareY==True and i==0:
                        # if y-axes are shared on all subplots, show only y-labels on left-most plots
                        ax.set_ylabel(data.yLabel,fontsize=data.axisLabelFontSize)
                    
                # implement equal aspect ratio on plot if requested.  
                if data.aspect == "equal" or data.aspect == "Equal": 
                    # TODO (John) Not working at 100%.  Seems to have issues with xlim.  Fix.  Not sure how...
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
                    ax.legend(loc=data.legendLoc, prop={'size': data.legendFontSize})
                    
                ## hide the top y tick label on each subplot
                if j!=0: 
                    ax.get_yticklabels()[-1].set_visible(False)
                
                ## create x label 
                if data.xLabel != []:
                    if self.shareX==True: # on bottom most plot only
                        if j==m-1:
                            ax.set_xlabel(data.xLabel,fontsize=data.axisLabelFontSize);
                    else: # on every subplot 
                        ax.set_xlabel(data.xLabel,fontsize=data.axisLabelFontSize);
                          
                ## create subtitle
                if not (data.subtitle==None or data.subtitle==''):
                    ax.annotate(data.subtitle, xy=(0.5, 0.98), 
                                xycoords='axes fraction', 
                                fontsize=data.subtitleFontSize,
                                horizontalalignment='center', 
                                verticalalignment='top', 
                                bbox=dict(pad=5.0, facecolor="w"))
                    
                ## create title on top most plot only
                if i==0 and j==0:
                    if not (data.title==None or data.title==''):
                        fig.suptitle(data.title, #verticalalignment='bottom',
                                     fontsize=data.titleFontSize)

        ## save figure
        if self.fileName != '':
            fig.savefig(self.fileName+'.png')   
                
        # plot
        if plotMe==True:
            _plt.show()
 
          
            
def _red_green_colormap():
    '''A colormap with a quick red-green transition at 0.5
    Default HBTEP colormap
    
    Example use (use explicit norm to center around 0):
    
    sig_range = _np.abs(signals).max()
    _plt.contourf(times*1e3, thetas, signals, 50, cmap=red_green_colormap(), 
                norm=Normalize(vmin=-sig_range, vmax=sig_range))
    '''
    
    cdict = {'red':   [(0.0,  0.0, 0.0),
                       (0.25, 0.0, 0.0),
                       (0.5,  0.0, 0.2),
                       (0.75, 1.0, 1.0),
                       (1.0,  1.0, 0.0)],
             
             'green': [(0.0,  0.0, 0.3),
                       (0.25, 1.0, 1.0),
                       (0.5,  0.2, 0.0),
                       (0.75, 0.0, 0.0),
                       (1.0,  1.0, 0.0)],
             
             'blue':  [(0.0,  0.0, 1.0),
                       (0.25, 0.0, 0.0),
                       (0.5,  0.0, 0.0),
                       (0.75, 0.0, 0.0),
                       (1.0,  0.12, 0.0)]} 

    return _lsc('Red-Green', cdict)

        
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

