import numpy as _np
import matplotlib.pyplot as _plt
#import _processData as _process
#from matplotlib.colors import LinearSegmentedColormap as _lsc


# color sequence for plotting.  add more colors if you need more than 8.  
_cSequence=['b', 'r', 'g', 'k', 'm', 'c', 'y', 'w']


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
        the z-data of the plot.  must be stored as a list of arrays.  multiple
        entries isn't supported because having multiple 3D plots on a single
        plot doesn't make sense.  
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
        
    # TODO(John) finish entries...
        
        
    Example
    -------
    # import library correctly.  I'm using hbt.plot
    fig=hbt.plot.plot()
    x = np.arange(0,10,.1)
    fig.xData=[x, x]
    fig.yData=[np.cos(x),np.sin(x)]
    fig.yLegendLabel=['cos(x)','sin(x)']
    fig.xLabel = 'x'
    fig.yLabel = 'y'
    fig.plot()
        
    """
    
    def __init__(self):
        self.title = ''
        self.subtitle = ''
        self.xLabel = ''
        self.yLabel = ''
        self.zLabel = ''
        self.xData = [] # possibility of multiple data arrays.  stored as lists.
        self.yData = [] # possibility of multiple data arrays.  stored as lists.
        self.zData = [] # for contour plot.  must be 2D
        self.yLegendLabel = [] # possibility of multiple data arrays.  stored as lists.
        self.linestyle=[] # '-' is default
        self.linewidth=[]
        self.marker=[] # Line2D.markers for list of markers.  '.' should be default ??
        self.xLim = []
        self.yLim = []
        self.legendLoc='upper right'; #'bottom left' 'center right' etc...
        self.showGrid = True
        self.yErData=[[]]
        self.axvspan=[]  # http://stackoverflow.com/questions/8270981/in-a-matplotlib-plot-can-i-highlight-specific-x-value-ranges
        self.axvspanColor=[[]]
        self.color=[]
        self.plotType='' # 'standard', 'errorbar', 'scatter', 'contour'
        self.colorData=[];  # scatter plot color
        self.alpha=[]#[1.0]
        self.fileName=''
        self.aspect=None  # "equal"
        self.cmap='nipy_spectral'
        self.shotno=[]
            
        
    def plot(self):
        subPlot([self])
        
    # TODO(John) add the other subfunctions from the previous prePlot class
        
        
        
        


class subPlot:
    
    def __init__(self,subPlots,fileName='',plot=True):
        
        # show x label on only bottom-most subplot
        self._showOnlyBottomXLabel=True
        
        self.fileName=fileName;
        self.subPlots=subPlots;
        
#        # show title on only top-most title
#        self.showOnlyTopTitle==True
        
        if plot==True:
            self.plot()
        
    def plot(self):
        
        # make sure x and y are lists
        if type(self.subPlots) is not list:
            self.subPlots=[self.subPlots];
            
        # determine the number of plots (rows and columns)
        m=len(self.subPlots);
        if type(self.subPlots[0]) is list:
            n=len(self.subPlots[0]);
        else:
            n=1;
            self.subPlots=[self.subPlots];
                
        # initialize subplot
        fig, axarr = _plt.subplots(nrows=m,
                                   ncols=n, 
                                   sharex=True,#True,
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
        for i in range(0,m):
            for j in range(0,n):
                
                # plot instance
                data = self.subPlots[j][i];
                
                # axis handle.  plt.subplots has trouble indexing when the 
                # subplot changes from 0D to 1D to 2D.  these next lines take
                # care of this
                if n==1 and m==1:
                    ax=axarr;
                elif n==1 and m!=1:
                    ax=axarr[i];
                else:
                    ax=axarr[j][i];
                    
                # check xData formatting
                if type(data.xData) is _np.ndarray:
                    data.xData=[data.xData]
                    data.yData=[data.yData]
                    
                # iterate through data within the same plot
                for k in range(0,len(data.xData)):
                    
#                    # set xData
#                    try:
#                        xData=data.xData[k]
#                    except IndexError:
#                        xData=data.xData;
#                        
#                    # set yData
#                    try:
#                        yData=data.yData[k]
#                    except IndexError:
#                        yData=data.yData;
#                    
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
#                        try:
                        alpha=data.alpha
#                        except TypeError:
#                            alpha
                        
                    # set color
                    if type(data.color) is str:
                        data.color=[data.color]
                    try:
                        color=data.color[k];
                    except IndexError:
#                        print k
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
    #                    lineWidth=0.5;
                        lineWidth=0.1;
                        
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
                        X, Y = _np.meshgrid(data.xData[k], data.yData[k])
                        X=_np.transpose(X)
                        Y=_np.transpose(Y)
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
        _plt.show()