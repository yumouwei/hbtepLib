"""
I AM TRYING TO PHASE THIS FUNCTION OUT!!!  STOP USING !!!
"""

#!/usr/bin/env python
import numpy as _np
import matplotlib.pyplot as _plt
import _processData as _pd
from matplotlib.colors import LinearSegmentedColormap as _lsc


# color sequence for plotting.  add more colors if you need more than 8.  
_cSequence=['b', 'r', 'g','m', 'k', 'c', 'y', 'w']



def scatterPlot(x,y, color, cmap='BuPu',markerSize=35, lineWidth=0.5):
    """
    Simple scatter plot with colors
    
    example: 
    x=_np.linspace(0,4*3.14159,30)
    y=_np.sin(x)
    pdt.scatterPlot(x,y,color=y)
    """
    # TODO(John) Implement this code with preplot below
    cm = _plt.cm.get_cmap(cmap)
    sc = _plt.scatter(x, y, c=color, s=markerSize, cmap=cm,lw=lineWidth) #vmin=0, vmax=20, 
    _plt.colorbar(sc)
    _plt.show()

    
class prePlot:
    """
    Data class that contains everything needed to plot itself
    """
    def __init__(self):
        self.title = None
        self.subtitle = None
        self.xLabel = None
        self.yLabel = None 
        self.zLabel = None
        self.xData = [] # possibility of multiple data arrays.  stored as lists.
        self.yData = [] # possibility of multiple data arrays.  stored as lists.
        self.zData = [] # for contour plot.  must be 2D
        self.yLegendLabel = [] # possibility of multiple data arrays.  stored as lists.
        self.linestyle=[] # '-' is default
        self.linewidth=[]
        self.marker=[] # Line2D.markers for list of markers.  '.' should be default ??
        self.xLim = None
        self.yLim = None
        self.legendLoc='upper right';
        self.showGrid = True
        self.yErData=[[]]
        self.axvspan=[]  # http://stackoverflow.com/questions/8270981/in-a-matplotlib-plot-can-i-highlight-specific-x-value-ranges
        self.axvspanColor=[[]]
        self.color=[]
        self.plotType='' # 'standard', 'errorbar', 'scatter', 'contour'
        self.colorData=[];  # scatter plot color
        self.alpha=[]#[1.0]
        self.fileName=''
        self.titleFontSize=12
        self.aspect=None  # "equal"
        self.cmap='nipy_spectral'
        self.shotno=[]
    
    def plot(self):
        """
        This subfunction causes the data class to plot itself
        """
        _plt.figure(num=None, figsize=(10, 5), dpi=80, facecolor='w', edgecolor='k')
    
        if isinstance(self.xData,list):
            m=len(self.xData);
        else:
            m=1;
            self.xData=[self.xData];
            self.yData=[self.yData];
#            self.yLegendLabel=[self.yLegendLabel];    
#            if self.yLegendLabel == []:
#                for j in range(0,m):
#                    self.yLegendLabel.append('')

        
        for i in range(0,m):
            
            # make sure that yLegendLabel is formatted correctly
            if self.yLegendLabel==[]:
                yLegendLabel='';
            elif isinstance(self.yLegendLabel,list):
                yLegendLabel = self.yLegendLabel[i];
            else:
                yLegendLabel = self.yLegendLabel;
                
            # make sure that alpha is formatted correctly
            if self.alpha==[]:
                alpha=1.0;
            elif isinstance(self.alpha,list):
                alpha = self.alpha[i];
            else:
                alpha=self.alpha;
            
            # make sure that linestyle is formatted correctly
            if self.linestyle == []:
                linestyle = '-';
            elif isinstance(self.linestyle,list):
                linestyle = self.linestyle[i];
            else:
                linestyle = self.linestyle;

            # make sure that marker is formatted correctly
            if self.marker == []:
                marker='None';
            else:
                marker=self.marker[i];
                
            # make sure that color is formatted correctly
            if self.color == []:
                color=_cSequence[i]
            else:
                color=self.color[i];

            # standard plot
            if (self.plotType == '') or (self.plotType == 'standard'):
                _plt.plot(self.xData[i],self.yData[i],label=yLegendLabel, marker=marker, linestyle=linestyle, color=color,alpha=alpha)
    
            # standard plot with error bars
            elif (self.plotType == 'errorbar') or (self.plotType == 'errorBar'):
                print "not currently working..."
                # TODO(John) Implement
                
            # scatter plot
            elif (self.plotType == 'scatter'):
                #cmap='nipy_spectral'#'gist_rainbow'#'nipy_spectral'#'BuPu' # 'brg'
                markerSize=25; #35
                lineWidth=0.0; #0.05
                
                
                if isinstance(self.cmap,basestring):
                    cm = _plt.cm.get_cmap(self.cmap)
                else:
                    cm=self.cmap
                sc = _plt.scatter(self.xData[i], self.yData[i], c=self.colorData[i], s=markerSize, cmap=cm,lw=lineWidth,alpha=alpha) #vmin=0, vmax=20, 
                cb=_plt.colorbar(sc)
                _plt.show()
#                cb=_plt.gcf().colorbar(sc)
                cb.set_label(self.zLabel)
            
            # contour plot
            elif (self.plotType == 'contour'):                
                X, Y = _np.meshgrid(self.xData[i], self.yData[i])
                X=_np.transpose(X)
                Y=_np.transpose(Y)
                Z=self.zData[i]
                # Z[Z==0]=_np.nan  ## any value with z=0 is turned white.  
                cp = _plt.contourf(X, Y, Z,100)
                _plt.colorbar(cp)
                _plt.show()
                
        # implement equal aspect ratio on plot if requested.  
        if self.aspect == "equal" or self.aspect == "Equal": 
            # TODO(John) Not working at 100%.  Seems to have issues with xlim.
            _plt.axes().set_aspect('equal',adjustable='box')  
            print "setting aspect to equal"
          
        # redundant with below (I think)
#        # create legend           
#        if m > 1 and self.yLegendLabel != None:
#            _plt.legend(loc=self.legendLoc)
        
        # add x and y labels
        if self.xLabel != None:
            _plt.xlabel(self.xLabel,fontsize=16)
        if self.yLabel != None:
            _plt.ylabel(self.yLabel,fontsize=16)
#        if self.zLabel != None:
#            print 'I havent figured out how to do this...  '   
            
        # add x and y limits
        if self.xLim != None:
            print self.xLim
            _plt.xlim(self.xLim)         
        if self.yLim != None:
            _plt.ylim(self.yLim)   
            
        # add origin-axes if applicable
        # TODO(John) This code is returning a "'AxesSubplot' object has no 
            # attribute 'get_yLim'" error which I don't understand.  Fix.
#        axes=_plt.gca()
##        return axes
##        import time
##        time.sleep(1)
#        print axes
#        print axes.get_yLim()
#        print axes.get_xlim()
#        (yLimLower,yLimUpper)=axes.get_yLim()
#        (xLimLower,xLimUpper)=axes.get_xlim()
#        if yLimLower < 0 and yLimUpper > 0:
#            _plt.plot([0,0],[yLimLower,yLimUpper])
#        if xLimLower < 0 and xLimUpper > 0:
#            _plt.plot([xLimLower,xLimUpper],[0,0])

            
        # add title
        if self.title != None:
            _plt.title(self.title,fontsize=self.titleFontSize)

        # create subtitle
        if self.subtitle != None:
            _plt.annotate(self.subtitle, xy=(0.5, 0.98), xycoords='axes fraction', fontsize=14,
                horizontalalignment='center', verticalalignment='top', bbox=dict(pad=5.0, facecolor="w"))
           
        # create legend
        if self.yLegendLabel!=[]:
            _plt.legend(loc=self.legendLoc, fontsize=10)

        # save figure
        if self.fileName != '':
            _plt.savefig(self.fileName+'.png')
                
        # add grid
        if self.showGrid == True:
            _plt.grid()
            
        # add shot number(s) to bottom right corner
        if  self.shotno != []:
            name="";
            for k in range (0,len(self.shotno)):
                if name == "":
                    name = str(self.shotno[k]);
                else:
                    name += ", " + str(self.shotno[k]);
    #        print("name: " + name)
            _plt.annotate(name, xy=(0.999,0.01), xycoords='axes fraction', fontsize=6,
                    horizontalalignment='right', verticalalignment='bottom')
            
        # draw figure.  not sure if this is needed.
        _plt.draw()
            

    def combinePlot(self,newPlot):
         self.combinePlot2(newPlot);
#        """
#        Allows for a plot to easily be combined with another plot
#        
#        TODO(JOHN) this code still needs a lot of work.  it should automatically define alpha, etc if not previously defined
#        """
##        
##        n=len(self.xData)
##        for i in range(0,n):
##            if self.marker[i]==[]
##        
#        m=len(newPlot.xData);
#        print "self.alpha="
#        print self.alpha
#        print "newPlot.alpha="
#        print newPlot.alpha
#        for i in range(0,m):
#            self.xData.append(newPlot.xData[i])
#            self.yData.append(newPlot.yData[i])
#            
#            if newPlot.marker[i] == []:
#                self.marker.append('')
#            else:
#                self.marker.append(newPlot.marker[i])
##            if newPlot.alpha[i] == []:
##                self.marker.append('')
##            else:
##                self.marker.append(newPlot.marker[i])
#                
#            self.linestyle.append(newPlot.linestyle[i])
##            if newPlot.color == []:
##                self.color.append(newPlot.color[i])
##            print self.alpha
##            print newPlot.alpha
#            self.alpha.append(newPlot.alpha)
#            self.yLegendLabel.append(newPlot.yLegendLabel[i])
            
    def keepEntries(self,indices=[1,2]):
        """
        self may contain multiple plots, but you don't wish to plot all of them.  
        this function returns a new prePlot class that contains only the plots
        you specify in indices
        """
        import copy
        a2=copy.copy(self)            
        a2.xData=[]
        a2.yData=[]
        a2.alpha=[]
        a2.marker=[]
        a2.yLegendLabel=[]
        a2.linestyle=[]
        a2.color=[]
        
        for i in range(0,len(indices)):
            a2.xData.append(self.xData[indices[i]])
            a2.yData.append(self.yData[indices[i]])
            if self.alpha!=[]:
                a2.alpha.append(self.alpha[indices[i]])
            if self.marker!=[]:
                a2.marker.append(self.marker[indices[i]])
            if self.linestyle!=[]:
                a2.linestyle.append(self.linestyle[indices[i]])
            if self.color!=[]:
                a2.color.append(self.color[indices[i]])
            if self.yLegendLabel!=[]:
                a2.yLegendLabel.append(self.yLegendLabel[indices[i]])
        return a2      

            
    def combinePlot2(self,newPlot):
        """
        Allows for a plot to easily be combined with another plot (prePlot)
        
        this is my attempt to improve the previous combinePlot subfunction
        """

        m=len(newPlot.xData);        
        
        for i in range(0,m):
            print i
            if self.alpha!=[]:
                if newPlot.alpha!=[]:
                    self.alpha.append(newPlot.alpha[i])
                else:
                    self.alpha=[];
                
            if self.color!=[]:
                if newPlot.color!=[]:
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
                
            self.xData.append(newPlot.xData[i])
            self.yData.append(newPlot.yData[i])
            

 
 
class preSubPlot:
    """
    Data class for sub plots.  This structure effective contains a number of
    prePlot classes and plots each as a subplot.  
    """
    # TODO(John) This class needs an overhaul
    def __init__(self,data=[],numRows=2,shareXAxis="all",hSpace=0, fileName='',shotNo=[],plot=False,yScaleSize=2.5):
        
        if shareXAxis=='none':
            self.hSpace=0.15; #0.05;
            self.showOnlyBottomXLabel=False
        else:
            self.hSpace=0.; #0.05;
            self.showOnlyBottomXLabel=True;
        self.data=data; # multiple prePlot structures are stored here
        if len(data)!=0:
            self.numRows=len(data);
        else:
            self.numRows=numRows;
        self.showOnlyTopTitle=True;
        self.shareXAxis=shareXAxis;  #sharex [] must be one of [u'all', u'row', u'col', u'none']
        self.numCols=1;
#        self.shotNo=[shotNo];
        self.fileName=fileName
        self.yScaleSize=yScaleSize
        
        if plot==True:
            self.plot()

        
    def plot(self):
        m=self.numRows;
        n=self.numCols;
        f, axarr = _plt.subplots(m,n, sharex=self.shareXAxis,figsize=(15*n, self.yScaleSize*m), dpi=80, facecolor='w', edgecolor='k')

        
        ## iterate through subfigs
        for i in range(0,m):
            data=self.data[i];

            



            if data.yLegendLabel == []:
                for j in range(0,m):
                    data.yLegendLabel.append('')
            ## for each subfig, iterate through each line of data
            for k in range(0,len(data.xData)):
                if data.linestyle == []:
                    linestyle = '-';
                else:
                    linestyle = data.linestyle[k];
                    
                if data.marker == []:
                    marker='None';
                else:
                    marker=data.marker[k];

                if data.alpha == []:
                    alpha=1.0;
                else:
                    print k
                    alpha=data.alpha[k];

                # make sure that color is formatted correctly
                if data.color == []:
                    color=_cSequence[k]
                else:
                    color=data.color[k];
                #            if data.color!=[]:
                #                axarr[i].set_color_cycle(data.color);
                
                if (data.plotType == '') or (data.plotType == 'standard'):
                    axarr[i].plot(data.xData[k], data.yData[k], marker=marker, linestyle=linestyle,label=data.yLegendLabel[k],alpha=alpha,color=color) # , 

                
                elif (data.plotType == 'errorbar') or (data.plotType == 'errorBar'):
                    axarr[i].errorbar(data.xData[k], data.yData[k], yerr=data.yErData[k], marker=marker, linestyle=linestyle,label=data.yLegendLabel[k]) # 
                    
                elif (data.plotType == 'errorribbon') or (data.plotType == 'errorRibbon'):
                    axarr[i].fill_between(data.xData[k], data.yData[k]-data.yErData[k],data.yData[k]+data.yErData[k],alpha=0.3,facecolor=data.color[k]) # 
                    axarr[i].plot(data.xData[k], data.yData[k], marker=marker, linestyle=linestyle,label=data.yLegendLabel[k],color=data.color[k]) # , 
                    
                elif (data.plotType == 'scatter'):
                    cmap='BuPu';
                    markerSize=35;
#                    lineWidth=0.5;
                    lineWidth=0.1;
                    
                    cm = _plt.cm.get_cmap(cmap)
                    sc = _plt.scatter(data.xData[k], data.yData[k], c=data.colorData[k], s=markerSize, cmap=cm,lw=lineWidth,alpha=data.alpha) #vmin=0, vmax=20, 
#                    f.colorbar(sc, axarr[i])
                    
                    #cax = f.add_axes([0.6, 0.05, 0.3, 0.02])  # [left, bottom, width, height]
                    #f.colorbar(sc,cax,orientation='horizontal')
                    a=axarr[i].get_position()
                    cax = f.add_axes([a.x0+a.width+0.005, a.y0, 0.01, a.height*1.13])  # [left, bottom, width, height]
                    b=f.colorbar(sc,cax)
                    if data.zLabel!=None:
                        b.set_label(data.zLabel)
                    
                elif (data.plotType == 'contour'):                
                    X, Y = _np.meshgrid(data.xData[k], data.yData[k])
                    X=_np.transpose(X)
                    Y=_np.transpose(Y)
                    Z=data.zData[k]
                    # Z[Z==0]=_np.nan  ## any value with z=0 is turned white.  
                    cp = _plt.contourf(X, Y, Z,100)
#                    _plt.colorbar(cp)
#                    f.colorbar(sc, axarr[i])
                    _plt.show()
            
            ## annotate shot number
            if  data.shotno != [] and i==0:
                name="";
                for k in range (0,len(data.shotno)):
                    if name == "":
                        name = str(data.shotno[k]);
                    else:
                        name += ", " + str(data.shotno[k]);
        #        print("name: " + name)
                axarr[0].annotate(name, xy=(0.999,0.01), xycoords='axes fraction', fontsize=6,
                        horizontalalignment='right', verticalalignment='bottom')
                                        
            ## create shaded regions
            if len(data.axvspan)!= 0:
                for j in range(0,len(data.axvspanColor[0])):
                    axarr[i].axvspan(data.axvspan[j*2],data.axvspan[j*2+1], color=data.axvspanColor[0][j], alpha=0.25)
            
            ## create grid    
            if data.showGrid==True:
                axarr[i].grid()
                
                
            ## create yaxis
            if data.yLabel!=None:
                axarr[i].set_ylabel(data.yLabel)

            ## create y limit    
            if data.yLim!=[]:
                axarr[i].set_ylim(data.yLim)
                
            ## create x limit    
            if data.xLim!=[]:
                axarr[i].set_xlim(data.xLim)
                    
            ## create legend
            if data.yLegendLabel!=[]:
                axarr[i].legend(loc=data.legendLoc, fontsize=10)
                
            ## hide the top y tick label on each subplot
            axarr[i].get_yticklabels()[-1].set_visible(False)
            
            ## create x label
            if data.xLabel != []:
                if self.showOnlyBottomXLabel==True:
                    if i==m-1:
                        axarr[i].set_xlabel(data.xLabel);
                else:
                    axarr[i].set_xlabel(data.xLabel);
            
            ## create title
            if self.showOnlyTopTitle==True:
                if self.data[0].title != None:
                    axarr[0].set_title(self.data[0].title, verticalalignment='bottom');
                
            ## create subtitle
            if data.subtitle != None:
                axarr[i].annotate(data.subtitle, xy=(0.5, 0.98), xycoords='axes fraction', fontsize=14,
                    horizontalalignment='center', verticalalignment='top', bbox=dict(pad=5.0, facecolor="w"))
                
                
        
        
        ## horizontal spacing between subfigures
        f.subplots_adjust(hspace=self.hSpace)
        
        ## save figure
        if self.fileName != '':
            f.savefig(self.fileName+'.png')
        
            
            
            
def sPlot(xData,yData,xLim=None,yLim=None,xLabel=None,yLabel=None,title=None,marker=None,linestyle=None,color=None,legendLabel=None):
    """
    simple plot.  simple wrapper for the prePlot class that also plots itself
    """
    
    p1=prePlot();
    p1.yData=[yData]; 
    p1.xData=[xData]
    if xLim!=None:
        p1.xLim=xLim
    if yLim!=None:
        p1.yLim=xLim
    if xLabel!=None:
        p1.xLabel=xLabel
    if yLabel!=None:
        p1.yLabel=yLabel
    if title!=None:
        p1.title=title
    if marker!=None:
        p1.marker=marker
    if linestyle!=None:
        p1.linestyle=linestyle
    if color!=None:
        p1.color=color
    if legendLabel!=None:
        p1.yLegendLabel=legendLabel
    p1.plot()
    return p1

 
    
def _plotTAPStripey(TAPData):#,figNum,sbNum,ls, mk,color, rightIndex):
    """
    outdated but likely still works
    """
    Z= TAPData[0];
    [a,b]=_np.shape(Z);
    for i in range(0,a):
        Z[i,:]=Z[i,:]-_pd.boxCarSmoothing(Z[i,:],'flat',50)
    _plt.figure()
    theta=TAPData[3];
    t=TAPData[1][0, 1999:3000]*1000;
#    t=TAPData[1][0, :]*1000
    X, Y = _np.meshgrid(t, theta)
    Z= Z[:, 1999:3000];
    cp = _plt.contourf(X, Y, Z,100,cmap=_red_green_colormap())
    _plt.colorbar(cp)
#    _plt.title('Filled Contours Plot')
    _plt.xlabel('time [ms]')
    _plt.ylabel('TAp')
    _plt.show()
    
def _plotPA1Stripey(PA1Data):
    """
    outdated but likely still works
    """
    Z= PA1Data[0];
    [a,b]=_np.shape(Z);
    for i in range(0,a):
        Z[i,:]=Z[i,:]-_pd.boxCarSmoothing(Z[i,:],'flat',50)
    _plt.figure()
    theta=PA1Data[3];
    t=PA1Data[1][0, 1999:3000]*1000;
#    t=TAPData[1][0, :]*1000
    X, Y = _np.meshgrid(t, theta)
    Z= Z[:, 1999:3000];
    cp = _plt.contourf(X, Y, Z,100,cmap=_red_green_colormap())
    _plt.colorbar(cp)
#    _plt.title('Filled Contours Plot')
    _plt.xlabel('time [ms]')
    _plt.ylabel('PA1')
    _plt.show()
    
def _red_green_colormap():
    '''A colormap with a quick red-green transition at 0.5
    
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






def singleColorMapWithLowerAndUpperCutoffs(upperColor=(0.93, 0.51, 0.93),lowerColor=(1.,1.,1.),lowerCutoff=0.2,upperCutoff=0.8,plotExample=False):
    """
    Creates a simple linear colormap between two colors with a lower and upper cutoff.
    Colors in tuple format with values between 0.0 and 1.0. 
    Cutoffs also between 0.0 and 1.0.  
    
    Note: easy color creation
    c = mcolors.ColorConverter().to_rgb
    c('red')
    c('violet')
    c('white')
    # etc...
    
    # lowerColor=(0.98,0.98,0.98) # very light grey
    """
    import matplotlib.colors as mcolors

    def make_colormap(seq):
        """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
        seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)

    c = mcolors.ColorConverter().to_rgb
    vlg=(0.98,0.98,0.98) # very light grey
    rvb = make_colormap(
        [vlg, vlg, lowerCutoff, c('lightblue'), c('darkblue'),upperCutoff,c('black'),c('black')])
    if plotExample==True:
        N = 1000
        array_dg = _np.random.uniform(0, 10, size=(N, 2))
        colors = _np.random.uniform(0, 100, size=(N,))
        _plt.scatter(array_dg[:, 0], array_dg[:, 1], c=colors, cmap=rvb)
        _plt.colorbar()
        _plt.show()
    return rvb