"""
this library contains a number of useful functions that are general purpose 
data processing functions.
"""
            
###############################################################################
### import libraries
#from __init__ import (_np,_copy,_pd,_plt,_plot,_math)

    
### import libraries
            
# common libraries
import numpy as _np
import matplotlib.pyplot as _plt
import copy as _copy
import math as _math

# hbtepLib libraries
import _plotTools as _plot


            
###############################################################################
### misc functions
            
def convertDataToStairstepData(x,y):
    """
    When dealting with discrete data, often it makes sense to show the data 
    look like stair-steps insteady of a the more typical smooth plot.  This
    function makes this happen #TODO(John) there is a better way to explain 
    this...
    
    Example
    -------
    x=np.arange(0,1,0.01)
    y=np.sin(2*np.pi*12.345*x)
    (xOut,yOut)=hbt.process.convertDataToStaircaseDate(x,y)
    plt.plot(xOut,yOut)

    """
    
    xOut=_np.zeros(len(x)*2)
    yOut=_np.zeros(len(x)*2)
    dx=x[1]-x[0];
    for i in range(0,len(x)):
#        print yOut[i*2]
#        print y[i]
        xOut[2*i]=x[i]
        if i == len(x)-1:
            xOut[2*i+1]=xOut[2*i]+dx
        else:
            xOut[2*i+1]=x[i+1]
        yOut[2*i]=y[i]
        yOut[2*i+1]=y[i]
    return (xOut,yOut)
    
            
def findNearest(array,value):
    """
    search through `array` and returns the `index` of the cell closest to the 
    `value`.   `array` should be sorted in ascending order
    
    Parameters
    ----------
    array : numpy.array
        data array to search through
    value : float (or int)
        value to look for in array
        
    Return
    ------
    index : int
        index of value in array that is closest to value
    
    References
    ----------
    http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    index = (_np.abs(array-value)).argmin()
    # value = array[index] 
    return index 
    # return index, value # uncomment to return both the index AND the value
    
    
def rmse(data, targets=0):
    """
    Root mean square error function.  
    If targets = 0, this is also the root mean square function.  
    
    Parameters
    ----------
    data : numpy.array 
        Data being considered
    target : numpy.array of floats
        values being compared against
        
    Return
    ------
    : numpy.array 
        root mean square error of data
    
    References
    ----------
    # http://stackoverflow.com/questions/17197492/root-mean-square-error-in-python
    # http://statweb.stanford.edu/~susan/courses/s60/split/node60.html
    """
    return _np.sqrt(((data - targets) ** 2).mean())

    
def rms(data):
    """
    Root mean square function.   
    
    Parameters
    ----------
    data : numpy.array 
        Data to be processed
        
    Return
    ------
    : numpy.array 
        root mean square of data
    
    References
    ----------
    # http://stackoverflow.com/questions/17197492/root-mean-square-error-in-python
    # http://statweb.stanford.edu/~susan/courses/s60/split/node60.html
    """
    return _np.sqrt(((data - 0) ** 2).mean())
    
    
def rejectOutliers(data, sigma=2):
	"""
	remove outliers from set of data
	
	Parameters
	----------
	data : numpy.ndarray
		data array being considered
	sigma : int
		the number of std. devs. about which to reject data.  E.g. sigma=2 
		rejects outliers outside of +-2*sigma
		
	Return
	------
	 : numpy.ndarray 
		Same as databut missing entires considered as outliers
	indicesToKeep : numpy.ndarray (of bool)
		Boolean indices associated with entries of data that are kept
	
	References
	----------
	http://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
	"""
	indicesToKeep=abs(data - _np.mean(data)) < sigma* _np.std(data)
	return data[indicesToKeep],indicesToKeep
				    
def rmPhaseJumps(time,data,cut=_np.pi):
    """
    Removed phase jumps of greater than "cut", defaults to pi
    Assumes phase is wrapped.
    Inserts NaN at jump points in data and time vector 
    """
    indx=_np.argwhere(_np.abs(_np.diff(data))>=cut)+1
    if indx[-1]==len(data):indx[-1]-=1# in case last point jumps
    data=list(data)
    time=list(time)
    # Iterate over jump points, insert NaNs
    for i in indx:
        data.insert(i,_np.NaN)
        time.insert(i,_np.NaN)
        indx += 1
    return _np.array(time),_np.array(data)


def wrapPhase(data): 
    """
    Wraps phase data so that it repeats every 2pi.
    This is important for phase data when you want it to all fit nicely on a 
    plot.  
    
    Parameters
    ----------
    data : numpy.ndarray
        data being wrapped
        
    Return
    ------
    outData : numpy.ndarray
        wrapped data array
        
    """
    inData=data*1.0
    outData=_np.zeros(inData.size);
    inData-=_np.pi;
    for i in range(0,_np.size(inData)-1):
        a=_np.floor(inData[i]/(_np.pi*2))
        outData[i]=inData[i]-(a+1)*2*_np.pi+_np.pi
    return outData
    
    
def unwrapPhase(inData):
    """
    Takes in phase array (in radians).  I think it needs to be centered about 0.
    Unwraps phase data so that it is continuous.
    This is important for phase data when you want to take it's derivative to
    get frequency.  
    
    Parameters
    ----------
    data : numpy.ndarray
        data being unwrapped
        
    Return
    ------
    outData : numpy.ndarray
        unwrapped data array
        
    """
    outData=_np.zeros(_np.size(inData));
    offset=0;
    outData[0]=inData[0];
    for i in range(0,_np.size(inData)-1):
        if inData[i] > _np.pi/4 and inData[i+1] < -_np.pi/4:
            offset=offset+2*_np.pi;
        elif inData[i] < -_np.pi/4 and inData[i+1] > _np.pi/4:
            offset=offset-2*_np.pi;
        outData[i+1]=inData[i+1]+offset;
    return outData


def hasNan(inArray):
    """
    searches array, inArray, for any occurances of NaN.  returns True if
    yes, returns False if no.
    
    Parameters
    ----------
    inArray : numpy.ndarray
        data array being considered for NaN entries
        
    Return
    ------
        : bool
        True if NaNs are in array, False otherwise
    """
    count = 0;
    for i in range(0,len(inArray)):
        if _math.isnan(inArray[i]):
            count+=1;
            
    print("There was/were %d instances of NaN" % count)
    
    if count == 0:
        return False
    if count != 0:
        return True
        
        
def sort2Arrays(array1, array2):
    """
    sorts array1 and array2 in ascending order of array1
    
    outdated:  replaced by sortArrays() below
    """
    array1, array2 = zip(*sorted(zip(array1, array2)));
    array1=_np.array(array1);
    array2=_np.array(array2);
    return array1, array2
    
    
def sortArrays(arrays,sortIndex):
    """
    sorts a list of n arrays.  sortIndex is the index of the array to sort all arrays.
    
    example use:  [V, I, phi]=sortArrays([V,I,phi],1) to sort all arrays by ascending I
    
    reference:  https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
    """
    sortedIndices=arrays[sortIndex].argsort()
    for i in range(len(arrays)):
        arrays[i]=arrays[i][sortedIndices]
    return arrays
    
    
def downSampleData(downX,upX,data):
    """
    Down samples data by finding the nearest x values on data that match the 
    downselecting x (downX)
    
    Sometimes, you want to compare two different sets of data that do not have
    the same time basis.  This function down-samples the data set with more 
    data points so that its x-data matches the other
    
    Parameters
    ----------
    downX : numpy.ndarray
        the x-data that will provide the down smpaled reference to the 
        upsampled x-data
    upX : numpy.ndarray
        the upsampled x-data that will be downsampled
    data : list (of numpy.array)
        the upsampled y-data that will be downsampled
        
    Returns
    -------
    out : list (of np.ndarray)
        list of trimmed y-data
        
    Notes
    -----
    upX is not actually trimmed in this instance.  it is assumed that you user
    will use downX as their new time basis
    """
    if type(data) is not list:
        data=[data]
    
    m=len(downX)
    indices=_np.zeros(m,dtype=_np.int16)
    out = []
    
    for i in range(0,m):
        indices[i]=int(findNearest(upX,downX[i]))

    for i in range(0,len(data)):
        out.append(data[i][indices])
    return out
    
    
def upSampleData(upX,downX,data):
    """
    Up samples data by linear interpolating 
    
    Similar to downSampleData() but up-samples instead
    
    Parameters
    ----------
    upX : np.ndarray
        the x-data that will be used in up-sampling the under-sampled data
    downX : np.ndarray
        the undersampled x-data
    data : list (of np.ndarray)
        the y-data to be up-samples.  downX is its time-base before 
        up-sampling.  upX will be its time-base after up-sampling
    
    Returns
    -------
    out : list (of np.ndarray)
        list of up-sampled y-data
    """
    out = _np.interp(upX,downX,data)
    return out
    

    
    
def linearizeDataMatrix(data):
    """
    data is assumed to be a list of arrays
    
    this function converts the data to a single, appended array
    """
    m=len(data);
    temp=_np.array([])
    for i in range(0,m):
        temp=_np.append(temp,data[i])
    return temp
    
    

    
    
###############################################################################
### filters and smoothing algorithms
    

def downSample(data,time,new_dT,shftStart=False):
    if new_dT%_np.mean(_np.diff(time))>=1e-7:
        raise SyntaxError("Sampling rate must be integer multiple of original: %e vs %e"%(_np.mean(_np.diff(time)),new_dT))
    dS=(new_dT/_np.mean(_np.diff(time))).astype(int)
    
    # Verify data shape
    if data.shape[0] != len(time):data=data.T
    
    # Make smoothing kernel
    boxMat=_np.zeros((len(time)/dS,len(time)))
    for i in range(len(boxMat)):boxMat[i,i*dS:(i+1)*dS]=1./dS
    
    return _np.matmul(boxMat,data),_np.matmul(boxMat,time)-shftStart*new_dT*(3./2)#(dS/2)#*_np.mean(_np.diff(time))
    
def nPoleFilter(data,xData=None,numPoles=1,alpha=0.0625,filterType='lowPass',plot=False):
    """
    n-pole filter.  
    
    Parameters
    ----------
    data : numpy.ndarray
        data to be smoothed
    xData : numpy.ndarray or NoneType
        (optional) array of x-data
    numPoles : int
        number of poles for the filter
    alpha : float
        weight of the filter, float between 0 and 1.  Close to zero for a low
        pass and close to 1 for a high pass
    filterType : str
        'lowPass' - Low pass filter
        'highPhass' - High pass filter
    plot : bool
        plots results if true
    
    Returns
    -------
    filteredData : 2D numpy.ndarray
        filtered data
    
    References
    ----------
    http://techteach.no/simview/lowpass_filter/doc/filter_algorithm.pdf
    https://en.wikipedia.org/wiki/Low-pass_filter#Discrete-time_realization
    https://en.wikipedia.org/wiki/High-pass_filter#Discrete-time_realization
    
    Notes
    -----
    this method is pulled from Qian Peng's 2016 GPU code.  
    his highpass filter does NOT follow this code
    
    Example #1
    ----------
    t=np.arange(0,.01,6e-8);
    y2=np.sin(2*np.pi*300*t)+np.sin(2*np.pi*30000*t)
    alpha=0.00625;
    hbt.process.nPoleFilter(y2,t,numPoles=2,filterType='lowPass',plot=True,alpha=alpha)
    hbt.process.nPoleFilter(y2,t,numPoles=2,filterType='highPass',plot=True,alpha=1.0-0.000625)
    
    Example #2
    ----------
    t=np.arange(0,.01,6e-6);
    from scipy.signal import square
    y=square(t,0.001)
    hbt.process.nPoleFilter(y,t,numPoles=2,filterType='lowPass',plot=True)
    hbt.process.nPoleFilter(y,t,numPoles=2,filterType='highPass',plot=True)

    """
    
    # initialize data arrays
    procData=_np.zeros((numPoles+1,len(data)));
    procData[0,:]=data;
    
    # filter.  for loop controls the number of poles
    for i in range(1,numPoles+1):
        
        if filterType=='lowPass':
            for j in range(0,len(data)-1):
                procData[i,j+1]=procData[i,j]+alpha*(procData[i-1,j]-procData[i,j])
        elif filterType=='highPass':
            for j in range(0,len(data)-1):
                procData[i,j+1]=alpha*(procData[i,j]+procData[i-1,j+1]-procData[i-1,j])
                    
    # plot results
    if plot==True:
        p1=_plot.plot(title=str(numPoles)+" pole "+filterType+" filter, alpha="+str(alpha),
                   xLabel="x-axis",yLabel='y-axis')
                   
        if type(xData)==type(None):
            xData=_np.arange(0,len(data));
            
        p1.addTrace(xData=xData,yData=data,yLegendLabel='raw')

        for i in range(1,numPoles+1):
            p1.addTrace(xData=xData,yData=procData[i,:],yLegendLabel=str(i))
            
        p1.plot()

    # return results
    return procData
   
           
def savgolFilter(data,numPoints,polynomialOrder,plot=False):
    """
    Ssavitzky-Golay moving average smoothing filter.  Applies a nth order 
    polynomial to a moving window of data.
    
    Parameters
    ----------
    data : numpy.ndarray
        data to be smoothed
    numPoints : int
        number of points for smoothing
    method : str
        method to use
        'box' - box car type of smoothing
        'gaussian' - guassian or normal smoothing
    plot : bool
        plots results if true
    
    Returns
    -------
    smoothedData : numpy.ndarray
        smoothed data
    
    References
    ----------
    https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.signal.savgol_filter.html
    https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
    
    Notes
    -----
    this is a wrapper for scipy.signal.savgol_filter
    
    I think I prefer the gaussian convolution algorithm over this one (john)
    """
    
    def plotResults():
        """
        plots results (before and after smoothing)
        """
        p1=_plot.plot()
        x=_np.arange(0,len(data))
        p1.xData=[x,x]
        p1.yData=[data,smoothedData]
        p1.yLegendLabel=['raw data','smoothed data']
        p1.marker=['.','']
        p1.linestyle=['','-']
        p1.title='%d point, Savitzky-Golay smoothing of order %s'% (numPoints,polynomialOrder)
        p1.plot()
        
    # import savgol package
    from scipy.signal import savgol_filter
    
    # perform filter
    smoothedData=savgol_filter(data,numPoints,polynomialOrder)
        
    # plot results if requested
    if plot==True or plot=='all':
        plotResults()
        
    return smoothedData
        
    
    
def convolutionSmoothing(data,numPoints,method='gaussian',plot=False):
    """
    Convolution moving average smoothing filter
    
    Parameters
    ----------
    data : numpy.ndarray
        data to be smoothed
    numPoints : int
        number of points for smoothing.  should be an odd number
    method : str
        method to use
        'box' - box car type of smoothing.  keywords: boxcar
        'gaussian' - guassian or normal smoothing
    plot : bool
        plots results if true
    
    Returns
    -------
    smoothedData : numpy.ndarray
        smoothed data
    
    References
    ----------
    https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
    
    Example 1
    ---------
    # randon noise on sine way
    x = np.arange(0,2*np.pi,.1)
    y1 = np.sin(x) + np.random.random(len(x)) * 0.8
    convolutionSmoothing(y,numPoints,method='box',plot='all')
        
    Example 2
    ---------
    # step function
    x = np.arange(0,2*np.pi,.1)
    y2=np.zeros(len(x))
    y2[np.where(x>np.pi)[0]]=1
    convolutionSmoothing(y,numPoints,method='box',plot='all')

    Notes
    -----
    I want to justiy some of my weird code below.  When numPoints is an even 
    number, the filter centers itself 1.5 points ahead of itself in the window 
    which results in a small time/phase shift in the filtered data.  When it's 
    an odd number, it centers itself 1.0 points ahead in the window and 
    smaller time/phase shift.  By requiring that numPoints is odd and 
    temporarily removing the first point in the data (I put it back in the 
    end), I center the filter in the window and have no time/phase shift.  
    
    The convolution filters "mess up" the last set of points (on the order of
    numPoints).  If you want to get around this, you should give it more data 
    than you actually want and trim the end off.
    """
    def plotResults():
        """
        plots before and after of smoothing
        """
        p1=_plot.plot()
        x=_np.arange(0,len(data))
        p1.xData=[x,x]
        p1.yData=[data,smoothedData]
        p1.yLegendLabel=['raw data','smoothed data']
        p1.marker=['.','']
        p1.linestyle=['','-']
        p1.title='%d point, %s smoothing'% (numPoints,method)
        p1.plot()
        
    def plotSmoothingFunction():
        """
        plots smoothing function
        """
        p1=_plot.plot()
        x=_np.arange(0,len(smoothingFunction))
        p1.xData=[x]
        p1.yData=[smoothingFunction]
        p1.yLegendLabel=['smoothing function']
        p1.marker=['o']
        p1.linestyle=['-']
        p1.title='%d point, %s smoothing function'% (numPoints,method)
        p1.plot()
        
    # make sure numPoints is an odd number
    if (numPoints % 2 == 0): # is an even number
        numPoints+=1;
        print("Warning: numPoints was not an odd number.  +1 was added.  " 
              "numPoints is now %d" % numPoints)
        
    # temporarily remove first point (see Notes)
    tempPoint=data[0]
    data=data[1:]   
    
    # box smoothing
    if method=='box':
        smoothingFunction=_np.ones(numPoints)/numPoints
        
    # gaussian smoothing
    elif method=='gaussian' or method == 'normal':        
        temp=_np.arange(numPoints)
        sigma=numPoints/(2*_np.pi);
        smoothingFunction=_np.exp(-((temp-(numPoints-1)/2.)/sigma)**2/2)
        smoothingFunction/=_np.sum(smoothingFunction)
        
    # perform smoothing
    smoothedData = _np.convolve(data, smoothingFunction, mode='same')
    
    # plot if requested
    if plot==True:
        plotResults()
    if plot=='all':
        plotSmoothingFunction()
        plotResults()
        
    # add point back to smoothedData (see Notes
    smoothedData=_np.append(tempPoint,smoothedData)
        
    return smoothedData
    

def gaussianFilter(t,y,timeFWHM,filterType='high',plot=False,plotGaussian=False):
	"""
	Low and pass filters using scipy's gaussian convolution filter
	
	Parameters
	----------
	t : numpy.array
		time
	y : numpy.array
		time dependent data
	timeFWHM : float
		full width at half maximum of the gaussian with units in time.  this
		effectively sets the corner frequency of the filter
	filterType : str
		'high' - high-pass filter
		'low' - low-pass filter
	plot : bool
		plots the results
	plotGaussian : bool
		plots the gaussian distribution used for the filter
		
	Returns
	-------
	yFiltered : numpy.array
		filtered time dependent data
		
	References
	----------
	https://en.wikipedia.org/wiki/Full_width_at_half_maximum
	https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.signal.gaussian.html
	https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html
	"""
	
	dt=t[1]-t[0]
	
#	from scipy import signal
	
	from scipy.ndimage import gaussian_filter1d
	
	def fwhmToGaussFilterStd(fwhm,dt):
		
		std=1.0/_np.sqrt(8*_np.log(2))*fwhm/dt
		return std
	
	
	std=fwhmToGaussFilterStd(timeFWHM,dt)
#	yFiltered=signal.gaussian(len(t), std=std)
	
#	if filterType=='low':
	yFiltered=gaussian_filter1d(y*1.0,std,mode='nearest')
#	elif filterType=='high':
#		yFiltered=y-gaussian_filter1d(y*1.0,std)
	
	if plot==True:
		
		_plt.figure()
		_plt.plot(t,y,label='Raw')
		_plt.plot(t,yFiltered,label='Low-pass')
		_plt.plot(t,y-yFiltered,label='High-pass')
		_plt.legend()

	if plotGaussian==True:
		
		from scipy import signal
		_plt.figure()
		_plt.plot(t,signal.gaussian(len(t), std=std),label='gaussian')
		_plt.legend()
		
	if filterType=='low':
		return yFiltered
	else:
		return y-yFiltered

def gaussianLowPassFilter(y,t,timeWidth=1./20000,plot=False,plotGaussian=False):
	"""
	Low pass filter using scipy's gaussian filters
	
	Parameters
	----------
	y : numpy.array
		time dependent data
	t : numpy.array
		time
	timeWidth : float
		full width at half maximum of the gaussian with units in time.  this
		effectively sets the corner frequency of the filter
		(f_{corner} \approx 1/timeWidth)
	plot : bool
		plots the results
	plotGaussian : bool
		plots the gaussian distribution used for the filter
		
	Returns
	-------
	yFiltered : numpy.array
		filtered time dependent data
		
	References
	----------
	https://en.wikipedia.org/wiki/Full_width_at_half_maximum
	https://docs.scipy.org/doc/8scipy-0.19.0/reference/generated/scipy.signal.gaussian.html
	https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html
	
	Example
	-------
	import numpy as np
	t=np.arange(0,10e-3,2e-6)
	y = np.random.randn(len(t)).cumsum()
	y+=np.sin(2*np.pi*1000+np.pi*2*np.random.rand())
	y+=np.sin(2*np.pi*3300+np.pi*2*np.random.rand())
	y+=np.sin(2*np.pi*10000+np.pi*2*np.random.rand())
	y+=np.sin(2*np.pi*33000+np.pi*2*np.random.rand())
	gaussianLowPassFilter(y,t,timeWidth=1e-4,plot=True,plotGaussian=True)
	"""
	
	
	from scipy.ndimage import gaussian_filter1d

	dt=t[1]-t[0]
    #1/(dt*timeWidth*2*_np.pi)#
	sigma= (1./(2*_np.pi))*timeWidth/dt#2.355*timeWidth/dt#  #TODO(John)  This equation is wrong.  Should be dividing by 2.355, not multiplying.  Fix here and with all dependencies
	yFiltered=gaussian_filter1d(y,sigma)
	
	if plot==True:
		
		_plt.figure()
		_plt.plot(t,y,label='Raw')
		_plt.plot(t,yFiltered,label='Filtered')
		_plt.legend()
        _plt.grid()

	if plotGaussian==True:
		
		from scipy import signal
		_plt.figure()
		_plt.plot(t,signal.gaussian(len(t), std=sigma),label='gaussian')
		_plt.legend()
		
	return yFiltered


def gaussianHighPassFilter(y,t,timeWidth=1./20000,plot=False,plotGaussian=False):
	"""
	High pass filter using scipy's gaussian filters
	
	Parameters
	----------
	y : numpy.array
		time dependent data
	t : numpy.array
		time
	timeWidth : float
		full width at half maximum of the gaussian with units in time.  this
		effectively sets the corner frequency of the filter
		(f_{corner} \approx 1/timeWidth)
	plot : bool
		plots the results
	plotGaussian : bool
		plots the gaussian distribution used for the filter
		
	Returns
	-------
	yFiltered : numpy.array
		filtered time dependent data
		
	References
	----------
	https://en.wikipedia.org/wiki/Full_width_at_half_maximum
	https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.signal.gaussian.html
	https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html
	
	Example
	-------
	import numpy as np
	t=np.arange(0,10e-3,2e-6)
	y = np.random.randn(len(t)).cumsum()
	y+=np.sin(2*np.pi*1000+np.pi*2*np.random.rand())
	y+=np.sin(2*np.pi*3300+np.pi*2*np.random.rand())
	y+=np.sin(2*np.pi*10000+np.pi*2*np.random.rand())
	y+=np.sin(2*np.pi*33000+np.pi*2*np.random.rand())
	gaussianHighPassFilter(y,t,timeWidth=1./20000,plot=True,plotGaussian=True)
	"""
	fit=gaussianLowPassFilter(y,t,timeWidth,plot=False,plotGaussian=plotGaussian)
	yFiltered= y-fit
	
	if plot==True:
		
		_plt.figure()
		_plt.plot(t,y,label='Raw')
		_plt.plot(t,fit,label='Fit')
		_plt.plot(t,yFiltered,label='Filtered')
		_plt.legend()
        _plt.grid()
		
	return yFiltered, fit


    
def butterworthFilter(y, x,filterOrder=2, samplingRate=1/(2*1e-6), 
                      cutoffFreq=20*1e3, filterType='low',plot=False):
    """
    Apply a digital butterworth filter on your data
    
    Parameters
    ----------
    y : numpy.ndarray
        unfiltered dependent data
    x : numpy.ndarray
        independent data
    filterOrder : int
        Butterworth filter order
    samplingRate : float
        Data sampling rate.  
    cutoffFreq : float
        cutoff frequency for the filter
    filterType : str
        filter type.  'low' is lowpass filter
    plot : bool or str
        - True - plots filter results. 
        - 'all'- plots filter results and filter response (psuedo-BODE plot)
        
    Returns
    -------
    filteredData : numpy.ndarray
        Filtered dependent data
        
    References
    ----------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqz.html
    https://stackoverflow.com/questions/25191620/creating-lowpass-filter-in-scipy-understanding-methods-and-units
    """    
    
    
    from scipy.signal import butter, lfilter, freqz
    
    def butter_lowpass(cutoff, fs, order=5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = butter(order, normal_cutoff, btype=filterType, analog=False)
        return b, a
    
    def butter_lowpass_filter(data, cutoff, fs, order=5):
        b, a = butter_lowpass(cutoff, fs, order=order)
        y = lfilter(b, a, data)
        return y
        
    def plotOfFreqResponse():
        
        # Get the filter coefficients so we can check its frequency response.
        b, a = butter_lowpass(cutoffFreq, samplingRate, filterOrder)
        
        # calc frequency response 
        w, h = freqz(b, a, worN=8000)
        gain=_np.abs(h)
        phase=_np.unwrap(_np.angle(h))*180/_np.pi
        
        # generate gain plot        
        p1=_plot.plot()
        p1.xLim=[0, 0.5*samplingRate]
        p1.xLabel='Frequency [Hz]'
        p1.subtitle="Gain Response"
        p1.yLabel='Gain'
        p1.yLim=[0,1.]
        p1.title='Butterworth filter. Order=%d. Cutoff Freq=%.1f. Sampling Rate = %.1f Hz. Type = %s.' % (filterOrder, cutoffFreq, samplingRate, filterType)
        
        p1.xData.append(0.5*samplingRate*w/_np.pi)
        p1.yData.append(gain)
        p1.yLegendLabel.append('Frequency Response')
        
        p1.xData.append(_np.array([cutoffFreq,cutoffFreq]))
        p1.yData.append(_np.array([0,1])) #0.5*_np.sqrt(2)
        p1.yLegendLabel.append('Cuttoff Frequency')
        
        # generate phase plot        
        p2=_plot.plot()
        p2.xLim=[0, 0.5*samplingRate]
        p2.xLabel='Frequency [Hz]'
        p2.subtitle="Phase Response"
        p2.yLabel='Phase [Degrees]'
        p2.yLim=[_np.min(phase),0]
        
        p2.xData.append(0.5*samplingRate*w/_np.pi)
        p2.yData.append(phase)
        p2.yLegendLabel.append('Frequency Response')
        
        p2.xData.append(_np.array([cutoffFreq,cutoffFreq]))
        p2.yData.append(_np.array([_np.min(phase),0])) #0.5*_np.sqrt(2)
        p2.yLegendLabel.append('Cuttoff Frequency')
        
        # combine into subplot
        sp1=_plot.subPlot([p1,p2],plot=False)

        return sp1
        
    def plotOfResults():
        p1=_plot.plot()
        
        p1.xData.append(x)
        p1.yData.append(y)
        p1.yLegendLabel.append('Unfiltered Data')
        
        p1.xData.append(x)
        p1.yData.append(filteredData)
        p1.yLegendLabel.append('Filtered Data')
        
        return p1
        
        
    filteredData=butter_lowpass_filter(y, cutoffFreq,
                                       samplingRate, filterOrder)
                                       
    if plot==True:
        plotOfResults().plot()
       
    if plot == 'all':
        plotOfResults().plot()
        plotOfFreqResponse().plot()
        
    return filteredData
        
                
    
    
###############################################################################
### fitting functions and related
    
class polyFitData:
    """ 
    Polynomial fit function.  
    
    Parameters
    ----------
    yData : 'numpy.array'
        dependent variable
    xData : 'numpy.array'
        independent variable
    order : int
        order of polynomial fit.  1 = linear, 2 = quadratic, etc.
    plot : bool
        Causes a plot of the fit to be generated

    Attributes
    ----------
    coefs : 'numpy.array'
        array of fit coefficients, starts at highest order.  
    ffit : 'numpy.lib.polynomial.poly1d'
        function that returns yFit data given ANY numpy.array of x values
    fitData : 'numpy.array'
        yFit data corresponding to xData
    plotOfFit : 
        custom plot class of data. 
        
    Notes
    -----
    This is merely a wrapper function for the numpy.polyfit function.  However,
    it also plots the result automatically.  
    
    output:
    fitData is the y fit data.   
    """
    
    def __init__(self, yData, xData,order=2, plot=True):
        title = str(order)+' order Polynomial fit'
    
        ### do fit
        self.coefs=_np.polyfit(xData, yData, order)
        self.ffit = _np.poly1d(self.coefs)
        self.fitData=self.ffit(xData)
        
        ### generate plot        
        self.plotOfFit=_plot.plot();
        self.plotOfFit.xLabel='x'
        self.plotOfFit.yLabel='y'
        self.plotOfFit.title=title;
        
        ### raw data
        self.plotOfFit.xData.append(xData)
        self.plotOfFit.yData.append(yData)
        self.plotOfFit.marker.append('.')
        self.plotOfFit.linestyle.append('')
        self.plotOfFit.alpha.append(.15)
        self.plotOfFit.yLegendLabel.append('raw data')
        
        ### fit
        x=_np.linspace(_np.min(xData),_np.max(xData),1000);
        self.plotOfFit.xData.append(x)
        self.plotOfFit.yData.append(self.ffit(x))
        self.plotOfFit.marker.append('')
        self.plotOfFit.linestyle.append('-')
        self.plotOfFit.alpha.append(1.)
        self.plotOfFit.yLegendLabel.append('poly fit order %d'%order)
        
        if plot==True:
            self.plotOfFit.plot()
            
    
class genericCurveFit:
    """
    generic curve fitting function that uses scipy.optimize.curve_fit solution
    
    I "think" I like the genericLeastSquaresFit code better than this function.  See below.
    
    func = fit function
    indepVars = independent variables.  for multivariable, use indepVars = (x,y) etc.
    depVars = dependent variable.  this is the single dependent variable that we are trying to model
    guess = guess parameters.  use: guess = 8., 2., 7. etc.  
    
    note:  this function CAN be upgraded to include bounds
    
    references:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
    
    Example use
    -----------
    def dumbModel(variables,a,b,c):
        x,y=variables
        return _np.log(a) + b*_np.log(x) + c*_np.log(y)
    
    # some artificially noisy data to fit
    x = _np.linspace(0.1,1.1,101)
    y = _np.linspace(1.,2., 101)
    a, b, c = 10., 4., 6.
    z = dumbModel((x,y), a, b, c) * 1 + _np.random.random(101) / 100
    
    # initial guesses for a,b,c:
    p0 = 8., 2., 7.
    
    # solve
    indepVars=(x,y)
    a=genericCurveFit(dumbModel, indepVars, z, p0)

    """    

    def __init__(self,func, indepVars,depVar,guess, bounds=None,plot=True , maxNumIterations = None, fileName=''):
        self.indepVars=indepVars;
        self.depVars=depVar
        self.fileName=fileName
        
        
        from scipy.optimize import curve_fit
    
        if bounds == None:
            self.fitParams,self.covMatrix = curve_fit(func, indepVars, depVar, guess) #, max_nfev = 1e2 , max_nfev=maxNumIterations
        else:
            self.fitParams,self.covMatrix = curve_fit(func, indepVars, depVar, guess, bounds=bounds) #, max_nfev = 1e2 , max_nfev=maxNumIterations
#        print self.fitParams
#        print self.covMatrix
        self.fitSoln=func(indepVars, *self.fitParams)
#        print fitSoln

        self.R2=rSquared(indepVars, self.fitSoln)
        print(r"$R^2$ = %.3E" % self.R2)

        if plot==True:
            self.plot();
    
    def plot(self):
        _plt.figure()
        _plt.plot(self.fitSoln, self.depVars, '.', 
                  label='scipy.optimize.curve_fit solution',alpha=0.10)
#        _plt.plot(self.solnRobust, self.depData, 'x', label='robust least squares soln')
        _plt.plot(_np.array([_np.min(self.depVars),_np.max(self.depVars)]),
                  _np.array([_np.min(self.depVars),_np.max(self.depVars)]),
                  color='r',linestyle='--',linewidth=5)
        _plt.xlabel('Fit Solution')
        _plt.ylabel('Dependent Data')
        _plt.legend()    
        _plt.grid()
#        alphaB, alphaV, alphaR, phi0, C, C2 = self.fitParams
        #        _plt.ylim([-5,30])
        _plt.axes().set_aspect('equal') #, 'datalim'
        if self.fileName != '':
            _plt.savefig(self.fileName+'.png')


def _expFunction(x, a,b,c):
    """
    Basic exponential function.  Used primarily with fitting functions. 
    Output = a*_np.exp(x/b)+c
    
    Parameters
    ----------
    x : numpy.ndarray
        Independent variable
    a : float
        Fitting parameter.  Output = a*_np.exp(x/b)+c
    b : float
        Fitting parameter.  Output = a*_np.exp(x/b)+c
    c : float
        Fitting parameter.  Output = a*_np.exp(x/b)+c
        
    Returns
    -------
    : numpy.ndarray
        Output = a*_np.exp(x/b)+c
    
    """
    return a*_np.exp(x/b)+c
    
    
def _cosFunction(x, a,b,c,d):
    """
    Cosine function.  Used primarily with fitting functions.  
    Output = a*_np.cos(x*d*2*_np.pi+b)+c
    
    Parameters
    ----------
    x : numpy.ndarray
        Independent variable
    a : float
        Fitting parameter.  Output = a*_np.cos(x*d*2*_np.pi+b)+c
    b : float
        Fitting parameter.  Output = a*_np.cos(x*d*2*_np.pi+b)+c
    c : float
        Fitting parameter.  Output = a*_np.cos(x*d*2*_np.pi+b)+c
    d : float
        Fitting parameter.  Output = a*_np.cos(x*d*2*_np.pi+b)+c
        
    Returns
    -------
    : numpy.ndarray
        Output = a*_np.cos(x*d*2*_np.pi+b)+c
    """
    return a*_np.cos(x*d*2*_np.pi+b)+c
    

def singlePowerTerm(x, a,b,c):
    """
    Basic power term.  Used primarily with fitting functions.  
    Output = a*(x)**b+c
    
    Parameters
    ----------
    x : numpy.ndarray
        Independent variable
    a : float
        Fitting parameter.  Output = a*(x)**b+c
    b : float
        Fitting parameter.  Output = a*(x)**b+c
    c : float
        Fitting parameter.  Output = a*(x)**b+c
        
    Returns
    -------
    : numpy.ndarray
        Output = a*(x)**b+c
    """
    return a*(x)**b+c
    
    
class cosFit:
    """
    Cos fit function.  a*_np.cos(x*d*2*_np.pi+b)+c

    Parameters
    ----------
    y : numpy.ndarray
        dependent data
    x : numpy.ndarray
        independent array
    guess : list
        list of four floats [a, b, c, d]=[amplitude, phase offset, amplitude offest, linear frequency].  these are the guess values.
    plot : bool
        causes the results to be plotted
        
    Attributes
    ----------
    fit : genericLeastSquaresFit
    
    Notes
    -----
    if you are receiving the error: "ValueError: Residuals are not finite in 
    the initial point.", most likely, you need to play with your initial 
    conditions to get them closer to the right answer before the fit will work
    correctly.  
    
    Example use
    -----------
    # import library first.  I set it as hbt.pd 
    
    >>> y=np.array([11.622967, 12.006081, 11.760928, 12.246830, 12.052126, 12.346154, 12.039262, 12.362163, 12.009269, 11.260743, 10.950483, 10.522091,  9.346292,  7.014578,  6.981853,  7.197708,  7.035624,  6.785289, 7.134426,  8.338514,  8.723832, 10.276473, 10.602792, 11.031908, 11.364901, 11.687638, 11.947783, 12.228909, 11.918379, 12.343574, 12.046851, 12.316508, 12.147746, 12.136446, 11.744371,  8.317413, 8.790837, 10.139807,  7.019035,  7.541484,  7.199672,  9.090377,  7.532161,  8.156842,  9.329572, 9.991522, 10.036448, 10.797905])
    >>> x=np.linspace(0,2*np.pi,48)
    >>> c=hbt.process.cosFit(y,x,guess=[2,0,10,.3])
    # note that this example took me quite a bit of guessing with the guess 
    # values before everything fit correctly.

    """
    def __init__(self,y,x,guess,plot=True):
        self.fit=genericLeastSquaresFit(x=x,paramsGuess=guess,y=y, 
                                        function=_cosFunction,plot=plot)


class expFit:
    """
    Exponential fit function

    Parameters
    ----------
    y : numpy.ndarray
        dependent data
    x : numpy.ndarray
        independent array
    guess : list
        list of three floats.  these are the guess values.
    plot : bool
        causes the results to be plotted
        
    Attributes
    ----------
    fit : genericLeastSquaresFit
    
    Notes
    -----
    if you are receiving the error: "ValueError: Residuals are not finite in 
    the initial point.", most likely, you need to play with your initial 
    conditions to get them closer to the right answer before the fit will work
    correctly.  
    
    Example use
    -----------
    # import library first.  I set it as hbt.pd 
    >>> x = np.array([399.75, 989.25, 1578.75, 2168.25, 2757.75, 3347.25, 3936.75, 4526.25, 5115.75, 5705.25])
    >>> y = np.array([109,62,39,13,10,4,2,0,1,2])
    >>> hbt.pd.expFit(y,x,guess=[10,-100,1])

    """
    def __init__(self,y,x,guess=[1,1,1], plot=True):
        self.fit=genericLeastSquaresFit(x=x,paramsGuess=guess,y=y, 
                                        function=_expFunction,plot=plot)
        
    
class genericLeastSquaresFit:
    """
    Least squares fitting function(class)
    This is a wrapper for scipy.optimize.least_squares
    
    Parameters
    ----------
    y : numpy.ndarray
        dependent data
    x : numpy.ndarray (multidimensional)
        independent array(s).  multiple ind. variables are supported.  
    paramsGuess : list (of floats)
        guess values for the fit parameters
    function : function
        function to attempt to fit the data to.  see exdamples _cosFunction and 
        _expFunction to see how this function should be constructed
    yTrue : (optional) numpy.ndarray 
        if you know what the actual fit should be (for example when testing 
        this code against a known), include it here, and it will plot 
        alongside the other data.  also useful for debugging.  
    plot : bool
        causes the results to be plotted
        
    Attributes
    ----------
    fitParams : numpy.ndarray
        array of fit parameters, in the same order as the guess
    plotOfFit :
    plotOfFitDep :
        custom plot function of fit data plotted against dependent data.  this 
        is important if there is more than 1 dependent data array.
    rSquared : float
        r^2 result of the fit
    res : 
        fit output from scipy.optimize.least_squares
    yFit : numpy.ndarray
        y-fit data that corresponds with the independent data
    
    Notes
    -----
    i've found that the guess values often NEED to be somewhat close to the 
    actual values for the solution to converge
    
    i've implemented several specific functions that implement this function.
    see expFit and cosFit
    
    Example use
    ------------
    # define expoential function
    def _expFunction(x, a,b,c):
        return a*_np.exp(x/b)+c
    # generate noisy exponential signal
    x1=_np.linspace(-1,1,100);
    a=_np.zeros(len(x1));
    b=_np.zeros(len(x1));
    c=_np.zeros(len(x1));
    for i in range(0,len(x1)):
        a[i]=(random.random()-0.5)/4. + 1.
        b[i]=(random.random()-0.5)/4. + 1.
        c[i]=(random.random()-0.5)/4. + 1.
    y1=1+_np.pi*_np.exp(x1/_np.sqrt(2)) # actual solution
    y2=1*a+_np.pi*b*_np.exp(x1/_np.sqrt(2)/c) # noisy solution
    # perform fit
    d=genericLeastSquaresFit(x1,[1,1,1],y2, _expFunction, y1,plot=True)
    """           
    
    def __init__(self, x, paramsGuess, y, function,yTrue=[],plot=True ):

        def fit_fun(paramsGuess,x,y):
            return function(x, *paramsGuess) - y

        from scipy.optimize import least_squares
        
        self.x=x
        self.y=y
        self.yTrue=yTrue
        
        # perform least squares fit and record results
        self.res=least_squares(fit_fun, paramsGuess, args=[x,y])  #args=(y)
        self.yFit=function(x, *self.res.x) 
        self.fitParams=self.res.x;
        
        # calculate r^2
        self.rSquared=rSquared(y,self.yFit)
        
        # print results to screen
        print(r'R2 =  %.5E' % self.rSquared)
        print('fit parameters')
        print('\n'.join('{}: {}'.format(*k) for k in enumerate(self.res.x)))
        
        # plot data
        if plot==True:
            if type(x) is _np.ndarray or len(x)==1:
                self.plotOfFit().plot()
            self.plotOfFitDep().plot()
        
  
    
    def plotOfFit(self):
        """
        plots raw and fit data vs. its indep. variable.  
        
        Notes
        -----
        this function only works if there is a single indep. variable
        """            
        # make sure that there is only a single indep. variable
        if type(self.x) is _np.ndarray or len(self.x)==1:
            p1=_plot.plot();
            p1.yLabel='y'
            p1.xLabel='x'
            p1.title=r'Fit results.  R$^2$ = %.5f' % self.rSquared
            
            if isinstance(self.x,list):
                x=self.x[0]
            else:
                x=self.x
                
            # raw data
            p1.xData.append(x)
            p1.yData.append(self.y)
            p1.yLegendLabel.append('raw data')
            p1.marker.append('.')
            p1.linestyle.append('')
            p1.color.append('b')
            p1.alpha.append(0.3)
            
            # fit data
            p1.xData.append(x)
            p1.yData.append(self.yFit)
            p1.yLegendLabel.append('fit')
            p1.marker.append('')
            p1.linestyle.append('-')
            p1.color.append('r')
            p1.alpha.append(1.)
            
            # the true data (if applicable)
            if self.yTrue!=[]:
                p1.xData.append(x)
                p1.yData.append(self.yTrue)
                p1.yLegendLabel.append('True Soln')
                p1.marker.append('')
                p1.linestyle.append('-')
                p1.plotOfFit.color.append('k')
                p1.plotOfFit.alpha.append(1.)
                
            return p1
            
    def plotOfFitDep(self):
        """ 
        plot of fit vs dependent data.  important if there are multiple
        independent variables.
        """
        p1=_plot.plot();
        p1.yLabel='fit data'
        p1.xLabel='raw data'
        p1.aspect="equal"
        
        p1.xData.append(self.y)
        p1.yData.append(self.yFit)
        p1.yLegendLabel.append('actual fit')
        p1.marker.append('.')
        p1.linestyle.append('')
        p1.color.append('b')
        p1.alpha.append(0.3)
        
        p1.xData.append(_np.array([_np.min(self.y),_np.max(self.y)]))
        p1.yData.append(_np.array([_np.min(self.y),_np.max(self.y)]))
        p1.yLegendLabel.append('ideal fit line')
        p1.marker.append('')
        p1.linestyle.append('-')
        p1.color.append('r')
        p1.alpha.append(1.)
        p1.legendLoc=  'upper left'
        p1.title=r'Fit quality.  R$^2$ = %.5f' % self.rSquared
        
        return p1
        
        
def rSquared(y,f):
    """
    calculates R^2 of data fit
        
    Parameters
    ----------
    y : numpy.ndarray
        data being fit to, the dependent variable (NOT THE INDEPENDENT VARIABLE).  y is a functino of x, i.e. y=y(x)
    f : float
        fit data
        
    Returns
    -------
    : float 
        R^2 = 1 - \frac{\sum (f-y)^2 }{\sum (y-<y>)^2 }
    
    Reference
    ---------
    https://en.wikipedia.org/wiki/Coefficient_of_determination
    """
    yAve=_np.average(y);
    SSres = _np.sum( (y-f)**2 )
    SStot = _np.sum( (y-yAve)**2 )
    return 1-SSres/SStot
    
    
###############################################################################
### data management related

def listArrayToNumpyArray(inData):
    """
    Converts data from format a list of numpy.ndarrays to a 2D numpy.ndarray
    e.g. inData=[_np.array, _np.array, _np.array] to outData=_np.array([3,:])
    
    Parameters
    ----------
    inData : list (of numpy.ndarray)
        e.g. inData=[array([ 10.,  10.,  10.,  10.]),
                     array([ 15.,  15.,  15.,  15.]),
                     array([ 2.,  2.,  2.,  2.])]

    Returns
    -------
    outData : numpy.ndarray (2D)
        e.g. outData=   array([[ 10.,  15.,   2.],
                               [ 10.,  15.,   2.],
                               [ 10.,  15.,   2.],
                               [ 10.,  15.,   2.]])

    Notes
    -----
    -Note that this code requires that all arrays in the list have the same
    length
    -This code is so obvious that having a wrapper for it is kinda dumb...
    
    """
    outData = _np.array(inData)
    return outData
    
    
###############################################################################
### string manipulation related
    
def extractIntsFromStr(string):
    """
    Extracts all integers from a string.
    
    Parameters
    ----------
    string : str
        str with numbers embedded
        
    Returns
    -------
    numbers : list (of int)
        list of numbers that were within string
        
    Example
    -------
    print(extractNumsFromStr("123HelloMy65Is23"))
    
    Notes
    -----
    Does not work with decimal points.  Integers only.
    """
    import re
    
    # get list of numbers
    numbers=re.findall(r'\d+',string)
    
    # convert to integers
    for i in range(0,len(numbers)):
        numbers[i]=int(numbers[i])
        
    return numbers
     
