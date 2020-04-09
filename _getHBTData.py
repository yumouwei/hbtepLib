"""
load HBT data, both processed and unprocessed

NOTES
-----
The convention for units is to maintain everything in SI until they are
plotted.  

A few of these functions are merely wrappers for other people's code

In most of the functions below, I've specified default shotno's.  This is 
largely to make debugging easier as there is nothing special about the 
provided shotnos.  
"""

###############################################################################
### import libraries

#from __init__ import (_np,_mds,_copy,_sys,_socket,os,_pd,_plt,_plot,_time,_process,_pref)

### import libraries

# common libraries 
import numpy as _np
import MDSplus as _mds
from copy import copy as _copy
import sys as _sys
import _socket
import matplotlib.pyplot as _plt
import time as _time
import pandas as _pd

# hbtepLib libraries
import _processData as _process
import _plotTools as _plot
try:
    import _hbtPreferences as _pref
except ImportError:
    _sys.exit("Code hault: _hbtPreferences.py file not found.  See readme.md" +
    " concerning the creation of _hbtPreferences.py")

###############################################################################
### constants
#_REMOTE_DATA_PATH='/opt/hbt/data/control'  
if _socket.gethostname()==_pref._HBT_SERVER_NAME:
    _ON_HBTEP_SERVER=True;
else:
    _ON_HBTEP_SERVER=False;


###############################################################################
### global variables

# default time limits for data.  units in seconds
_TSTART = 0*1e-3;
_TSTOP  = 10*1e-3;

# list of known bad sensors.   likely VERY outdated.  also note that the code does not YET do anything with this info...
_SENSORBLACKLIST = ['PA1_S29R', 'PA1_S16P', 'PA2_S13R', 'PA2_S14P', 'PA2_S27P', 'FB03_S4R', 'FB04_S4R', 'FB08_S3P', 'FB10_S1P', 'FB04_S3P', 'FB06_S2P', 'FB08_S4P', 'TA07_S1P', 'TA02_S1P', 'TA02_S2P'];

# directory where unprocessed or minimally processed data is written locally.   
#_FILEDIR='/home/john/shotData/'


###############################################################################
### decorators
    
def _prepShotno(func,debug=False):
    """
    If no shot number is None or -1, -2, -3, etc, this decorator grabs the 
    latest shotnumber for whatever function is called
    
    References
    ----------
    https://stackoverflow.com/questions/2536307/decorators-in-the-python-standard-lib-deprecated-specifically/30253848#30253848
    
    Notes
    -----
    # TODO(John) Add a try/except error handling for bad shot numbers or 
    # missing data
    """ 
    from functools import wraps
    
    @wraps(func) # allows doc-string to be visible through the decorator function
    def inner1(*args, **kwargs):
        
        # check to see if shotno is an arg or kwarg.  if kwarg, effectively
        # move it to be an arg and delete the redundant kwarg key
        if len(args)>0:
            shotno=args[0]
        else:
            shotno=kwargs.get('shotno')
            del(kwargs['shotno'])
            
        # if shotno == None
        if debug==True:
            print('args = ')
            print(shotno)
            print(type(shotno))
            
        # check to see if it is a number (float,int,etc)
        if _np.issubdtype(type(shotno),_np.integer):
            # shotno = number if int(args[0]) does not throw an error
            int(shotno)
                        
            # if less than 0, use python reverse indexing notation to return the most recent shots
            if shotno<0:
                args=(latestShotNumber()+shotno+1,)+args[1:]
                waitUntilLatestShotNumber(args[0]) 
                return func(*args, **kwargs)
            
            # if a standard shot number (default case)
            else:
                # make sure the value is an integer
                args=(int(shotno),)+args[1:]
                waitUntilLatestShotNumber(int(shotno)) 
                return func(*args, **kwargs)
                
#        except ValueError:
#            # it must be a string
#            print("string... skipping...")
            
        else:
            # it might be None, a list, or an array
            
#            try:
                # try if it has a length (ie, it's either an array or list)
            n=len(shotno)
            
            out=[]
            for i in range(n):
                
                # if less than 0
                if shotno[i]<0:
                    arg=(latestShotNumber()+shotno[i]+1,)+args[1:]
                    waitUntilLatestShotNumber(int(arg[0])) 
                    out.append(func(*arg, **kwargs))
                    
                # if a standard shot number
                else:
                    arg=(shotno[i],)+args[1:]
                    waitUntilLatestShotNumber(int(arg[0])) 
                    out.append(func(*arg, **kwargs))
            return out
            
#            except TypeError:
#                # it must be None
#                
#                args=(latestShotNumber(),)+args[1:]
#                waitUntilLatestShotNumber(int(shotno)) 
#                return func(*args, **kwargs)
        
        
    return inner1


###############################################################################
### MDSplus tree data collection and misc. related functions
def _trimTime(time,data,tStart,tStop):
    """
    Trims list of data arrays down to desired time
    
    Parameters
    ----------
    time : numpy.ndarray
        time array
    data : list (of 1D numpy.ndarrays)
        list of data arrays to be trimmed
    tStart : float
        trims data before start time
    tStop : float
        trims data after stop time
        
    Returns
    -------
    time : numpy.ndarray
        trimmed time array
    data : list (of numpy.ndarrays)
        trimmed data
        
    Notes
    -----
    This function does not concern itself with units (e.g. s or ms). Instead, 
    it is assumed that tStart and tStop have the same units as the variable, time.  
    """    
    if tStart is None:
        iStart=0;
        iStop=len(time);
    else:
        # determine indices of cutoff regions
        iStart=_process.findNearest(time,tStart);   # index of lower cutoff
        iStop=_process.findNearest(time,tStop);     # index of higher cutoff
        
    # trim time
    time=time[iStart:iStop];
    
    # trim data
    if type(data) is not list:
        data=[data];
    for i in range(0,len(data)):
        data[i]=data[i][iStart:iStop];
        
    return time, data
    
    
def _initRemoteMDSConnection(shotno):
    """
    Initiate remote connection with MDSplus HBT-EP tree
    
    Parameters
    ----------
    shotno : int
    
    Returns
    -------
    conn : MDSplus.connection
        connection class to mdsplus tree
    """
    conn = _mds.Connection(_pref._HBT_SERVER_ADDRESS+':8003');
    conn.openTree('hbtep2', shotno);
    return conn
    
    
def latestShotNumber():
    """
    Gets the latest shot number from the tree
    
    Parameters
    ----------
    
    Returns
    -------
    shot_num : int
        latest shot number
    """
    conn = _mds.Connection(_pref._HBT_SERVER_ADDRESS+':8003');
    shot_num = conn.get('current_shot("hbtep2")')
    return int(shot_num)


def waitUntilLatestShotNumber(shotno,debug=False):
    """
    If the shotno that you are interested is the latest shotno,
    this code checks to see if all of the data has finished recording.  
    If not, the code pauses until it has.  Then it exits the code.
    
    This code is useful to include in anything where you want to make sure
    you aren't trying to get data from a shot number that hasn't finished recording yet.
    
    
    Parameters
    ----------
    shotno : int
        shot number
    debug : bool
        default False.  
        prints text to screen to help with debugging.  
    """

    latestShotno=latestShotNumber()
    if debug==True:
        print("latest shot number : %d"%latestShotno)
        print("shot number in question : %d"%shotno)
    
    # if you haven't made it to the latest number yet
    if shotno<latestShotno:
        return
    
    # if you are trying to access a number that hasn't even been created yet
    elif shotno>latestShotno:
        print("This shot number hasn't even been created yet.  Waiting...")
        stop=True
        while(stop==True):
            if shotno==latestShotNumber():
                stop=False
            else:
                pass
    
            _time.sleep(2)
        
    # if you are trying to access a number that has been created but not finished        
#    if shotno==latestShotno:
#    print("shotno==latestShotno")
    stop=False
    try:
        mdsData(shotno,dataAddress=['\HBTEP2::TOP.DEVICES.WEST_RACK:CPCI:INPUT_96'],
                    tStart=[],tStop=[])
    except _mds.TreeNODATA:
        stop=True
        print("Shot number has not finished.  Waiting...")
        
    while(stop==True):
        
        _time.sleep(2)
        try:
            mdsData(shotno,
                        dataAddress=['\HBTEP2::TOP.DEVICES.WEST_RACK:CPCI:INPUT_96'],
                        tStart=[],tStop=[])
            stop=False
            _time.sleep(1)
        except _mds.TreeNODATA:
            pass
                
    return

        
def mdsData(shotno=None,
            dataAddress=['\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_94',
                         '\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_95'],
            tStart=[],tStop=[]):
    """
    Get data and optionally associated time from MDSplus tree
    
    Parameters
    ----------
    shotno : int
        shotno of data.  this function will establish its own mdsConn of this 
        shotno
    dataAddress : list (of strings)
        address of desired data on MDSplus tree
    tStart : float
        trims data before this time
    tStop : float
        trims data after this time
    
    Returns
    -------
    data : list (of numpy.ndarray)
        requested data
    time : numpy.ndarray
        time associated with data array
    """            
        
    # convert dataAddress to a list if it not one originally 
    if type(dataAddress) is not list:
        dataAddress=[dataAddress];
        
#    # if shotno == -1, use the latest shot number
#    if shotno==-1:
#        shotno=latestShotNumber()
        
    # init arrays
    time = []
    data = []
        
    # check if computer is located locally or remotely.  The way it connects to spitzer remotely can only use one method, but locally, either method can be used.  
    if _ON_HBTEP_SERVER==True: # if operating local to the tree
        # converted from Ian's code
        
        tree = _mds.Tree('hbtep2', shotno)  
        for i in range(0,len(dataAddress)):

            node = tree.getNode(dataAddress[i])            #Get the proper node    
            data.append(node.data())                      #Get the data from this node 
        if type(data[0]) is _np.ndarray: # if node is an array, return data and time
            time = node.dim_of().data()        

    
    else: # operaeting remotely
    
        # if shotno is specified, this function gets its own mdsConn
        if type(shotno) is float or type(shotno) is int or type(shotno) is _np.int64:
            mdsConn=_initRemoteMDSConnection(shotno);

        for i in range(0,len(dataAddress)):
            data.append(mdsConn.get(dataAddress[i]).data())
        
        # if data is an array, also get time
        if type(data[0]) is _np.ndarray:
    
            time = mdsConn.get('dim_of('+dataAddress[0]+')').data();  # time assocated with data
            
        mdsConn.closeTree('hbtep2', shotno)
        mdsConn.disconnect()
 
    if time != [] and type(tStop)!=list:
        # trim time and data
        time,data= _trimTime(time,data,tStart,tStop)
        
    if time != []:
        return data, time
    else: 
        return data
    
    
###############################################################################
### get device specific data
    
@_prepShotno
class ipData:
    """
    Gets plasma current (I_p) data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    ip : numpy.ndarray
        plasma current data
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plotOfIP : 
        returns the plot of IP vs time
    plot :
        Plots all relevant plots
    
    """
    def __init__(    self,
                shotno=96530,
                tStart=_TSTART,
                tStop=_TSTOP,
                plot=False,
                findDisruption=True,
                verbose=False,
                paIntegrate=False):
        
        self.shotno = shotno
        self.title = "%d, Ip Data" % shotno
        
        # use the IP Rogowski coil to get IP
        data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.ROGOWSKIS:IP'],
                              tStart=tStart, tStop=tStop)
        
        self.ip=data[0];
        self.time=time;
            
        if paIntegrate==True:
            # integrate PA1 sensor data to get IP
            dfPA=paData(shotno,tStart,tStop).dfDataRaw
            dfPA=dfPA.drop(columns=['PA2_S27P','PA2_S14P','PA1_S16P'])
            dfPA1=dfPA.iloc[:,dfPA.columns.str.contains('PA1')]
            dfPA2=dfPA.iloc[:,dfPA.columns.str.contains('PA2')]
            
            mu0=4*_np.pi*1e-7
            minorRadius=0.16
            self.ipPA1Integration=_np.array(dfPA1.sum(axis=1)*1.0/dfPA1.shape[1]*2*_np.pi*minorRadius/mu0)
            self.ipPA2Integration=_np.array(dfPA2.sum(axis=1)*1.0/dfPA2.shape[1]*2*_np.pi*minorRadius/mu0)


        
        if findDisruption==True:
            
            try:
                time=self.time
                data=self.ip
                
                # only look at data after breakdown
                iStart=_process.findNearest(time,1.5e-3)
                ipTime=time[iStart:]
                ip=data[iStart:]
                
                # filter data to remove low-frequency offset
                ip2,temp=_process.gaussianHighPassFilter(ip,ipTime,
                                                timeWidth=1./20e3,plot=verbose)
                
                # find time derivative of ip2
                dip2dt=_np.gradient(ip2)                
                
                # find the first large rise in d(ip2)/dt
                threshold=40.0
                index=_np.where(dip2dt>threshold)[0][0]
                
                # debugging feature
                if verbose==True:
                    _plt.figure()
                    _plt.plot(ipTime,dip2dt,label=r'$\frac{d(ip)}{dt}$')
                    _plt.plot([ipTime[0],ipTime[-1]],[threshold,threshold],
                              label='threshold')
                    _plt.legend()
                
                # find the max value of ip immediately after the disrup. onset
                while(ip[index]<ip[index+1]):
                    index+=1
                self.timeDisruption=ipTime[index]
                
            except:
                print("time of disruption could not be found")
                self.timeDisruption=None
        
        if plot == True or plot=='all':
            self.plot()
        
        
    def plotOfIP(self):
        """
        returns the plot of IP vs time
        """
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.ip*1e-3,label='IP Rogowski')
        try:
            p1.plot(self.time*1e3,self.ipPA1Integration*1e-3,label='PA1')
            p1.plot(self.time*1e3,self.ipPA2Integration*1e-3,label='PA2')
        except:
            pass
        try:
            p1.plot(self.timeDisruption*1e3,self.ip[self.time==self.timeDisruption]*1e-3,label='Disruption',marker='x',linestyle='')
        except:
            "Time of disruption not available to plot"
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='Plasma Current (kA)')
        _plot.finalizeFigure(fig,title=self.title)
        
        return p1
        
            
    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfIP().plot()
        
        
@_prepShotno
def ipData_df(    shotno=96530,
                tStart=_TSTART,
                tStop=_TSTOP,
                plot=False,
                findDisruption=True,
                verbose=False,
                paIntegrate=True):
    """
    Gets plasma current (I_p) data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    ip : numpy.ndarray
        plasma current data
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plotOfIP : 
        returns the plot of IP vs time
    plot :
        Plots all relevant plots
    
    """
    
    shotno = shotno
    title = "%d, Ip Data" % shotno
    
    # use the IP Rogowski coil to get IP
    data, time=mdsData(shotno=shotno,
                          dataAddress=['\HBTEP2::TOP.SENSORS.ROGOWSKIS:IP'],
                          tStart=tStart, tStop=tStop)
    
    ip=data[0];
    time=time;
    
    dfData=_pd.DataFrame(index=time)
    dfData['ip']=ip*1.0
    
#     dfData.plot()
    
    if paIntegrate==True:
        # integrate PA1 sensor data to get IP
        dfPA=paData(shotno,tStart,tStop).dfDataRaw
        key='PA1'
        dfPA=dfPA.iloc[:,dfPA.columns.str.contains(key)]
        
        mu0=4*_np.pi*1e-7
        minorRadius=0.16
        ipPAIntegration=_np.array(dfPA.sum(axis=1)*1.0/dfPA.shape[1]*2*_np.pi*minorRadius/mu0)
        dfData['ip%s'%key]=ipPAIntegration

    
    if findDisruption==True:
        
        try:
            dfTemp=dfData[dfData.index>1.5e-3]
            
            # filter data to remove low-frequency offset
            ipSmooth,temp=_process.gaussianHighPassFilter(dfTemp.ip.to_numpy()*1.0,dfTemp.index.to_numpy(),
                                            timeWidth=1./20e3,plot=verbose)
            
            # find time derivative of smoothed ip
            dip2dt=_np.gradient(ipSmooth)                
            
            # find the first large rise in d(ip2)/dt
            threshold=15.0
            index=_np.where(dip2dt>threshold)[0][0]
            
            # debugging feature
            if verbose==True:
                _plt.figure()
                t=dfTemp.index.to_numpy()
                _plt.plot(t,dip2dt,label=r'$\frac{d(ip)}{dt}$')
                _plt.plot([t[0],t[-1]],[threshold,threshold],
                          label='threshold')
                _plt.legend()
            
            # find the max value of ip immediately after the disrup. onset
            while(dfTemp.ip.to_numpy()[index]<dfTemp.ip.to_numpy()[index+1]):
                index+=1
            tDisrupt=dfTemp.iloc[index].name
            timeDisruption=_np.zeros(dfData.shape[0])
            dfData['timeDisruption']=timeDisruption
            dfData['timeDisruption'].at[tDisrupt]=dfData['ip'].at[tDisrupt]
#             plt.plot(timeDisruption)
            
        except:
            print("time of disruption could not be found")
            timeDisruption=None
    
        
    def plot():
        """
        returns the plot of IP vs time
        """
        fig,p1=_plt.subplots()
        p1.plot(time*1e3,ip*1e-3,label='IP Rogowski')
        try:
            p1.plot(time*1e3,ipPAIntegration*1e-3,label='PA')
        except:
            pass
        try:
            p1.plot(timeDisruption*1e3,ip[time==timeDisruption]*1e-3,label='Disruption',marker='x',linestyle='')
        except:
            "Time of disruption not available to plot"
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='Plasma Current (kA)')
        _plot.finalizeFigure(fig,title=title)
        
        return p1
    
    if plot == True or plot=='all':
        plot()
        
    return dfData
        
            
        
        
@_prepShotno
class egunData:
    """
    Gets egun data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    heatingCurrent : numpy.ndarray
        egun heating current
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plot :
        Plots all relevant plots
    
    """
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno = shotno
        self.title = "%d, egun Data" % shotno
        
        # get data
        data,time=mdsData(101169,
              dataAddress=['\HBTEP2::TOP.OPER_DIAGS.E_GUN:I_EMIS', 
                        '\HBTEP2::TOP.OPER_DIAGS.E_GUN:I_HEAT',
                        '\HBTEP2::TOP.OPER_DIAGS.E_GUN:V_BIAS',])
        self.I_EMIS=data[0]
        self.heatingCurrent=data[1]
        self.biasVoltage=data[2]
        self.heatingCurrentRMS=_np.sqrt(_np.average((self.heatingCurrent-_np.average(self.heatingCurrent))**2));
        self.time=time;
        
        if plot == True or plot=='all':
            self.plot()
        
        
    def plotOfHeatingCurrent(self):
        """
        returns the plot of heating current vs time
        """
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.heatingCurrent)
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='Current (A)')
        _plot.finalizeFigure(fig,title=self.title)
        
        return p1
            
    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfHeatingCurrent().plot()
        
        
@_prepShotno
class cos1RogowskiData:
    """
    Gets cos 1 rogowski data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    cos1 : numpy.ndarray
        cos1 data
    time : numpy.ndarray
        time data
    cos2Raw : numpy.ndarray
        raw cos1 data
        
    Subfunctions
    ------------
    plotOfIP : 
        returns the plot of IP vs time
    plot :
        Plots all relevant plots
    
    Notes
    -----
    This function initially grabs data starting at -1 ms.  This is because it 
    needs time data before 0 ms to calculate the cos1RawOffset value.  After 
    calculating this value, the code trims off the time before tStart.
    """
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno = shotno
        self.title = "%d, Cos1 Rog. Data" % shotno
        
        # get data.  need early time data for offset subtraction
        data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.ROGOWSKIS:COS_1',
                                           '\HBTEP2::TOP.SENSORS.ROGOWSKIS:COS_1:RAW'],
                              tStart=-1*1e-3, tStop=tStop)
        
        # calculate offset
        self.cos1Raw=data[1]
        indices=time<0.0*1e-3
        self.cos1RawOffset=_np.mean(self.cos1Raw[indices])
        
        # trime time before tStart
        iStart=_process.findNearest(time,tStart)
        self.cos1=data[0][iStart:];
        self.time=time[iStart:];
        self.cos1Raw=self.cos1Raw[iStart:]
        
        if plot == True or plot=='all':
            self.plot()
        
        
    def plotOfCos1(self):
        """
        returns the plot of cos1 rog vs time
        """
#        p1=_plot.plot(yLabel='',xLabel='time [ms]',title=self.title,
#                      subtitle='Cos1 Rogowski',shotno=[self.shotno])
#        p1.addTrace(xData=self.time*1000,yData=self.cos1)
#        
#        return p1
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.cos1)
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='')
        _plot.finalizeFigure(fig,title=self.title)
        
        return p1
            
    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfCos1().plot()
        
  
@_prepShotno
class bpData:
    """
    Downloads bias probe data from both probes.  
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before \line
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    ip : numpy.ndarray
        plasma current data
    time : numpy.ndarray
        time data
    title : str
        title of all included figures
    bps9Voltage : numpy.ndarray
        bps9 voltage
    bps9Current : numpy.ndarray
        bps9 current
    bps5Voltage : numpy.ndarray
        bps9 voltage
    bps5Current : numpy.ndarray
        bps9 current
    bps9GPURequestVoltage : numpy.ndarray
        CPCI measurement of pre-amp voltage, out from the GPU, and going to 
        control bps9
        
    Subfunctions
    ------------
    plotOfGPUVoltageRequest : 
        Plot of gpu request voltage (as measured by the CPCI)
    plotOfVoltage : 
        Plot of both bias probe voltages
    plotOfCurrent : 
        Plot of both bias probe currents
    plotOfBPS9Voltage : 
        Plot of bps9 voltage only
    plotOfBPS9Current : 
        Plot of bps9 current only
    plot :
        Plots all relevant plots
        
    Notes
    -----
    BPS5 was moved to section 2 (now BPS2) summer of 2017.  Instrumented on May 22, 2018.  
    
    """
    # TODO(John) Time should likely be split into s2 and s9 because different 
    # racks often have slightly different time bases
    # TODO(John) Determine when the Tree node for the BP was added, and have
    # this code automatically determine whether to use the old or new loading
    # method
    # TODO(John) Also, one of the BPs was moved recently.  Need to figure out
    # how to handle this
    # TODO(John) these probes have been periodically moved to different nodes.  
    # implement if lowerbound < shotno < upperbound conditions to handle these cases
    def __init__(self,shotno=98147,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno = shotno
        self.title = "%s, BP Data." % shotno
        
        if shotno > 99035 or shotno==-1:
            # BPS5 was moved to section 2
            
            # get voltage data
            data, time=mdsData(shotno=shotno,
                               dataAddress=['\HBTEP2::TOP.SENSORS.BIAS_PROBE_9:VOLTAGE',
                                            '\HBTEP2::TOP.SENSORS.BIAS_PROBE_9:CURRENT'],
                               tStart=tStart, tStop=tStop)
            self.bps9Voltage=data[0];
            self.bps9Current=data[1];#r*-1; # signs are flipped somewhere
            self.time=time;
            
            # get current data
            try:
                data, time=mdsData(shotno=shotno,
                                   dataAddress=['\HBTEP2::TOP.SENSORS.BIAS_PROBE_2:VOLTAGE',
                                                '\HBTEP2::TOP.SENSORS.BIAS_PROBE_2:CURRENT'],
    #                               dataAddress=['\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_85',
    #                                            '\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_84'],
                                   tStart=tStart, tStop=tStop)
                self.bps2Voltage=data[0]#*100; #TODO get actual voltage divider info
                self.bps2Current=data[1]*-1#/0.05;  
            except:
                "no bp2"
                

        elif shotno > 96000 and shotno < 99035 :
            #TODO(determine when this probe was rewired or moved)
        
            # get voltage data
            data, time=mdsData(shotno=shotno,
                               dataAddress=['\HBTEP2::TOP.SENSORS.BIAS_PROBE_9:VOLTAGE',
                                            '\HBTEP2::TOP.SENSORS.BIAS_PROBE_9:CURRENT'],
                               tStart=tStart, tStop=tStop)
            self.bps9Voltage=data[0];
            self.bps9Current=data[1]*-1; # signs are flipped somewhere
            self.time=time;
            
            # get current data
            data, time=mdsData(shotno=shotno,
                               dataAddress=['\HBTEP2::TOP.SENSORS.BIAS_PROBE_5:VOLTAGE',
                                            '\HBTEP2::TOP.SENSORS.BIAS_PROBE_5:CURRENT'],
                               tStart=tStart, tStop=tStop)
            self.bps5Voltage=data[0];
            self.bps5Current=data[1];
            
            ## previous BP addresses.  do not delete this until implemented
            # if probe == 'BPS5' or probe == 'both':
            #     self.currentBPS5=conn.get('\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_83').data()/.01/5;
            #     self.voltageBPS5 = conn.get('\HBTEP2::TOP.DEVICES.NORTH_RACK:CPCI:INPUT_82').data()*80;
            # if probe == 'BPS9' or probe == 'both':
            #     self.timeBPS9 = conn.get('dim_of(\TOP.DEVICES.SOUTH_RACK:A14:INPUT_3)').data();
            #     self.voltageBPS9 = (-1.)*conn.get('\TOP.DEVICES.SOUTH_RACK:A14:INPUT_4').data()/.00971534052268532 / 1.5

        else:
                        # get voltage data
            data, time=mdsData(shotno=shotno,
                               dataAddress=['\HBTEP2::TOP.SENSORS.BIAS_PROBE:VOLTAGE',
                                            '\HBTEP2::TOP.SENSORS.BIAS_PROBE:CURRENT'],
                               tStart=tStart, tStop=tStop)
            self.bps9Voltage=data[0];
            self.bps9Current=data[1]*-1; # signs are flipped somewhere
            self.time=time;
            
            # get current data
            try:
                data, time=mdsData(shotno=shotno,
                                   dataAddress=['\HBTEP2::TOP.SENSORS.BIAS_PROBE_2:VOLTAGE',
                                                '\HBTEP2::TOP.SENSORS.BIAS_PROBE_2:CURRENT'],
                                   tStart=tStart, tStop=tStop)
                self.bps5Voltage=data[0];
                self.bps5Current=data[1];
            except:
                print("skipping bps5")
            
        # transformer primary voltage.  first setup for shot 100505 and on.  
        [primaryVoltage,primaryCurrent], time=mdsData(shotno=shotno,
                                                     dataAddress=['\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_86',
                                                                '\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_87'],
                                                      tStart=tStart, tStop=tStop)
        self.primaryVoltage=primaryVoltage*(0.745/(110+.745))**(-1) # correct for voltage divider
        self.primaryCurrent=primaryCurrent*0.01**(-1) # correct for Pearson correction factor
#        self.primaryCurrent*=1; #the sign is wrong.  
#        self.primaryVoltage*=-1; #the sign is wrong.  
        
        # get gpu request voltage (for when the BP is under feedforward or feedback control)
        data, time=mdsData(shotno=shotno,
                           dataAddress=['\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_93'],
                           tStart=tStart, tStop=tStop)
        self.bps9GPURequestVoltage=data[0];
        
        if plot==True:
            self.plot()
        if plot=='all':
            self.plot(True)

    def plotOfGPUVoltageRequest(self):
        """ 
        returns plot of gpu voltage request 
        (Preamp signal out from caliban)            
        """
        p1=_plot.plot(title=self.title,xLabel='ms',yLabel='V',
                      subtitle='Voltage Request from GPU (pre-amplifier)',
                      shotno=self.shotno);
        p1.addTrace(xData=self.time*1000,yData=self.bps9GPURequestVoltage,
                    yLegendLabel='BPS9')
        return p1
    
#    def plotOfPrimaryVoltage(self):
#        """ 
#        returns plot of transformer primary values
#        (Preamp signal out from caliban)            
#        """
#        p1=_plot.plot(title=self.title,xLabel='ms',yLabel='V',
#                      subtitle='Primary voltage',
#                      shotno=self.shotno);
#        p1.addTrace(xData=self.time*1000,yData=self.primaryVoltage,
#                    yLegendLabel='')
#        return p1
#    
#    def plotOfPrimaryCurrent(self):
#        """ 
#        returns plot of transformer primary values
#        (Preamp signal out from caliban)            
#        """
#        p1=_plot.plot(title=self.title,xLabel='ms',yLabel='A',
#                      subtitle='Primary current',
#                      shotno=self.shotno);
#        p1.addTrace(xData=self.time*1000,yData=self.primaryCurrent,
#                    yLegendLabel='')
#        return p1
        
        
    def plotOfVoltage(self,primary=False):
        """ 
        returns plot of BP voltage          
        """
        p1=_plot.plot(title=self.title,yLabel='V', #yLim=[-200,200]
                      xLabel='Time [ms]',subtitle='BP Voltage',
                      shotno=[self.shotno]);
        p1.addTrace(xData=self.time*1000,yData=self.bps9Voltage,
                    yLegendLabel='BPS9')
        try:
            p1.addTrace(xData=self.time*1000,yData=self.bps2Voltage,
                    yLegendLabel='BPS2')
        except:
            p1.addTrace(xData=self.time*1000,yData=self.bps5Voltage,
                    yLegendLabel='BPS5')
        if primary==True:
            p1.addTrace(xData=self.time*1000,yData=self.primaryVoltage,
                    yLegendLabel='Primary')
        
        return p1
        
    def plotOfCurrent(self,primary=False):
        """ 
        returns plot of BP current         
        """
        p1=_plot.plot(title=self.title,yLabel='A',
                      xLabel='Time [ms]',subtitle='BP Current',
                      shotno=[self.shotno])
        p1.addTrace(xData=self.time*1000,yData=self.bps9Current,
                    yLegendLabel='BPS9')
        try:
            p1.addTrace(xData=self.time*1000,yData=self.bps2Current,
                        yLegendLabel='BPS2')
        except:
            p1.addTrace(xData=self.time*1000,yData=self.bps5Current,
                        yLegendLabel='BPS5')
        if primary==True:
            p1.addTrace(xData=self.time*1000,yData=self.primaryCurrent,
                    yLegendLabel='Primary')
        
        return p1
            
    def plotOfBPS9Voltage(self):
        """ 
        returns plot of BPS9 voltage         
        """
        p1=_plot.plot(title=self.title,yLabel='V',yLim=[-200,200],
                      xLabel='Time [ms]',subtitle='BP Voltage',
                      shotno=[self.shotno])
        p1.addTrace(xData=self.time*1000,yData=self.bps9Voltage,
                    yLegendLabel='BPS9')
        
        return p1
    
    def plotOfBPS9Current(self):
        """ 
        returns plot of BPS9 current         
        """
        p1=_plot.plot(title=self.title,yLabel='A',
                      xLabel='Time [ms]',subtitle='BP Current',
                      shotno=[self.shotno])
        p1.addTrace(xData=self.time*1000,yData=self.bps9Current,
                    yLegendLabel='BPS9')
        
        return p1
    
        # TODO(john) also make plots for BPS5 only

    def plot(self,plotAll=False):
        """ Plot relevant plots """
        if plotAll==False:
            sp1=_plot.subPlot([self.plotOfVoltage(),self.plotOfCurrent()])
        else:
            sp1=_plot.subPlot([self.plotOfVoltage(True),self.plotOfCurrent(True),self.plotOfGPUVoltageRequest()])
        return sp1
    
    
@_prepShotno
class dpData:
    """
    Downloads directional (double) probe data f
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before \line
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    ip : numpy.ndarray
        plasma current data
    time : numpy.ndarray
        time data
    title : str
        title of all included figures
    dp1Voltage : numpy.ndarray
        double probe 1 voltage
    dp1Current : numpy.ndarray
        double probe 1 current
    gpuRequestVoltage : numpy.ndarray
        CPCI measurement of pre-amp voltage, out from the GPU, and going to 
        the probe
        
    Subfunctions
    ------------
    #TODO
        
    Notes
    -----  
    
    """
    def __init__(self,shotno,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno = shotno
        self.title = "%s, DP Data." % shotno
                
        # get probe results
        data, time=mdsData(shotno=shotno,
                               dataAddress=['\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_85',
                                            '\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_84',
                                            '\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_83'],
                               tStart=tStart, tStop=tStop)
        self.dp1VoltageLeft=data[0]*(469.7/(469.7+100000))**(-1);
        self.dp1Current=data[1]*0.1**(-1)
        self.dp1VoltageRight=data[2]*(470./(470+99800))**(-1);
        self.dp1VoltageDiff=self.dp1VoltageLeft-self.dp1VoltageRight
        self.time=time;
#        self.dp1Current*=-1;
                    
        # transformer primary voltage.  first setup for shot 100505 and on.  
        [primaryVoltage,primaryCurrent,secondaryVoltage], time=mdsData(shotno=shotno,
                                                     dataAddress=['\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_86',
                                                                '\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_87',
                                                                '\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_88'],
                                                      tStart=tStart, tStop=tStop)
#        self.primaryVoltage=primaryVoltage*(0.745/(110+.745))**(-1) # correct for voltage divider
        self.primaryVoltage=primaryVoltage*(507.2/102900.)**(-1) # correct for voltage divider
        self.primaryCurrent=primaryCurrent*0.01**(-1) # correct for Pearson correction factor
#        self.primaryCurrent*=-1
        self.secondaryVoltage=secondaryVoltage*(514.6/102500)**(-1)
        
        # get gpu request voltage (for when the BP is under feedforward or feedback control)
        data, time=mdsData(shotno=shotno,
                           dataAddress=['\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_93'],
                           tStart=tStart, tStop=tStop)
        self.gpuRequestVoltage=data[0];
        
        if plot==True:
            self.plot()
        if plot=='all':
            self.plot(True)

    def plotOfGPUVoltageRequest(self):
        """ 
        returns plot of gpu voltage request 
        (Preamp signal out from caliban)            
        """
        p1=_plot.plot(title=self.title,xLabel='ms',yLabel='V',
                      subtitle='Voltage Request from GPU (pre-amplifier)',
                      shotno=self.shotno);
        p1.addTrace(xData=self.time*1000,yData=self.gpuRequestVoltage)#,
                    #yLegendLabel='BPS9')
        return p1
    
#    def plotOfPrimaryVoltage(self):
#        """ 
#        returns plot of transformer primary values
#        (Preamp signal out from caliban)            
#        """
#        p1=_plot.plot(title=self.title,xLabel='ms',yLabel='V',
#                      subtitle='Primary voltage',
#                      shotno=self.shotno);
#        p1.addTrace(xData=self.time*1000,yData=self.primaryVoltage,
#                    yLegendLabel='')
#        return p1
#    
#    def plotOfPrimaryCurrent(self):
#        """ 
#        returns plot of transformer primary values
#        (Preamp signal out from caliban)            
#        """
#        p1=_plot.plot(title=self.title,xLabel='ms',yLabel='A',
#                      subtitle='Primary current',
#                      shotno=self.shotno);
#        p1.addTrace(xData=self.time*1000,yData=self.primaryCurrent,
#                    yLegendLabel='')
#        return p1
        
        
    def plotOfVoltage(self,primary=False):
        """ 
        returns plot of DP voltage          
        """
        p1=_plot.plot(title=self.title,yLabel='V', #yLim=[-200,200]
                      xLabel='Time [ms]',subtitle='DP Voltage',
                      shotno=[self.shotno]);
        p1.addTrace(xData=self.time*1000,yData=self.dp1VoltageRight,
                    yLegendLabel='DP1 Right')
        p1.addTrace(xData=self.time*1000,yData=self.dp1VoltageLeft,
                    yLegendLabel='DP1 Left')
        p1.addTrace(xData=self.time*1000,yData=self.dp1VoltageDiff,
                    yLegendLabel='DP1 Difference')
        if primary==True:
            p1.addTrace(xData=self.time*1000,yData=self.primaryVoltage,
                    yLegendLabel='Primary')
            p1.addTrace(xData=self.time*1000,yData=self.secondaryVoltage,
                    yLegendLabel='Secondary')
        
        return p1
        
    def plotOfCurrent(self,primary=False):
        """ 
        returns plot of DP current         
        """
        p1=_plot.plot(title=self.title,yLabel='A',
                      xLabel='Time [ms]',subtitle='DP Current',
                      shotno=[self.shotno])
        p1.addTrace(xData=self.time*1000,yData=self.dp1Current,
                    yLegendLabel='DP1')
        if primary==True:
            p1.addTrace(xData=self.time*1000,yData=self.primaryCurrent,
                    yLegendLabel='Primary')
        
        return p1
    
    def plot(self,plotAll=False):
        """ Plot relevant plots """
        if plotAll==False:
            sp1=_plot.subPlot([self.plotOfVoltage(),self.plotOfCurrent()])
        else:
            sp1=_plot.subPlot([self.plotOfVoltage(True),self.plotOfCurrent(True),self.plotOfGPUVoltageRequest()])
        return sp1
        
    
@_prepShotno
class tpData:
    """
    Triple probe data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
    probes : str
        This parameter allows the user to specify which probe from which to 
        load data.  There are two triple probes: tps5 (triple probe section 5) 
        and tps8.  This str can be 'tps5', 'tps8', or 'both'.  
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title of all included figures
    self.tps5TipA : numpy.ndarray
        tps5 tip A voltage data.  (note that this channel is typically 
        disconnected)
    tps5TipB : numpy.ndarray
        tps5 tip B voltage data.  
    tps5TipC : numpy.ndarray
        tps5 tip C voltage data.  
    tps5Time : numpy.ndarray
        tps5 time data 
    tps5Current : numpy.ndarray
        tps5 current data
    tps5Temp : numpy.ndarray    
        tps5 temperature data.  
    tps5VFloat : numpy.ndarray
        tps5 floating voltage data 
    tps5Density : numpy.ndarray
        tps5 density data
    tps8TipA : numpy.ndarray
        tps8 tip A voltage data.  (note that this channel is typically 
        disconnected)
    tps8TipB : numpy.ndarray
        tps8 tip B voltage data.  
    tps8TipC : numpy.ndarray
        tps8 tip C voltage data.  
    tps8Time : numpy.ndarray
        tps8 time data 
    tps8Current : numpy.ndarray
        tps8 current data
    tps8Temp : numpy.ndarray    
        tps8 temperature data.  
    tps8VFloat : numpy.ndarray
        tps8 floating voltage data 
    tps8Density : numpy.ndarray
        tps8 density data
        
    Subfunctions
    ------------
    plotOfISat : _plotTools.plot
        ion saturation current
    plotOfTipC : _plotTools.plot
        tip C voltages
    plotOfTipB : _plotTools.plot
        tip B voltages
    plotOfTipA : _plotTools.plot
        tip A voltages
    plotOfVf : _plotTools.plot
        floating voltages
    plotOfNe : _plotTools.plot
        density
    plotOfKTe : _plotTools.plot
        temperature
    plot :
        plots all relevant plots
        
    Notes
    -----
    I am not using the same time array for the section 5 or the section 8 
    triple probes.  I do this because different data acq. systems (e.g. north
    rack CPCI vs south rack CPCI) doesn't always return the EXACT same array.  
    I've run into issues where the length of the time arrays weren't the same
    length which causes problems during plotting. 
    
    TPS2 was moved to section 5 (now called TPS5) during the summer of 2017.
    This may cause variable naming issues.  Be warned.  
    
    TODO The other cases for shotnumbers need to finalized so that legacy data
    can still be loaded.  
    """
    
    def __init__(self,shotno=95996,tStart=_TSTART,tStop=_TSTOP,plot=False,probes='both'):  #sectionNum=2,
        
        self.shotno = shotno
        self.title = '%s, triple probe data' % shotno
        self.probes=probes
        
        # enforce probes naming convetion
        if probes=='5':
            probes = 'tps5'
        if probes=='8':
            probes = 'tps8'

        # constants
#        A=1.5904e-5 #(1.5mm)^2/4*pi + pi*(3.0mm)*(1.5mm), probe area
        A = 2.0*(1.5/1000.0*3.0/1000.0)
        e=1.602e-19; # fundamental charge
        eV=1.60218e-19; # 1 eV = 1.60218e-19 joules
        M=2.014102*1.66054e-27;  # mass of ion, approx 2 amu converted to kg
        me=9.109e-31; # mass of an electron
      
        ## Grab data
        if shotno > 95000: # Shotno after 2017 summer upgrade = 97239.  TPS2 was moved to section 5.  Now, it's TPS5.
            if probes=='both' or probes=='tps5' or probes=='tps2':
                
                # get data   
                try:         
                    data, time=mdsData(shotno=shotno,
                                  # TODO these addresses need to be updated to section 5 in the tree before they can be updated here
                                  dataAddress=['\HBTEP2::TOP.SENSORS.TRI_PROBE_S5.V_ION',
                                               '\HBTEP2::TOP.SENSORS.TRI_PROBE_S5.V_ELEC',
                                               '\HBTEP2::TOP.SENSORS.TRI_PROBE_S5.V_FLOAT',
                                               '\HBTEP2::TOP.SENSORS.TRI_PROBE_S5.I_SAT'],
                                  tStart=tStart, tStop=tStop)
                except:
                    data, time=mdsData(shotno=shotno,
                                  # TODO these addresses need to be updated to section 5 in the tree before they can be updated here
                                  dataAddress=['\HBTEP2::TOP.SENSORS.TRI_PROBE_S2.V_ION',
                                               '\HBTEP2::TOP.SENSORS.TRI_PROBE_S2.V_ELEC',
                                               '\HBTEP2::TOP.SENSORS.TRI_PROBE_S2.V_FLOAT',
                                               '\HBTEP2::TOP.SENSORS.TRI_PROBE_S2.I_SAT'],
                                  tStart=tStart, tStop=tStop)
                      
                # raw TPS5 Data
                self.tps5TipA = data[0] # the 180 is a ballparked number.  needs "actual" calibration       
                self.tps5TipB = data[1]
                self.tps5TipC = data[2]
                self.tps5Current=data[3]
                self.tps5Time = time
                
                # processed TPS5 Data
                self.tps5VFloat=self.tps5TipC;
                self.tps5Temp=(self.tps5TipB-self.tps5TipC)/.693;
                self.tps5Temp[self.tps5Temp>=200]=0; # trim data over 200eV.  I trim this data because there are a few VERY high temperature points that throw off the autoscaling
                tps5Temp=_copy(self.tps5Temp);
                tps5Temp[tps5Temp<=0]=1e6; # i trim here to avoid imaginary numbers when I take the square root below
                self.tps5Density=self.tps5Current/(0.605*e*_np.sqrt(tps5Temp*eV/(M))*A);
                self.tps5PlasmaPotential=self.tps5VFloat-self.tps5Temp/e*_np.log(0.6*_np.sqrt(2*_np.pi*me/M))
    
            if probes=='both' or probes=='tps8':
                
                # get data                 
                data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.TRI_PROBE_S8.V_ION',
                                           '\HBTEP2::TOP.SENSORS.TRI_PROBE_S8.V_ELEC',
                                           '\HBTEP2::TOP.SENSORS.TRI_PROBE_S8.V_FLOAT',
                                           '\HBTEP2::TOP.SENSORS.TRI_PROBE_S8.I_SAT'],
                              tStart=tStart, tStop=tStop)
                  
                # raw TPS8 Data
                self.tps8TipA = data[0] # the 180 is a ballparked number.  needs "actual" calibration       
                self.tps8TipB = data[1]
                self.tps8TipC = data[2]
                self.tps8Current=data[3]
                self.tps8Time = time
                
                # processed TPS8 Data
                self.tps8VFloat=self.tps8TipC;
                self.tps8Temp=(self.tps8TipB-self.tps8TipC)/.693;
                self.tps8Temp[self.tps8Temp>=200]=0; # trim data over 200eV.  I trim this data because there are a few VERY high temperature points that throw off the autoscaling
                tps8Temp=_copy(self.tps8Temp);
                tps8Temp[tps8Temp<=0]=1e6; # i trim here to avoid imaginary numbers when I take the square root below
                self.tps8Density=self.tps8Current/(0.605*e*_np.sqrt(tps8Temp*eV/(M))*A);
                self.tps8PlasmaPotential=self.tps8VFloat-self.tps8Temp/e*_np.log(0.6*_np.sqrt(2*_np.pi*me/M))
                
        else: # Shotno after 2017 summer upgrade = 97239.  TPS2 was moved to section 5.  Now, it's TPS5.
            if probes=='both' or probes=='tps5' or probes=='tps2':
                
                # get data                
                data, time=mdsData(shotno=shotno,
                              # TODO these addresses need to be updated to section 5 in the tree before they can be updated here
                              dataAddress=['\HBTEP2::TOP.SENSORS.TRI_PROBE_1.V_ION',
                                           '\HBTEP2::TOP.SENSORS.TRI_PROBE_1.V_ELEC',
                                           '\HBTEP2::TOP.SENSORS.TRI_PROBE_1.V_FLOAT',
                                           '\HBTEP2::TOP.SENSORS.TRI_PROBE_1.I_SAT'],
                              tStart=tStart, tStop=tStop)
                  
                # raw TPS5 Data
                self.tps5TipA = data[0] # the 180 is a ballparked number.  needs "actual" calibration       
                self.tps5TipB = data[1]
                self.tps5TipC = data[2]
                self.tps5Current=data[3]
                self.tps5Time = time
                
                # processed TPS5 Data
                self.tps5VFloat=self.tps5TipC;
                self.tps5Temp=(self.tps5TipB-self.tps5TipC)/.693;
                self.tps5Temp[self.tps5Temp>=200]=0; # trim data over 200eV.  I trim this data because there are a few VERY high temperature points that throw off the autoscaling
                tps5Temp=_copy(self.tps5Temp);
                tps5Temp[tps5Temp<=0]=1e6; # i trim here to avoid imaginary numbers when I take the square root below
                self.tps5Density=self.tps5Current/(0.605*e*_np.sqrt(tps5Temp*eV/(M))*A);
                self.tps5PlasmaPotential=self.tps5VFloat-self.tps5Temp/e*_np.log(0.6*_np.sqrt(2*_np.pi*me/M))

#                
#        else:
#            _sys.exit("Requested shot number range not supported yet.  Update code.")
        
        if plot==True:
            self.plot();
        elif plot=='all':
            self.plot(True)
        
    def plotOfKTe(self):
        p1=_plot.plot(yLabel='eV',subtitle='Electron Temperature',
                      title=self.title,yLim=[-50, 100],
                      shotno=self.shotno,xLabel='time [ms]');
        if self.probes=='both' or self.probes=='tps5': 
            p1.addTrace(xData=self.tps5Time*1000,yData=self.tps5Temp,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8Temp,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        return p1
        
    def plotOfNe(self):
        p1=_plot.plot(yLabel=r'$m^{-3}$ $10^{18}$',subtitle='Density',
                      yLim=[-1, 4.5],shotno=self.shotno,xLabel='time [ms]');
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5Density/1e18,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8Density/1e18,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        return p1
            
    def plotOfVf(self):
        p1=_plot.plot(yLabel='V',subtitle='Floating Potential',
                      xLabel='time [ms]',yLim=[-150, 75],shotno=[self.shotno]);
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5VFloat,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8VFloat,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        return p1
            
    def plotOfTipA(self):
        # initialize tip A potential plot
        p1=_plot.plot(yLabel='V',subtitle=r'Tip A, V$_{-}$',xLabel='time [ms]',
                      shotno=[self.shotno],title=self.title);
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5TipA,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8TipA,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        return p1
            
    def plotOfTipB(self):
        # initialize tip B potential plot
        p1=_plot.plot(yLabel='V',subtitle=r'Tip B, V$_{+}$',xLabel='time [ms]',
                      shotno=[self.shotno]);
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5TipB,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8TipB,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        return p1
            
    def plotOfTipC(self):            
        # initialize tip C potential plot
        p1=_plot.plot(yLabel='V',subtitle=r'Tip C, V$_{f}$',xLabel='time [ms]',
                      shotno=[self.shotno]);
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5TipC,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8TipC,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        return p1
        
    def plotOfISat(self):            
        # initialize ion sat current
        p1=_plot.plot(yLabel='A',xLabel='time [ms]',
                      subtitle='Ion Sat. Current',shotno=[self.shotno]);
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5Current,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8Current,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')   
        return p1
            
    def plot(self,plotAll=False):
        if plotAll == False:
            _plot.subPlot([self.plotOfKTe(),self.plotOfNe(),self.plotOfVf()]);
        else:
            _plot.subPlot([self.plotOfKTe(),self.plotOfNe(),self.plotOfVf()]);
#            _plot.subPlot([self.plotOfTipA(),self.plotOfTipB(),self.plotOfTipC(),
#                         self.plotOfISat()]);          
            _plot.subPlot([self.plotOfTipB(),self.plotOfTipC(),
                         self.plotOfISat()]);              
  
    
@_prepShotno
class paData:
    """
    Downloads poloidal array sensor data.  Presently, only poloidal
    measurements as the radial sensors are not yet implemeted.  
    
    Parameters
    ----------

    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        True - plots all 64 sensors
        'sample' - plots one of each (PA1 and PA2)
        'all' - same as True
    smoothingAlgorithm : str
        informs function as to which smoothing algorithm to use on each PA 
        sensor
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to put at the top of figures
    theta : numpy.ndarray
        poloidal location of sensors.  
    namesPA1 : numpy.ndarray
        1D array of all PA1 sensor names
    namesPA2 : numpy.ndarray
        1D array of all PA2 sensor names
    pa1Raw : list (of numpy.ndarray)
        raw PA1 sensor data
    pa2Raw : list (of numpy.ndarray)
        raw PA2 sensor data
    pa1Data : list (of numpy.ndarray)
        PA1 sensor data, processed
    pa2Data : list (of numpy.ndarray)
        PA2 sensor data, processed
    pa1RawFit : list (of numpy.ndarray)
        fit applied to raw data
    pa2RawFit : list (of numpy.ndarray)
        fit applied to raw data
        
    Subfunctions
    ------------
    plotOfPA1 : 
        returns plot of PA1 sensor based on the provided index
    plotOfPA2 : 
        returns plot of PA2 sensor based on the provided index
    plot :
        plots all relevant plots    
        
        
    Notes
    -----
    'PA2_S14P' is a known bad sensor
    pa1_s16 ???

    """
    
    def __init__(self,shotno=98170,tStart=_TSTART,tStop=_TSTOP,plot=False,
              removeBadSensors=True,correctTheta=False):
        self.shotno = shotno
        self.title1 = '%d, PA1 sensors' % shotno
        self.title2 = '%d, PA2 sensors' % shotno
        self.badSensors=['PA2_S14P','PA2_S27P']
        self.correctTheta=correctTheta
        self.tStart = tStart
        self.tStop = tStop
        
        # poloidal location
        self.thetaPA1 =_np.array([-174.74778518, -164.23392461, -153.66901098, -143.01895411,       -132.24974382, -121.3277924 , -110.22067715,  -98.93591492,        -87.23999699,  -75.60839722,  -63.97679673,  -52.34519359,        -40.71359604,  -29.08199717,  -17.45039318,   -5.81879416,          5.81280487,   17.44440438,   29.07600466,   40.70760263,         52.33920936,   63.97080017,   75.60240749,   87.23400093,         98.93591492,  110.22067715,  121.3277924 ,  132.24974382,        143.01895411,  153.66901098,  164.23392461,  174.74778518])*_np.pi/180.
        self.thetaPA2 =_np.array([-174.74778518, -164.23392461, -153.66901098, -143.01895411,       -132.24974382, -121.3277924 , -110.22067715,  -98.93591492,        -87.23999699,  -75.60839722,  -63.97679673,  -52.34519359,        -40.71359604,  -29.08199717,  -17.45039318,   -5.81879416,          5.81280487,   17.44440438,   29.07600466,   40.70760263,         52.33920936,   63.97080017,   75.60240749,   87.23400093,         98.93591492,  110.22067715,  121.3277924 ,  132.24974382,        143.01895411,  153.66901098,  164.23392461,  174.74778518])*_np.pi/180.
        
        # toroidal location
        self.phiPA1=_np.ones(len(self.thetaPA1))*317.5*_np.pi/180.
        self.phiPA2=_np.ones(len(self.thetaPA2))*137.5*_np.pi/180.
        
        # sensor names
        self.namesPA1=_np.array([   'PA1_S01P', 'PA1_S02P', 'PA1_S03P', 'PA1_S04P', 'PA1_S05P', 'PA1_S06P', 'PA1_S07P', 'PA1_S08P', 'PA1_S09P', 'PA1_S10P', 'PA1_S11P', 'PA1_S12P', 'PA1_S13P', 'PA1_S14P', 'PA1_S15P', 'PA1_S16P', 'PA1_S17P', 'PA1_S18P', 'PA1_S19P', 'PA1_S20P', 'PA1_S21P', 'PA1_S22P', 'PA1_S23P', 'PA1_S24P', 'PA1_S25P', 'PA1_S26P', 'PA1_S27P', 'PA1_S28P', 'PA1_S29P', 'PA1_S30P', 'PA1_S31P', 'PA1_S32P'])
        self.namesPA2=_np.array([   'PA2_S01P', 'PA2_S02P', 'PA2_S03P', 'PA2_S04P', 'PA2_S05P', 'PA2_S06P', 'PA2_S07P', 'PA2_S08P', 'PA2_S09P', 'PA2_S10P', 'PA2_S11P', 'PA2_S12P', 'PA2_S13P', 'PA2_S14P', 'PA2_S15P', 'PA2_S16P', 'PA2_S17P', 'PA2_S18P', 'PA2_S19P', 'PA2_S20P', 'PA2_S21P', 'PA2_S22P', 'PA2_S23P', 'PA2_S24P', 'PA2_S25P', 'PA2_S26P', 'PA2_S27P', 'PA2_S28P', 'PA2_S29P', 'PA2_S30P', 'PA2_S31P', 'PA2_S32P'])

                    

        if self.correctTheta:
            self.thetaPA1=processPlasma.thetaCorrection(self.shotno,self.thetaPA1,\
                    self.tStart,self.tStop)[1]
            self.thetaPA2=processPlasma.thetaCorrection(self.shotno,self.thetaPA2,\
                    self.tStart,self.tStop)[1]
        # compile full sensor addresses names
        pa1SensorAddresses=[]
        pa2SensorAddresses=[]        
        rootAddress='\HBTEP2::TOP.SENSORS.MAGNETIC:';
        for i in range(0,len(self.namesPA1)):
            pa1SensorAddresses.append(rootAddress+self.namesPA1[i])
        for i in range(0,len(self.namesPA2)):
            pa2SensorAddresses.append(rootAddress+self.namesPA2[i])
            
        # get raw data
        self.pa1Raw,self.pa1Time=mdsData(shotno,pa1SensorAddresses, tStart, tStop)
        self.pa2Raw,self.pa2Time=mdsData(shotno,pa2SensorAddresses, tStart, tStop)
        
        # data smoothing algorithm
        self.pa1Data=[]
        self.pa1RawFit=[]
        self.pa2Data=[]
        self.pa2RawFit=[]
        
        # gaussian offset subtraction
        for i in range(0,len(self.namesPA1)):
            temp,temp2=_process.gaussianHighPassFilter(self.pa1Raw[i][:],self.pa1Time,timeWidth=1./20000)
            self.pa1RawFit.append(temp2)
            self.pa1Data.append(temp)
        for i in range(0,len(self.namesPA2)):
            temp,temp2=_process.gaussianHighPassFilter(self.pa2Raw[i][:],self.pa2Time,timeWidth=1./20000)
            self.pa2RawFit.append(temp2)
            self.pa2Data.append(temp)
            
        # pandas dataframes 
        #TODO(John) rewrite entire class.  start with dataframes instead of lists
        self.dfData=_pd.DataFrame(     data=_np.append(_np.array([self.pa1Time]).transpose(),_np.array(self.pa1Data+self.pa2Data).transpose(),axis=1),
                            columns=_np.append(_np.append(_np.array(['Time']),self.namesPA1),self.namesPA2)).set_index('Time')
        self.dfDataRaw=_pd.DataFrame(     data=_np.append(_np.array([self.pa1Time]).transpose(),_np.array(self.pa1Raw+self.pa2Raw).transpose(),axis=1),
                            columns=_np.append(_np.append(_np.array(['Time']),self.namesPA1),self.namesPA2)).set_index('Time')
        self.dfMeta=_pd.DataFrame(     data={'SensorNames':_np.append(self.namesPA1,self.namesPA2),
                                    'Phi':_np.append(self.phiPA1,self.phiPA2),
                                    'Theta':_np.append(self.thetaPA1,self.thetaPA2),
                                    'Address':pa1SensorAddresses+pa2SensorAddresses,
                                    'SectionNum':_np.append([3]*32,[8]*32),
                                    },
                                columns=['SensorNames','Phi','Theta','Address','SectionNum']).set_index('SensorNames')#,'Address','SectionNum'])
        
        
        if removeBadSensors==True:
            self.dfData=self.dfData.drop(columns=self.badSensors)
            self.dfMeta=self.dfMeta.drop(index=self.badSensors)
            
            # Non-Pandas solution
            
            for i in self.badSensors:
                iBad=_np.where(self.namesPA1==i)
                
                self.namesPA1=_np.delete(self.namesPA1,iBad)
                self.thetaPA1=_np.delete(self.thetaPA1,iBad)
                self.phiPA2=_np.delete(self.phiPA2,iBad)
                
                iBad=_np.where(self.namesPA2==i)
                self.namesPA2=_np.delete(self.namesPA2,iBad)
                self.thetaPA2=_np.delete(self.thetaPA2,iBad)
                self.phiPA2=_np.delete(self.phiPA2,iBad)
                '''
                self.pa1Raw=_np.delete(self.pa1Raw,iBad)
                self.pa1RawFit=_np.delete(self.pa1RawFit,iBad)
                self.pa1Data=_np.delete(self.pa1Data,iBad) 
                
                
                
                self.pa2Raw=_np.delete(self.pa2Raw,iBad)
                self.pa2RawFit=_np.delete(self.pa2RawFit,iBad)
                self.pa2Data=_np.delete(self.pa2Data,iBad)
                '''
                
        if plot==True or plot=='all':
            self.plot(True)
        if plot=='sample':
            self.plotOfPA1().plot();
            self.plotOfPA2().plot();
            
    def plotOfPA1Stripey(self,tStart=2e-3,tStop=4e-3):
        iStart=_process.findNearest(self.pa1Time,tStart)
        iStop=_process.findNearest(self.pa1Time,tStop)
        p1=_plot.plot(title=self.title1,subtitle='PA1 Sensors',
                      xLabel='Time [ms]', yLabel='theta [rad]',zLabel='Gauss',
                      plotType='contour',colorMap=_plot._red_green_colormap(),
                      centerColorMapAroundZero=True)
        data=self.pa1Data[0:len(self.namesPA1)]
        for i in range(0,len(data)):
            data[i]=data[i][iStart:iStop]*1e4
        p1.addTrace(self.pa1Time[iStart:iStop]*1e3,self.thetaPA1,
                    _np.array(data))
        return p1
        
    def plotOfPA2Stripey(self,tStart=2e-3,tStop=4e-3):
        iStart=_process.findNearest(self.pa2Time,tStart)
        iStop=_process.findNearest(self.pa2Time,tStop)
        p1=_plot.plot(title=self.title2,subtitle='PA2 Sensors',
                      xLabel='Time [ms]', yLabel='theta [rad]',zLabel='Gauss',
                      plotType='contour',colorMap=_plot._red_green_colormap(),
                      centerColorMapAroundZero=True)
        data=self.pa2Data[0:len(self.namesPA2)]
        data2=[]
        theta=[]
        for i in range(0,len(data)):
            if self.namesPA2[i] in ['PA2_S14P', 'PA2_S27P']:
                pass
            else:
                data2.append(data[i][iStart:iStop]*1e4)
                theta.append(self.thetaPA2[i])
        p1.addTrace(self.pa2Time[iStart:iStop]*1e3,theta,
                    data2)
        return p1

    def plotOfPA1(self, i=0, alsoPlotRawAndFit=True):
        """ Plot one of the PA1 plots.  based on index, i. """
        p1=_plot.plot(xLabel='time [ms]',yLabel=r'Gauss',title=self.title1,
                      shotno=[self.shotno],subtitle=self.namesPA1[i]);
        
        # smoothed data
        p1.addTrace(yData=self.pa1Data[i],xData=self.pa1Time*1000,
                    yLegendLabel='smoothed')   

        if alsoPlotRawAndFit==True:
            # raw data
            p1.addTrace(yData=self.pa1Raw[i],xData=self.pa1Time*1000,
                        yLegendLabel='raw')   
            
            # fit data (which is subtracted from raw)
            p1.addTrace(yData=self.pa1RawFit[i],xData=self.pa1Time*1000,
                        yLegendLabel='fit')   
            
        return p1
        
    def plotOfPA2(self, i=0, alsoPlotRawAndFit=True):
        """ Plot one of the PA2 plots.  based on index, i. """
        p1=_plot.plot(xLabel='time [ms]',yLabel='Gauss',
                      title=self.title2,subtitle=self.namesPA2[i],
                      shotno=self.shotno);
        
        # smoothed data
        p1.addTrace(yData=self.pa2Data[i],xData=self.pa2Time*1000,
                    yLegendLabel='smoothed')   

        if alsoPlotRawAndFit==True:
            # raw data
            p1.addTrace(yData=self.pa2Raw[i],xData=self.pa2Time*1000,
                        yLegendLabel='raw')   
            
            # fit data (which is subtracted from raw)
            p1.addTrace(yData=self.pa2RawFit[i],xData=self.pa2Time*1000,
                        yLegendLabel='fit')   

        return p1
        
    def plot(self,plotAll=False):
        sp1=[[],[],[],[]]
        sp2=[[],[],[],[]]
        count=0
        for i in range(0,4):
            for j in range(0,8):
                if plotAll==True:
                    newPlot=self.plotOfPA1(count,alsoPlotRawAndFit=True)
                else:
                    newPlot=self.plotOfPA1(count,alsoPlotRawAndFit=False)
                newPlot.subtitle=self.namesPA1[count]
                newPlot.yLegendLabel=[]
                sp1[i].append(newPlot)
                count+=1;

        k=0
        count=0
        for i in range(0,4):
            for j in range(0,8):
                k=i*8+j*1
                print("i %d, j %d, k %d"%(i,j,k))
                
                # check to see if all 32 sensors are present
                if k>=len(self.namesPA2):
                    # and create an empty plot if not
                    newPlot=_plot.plot()
                else:
                    print("%d" %k)
                    if plotAll==True:
                        newPlot=self.plotOfPA2(count,alsoPlotRawAndFit=True)
                    else:
                        newPlot=self.plotOfPA2(count,alsoPlotRawAndFit=False)
                    newPlot.subtitle=self.namesPA2[count]
                    newPlot.yLegendLabel=[]
                    sp2[i].append(newPlot)
                    count+=1;
        sp1[0][0].title=self.title1
        sp2[0][0].title=self.title2
        sp1=_plot.subPlot(sp1,plot=False)
        sp2=_plot.subPlot(sp2,plot=False)
        # sp1.shareY=True;
        sp1.plot()
        sp2.plot()
    
    
    # Build polar plot a la JEff
    def plotPolar(self, tPoint = 1e-3):
        tPoint = _process.findNearest(self.pa1Time,tPoint) # time index
        data = _np.array(self.pa1Data)[:,tPoint]
        data = _np.hstack( (data,data[0])) # wrap data
        offset = _np.min(data)
        if self.correctTheta:
            self.thetaPA1=processPlasma.thetaCorrection(self.shotno,self.thetaPA1,\
                    self.tStart,self.tStop)[0]
            self.thetaPA1=self.thetaPA1[tPoint,:]
        # Build plot 
        _plt.figure()
        ax = _plt.subplot(111,projection='polar')
        ax.plot(_np.hstack((self.thetaPA1,self.thetaPA1[0])),\
                data+1.5*offset,'-*')
        _plt.show()
        
        
@_prepShotno
def paData_df(shotno=98170,tStart=_TSTART,tStop=_TSTOP,plot=False,
              removeBadSensors=True,correctTheta=False):
    """
    Downloads poloidal array sensor data.  Presently, only poloidal
    measurements as the radial sensors are not yet implemeted.  
    
    Parameters
    ----------

    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        True - plots all 64 sensors
        'sample' - plots one of each (PA1 and PA2)
        'all' - same as True
    smoothingAlgorithm : str
        informs function as to which smoothing algorithm to use on each PA 
        sensor
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to put at the top of figures
    theta : numpy.ndarray
        poloidal location of sensors.  
    namesPA1 : numpy.ndarray
        1D array of all PA1 sensor names
    namesPA2 : numpy.ndarray
        1D array of all PA2 sensor names
    pa1Raw : list (of numpy.ndarray)
        raw PA1 sensor data
    pa2Raw : list (of numpy.ndarray)
        raw PA2 sensor data
    pa1Data : list (of numpy.ndarray)
        PA1 sensor data, processed
    pa2Data : list (of numpy.ndarray)
        PA2 sensor data, processed
    pa1RawFit : list (of numpy.ndarray)
        fit applied to raw data
    pa2RawFit : list (of numpy.ndarray)
        fit applied to raw data
        
    Subfunctions
    ------------
    plotOfPA1 : 
        returns plot of PA1 sensor based on the provided index
    plotOfPA2 : 
        returns plot of PA2 sensor based on the provided index
    plot :
        plots all relevant plots    
        
        
    Notes
    -----
    'PA2_S14P' is a known bad sensor
    pa1_s16 ???

    """

    shotno = shotno
    title1 = '%d, PA1 sensors' % shotno
    title2 = '%d, PA2 sensors' % shotno
    badSensors=['PA2_S14P','PA2_S27P']
    correctTheta=correctTheta
    tStart = tStart
    tStop = tStop
    
    # poloidal location
    thetaPA1 =_np.array([-174.74778518, -164.23392461, -153.66901098, -143.01895411,       -132.24974382, -121.3277924 , -110.22067715,  -98.93591492,        -87.23999699,  -75.60839722,  -63.97679673,  -52.34519359,        -40.71359604,  -29.08199717,  -17.45039318,   -5.81879416,          5.81280487,   17.44440438,   29.07600466,   40.70760263,         52.33920936,   63.97080017,   75.60240749,   87.23400093,         98.93591492,  110.22067715,  121.3277924 ,  132.24974382,        143.01895411,  153.66901098,  164.23392461,  174.74778518])*_np.pi/180.
    thetaPA2 =_np.array([-174.74778518, -164.23392461, -153.66901098, -143.01895411,       -132.24974382, -121.3277924 , -110.22067715,  -98.93591492,        -87.23999699,  -75.60839722,  -63.97679673,  -52.34519359,        -40.71359604,  -29.08199717,  -17.45039318,   -5.81879416,          5.81280487,   17.44440438,   29.07600466,   40.70760263,         52.33920936,   63.97080017,   75.60240749,   87.23400093,         98.93591492,  110.22067715,  121.3277924 ,  132.24974382,        143.01895411,  153.66901098,  164.23392461,  174.74778518])*_np.pi/180.
    
    # toroidal location
    phiPA1=_np.ones(len(thetaPA1))*317.5*_np.pi/180.
    phiPA2=_np.ones(len(thetaPA2))*137.5*_np.pi/180.
    
    # sensor names
    namesPA1=_np.array([   'PA1_S01P', 'PA1_S02P', 'PA1_S03P', 'PA1_S04P', 'PA1_S05P', 'PA1_S06P', 'PA1_S07P', 'PA1_S08P', 'PA1_S09P', 'PA1_S10P', 'PA1_S11P', 'PA1_S12P', 'PA1_S13P', 'PA1_S14P', 'PA1_S15P', 'PA1_S16P', 'PA1_S17P', 'PA1_S18P', 'PA1_S19P', 'PA1_S20P', 'PA1_S21P', 'PA1_S22P', 'PA1_S23P', 'PA1_S24P', 'PA1_S25P', 'PA1_S26P', 'PA1_S27P', 'PA1_S28P', 'PA1_S29P', 'PA1_S30P', 'PA1_S31P', 'PA1_S32P'])
    namesPA2=_np.array([   'PA2_S01P', 'PA2_S02P', 'PA2_S03P', 'PA2_S04P', 'PA2_S05P', 'PA2_S06P', 'PA2_S07P', 'PA2_S08P', 'PA2_S09P', 'PA2_S10P', 'PA2_S11P', 'PA2_S12P', 'PA2_S13P', 'PA2_S14P', 'PA2_S15P', 'PA2_S16P', 'PA2_S17P', 'PA2_S18P', 'PA2_S19P', 'PA2_S20P', 'PA2_S21P', 'PA2_S22P', 'PA2_S23P', 'PA2_S24P', 'PA2_S25P', 'PA2_S26P', 'PA2_S27P', 'PA2_S28P', 'PA2_S29P', 'PA2_S30P', 'PA2_S31P', 'PA2_S32P'])

                

#    if correctTheta:
#        thetaPA1=processPlasma.thetaCorrection(shotno,thetaPA1,\
#                tStart,tStop)[1]
#        thetaPA2=processPlasma.thetaCorrection(shotno,thetaPA2,\
#                tStart,tStop)[1]
    # compile full sensor addresses names
    pa1SensorAddresses=[]
    pa2SensorAddresses=[]        
    rootAddress='\HBTEP2::TOP.SENSORS.MAGNETIC:';
    for i in range(0,len(namesPA1)):
        pa1SensorAddresses.append(rootAddress+namesPA1[i])
    for i in range(0,len(namesPA2)):
        pa2SensorAddresses.append(rootAddress+namesPA2[i])
        
    # get raw data
    pa1Raw,pa1Time=mdsData(shotno,pa1SensorAddresses, tStart, tStop)
    pa2Raw,pa2Time=mdsData(shotno,pa2SensorAddresses, tStart, tStop)
    
    # data smoothing algorithm
    pa1Data=[]
    pa1RawFit=[]
    pa2Data=[]
    pa2RawFit=[]
    
    # gaussian offset subtraction
    for i in range(0,len(namesPA1)):
        temp,temp2=_process.gaussianHighPassFilter(pa1Raw[i][:],pa1Time,timeWidth=1./20000)
        pa1RawFit.append(temp2)
        pa1Data.append(temp)
    for i in range(0,len(namesPA2)):
        temp,temp2=_process.gaussianHighPassFilter(pa2Raw[i][:],pa2Time,timeWidth=1./20000)
        pa2RawFit.append(temp2)
        pa2Data.append(temp)
        
    # pandas dataframes 
    dfData=_pd.DataFrame(     data=_np.append(_np.array([pa1Time]).transpose(),_np.array(pa1Data+pa2Data).transpose(),axis=1),
                        columns=_np.append(_np.append(_np.array(['Time']),namesPA1),namesPA2)).set_index('Time')
    dfDataRaw=_pd.DataFrame(     data=_np.append(_np.array([pa1Time]).transpose(),_np.array(pa1Raw+pa2Raw).transpose(),axis=1),
                        columns=_np.append(_np.append(_np.array(['Time']),namesPA1),namesPA2)).set_index('Time')
    dfMeta=_pd.DataFrame(     data={'SensorNames':_np.append(namesPA1,namesPA2),
                                'Phi':_np.append(phiPA1,phiPA2),
                                'Theta':_np.append(thetaPA1,thetaPA2),
                                'Address':pa1SensorAddresses+pa2SensorAddresses,
                                'SectionNum':_np.append([3]*32,[8]*32),
                                },
                            columns=['SensorNames','Phi','Theta','Address','SectionNum']).set_index('SensorNames')#,'Address','SectionNum'])
    
    if removeBadSensors==True:
        dfData=dfData.drop(columns=badSensors)
        dfDataRaw=dfDataRaw.drop(columns=badSensors)
        dfMeta=dfMeta.drop(index=badSensors)
        
#    return dfData,dfDataRaw,dfMeta
    return dfData,dfDataRaw,dfMeta
        
        
        
@_prepShotno
class sxrData:
    """
    Downloads (and optionally plots) soft xray sensor data.   
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        default is False
        True - plots far array of all 11 (of 16) channels
        'all' - plots 
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to put at the top of figures
    data : list (of numpy.ndarray)
        list of 11 (of 16) data arrays, one for each channel
        
    Subfunctions
    ------------
    plotAll :
        plots all relevant plots     
    plotOfSXRStripey : 
        returns a stripey plot of the sensors
    plotOfOneChannel :
        returns a plot of a single channel based on the provided index, i
        
    Notes
    -----
    Only 11 (of 16) of the SXR sensors are included in the data below.  Some of the 
    missing sensors are broken and others include anamolous or attenuated 
    results.  
    

    """
    
    def __init__(self,shotno=98170,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno = shotno
        self.title = '%d, SXR sensors' % shotno
        
        # note that sensors 5, 8, 10, 13 and 15 are not included.  
        self.sensor_num=_np.array([  1,  2,  3,  4,  6,   9, 11, 12, 14, 16])
#        self.sensor_num=_np.array([ 0,  1,  2,  3,  4,  5,  6,  7, 8,  9,  10, 11, 12, 13, 14, 15])
        channels=self.sensor_num+75;

        # compile full sensor addresses names
        addresses=[]
        self.sensorNames=[]
        mdsAddressRaw = []
        mdsAddressR = []
        mdsAddressZ = []
        mdsAddressMid = []
        for i in self.sensor_num:
            self.sensorNames.append('CHANNEL_%02d'%i)
            #sensorAddresses.append('\HBTEP2::TOP.DEVICES.WEST_RACK:CPCI:INPUT_%02d' %channels[i])
            address = '\HBTEP2::TOP.SENSORS.SXR_FAN.CHANNEL_%02d' %i                
            mdsAddressRaw.append(address + ':RAW')
            mdsAddressR.append(address+':R')
            mdsAddressZ.append(address+':Z')
            mdsAddressMid.append(address+':MIDPLANE')
        # get data
        #self.data,self.time=mdsData(shotno,sensorAddresses, tStart, tStop)
        self.data, self.time = mdsData(shotno,mdsAddressRaw,tStart,tStop)
        self.R = mdsData(shotno,mdsAddressR)
        self.Z = mdsData(shotno,mdsAddressZ)
        self.Midplane = mdsData(shotno,mdsAddressMid)
        # aperture location
        self.det_ap_R = [.9654660224914551]
        self.det_ap_Z = [.2211250066757202]
        self.det_ap_Pol = [.001270000007934868]
        self.det_ap_Tor =[.02539999969303608]
        # plot 
        if plot==True:
            self.plotOfSXRStripey(tStart,tStop).plot()
        elif plot=='all':
            self.plotAll()
            self.plotOfSXRStripey(tStart,tStop).plot()
            
            
    def plotOfSXRStripey(self,tStart=1e-3,tStop=10e-3):
        iStart=_process.findNearest(self.time,tStart)
        iStop=_process.findNearest(self.time,tStop)
        p1=_plot.plot(title=self.title,subtitle='SXR Fan Array',
                      xLabel='Time [ms]', yLabel='Sensor Number',zLabel='a.u.',
                      plotType='contour',#colorMap=_plot._red_green_colormap(),
                      centerColorMapAroundZero=True)
        data=self.data;
        for i in range(0,len(data)):
            data[i]=data[i][iStart:iStop]
        p1.addTrace(self.time[iStart:iStop]*1e3,self.sensor_num,
                    _np.array(data))
        return p1
        
            
    def plotOfOneChannel(self, i=0):
        """ Plot one of the SXR chanenl.  based on index, i. """
        p1=_plot.plot(xLabel='time [ms]',yLabel=r'a.u.',title=self.title,
                      shotno=[self.shotno],subtitle=self.sensorNames[i]);
        
        # smoothed data
        p1.addTrace(yData=self.data[i],xData=self.time*1000,
                    yLegendLabel=self.sensorNames[i])   
            
        return p1
        

    def plotAll(self):
        sp1=[]
        count=0
        for i in range(0,len(self.data)):
                newPlot=self.plotOfOneChannel(count)
                newPlot.subtitle=self.sensorNames[count]
                newPlot.yLegendLabel=[]
                newPlot.plot()
                sp1.append(newPlot)
 
                count+=1;
                
        sp1[0].title=self.title
        sp1=_plot.subPlot(sp1,plot=False)
        # sp1.shareY=True;
#        sp1.plot()
        
        return sp1
            
    
def fbData_df(shotno=98170,tStart=_TSTART,tStop=_TSTOP,plot=False,badSensors=['FB03_S1P','FB06_S2P','FB08_S3P'],poloidalOnly=True):
    
    rootAddress='\HBTEP2::TOP.SENSORS.MAGNETIC:';
    
    names=[]
    addresses=[]
    
    if poloidalOnly==True:
        sensors=["P"]
    else:
        sensors=["P","R"]
    for k in sensors:
        for i in range(4):
#            for j in range(10):
            for j in _np.array([5,6,7,8,9,10,1,2,3,4])-1:
                names.append('FB%1.2d_S%d%s'%(j+1,i+1,k))
                addresses.append('%s%s'%(rootAddress,names[-1]))
    
    dataRaw, time =mdsData(shotno,addresses, tStart, tStop)
#     dataRaw
    dfData=_pd.DataFrame(data=_np.array(dataRaw).transpose(),columns=["%s_RAW"%i for i in names])
    dfData['time']=time
    dfData=dfData.set_index('time')
    
#    phi=_np.pi/180.*_np.array([241,277,313,349,25,61, 97,133,169,205]*4)
    phi=_np.pi/180.*_np.array([25,61, 97,133,169,205,241,277,313,349]*4)
    theta=_np.pi/180.*_np.array([_np.ones(10)*(-83.4),_np.ones(10)*(-29.3),_np.ones(10)*29.3,_np.ones(10)*83.4]).flatten()
    dfMeta=_pd.DataFrame(index=names)
    dfMeta['Phi']=phi
    dfMeta['Theta']=theta
    dfMeta['address']=addresses
    
    dfMeta=dfMeta.drop(index=badSensors)
    badSensors=['%s_RAW'%i for i in badSensors]
    dfData=dfData.drop(columns=badSensors)
    
    for i, (key, vals) in enumerate(dfData.iteritems()):
        t=dfData.index.to_numpy()
        dfData[key[0:8]]=_process.gaussianFilter(t,dfData[key],timeFWHM=5e-4,filterType='high',plot=False)

    return dfData,dfMeta


        
@_prepShotno
class fbData:
    """
    Downloads feedback (FB) array sensor data.   
    
    Parameters
    ----------

    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        True - Plots a sample of each FB poloidal and radial data
        'sample'- same as True
        'all' - Plots all 80 sensor data
    smoothingAlgorithm : str
        informs function as to which smoothing algorithm to use on each PA 
        sensor
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to be added to each plot
    fbPolNames : 2D list (of str)
        name of every poloidal FB sensor
    fbRadNames : 2D list (of str)
        name of every radial FB sensor
    phi : 2D numpy.ndarray
        toroidal locations for all sensors.  units in radians.
    theta : 2D numpy.ndarray
        poloidal locations for all sensors.  units in radians.
    fbPolRaw : 2D list (of numpy.ndarray)
        raw FB-poloidal data
    fbRadRaw : 2D list (of numpy.ndarray)
        raw FB-radial data
    fbPolData : 2D list (of numpy.ndarray)
        FB-poloidal data, processed
    fbRadData : 2D list (of numpy.ndarray)
        FB-radial data, processed
    fbPolRawFit : 2D list (of numpy.ndarray)
        smoothed fit of raw poloidal data.  subtracted from data to get 
        fbPolData
    fbRadRawFit : 2D list (of numpy.ndarray)
        smoothed fit of raw radial data.  subtracted from data to get fbRadData
        
    Subfunctions
    ------------
    plotOfSinglePol :
        returns plot of a specified poloidal sensor
    plotOfSingleRad :
        returns plot of a specified radial sensor
    plot :
        plots all relevant data
    
    Notes
    -----
    Known bad sensors: 'FB03_S1P','FB06_S2P','FB08_S3P'
    The S4P array has no broken sensors at present.  
    """
    def __init__(self,shotno=98170,tStart=_TSTART,tStop=_TSTOP,plot=False,removeBadSensors=True,invertNegSignals=True):
        self.shotno = shotno
        self.title = "%d, FB sensors" % shotno
#        self.badSensors=['FB03_S1P','FB06_S2P','FB08_S3P'] # some sensors appear to be broken

        # sensor names
        fbPolNames=[['FB01_S1P', 'FB02_S1P', 'FB03_S1P', 'FB04_S1P', 'FB05_S1P', 'FB06_S1P', 'FB07_S1P', 'FB08_S1P', 'FB09_S1P', 'FB10_S1P'], ['FB01_S2P', 'FB02_S2P', 'FB03_S2P', 'FB04_S2P', 'FB05_S2P', 'FB06_S2P', 'FB07_S2P', 'FB08_S2P', 'FB09_S2P', 'FB10_S2P'], ['FB01_S3P', 'FB02_S3P', 'FB03_S3P', 'FB04_S3P', 'FB05_S3P', 'FB06_S3P', 'FB07_S3P', 'FB08_S3P', 'FB09_S3P', 'FB10_S3P'], ['FB01_S4P', 'FB02_S4P', 'FB03_S4P', 'FB04_S4P', 'FB05_S4P', 'FB06_S4P', 'FB07_S4P', 'FB08_S4P', 'FB09_S4P', 'FB10_S4P']]
        self.fbRadNames=[['FB01_S1R', 'FB02_S1R', 'FB03_S1R', 'FB04_S1R', 'FB05_S1R', 'FB06_S1R', 'FB07_S1R', 'FB08_S1R', 'FB09_S1R', 'FB10_S1R'], ['FB01_S2R', 'FB02_S2R', 'FB03_S2R', 'FB04_S2R', 'FB05_S2R', 'FB06_S2R', 'FB07_S2R', 'FB08_S2R', 'FB09_S2R', 'FB10_S2R'], ['FB01_S3R', 'FB02_S3R', 'FB03_S3R', 'FB04_S3R', 'FB05_S3R', 'FB06_S3R', 'FB07_S3R', 'FB08_S3R', 'FB09_S3R', 'FB10_S3R'], ['FB01_S4R', 'FB02_S4R', 'FB03_S4R', 'FB04_S4R', 'FB05_S4R', 'FB06_S4R', 'FB07_S4R', 'FB08_S4R', 'FB09_S4R', 'FB10_S4R']]
        
        # sensor, toroidal location
#        self.phi=_np.pi/180.*_np.array([242.5-360, 278.5-360, 314.5-360, 350.5-360, 26.5, 62.5, 98.5, 134.5, 170.5, 206.5]);#*_np.pi/180.
        phi=_np.pi/180.*_np.array([241,277,313,349,25,61, 97,133,169,205])
        phi=[phi,phi,phi,phi]
        theta=_np.pi/180.*_np.array([_np.ones(10)*(-83.4),_np.ones(10)*(-29.3),_np.ones(10)*29.3,_np.ones(10)*83.4])
        theta=[theta[0,:],theta[1,:],theta[2,:],theta[3,:]]
        
        ## construct full sensor addresses 
        fbPolSensorAddresses=[[],[],[],[]]
        fbRadSensorAddresses=[[],[],[],[]]   
        rootAddress='\HBTEP2::TOP.SENSORS.MAGNETIC:';
        for j in range(0,4):
            for i in range(0,len(fbPolNames[j])):
                fbPolSensorAddresses[j].append(rootAddress+fbPolNames[j][i])
                fbRadSensorAddresses[j].append(rootAddress+self.fbRadNames[j][i])
        
        # get raw data
        fbPolRaw=[[],[],[],[]];
        fbPolRaw[0], self.fbPolTime =mdsData(shotno,fbPolSensorAddresses[0], tStart, tStop)
        fbPolRaw[1], self.fbPolTime =mdsData(shotno,fbPolSensorAddresses[1], tStart, tStop)
        fbPolRaw[2], self.fbPolTime =mdsData(shotno,fbPolSensorAddresses[2], tStart, tStop)
        fbPolRaw[3], self.fbPolTime =mdsData(shotno,fbPolSensorAddresses[3], tStart, tStop)

        self.fbRadRaw=[[],[],[],[]];
        self.fbRadRaw[0], self.fbRadTime =mdsData(shotno,fbRadSensorAddresses[0], tStart, tStop)
        self.fbRadRaw[1], self.fbRadTime =mdsData(shotno,fbRadSensorAddresses[1], tStart, tStop)
        self.fbRadRaw[2], self.fbRadTime =mdsData(shotno,fbRadSensorAddresses[2], tStart, tStop)
        self.fbRadRaw[3], self.fbRadTime =mdsData(shotno,fbRadSensorAddresses[3], tStart, tStop)        
               
        # remove bad/broken sensors using a sigma=1 outlier rejection method
        if removeBadSensors==True:
            linArrayOfSensorData=[]
            linListOfSensorNames=[]
            for i in range(0,4):
                for j in range(0,len(fbPolNames[i])):
                    linArrayOfSensorData.append(_np.average(_np.abs(fbPolRaw[i][j])))
                    linListOfSensorNames.append(fbPolNames[i][j])
            temp,indicesOfGoodSensors=_process.rejectOutliers(_np.array(linArrayOfSensorData),sigma=1.0)
            indicesOfBadSensors = indicesOfGoodSensors==False
            self.badSensors=_np.array(linListOfSensorNames)[indicesOfBadSensors]
            self.fbPolNames=[[],[],[],[]]
            self.fbPolRaw=[[],[],[],[]];
            self.phi=[[],[],[],[]];
            self.theta=[[],[],[],[]];
            for j in range(0,4):
                for i in range(0,len(fbPolNames[j])):
                    if fbPolNames[j][i] not in self.badSensors:
                        self.fbPolRaw[j].append(fbPolRaw[j][i])
                        self.fbPolNames[j].append(fbPolNames[j][i])
                        self.theta[j].append(theta[j][i])
                        self.phi[j].append(phi[j][i])
                    else:
                        print("Removing broken signal: %s" % fbPolNames[j][i])
        else:
            self.fbPolNames=fbPolNames
            self.fbPolRaw=fbPolRaw
            self.theta=theta
            self.phi=phi

        # make sure the signals are not inverted
        if invertNegSignals==True:
            for i in range(0,4):
                for j in range(0,len(self.fbPolRaw[i])):
                    if _np.average(self.fbPolRaw[i][j])<0:
                        self.fbPolRaw[i][j]*=-1
                        print("inverting signal %s"%self.fbPolNames[i][j])
        
        # remove low-frequency offset (we are only interested in high-freq data)
        self.fbPolData=[[],[],[],[]]
        self.fbPolRawFit=[[],[],[],[]]
        self.fbRadData=[[],[],[],[]]
        self.fbRadRawFit=[[],[],[],[]]
        for j in range(0,4):
            for i in range(0,len(self.fbPolNames[j])):
                temp,temp2=_process.gaussianHighPassFilter(self.fbPolRaw[j][i][:],self.fbPolTime,timeWidth=1./20000*1.0,plot=False) 
                self.fbPolRawFit[j].append(temp2)
                self.fbPolData[j].append(temp)
                
        # pandas dataframes 
        #TODO(John) rewrite entire class.  start with dataframes instead of lists
        flatten = lambda l: [item for sublist in l for item in sublist]
        self.dfData=_pd.DataFrame(     data=_np.append(_np.array([self.fbPolTime]).transpose(),_np.array(flatten(self.fbPolData)).transpose(),axis=1),
                                columns=['Time']+flatten(self.fbPolNames)).set_index('Time')
        self.dfDataRaw=_pd.DataFrame(     data=_np.append(_np.array([self.fbPolTime]).transpose(),_np.array(flatten(self.fbPolRaw)).transpose(),axis=1),
                                columns=['Time']+flatten(self.fbPolNames)).set_index('Time')
        self.dfMeta=_pd.DataFrame(     data={'SensorNames':flatten(self.fbPolNames),
                                    'Phi':flatten(self.phi),
                                    'Theta':flatten(self.theta),
#                                    'Address':dataAddress,
#                                    'SectionNum':[ ],
#                                    'RowNum':[1,1,1....2,2,2....,4,4,4,4],
                                    },
                                columns=['SensorNames','Phi','Theta']).set_index('SensorNames')#,'Address','SectionNum'])
        

        # plot
        if plot=='sample':
            self.plotOfSinglePol().plot();
            self.plotOfSingleRad().plot();
        elif plot == True or plot=='all':
            self.plot(True)
            
    def plotOfFBPolStripey(self,tStart=2e-3,tStop=4e-3,sensorArray='S4P'):
        # grab and trim data to desired time rane
        iStart=_process.findNearest(self.fbPolTime,tStart)
        iStop=_process.findNearest(self.fbPolTime,tStop)
        data=self.fbPolData[int(sensorArray[1])-1]*1
        for i in range(0,len(data)):
            data[i]=data[i][iStart:iStop]*1e4
            
        # create and return plot
        p1=_plot.plot(title=self.title,subtitle="FB "+sensorArray+' Sensors',
                      xLabel='Time [ms]', yLabel='phi [rad]',zLabel='Gauss',
                      plotType='contour',colorMap=_plot._red_green_colormap(),
                      centerColorMapAroundZero=True)
        phi=_np.array(self.phi[int(sensorArray[1])-1])
        index=_np.where(phi==_np.min(phi))[0][0]
        phi[0:index]-=2*_np.pi
        p1.addTrace(self.fbPolTime[iStart:iStop]*1e3,phi,
                    zData=_np.array(data))
        return p1

    def plotOfSinglePol(self, row=0, col=0,plot=True,alsoPlotRawAndFit=True):
        """
        Plots poloidal data from FB sensors
        """
        i=col; j=row;
        
        # initialize plot
        p1=_plot.plot(xLabel='time [ms]',yLabel=r'Gauss',
                      subtitle=self.fbPolNames[j][i],shotno=[self.shotno],
                      title=self.title);
                      
        # smoothed data
        p1.addTrace(yData=self.fbPolData[j][i],xData=self.fbPolTime*1000,
                    yLegendLabel='smoothed')   
        
        if alsoPlotRawAndFit==True:
            # raw data
            p1.addTrace(yData=self.fbPolRaw[j][i],xData=self.fbPolTime*1000,
                        yLegendLabel='Raw')   
            
            # fit data (which is subtracted from raw)
            p1.addTrace(yData=self.fbPolRawFit[j][i],xData=self.fbPolTime*1000,
                        yLegendLabel='Fit')  
            
        return p1
        
    def plotOfSingleRad(self, row=0, col=0,plot=True,alsoPlotRawAndFit=True):
        """
        Plots radial data from FB sensors
        """
        i=col; j=row;
        
        # initialize plot
        p1=_plot.plot(xLabel='time [ms]',yLabel=r'Gauss',
                      subtitle=self.fbRadNames[j][i],shotno=[self.shotno],
                      title=self.title,yLim=[-0.01,0.05]);
        
        # smoothed data
        p1.addTrace(yData=self.fbRadData[j][i],xData=self.fbRadTime*1000,
                    yLegendLabel='smoothed')   
        
        if alsoPlotRawAndFit==True:
            # raw data
            p1.addTrace(yData=self.fbRadRaw[j][i],xData=self.fbRadTime*1000,
                        yLegendLabel='Raw')   
            
            # fit data (which is subtracted from raw)
            p1.addTrace(yData=self.fbRadRawFit[j][i],xData=self.fbRadTime*1000,
                        yLegendLabel='Fit')  
        
        return p1
    
    def plot(self,plotAll=True):
        """
        Plots all 40 poloidal FB sensors
        """
        for i in range(0,4):
            for j in range(0,len(self.fbPolNames[i])):
                
                if plotAll==True:
                    newPlot=self.plotOfSinglePol(i,j,alsoPlotRawAndFit=True)
                else:
                    newPlot=self.plotOfSinglePol(i,j,alsoPlotRawAndFit=False)
                newPlot.plot()

    
@_prepShotno
class taData:
    """
    Downloads toroidal array (TA) sensor data.  Presently, only poloidal
    measurements as the radial sensors are not yet implemeted.  
    
    Parameters
    ----------

    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        True - Plots a sample of each FB poloidal and radial data
        'sample'- same as True
        'all' - Plots all 80 sensor data
#    smoothingAlgorithm : str
#        informs function as to which smoothing algorithm to use on each PA 
#        sensor
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    namesTAPol : list (of str)
        names of poloidal-TA sensors
    namesTARad : list (of str)
        names of radial-TA sensors
    phi : numpy.ndarray
        toroidal location of each poloidal-TA sensor.  units in radians.
    phiR : numpy.ndarray
        toroidal location of each raidal-TA sensor.  units in radians.
    taPolRaw : list (of numpy.ndarray)
        raw poloidal-TA sensor data
    taRadRaw : list (of numpy.ndarray)
        raw radial-TA sensor data
    taPolTime : numpy.ndarray
        time data associated with poloidal-TA sensor data
    taRadTime : numpy.ndarray
        time data associated with radial-TA sensor data
    taPolData : list (of numpy.ndarray)
        poloidal-TA sensor data
    taRadData : list (of numpy.ndarray)
        radial-TA sensor data
    taPolRawFit : list (of numpy.ndarray)
        fit of raw poloidal-TA sensor data.  subtract this from taPolRaw to get
        taPolData
    taRadRawFit : list (of numpy.ndarray)
        fit of raw radial-TA sensor data.  subtract this from taRadRaw to get
        taRadData
        
    Subfunctions
    ------------
    plotOfSinglePol :
        returns plot function of a single poloidal sensor, based on provided 
        index
    plot :
        plots all relevant datas

    """
        
    def __init__(    self,
                shotno=98173,
                tStart=_TSTART,
                tStop=_TSTOP,
                plot=False,
                removeBadSensors=False):
        self.shotno = shotno
        self.title = "%d, TA sensor data." % shotno
        self.badSensors=[] # no bad sensors as of present
        
        # names of poloidal and radial sensors
        self.namesTAPol=['TA01_S1P', 'TA01_S2P', 'TA01_S3P', 'TA02_S1P', 'TA02_S2P', 'TA02_S3P', 'TA03_S1P', 'TA03_S2P', 'TA03_S3P', 'TA04_S1P', 'TA04_S2P', 'TA04_S3P', 'TA05_S1P', 'TA05_S2P', 'TA05_S3P', 'TA06_S1P', 'TA06_S2P', 'TA06_S3P', 'TA07_S1P', 'TA07_S2P', 'TA07_S3P', 'TA08_S1P', 'TA08_S2P', 'TA08_S3P', 'TA09_S1P', 'TA09_S2P', 'TA09_S3P', 'TA10_S1P', 'TA10_S2P', 'TA10_S3P'];
        self.namesTARad=['TA01_S2R', 'TA02_S2R', 'TA03_S2R', 'TA04_S2R', 'TA05_S2R', 'TA06_S2R', 'TA07_S2R', 'TA08_S2R', 'TA09_S2R', 'TA10_S2R']
        
        # toroidal locations for the poloidal measurements
        self.phi=_np.pi/180.*_np.array([241.5,250.5,259.5,277.5,286.5,295.5,313.5,322.5,331.5,349.5,358.5,7.5,25.5,34.5,43.5,61.5,70.5,79.5,97.5,106.5,115.5,133.5,142.5,151.5,169.5,178.5,187.5,205.5,214.5,223.5])    
        
        # poloidal locations of sensors
        self.theta=_np.ones(len(self.phi))*(189-360)*_np.pi/180
        
#        # toroidal locations for the radial measurements
#        self.phiR=_np.pi/180.*_np.array([-108.,  -72.,  -36.,    0.,   36.,   72.,  108.,  144.,  180.,  216.])
        
        # compile full sensor addresses names
        taPolSensorAddresses=[]
        taRadSensorAddresses=[]    
        rootAddress='\HBTEP2::TOP.SENSORS.MAGNETIC:';
        for i in range(0,30):
            taPolSensorAddresses.append(rootAddress+self.namesTAPol[i])
            if i < 10:
                taRadSensorAddresses.append(rootAddress+self.namesTARad[i])
                
        # get raw data
        self.taPolRaw,self.taPolTime=mdsData(shotno,taPolSensorAddresses, tStart, tStop)
        self.taRadRaw,self.taRadTime=mdsData(shotno,taRadSensorAddresses, tStart, tStop)
          
        # data smoothing algorithm
        self.taPolData=[]
        self.taPolRawFit=[]
        self.taRadData=[]
        self.taRadRawFit=[]
        
        # high pass filter the measurements
        for i in range(0,30):
            temp,temp2=_process.gaussianHighPassFilter(self.taPolRaw[i][:],self.taPolTime,timeWidth=1./20000)
            self.taPolData.append(temp)
            self.taPolRawFit.append(temp2)
            
        # pandas dataframes 
        #TODO(John) rewrite entire class.  start with dataframes instead of lists
        self.dfData=_pd.DataFrame(     data=_np.append(_np.array([self.taPolTime]).transpose(),_np.array(self.taPolData).transpose(),axis=1),
                                columns=['Time']+self.namesTAPol).set_index('Time')
        self.dfDataRaw=_pd.DataFrame(     data=_np.append(_np.array([self.taPolTime]).transpose(),_np.array(self.taPolRaw).transpose(),axis=1),
                                columns=['Time']+self.namesTAPol).set_index('Time')
        self.dfMeta=_pd.DataFrame(     data={'SensorNames':self.namesTAPol,
                                    'Phi':self.phi,
                                    'Theta':self.theta,
                                    'Address':taPolSensorAddresses,
                                    'SectionNum':_np.array([ 1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,  5,  6,  6,
                                                           6,  7,  7,  7,  8,  8,  8,  9,  9,  9, 10, 10, 10]),
                                    'CenterSensor':[False,  True, False, False,  True, False, False,  True, False,
                                                   False,  True, False, False,  True, False, False,  True, False,
                                                   False,  True, False, False,  True, False, False,  True, False,
                                                   False,  True, False],
                                    },
                                columns=['SensorNames','Phi','Theta','Address','SectionNum','CenterSensor']).set_index('SensorNames')#,'Address','SectionNum'])
        

        if removeBadSensors==True:
            self.dfData=self.dfData.drop(columns=self.badSensors)
            self.dfDataRaw=self.dfDataRaw.drop(columns=self.badSensors)
            self.dfMeta=self.dfMeta.drop(index=self.badSensors)

        # plot
        if plot=='sample':
            self.plotOfSinglePol().plot();
        elif plot==True or plot=='all':
            self.plot(True);
            
    # TODO Add plotOfSingleRad function
            
    def plotOfTAStripey(self,tStart=2e-3,tStop=4e-3):
        iStart=_process.findNearest(self.taPolTime,tStart)
        iStop=_process.findNearest(self.taPolTime,tStop)
        p1=_plot.plot(title=self.title,subtitle='TA Sensors',
                      xLabel='Time [ms]', yLabel='phi [rad]',zLabel='Gauss',
                      plotType='contour',colorMap=_plot._red_green_colormap(),
                      centerColorMapAroundZero=True)
        
        # phi is not in order.  this corrects for that
        maxValInd=_np.where(self.phi==self.phi.min())[0][0]
#        print(type(self.taPolData[0:30]))
        data=list(_np.concatenate((self.taPolData[maxValInd:30],self.taPolData[0:maxValInd]),axis=0))
        phi=_np.concatenate((self.phi[maxValInd:30],self.phi[0:maxValInd]))

        # truncate data with time
        for i in range(0,len(data)):
            data[i]=data[i][iStart:iStop]*1e4
        p1.addTrace(self.taPolTime[iStart:iStop]*1e3,phi,
                    _np.array(data))
        return p1
            
    def plotOfSinglePol(self, i=0, alsoPlotRawAndFit=True):
        """ Plot one of the PA1 plots.  based on index, i. """
        p1=_plot.plot(xLabel='time [ms]',yLabel=r'Gauss',shotno=[self.shotno],
                      title=self.title,subtitle=self.namesTAPol[i]);
        
        # smoothed data
        p1.addTrace(yData=self.taPolData[i],xData=self.taPolTime*1000,
                    yLegendLabel='smoothed')  

        if alsoPlotRawAndFit==True:
            # raw data
            p1.addTrace(yData=self.taPolRaw[i],xData=self.taPolTime*1000,
                        yLegendLabel='raw') 
            
            # fit data (which is subtracted from raw)
            p1.addTrace(yData=self.taPolRawFit[i],xData=self.taPolTime*1000,
                        yLegendLabel='fit') 

        return p1
        
    def plot(self,plotAll=True):
        """
        Plots poloidal sensor data for all 40 sensors
        
        Warning, 40 plots is tough on memory.  
        """
        # TODO(john) update this so that all poloidal data is on a single 
        # window. Same with radial
        
        sp1=[[],[],[],[],[]]
        count=0
        for i in range(0,5):
            for j in range(0,6):
                if plotAll==True:
                    newPlot=self.plotOfSinglePol(count,alsoPlotRawAndFit=True);
                else:
                    newPlot=self.plotOfSinglePol(count,alsoPlotRawAndFit=False);
                newPlot.subtitle=self.namesTAPol[count]
                newPlot.yLegendLabel=[]
                sp1[i].append(newPlot)
                count+=1;
        sp1[0][0].title=self.title
        sp1=_plot.subPlot(sp1,plot=False)
        sp1.shareY=True;
        sp1.plot()
  

@_prepShotno
def taData_df(shotno=98173,
            tStart=_TSTART,
            tStop=_TSTOP,
            plot=False,
            removeBadSensors=False):
    """
    Downloads toroidal array (TA) sensor data.  Presently, only poloidal
    measurements as the radial sensors are not yet implemeted.  
    
    Parameters
    ----------

    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        True - Plots a sample of each FB poloidal and radial data
        'sample'- same as True
        'all' - Plots all 80 sensor data
#    smoothingAlgorithm : str
#        informs function as to which smoothing algorithm to use on each PA 
#        sensor
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    namesTAPol : list (of str)
        names of poloidal-TA sensors
    namesTARad : list (of str)
        names of radial-TA sensors
    phi : numpy.ndarray
        toroidal location of each poloidal-TA sensor.  units in radians.
    phiR : numpy.ndarray
        toroidal location of each raidal-TA sensor.  units in radians.
    taPolRaw : list (of numpy.ndarray)
        raw poloidal-TA sensor data
    taRadRaw : list (of numpy.ndarray)
        raw radial-TA sensor data
    taPolTime : numpy.ndarray
        time data associated with poloidal-TA sensor data
    taRadTime : numpy.ndarray
        time data associated with radial-TA sensor data
    taPolData : list (of numpy.ndarray)
        poloidal-TA sensor data
    taRadData : list (of numpy.ndarray)
        radial-TA sensor data
    taPolRawFit : list (of numpy.ndarray)
        fit of raw poloidal-TA sensor data.  subtract this from taPolRaw to get
        taPolData
    taRadRawFit : list (of numpy.ndarray)
        fit of raw radial-TA sensor data.  subtract this from taRadRaw to get
        taRadData
        
    Subfunctions
    ------------
    plotOfSinglePol :
        returns plot function of a single poloidal sensor, based on provided 
        index
    plot :
        plots all relevant datas

    """
    
    shotno = shotno
#    title = "%d, TA sensor data." % shotno
    badSensors=[] # no bad sensors as of present
    
    # names of poloidal and radial sensors
    namesTAPol=['TA04_S3P', 'TA05_S1P', 'TA05_S2P', 'TA05_S3P', 'TA06_S1P', 'TA06_S2P', 'TA06_S3P', 'TA07_S1P', 'TA07_S2P', 'TA07_S3P', 'TA08_S1P', 'TA08_S2P', 'TA08_S3P', 'TA09_S1P', 'TA09_S2P', 'TA09_S3P', 'TA10_S1P', 'TA10_S2P', 'TA10_S3P','TA01_S1P', 'TA01_S2P', 'TA01_S3P', 'TA02_S1P', 'TA02_S2P', 'TA02_S3P', 'TA03_S1P', 'TA03_S2P', 'TA03_S3P', 'TA04_S1P', 'TA04_S2P'];
    namesTARad=['TA01_S2R', 'TA02_S2R', 'TA03_S2R', 'TA04_S2R', 'TA05_S2R', 'TA06_S2R', 'TA07_S2R', 'TA08_S2R', 'TA09_S2R', 'TA10_S2R']
    
    # toroidal locations for the poloidal measurements
    phi=_np.pi/180.*_np.array([7.5,25.5,34.5,43.5,61.5,70.5,79.5,97.5,106.5,115.5,133.5,142.5,151.5,169.5,178.5,187.5,205.5,214.5,223.5,241.5,250.5,259.5,277.5,286.5,295.5,313.5,322.5,331.5,349.5,358.5])    
    
    # poloidal locations of sensors
    theta=_np.ones(len(phi))*(189-360)*_np.pi/180
    
#        # toroidal locations for the radial measurements
#        phiR=_np.pi/180.*_np.array([-108.,  -72.,  -36.,    0.,   36.,   72.,  108.,  144.,  180.,  216.])
    
    # compile full sensor addresses names
    taPolSensorAddresses=[]
    taRadSensorAddresses=[]    
    rootAddress='\HBTEP2::TOP.SENSORS.MAGNETIC:';
    for i in range(0,30):
        taPolSensorAddresses.append(rootAddress+namesTAPol[i])
        if i < 10:
            taRadSensorAddresses.append(rootAddress+namesTARad[i])
            
    # get raw data
    taPolRaw,taPolTime=mdsData(shotno,taPolSensorAddresses, tStart, tStop)
    taRadRaw,taRadTime=mdsData(shotno,taRadSensorAddresses, tStart, tStop)
      
    # data smoothing algorithm
    taPolData=[]
    taPolRawFit=[]
    taRadData=[]
    taRadRawFit=[]
    
    # high pass filter the measurements
    for i in range(0,30):
        temp,temp2=_process.gaussianHighPassFilter(taPolRaw[i][:],taPolTime,timeWidth=1./20000)
        taPolData.append(temp)
        taPolRawFit.append(temp2)
        
    # pandas dataframes 
    #TODO(John) rewrite entire class.  start with dataframes instead of lists
    dfData=_pd.DataFrame(     data=_np.append(_np.array([taPolTime]).transpose(),_np.array(taPolData).transpose(),axis=1),
                            columns=['Time']+namesTAPol).set_index('Time')
    dfDataRaw=_pd.DataFrame(     data=_np.append(_np.array([taPolTime]).transpose(),_np.array(taPolRaw).transpose(),axis=1),
                            columns=['Time']+namesTAPol).set_index('Time')
    dfMeta=_pd.DataFrame(     data={'SensorNames':namesTAPol,
                                'Phi':phi,
                                'Theta':theta,
                                'Address':taPolSensorAddresses,
                                'SectionNum':_np.array([ 1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,  5,  6,  6,
                                                       6,  7,  7,  7,  8,  8,  8,  9,  9,  9, 10, 10, 10]),
                                'CenterSensor':[False,  True, False, False,  True, False, False,  True, False,
                                               False,  True, False, False,  True, False, False,  True, False,
                                               False,  True, False, False,  True, False, False,  True, False,
                                               False,  True, False],
                                },
                            columns=['SensorNames','Phi','Theta','Address','SectionNum','CenterSensor']).set_index('SensorNames')#,'Address','SectionNum'])
    

    if removeBadSensors==True:
        dfData=dfData.drop(columns=badSensors)
        dfDataRaw=dfDataRaw.drop(columns=badSensors)
        dfMeta=dfMeta.drop(index=badSensors)

    return dfData,dfDataRaw,dfMeta


@_prepShotno
class groundCurrentData:
    """
    Ground current flowing through the west and north racks to the grounding bus.
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    wRackCurrent : numpy.ndarray
        west rack current to grounding bus
    nRackCurrent : numpy.ndarray
        north rack current to grounding bus
    wRackTime : numpy.ndarray
        time data
    nRackTime : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plot :
        plots data
    
    """
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno = shotno
        self.title = "%d, Ext. Rogowski Data" % shotno
        
        # get north rack data
        data, time=mdsData(shotno=shotno,
                           dataAddress=['\HBTEP2::TOP.DEVICES.NORTH_RACK:CPCI:INPUT_96'],
                           tStart=tStart, tStop=tStop)
        self.nRackCurrent=data[0];
        self.nRackTime=time;
        
        # get west rack data
        data, time=mdsData(shotno=shotno,
                           dataAddress=['\HBTEP2::TOP.DEVICES.WEST_RACK:CPCI:INPUT_96'],
                           tStart=tStart, tStop=tStop)
        self.wRackCurrent=data[0];
        self.wRackTime=time;
        
        if plot == True:
            self.plot().plot()
        
    def plot(self):
        """ Plot all relevant plots """
    
        fig,p1=_plt.subplots()
        p1.plot(self.nRackTime*1e3,self.nRackCurrent,label='North Rack Ground Current')
        p1.plot(self.wRackTime*1e3,self.wRackCurrent,label='West Rack Ground Current')
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='Current (A)')
        _plot.finalizeFigure(fig,title=self.title)
        
        return p1
        

@_prepShotno
class quartzJumperData:
    """
    External rogowski data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    eRogA : numpy.ndarray
        external rogowski A data
    eRogB : numpy.ndarray
        external rogowski B data
    eRogC : numpy.ndarray
        external rogowski C data
    eRogD : numpy.ndarray
        external rogowski D data
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plotOfERogA : _plotTools.plot
        plot of external rogowski A data
    plotOfERogB : _plotTools.plot
        plot of external rogowski B data
    plotOfERogC : _plotTools.plot
        plot of external rogowski C data
    plotOfERogD : _plotTools.plot
        plot of external rogowski D data
    plotOfERogAll : _plotTools.plot
        plot of all 4 external rogowskis
    plot :
        plots plotOfERogAll()
        
    Notes
    -----
    Rog. D is permanently off for the time being
    Rog. C is permanently off for the time being
    
    """
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno = shotno
        self.title = "%d, Ext. Rogowski Data" % shotno
        
        # get data
        dataAddress=['\HBTEP2::TOP.SENSORS.EXT_ROGS:EX_ROG_A',
                                        '\HBTEP2::TOP.SENSORS.EXT_ROGS:EX_ROG_B',
                                        '\HBTEP2::TOP.SENSORS.EXT_ROGS:EX_ROG_C',
                                        '\HBTEP2::TOP.SENSORS.EXT_ROGS:EX_ROG_D',
                                        '\HBTEP2::TOP.DEVICES.WEST_RACK:CPCI:INPUT_96 ',
                                        '\HBTEP2::TOP.DEVICES.NORTH_RACK:CPCI:INPUT_96 ']
        data, time=mdsData(shotno=shotno,
                           dataAddress=dataAddress,
                           tStart=tStart, tStop=tStop)
        self.eRogA=data[0];
        self.eRogB=data[1];
        self.eRogC=data[2];
        self.eRogD=data[3];
        data[4]*=100
        data[5]*=100
        self.westRackGround=data[4]*100
        self.northRackGround=data[5]*100
        self.time=time;
        self.sensorNames=['A. Section 9-10','B. Section 3-4','C. Section 10-1','D. Section 5-6']
        self.phi=_np.array([198,342,234,54])*_np.pi/180.
        self.theta=_np.array([0,0,0,0])
        
        # pandas dataframes
        self.dfData=_pd.DataFrame(     data=_np.array((time,data[0],data[1],data[2],data[3],data[4],data[5])).transpose(),
                                columns=['Time','Jumper9_10','Jumper3_4','Jumper10_1','Jumper5_6','WestRackGround','NorthRackGround']).set_index('Time')
        self.dfMeta=_pd.DataFrame(     data={'SensorNames':['Jumper9_10','Jumper3_4','Jumper10_1','Jumper5_6','WestRackGround','NorthRackGround'],
                                    'Phi':_np.array([198,342,234,54,0,0])*_np.pi/180.,
                                    'Theta':_np.array([0,0,0,0,0,0]),
                                    'Address':dataAddress,
                                    'SectionNum':[9.5,3.5,0.5,5.5,0,0],
                                    'JumperLetter':['A','B','C','D','na','na']},
                                columns=['SensorNames','Phi','Theta','Address','SectionNum','JumperLetter']).set_index('SensorNames')
        
        if plot == True:
            self.plot()
        
    def plot(self):
        """ Plot all relevant plots """
            
        fig,p1=_plt.subplots(2,sharex=True)
        p1[0].plot(self.time*1e3,self.eRogA,label='Rogowski A')
        p1[0].plot(self.time*1e3,self.eRogB,label='Rogowski B')
        p1[0].plot(self.time*1e3,self.eRogC,label='Rogowski C')
        p1[0].plot(self.time*1e3,self.eRogD,label='Rogowski D')
        p1[1].plot(self.time*1e3,self.westRackGround,label='West rack ground')
        p1[1].plot(self.time*1e3,self.northRackGround,label='North rack ground')
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='Current (A)')
        _plot.finalizeFigure(fig,title=self.title)
        
        return p1
        
@_prepShotno
class spectrometerData:
    """
    Spectrometer data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    spect : numpy.ndarray
        spectrometer current data
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plotOfSpect : _plotTools.plot
        plot of sprectrometer data
    plot :
        plots all relevant data
    
    """
    def __init__(self,shotno=98030,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno = shotno
        self.title = "%d, Spectrometer Data" % shotno
        
        # get data
        data, self.time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.SPECTROMETER'],
                              tStart=tStart, tStop=tStop)
        self.spect=data[0];
        
        if plot == True or plot=='all':
            self.plot()
        
    def plotOfSpect(self):
        # generate plot
        
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.spect,label='Spectrometer Intensity')
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='Voltage (V)')
        _plot.finalizeFigure(fig,title=self.title)
        
        return p1
    
    def plot(self):
        """ Plot all relevant plots """
        self.plotOfSpect().plot()
        
        
@_prepShotno    
class usbSpectrometerData:
    """
    USB spectrometer data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    spect : numpy.ndarray
        spectrometer current data
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plotOfSpect : _plotTools.plot
        plot of usb sprectrometer data
    plotOfStripey : _plotTools.plot
        stripey plot of usb spectrometer data
    plot :
        plots all relevant data
    
    """
    def __init__(self,shotno=98415,plot=False):
        self.shotno = shotno
        self.title = "%d, USB Spectrometer Data" % shotno
        
        # get data
        dataAddressRoot = '\HBTEP2::TOP.SENSORS.USB_SPECTROM:SPECTRUM_'
        self.spectrometerArrayNumber=[]
        self.spectrometerData=[]
        for i in range(1,11):
            if i < 10:
                dataAddress='%s0%d' % (dataAddressRoot, i)
            else:
                dataAddress='%s%d' % (dataAddressRoot, i)
            try:
                data, xData=mdsData(shotno=shotno,
                                    dataAddress=dataAddress)
                self.spectrometerArrayNumber.append(i)
                self.spectrometerData.append(data[0])
            except:# _mds.MdsIpException:
                print("usb spectrometer channel %d data does not exist for shot number %d" % (i, shotno))
        
        self.spectrometerArrayNumber=_np.array(self.spectrometerArrayNumber)
        # get wavelength
        yData, xData=mdsData(shotno=shotno,
                            dataAddress='\HBTEP2::TOP.SENSORS.USB_SPECTROM:WAVELENGTH')
        self.wavelength=yData[0]
               
        # plot if requested
        if plot == True:
            self.plotOfSpect()
        if plot == 'all':
            self.plotOfSpect()
            self.plotOfStripey().plot()
            
    def plotOfSpect(self):
        # generate subplot of data
        
        fig,p1=_plt.subplots(len(self.spectrometerArrayNumber),sharex=True)
        for i in range(0,len(self.spectrometerArrayNumber)):
            p1[i].plot(self.wavelength,self.spectrometerData[i],label='Time slice %d'%(self.spectrometerArrayNumber[i]))
        _plot.finalizeSubplot(p1,xlabel='Wavelength [nm]',ylabel='Intensity')
        _plot.finalizeFigure(fig,title=self.title)
        
        return fig

    def plotOfStripey(self):
        # generate stripey plot of data
        p1=_plot.plot(yLabel='Channel',xLabel='Wavelength [nm]',zLabel='Intensity',
                     plotType='contour', shotno=self.shotno,title=self.title)
                     
        p1.addTrace(zData=_np.array(self.spectrometerData),xData=self.wavelength,yData=_np.array(self.spectrometerArrayNumber))
        return p1
        
    def plot(self):
        """ Plot all relevant plots """
        self.plotOfSpect()
        
    
@_prepShotno            
class solData:
    """
    SOL tile sensor data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
    numPointsForSmothing : int
        number of points to be used in removing the offset.  Note that 
        numPointsForSmothing is effectively a high pass filter. the smaller the 
        value, the more aggressive the filter is on high frequencies.
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    sensorNames : list (of str)
        names of each SOL sensor
    solData : list (of numpy.ndarray)
        SOL sensor data with offset subtracted
    solFitData : list (of numpy.ndarray)
        Fits that was used to remove the offset
    solDataRaw : list (of numpy.ndarray)
        Raw SOL data (prior to offset subtraction)
    time : numpy.ndarray
        time data
        
    Subfunctions 
    ------------
    plotOfSingleSensor : _plotTools.plot
        returns plot of single sensor
    plot :
        plots all SOL data
    
    """
    def __init__(self,shotno=98030,tStart=_TSTART,tStop=_TSTOP,plot=False,
            numPointsForSmothing=201):
        # note that numPointsForSmothing is effectively a high pass filter.
        # the larger the value, the more aggressive the filter is on high frequencies.

        # initialize
        self.shotno = shotno
        self.title = "%d, SOL Data" % shotno
        self.sensorNames = ['LFS01_S1', 'LFS01_S2', 'LFS01_S3', 'LFS01_S4', 'LFS01_S5', 'LFS01_S6', 'LFS01_S7', 'LFS01_S8', 'LFS04_S1', 'LFS04_S2', 'LFS04_S3', 'LFS04_S4', 'LFS08_S1', 'LFS08_S2', 'LFS08_S3', 'LFS08_S4', 'LFS08_S5', 'LFS08_S6', 'LFS08_S7', 'LFS08_S8']
        self.phis=_np.array([234.8, 234.8, 234.8, 234.8, 234.8, 234.8, 234.8, 234.8, 342.8, 342.8, 342.8, 342.8, 126.8, 126.8, 126.8, 126.8, 126.8, 126.8, 126.8, 126.8])
        self.thetas=_np.array([-70. , -50. , -30. , -10. ,  10. ,  30. ,  50. ,  70. , -83. ,       -28.2,  28.2,  83. , -70. , -50. , -30. , -10. ,  10. ,  30. ,  50. ,  70. ])
        sensorPathRoot='\HBTEP2::TOP.SENSORS.SOL:'
        
        # compile list of sensor addresses for all 20 SOL tiles
        sensorAddress=[]
        for i in range(0,len(self.sensorNames)):
            sensorAddress.append(sensorPathRoot+'%s' % self.sensorNames[i]) 
            
        # get raw data from the tree
        self.solDataRaw, self.time=mdsData(shotno=shotno,
                              dataAddress=sensorAddress,
                              tStart=tStart, tStop=tStop)
                              
        # subtract offset from sensors
        self.solDataFit=[]
        self.solData=[]
        for i in range(0,len(self.sensorNames)):
            temp,temp2=_process.gaussianHighPassFilter(self.solDataRaw[i],self.time,timeWidth=1./20000,plot=False)
            self.solData.append(temp)
            self.solDataFit.append(temp2)
            
        # pandas dataframes
        self.dfData=_pd.DataFrame(     data=_np.append(_np.array([self.time]).transpose(),_np.array(self.solData).transpose(),axis=1),
                                columns=['Time']+self.sensorNames).set_index('Time')
        self.dfMeta=_pd.DataFrame(     data={'SensorNames':self.sensorNames,
                                    'Phi':self.phis*_np.pi/180,
                                    'Theta':self.thetas*_np.pi/180,
                                    'Address':sensorAddress},
                                columns=['SensorNames','Phi','Theta','Address']).set_index('SensorNames')
        self.dfData['AllTotal']=self.dfData.sum(axis=1)
        self.dfData['S01Total']=self.dfData.iloc[:,self.dfData.columns.str.contains('S01')].sum(axis=1)
        self.dfData['S04Total']=self.dfData.iloc[:,self.dfData.columns.str.contains('S04')].sum(axis=1)
        self.dfData['S08Total']=self.dfData.iloc[:,self.dfData.columns.str.contains('S08')].sum(axis=1)
        
        # optional plotting    
        if plot == True:
            self.plot()
        if plot == 'all':
            self.plot('all')
                              
    def plotOfSingleSensor(self,index,plot='all'): #name='LFS01_S1'
        """ 
        Returns plot of a single sol sensor.  Plots raw, fit, and smoothed
        """
        p1=_plot.plot(yLabel='V',xLabel='time [ms]',
                      subtitle=self.sensorNames[index],title=self.title,
                      shotno=self.shotno)
        if plot=='all' or plot=='raw':
            p1.addTrace(yData=self.solDataRaw[index],xData=self.time*1000,
                        yLegendLabel=self.sensorNames[index]+' Raw')
        if plot=='all' or plot=='fit': 
            p1.addTrace(yData=self.solDataFit[index],xData=self.time*1000,
                        yLegendLabel=self.sensorNames[index]+' Fit') 
        if plot=='all' or plot=='smoothed' or plot=='smoothedOnly': 
            p1.addTrace(yData=self.solData[index],xData=self.time*1000,
                        yLegendLabel=self.sensorNames[index]+' Without Offset') 
        return p1
                
    def plotOfContour(self,tStart=2e-3,tStop=4e-3,section='LFS01'):
        """ 
        contour plot of LFS01 Data
        """
        iStart=_process.findNearest(self.time,tStart)
        iStop=_process.findNearest(self.time,tStop)
        p1=_plot.plot(title=self.title,subtitle=section+' SOL sensors',
                      xLabel='Time [ms]', yLabel='phi [rad]',zLabel='A',
                      plotType='contour')
        if section=='LFS01':
            data=self.solData[0:8]
        elif section=='LFS04':
            data=self.solData[8:12]
        elif section=='LFS08':
            data=self.solData[12:20]
        elif section=='all':
            data=self.solData[0:20]
            
        for i in range(0,len(data)):
            data[i]=data[i][iStart:iStop]
#        return data
        p1.addTrace(self.time[iStart:iStop]*1e3,_np.arange(0,8),_np.array(data))
        return p1
        
    def plot(self,plot='smoothedOnly',includeBP=True):
        """ plots all 20 sol sensor currents on three plots """

#        if plot=='all':
#            for j in range(0,20):
#                    p1=self.plotOfSingleSensor(j,'all').plot()
 
        if True:
            fig,ax=_plt.subplots(4,sharex=True)
            
            time=self.dfData.index.to_numpy()*1e3
            dfS01=self.dfData.iloc[:,self.dfData.columns.str.contains('S01')]
            for key,values in dfS01.iteritems():
                ax[0].plot(time,dfS01[key],label=key)
            _plot.finalizeSubplot( ax[0],
                                     ylabel='A',
                                     subtitle='S01')
            
            dfS04=self.dfData.iloc[:,self.dfData.columns.str.contains('S04')]
            for key,values in dfS04.iteritems():
                ax[1].plot(time,dfS04[key],label=key)
            _plot.finalizeSubplot( ax[1],
                                     ylabel='A',
                                     subtitle='S04')
            
            dfS08=self.dfData.iloc[:,self.dfData.columns.str.contains('S08')]
            for key,values in dfS08.iteritems():
                ax[2].plot(time,dfS08[key],label=key)
            _plot.finalizeSubplot( ax[2],
                                     ylabel='A',
                                     subtitle='S08')
            
            dfTotal=self.dfData.iloc[:,self.dfData.columns.str.contains('Total')]
            for key,values in dfTotal.iteritems():
                ax[3].plot(time,dfTotal[key],label=key)
            _plot.finalizeSubplot( ax[3],
                                     ylabel='A',
                                     subtitle='Total')
#            for j in range(0,8):
#                if j==0:
#                    p1=self.plotOfSingleSensor(j,plot) 
#                    p3=self.plotOfSingleSensor(12+j,plot) 
#                    if j<4:
#                        p2=self.plotOfSingleSensor(8+j,plot) 
#                else:
#                    p1.mergePlots(self.plotOfSingleSensor(j,plot))
#                    p3.mergePlots(self.plotOfSingleSensor(12+j,plot))
#                    if j<4:
#                        p2.mergePlots(self.plotOfSingleSensor(8+j,plot))     
#            p1.subtitle='Section 1 SOL Sensors'    
#            p2.subtitle='Section 4 SOL Sensors'    
#            p3.subtitle='Section 8 SOL Sensors'            
#            return _plot.subPlot([p1,p2,p3],plot=True)
        

@_prepShotno
class loopVoltageData:
    """
    loo voltage data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots        
    loopVoltage : numpy.ndarray
        SOL sensor data
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plotOfLV : _plotTools.plot
        returns plot of loop voltage data
    plot : 
        plots loop voltage data
    
    """
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno = shotno
        self.title = "%d, Loop voltage data." % shotno
        
        # get data
        data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.LOOP_VOlTAGE'],
                              tStart=tStart, tStop=tStop)   
        self.loopVoltage=data[0];
        self.time=time;
        
        if plot == True or plot=='all':
            self.plot()
        
    def plotOfLoopVoltage(self):
        # generate plot
        
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.loopVoltage)
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='Voltage (V)',ylim=[0,15])
        _plot.finalizeFigure(fig,title=self.title)
        
        return fig
            
    def plot(self):
        """ Plot all relevant plots """
        self.plotOfLoopVoltage()
        
        
@_prepShotno
class tfData:
    """
    Toroidal field data  
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots        
    tfBankField : numpy.ndarray
        Toroidal mangetic field data
    time : numpy.ndarray
        Toroidla field time data
        
    Subfunctions
    ------------
    plotOfTF : _plotTools.plot
        returns plot of TF data
    plot :
        plots all relevant data
    upSample :
        upsamples TF's normal A14 time base (1e-5 s period) to the CPCI time 
        base (2e-6 s) using a linear interpolation method
    
    Notes
    -----
    note that the TF field data is recorded on an A14 where most of HBTEP data
    is stored with the CPCI.  Because the A14 has a slower sampling rate, this
    means that the TF data has fewer points than the rest of the HBTEP data, 
    and this makes comparing data difficult.  Therefore by default, I up-sample
    the data to match the CPCI sampling rate.  
    """
    def __init__(self,shotno=96530,tStart=None,tStop=None,plot=False,
                 upSample=True):
        self.shotno = shotno
        self.title = "%d, TF Field Data" % shotno
        
        # get tf data
        data, self.time=mdsData(shotno=shotno,
                                  dataAddress=['\HBTEP2::TOP.SENSORS.TF_PROBE'],
                                  tStart=tStart, tStop=tStop) 
        self.tfBankField=data[0];

        if upSample==True:
            self.upSample()
        if plot == True:
            self.plot()
            
    def upSample(self):
        
        # time step sizes
        dtUp=2*1e-6 # CPCI sampling period
        dtDown=self.time[-1]-self.time[-2] # A14 sampling period
        
        # reconstruct CPCI time base
        upTime=_np.arange(self.time[0],self.time[-1]+dtDown-dtUp,dtUp) # note that there is some trickery here with reconstructing the CPCI time base.  
        
        # upsample data
        self.tfBankField=_process.upSampleData(upTime,self.time,self.tfBankField)
        self.time=upTime
            
               
    def plotOfTF(self,tStart=None,tStop=None):
        # generate tf plot
        p1=_plot.plot(yLabel='T',xLabel='time [ms]',subtitle='TF Bank Field',
                      title=self.title,shotno=self.shotno)
        p1.addTrace(yData=self.tfBankField,xData=self.time*1000) 
        return p1
            
    def plot(self):
        """ Plot all relevant plots """
        self.plotOfTF().plot()
        
    
@_prepShotno    
class capBankData:
    """
    Capacitor bank data.  Currents.  
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots       
        Toroidla field time data
    vfBankCurrent : numpy.ndarray
        Vertical Field (VF) bank current data
    vfTime : numpy.ndarray
        VF time data
    ohBankCurrent : numpy.ndarray
        Ohmic heating (OH) bank current data
    ohTime : numpy.ndarray
        OH time data
    shBankCurrent : numpy.ndarray
        SHaping (SH) bank current data
    shTime : numpy.ndarray
        SH time data
        
    Subfunctions
    ------------
    plotOfVF : _plotTools.plot
        returns plot of VF data
    plotOfOH : _plotTools.plot
        returns plot of OH data
    plotOfSH : _plotTools.plot
        returns plot of SH data
    plot :
        plots all relevant data
    
    Notes
    -----
    Note that all 3 banks have their own time array.  This is because the 
    data doesn't always have the same length and therefore must have their own
    time array. 
    
    Note that tStart and tStop are intentionally left as None because the TF
    data is so incredibly long next to the other data.
    """
    
    def __init__(self,shotno=96530,tStart=None,tStop=None,plot=False):
        self.shotno = shotno
        self.title = "%d, Capacitor Bank Data" % shotno
        
        # get vf data
        data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.VF_CURRENT'],
                              tStart=tStart, tStop=tStop) 
        self.vfBankCurrent=data[0];    
        self.vfTime=time;
        
        # get oh data
        data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.OH_CURRENT'],
                              tStart=tStart, tStop=tStop) 
        self.ohBankCurrent=data[0];    
        self.ohTime=time;
        
        # get sh data
        data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.SH_CURRENT'],
                              tStart=tStart, tStop=tStop) 
        self.shBankCurrent=data[0];    
        self.shTime=time;

        if plot == True:
            self.plot()
            
    def plot(self):
        """ Plot all relevant plots """
        
        tf=tfData(self.shotno,tStart=None,tStop=None)
        
        _plt.figure()
        ax1 = _plt.subplot2grid((3,2), (0,1), rowspan=3)  #tf
        ax2 = _plt.subplot2grid((3,2), (0,0)) #vf
        ax3 = _plt.subplot2grid((3,2), (1,0),sharex=ax2) #oh
        ax4 = _plt.subplot2grid((3,2), (2, 0),sharex=ax2) #sh
        fig=_plt.gcf()
        fig.set_size_inches(10,5)
                
        tStart=-2
        tStop=20
        
        ax1.plot(tf.time*1e3,tf.tfBankField)
        ax1.axvspan(tStart,tStop,color='r',alpha=0.3)
        _plot.finalizeSubplot(ax1,xlabel='Time (s)',xlim=[-150,450],ylabel='TF Field (T)')#,title=self.title
        
        ax2.plot(self.vfTime*1e3,self.vfBankCurrent*1e-3)
        _plot.finalizeSubplot(ax2,ylabel='VF Current\n(kA)')
        
        ax3.plot(self.ohTime*1e3,self.ohBankCurrent*1e-3)
        _plot.finalizeSubplot(ax3,ylim=[-20,30],ylabel='OH Current\n(kA)')
        
        ax4.plot(self.shTime*1e3,self.shBankCurrent*1e-3)
        _plot.finalizeSubplot(ax4,ylim=[tStart,tStop],xlabel='Time (s)',ylabel='SH Current\n(kA)')
        
        _plot.finalizeFigure(fig,title=self.title)
#        fig.set_tight_layout(True)
        
        return fig
        
        
#####################################################
        
@_prepShotno
class plasmaRadiusData:
    """
    Calculate the major and minor radius.
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    majorRadius : numpy.ndarray
        plasma major radius in meters
    minorRadius : numpy.ndarray
        plasma minor radius in meters
    time : numpy.ndarray
        time (in seconds) associated with data
        
    Subfunctions
    ------------
    plotOfMajorRadius : 
        returns the plot of major radius vs time
    plotOfMinorRadius : 
        returns the plot of major radius vs time
    plot :
        Plots all relevant plots
        
    Notes
    -----
    The radius calculations below are pulled from Paul Hughes's 
    pauls_MDSplus_toolbox.py code.  In that code, he attributes Niko Rath for 
    its implementation
    
    """
    
    def __init__(self,shotno=95782,tStart=_TSTART,tStop=_TSTOP, plot=False, probeRadius=[]):
        self.shotno=shotno;
        self.title = "%d, plasma radius" % shotno
        
        # Determined by Daisuke during copper plasma calibration
        a=.00643005
        b=-1.10423
        c=48.2567
        
        # Calculated by Jeff, but still has errors
        vf_pickup = 0.0046315133 * -1e-3
        oh_pickup = 7.0723416e-08
        
        # get vf and oh data
        capBank=capBankData(shotno,tStart=tStart,tStop=tStop)
        vf=capBank.vfBankCurrent
        oh=capBank.ohBankCurrent
        self.time=capBank.vfTime
        
        # get plasma current
        ip=ipData(shotno,tStart=tStart,tStop=tStop)
        ip=ip.ip*1212.3*1e-9  # ip gain
        
        # get cos-1 raw data
        cos1=cos1RogowskiData(shotno,tStart=tStart,tStop=tStop+2e-06) # note that the cumtrapz function below loses a data point.  by adding 2e-06 to the time, i start with an additional point that it's ok to lose
        # subtract offset
        cos1Raw=cos1.cos1Raw-cos1.cos1RawOffset        
        
        # integrate cos-1 raw 
        from scipy.integrate import cumtrapz
        cos1 = cumtrapz(cos1Raw,cos1.time) + cos1Raw[:-1]*.004571
        
        # r-major calculations
        pickup = vf * vf_pickup + oh * oh_pickup
        ratio = ip / (cos1 - pickup)
        arg = b**2 - 4 * a * (c-ratio)
        arg[arg < 0] = 0
        r_major = (-b + _np.sqrt(arg)) / (2*a)
        self.majorRadius  = r_major / 100 # Convert to meters
#        self.majorRadius -= 0.45/100
        
        # r-minor calculations
        self.minorRadius=_np.ones(len(self.majorRadius))*0.15
        outwardLimitedIndices=self.majorRadius > (0.92)
        self.minorRadius[outwardLimitedIndices] = 1.07 - self.majorRadius[outwardLimitedIndices] # Outboard limited
        inwardLimitedIndices=self.majorRadius < (0.92 - 0.01704)   
        self.minorRadius[inwardLimitedIndices] = self.majorRadius[inwardLimitedIndices] - 0.75296 # inward limited
        
        if plot==True:
            self.plot();

    def plot(self,plotAll=False):
        fig,p1=_plt.subplots(2,sharex=True)
        p1[0].plot(self.time*1e3,self.majorRadius*1e2,label='Major Radius')
        p1[0].plot(self.time*1e3,_np.ones(_np.shape(self.time*1e3))*92)
        p1[0].plot(self.time*1e3,_np.ones(_np.shape(self.time*1e3))*90)
        p1[1].plot(self.time*1e3,self.minorRadius*1e2,label='Minor Radius')
        _plot.finalizeSubplot(p1[0],xlabel='Time (ms)',ylabel='Major Radius (cm)',ylim=[89,95])
        _plot.finalizeSubplot(p1[1],ylabel='Minor Radius (cm)',ylim=[10,16])
        _plot.finalizeFigure(fig,title=self.title)
        
        return p1
    
def plasmaRadiusData_df(shotno=95782,tStart=_TSTART,tStop=_TSTOP, plot=False, probeRadius=[]):
    """
    Calculate the major and minor radius.
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    majorRadius : numpy.ndarray
        plasma major radius in meters
    minorRadius : numpy.ndarray
        plasma minor radius in meters
    time : numpy.ndarray
        time (in seconds) associated with data
        
    Subfunctions
    ------------
    plotOfMajorRadius : 
        returns the plot of major radius vs time
    plotOfMinorRadius : 
        returns the plot of major radius vs time
    plot :
        Plots all relevant plots
        
    Notes
    -----
    The radius calculations below are pulled from Paul Hughes's 
    pauls_MDSplus_toolbox.py code.  In that code, he attributes Niko Rath for 
    its implementation
    
    """
    
#    def __init__(self,
#        self.shotno=shotno;
#        self.title = "%d, plasma radius" % shotno
        
    # Determined by Daisuke during copper plasma calibration
    a=.00643005
    b=-1.10423
    c=48.2567
    
    # Calculated by Jeff, but still has errors
    vf_pickup = 0.0046315133 * -1e-3
    oh_pickup = 7.0723416e-08
    
    # get vf and oh data
    capBank=capBankData(shotno,tStart=tStart,tStop=tStop)
    vf=capBank.vfBankCurrent
    oh=capBank.ohBankCurrent
    time=capBank.vfTime
    
    # get plasma current
    ip=ipData(shotno,tStart=tStart,tStop=tStop)
    ip=ip.ip*1212.3*1e-9  # ip gain
    
    # get cos-1 raw data
    cos1=cos1RogowskiData(shotno,tStart=tStart,tStop=tStop+2e-06) # note that the cumtrapz function below loses a data point.  by adding 2e-06 to the time, i start with an additional point that it's ok to lose
    # subtract offset
    cos1Raw=cos1.cos1Raw-cos1.cos1RawOffset        
    
    # integrate cos-1 raw 
    from scipy.integrate import cumtrapz
    cos1 = cumtrapz(cos1Raw,cos1.time) + cos1Raw[:-1]*.004571
    
    # r-major calculations
    pickup = vf * vf_pickup + oh * oh_pickup
    ratio = ip / (cos1 - pickup)
    arg = b**2 - 4 * a * (c-ratio)
    arg[arg < 0] = 0
    r_major = (-b + _np.sqrt(arg)) / (2*a)
    majorRadius  = r_major / 100 # Convert to meters
#        self.majorRadius -= 0.45/100
        
    # r-minor calculations
    minorRadius=_np.ones(len(majorRadius))*0.15
    outwardLimitedIndices=majorRadius > (0.92)
    minorRadius[outwardLimitedIndices] = 1.07 - majorRadius[outwardLimitedIndices] # Outboard limited
    inwardLimitedIndices=majorRadius < (0.92 - 0.01704)   
    minorRadius[inwardLimitedIndices] = majorRadius[inwardLimitedIndices] - 0.75296 # inward limited
    
    dfData=_pd.DataFrame()
    dfData['time']=time
    dfData['minorRadius']=minorRadius
    dfData['majorRadius']=majorRadius
    dfData=dfData.set_index('time')
    
    return dfData
#    
#        if plot==True:
#            self.plot();
#
#    def plot(self,plotAll=False):
#        fig,p1=_plt.subplots(2,sharex=True)
#        p1[0].plot(self.time*1e3,self.majorRadius*1e2,label='Major Radius')
#        p1[0].plot(self.time*1e3,_np.ones(_np.shape(self.time*1e3))*92)
#        p1[0].plot(self.time*1e3,_np.ones(_np.shape(self.time*1e3))*90)
#        p1[1].plot(self.time*1e3,self.minorRadius*1e2,label='Minor Radius')
#        _plot.finalizeSubplot(p1[0],xlabel='Time (ms)',ylabel='Major Radius (cm)',ylim=[89,95])
#        _plot.finalizeSubplot(p1[1],ylabel='Minor Radius (cm)',ylim=[10,16])
#        _plot.finalizeFigure(fig,title=self.title)
#        
#        return p1


@_prepShotno
class qStarData:
    """
    Gets qstar data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    qStar : numpy.ndarray
        plasma current data
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plotOfQStar : 
        returns the plot of IP vs time
    plot :
        Plots all relevant plots
    """
    
    def __init__(self,shotno=96496, tStart=_TSTART, tStop=_TSTOP, plot=False):
        self.shotno = shotno
        self.title = r"%d, q$^*$ Data" % shotno
        
        # get data
        ip=ipData(shotno,tStart=tStart,tStop=tStop)
        plasmaRadius=plasmaRadiusData(shotno,tStart=tStart,tStop=tStop)
#        tfProbeData,tfProbeTime=mdsData(shotno=96496,
        tfProbeData,tfProbeTime=mdsData(shotno,
                                        dataAddress=['\HBTEP2::TOP.SENSORS.TF_PROBE'],
                                        tStart=tStart,tStop=tStop)
                                        
        # upsample tfprobe data  (its recorded on the A14)      
        data=_process.upSampleData(ip.time,tfProbeTime,tfProbeData[0])
        
        # more tf calculations
        tfProbeData=data*1.23/plasmaRadius.majorRadius
        
        # calc q star
        self.qStar= plasmaRadius.minorRadius**2 * tfProbeData / (2e-7 * ip.ip * plasmaRadius.majorRadius)
        self.qStarCorrected=self.qStar*(1.15) # 15% correction factor.  jeff believes our qstar measurement might be about 15% to 20% too low.  
        self.time=ip.time
        
        if plot == True:
            self.plot()
            
    def plot(self):
        """ 
        Plot all relevant plots 
        """
        
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.qStar,label=r'q$^*$')
        p1.plot(self.time*1e3,self.qStarCorrected,label=r'q$^* * 1.15$')
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel=r'q$^*$',ylim=[1,5])
        _plot.finalizeFigure(fig,title=self.title)
       
       
###############################################################################
### sensor black list data.  presently not used anywhere
    
def checkBlackList_depricated(inData,inName):
    # TODO(John) this needs an overhaul
    """
    Takes in data and sensor name.  Checks sensor name against blacklist.  
    If bad sensor, return all zeros # and a true boolean.  
    Otherwise, returns original data # and a false boolean.
    
    """
    if inName in _SENSORBLACKLIST==True:
        outData=_np.zeros(inData.size);
    else:
        outData=inData;
        
    return outData
    
    
       
  
###############################################################################
### Processed data from HBTEP
    
@_prepShotno
class nModeData:
    """
    This function performs n-mode (toroidal) mode analysis on the plasma.
    Provides mode amplitude, phase, and frequency
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        True - plots relevant data
        'all' - plots all data
    nModeSensor : str
        sensors to be used to calculate the modes
        'FB' - feedback sensors
        'TA' - toroidal array sensors
    method : str
        method to calculate mode analysis
        'leastSquares' - performs a matrix least squares analysis
        
    Attributes
    ----------
    shotno : int
        data shot number
    title : str
        title to be placed on each plot
    nModeSensor : str
        sensors to be used to calculate the modes
    n1Amp : numpy.ndarray
        n=1 mode amplitude data
    n2Amp : numpy.ndarray
        n=2 mode amplitude data
    n1Phase : numpy.ndarray
        filtered n=1 mode phase data
    n1PhaseRaw : numpy.ndarray
        raw n=1 mode phase data
    n1Freq : numpy.ndarray
        filtered n=1 mode frequency data
    n1FreqRaw : numpy.ndarray
        raw n=1 mode frequency data
        
    Subfunctions
    ------------
    plot :
        plots relevant plots
    Bn1 : 
        Generates a pretend B_{n=1} signal at the toroidal location, phi0
        Not presently in use
    plotOfAmps : 
        returns a plot of n=1 and n=2 mode amplitudes
    plotOfN1Phase
        returns a plot of the n=1 mode phase
    self.plotOfN1Freq
        returns a plot of the n=1 mode frequency
    self.plotOfN1Amp
        returns a plot of the n=1 mode amplitude
        
    Notes
    -----
    The convolution filters used with the phase and frequency "mess up" the 
    last tenth of a millisecond of data
    
    """    

    def Bn1(self,phi0=0):
        """
        Generates a pretend B_{n=1} signal at the toroidal location, phi0
        Not presently in use
        """
        return self.x[1,:]*_np.sin(self.phi0)+self.x[2,:]*_np.cos(self.phi0)
        
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False,
                 nModeSensor='FB',method='leastSquares',phaseFilter='gaussian',
                 phaseFilterTimeConstant=1.0/100e3):
        
        self.shotno=shotno
        self.title = '%d.  %s sensor.  n mode analysis' % (shotno,nModeSensor)
        self.nModeSensor=nModeSensor        
        
        # load data from requested sensor array
        if nModeSensor=='TA':
            ## load TA data
            temp=taData(self.shotno,tStart,tStop+0.5e-3);  # asking for an extra half millisecond (see Notes above) 
            data=temp.taPolData
            self.time=temp.taPolTime
            phi=temp.phi
            [n,m]=_np.shape(data)
        elif nModeSensor=='FB' or nModeSensor=='FB_S4':
            ## load FB data
            array=3 # the 4th array (4-1=3) is the top most FB array and has no broken sensors
            temp=fbData(self.shotno,tStart=tStart,tStop=tStop+0.5e-3);  # asking for an extra half millisecond (see Notes above) 
            data=temp.fbPolData[array]  ## top toroidal array = 0, bottom = 3
            self.time=temp.fbPolTime
            phi=_np.array(temp.phi[array])
            self._theta=_np.array(temp.theta[array])
            [n,m]=_np.shape(data)
        self._data=data
        self._phi=phi

        if method=='leastSquares':
            ## Construct A matrix and its inversion.  Only interested in n=1 and n=2 mode fits at present
            A=_np.zeros((n,5))
            A[:,0]=_np.ones(n);
            A[:,1]=_np.sin(phi)
            A[:,2]=_np.cos(phi) 
            A[:,3]=_np.sin(2*phi)
            A[:,4]=_np.cos(2*phi)
            Ainv=_np.linalg.pinv(A)
            
            ## Solve for coefficients, x, for every time step and assign values to appropriate arrays 
            x=_np.zeros([5,m]);
            self.n1Amp=_np.zeros(m)
            self.n1PhaseRaw=_np.zeros(m)
            self.n2Amp=_np.zeros(m)
            
            x=Ainv.dot(data)
            self.n1Amp=_np.sqrt(x[1,:]**2+x[2,:]**2)*1e4
            self.n2Amp=_np.sqrt(x[3,:]**2+x[4,:]**2)*1e4
            self.n1PhaseRaw=_np.arctan2(x[1,:],x[2,:])
            self.n2PhaseRaw=_np.arctan2(x[3,:],x[4,:])
#            for j in range(0,m):
#                y=_np.zeros(n);
#                for i in range(0,n):
#                    y[i]=data[i][j]*1e4
#                x[:,j]=Ainv.dot(y)
#                self.n1Amp[j]=_np.sqrt(x[1,j]**2+x[2,j]**2)
#                self.n2Amp[j]=_np.sqrt(x[3,j]**2+x[4,j]**2)
#                self.n1PhaseRaw[j]=_np.arctan2(x[1,j],x[2,j])
            self._x=x
            self.n1PhaseRaw*=-1  # for some reason, the slope of phase had the wrong sign.  this corrects that.
            self.n2PhaseRaw*=-1  # for some reason, the slope of phase had the wrong sign.  this corrects that.


        else:
            _sys.exit("Invalid mode analysis method requested.")
            
        # filter phase.  (this is necessary to get a clean frequency)
        self.n1Phase=_np.zeros(len(self.n1PhaseRaw))   
        if phaseFilter == 'gaussian':
            self.n1Phase=_process.wrapPhase(
                    _process.gaussianLowPassFilter(
                            _process.unwrapPhase(self.n1PhaseRaw),
                            self.time,
                            timeWidth=phaseFilterTimeConstant))
        else:
            _sys.exit("Invalid phase filter requested.")
                    
        ## Calculate frequency (in Hz) using second order deriv 
        self.n1Freq=_np.gradient(_process.unwrapPhase(self.n1Phase))/_np.gradient(self.time)/(2*_np.pi)        
            
        # trim off extra half millisecond (see Notes)
        self.time, temp=_trimTime(self.time,
                                  [self.n1Amp,self.n2Amp,self.n1Phase,
                                   self.n1PhaseRaw,self.n1Freq],
                                  tStart,tStop)
        self.n1Amp=temp[0]
        self.n2Amp=temp[1]
        self.n1Phase=temp[2]
        self.n1PhaseRaw=temp[3]
        self.n1Freq=temp[4]
        
        self.dfData=_pd.DataFrame(     data=_np.concatenate((    [_np.array(self.time)],
                                                    [_np.array(self.n1Amp).transpose()],
                                                    [_np.array(self.n1Phase).transpose()],
                                                    [_np.array(self.n1PhaseRaw).transpose()],
                                                    [_np.array(self.n1Freq).transpose()])).transpose(),
                                columns=['Time','n1Amp','n1Phase','n1PhaseRaw','n1Freq']).set_index('Time')
        self.dfMeta=_pd.DataFrame() # intentionally left empty

        
        ## plot data
        if plot==True:
            self.plot()
            
        elif plot == 'all':
            self.plotOfSlice(index=int(m/4)).plot();
            self.plotOfSlice(index=int(m/2)).plot();
            self.plot(True)
            
    def plot(self,plotAll=False):
        """
        plots and returns a subplot of n=1 mode amplitude, phase, and frequency
        
        Parameters
        ----------
        includeRaw : bool
            if True, also plots the raw (unfiltered) phase and frequency
        """
        fig,p1=_plt.subplots(3,sharex=True)
        p1[0].plot(self.time*1e3,self.n1Amp,label=r'$\left| \delta B\right|_{n=1}$')
        p1[1].plot(self.time*1e3,self.n1Phase,label=r'$\angle \delta B_{n=1}$',marker='.',markersize=2,linestyle='')
        p1[2].plot(self.time*1e3,self.n1Freq*1e-3,label='Frequency (kHz)')
        if plotAll==True:
            p1[0].plot(self.time*1e3,self.n2Amp,label=r'$\left| \delta B\right|_{n=2}$')
            p1[1].plot(self.time*1e3,self.n1PhaseRaw,label=r'$\angle \delta B_{n=1}$ unfiltered',marker='.',markersize=2,linestyle='')
        _plot.finalizeSubplot(p1[0],ylabel='Mode amplitude (G)',ylim=[-0.5,12])#,title=self.title)#xlabel='Time (ms)',
        _plot.finalizeSubplot(p1[1],ylabel='Mode phase (rad)',ylim=[-_np.pi,_np.pi]) #xlabel='Time (ms)',
        _plot.finalizeSubplot(p1[2],ylabel='Frequency (kHz)',xlabel='Time (ms)')#,ylim=[-_np.pi,_np.pi]
        _plot.finalizeFigure(fig,self.title)
        
#        
#    def plotOfPhaseAmp(self):
#                   
#        # hybrid plot of phase AND amplitude
#        # TODO(John) implement in new plot function
#                   # TODO(john) subfunction needs overhaul
#        p1=_plot.plot() 
#        p1.yData=[self.n1Phase]
#        p1.xData=[self.time*1000]
#        p1.colorData=[self.n1Amp]#[self.n1Amp]
#        p1.linestyle=['']
#        p1.marker=['.']
#        p1.subtitle='n=1 Phase and Filtered Amplitude'
#        p1.title=self.title
#        p1.shotno=[self.shotno]
#        p1.xLabel='ms'
#        p1.yLabel=r'$\phi$'
#        p1.zLabel='Gauss'
#        p1.yLegendLabel=['TA sensors']
#        p1.plotType='scatter'
#        p1.yLim=[-_np.pi,_np.pi]
#        return p1      
##        mx=_np.max(self.n1AmpFiltered)
##        lCutoff=2.5
##        uCutoff=8.
##        cm = _processt.singleColorMapWithLowerAndUpperCutoffs(lowerCutoff=lCutoff/mx,upperCutoff=uCutoff/mx)
##        self.plotOfPhaseAmp.cmap=cm
                        
    def plotOfSlice(self,index=0):
        """
        Plots fit data for a single time value
        """
        j=index;
        [n,m]=_np.shape(self._data)
        y=_np.zeros(n);
        for i in range(0,n):
                y[i]=self._data[i][j]*1e4
        p1=_plot.plot(shotno=[self.shotno],
                      title=self.title+', t='+str(self.time[j]*1000)+'ms.')
        phi=_np.linspace(0,_np.pi*2,100)
        n1Fit=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)
        n2Fit=self._x[0,j]+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)
        fitTotal=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)

        # plot
        p1.addTrace(yData=y,xData=self._phi,
                    marker='x',linestyle='',yLegendLabel='raw') 
        p1.addTrace(yData=n1Fit,xData=phi,
                    yLegendLabel='n=1') 
        p1.addTrace(yData=n2Fit,xData=phi,
                    yLegendLabel='n=2') 
        p1.addTrace(yData=fitTotal,xData=phi,
                    yLegendLabel='Superposition') 
        return p1
    
    

# @_prepShotno
# class nModeData_df:
#     """
#     This function performs n-mode (toroidal) mode analysis on the plasma.
#     Provides mode amplitude, phase, and frequency
#     
#     Parameters
#     ----------
#     shotno : int
#         shot number of desired data
#     tStart : float
#         time (in seconds) to trim data before
#         default is 0 ms
#     tStop : float
#         time (in seconds) to trim data after
#         default is 10 ms
#     plot : bool
#         plots all relevant plots if true
#         default is False
#         True - plots relevant data
#         'all' - plots all data
#     nModeSensor : str
#         sensors to be used to calculate the modes
#         'FB' - feedback sensors
#         'TA' - toroidal array sensors
#     method : str
#         method to calculate mode analysis
#         'leastSquares' - performs a matrix least squares analysis
#         
#     Attributes
#     ----------
#     shotno : int
#         data shot number
#     title : str
#         title to be placed on each plot
#     nModeSensor : str
#         sensors to be used to calculate the modes
#     n1Amp : numpy.ndarray
#         n=1 mode amplitude data
#     n2Amp : numpy.ndarray
#         n=2 mode amplitude data
#     n1Phase : numpy.ndarray
#         filtered n=1 mode phase data
#     n1PhaseRaw : numpy.ndarray
#         raw n=1 mode phase data
#     n1Freq : numpy.ndarray
#         filtered n=1 mode frequency data
#     n1FreqRaw : numpy.ndarray
#         raw n=1 mode frequency data
#         
#     Subfunctions
#     ------------
#     plot :
#         plots relevant plots
#     Bn1 : 
#         Generates a pretend B_{n=1} signal at the toroidal location, phi0
#         Not presently in use
#     plotOfAmps : 
#         returns a plot of n=1 and n=2 mode amplitudes
#     plotOfN1Phase
#         returns a plot of the n=1 mode phase
#     self.plotOfN1Freq
#         returns a plot of the n=1 mode frequency
#     self.plotOfN1Amp
#         returns a plot of the n=1 mode amplitude
#         
#     Notes
#     -----
#     The convolution filters used with the phase and frequency "mess up" the 
#     last tenth of a millisecond of data
#     
#     """    

#     def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False,
#                  nModeSensor='FB',method='leastSquares',phaseFilter='gaussian',
#                  phaseFilterTimeConstant=1.0/100e3):
#         
#         self.shotno=shotno
#         self.title = '%d.  %s_4 sensor array.  n mode analysis' % (shotno,nModeSensor)
#         self.nModeSensor=nModeSensor        
#         
#         # load data from requested sensor array
#         if nModeSensor=='TA':
#             ## load TA data
#             temp=taData(self.shotno,tStart,tStop+0.5e-3);  # asking for an extra half millisecond (see Notes above) 
#             b=temp.taPolData
#             time=temp.taPolTime
#             phi=temp.phi
# #            [n,m]=_np.shape(data)
#         elif nModeSensor=='FB' or nModeSensor=='FB_S4':
#             fb=fbData(shotno,tStart,tStop)
#             time=fb.dfData.index.to_numpy()
#             b=fb.dfData.iloc[:,fb.dfData.columns.str.contains('S4P')].to_numpy()
#             phi=fb.dfMeta[fb.dfMeta.index.str.contains('S4P')].Phi.to_numpy()
#         
#         if method=='leastSquares':
#             n=len(phi)
#             A=_np.zeros((n,5))
#             A[:,0]=_np.ones(n);
#             A[:,1]=_np.sin(-phi)
#             A[:,2]=_np.cos(-phi) 
#             A[:,3]=_np.sin(-2*phi)
#             A[:,4]=_np.cos(-2*phi)
#             Ainv=_np.linalg.pinv(A)
#             x=Ainv.dot(b.transpose())
#             dfResults=_pd.DataFrame(data=x.transpose(),index=time,columns=['n0','n1Sin','n1Cos','n2Sin','n2Cos'])
#             dfResults['X1']=1j*dfResults['n1Sin']+dfResults['n1Cos']
#             dfResults['X2']=1j*dfResults['n2Sin']+dfResults['n2Cos']
#             dfResults['n1Amp']=_np.sqrt(dfResults['n1Sin']**2+dfResults['n1Cos']**2)
#             dfResults['n1Phase']=_np.arctan2(dfResults['n1Sin'],dfResults['n1Cos'])
#             dfResults['n2Amp']=_np.sqrt(dfResults['n2Sin']**2+dfResults['n2Cos']**2)
#             dfResults['n2Phase']=_np.arctan2(dfResults['n2Sin'],dfResults['n2Cos'])
#             dfResults['n1PhaseFilt']=_process.gaussianFilter(    dfResults.index.to_numpy(),
#                                                                 _process.unwrapPhase(dfResults['n1Phase'].to_numpy()),
#                                                                 timeFWHM=0.5e-4,
#                                                                 plot=False,
#                                                                 filterType='low')
#             dfResults['n2PhaseFilt']=_process.gaussianFilter(    dfResults.index.to_numpy(),
#                                                                 _process.unwrapPhase(dfResults['n2Phase'].to_numpy()),
#                                                                 timeFWHM=0.5e-4,
#                                                                 plot=False,
#                                                                 filterType='low')
#             dfResults['n1Freq']=_np.gradient(dfResults['n1PhaseFilt'])/_np.gradient(dfResults.index.to_numpy())/(_np.pi*2)
#             dfResults['n2Freq']=_np.gradient(dfResults['n2PhaseFilt'])/_np.gradient(dfResults.index.to_numpy())/(_np.pi*2)
#             dfResults['n1PhaseFilt']=_process.wrapPhase(dfResults['n1PhaseFilt'].to_numpy())
#             dfResults['n2PhaseFilt']=_process.wrapPhase(dfResults['n2PhaseFilt'].to_numpy())
#             
#         self.dfResults=dfResults
#         
#         ## plot data
#         if plot==True:
#             self.plot()
#             
#         elif plot == 'all':
#             self.plotOfSlice(index=int(n/4)).plot();
#             self.plotOfSlice(index=int(n/2)).plot();
#             self.plot(True)
#             
#     def plot(self,plotAll=False):
#         """
#         plots
#         """
#         fig,p1=_plt.subplots(4,sharex=True,figsize=[6,4])
#         p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
#                     self.dfResults['n1Cos']*1e4,
#                     label=r'n=1 Cos.',
#                     linewidth=0.8)
#         p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
#                     self.dfResults['n1Sin']*1e4,
#                     label=r'n=1 Sin.',
#                     linewidth=0.8)
#         p1[1].plot(    self.dfResults.index.to_numpy()*1e3,
#                     self.dfResults['n1Amp']*1e4,
#                     label=r'$\left| \delta B\right|_{n=1}$')
#         p1[2].plot(    self.dfResults.index.to_numpy()*1e3,
#                     self.dfResults['n1Phase'],
#                     label=r'$\angle \delta B_{n=1}$',
#                     marker='.',
#                     markersize=2,
#                     linestyle='')
# #        p1[2].plot(    self.dfResults.index.to_numpy()*1e3,
# #                    self.dfResults['n1PhaseFilt'],
# #                    label=r'$\angle \delta B_{n=1}$',
# #                    marker='.',
# #                    markersize=2,
# #                    linestyle='')
#         p1[3].plot(    self.dfResults.index.to_numpy()*1e3,
#                     self.dfResults['n1Freq']*1e-3,
#                     label='Frequency (kHz)')
#         if plotAll==True:
#             p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
#                         self.dfResults['n2Cos']*1e4,
#                         label=r'n=2 Cos.')
#             p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
#                         self.dfResults['n2Sin']*1e4,
#                         label=r'n=2 Sin.')
#             p1[1].plot(    self.dfResults.index.to_numpy()*1e3,
#                         self.dfResults['n2Amp']*1e4,
#                         label=r'$\left| \delta B\right|_{n=2}$')
#             p1[2].plot(    self.dfResults.index.to_numpy()*1e3,
#                         self.dfResults['n2Phase'],
#                         label=r'$\angle \delta B_{n=2}$',
#                         marker='.',
#                         markersize=2,
#                         linestyle='')
#             p1[3].plot(    self.dfResults.index.to_numpy()*1e3,
#                         self.dfResults['n2Freq']*1e-3,
#                         label='Frequency (kHz)')
#         _plot.finalizeSubplot(    p1[0],
#                                 ylabel='G',
#                                 ylim=[-12,12],
#                                 legendLoc='upper left',
#                                 fontSizeStandard=8,
#                                 fontSizeTitle=8,
#                                 title=self.title)
#         _plot.finalizeSubplot(    p1[1],
#                                 ylabel='G',
#                                 ylim=[-0.5,12],
#                                 legendLoc='upper left',
#                                 fontSizeStandard=8,
#                                 fontSizeTitle=8)
#         _plot.finalizeSubplot(    p1[2],
#                                 ylabel='rad.',
#                                 ylim=[-_np.pi,_np.pi],
#                                 legendLoc='upper left',
#                                 fontSizeStandard=8,
#                                 fontSizeTitle=8,
# #                                yticks=_np.array([-_np.pi,-_np.pi/2,0,_np.pi/2,_np.pi]),
#                                 )
#         _plot.finalizeSubplot(    p1[3],
#                                 ylabel='kHz',
#                                 xlabel='Time (ms)',
#                                 ylim=[-5,25],
#                                 legendLoc='upper left',
#                                 fontSizeStandard=8,
#                                 fontSizeTitle=8)
#         _plot.finalizeFigure(fig)
#         
#         
#     def plotOfSlice(self,index=0):
#         """
#         Plots fit data for a single time value
#         """
#         j=index;
#         [n,m]=_np.shape(self._data)
#         y=_np.zeros(n);
#         for i in range(0,n):
#                 y[i]=self._data[i][j]*1e4
#         p1=_plot.plot(shotno=[self.shotno],
#                       title=self.title+', t='+str(self.time[j]*1000)+'ms.')
#         phi=_np.linspace(0,_np.pi*2,100)
#         n1Fit=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)
#         n2Fit=self._x[0,j]+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)
#         fitTotal=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)

#         # plot
#         p1.addTrace(yData=y,xData=self._phi,
#                     marker='x',linestyle='',yLegendLabel='raw') 
#         p1.addTrace(yData=n1Fit,xData=phi,
#                     yLegendLabel='n=1') 
#         p1.addTrace(yData=n2Fit,xData=phi,
#                     yLegendLabel='n=2') 
#         p1.addTrace(yData=fitTotal,xData=phi,
#                     yLegendLabel='Superposition') 
#         return p1
#     
    
@_prepShotno
def nModeData_df(shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False,
                 nModeSensor='FB',method='leastSquares',phaseFilter='gaussian',
                 phaseFilterTimeConstant=0.5e-4):
    """
    This function performs n-mode (toroidal) mode analysis on the plasma.
    Provides mode amplitude, phase, and frequency
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        True - plots relevant data
        'all' - plots all data
    nModeSensor : str
        sensors to be used to calculate the modes
        'FB' - feedback sensors
        'TA' - toroidal array sensors
    method : str
        method to calculate mode analysis
        'leastSquares' - performs a matrix least squares analysis
        
    Attributes
    ----------
    shotno : int
        data shot number
    title : str
        title to be placed on each plot
    nModeSensor : str
        sensors to be used to calculate the modes
    n1Amp : numpy.ndarray
        n=1 mode amplitude data
    n2Amp : numpy.ndarray
        n=2 mode amplitude data
    n1Phase : numpy.ndarray
        filtered n=1 mode phase data
    n1PhaseRaw : numpy.ndarray
        raw n=1 mode phase data
    n1Freq : numpy.ndarray
        filtered n=1 mode frequency data
    n1FreqRaw : numpy.ndarray
        raw n=1 mode frequency data
        
    Subfunctions
    ------------
    plot :
        plots relevant plots
    Bn1 : 
        Generates a pretend B_{n=1} signal at the toroidal location, phi0
        Not presently in use
    plotOfAmps : 
        returns a plot of n=1 and n=2 mode amplitudes
    plotOfN1Phase
        returns a plot of the n=1 mode phase
    self.plotOfN1Freq
        returns a plot of the n=1 mode frequency
    self.plotOfN1Amp
        returns a plot of the n=1 mode amplitude
        
    Notes
    -----
    The convolution filters used with the phase and frequency "mess up" the 
    last tenth of a millisecond of data
    
    """    
    
    shotno=shotno
    title = '%d.  %s_4 sensor array.  n mode analysis' % (shotno,nModeSensor)
    nModeSensor=nModeSensor        
    
    # load data from requested sensor array
    if nModeSensor=='TA':
        ## load TA data
        temp=taData(shotno,tStart,tStop);  # asking for an extra half millisecond (see Notes above) 
        b=_np.array(temp.taPolData)
        time=temp.taPolTime
        phi=temp.phi
#            [n,m]=_np.shape(data)
    elif nModeSensor=='FB' or nModeSensor=='FB_S4':
        fb=fbData(shotno,tStart,tStop)
        time=fb.dfData.index.to_numpy()
        b=fb.dfData.iloc[:,fb.dfData.columns.str.contains('S4P')].to_numpy()
        b=b.transpose()
        phi=fb.dfMeta[fb.dfMeta.index.str.contains('S4P')].Phi.to_numpy()
    
    if method=='leastSquares':
        n=len(phi)
        A=_np.zeros((n,5))
        A[:,0]=_np.ones(n);
        A[:,1]=_np.sin(-phi)
        A[:,2]=_np.cos(-phi) 
        A[:,3]=_np.sin(-2*phi)
        A[:,4]=_np.cos(-2*phi)
        Ainv=_np.linalg.pinv(A)
        x=Ainv.dot(b)
        dfResults=_pd.DataFrame(data=x.transpose(),index=time,columns=['n0','n1Sin','n1Cos','n2Sin','n2Cos'])
        dfResults['X1']=1j*dfResults['n1Sin']+dfResults['n1Cos']
        dfResults['X2']=1j*dfResults['n2Sin']+dfResults['n2Cos']
        dfResults['n1Amp']=_np.sqrt(dfResults['n1Sin']**2+dfResults['n1Cos']**2)
        dfResults['n1Phase']=_np.arctan2(dfResults['n1Sin'],dfResults['n1Cos'])
        dfResults['n2Amp']=_np.sqrt(dfResults['n2Sin']**2+dfResults['n2Cos']**2)
        dfResults['n2Phase']=_np.arctan2(dfResults['n2Sin'],dfResults['n2Cos'])
        dfResults['n1PhaseFilt']=_process.gaussianFilter(    dfResults.index.to_numpy(),
                                                            _process.unwrapPhase(dfResults['n1Phase'].to_numpy()),
                                                            timeFWHM=phaseFilterTimeConstant,
                                                            plot=False,
                                                            filterType='low')
        dfResults['n2PhaseFilt']=_process.gaussianFilter(    dfResults.index.to_numpy(),
                                                            _process.unwrapPhase(dfResults['n2Phase'].to_numpy()),
                                                            timeFWHM=phaseFilterTimeConstant,
                                                            plot=False,
                                                            filterType='low')
        dfResults['n1Freq']=_np.gradient(dfResults['n1PhaseFilt'])/_np.gradient(dfResults.index.to_numpy())/(_np.pi*2)
        dfResults['n2Freq']=_np.gradient(dfResults['n2PhaseFilt'])/_np.gradient(dfResults.index.to_numpy())/(_np.pi*2)
        dfResults['n1PhaseFilt']=_process.wrapPhase(dfResults['n1PhaseFilt'].to_numpy())
        dfResults['n2PhaseFilt']=_process.wrapPhase(dfResults['n2PhaseFilt'].to_numpy())
        
    dfResults=dfResults
    
            
    def Plot(plotAll=False):
        """
        plots
        """
        fig,p1=_plt.subplots(4,sharex=True,figsize=[6,4])
        p1[0].plot(    dfResults.index.to_numpy()*1e3,
                    dfResults['n1Cos']*1e4,
                    label=r'n=1 Cos.',
                    linewidth=0.8)
        p1[0].plot(    dfResults.index.to_numpy()*1e3,
                    dfResults['n1Sin']*1e4,
                    label=r'n=1 Sin.',
                    linewidth=0.8)
        p1[1].plot(    dfResults.index.to_numpy()*1e3,
                    dfResults['n1Amp']*1e4,
                    label=r'$\left| \delta B\right|_{n=1}$')
        p1[2].plot(    dfResults.index.to_numpy()*1e3,
                    dfResults['n1Phase'],
                    label=r'$\angle \delta B_{n=1}$',
                    marker='.',
                    markersize=2,
                    linestyle='')
    #        p1[2].plot(    dfResults.index.to_numpy()*1e3,
    #                    dfResults['n1PhaseFilt'],
    #                    label=r'$\angle \delta B_{n=1}$',
    #                    marker='.',
    #                    markersize=2,
    #                    linestyle='')
        p1[3].plot(    dfResults.index.to_numpy()*1e3,
                    dfResults['n1Freq']*1e-3,
                    label='Frequency\n(kHz)')
        if plotAll==True:
            p1[0].plot(    dfResults.index.to_numpy()*1e3,
                        dfResults['n2Cos']*1e4,
                        label=r'n=2 Cos.')
            p1[0].plot(    dfResults.index.to_numpy()*1e3,
                        dfResults['n2Sin']*1e4,
                        label=r'n=2 Sin.')
            p1[1].plot(    dfResults.index.to_numpy()*1e3,
                        dfResults['n2Amp']*1e4,
                        label=r'$\left| \delta B\right|_{n=2}$')
            p1[2].plot(    dfResults.index.to_numpy()*1e3,
                        dfResults['n2Phase'],
                        label=r'$\angle \delta B_{n=2}$',
                        marker='.',
                        markersize=2,
                        linestyle='')
            p1[3].plot(    dfResults.index.to_numpy()*1e3,
                        dfResults['n2Freq']*1e-3,
                        label='Frequency (kHz)')
        _plot.finalizeSubplot(    p1[0],
                                ylabel='G',
                                ylim=[-12,12],
                                legendLoc='upper left',
                                fontSizeStandard=8,
                                fontSizeTitle=8,
                                title=title)
        _plot.finalizeSubplot(    p1[1],
                                ylabel='G',
                                ylim=[-0.5,12],
                                legendLoc='upper left',
                                fontSizeStandard=8,
                                fontSizeTitle=8)
        _plot.finalizeSubplot(    p1[2],
                                ylabel='rad.',
                                ylim=[-_np.pi,_np.pi],
                                legendLoc='upper left',
                                fontSizeStandard=8,
                                fontSizeTitle=8,
    #                                yticks=_np.array([-_np.pi,-_np.pi/2,0,_np.pi/2,_np.pi]),
                                )
        _plot.finalizeSubplot(    p1[3],
                                ylabel='kHz',
                                xlabel='Time (ms)',
                                ylim=[-5,25],
                                legendLoc='upper left',
                                fontSizeStandard=8,
                                fontSizeTitle=8)
        _plot.finalizeFigure(fig)
        
        print('saf')
            
        
    ## plot data
    if plot==True:
        Plot()
        
    return dfResults


def loadAllMagData(shotno):
    """
    Downloads all TA, FB, and PA data into a single dataframe for data and another
    for meta data
    
    #TODO
    -----
    taData_df, fbData_df, paData_df need to operating the same way.
    """
    
    taData,_,taMeta=taData_df(shotno)
    fbData,fbMeta=fbData_df(shotno)
    fbData=fbData.iloc[:,_np.invert(fbData.columns.str.contains('RAW'))]
    paData,_,paMeta=paData_df(shotno)
    
    dfDataAll=_pd.concat((taData,paData,fbData),axis=1)
    dfMetaAll=_pd.concat((taMeta,paMeta,fbMeta))
    
    for i,(key,val) in enumerate(dfMetaAll.iterrows()):
        if 'FB' in val.name:
            dfMetaAll.at[key,'marker']=u'^'
        elif 'TA' in val.name:
            dfMetaAll.at[key,'marker']=u'v'
        elif 'PA1' in val.name:
            dfMetaAll.at[key,'marker']=u'>'
        elif 'PA2' in val.name:
            dfMetaAll.at[key,'marker']=u'<'
            
    return dfDataAll,dfMetaAll



def loadAllRawMagData(shotno):
    """
    Downloads all TA, FB, and PA data into a single dataframe for data and another
    for meta data
    
    #TODO
    -----
    taData_df, fbData_df, paData_df need to operating the same way.
    """
    
    taData,taDataRaw,taMeta=taData_df(shotno)
    fbDataAll,fbMeta=fbData_df(shotno)
    fbData=fbDataAll.iloc[:,_np.invert(fbDataAll.columns.str.contains('RAW'))]
    fbDataRaw=fbDataAll.iloc[:,fbDataAll.columns.str.contains('RAW')]
    fbDataRaw.columns =fbData.columns.to_numpy().tolist()
    paData,paDataRaw,paMeta=paData_df(shotno)
    
    dfDataAll=_pd.concat((taDataRaw,paDataRaw,fbDataRaw),axis=1)
    dfMetaAll=_pd.concat((taMeta,paMeta,fbMeta))
    
    for i,(key,val) in enumerate(dfMetaAll.iterrows()):
        if 'FB' in val.name:
            dfMetaAll.at[key,'marker']=u'^'
        elif 'TA' in val.name:
            dfMetaAll.at[key,'marker']=u'v'
        elif 'PA1' in val.name:
            dfMetaAll.at[key,'marker']=u'>'
        elif 'PA2' in val.name:
            dfMetaAll.at[key,'marker']=u'<'
            
    return dfDataAll,dfMetaAll




@_prepShotno
class mModeData_df:
    """
    This function performs m-mode (poloidal) mode analysis on the plasma.
    Provides mode amplitude, phase, and frequency
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        True - plots relevant data
        'all' - plots all data
    mModeSensor : str
        sensors to be used to calculate the modes
        'PA1' or 'PA2' 
    method : str
        method to calculate mode analysis
        'leastSquares' - performs a matrix least squares analysis
        
    Attributes
    ----------
    shotno : int
        data shot number
    title : str
        title to be placed on each plot
    TODO add more
        
    Subfunctions
    ------------
    plot :
        plots relevant plots
    TODO add more
        
    Notes
    -----
    
    
    """    

    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False,
                 mModeSensor='PA1',method='leastSquares',phaseFilter='gaussian',
                 Lambda=0.0):
        
        self.shotno=shotno
        self.title = '%d.  %s sensor array.  n mode analysis' % (shotno,mModeSensor)
#        self.nModeSensor=nModeSensor    
        
        data=paData(shotno,tStart=tStart,tStop=tStop)
        dfPA=data.dfData
        dfMeta=data.dfMeta
        
        dfPA=dfPA.iloc[:,dfPA.columns.str.contains(mModeSensor)]
        dfMeta=dfMeta.iloc[dfMeta.index.str.contains(mModeSensor)]
        theta=dfMeta['Theta'].to_numpy()
#        Lambda=0.1
        
        def thetaStarCalc(theta, L=0):
            """
            theta correction in Jeff's thesis.  Eq. 4.4.
            L=0 does not alter theta
            """
            return theta-L*_np.sin(theta)
        
        theta=thetaStarCalc(theta,Lambda)
        
        b=dfPA.to_numpy()
        self._data=b
        
        time=dfPA.index.to_numpy()
        
        if method=='leastSquares':
            n=len(theta)
            A=_np.zeros((n,11))
            A[:,0]=_np.ones(n);
            A[:,1]=_np.sin(1*theta)
            A[:,2]=_np.cos(1*theta) 
            A[:,3]=_np.sin(2*theta)
            A[:,4]=_np.cos(2*theta) 
            A[:,5]=_np.sin(3*theta)
            A[:,6]=_np.cos(3*theta) 
            A[:,7]=_np.sin(4*theta)
            A[:,8]=_np.cos(4*theta) 
            A[:,9]=_np.sin(5*theta)
            A[:,10]=_np.cos(5*theta) 
            Ainv=_np.linalg.pinv(A)
            x=Ainv.dot(b.transpose())
            self._x=x
            dfResults=_pd.DataFrame(data=x.transpose(),index=time,columns=['n0','m1Sin','m1Cos','m2Sin','m2Cos','m3Sin','m3Cos','m4Sin','m4Cos','m5Sin','m5Cos'])
            dfResults['m1Amp']=_np.sqrt(dfResults['m1Sin']**2+dfResults['m1Cos']**2)
            dfResults['m1Phase']=_np.arctan2(dfResults['m1Sin'],dfResults['m1Cos'])
            dfResults['m2Amp']=_np.sqrt(dfResults['m2Sin']**2+dfResults['m2Cos']**2)
            dfResults['m2Phase']=_np.arctan2(dfResults['m2Sin'],dfResults['m2Cos'])
            dfResults['m3Amp']=_np.sqrt(dfResults['m3Sin']**2+dfResults['m3Cos']**2)
            dfResults['m3Phase']=_np.arctan2(dfResults['m3Sin'],dfResults['m3Cos'])
            dfResults['m4Amp']=_np.sqrt(dfResults['m4Sin']**2+dfResults['m4Cos']**2)
            dfResults['m4Phase']=_np.arctan2(dfResults['m4Sin'],dfResults['m4Cos'])
            dfResults['m5Amp']=_np.sqrt(dfResults['m5Sin']**2+dfResults['m5Cos']**2)
            dfResults['m5Phase']=_np.arctan2(dfResults['m5Sin'],dfResults['m5Cos'])
            dfResults['m1PhaseFilt']=_process.gaussianFilter(    dfResults.index.to_numpy(),
                                                                _process.unwrapPhase(dfResults['m1Phase'].to_numpy()),
                                                                timeFWHM=0.5e-4,
                                                                plot=False,
                                                                filterType='low')
            dfResults['m2PhaseFilt']=_process.gaussianFilter(    dfResults.index.to_numpy(),
                                                                _process.unwrapPhase(dfResults['m2Phase'].to_numpy()),
                                                                timeFWHM=0.5e-4,
                                                                plot=False,
                                                                filterType='low')
            dfResults['m3PhaseFilt']=_process.gaussianFilter(    dfResults.index.to_numpy(),
                                                                _process.unwrapPhase(dfResults['m3Phase'].to_numpy()),
                                                                timeFWHM=0.5e-4,
                                                                plot=False,
                                                                filterType='low')
            dfResults['m4PhaseFilt']=_process.gaussianFilter(    dfResults.index.to_numpy(),
                                                                _process.unwrapPhase(dfResults['m4Phase'].to_numpy()),
                                                                timeFWHM=0.5e-4,
                                                                plot=False,
                                                                filterType='low')
            dfResults['m5PhaseFilt']=_process.gaussianFilter(    dfResults.index.to_numpy(),
                                                                _process.unwrapPhase(dfResults['m5Phase'].to_numpy()),
                                                                timeFWHM=0.5e-4,
                                                                plot=False,
                                                                filterType='low')
            
            dfResults['m1Freq']=_np.gradient(dfResults['m1PhaseFilt'])/_np.gradient(dfResults.index.to_numpy())/(_np.pi*2)
            dfResults['m2Freq']=_np.gradient(dfResults['m2PhaseFilt'])/_np.gradient(dfResults.index.to_numpy())/(_np.pi*2)
            dfResults['m3Freq']=_np.gradient(dfResults['m3PhaseFilt'])/_np.gradient(dfResults.index.to_numpy())/(_np.pi*2)
            dfResults['m4Freq']=_np.gradient(dfResults['m4PhaseFilt'])/_np.gradient(dfResults.index.to_numpy())/(_np.pi*2)
            dfResults['m5Freq']=_np.gradient(dfResults['m5PhaseFilt'])/_np.gradient(dfResults.index.to_numpy())/(_np.pi*2)
            dfResults['m1PhaseFilt']=_process.wrapPhase(dfResults['m1PhaseFilt'].to_numpy())
            dfResults['m2PhaseFilt']=_process.wrapPhase(dfResults['m2PhaseFilt'].to_numpy())
            dfResults['m3PhaseFilt']=_process.wrapPhase(dfResults['m3PhaseFilt'].to_numpy())
            dfResults['m4PhaseFilt']=_process.wrapPhase(dfResults['m4PhaseFilt'].to_numpy())
            dfResults['m5PhaseFilt']=_process.wrapPhase(dfResults['m5PhaseFilt'].to_numpy())
            
        self.dfResults=dfResults
        
        ## plot data
        if plot==True:
            self.plot()
            
#        elif plot == 'all':
#            self.plotOfSlice(index=int(m/4)).plot();
#            self.plotOfSlice(index=int(m/2)).plot();
#            self.plot(True)
#            
    def plot(self,plotAll=False):
        """
        plots
        """
        fig,p1=_plt.subplots(4,sharex=True,figsize=[6,4])
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m1Cos']*1e4,
                    label=r'm=1 Cos.',
                    linewidth=0.8)
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m1Sin']*1e4,
                    label=r'm=1 Sin.',
                    linewidth=0.8)
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m2Cos']*1e4,
                    label=r'm=2 Cos.',
                    linewidth=0.8)
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m2Sin']*1e4,
                    label=r'm=2 Sin.',
                    linewidth=0.8)
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m3Cos']*1e4,
                    label=r'm=3 Cos.',
                    linewidth=0.8)
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m3Sin']*1e4,
                    label=r'm=3 Sin.',
                    linewidth=0.8)
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m4Cos']*1e4,
                    label=r'm=4 Cos.',
                    linewidth=0.8)
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m4Sin']*1e4,
                    label=r'm=4 Sin.',
                    linewidth=0.8)
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m5Cos']*1e4,
                    label=r'm=5 Cos.',
                    linewidth=0.8)
        p1[0].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m5Sin']*1e4,
                    label=r'm=5 Sin.',
                    linewidth=0.8)
        p1[1].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m1Amp']*1e4,
                    label=r'$\left| \delta B\right|_{m=1}$')
        p1[1].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m2Amp']*1e4,
                    label=r'$\left| \delta B\right|_{m=2}$')
        p1[1].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m3Amp']*1e4,
                    label=r'$\left| \delta B\right|_{m=3}$')
        p1[1].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m4Amp']*1e4,
                    label=r'$\left| \delta B\right|_{m=4}$')
        p1[1].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m5Amp']*1e4,
                    label=r'$\left| \delta B\right|_{m=5}$')
        p1[2].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m1Phase'],
                    label=r'$\angle \delta B_{m=1}$',
                    marker='.',
                    markersize=2,
                    linestyle='')
        p1[2].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m2Phase'],
                    label=r'$\angle \delta B_{m=2}$',
                    marker='.',
                    markersize=2,
                    linestyle='')
        p1[2].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m3Phase'],
                    label=r'$\angle \delta B_{m=3}$',
                    marker='.',
                    markersize=2,
                    linestyle='')
        p1[2].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m4Phase'],
                    label=r'$\angle \delta B_{m=4}$',
                    marker='.',
                    markersize=2,
                    linestyle='')
        p1[2].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m5Phase'],
                    label=r'$\angle \delta B_{m=5}$',
                    marker='.',
                    markersize=2,
                    linestyle='')
#        p1[2].plot(    self.dfResults.index.to_numpy()*1e3,
#                    self.dfResults['n1PhaseFilt'],
#                    label=r'$\angle \delta B_{n=1}$',
#                    marker='.',
#                    markersize=2,
#                    linestyle='')
        p1[3].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m1Freq']*1e-3,
                    label='Frequency')
        p1[3].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m2Freq']*1e-3,
                    label='Frequency')
        p1[3].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m3Freq']*1e-3,
                    label='Frequency')
        p1[3].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m4Freq']*1e-3,
                    label='Frequency')
        p1[3].plot(    self.dfResults.index.to_numpy()*1e3,
                    self.dfResults['m5Freq']*1e-3,
                    label='Frequency')
        _plot.finalizeSubplot(    p1[0],
                                ylabel='G',
                                ylim=[-12,12],
                                legendLoc='upper left',
                                fontSizeStandard=8,
                                fontSizeTitle=8,
                                title=self.title,
                                ncol=2)
        _plot.finalizeSubplot(    p1[1],
                                ylabel='G',
                                ylim=[-0.5,12],
                                legendLoc='upper left',
                                fontSizeStandard=8,
                                fontSizeTitle=8,
                                ncol=2)
        _plot.finalizeSubplot(    p1[2],
                                ylabel='rad.',
                                ylim=[-_np.pi,_np.pi],
                                legendLoc='upper left',
                                fontSizeStandard=8,
                                fontSizeTitle=8,
                                ncol=2
#                                yticks=_np.array([-_np.pi,-_np.pi/2,0,_np.pi/2,_np.pi]),
                                )
        _plot.finalizeSubplot(    p1[3],
                                ylabel='kHz',
                                xlabel='Time (ms)',
                                ylim=[-5,25],
                                legendLoc='upper left',
                                fontSizeStandard=8,
                                fontSizeTitle=8,
                                ncol=2)
        _plot.finalizeFigure(fig)
        
        
    def plotOfSlice(self,index=0):
        """
        Plots fit data for a single time value
        """
        j=index;
        [n,m]=_np.shape(self._data)
        y=_np.zeros(n);
        for i in range(0,n):
                y[i]=self._data[i][j]*1e4
        p1=_plot.plot(shotno=[self.shotno],
                      title=self.title+', t='+str(self.time[j]*1000)+'ms.')
        phi=_np.linspace(0,_np.pi*2,100)
        n1Fit=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)
        n2Fit=self._x[0,j]+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)
        fitTotal=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)

        # plot
        p1.addTrace(yData=y,xData=self._phi,
                    marker='x',linestyle='',yLegendLabel='raw') 
        p1.addTrace(yData=n1Fit,xData=phi,
                    yLegendLabel='n=1') 
        p1.addTrace(yData=n2Fit,xData=phi,
                    yLegendLabel='n=2') 
        p1.addTrace(yData=fitTotal,xData=phi,
                    yLegendLabel='Superposition') 
        return p1


@_prepShotno
class mModeData:
    """
    This function performs a least squares fit to a poloidal array of sensors and analyzes m=2,3 and 4 modes.  Mode amplitude, phase, and phase velocity. 
    In addtion, this code generates a perturbed B_pol(t) measurement as observed by a sensor at location, theta0
    Function uses either 32 poloidal PA1 or PA2 sensors
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    
    """  
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False,theta0=0,\
                 sensor='PA1',phaseFilter = 'gaussian',correctTheta=False):
        self.shotno=shotno
        self.title= '%d.  sensor = %s.  m mode analysis' % (shotno, sensor)
        
        data=paData(self.shotno,tStart=tStart,tStop=tStop,removeBadSensors=True,correctTheta=correctTheta);
        if sensor=='PA1':
            #self._data=data.pa1Data # Non pandas
            self._data=data.dfData.filter(regex="PA1_*").to_numpy().transpose()
            self.time=data.pa1Time
            self._theta=data.thetaPA1
            self._phi=data.phiPA1 # Account for torroidal location, implicitly assume n = 1 dominant structure
            [n,m]=_np.shape(self._data)
        if sensor=='PA2':
            #self._data=data.pa2Data
            self._data=data.dfData.filter(regex="PA2_*").to_numpy().transpose()
            self.time=data.pa2Time
            self._theta=data.thetaPA2
            self._phi=data.phiPA2
            [n,m]=_np.shape(self._data)

        ## Construct A matrix and its inversion
        A=_np.zeros((n,11))
        A[:,0]=_np.ones(n);
        A[:,1]=_np.sin(self._theta - self._phi)
        A[:,2]=_np.cos(self._theta - self._phi)
        A[:,3]=_np.sin(2*self._theta - self._phi)
        A[:,4]=_np.cos(2*self._theta - self._phi)
        A[:,5]=_np.sin(3*self._theta - self._phi)
        A[:,6]=_np.cos(3*self._theta - self._phi)
        A[:,7]=_np.sin(4*self._theta - self._phi)
        A[:,8]=_np.cos(4*self._theta - self._phi)
        A[:,9]=_np.sin(5*self._theta - self._phi)
        A[:,10]=_np.cos(5*self._theta - self._phi)
        Ainv=_np.linalg.pinv(A)
        self.A=A# Store for later use
        
        ## Solve for coefficients, x, for every time step and assign values to appropriate arrays 
        self._x=_np.zeros([11,m]);
        # self.m0Offset=_np.zeros(m)
        self.m0Amp=_np.zeros(m)
        self.m1Amp=_np.zeros(m)
        self.m1PhaseRaw=_np.zeros(m)
        self.m2Amp=_np.zeros(m)
        self.m2PhaseRaw=_np.zeros(m)   
        self.m3Amp=_np.zeros(m)
        self.m3PhaseRaw=_np.zeros(m)     
        self.m4Amp=_np.zeros(m)
        self.m4PhaseRaw=_np.zeros(m)     
        self.m5Amp=_np.zeros(m)
        self.m5PhaseRaw=_np.zeros(m)          
        for j in range(0,m):
            y=_np.zeros(n);
            for i in range(0,n):
                y[i]=self._data[i][j]*1e4
            self._x[:,j]=Ainv.dot(y)
            # self.m0Offset=self.x[0,j]
            self.m0Amp[j]=self._x[0,j]**2
            self.m1Amp[j]=_np.sqrt(self._x[1,j]**2+self._x[2,j]**2)
            self.m2Amp[j]=_np.sqrt(self._x[3,j]**2+self._x[4,j]**2)
            self.m3Amp[j]=_np.sqrt(self._x[5,j]**2+self._x[6,j]**2)
            self.m4Amp[j]=_np.sqrt(self._x[7,j]**2+self._x[8,j]**2)
            self.m5Amp[j]=_np.sqrt(self._x[9,j]**2+self._x[10,j]**2)
            self.m1PhaseRaw[j]=_np.arctan2(self._x[1,j],self._x[2,j])
            self.m2PhaseRaw[j]=_np.arctan2(self._x[3,j],self._x[4,j])
            self.m3PhaseRaw[j]=_np.arctan2(self._x[5,j],self._x[6,j])
            self.m4PhaseRaw[j]=_np.arctan2(self._x[7,j],self._x[8,j])
            self.m5PhaseRaw[j]=_np.arctan2(self._x[9,j],self._x[10,j])

        if phaseFilter == 'gaussian':
            self.m1Phase=_process.wrapPhase(
                            _process.gaussianLowPassFilter(
                                _process.unwrapPhase(self.m1PhaseRaw),
                                self.time,timeWidth=1./20e3))
            self.m2Phase=_process.wrapPhase(
                            _process.gaussianLowPassFilter(
                                _process.unwrapPhase(self.m2PhaseRaw),
                                self.time,timeWidth=1./20e3))
            self.m3Phase=_process.wrapPhase(
                            _process.gaussianLowPassFilter(
                                _process.unwrapPhase(self.m3PhaseRaw),
                                self.time,timeWidth=1./20e3))
            self.m4Phase=_process.wrapPhase(
                            _process.gaussianLowPassFilter(
                                _process.unwrapPhase(self.m4PhaseRaw),
                                self.time,timeWidth=1./20e3))
            self.m5Phase=_process.wrapPhase(
                            _process.gaussianLowPassFilter(
                                _process.unwrapPhase(self.m5PhaseRaw),
                                self.time,timeWidth=1./20e3))
            
        else:
            self.m1Phase=_np.zeros(len(self.m1PhaseRaw))  
            self.m1Phase[:]=self.m1PhaseRaw[:]
            self.m2Phase=_np.zeros(len(self.m2PhaseRaw))  
            self.m2Phase[:]=self.m2PhaseRaw[:]
            self.m3Phase=_np.zeros(len(self.m3PhaseRaw))  
            self.m3Phase[:]=self.m3PhaseRaw[:]
            self.m4Phase=_np.zeros(len(self.m4PhaseRaw))  
            self.m4Phase[:]=self.m4PhaseRaw[:]
            self.m5Phase=_np.zeros(len(self.m5PhaseRaw))  
            self.m5Phase[:]=self.m5PhaseRaw[:]
            

        self.m1Freq=_np.gradient(_process.unwrapPhase(self.m1Phase))/_np.gradient(self.time)/(2*_np.pi)        
        self.m2Freq=_np.gradient(_process.unwrapPhase(self.m2Phase))/_np.gradient(self.time)/(2*_np.pi)        
        self.m3Freq=_np.gradient(_process.unwrapPhase(self.m3Phase))/_np.gradient(self.time)/(2*_np.pi)        
        self.m4Freq=_np.gradient(_process.unwrapPhase(self.m4Phase))/_np.gradient(self.time)/(2*_np.pi)        
        self.m5Freq=_np.gradient(_process.unwrapPhase(self.m5Phase))/_np.gradient(self.time)/(2*_np.pi)        
        
            
        if plot == True:

            self.plot()
            
        elif plot == 'all':
            self.plotOfSlice(index=int(m/4)).plot();
            self.plotOfSlice(index=int(m/2)).plot();
            self.plotOfAmplitudes().plot()
            
        
    def plotOfAmplitudes(self):
        # plot amplitudes
        p1=_plot.plot(title=self.title,shotno=self.shotno,
                      xLabel='ms',yLabel='G',yLim=[0,20])
        p1.addTrace(yData=self.m1Amp,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=1}|$') 
        p1.addTrace(yData=self.m2Amp,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=2}|$') 
        p1.addTrace(yData=self.m3Amp,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=3}|$') 
        p1.addTrace(yData=self.m4Amp,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=4}|$') 
        p1.addTrace(yData=self.m5Amp,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=5}|$') 
        return p1
    
    def plotOfPhases(self):
        # plot amplitudes
        p1=_plot.plot(title=self.title,shotno=self.shotno,
                      xLabel='ms',yLabel='rad')
        p1.addTrace(yData=self.m1Phase,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=1}|$',linestyle='',marker='.')
        p1.addTrace(yData=self.m2Phase,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=2}|$',linestyle='',marker='.')
        p1.addTrace(yData=self.m3Phase,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=3}|$',linestyle='',marker='.')
        p1.addTrace(yData=self.m4Phase,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=4}|$',linestyle='',marker='.')
        p1.addTrace(yData=self.m5Phase,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=5}|$',linestyle='',marker='.')
        return p1
        
    def plotOfFreqs(self):
        # plot amplitudes
        p1=_plot.plot(title=self.title,shotno=self.shotno,
                      xLabel='ms',yLabel='kHz',yLim=[-20,20])
        p1.addTrace(yData=self.m1Freq*1e-3,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=1}|$')
        p1.addTrace(yData=self.m2Freq*1e-3,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=2}|$')
        p1.addTrace(yData=self.m3Freq*1e-3,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=3}|$')
        p1.addTrace(yData=self.m4Freq*1e-3,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=4}|$')
        p1.addTrace(yData=self.m5Freq*1e-3,xData=self.time*1000,
                   yLegendLabel=r'$|B_{pol, m=5}|$')
        return p1
    
    def plot(self):
        sp1=_plot.subPlot([self.plotOfAmplitudes(),self.plotOfPhases(),
                           self.plotOfFreqs()],plot=False)
        sp1.plot()
        return sp1
        
    # TODO:  add other plots
        
        
    def plotOfSlice(self,index=0):
        """
        Plot fits for a single instant in time
        """
        j=index;
        [n,m]=_np.shape(self._data)
        y=_np.zeros(n);
        for i in range(0,n):
            y[i]=self._data[i][j]*1e4
        p1=_plot.plot(title='t=%.3f ms. %s ' % (self.time[j]*1000, self.title),
                      shotno=self.shotno)
        theta=_np.linspace(self._theta[0],self._theta[-1],100)
#        m0Fit=self._x[0,j]
        m1Fit=self._x[0,j]+self._x[1,j]*_np.sin(theta)+self._x[2,j]*_np.cos(theta)
        m2Fit=self._x[0,j]+self._x[3,j]*_np.sin(2*theta)+self._x[4,j]*_np.cos(2*theta)
        m3Fit=self._x[0,j]+self._x[5,j]*_np.sin(3*theta)+self._x[6,j]*_np.cos(3*theta)
        m4Fit=self._x[0,j]+self._x[7,j]*_np.sin(4*theta)+self._x[8,j]*_np.cos(4*theta)
        m5Fit=self._x[0,j]+self._x[9,j]*_np.sin(5*theta)+self._x[10,j]*_np.cos(5*theta)
        fitTotal=(-4.)*self._x[0,j]+m1Fit+m2Fit+m3Fit+m4Fit+m5Fit  # the -4 corrects for the 4 extra offsets added from the preview 5 fits
        
        p1.addTrace(yData=y,xData=self._theta,
                    linestyle='',marker='.',yLegendLabel='raw')
        p1.addTrace(yData=m1Fit,xData=theta,
                    yLegendLabel='m=1')
        p1.addTrace(yData=m2Fit,xData=theta,
                    yLegendLabel='m=2')
        p1.addTrace(yData=m3Fit,xData=theta,
                    yLegendLabel='m=3')
        p1.addTrace(yData=m4Fit,xData=theta,
                    yLegendLabel='m=4')
        p1.addTrace(yData=m5Fit,xData=theta,
                    yLegendLabel='m=5')
        p1.addTrace(yData=fitTotal,xData=theta,
                    yLegendLabel='m=1-5')
        return p1
    
     # Build polar plot a la Jeff
    def plotPolar(self, tPoint = 1e-3):
        tPoint = _process.findNearest(self.time,tPoint) # time index
        data = _np.array(self._data)[:,tPoint]
        #data = _np.hstack( (data,data[0])) # wrap data
        #data = _np.hstack( (data,0)) # wrap data
        offset = _np.min(data)
        '''
        if self.correctTheta:
            self.thetaPA1=processPlasma.thetaCorrection(self.shotno,self.thetaPA1,\
                    self.tStart,self.tStop)[0]
            self.thetaPA1=self.thetaPA1[tPoint,:]
        '''
        # Build plot 
        _plt.figure()
        ax = _plt.subplot(111,projection='polar')
        # Plot Raw Data
        ax.plot(_np.hstack((self._theta,self._theta[0])),\
                _np.hstack( (data,data[0]))+1.5*offset,'*',ms=7,label='Raw')
        # Plot LSQ Fit
        theta = _np.hstack( (self._theta,self._theta[0]))
        r = _np.hstack( (self.A.dot(self._x[:,tPoint]),self.A.dot(self._x[:,tPoint])[0] ))*1e-4
        ax.plot(theta,r+1.5*offset,'-',lw=3,label='Sum')
        r = self.A[:,3:5].dot(self._x[3:5,tPoint])*1e-4
        r = _np.hstack((r,r[0]))
        ax.plot(theta,r+1.5*offset,'-',label='m=2: %2.1f'%\
                _np.sqrt(self._x[3,tPoint]**2+self._x[4,tPoint]**2))
        r = self.A[:,5:7].dot(self._x[5:7,tPoint])*1e-4
        r = _np.hstack((r,r[0]))        
        ax.plot(theta,r+1.5*offset,'-',label='m=3: %2.1f'%\
                _np.sqrt(self._x[5,tPoint]**2+self._x[7,tPoint]**2))
        r = self.A[:,7:9].dot(self._x[7:9,tPoint])*1e-4
        r = _np.hstack((r,r[0]))
        ax.plot(theta,r+1.5*offset,'-',label='m=4: %2.1f'%\
                _np.sqrt(self._x[7,tPoint]**2+self._x[8,tPoint]**2))
        r = self.A[:,9:11].dot(self._x[9:19,tPoint])*1e-4
        r = _np.hstack((r,r[0]))
        ax.plot(theta,r+1.5*offset,'-',label='m=4: %2.1f'%\
                _np.sqrt(self._x[7,tPoint]**2+self._x[10,tPoint]**2))
        
       
        _plt.legend()
        
        ax.set_rticks(_np.arange(-0.002,0,.0005))
        
        # Calculation of goodness of fit
        chiSq = _np.sum( (data-self.A.dot(self._x[:,tPoint])*1e-4)**2 /\
                        (self.A.dot(self._x[:,tPoint])*1e-4)**2 )
        ax.set_title("$\chi^2$ = %2.3e"%chiSq)
        '''
        tick = [1.05*ax.get_rmax(),ax.get_rmax()*0.97]
        for t  in _np.deg2rad(_np.arange(0,360,5)):
            ax.plot([t,t], tick, lw=0.72, color="k")

        '''
        _plt.show()
        return ax
        


@_prepShotno
class euvData:
    """""
    Pull EUV array data, similar to SXR format
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        default is False
        True - plots far array of all 11 (of 16) channels
        'all' - plots 
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to put at the top of figures
    data : list (of numpy.ndarray)
        list of 11 (of 16) data arrays, one for each channel
        
    Subfunctions
    ------------
    plotAll :
        plots all relevant plots     
    plotOfEUVStripey : 
        returns a stripey plot of the sensors
    plotOfOneChannel :
        returns a plot of a single channel based on the provided index, i
        
    Notes
    -----
    
    """""
    
    def __init__(self,shotno=101393,tStart=_TSTART,tStop=_TSTOP,plot=False):
        self.shotno=shotno
        self.title = '%d, EUV Array' % shotno
        
        self.detectors = [0,25,90,270]
        #".sensors.euv.pol.det000:channel_01:raw
        mdsAddressRaw = []
        mdsAddressR = []
        mdsAddressZ = []
        mdsAddressGain = []
        for det in self.detectors:
            for i in range(16):
                address = '\HBTEP2::TOP.SENSORS.EUV.POL.DET%03d' %det + '.CHANNEL_%02d' %(i+1)                
                mdsAddressRaw.append(address + ':RAW')
                mdsAddressR.append(address+':R')
                mdsAddressZ.append(address+':Z')
                mdsAddressGain.append(address+':GAIN')
        # Pull data
        self.data, self.time = mdsData(shotno,mdsAddressRaw,tStart,tStop)
        self.R = mdsData(shotno,mdsAddressR)
        self.Z = mdsData(shotno,mdsAddressZ)
        self.Gain = mdsData(shotno,mdsAddressGain)
        # aperture location
        self.det_ap_R = _np.array([1.160508, 1.090508, 0.907057, 0.929746]) 
        self.det_ap_Z = _np.array([0.000000, 0.099975, 0.166773,-0.173440])
        self.det_ap_Pol = _np.array([0.000635,.000635,.000635,.000635])
        self.det_ap_Tor = _np.array([0.000635,.0254,.0254,.0244])
        self.orient = _np.array([90,90,150,145]) # in degrees
        self.impactParameters=_np.array([-15.20239282, -13.88379581, -12.34664158, -10.57283874,\
        -8.55643016,  -6.31046679,  -3.87255431,  -1.30599562,\
         1.30599562,   3.87255431,   6.31046679,   8.55643016,\
        10.57283874,  12.34664158,  13.88379581,  15.20239282,\
         4.30322029,   5.05782374,   5.82121675,   6.58965689,\
         7.3590901 ,   8.12525048,   8.88376096,   9.63024278,\
        10.36044694,  11.07036625,  11.75634593,  12.41518118,\
        13.04417992,  13.64121203,  14.20472504,  14.73373985,\
         0.48160587,  -0.45909578,  -1.48806338,  -2.60489231,\
        -3.80409537,  -5.0735163 ,  -6.39312401,  -7.73487081,\
        -9.06416928, -10.34320912, -11.53564661, -12.61149041,\
       -13.55074303, -14.34481908, -14.99565891, -15.51324615,\
        14.72875124,  14.0599001 ,  13.28190271,  12.39369218,\
        11.39994292,  10.31173957,   9.14636266,   7.92605538,\
         6.67595335,   5.42160618,   4.18662977,   2.99095252,\
         1.84984336,   0.7737403 ,  -0.23139142,  -1.16322417])
        # plot 
        if plot==True:
            self.plotOfEUVStripey(tStart,tStop).plot()
        elif plot=='all':
            self.plotAll()
            self.plotOfEUVStripey(tStart,tStop).plot()
        
        
    def plotOfEUVStripey(self,tStart=1e-3,tStop=10e-3):
        iStart=_process.findNearest(self.time,tStart)
        iStop=_process.findNearest(self.time,tStop)
        p1=_plot.plot(title=self.title,subtitle='EUV Fan Array',
                      xLabel='Time [ms]', yLabel='Sensor Number',zLabel='a.u.',
                      plotType='contour',colorMap=_plot._red_green_colormap(),
                      centerColorMapAroundZero=True)
        #_plt.set_cmap('plasma')
        data=self.data;
        print(_np.shape(data))
        print(iStart)
        print(iStop)
        print(_np.shape(self.time))
        for i in range(0,len(data)):
            data[i]=data[i][iStart:iStop]
        p1.addTrace(self.time[iStart:iStop]*1e3,_np.arange(64),
                    _np.array(data))
        return p1
    
    
    def plotOfOneChannel(self, i=0):
        """ Plot one of the EUV chanenl.  based on index, i. """
        p1=_plot.plot(xLabel='time [ms]',yLabel=r'a.u.',title=self.title,
                      shotno=[self.shotno],subtitle='%d' %i);
        
        # smoothed data
        p1.addTrace(yData=self.data[i],xData=self.time*1000,
                    yLegendLabel=str(i))   
            
        return p1
    
    def plotAll(self):
        sp1=[]
        count=0
        for i in range(0,len(self.data)):
            newPlot=self.plotOfOneChannel(count)
            newPlot.subtitle=str(count)
            newPlot.yLegendLabel=[]
            newPlot.plot()
            sp1.append(newPlot)
 
            count+=1;
                
        sp1[0].title=self.title
        sp1=_plot.subPlot(sp1,plot=False)
        # sp1.shareY=True;
#        sp1.plot()
        
        return sp1
#@_prepShotno
#def _hbtPlot(shotnos=_np.array([98147, 98148]),plot=True,bp=False,tZoom=[2e-3,4e-3],saveFig=True):
#    """
#    This function acts similarly to hbtplot.py
#    Still under development
#    """
#    try:
#        len(shotnos)
#    except:
#        shotnos=_np.array([shotnos])
#        
#    for i in range(0,len(shotnos)):
#        
#        shotno=shotnos[i]
#        print(str(shotno))
#        
#        subplots=[]
#        
#        
#        ip=ipData(shotno)
#        subplots.append(ip.plotOfIP())
#        
#        q=qStarData(shotno)
#        subplots.append(q.plotOfQStar())
#        
#        rad=plasmaRadiusData(shotno)
#        subplots.append(rad.plotOfMajorRadius())
#        
#        pol=paData(shotno)
#        subplots.append(pol.plotOfPA1Stripey(tZoom[0],tZoom[1]))
#        
#        mode=nModeData(shotno)
#        subplots.append(mode.plotOfN1Amp())
#        n1freq=mode.plotOfN1Freq()
#        n1freq.yLim=[-10,20]
#        subplots.append(n1freq)
#        
#        if bp==True:
#            bpIn=bpData(shotno)
#            bpIn.bps9Current*=-1
#            bpIn.bps9Voltage*=-1
#            subplots.append(bpIn.plotOfBPS9Voltage())
#            subplots.append(bpIn.plotOfBPS9Current())
#            
#        _plot.subPlot(subplots)
#    
#        if saveFig==True:
#            p=_plt.gcf()
#            p.savefig(str(shotno)+'.png')
#            _plt.close('all')
#            


@_prepShotno
class xrayData:
    """
    Gets x-ray sensor data
    X-ray sensor is currently connected at: devices.north_rack:cpci:input_96 
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    xray : numpy.ndarray
        xray data
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plotOfXray : 
        returns the plot of x-ray intensity vs time
    plot :
        Plots all relevant plots
    
    """
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False,verbose=False):
        
        self.shotno = shotno
        self.title = "%d, X-ray Data" % shotno
        
        # get data
        
        if shotno>=103281:
            data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_12:INPUT_05 '],
                              tStart=tStart, tStop=tStop)
        else:
            data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.DEVICES.NORTH_RACK:CPCI:INPUT_96 '],
                              tStart=tStart, tStop=tStop)
        self.xray=data[0];
        self.time=time;
        
        if plot == True or plot=='all':
            self.plot()
        
        
    def plotOfXray(self):
        """
        returns the plot of xray data vs time
        """
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.xray)
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='Intensity')
        _plot.finalizeFigure(fig,title=self.title)
        
        return p1
        
            
    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfXray().plot()


class sxrMidplaneData:
    """
    Gets Soft X-ray midplane sensor data at: devices.north_rack:cpci:input_74 
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    sxr : numpy.ndarray
        sxr data
    time : numpy.ndarray
        time data
        
    Subfunctions
    ------------
    plotOfSXR : 
        returns the plot of sxr vs time
    plot :
        Plots all relevant plots
    
    """
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False,verbose=False):
        
        self.shotno = shotno
        self.title = "%d, SXR Midplane Data" % shotno
        
        # get data
        
        data, time=mdsData(shotno=shotno,
                           dataAddress=['\HBTEP2::TOP.DEVICES.NORTH_RACK:CPCI:INPUT_74 '],
                           tStart=tStart, tStop=tStop)

        self.sxr=-1*data[0];
        self.time=time;
        
        if plot == True or plot=='all':
            self.plot()
        
        
    def plotOfSXR(self):
        """
        returns the plot of xray data vs time
        """
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.sxr)
        _plot.finalizeSubplot(p1,xlabel='Time (ms)',ylabel='Intensity')
        _plot.finalizeFigure(fig,title=self.title)
        
        return p1
        
            
    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfSXR().plot()



class polBetaLi:
    """
    Calculate the sum beta_p + li/2 using Friedberg eq 6.90 (p 150)
    betap + li/2 = 4*pi*R*Bv/(mu0*Ip) + 3/2 - ln(8*R/a)
    Coefficients for calculating Bv are found by Jeff (see get_lambda.pro)
    """
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False,verbose=False):
        
        self.shotno = shotno
        self.title = "%d, Beta_p + li/2" % shotno
        
        # get data
        mu0 = 4E-7 * _np.pi 
        
        ip = ipData(shotno, tStart, tStop, False, findDisruption=False)
        Ip = ip.ip
        
        plasmaRadius = plasmaRadiusData(shotno, tStart, tStop, False)
        R = plasmaRadius.majorRadius
        a = plasmaRadius.minorRadius
        self.time = plasmaRadius.time
        
        bankData = capBankData(shotno, tStart, tStop, False)
        OH_cur = bankData.ohBankCurrent
        VF_cur = bankData.vfBankCurrent
        Bv = VF_cur * 2.6602839e-06 - OH_cur * 1.9580808e-08
        
        #Calculate sum
        self.polBetaLi = (4 * _np.pi * R * Bv)/(mu0 * Ip) + 1.5 - _np.log(8 * R / a)
        
        if plot == True or plot=='all':
            self.plot()
            
    def plotOfpolBetaLi(self):
        """
        returns the plot of xray data vs time
        """
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.polBetaLi)
        _plot.finalizeSubplot(p1,xlabel='Time (ms)', ylabel='sum')
        _plot.finalizeFigure(fig,title=self.title)
        _plt.ylim(0,3)
        #_plot.yLim(0,10)
        return p1

    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfpolBetaLi().plot()        
        
        return

class nGreenwald:
    """
    Calculate Greenwald's density limit
    ng = 1e14 * elong * IP(kA)*1000. / (!pi*minor_r^2) ; in m^-3
    """
    def __init__(self,shotno=96530,tStart=_TSTART,tStop=_TSTOP,plot=False,verbose=False):
        
        self.shotno = shotno
        self.title = "%d, Greenwald limit" % shotno
        # get data
        elong = 1
        ip = ipData(shotno, tStart, tStop, False, findDisruption=False)
        Ip = ip.ip
        self.time = ip.time
        plasmaRadius = plasmaRadiusData(shotno, tStart, tStop, False)
        a = plasmaRadius.minorRadius
        
        self.nGreenwald = (1e14*elong*Ip)/(_np.pi*a**2) 
        if plot == True or plot=='all':
            self.plot()
            
    def plotOfNGreenwald(self):
        """
        returns the plot of xray data vs time
        """
        fig,p1=_plt.subplots()
        p1.plot(self.time*1e3,self.nGreenwald)
        _plot.finalizeSubplot(p1,xlabel='Time (ms)', ylabel='m^-3')
        _plot.finalizeFigure(fig,title=self.title)
        _plt.ylim(0,5e19)
        #_plot.yLim(0,10)
        return p1

    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfNGreenwald().plot()        
        return

class metadata:
    """
    Retrive shot metadata
    """
    def __init__(self,shotno=96530, display=False):
        
        self.shotno = shotno
        try:
            self.comment = mdsData(shotno=shotno, dataAddress=['\HBTEP2::TOP.METADATA:COMMENT'])[0]
        except:
            self.comment = ''
        try:
            self.date = mdsData(shotno=shotno, dataAddress=['\HBTEP2::TOP.METADATA:DATE'])[0]
        except:
            self.date = ''
        try:
            self.operator = mdsData(shotno=shotno, dataAddress=['\HBTEP2::TOP.METADATA:OPERATOR'])[0]
        except:
            self.operator = ''
        try:    
            self.post_comment = mdsData(shotno=shotno, dataAddress=['\HBTEP2::TOP.METADATA:POST_COMMENT'])[0]
        except:
            self.post_comment = ''
            
        if display == True:
            self.show()
            
    def printMetadata(self):
        """
        Print all metadata in console
        """
        print('--COMMENT--')
        print(self.comment)
        print('--DATE--')
        print(self.date)
        print('--OPERATOR--')
        print(self.operator)
        print('--POST COMMENT--')
        print(self.post_comment)
        return 

    def show(self):
        """ 
        Print
        """
        self.printMetadata()    
        return
    
    def toDataFrame(self, include_shotno=False):
        '''
        Convert metadata to a 1-row pd.dataframe
        '''
        if include_shotno == False:
            df = _pd.DataFrame(columns=['comment', 'date', 'operator', 'post_comment'])
            df = df.append({'comment':self.comment, 'date':self.date, 
                            'operator':self.operator, 'post_comment':self.post_comment}, 
                           ignore_index=True)
        else:
            df = _pd.DataFrame(columns=['shotno', 'comment', 'date', 'operator', 'post_comment'])
            df = df.append({'shotno':self.shotno, 'comment':self.comment, 'date':self.date, 
                            'operator':self.operator, 'post_comment':self.post_comment}, 
                           ignore_index=True)               
        return df
###############################################################################
### debugging code

def _debugPlotExamplesOfAll():
    """ 
    This code executes each function above and plots an example of most every 
    plot.  Effectively, this allows the testing of most every function all in 
    one go.
    
    TODO(John):  This needs to be updated in light of _prepShotno().
    TODO(John):  Also, provide default shot numbers here.
    """
    bpData(plot=True)
    capBankData(plot=True)
    cos1RogowskiData(plot=True)
    quartzJumperData(plot=True)
    fbData(plot=True)
    ipData(plot=True)
    loopVoltageData(plot=True)
    mModeData(plot=True)
    nModeData(plot=True)
    paData(plot=True)
    plasmaRadiusData(plot=True)
    qStarData(plot=True)
    solData(plot=True)
    spectrometerData(plot=True)
    taData(plot=True)
    tfData(plot=True)
    tpData(plot=True)
    sxrData(plot=True)
    

