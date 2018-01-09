"""
_rwHBTData.py - load HBT data, both processed and unprocessed

NOTES
-----
The convention for units is to maintain everything in SI until they are
plotted.  

A few of these functions are merely wrappers for other people's code

In most of the functions below, I've specified default shotno's.  This is 
largely to make bebugging easier as there is nothing special about the 
provided shotnos.  
"""

###############################################################################
### import libraries

# common libraries 
import numpy as _np
import MDSplus as _mds
from copy import copy as _copy
#import pickle as _pk
import sys as _sys
import _socket
import os as _os

# hbtepLib libraries
import _rwDataTools as _rwData
import _processData as _process
import _plotTools as _plot
try:
    import _hbtPreferences as _pref
except ImportError:
    _sys.exit("Code hault: _hbtPreferences.py file not found.  See readme.md" +
    " concerning the creation of hbtPreferences.py")
        
        
###############################################################################
### constants
_REMOTE_DATA_PATH='/opt/hbt/data/control'  
if _socket.gethostname()==_pref._HBT_SERVER_NAME:
    _ON_HBTEP_SERVER=True;
else:
    _ON_HBTEP_SERVER=False;


###############################################################################
### global variables

# default time limits for data.  units in seconds
_TSTART = 0*1e-3;
_TSTOP  = 10*1e-3;

# list of known bad sensors.   outdated???  also note that i'm not YET doing anything with this info...
SENSORBLACKLIST = ['PA1_S29R', 'PA1_S16P', 'PA2_S13R', 'PA2_S14P', 'PA2_S27P', 'FB03_S4R', 'FB04_S4R', 'FB08_S3P', 'FB10_S1P', 'FB04_S3P', 'FB06_S2P', 'FB08_S4P', 'TA07_S1P', 'TA02_S1P', 'TA02_S2P'];

# directory where unprocessed or minimally processed data is written locally.   
# TODO:  move to hbtPreferences.py
_FILEDIR='/home/john/shotData/'



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
    it is assumed that tStart and tStop have the same units as time.  
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
    
    
def _initMDSConnection(shotno):
    """
    Initiate connection with MDSplus HBT-EP tree
    
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
        
        
def mdsData(shotno=None,
            dataAddress=['\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_95',
                         '\HBTEP2::TOP.DEVICES.SOUTH_RACK:CPCI_10:INPUT_96'],
            tStart=0*1e-3,tStop=10*1e-3):
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
    # if shotno is specified, this function gets its own mdsConn
    if type(shotno) is float or type(shotno) is int or type(shotno) is _np.int64:
        mdsConn=_initMDSConnection(shotno);
        
#    return mdsConn
    if type(dataAddress) is not list:
        dataAddress=[dataAddress];
    
    # get data and time
    data = [];
    for i in range(0,len(dataAddress)):
        data.append(mdsConn.get(dataAddress[i]).data())
    time = mdsConn.get('dim_of('+dataAddress[0]+')').data();  # time assocated with data
        
    # check to see if units are in seconds and NOT in milliseconds
    if tStop > 1:
        tStart=tStart*1e-3;
        tStop=tStop*1e-3;
        
    # trim time and data
    time,data= _trimTime(time,data,tStart,tStop)
        
    # return data and time
    return data, time
    
    
###############################################################################
### get device specific data
    
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
    def __init__(self,shotno=96530,tStart=0*1e-3,tStop=10*1e-3,plot=False):
        self.shotno = shotno
        self.title = "shotno = %d, Ip Data" % shotno
        
        # get data
        data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.ROGOWSKIS:IP'],
                              tStart=tStart, tStop=tStop)
        self.ip=data[0];
        self.time=time;
        
        if plot == True:
            self.plot()
        
        
    def plotOfIP(self):
        """
        returns the plot of IP vs time
        """
        p1=_plot.plot()
        p1.yLabel='kA'
        p1.xLabel='time [ms]'
        p1.subtitle='Plasma Current'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(xData=self.time*1000,yData=self.ip*1e-3)
        
        return p1
        
            
    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfIP().plot()
        
        
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
    def __init__(self,shotno=96530,tStart=0*1e-3,tStop=10*1e-3,plot=False):
        self.shotno = shotno
        self.title = "shotno = %d, Cos1 Rog. Data" % shotno
        
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
        
        if plot == True:
            self.plot()
        
        
    def plotOfCos1(self):
        """
        returns the plot of cos1 rog vs time
        """
        p1=_plot.plot()
        p1.yLabel=''
        p1.xLabel='time [ms]'
        p1.subtitle='Cos1 Rogowski'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(xData=self.time*1000,yData=self.cos1)
        
        return p1
        
            
    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfCos1().plot()
        
  
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
    BPS5 was moved to section 2 (now BPS2) summer of 2017.  
    
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
    def __init__(self,shotno=98147,tStart=0*1e-3,tStop=10*1e-3,plot=False):
        self.shotno = shotno
        self.title = "%s, BP Data." % shotno

        if shotno > 92000:
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
        
    def plotOfVoltage(self):
        """ 
        returns plot of BP voltage          
        """
        p1=_plot.plot(title=self.title,yLim=[-200,200],yLabel='V',
                      xLabel='Time [ms]',subtitle='BP Voltage',
                      shotno=[self.shotno]);
        p1.addTrace(xData=self.time*1000,yData=self.bps9Voltage,
                    yLegendLabel='BPS9')
        p1.addTrace(xData=self.time*1000,yData=self.bps5Voltage,
                    yLegendLabel='BPS5')
        
        return p1
        
    def plotOfCurrent(self):
        """ 
        returns plot of BP current         
        """
        p1=_plot.plot(title=self.title,yLabel='A',
                      xLabel='Time [ms]',subtitle='BP Current',
                      shotno=[self.shotno])
        p1.addTrace(xData=self.time*1000,yData=self.bps9Current,
                    yLegendLabel='BPS9')
        p1.addTrace(xData=self.time*1000,yData=self.bps5Current,
                    yLegendLabel='BPS5')
        
        return p1
            
    def plotOfBPS9Voltage(self):
        """ 
        returns plot of BPS9 voltage         
        """
        p1=_plot.plot(title=self.title,yLabel='A',yLim=[-200,200],
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
            sp1=_plot.subPlot([self.plotOfVoltage(),self.plotOfCurrent(),self.plotOfGPUVoltageRequest()])
        return sp1
        
    
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
    
    def __init__(self,shotno=95996,tStart=0*1e-3,tStop=10*1e-3,plot=False,probes='both'):  #sectionNum=2,
        
        self.shotno = shotno
        self.title = 'shotno = %s, triple probe data' % shotno
        self.probes=probes
        
        # enforce probes naming convetion
        if probes=='5':
            probes = 'tps5'
        if probes=='8':
            probes = 'tps8'

        # constants
        A=1.5904e-5 #(1.5mm)^2/4*pi + pi*(3.0mm)*(1.5mm), probe area
        e=1.602e-19; # fundamental charge
        eV=1.60218e-19; # 1 eV = 1.60218e-19 joules
        M=2.014102*1.66054e-27;  # approx 2 amu converted to kg
      
        ## Grab data
        if shotno > 90000: # Shotno after 2017 summer upgrade = 97239.  TPS2 was moved to section 5.  Now, it's TPS5.
            if probes=='both' or probes=='tps5' or probes=='tps2':
                
                # get data                
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
                self.tps5Density=self.tps5Current/(e*_np.sqrt(tps5Temp*eV/(M))*A);
    
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
                tps8Temp[tps5Temp<=0]=1e6; # i trim here to avoid imaginary numbers when I take the square root below
                self.tps8Density=self.tps8Current/(e*_np.sqrt(tps8Temp*eV/(M))*A);
                
        else:
            _sys.exit("Requested shot number range not supported yet.  Update code.")
        
        if plot==True:
            self.plot();
        elif plot=='all':
            self.plot(True)
        
    def plotOfKTe(self):
        p1=_plot.plot();
        p1.yLabel='eV'
        p1.subtitle='Electron Temperature'
        p1.title=self.title
        p1.yLim=[-50, 100]       
        if self.probes=='both' or self.probes=='tps5': 
            p1.addTrace(xData=self.tps5Time*1000,yData=self.tps5Temp,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8Temp,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        p1.shotno=[self.shotno]
        return p1
        
    def plotOfNe(self):
        p1=_plot.plot();
        p1.yLabel=r'$m^{-3}$ $10^{18}$'
        p1.subtitle='Density'
        p1.yLim=[-1, 4.5]
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5Density/1e18,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8Density/1e18,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        p1.shotno=[self.shotno]
        return p1
            
    def plotOfVf(self):
        p1=_plot.plot();
        p1.yLabel='V'
        p1.subtitle='Floating Potential'
        p1.xLabel='time [ms]'
        p1.yLim=[-150, 75]
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5VFloat,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8VFloat,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        p1.shotno=[self.shotno]
        return p1
            
    def plotOfTipA(self):
        # initialize tip A potential plot
        p1=_plot.plot();
        p1.yLabel='V'
        p1.subtitle='Tip A Potential'
        p1.xLabel='time [ms]'
#        self.plotOfTipA.yLim=[-100, 100]
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5TipA,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8TipA,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        p1.shotno=[self.shotno]
        return p1
            
    def plotOfTipB(self):
        # initialize tip B potential plot
        p1=_plot.plot();
        p1.yLabel='V'
        p1.subtitle='Tip B Potential'
        p1.xLabel='time [ms]'
#        self.plotOfTipA.yLim=[-100, 100]
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5TipB,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8TipB,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        p1.shotno=[self.shotno]
        return p1
            
    def plotOfTipC(self):            
        # initialize tip C potential plot
        p1=_plot.plot();
        p1.yLabel='V'
        p1.subtitle='Tip C Potential'
        p1.xLabel='time [ms]'
#        self.plotOfTipA.yLim=[-100, 100]
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5TipC,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8TipC,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')
        p1.shotno=[self.shotno]
        return p1
        
    def plotOfISat(self):            
        # initialize ion sat current
        p1=_plot.plot();
        p1.yLabel='A'
        p1.subtitle='Ion Sat. Current'
        p1.xLabel='time [ms]'
#        self.plotOfTipA.yLim=[-100, 100]
        if self.probes=='both' or self.probes=='tps5':
            p1.addTrace(yData=self.tps5Current,xData=self.tps5Time*1000,
                        yLegendLabel='TPS5')
        if self.probes=='both' or self.probes=='tps8':
            p1.addTrace(yData=self.tps8Current,xData=self.tps8Time*1000,
                        yLegendLabel='TPS8')   
        p1.shotno=[self.shotno]
        return p1
            
    def plot(self,plotAll=False):
        if plotAll == False:
            _plot.subPlot([self.plotOfKTe(),self.plotOfNe(),self.plotOfVf()]);
        else:
            _plot.subPlot([self.plotOfKTe(),self.plotOfNe(),self.plotOfVf()]);
            _plot.subPlot([self.plotOfTipA(),self.plotOfTipB(),self.plotOfTipC(),
                         self.plotOfISat()]);            
  
    
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

    """
    
    def __init__(self,shotno=95540,tStart=0*1e-3,tStop=10*1e-3,plot=False,
                 smoothingAlgorithm='tripleBoxCar'):
        self.shotno = shotno
        self.title = '%d, PA sensors' % shotno
        
        # poloidal location (in degrees)
        self.theta = _np.array([    5.625,      16.875,     28.125,     39.375,     50.625,     61.875,     73.125,     84.375,     95.625,     106.875,    118.125,    129.375,    140.625,    151.875,    163.125,    174.375,    185.625,    196.875,    208.125,    219.375,    230.625,    241.875,    253.125,    264.375,    275.625,    286.875,    298.125,    309.375,    320.625,    331.875,    343.125,    354.375])*_np.pi/180.
        
        # sensor names
        self.namesPA1=_np.array([   'PA1_S01P', 'PA1_S02P', 'PA1_S03P', 'PA1_S04P', 'PA1_S05P', 'PA1_S06P', 'PA1_S07P', 'PA1_S08P', 'PA1_S09P', 'PA1_S10P', 'PA1_S11P', 'PA1_S12P', 'PA1_S13P', 'PA1_S14P', 'PA1_S15P', 'PA1_S16P', 'PA1_S17P', 'PA1_S18P', 'PA1_S19P', 'PA1_S20P', 'PA1_S21P', 'PA1_S22P', 'PA1_S23P', 'PA1_S24P', 'PA1_S25P', 'PA1_S26P', 'PA1_S27P', 'PA1_S28P', 'PA1_S29P', 'PA1_S30P', 'PA1_S31P', 'PA1_S32P'])
        self.namesPA2=_np.array([   'PA2_S01P', 'PA2_S02P', 'PA2_S03P', 'PA2_S04P', 'PA2_S05P', 'PA2_S06P', 'PA2_S07P', 'PA2_S08P', 'PA2_S09P', 'PA2_S10P', 'PA2_S11P', 'PA2_S12P', 'PA2_S13P', 'PA2_S14P', 'PA2_S15P', 'PA2_S16P', 'PA2_S17P', 'PA2_S18P', 'PA2_S19P', 'PA2_S20P', 'PA2_S21P', 'PA2_S22P', 'PA2_S23P', 'PA2_S24P', 'PA2_S25P', 'PA2_S26P', 'PA2_S27P', 'PA2_S28P', 'PA2_S29P', 'PA2_S30P', 'PA2_S31P', 'PA2_S32P'])

        # compile full sensor addresses names
        pa1SensorAddresses=[]
        pa2SensorAddresses=[]        
        rootAddress='\HBTEP2::TOP.SENSORS.MAGNETIC:';
        for i in range(0,32):
            pa1SensorAddresses.append(rootAddress+self.namesPA1[i])
            pa2SensorAddresses.append(rootAddress+self.namesPA2[i])
            
        # get raw data
        self.pa1Raw,self.pa1Time=mdsData(shotno,pa1SensorAddresses, tStart, tStop)
        self.pa2Raw,self.pa2Time=mdsData(shotno,pa2SensorAddresses, tStart, tStop)
        
        # data smoothing algorithm
        self.pa1Data=[]
        self.pa1RawFit=[]
        self.pa2Data=[]
        self.pa2RawFit=[]
        if smoothingAlgorithm=='tripleBoxCar':
            # jeff's triple boxcar smoothing
            for i in range(0,32):
                temp=_process.convolutionSmoothing(self.pa1Raw[i][:],101,'box')
                temp=_process.convolutionSmoothing(temp,101,'box')
                temp=_process.convolutionSmoothing(temp,21,'box')
                self.pa1RawFit.append(temp)
                self.pa1Data.append(self.pa1Raw[i]-temp)
                
                temp=_process.convolutionSmoothing(self.pa2Raw[i][:],101,'box')
                temp=_process.convolutionSmoothing(temp,101,'box')
                temp=_process.convolutionSmoothing(temp,21,'box')
                self.pa2RawFit.append(temp)
                self.pa2Data.append(self.pa2Raw[i]-temp)
        else:
            _sys.exit("You must specify a correct smoothing algorithm.  Exiting code...")
        
        # plot 
        if plot==True or plot=='all':
            self.plot(True)
        if plot=='sample':
            self.plotOfPA1().plot();
            self.plotOfPA2().plot();
            

    def plotOfPA1(self, i=0, alsoPlotRawAndFit=True):
        """ Plot one of the PA1 plots.  based on index, i. """
        p1=_plot.plot();
        p1.xLabel='time [ms]'
        p1.yLabel=r'dB [G]'
        p1.title=str(self.shotno)+'. '+str(self.namesPA1[i])+' data'
        p1.shotno=[self.shotno]
        
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
        p1=_plot.plot();
        p1.xLabel='time [ms]'
        p1.yLabel=r'dB [G]'
        p1.title=str(self.shotno)+'. '+str(self.namesPA2[i])+' data'
        p1.shotno=[self.shotno]
        
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
                
                if plotAll==True:
                    newPlot=self.plotOfPA2(count,alsoPlotRawAndFit=True)
                else:
                    newPlot=self.plotOfPA2(count,alsoPlotRawAndFit=False)
                newPlot.subtitle=self.namesPA2[count]
                newPlot.yLegendLabel=[]
                sp2[i].append(newPlot)
                count+=1;
        sp1[0][0].title=self.title
        sp2[0][0].title=self.title
        sp1=_plot.subPlot(sp1,plot=False)
        sp2=_plot.subPlot(sp2,plot=False)
        # sp1.shareY=True;
        sp1.plot()
        sp2.plot()
            
    
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
    phi : numpy.ndarray
        toroidal locations for all sensors.  units in radians.
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
        
    """
    def __init__(self,shotno=95540,tStart=0*1e-3,tStop=10*1e-3,plot=False,smoothingAlgorithm='tripleBoxCar'):
        self.shotno = shotno
        self.title = "%d, FB sensors" % shotno

        # sensor names
        self.fbPolNames=[['FB01_S1P', 'FB02_S1P', 'FB03_S1P', 'FB04_S1P', 'FB05_S1P', 'FB06_S1P', 'FB07_S1P', 'FB08_S1P', 'FB09_S1P', 'FB10_S1P'], ['FB01_S2P', 'FB02_S2P', 'FB03_S2P', 'FB04_S2P', 'FB05_S2P', 'FB06_S2P', 'FB07_S2P', 'FB08_S2P', 'FB09_S2P', 'FB10_S2P'], ['FB01_S3P', 'FB02_S3P', 'FB03_S3P', 'FB04_S3P', 'FB05_S3P', 'FB06_S3P', 'FB07_S3P', 'FB08_S3P', 'FB09_S3P', 'FB10_S3P'], ['FB01_S4P', 'FB02_S4P', 'FB03_S4P', 'FB04_S4P', 'FB05_S4P', 'FB06_S4P', 'FB07_S4P', 'FB08_S4P', 'FB09_S4P', 'FB10_S4P']]
        self.fbRadNames=[['FB01_S1R', 'FB02_S1R', 'FB03_S1R', 'FB04_S1R', 'FB05_S1R', 'FB06_S1R', 'FB07_S1R', 'FB08_S1R', 'FB09_S1R', 'FB10_S1R'], ['FB01_S2R', 'FB02_S2R', 'FB03_S2R', 'FB04_S2R', 'FB05_S2R', 'FB06_S2R', 'FB07_S2R', 'FB08_S2R', 'FB09_S2R', 'FB10_S2R'], ['FB01_S3R', 'FB02_S3R', 'FB03_S3R', 'FB04_S3R', 'FB05_S3R', 'FB06_S3R', 'FB07_S3R', 'FB08_S3R', 'FB09_S3R', 'FB10_S3R'], ['FB01_S4R', 'FB02_S4R', 'FB03_S4R', 'FB04_S4R', 'FB05_S4R', 'FB06_S4R', 'FB07_S4R', 'FB08_S4R', 'FB09_S4R', 'FB10_S4R']]
        
        # sensor, toroidal location
        self.phi=_np.pi/180.*_np.array([242.5-360, 278.5-360, 314.5-360, 350.5-360, 26.5, 62.5, 98.5, 134.5, 170.5, 206.5]);#*_np.pi/180.
        
        ## compile full sensor addresses 
        fbPolSensorAddresses=[[],[],[],[]]
        fbRadSensorAddresses=[[],[],[],[]]   
        rootAddress='\HBTEP2::TOP.SENSORS.MAGNETIC:';
        for j in range(0,4):
            for i in range(0,10):
                
                fbPolSensorAddresses[j].append(rootAddress+self.fbPolNames[j][i])
                fbRadSensorAddresses[j].append(rootAddress+self.fbRadNames[j][i])
        
        # get raw data
        self.fbPolRaw=[[],[],[],[]];
        self.fbPolRaw[0], self.fbPolTime =mdsData(shotno,fbPolSensorAddresses[0], tStart, tStop)
        self.fbPolRaw[1], self.fbPolTime =mdsData(shotno,fbPolSensorAddresses[1], tStart, tStop)
        self.fbPolRaw[2], self.fbPolTime =mdsData(shotno,fbPolSensorAddresses[2], tStart, tStop)
        self.fbPolRaw[3], self.fbPolTime =mdsData(shotno,fbPolSensorAddresses[3], tStart, tStop)

        self.fbRadRaw=[[],[],[],[]];
        self.fbRadRaw[0], self.fbRadTime =mdsData(shotno,fbRadSensorAddresses[0], tStart, tStop)
        self.fbRadRaw[1], self.fbRadTime =mdsData(shotno,fbRadSensorAddresses[1], tStart, tStop)
        self.fbRadRaw[2], self.fbRadTime =mdsData(shotno,fbRadSensorAddresses[2], tStart, tStop)
        self.fbRadRaw[3], self.fbRadTime =mdsData(shotno,fbRadSensorAddresses[3], tStart, tStop)        
               
        # data smoothing algorithm
        self.fbPolData=[[],[],[],[]]
        self.fbPolRawFit=[[],[],[],[]]
        self.fbRadData=[[],[],[],[]]
        self.fbRadRawFit=[[],[],[],[]]
        if smoothingAlgorithm=='tripleBoxCar':
            # jeff's triple boxcar smoothing
            for j in range(0,4):
                for i in range(0,10):
                    temp=_process.convolutionSmoothing(self.fbPolRaw[j][i][:],101,'box')
                    temp=_process.convolutionSmoothing(temp,101,'box')
                    temp=_process.convolutionSmoothing(temp,21,'box')
                    self.fbPolRawFit[j].append(temp)
                    self.fbPolData[j].append(self.fbPolRaw[j][i]-temp)
                    
                    temp=_process.convolutionSmoothing(self.fbRadRaw[j][i][:],101,'box')
                    temp=_process.convolutionSmoothing(temp,101,'box')
                    temp=_process.convolutionSmoothing(temp,21,'box')
                    self.fbRadRawFit[j].append(temp)
                    self.fbRadData[j].append(self.fbRadRaw[j][i]-temp)
        else:
            _sys.exit("You must specify a correct smoothing algorithm.  Exiting code...")
     
        # plot
        if plot=='sample':
            self.plotOfSinglePol().plot();
            self.plotOfSingleRad().plot();
        elif plot == True or plot=='all':
            self.plot(True)
            

    def plotOfSinglePol(self, row=0, col=0,plot=True,alsoPlotRawAndFit=True):
        """
        Plots poloidal data from FB sensors
        """
        i=col; j=row;
        
        # initialize plot
        p1=_plot.plot();
        p1.xLabel='time [ms]'
        p1.yLabel=r'dB [G]'
        p1.title=str(self.shotno)+'. '+str(self.fbPolNames[j][i])+' data'
        p1.shotno=[self.shotno]
        
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
        p1=_plot.plot();
        p1.xLabel='time [ms]'
        p1.yLabel=r'dB [G]'
        p1.title=str(self.shotno)+'. '+str(self.fbRadNames[j][i])+' data'
        p1.shotno=[self.shotno]
        
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
        Plots all 80 poloidal and radial FB sensors
        """
        sp1=[[],[],[],[]]
        sp2=[[],[],[],[]]
        count=0
        for i in range(0,4):
            for j in range(0,10):
                if plotAll==True:
                    newPlot=self.plotOfSinglePol(i,j,alsoPlotRawAndFit=True)
                else:
                    newPlot=self.plotOfSinglePol(i,j,alsoPlotRawAndFit=False)                    
                newPlot.subtitle=self.fbPolNames[i][j]
                newPlot.yLegendLabel=[]
                sp1[i].append(newPlot)
                if plotAll==True:
                    newPlot=self.plotOfSingleRad(i,j,alsoPlotRawAndFit=True)
                else:
                    newPlot=self.plotOfSingleRad(i,j,alsoPlotRawAndFit=False)
                newPlot.subtitle=self.fbRadNames[i][j]
                newPlot.yLegendLabel=[]
                sp2[i].append(newPlot)
                count+=1;
        sp1[0][0].title=self.title
        sp2[0][0].title=self.title
        sp1=_plot.subPlot(sp1,plot=False)
        sp2=_plot.subPlot(sp2,plot=False)
        # sp1.shareY=True;
        sp1.plot()
        sp2.plot()    
    
    
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
    smoothingAlgorithm : str
        informs function as to which smoothing algorithm to use on each PA 
        sensor
        
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
        
    def __init__(self,shotno=95540,tStart=0*1e-3,tStop=10*1e-3,plot=False,
                 smoothingAlgorithm='tripleBoxCar'):
        self.shotno = shotno
        self.title = "%d, TA sensor data." % shotno
        
        # names of poloidal and radial sensors
        self.namesTAPol=['TA01_S1P', 'TA01_S2P', 'TA01_S3P', 'TA02_S1P', 'TA02_S2P', 'TA02_S3P', 'TA03_S1P', 'TA03_S2P', 'TA03_S3P', 'TA04_S1P', 'TA04_S2P', 'TA04_S3P', 'TA05_S1P', 'TA05_S2P', 'TA05_S3P', 'TA06_S1P', 'TA06_S2P', 'TA06_S3P', 'TA07_S1P', 'TA07_S2P', 'TA07_S3P', 'TA08_S1P', 'TA08_S2P', 'TA08_S3P', 'TA09_S1P', 'TA09_S2P', 'TA09_S3P', 'TA10_S1P', 'TA10_S2P', 'TA10_S3P'];
        self.namesTARad=['TA01_S2R', 'TA02_S2R', 'TA03_S2R', 'TA04_S2R', 'TA05_S2R', 'TA06_S2R', 'TA07_S2R', 'TA08_S2R', 'TA09_S2R', 'TA10_S2R']
        
        # toroidal locations for the poloidal measurements
        self.phi=_np.pi/180.*_np.array([-117., -108.,  -99.,  -81.,  -72.,  -63.,  -45.,  -36.,  -27.,  -9.,    0.,    9.,   27.,   36.,   45.,   63.,   72.,   81.,   99.,  108.,  117.,  135.,  144.,  153.,  171.,  180.,  189.,    207.,  216.,  225.])    
        
        # toroidal locations for the radial measurements
        self.phiR=_np.pi/180.*_np.array([-108.,  -72.,  -36.,    0.,   36.,   72.,  108.,  144.,  180.,  216.])
        
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
        if smoothingAlgorithm=='tripleBoxCar':
            # jeff's triple boxcar smoothing
            for i in range(0,30):
                temp=_process.convolutionSmoothing(self.taPolRaw[i][:],101,'box')
                temp=_process.convolutionSmoothing(temp,101,'box')
                temp=_process.convolutionSmoothing(temp,21,'box')
                self.taPolRawFit.append(temp)
                self.taPolData.append(self.taPolRaw[i]-temp)
                
                if i < 10:
                    temp=_process.convolutionSmoothing(self.taRadRaw[i][:],101,'box')
                    temp=_process.convolutionSmoothing(temp,101,'box')
                    temp=_process.convolutionSmoothing(temp,21,'box')
                    self.taRadRawFit.append(temp)
                    self.taRadData.append(self.taRadRaw[i]-temp)
        else:
            _sys.exit("You must specify a correct smoothing algorithm.  Exiting code...")
             
        # plot
        if plot=='sample':
            self.plotOfSinglePol().plot();
        elif plot==True or plot=='all':
            self.plot(True);
            
    # TODO Add plotOfSingleRad function
            
    def plotOfSinglePol(self, i=0, alsoPlotRawAndFit=True):
        """ Plot one of the PA1 plots.  based on index, i. """
        p1=_plot.plot();
        p1.xLabel='time [ms]'
        p1.yLabel=r'Gauss'
        p1.title=str(self.shotno)+'. '+str(self.namesTAPol[i])+' data'
        p1.shotno=[self.shotno]
        
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
  

class externalRogowskiData:
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
    
    """
    def __init__(self,shotno=96530,tStart=0*1e-3,tStop=10*1e-3,plot=False):
        self.shotno = shotno
        self.title = "shotno = %d, Ext. Rogowski Data" % shotno
        
        # get data
        data, time=mdsData(shotno=shotno,
                           dataAddress=['\HBTEP2::TOP.SENSORS.EXT_ROGS:EX_ROG_A',
                                        '\HBTEP2::TOP.SENSORS.EXT_ROGS:EX_ROG_B',
                                        '\HBTEP2::TOP.SENSORS.EXT_ROGS:EX_ROG_C',
                                        '\HBTEP2::TOP.SENSORS.EXT_ROGS:EX_ROG_D',],
                           tStart=tStart, tStop=tStop)
        self.eRogA=data[0];
        self.eRogB=data[1];
        self.eRogC=data[2];
        self.eRogD=data[3];
        self.time=time;
        
        if plot == True:
            self.plot()
        
    def plotOfERogA(self):
        # generate rog A plot
        p1=_plot.plot()
        p1.yLabel='A'
        p1.xLabel='time [ms]'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.eRogA,xData=self.time*1000,
                    yLegendLabel='Ext. Rog. A') 
        return p1
    
    def plotOfERogB(self):
        # generate rog B plot
        p1=_plot.plot()
        p1.yLabel='A'
        p1.xLabel='time [ms]'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.eRogB,xData=self.time*1000,
                    yLegendLabel='Ext. Rog. B') 
        return p1
    
    def plotOfERogC(self):
        # generate rog C plot
        p1=_plot.plot()
        p1.yLabel='A'
        p1.xLabel='time [ms]'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.eRogC,xData=self.time*1000,
                    yLegendLabel='Ext. Rog. C') 
        return p1
        
    def plotOfERogD(self):
        # generate rog D plot
        p1=_plot.plot()
        p1.yLabel='A'
        p1.xLabel='time [ms]'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.eRogD,xData=self.time*1000,
                    yLegendLabel='Ext. Rog. D') 
        return p1
        
    def plotOfERogAll(self):        
        # generate rog All plot
        p1=self.plotOfERogA()
        p1.mergePlots(self.plotOfERogB())
        p1.mergePlots(self.plotOfERogC())
        p1.mergePlots(self.plotOfERogD())
        return p1
        
    def plot(self):
        """ Plot all relevant plots """
        self.plotOfERogAll().plot()
        
        
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
    def __init__(self,shotno=98030,tStart=0*1e-3,tStop=10*1e-3,plot=False):
        self.shotno = shotno
        self.title = "shotno = %d, Spectrometer Data" % shotno
        
        # get data
        data, self.time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.SPECTROMETER'],
                              tStart=tStart, tStop=tStop)
        self.spect=data[0];
#        self.time=time;
        
        if plot == True:
            self.plot()
        
    def plotOfSpect(self):
        # generate plot
        p1=_plot.plot()
        p1.yLabel='V'
        p1.xLabel='time [ms]'
        p1.subtitle='Spectrometer Intensity'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.spect,xData=self.time*1000) 
        return p1
    
    def plot(self):
        """ Plot all relevant plots """
        self.plotOfSpect().plot()
        
        
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
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    title : str
        title to go on all plots
    sensorNames : list (of str)
        names of each SOL sensor
    solData : list (of numpy.ndarray)
        SOL sensor data
    time : numpy.ndarray
        time data
        
    Subfunctions 
    ------------
    plotOfSingleSensor : _plotTools.plot
        returns plot of single sensor
    plot :
        plots all SOL data
    
    """
    def __init__(self,shotno=98030,tStart=0*1e-3,tStop=10*1e-3,plot=False):
        self.shotno = shotno
        self.title = "shotno = %d, SOL Data" % shotno
        self.sensorNames = ['LFS01_S1', 'LFS01_S2', 'LFS01_S3', 'LFS01_S4', 'LFS01_S5', 'LFS01_S6', 'LFS01_S7', 'LFS01_S8', 'LFS04_S1', 'LFS04_S2', 'LFS04_S3', 'LFS04_S4', 'LFS08_S1', 'LFS08_S2', 'LFS08_S3', 'LFS08_S4', 'LFS08_S5', 'LFS08_S6', 'LFS08_S7', 'LFS08_S8']

        # compile list of sensor addresses for all 20 SOL tiles
        sensorAddress=[]
        for i in range(0,len(self.sensorNames)):
            sensorAddress.append('\HBTEP2::TOP.SENSORS.SOL:%s' % self.sensorNames[i]) 
            
        # get data
        self.solData, self.time=mdsData(shotno=shotno,
                              dataAddress=sensorAddress,
                              tStart=tStart, tStop=tStop)
                              
        if plot == True:
            self.plot()
                              
    def plotOfSingleSensor(self,index): #name='LFS01_S1'
        """ returns plot of a single sol sensor """
#        index = self.sensorNames.index(name)
        # generate plot
        p1=_plot.plot()
        p1.yLabel='V'
        p1.xLabel='time [ms]'
        p1.subtitle=self.sensorNames[index]
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.solData[index],xData=self.time*1000) 
            
        return p1
        
    def plot(self):
        """ plots all 20 sol sensor currents """
        plots=[[],[],[],[]]
        count=0
        for i in range(0,4):
            for j in range(0,5):
                newPlot=self.plotOfSingleSensor(count) #name=self.sensorNames[count]
                newPlot.subtitle=self.sensorNames[count]
                plots[i].append(newPlot)
                count+=1;
        sp1=_plot.subPlot(plots,plot=False)
        sp1.shareY=True;
        sp1.plot()
    
        
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
    def __init__(self,shotno=96530,tStart=0*1e-3,tStop=10*1e-3,plot=False):
        self.shotno = shotno
        self.title = "%d, Loop voltage data." % shotno
        
        # get data
        data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.LOOP_VOlTAGE'],
                              tStart=tStart, tStop=tStop)   
        self.loopVoltage=data[0];
        self.time=time;
        
        if plot == True:
            self.plot()
        
    def plotOfLoopVoltage(self):
        # generate plot
        p1=_plot.plot()
        p1.yLabel='V'
        p1.xLabel='time [ms]'
        p1.subtitle='Loop Voltage'
        p1.title=str(self.shotno);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.loopVoltage,xData=self.time*1000) 
        p1.yLim=[0,15]  # using same axis-limits as hbtplot.py
        return p1
            
    def plot(self):
        """ Plot all relevant plots """
        self.plotOfLoopVoltage().plot()
        
        
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
    means that the TF data has fewer poitns than the rest of the HBTEP data, 
    and this makes comparing data difficult.  
    """
    def __init__(self,shotno=96530,tStart=None,tStop=None,plot=False,upSample=False):
        self.shotno = shotno
        self.title = "shotno = %d, TF Field Data" % shotno
        
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
        p1=_plot.plot()
        p1.yLabel='T'
        p1.xLabel='time [ms]'
        p1.subtitle='TF Bank Field'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.tfBankField,xData=self.time*1000) 
        return p1
            
    def plot(self):
        """ Plot all relevant plots """
        self.plotOfTF().plot()
        
        
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
        self.title = "shotno = %d, Capacitor Bank Data" % shotno
        
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
        
    def plotOfVF(self):
        # generate vf plot
        p1=_plot.plot()
        p1.yLabel='kA'
        p1.xLabel='time [ms]'
        p1.subtitle='VF Bank Current'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.vfBankCurrent,xData=self.vfTime*1000) 
        return p1
        
    def plotOfOH(self):
        # generate oh plot
        p1=_plot.plot()
        p1.yLabel='kA'
        p1.xLabel='time [ms]'
        p1.subtitle='OH Bank Current'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.yLim=[-3.5e1,3.5e1]
        p1.addTrace(yData=self.ohBankCurrent,xData=self.ohTime*1000) 
        return p1
        
    def plotOfSH(self):
        # generate sh plot
        p1=_plot.plot()
        p1.yLabel='kA'
        p1.xLabel='time [ms]'
        p1.subtitle='SH Bank Current'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.addTrace(yData=self.shBankCurrent,xData=self.shTime*1000) 
        return p1
        
            
    def plot(self):
        """ Plot all relevant plots """
        # load tf data because it's always nice to compare it with the cap banks
        tf=tfData(shotno=self.shotno,tStart=None,tStop=None)
        
        # subplot of all 4 bank data
        # note: most of this code is determining and setting the appropriate
        # x and y limits
        sp1=_plot.subPlot([self.plotOfVF(),self.plotOfOH(),self.plotOfSH(),
                        tf.plotOfTF()],plot=False)
        xMin=_np.min(sp1.subPlots[0].xData)
        xMax=_np.max(sp1.subPlots[0].xData)
        sp1.subPlots[0].xLim=[xMin,xMax]
        subData=sp1.subPlots[3].yData[0]
        subTime=sp1.subPlots[3].xData[0]
        iMin=_process.findNearest(subTime,xMin)
        iMax=_process.findNearest(subTime,xMax)
        yMin=_np.min(subData[iMin:iMax])
        yMax=_np.max(subData[iMin:iMax])
        dY=yMax-yMin
        sp1.subPlots[3].yLim=[yMin,yMax+0.25*dY] #[yMin-0.25*dY,yMax+0.25*dY]
        sp1.plot()
        
        # plot of TF with shaded region representing x-limits in subplot
        p1=tf.plotOfTF()
        p1.axvspan=[self.ohTime[0]*1e3,self.ohTime[-1]*1e3]
        p1.axvspanColor=['r']
        p1.plot()
        
        
     
#####################################################
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
    
    def __init__(self,shotno=95782,tStart=0*1e-3,tStop=10*1e-3, plot=False, probeRadius=[]):
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
        capBank=capBankData(shotno=shotno,tStart=tStart,tStop=tStop)
        vf=capBank.vfBankCurrent
        oh=capBank.ohBankCurrent
        self.time=capBank.vfTime
        
        # get plasma current
        ip=ipData(shotno=shotno,tStart=tStart,tStop=tStop)
        ip=ip.ip*1212*1e-9  # ip gain
        
        # get cos-1 raw data
        cos1=cos1RogowskiData(shotno=shotno,tStart=tStart,tStop=tStop+2e-06) # note that the cumtrapz function below loses a data point.  by adding 2e-06 to the time, i start with an additional point that it's ok to lose
        
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

        
    def plotOfMajorRadius(self):
        # limited data
        innerLimiter=_np.array([0.90296,0.90296])*100
        innerLimiterTime=_np.array([self.time[0],self.time[-1]]);
        outerLimiter=_np.array([0.92,0.92])*100        
        outerLimiterTime=_np.array([self.time[0],self.time[-1]]);
        
        # plot
        p1=_plot.plot()
        p1.subtitle='major radius'
        p1.title=self.title
        p1.shotno=[self.shotno]
        p1.xLabel='time [ms]'
        p1.yLabel='cm'
        p1.yLim=[89, 95]
        p1.addTrace(yData=self.majorRadius*100,xData=self.time*1000,
                    yLegendLabel='major radius') 
        p1.addTrace(yData=innerLimiter,xData=innerLimiterTime*1000,
                    yLegendLabel='HFS limited') 
        p1.addTrace(yData=outerLimiter,xData=outerLimiterTime*1000,
                    yLegendLabel='LFS limited') 
        return p1
        
    def plotOfMinorRadius(self):
        p1=_plot.plot()
        p1.subtitle='minor radius'
        p1.title=self.title
        p1.shotno=[self.shotno]
        p1.xLabel='time [ms]'
        p1.yLabel='cm'
        p1.yLim=[10, 16]
        p1.addTrace(yData=self.minorRadius*100,xData=self.time*1000,
                    yLegendLabel='minor radius') 
        return p1

    def plot(self):
        self.p=_plot.subPlot([self.plotOfMajorRadius(),self.plotOfMinorRadius()]);


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
    
    def __init__(self,shotno=96496, tStart=0*1e-3, tStop=10*1e-3, plot=False):
        self.shotno = shotno
        self.title = r"shotno = %d, q$^*$ Data" % shotno
        
        # get data
        ip=ipData(shotno=shotno,tStart=tStart,tStop=tStop)
        plasmaRadius=plasmaRadiusData(shotno=shotno,tStart=tStart,tStop=tStop)
        tfProbeData,tfProbeTime=mdsData(shotno=96496,
                                        dataAddress=['\HBTEP2::TOP.SENSORS.TF_PROBE'],
                                        tStart=tStart,tStop=tStop)
                                        
        # upsample tfprobe data  (its recorded on the A14)      
        data=_process.upSampleData(ip.time,tfProbeTime,tfProbeData[0])
        
        # more tf calculations
        tfProbeData=data*1.23/plasmaRadius.majorRadius
        
        # calc q star
        self.qStar= plasmaRadius.minorRadius**2 * tfProbeData / (2e-7 * ip.ip * plasmaRadius.majorRadius)
        self.time=ip.time
        
                
        if plot == True:
            self.plot()
        
        
    def plotOfQStar(self):
        """
        returns the plot of IP vs time
        """
        p1=_plot.plot()
        p1.yLabel=''
        p1.xLabel='time [ms]'
        p1.subtitle=r'q$^*$'
        p1.title=str(self.title);
        p1.shotno=[self.shotno]
        p1.yLim=[2,5]
        p1.addTrace(yData=self.qStar,xData=self.time*1000) 
        
        return p1
        
            
    def plot(self):
        """ 
        Plot all relevant plots 
        """
        self.plotOfQStar().plot()
        
###############################################################################
### sensor black list data.  presently not used anywhere
    
def checkBlackList(inData,inName):
    # TODO(John) this needs an overhaul
    """
    Takes in data and sensor name.  Checks sensor name against blacklist.  
    If bad sensor, return all zeros # and a true boolean.  
    Otherwise, returns original data # and a false boolean.
    
    """
    if inName in SENSORBLACKLIST==True:
        outData=_np.zeros(inData.size);
    else:
        outData=inData;
        
    return outData
    
    
################################################################################
#### HBTEP shot data format
    
def shotData(shotnos=_np.array([98119,98120],dtype=int), #,
                     tStart=0.0*1e-3, tStop=8.0*1e-3, saveData=False, 
                     tags=['nMode','bp','tp','eRog','spect','sol','lv','ip','radius','qStar']):

    """
    Loads multiple shot data files.  
    
    Parameters
    ----------
    
    Returns 
    -------
    data : list (of _getHBTData.shotData)
        list of shot data structures
    """
                     
    # make sure shotnos is a np array
    if type(shotnos) is int or type(shotnos) is float:
        shotnos=_np.array([shotnos]);
    elif type(shotnos) is list:
        shotnos=_np.array(shotnos);
    m=len(shotnos)
        
    # make sure tStart and tStop are also lists (actually np arrays)
    if type(tStart) is int or type(tStart) is float: 
        tStart=_np.ones(m)*tStart
    if type(tStop) is int or type(tStop) is float: 
        tStop=_np.ones(m)*tStop
        
    # load each shotData
    data=[]
    for i in range(0,m):
        data.append(_getShotData(shotno=shotnos[i],
                                tStart=tStart[i],
                                tStop=tStop[i],
                                saveData=saveData,
                                tags=tags))
                                         
    # return
    return data
                 
    
class _getShotData:
    """
    Downloads HBTEP data into a single class for a single shotno.  
    If the file is present locally, it loads it from there.  If the file is NOT
    present at the local save directory, the data is loaded from the HBTEP 
    server
    
    Parameters
    ----------
    
    Attributes
    ----------
    
    """
    def __init__(self,shotno=98030,tStart=0.0*1e-3,tStop=8.0*1e-3,saveData=False,
                 tags=['nMode','bp','tp','eRog','spect','sol','lv','ip','radius','qStar']):
        self.shotno=shotno
        # download data
        
        # check if file is alreay downloaded.  download permanently otherwise.  also downloads file if specifically requested.
        
        localFilePath=_pref._LOCAL_DATA_DIR+'shotData_'+str(shotno)+'.pickle'
        filePresent =_os.path.isfile(localFilePath)

        # load from existing file        
        if filePresent==True: # and forceDownload!=True:
            saveData=False
            data=_rwData.loadFromPickle(localFilePath)
            try:
                if 'nMode' in tags:
                    self.nMode=data.nMode
                if 'bp' in tags:
                    self.bp=data.bp
                if 'tp' in tags:
                    self.tp=data.tp
                if 'eRog' in tags:
                    self.eRog=data.eRog
                if 'spect' in tags:
                    self.spectrometer=data.spectrometer
                if 'sol' in tags:
                    self.sol=data.sol
                if 'lv' in tags:
                    self.lv=data.lv
                if 'ip' in tags:
                    self.ip=data.ip
                if 'capBank' in tags:
                    self.capBank=data.capBank
                if 'radius' in tags:
                    self.radius=data.radius
                if 'qStar' in tags:
                    self.qStar=data.qStar
                if 'gpu' in tags:
                    self.gpu=data.gpu
                print "loading "+str(int(shotno))+" shot data from local storage"
                    
            except AttributeError:
                print str(int(shotno))+" shot data is incomplete. " + "Getting from HBTEP server instead and rewriting old data"
                filePresent=False
#                saveData=True
            
        # if file is not already downloaded locally, download from HBTEP serer
        if filePresent==False: # or forceDownload==True:
            saveData=True
            print "downloading "+str(int(shotno))+" shot data from HBTEP server"
            if 'nMode' in tags:
                self.nMode=nModeData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'bp' in tags:
                self.bp=bpData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'tp' in tags:
                self.tp=tpData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'eRog' in tags:
                self.eRog=externalRogowskiData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'spect' in tags:
                self.spectrometer=spectrometerData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'sol' in tags:
                self.sol=solData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'lv' in tags:
                self.lv=loopVoltageData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'ip' in tags:
                self.ip=ipData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'capBank' in tags:
                self.capBank=capBankData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'radius' in tags:
                self.radius=plasmaRadiusData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'qStar' in tags:
                self.qStar=qStarData(shotno=shotno,tStart=tStart,tStop=tStop)
            if 'gpu' in tags:
                import johnsUnderDevelopementToolbox as john
                self.gpu=john.gpuControlData(shotno=shotno,tStart=tStart,tStop=tStop,operatingMode='frequencyControl')
#                self.gpu
#                self.gpu=gpu.fb
        
        # save data to file
        if saveData==True:
            print "Writing to %s" % localFilePath
            _rwData.saveToPickle(self,localFilePath)
            
  
            
###############################################################################
### Processed data from HBTEP
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
        
    def __init__(self,shotno=96530,tStart=0*1e-3,tStop=10*1e-3,plot=False,nModeSensor='FB',method='leastSquares',phaseFilter='gaussian',frequencyFilter=''):
        
        self.shotno=shotno
        self.title = 'shotno = %d.  %s sensor.  n mode analysis' % (shotno,nModeSensor)
        self.nModeSensor=nModeSensor
#        self.frequencyFilter=frequencyFilter
#        self._phaseFilter=phaseFilter

        if nModeSensor=='TA':
            ## load TA data
            temp=taData(self.shotno,tStart,tStop+0.5e-3);  # asking for an extra half millisecond (see Notes above) 
            data=temp.taPolData
            self.time=temp.tbPolTime
            phi=temp.phi
            [n,m]=_np.shape(data)
        elif nModeSensor=='FB':
            ## load FB data
            temp=fbData(self.shotno,tStart=tStart,tStop=tStop+0.5e-3);  # asking for an extra half millisecond (see Notes above) 
            data=temp.fbPolData[0]  ## top toroidal array = 0
            self.time=temp.fbPolTime
            phi=temp.phi
            [n,m]=_np.shape(data)
        self._data=data
        self._phi=phi

        if method=='leastSquares':
            ## Construct A matrix and its inversion
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
            # TODO(John): remove for loop and convert into all matrix math 
            #             Should simplify code and make it run faster
            for j in range(0,m):
                y=_np.zeros(n);
                for i in range(0,n):
                    y[i]=data[i][j]*1e4
                x[:,j]=Ainv.dot(y)
                self.n1Amp[j]=_np.sqrt(x[1,j]**2+x[2,j]**2)
                self.n2Amp[j]=_np.sqrt(x[3,j]**2+x[4,j]**2)
                self.n1PhaseRaw[j]=_np.arctan2(x[1,j],x[2,j])
            self._x=x
            self.n1PhaseRaw*=-1  # for some reason, the slope of phase had the wrong sign.  this corrects that.

            
        else:
            _sys.exit("Invalid mode analysis method requested.")
            
        # filter phase
        self.n1Phase=_np.zeros(len(self.n1PhaseRaw))   
        if phaseFilter == 'gaussian':
            self.n1Phase=_process.wrapPhase(
                            _process.convolutionSmoothing(
                                _process.unwrapPhase(self.n1PhaseRaw),101,'gaussian'))
        elif phaseFilter == 'box' or phaseFilter == 'boxcar':
            self.n1Phase=_process.wrapPhase(
                            _process.convolutionSmoothing(
                                _process.unwrapPhase(self.n1PhaseRaw),101,'box'))
        elif phaseFilter == 'savgol':
            self.n1Phase=_process.wrapPhase(
                            _process.savgolFilter(
                                _process.unwrapPhase(self.n1PhaseRaw),41,1))
            
        elif phaseFilter == '' or phaseFilter == None:
            self.n1Phase=_copy(self.n1PhaseRaw)
        else:
            _sys.exit("Invalid phase filter requested.")
        
            
            
        ## Calculate frequency (in Hz) using second order deriv 
        self.n1FreqRaw=_np.gradient(_process.unwrapPhase(self.n1Phase))/_np.gradient(self.time)/(2*_np.pi)        
        
        # filter frequency
        self.n1Freq=_np.zeros(len(self.n1FreqRaw)) 
        if 'Butterworth' in frequencyFilter or 'butterworth' in frequencyFilter:
            
            filterOrder=1
            cutoffFreq=1000

            # implement filter
            self.n1Freq=_process.butterworthFilter(self.n1FreqRaw,
                                                           self.time,
                                                           filterOrder=filterOrder,
                                                           samplingRate=1./(2*1e-6),
                                                           cutoffFreq=cutoffFreq,
                                                           filterType='low')
                                                           
        elif frequencyFilter=='boxcar' or frequencyFilter=='boxCar':
            self.n1Freq=_process.convolutionSmoothing(self.n1FreqRaw,81,'box')
            
        elif frequencyFilter=='gaussian':
            self.n1Freq=_process.convolutionSmoothing(self.n1FreqRaw,101,'gaussian')
            
        elif frequencyFilter=='' or frequencyFilter==None:
            self.n1Freq=_copy(self.n1FreqRaw)
            
        else:
            _sys.exit("Invalid frequency filter requested.")
            
            
        # trim off extra half millisecond (see Notes)
        self.time, temp=_trimTime(self.time,
                                  [self.n1Amp,self.n2Amp,self.n1Phase,
                                   self.n1PhaseRaw,self.n1Freq,self.n1FreqRaw],
                                  tStart,tStop)
        self.n1Amp=temp[0]
        self.n2Amp=temp[1]
        self.n1Phase=temp[2]
        self.n1PhaseRaw=temp[3]
        self.n1Freq=temp[4]
        self.n1FreqRaw=temp[5]
        
        ## plot data
        if plot==True:
            self.plot(includeRaw=True)
            
        elif plot == 'all':
            self.plotOfSlice(index=int(m/4)).plot();
            self.plotOfSlice(index=int(m/2)).plot();
            self.plotOfAmps().plot()
            self.plot(includeRaw=True)
            
    def plot(self,includeRaw=True):
        """
        plots and returns a subplot of n=1 mode amplitude, phase, and frequency
        
        Parameters
        ----------
        includeRaw : bool
            if True, also plots the raw (unfiltered) phase and frequency
        """
        sp1=_plot.subPlot([self.plotOfN1Amp(),self.plotOfN1Phase(),
                           self.plotOfN1Freq()],plot=False)
        if includeRaw==True:
            # add phase raw data
            sp1.subPlots[1].addTrace(yData=self.n1PhaseRaw,
                                     xData=self.time*1000,
                                     linestyle='',
                                     marker='.',
                                     yLegendLabel='raw')
            
            # add frequency raw data
            sp1.subPlots[2].addTrace(yData=self.n1FreqRaw/1000.,
                                     xData=self.time*1000,
                                     yLegendLabel='raw')
            
        sp1.plot()
        return sp1
        
    def plotOfAmps(self):
        ## mode amplitude plots  
        p1=_plot.plot()
        p1.addTrace(yData=self.n1Amp,xData=self.time*1000,
                    yLegendLabel='n=1') 
        p1.addTrace(yData=self.n2Amp,xData=self.time*1000,
                    yLegendLabel='n=2') 
        p1.title=self.title
        p1.shotno=[self.shotno]
        p1.xLabel='ms'
        p1.yLabel='G'   
        return p1

    def plotOfN1Amp(self):
        # n=1 mode amplitude
        p1=_plot.plot()
        p1.subtitle='Mode amplitude, n=1'
        p1.yLim=[0,10]
        p1.addTrace(yData=self.n1Amp,xData=self.time*1000) 
        p1.title=str(self.title)
        p1.shotno=[self.shotno]
        p1.xLabel='ms'
        p1.yLabel='G'
        return p1

    def  plotOfN1Phase(self):
        # n=1 mode phase
        p1=_plot.plot()
        p1.subtitle='Mode phase, n=1'
        p1.yLim=[-_np.pi,_np.pi]
        p1.title=self.title
        p1.shotno=[self.shotno]
        p1.xLabel='ms'
        p1.yLabel='phi'
        p1.addTrace(yData=self.n1Phase,xData=self.time*1000,
                    marker='.',linestyle='',yLegendLabel='filtered') 
        
        return p1
        
    def plotOfN1Freq(self):
        # n=1 mode freq
        p1=_plot.plot()        
        p1.subtitle='Mode frequency, n=1'
        p1.title=self.title
        p1.shotno=[self.shotno]
        p1.xLabel='ms'
        p1.yLabel='kHz'
        p1.yLim=[-20,20]        
        p1.addTrace(yData=self.n1Freq/1000.,xData=self.time*1000,
                    yLegendLabel='filtered') 
        
        return p1
        
    def plotOfPhaseAmp(self):
                   
        # hybrid plot of phase AND amplitude
        # TODO(John) implement in new plot function
                   # TODO(john) subfunction needs overhaul
        p1=_plot.plot() 
        p1.yData=[self.n1Phase]
        p1.xData=[self.time*1000]
        p1.colorData=[self.n1Amp]#[self.n1Amp]
        p1.linestyle=['']
        p1.marker=['.']
        p1.subtitle='n=1 Phase and Filtered Amplitude'
        p1.title=self.title
        p1.shotno=[self.shotno]
        p1.xLabel='ms'
        p1.yLabel=r'$\phi$'
        p1.zLabel='Gauss'
        p1.yLegendLabel=['TA sensors']
        p1.plotType='scatter'
        p1.yLim=[-_np.pi,_np.pi]
        return p1      
#        mx=_np.max(self.n1AmpFiltered)
#        lCutoff=2.5
#        uCutoff=8.
#        cm = _processt.singleColorMapWithLowerAndUpperCutoffs(lowerCutoff=lCutoff/mx,upperCutoff=uCutoff/mx)
#        self.plotOfPhaseAmp.cmap=cm
                        
    def plotOfSlice(self,index=0):
        """
        Plots fit data for a single time value
        """
        j=index;
        [n,m]=_np.shape(self._data)
        y=_np.zeros(n);
        for i in range(0,n):
                y[i]=self._data[i][j]*1e4
        p1=_plot.plot()
        phi=_np.linspace(self._phi[0],self._phi[-1],100)
        n1Fit=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)
        n2Fit=self._x[0,j]+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)
        fitTotal=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)

        # plot
        p1.addTrace(yData=y,xData=self._phi,
                    marker='.',linestyle='',yLegendLabel='raw') 
        p1.addTrace(yData=n1Fit,xData=phi,
                    yLegendLabel='n=1') 
        p1.addTrace(yData=n2Fit,xData=phi,
                    yLegendLabel='n=2') 
        p1.addTrace(yData=fitTotal,xData=phi,
                    yLegendLabel='n=3') 
        p1.title='t='+str(self.time[j]*1000)+'ms.  '+str(self.shotno)
        p1.shotno=[self.shotno]
        return p1
        

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
    def __init__(self,shotno=96530,tStart=0*1e-3,tStop=10*1e-3,plot=False,theta0=0,sensor='PA1'):
        self.shotno=shotno
        self.title= 'shotno = %d.  sensor = %s.  m mode analysis' % (shotno, sensor)
#        self.tStart=tStart
#        self.tStop=tStop
#        self.theta0=theta0
        
        if sensor=='PA1':
            data=paData(self.shotno,tStart=tStart,tStop=tStop);
            self._data=data.pa1Data
            self.time=data.pa1Time
            self._theta=data.theta
            [n,m]=_np.shape(self._data)
        if sensor=='PA2':
            data=paData(self.shotno,tStart=tStart,tStop=tStop,sensor='PA2');
            self._data=data.pa2Data
            self.time=data.time
            self._theta=data.theta
            [n,m]=_np.shape(self._data)

        ## Construct A matrix and its inversion
        A=_np.zeros((n,11))
        A[:,0]=_np.ones(n);
        A[:,1]=_np.sin(self._theta)
        A[:,2]=_np.cos(self._theta)
        A[:,3]=_np.sin(2*self._theta)
        A[:,4]=_np.cos(2*self._theta)
        A[:,5]=_np.sin(3*self._theta)
        A[:,6]=_np.cos(3*self._theta)
        A[:,7]=_np.sin(4*self._theta)
        A[:,8]=_np.cos(4*self._theta)
        A[:,9]=_np.sin(5*self._theta)
        A[:,10]=_np.cos(5*self._theta)
        Ainv=_np.linalg.pinv(A)
        
        ## Solve for coefficients, x, for every time step and assign values to appropriate arrays 
        self._x=_np.zeros([11,m]);
        # self.m0Offset=_np.zeros(m)
        self.m1Amp=_np.zeros(m)
        self.m1Phase=_np.zeros(m)
        self.m2Amp=_np.zeros(m)
        self.m2Phase=_np.zeros(m)   
        self.m3Amp=_np.zeros(m)
        self.m3Phase=_np.zeros(m)     
        self.m4Amp=_np.zeros(m)
        self.m4Phase=_np.zeros(m)     
        self.m5Amp=_np.zeros(m)
        self.m5Phase=_np.zeros(m)          
        for j in range(0,m):
            y=_np.zeros(n);
            for i in range(0,n):
                y[i]=self._data[i][j]*1e4
            self._x[:,j]=Ainv.dot(y)
            # self.m0Offset=self.x[0,j]
            self.m1Amp[j]=_np.sqrt(self._x[1,j]**2+self._x[2,j]**2)
            self.m2Amp[j]=_np.sqrt(self._x[3,j]**2+self._x[4,j]**2)
            self.m3Amp[j]=_np.sqrt(self._x[5,j]**2+self._x[6,j]**2)
            self.m4Amp[j]=_np.sqrt(self._x[7,j]**2+self._x[8,j]**2)
            self.m5Amp[j]=_np.sqrt(self._x[9,j]**2+self._x[10,j]**2)
            self.m1Phase[j]=_np.arctan2(self._x[1,j],self._x[2,j])
            self.m2Phase[j]=_np.arctan2(self._x[3,j],self._x[4,j])
            self.m3Phase[j]=_np.arctan2(self._x[5,j],self._x[6,j])
            self.m4Phase[j]=_np.arctan2(self._x[7,j],self._x[8,j])
            self.m5Phase[j]=_np.arctan2(self._x[9,j],self._x[10,j])
#        m1PhaseUnwrapped=_process.unwrapPhase(self.m1Phase)
#        m2PhaseUnwrapped=_process.unwrapPhase(self.m2Phase)
#        m3PhaseUnwrapped=_process.unwrapPhase(self.m3Phase)
#        m4PhaseUnwrapped=_process.unwrapPhase(self.m4Phase)
#        m5PhaseUnwrapped=_process.unwrapPhase(self.m5Phase)
        
        # TODO:  add frequency and phase data, raw and filtered
        
        
        if plot == True:
            self.plotOfAmplitudes().plot()
        elif plot == 'all':
            self.plotOfSlice(index=int(m/4)).plot();
            self.plotOfSlice(index=int(m/2)).plot();
            self.plotOfAmplitudes().plot()
            
        
    def plotOfAmplitudes(self):
        # plot amplitudes
        p1=_plot.plot()
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
        p1.title=str(self.title)
        p1.shotno=[self.shotno]
        p1.xLabel='ms'
        p1.yLabel='G'
        return p1
        
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
        p1=_plot.plot()
        theta=_np.linspace(self._theta[0],self._theta[-1],100)
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
        p1.title='t=%.3f ms. %s ' % (self.time[j]*1000, self.title)
        p1.shotno=[self.shotno]
#        p1.plot()
        return p1
        

###############################################################################
### debugging code

def _debugPlotExamplesOfAll():
    """ 
    This code plots an example of most every function in this file.  
    Effectively, this allows the testing of most every function all in one go.
    """
    bpData(plot=True)
    capBankData(plot=True)
    cos1RogowskiData(plot=True)
    externalRogowskiData(plot=True)
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
    
