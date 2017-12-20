"""
_rwHBTData.py - load HBT data, both processed and unprocessed

NOTES
-----
-The convention for units is to maintain everything in SI until they are
plotted.  
-A few of these functions are merely wrappers for other people's code
"""

###############################################################################
### import libraries

# common libraries 
import numpy as _np
import MDSplus as _mds
from copy import copy as _copy
import pickle as _pk
import sys as _sys
import _socket

# hbtepLib libraries
import _rwDataTools as _rwData
import _processData as _process
import _plotTools as _plot
try:
    import _hbtPreferences as _pref
except ImportError:
    _sys.exit("Code hault: _hbtPreferences.py file not found.  See readme.md concerning the creation of hbtPreferences.py")
        
        
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
    # determine indices of cutoff regions
    iStart=_process.find_nearest(time,tStart);   # index of lower cutoff
    iStop=_process.find_nearest(time,tStop);     # index of higher cutoff
    
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
            tStart=None,tStop=None):
    """
    Get data and optionally associated time from MDSplus tree
    
    Parameters
    ----------
#    mdsConn : MDSplus.connection
#        connection class to mdsplus tree.  i.e. the output from
#        initMDSConnection 
    shotno : int
        if shotno is defined, this function will establish its own mdsConn
#        regardless of whether mdsConn was previously defined
    dataAddress : list (of strings)
        address of desired data on tree
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
    if type(shotno) is float or type(shotno) is int:
        mdsConn=_initMDSConnection(shotno);
        
#    return mdsConn
    if type(dataAddress) is not list:
        dataAddress=[dataAddress];
    
    # get data and time
    data = [];
    for i in range(0,len(dataAddress)):
        data.append(mdsConn.get(dataAddress[i]).data())
    time = mdsConn.get('dim_of('+dataAddress[0]+')').data();

    # if tStart is not defined, give it the default values
    if type(tStart) is not int and type(tStart) is not float:
        tStart = _TSTART;
        tStop  = _TSTOP;
        
    # check to see if units are in seconds and NOT in milliseconds
    if tStop > 1:
        tStart=tStart*1e-3;
        tStop=tStop*1e-3;
        
    # trim time and data
    time,data= _trimTime(time,data,tStart,tStop)
        
    # return data and time
    return data, time
    # TODO is there a way to make time an optional return??
    
    
###############################################################################
### get device specific data
    
class ipData:
    """
    Gets plasma current data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
    tStart : float
        time (in seconds) to drim data after
    plot : bool
        plots all relevant plots if true
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    ip : numpy.ndarray
        plasma current data
    time : numpy.ndarray
        time data
    plotOfIP : 
        custom plot function
    
    """
    def __init__(self,shotno=96530,tStart=None,tStop=None,plot=False):
        self.shotno = shotno
        
        # get data
        data, time=mdsData(shotno=shotno,
                              dataAddress=['\HBTEP2::TOP.SENSORS.ROGOWSKIS:IP'],
                              tStart=tStart, tStop=tStop)
        self.ip=data[0];
        self.time=time;
        
        # generate plot
        self.plotOfIP=_plot.plot()
        self.plotOfIP.yData=[self.ip*1e-3]
        self.plotOfIP.xData=[self.time*1000]
#        self.p1.xLim=[self.tStart,self.tStop]
        self.plotOfIP.yLabel='kA'
        self.plotOfIP.xLabel='time [ms]'
        self.plotOfIP.subtitle='Plasma Current'
        self.plotOfIP.title=str(self.shotno);
        
        if plot == True:
            self.plot()
            
    def plot(self):
        """ Plot all relevant plots """
        self.plotOfIP.plot()
        
  
class bpData:
    """
    Downloads bias probe data from both probes.  
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
    tStart : float
        time (in seconds) to drim data after
    plot : bool
        plots all relevant plots if true
        
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
    plotOfGPUVoltageRequest : _plotTools.plot 
        Plot of gpu request voltage (as measured by the CPCI)
    plotOfVoltage : _plotTools.plot 
        Plot of both bias probe voltages
    plotOfCurrent : _plotTools.plot 
        Plot of both bias probe currents
    plotOfBPS9Voltage : _plotTools.plot 
        Plot of bps9 voltage only
    plotOfBPS9Current : _plotTools.plot 
        Plot of bps9 current only
        
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
    def __init__(self,shotno=96530,tStart=None,tStop=None,plot=False):
        self.shotno = shotno
        self.title = "shotno=%s, BP Data." % shotno

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

        # initialize gpu voltage request plot
        self.plotOfGPUVoltageRequest=_plot.plot();
        self.plotOfGPUVoltageRequest.title=self.title
        self.plotOfGPUVoltageRequest.yData=[self.bps9GPURequestVoltage]
        self.plotOfGPUVoltageRequest.xData=[self.time*1000]
        self.plotOfGPUVoltageRequest.yLabel='V'
        self.plotOfGPUVoltageRequest.xLabel='ms'
        self.plotOfGPUVoltageRequest.subtitle='Voltage Request from GPU (pre-amplifier)'
        self.plotOfGPUVoltageRequest.title=str(self.shotno);
        self.plotOfGPUVoltageRequest.yLegendLabel='BPS9'

        # initialize voltage plot
        self.plotOfVoltage=_plot.plot();
        self.plotOfVoltage.title=self.title
        self.plotOfVoltage.yLim=[-200,200] # default [-200,200]
        self.plotOfVoltage.yLabel='V'
        self.plotOfVoltage.xLabel='Time [ms]'
        self.plotOfVoltage.subtitle='BP Voltage'
        self.plotOfVoltage.yData.append(self.bps9Voltage)
        self.plotOfVoltage.xData.append(self.time*1000)
        self.plotOfVoltage.yLegendLabel.append('BPS9')
        self.plotOfVoltage.yData.append(self.bps5Voltage)
        self.plotOfVoltage.xData.append(self.time*1000)
        self.plotOfVoltage.yLegendLabel.append('BPS5')
        
        # initialize current plot
        self.plotOfCurrent=_plot.plot();
        self.plotOfCurrent.title=self.title
        self.plotOfCurrent.yLim=[-5,60]  # default [-10,70]
        self.plotOfCurrent.yLabel='A'
        self.plotOfCurrent.subtitle='BP Current'
        self.plotOfCurrent.xLabel='Time [ms]'
        self.plotOfCurrent.yData.append(self.bps9Current)
        self.plotOfCurrent.xData.append(self.time*1000)
        self.plotOfCurrent.yLegendLabel.append('BPS9')
        self.plotOfCurrent.yData.append(self.bps5Current)
        self.plotOfCurrent.xData.append(self.time*1000)
        self.plotOfCurrent.yLegendLabel.append('BPS5')
                
        # initialize BPS9 voltage plot
        self.plotOfBPS9Voltage=_plot.plot();
        self.plotOfBPS9Voltage.yLim=[-200,200] # default [-200,200]
        self.plotOfBPS9Voltage.yLabel='V'
        self.plotOfBPS9Voltage.xLabel='Time [ms]'
        self.plotOfBPS9Voltage.subtitle='BP Voltage'
        self.plotOfBPS9Voltage.yData.append(self.bps9Voltage)
        self.plotOfBPS9Voltage.xData.append(self.time*1000)
        self.plotOfBPS9Voltage.yLegendLabel.append('BPS9')
        
        # initialize BPS9 current plot
        self.plotOfBPS9Current=_plot.plot();
        self.plotOfBPS9Current.yLim=[-5,60]  # default [-10,70]
        self.plotOfBPS9Current.yLabel='A'
        self.plotOfBPS9Current.subtitle='BP Current'
        self.plotOfBPS9Current.xLabel='Time [ms]'
        self.plotOfBPS9Current.yData.append(self.bps9Current)
        self.plotOfBPS9Current.xData.append(self.time*1000)
        self.plotOfBPS9Current.yLegendLabel.append('BPS9')
        
        # TODO(john) also make plots for BPS5 only

        if plot==True:
            self.plot()

    def plot(self):
        """ Plot all relevant plots """
        _plot.subPlot([self.plotOfVoltage,self.plotOfCurrent,self.plotOfGPUVoltageRequest])
        
    
class tpData:
    """
    Triple probe data
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
    tStart : float
        time (in seconds) to drim data after
    plot : bool
        plots all relevant plots if true
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
    plotOfISat : _plotTools.plot
        plot function of ion saturation current
    plotOfTipC : _plotTools.plot
        plot function of tip C voltages
    plotOfTipB : _plotTools.plot
        plot function of tip B voltages
    plotOfTipA : _plotTools.plot
        plot function of tip A voltages
    plotOfVf : _plotTools.plot
        plot function of floating voltages
    plotOfNe : _plotTools.plot
        plot function of density
    plotOfKTe : _plotTools.plot
        plot function of temperature
        
    Notes
    -----
    - I am not using the same time array for the section 5 or the section 8 
    triple probes.  I do this because different data acq. systems (e.g. north
    rack CPCI vs south rack CPCI) doesn't always return the EXACT same array.  
    I've run into issues where the length of the time arrays weren't the same
    length which causes problems during plotting. 
    - TPS2 was moved to section 5 (now called TPS5) during the summer of 2017.
    This may cause variable naming issues.  Be warned.  
    
    TODO The other cases for shotnumbers need to finalized so that legacy data
    can still be loaded.  
    """
    
    def __init__(self,shotno=95996,tStart=None,tStop=None,plot=False,probes='both'):  #sectionNum=2,
        
        self.shotno = shotno
        self.title = 'shotno = %s, triple probe data' % shotno
#        self.vTipA = None  # negative tip? TODO check
#        self.vTipB = None  # positive tip? TODO check
#        self.vTipC = None  # floating tip
        
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
        
        # initialize temperature plot
        self.plotOfKTe=_plot.plot();
        self.plotOfKTe.yLabel='eV'
        self.plotOfKTe.subtitle='Electron Temperature'
        self.plotOfKTe.title=self.title
        self.plotOfKTe.yLim=[-50, 100]       
        if probes=='both' or probes=='tps5': 
            self.plotOfKTe.yData.append(self.tps5Temp)
            self.plotOfKTe.xData.append(self.tps5Time*1000)
            self.plotOfKTe.yLegendLabel.append('TPS5')
        if probes=='both' or probes=='tps8':
            self.plotOfKTe.yData.append(self.tps8Temp)
            self.plotOfKTe.xData.append(self.tps8Time*1000)
            self.plotOfKTe.yLegendLabel.append('TPS8')
        
        # initialize density plot
        self.plotOfNe=_plot.plot();
        self.plotOfNe.yLabel=r'$m^{-3}$ $10^{18}$'
        self.plotOfNe.subtitle='Density'
        self.plotOfNe.yLim=[-1, 4.5]
        if probes=='both' or probes=='tps5':
            self.plotOfNe.yData.append(self.tps5Density/1e18)
            self.plotOfNe.xData.append(self.tps5Time*1000)
            self.plotOfNe.yLegendLabel.append('TPS5')
        if probes=='both' or probes=='tps8':
            self.plotOfNe.yData.append(self.tps8Density/1e18)
            self.plotOfNe.xData.append(self.tps8Time*1000)
            self.plotOfNe.yLegendLabel.append('TPS8')
            
        # initialize floating potential plot
        self.plotOfVf=_plot.plot();
        self.plotOfVf.yLabel='V'
        self.plotOfVf.subtitle='Floating Potential'
        self.plotOfVf.xLabel='time [ms]'
        self.plotOfVf.yLim=[-150, 75]
        if probes=='both' or probes=='tps5':
            self.plotOfVf.yData.append(self.tps5VFloat)
            self.plotOfVf.xData.append(self.tps5Time*1000)
            self.plotOfVf.yLegendLabel.append('TPS5')
        if probes=='both' or probes=='tps8':
            self.plotOfVf.yData.append(self.tps8VFloat)
            self.plotOfVf.xData.append(self.tps8Time*1000)
            self.plotOfVf.yLegendLabel.append('TPS8')
            
        # initialize tip A potential plot
        self.plotOfTipA=_plot.plot();
        self.plotOfTipA.yLabel='V'
        self.plotOfTipA.subtitle='Tip A Potential'
        self.plotOfTipA.xLabel='time [ms]'
#        self.plotOfTipA.yLim=[-100, 100]
        if probes=='both' or probes=='tps5':
            self.plotOfTipA.yData.append(self.tps5TipA)
            self.plotOfTipA.xData.append(self.tps5Time*1000)
            self.plotOfTipA.yLegendLabel.append('TPS5')
        if probes=='both' or probes=='tps8':
            self.plotOfTipA.yData.append(self.tps8TipA)
            self.plotOfTipA.xData.append(self.tps8Time*1000)
            self.plotOfTipA.yLegendLabel.append('TPS8')
            
        # initialize tip B potential plot
        self.plotOfTipB=_plot.plot();
        self.plotOfTipB.yLabel='V'
        self.plotOfTipB.subtitle='Tip B Potential'
        self.plotOfTipB.xLabel='time [ms]'
#        self.plotOfTipA.yLim=[-100, 100]
        if probes=='both' or probes=='tps5':
            self.plotOfTipB.yData.append(self.tps5TipB)
            self.plotOfTipB.xData.append(self.tps5Time*1000)
            self.plotOfTipB.yLegendLabel.append('TPS5')
        if probes=='both' or probes=='tps8':
            self.plotOfTipB.yData.append(self.tps8TipB)
            self.plotOfTipB.xData.append(self.tps8Time*1000)
            self.plotOfTipB.yLegendLabel.append('TPS8')
            
        # initialize tip C potential plot
        self.plotOfTipC=_plot.plot();
        self.plotOfTipC.yLabel='V'
        self.plotOfTipC.subtitle='Tip C Potential'
        self.plotOfTipC.xLabel='time [ms]'
#        self.plotOfTipA.yLim=[-100, 100]
        if probes=='both' or probes=='tps5':
            self.plotOfTipC.yData.append(self.tps5TipC)
            self.plotOfTipC.xData.append(self.tps5Time*1000)
            self.plotOfTipC.yLegendLabel.append('TPS5')
        if probes=='both' or probes=='tps8':
            self.plotOfTipC.yData.append(self.tps8TipC)
            self.plotOfTipC.xData.append(self.tps8Time*1000)
            self.plotOfTipC.yLegendLabel.append('TPS8')     
            
        # initialize ion sat current
        self.plotOfISat=_plot.plot();
        self.plotOfISat.yLabel='A'
        self.plotOfISat.subtitle='Ion Sat. Current'
        self.plotOfISat.xLabel='time [ms]'
#        self.plotOfTipA.yLim=[-100, 100]
        if probes=='both' or probes=='tps5':
            self.plotOfISat.yData.append(self.tps5Current)
            self.plotOfISat.xData.append(self.tps5Time*1000)
            self.plotOfISat.yLegendLabel.append('TPS5')
        if probes=='both' or probes=='tps8':
            self.plotOfISat.yData.append(self.tps8Current)
            self.plotOfISat.xData.append(self.tps8Time*1000)
            self.plotOfISat.yLegendLabel.append('TPS8')         
            
        if plot==True:
            self.plotRaw();
            self.plotProcessed();
            
    def plotRaw(self):
        _plot.subPlot([self.plotOfTipA,self.plotOfTipB,self.plotOfTipC,
                         self.plotOfISat]);            
            
    def plotProcessed(self):
        _plot.subPlot([self.plotOfKTe,self.plotOfNe,self.plotOfVf]);

    
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
    tStart : float
        time (in seconds) to drim data after
    plotSample : bool
        plots a single PA1 and PA2 sensor data
    plotAll : bool
        plots all 64 sensors.  Warning: this is hard on system memory
    smoothingAlgorithm : str
        informs function as to which smoothing algorithm to use on each PA 
        sensor
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    TODO add all remaining attributes here

    """
    
    def __init__(self,shotno=95540,tStart=None,tStop=None,plotSample=False, plotAll=False,
                 smoothingAlgorithm='tripleBoxCar'):
        self.shotno = shotno
        
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
                temp=_copy(_process.boxCar(data=self.pa1Raw[i][:],c=50))
                temp=_copy(_process.boxCar(data=temp,c=50))
                temp=_copy(_process.boxCar(data=temp,c=10))
                self.pa1RawFit.append(temp)
                self.pa1Data.append(self.pa1Raw[i]-temp)
                
                temp=_copy(_process.boxCar(data=self.pa2Raw[i][:],c=50))
                temp=_copy(_process.boxCar(data=temp,c=50))
                temp=_copy(_process.boxCar(data=temp,c=10))
                self.pa2RawFit.append(temp)
                self.pa2Data.append(self.pa2Raw[i]-temp)
        else:
            _sys.exit("You must specify a correct smoothing algorithm.  Exiting code...")
        
        # plot 
        if plotSample==True:
            self.plotPA1();
            self.plotPA2();
        if plotAll==True:
            self.plotAll()
            

    def plotPA1(self, i=0, alsoPlotRawAndFit=True):
        """ Plot one of the PA1 plots.  based on index, i. """
        p1=_plot.plot();
#        p1.xLim=[self.tStart,self.tStop]
        p1.xLabel='time [ms]'
        p1.yLabel=r'dB [G]'
        p1.title=str(self.shotno)+'. '+str(self.namesPA1[i])+' data'
        
        # smoothed data
        p1.yData.append(self.pa1Data[i]); #[self.taPDataRaw[i]]#
        p1.xData.append(self.pa1Time*1000);
        p1.yLegendLabel.append('smoothed')

        if alsoPlotRawAndFit==True:
            # raw data
            p1.yData.append(self.pa1Raw[i])
            p1.xData.append(self.pa1Time*1000);
            p1.yLegendLabel.append('raw')
            
            # fit data (which is subtracted from raw)
            p1.yData.append(self.pa1RawFit[i])
            p1.xData.append(self.pa1Time*1000);
            p1.yLegendLabel.append('fit')

        # plot
        p1.plot()
        
    def plotPA2(self, i=0, alsoPlotRawAndFit=True):
        """ Plot one of the PA2 plots.  based on index, i. """
        p1=_plot.plot();
#        p1.xLim=[self.tStart,self.tStop]
        p1.xLabel='time [ms]'
        p1.yLabel=r'dB [G]'
        p1.title=str(self.shotno)+'. '+str(self.namesPA2[i])+' data'
        
        # smoothed data
        p1.yData.append(self.pa2Data[i]); #[self.taPDataRaw[i]]#
        p1.xData.append(self.pa2Time*1000);
        p1.yLegendLabel.append('smoothed')

        if alsoPlotRawAndFit==True:
            # raw data
            p1.yData.append(self.pa2Raw[i])
            p1.xData.append(self.pa2Time*1000);
            p1.yLegendLabel.append('raw')
            
            # fit data (which is subtracted from raw)
            p1.yData.append(self.pa2RawFit[i])
            p1.xData.append(self.pa2Time*1000);
            p1.yLegendLabel.append('fit')

        # plot
        p1.plot()
        
    def plotAll(self,alsoPlotRawAndFit=True):
        """
        Plots poloidal sensor data for all 64 sensors
        
        Warning, 64 plots is tough on memory.  
        """
        for k in range(0,32):
            self.plotPA1(k,alsoPlotRawAndFit);
            self.plotPA2(k,alsoPlotRawAndFit);
            
    
class fbData:
    """
    Downloads feedback (FB) array sensor data.   
    
    Parameters
    ----------
    shotno : int
        shot number of desired data
    tStart : float
        time (in seconds) to trim data before
    tStart : float
        time (in seconds) to drim data after
    plotSample : bool
        plots results from only one sensor.  both poloidal and radial data
    plotAll : bool
        plots all 40 sensors.  both poloidal and radial data
    smoothingAlgorithm : str
        informs function as to which smoothing algorithm to use on each PA 
        sensor
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    phi : numpy.ndarray
        array of toroidal locations for all sensors.  units in degrees.
    TODO add all remaining attributes here
        
    """
    def __init__(self,shotno=95540,tStart=None,tStop=None,plotSample=False, plotAll=False,smoothingAlgorithm='tripleBoxCar'):
        self.shotno = shotno

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
                    temp=_copy(_process.boxCar(data=self.fbPolRaw[j][i][:],c=50))
                    temp=_copy(_process.boxCar(data=temp,c=50))
                    temp=_copy(_process.boxCar(data=temp,c=10))
                    self.fbPolRawFit[j].append(temp)
                    self.fbPolData[j].append(self.fbPolRaw[j][i]-temp)
                    
                    temp=_copy(_process.boxCar(data=self.fbRadRaw[j][i][:],c=50))
                    temp=_copy(_process.boxCar(data=temp,c=50))
                    temp=_copy(_process.boxCar(data=temp,c=10))
                    self.fbRadRawFit[j].append(temp)
                    self.fbRadData[j].append(self.fbRadRaw[j][i]-temp)
        else:
            _sys.exit("You must specify a correct smoothing algorithm.  Exiting code...")
     
        # plot
        if plotSample==True:
            self.plotPol();
            self.plotRad();
        if plotAll==True:
            self.plotAll()
            

    def plotPol(self, row=0, col=0,alsoPlotRawAndFit=True):
        """
        Plots poloidal data from FB sensors
        """
        i=col; j=row;
        
        # initialize plot
        p1=_plot.plot();
        p1.xLabel='time [ms]'
        p1.yLabel=r'dB [G]'
        p1.title=str(self.shotno)+'. '+str(self.fbPolNames[i][j])+' data'
        
        # smoothed data
        p1.yData.append(self.fbPolData[i][j]); #[self.taPDataRaw[i]]#
        p1.xData.append(self.fbPolTime*1000);
        p1.yLegendLabel.append('Smoothed')
        
        if alsoPlotRawAndFit==True:
            # raw data
            p1.yData.append(self.fbPolRaw[i][j])
            p1.xData.append(self.fbPolTime*1000);
            p1.yLegendLabel.append('Raw')
            
            # fit data (which is subtracted from raw)
            p1.yData.append(self.fbPolRawFit[i][j])
            p1.xData.append(self.fbPolTime*1000);
            p1.yLegendLabel.append('Fit')
            
        # plot
        p1.plot()
        
        
    def plotRad(self, row=0, col=0,alsoPlotRawAndFit=True):
        """
        Plots radial data from FB sensors
        """
        i=col; j=row;
        
        # initialize plot
        p1=_plot.plot();
        p1.xLabel='time [ms]'
        p1.yLabel=r'dB [G]'
        p1.title=str(self.shotno)+'. '+str(self.fbRadNames[i][j])+' data'
        
        # smoothed data
        p1.yData.append(self.fbRadData[i][j]); #[self.taPDataRaw[i]]#
        p1.xData.append(self.fbRadTime*1000);
        p1.yLegendLabel.append('Smoothed')
        
        if alsoPlotRawAndFit==True:
            # raw data
            p1.yData.append(self.fbRadRaw[i][j])
            p1.xData.append(self.fbRadTime*1000);
            p1.yLegendLabel.append('Raw')
            
            # fit data (which is subtracted from raw)
            p1.yData.append(self.fbRadRawFit[i][j])
            p1.xData.append(self.fbRadTime*1000);
            p1.yLegendLabel.append('Fit')
            
        # plot
        p1.plot()
        
    def plotAll(self,alsoPlotRawAndFit=True):
        """
        Plots poloidal and radial sensor data for all 40 sensors (80 in total)
        """
        for k in range(0,10):
            for l in range(0,4):
                self.plotPol(row=k,col=l,alsoPlotRawAndFit=alsoPlotRawAndFit);
                # self.plotRad(row=k,col=l,alsoPlotRawAndFit=alsoPlotRawAndFit);
    
    
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
    tStart : float
        time (in seconds) to drim data after
    plotSample : bool
        plots one PA1 and one PA2 sensor data
    plotAll : bool
        plots all 30 poloidal sensors
    smoothingAlgorithm : str
        informs function as to which smoothing algorithm to use on each PA 
        sensor
        
    Attributes
    ----------
    shotno : int
        shot number of desired data
    TODO add all remaining attributes here

    """
        
    def __init__(self,shotno=95540,tStart=1,tStop=8,plotSample=False,plotAll=False,smoothingAlgorithm='tripleBoxCar'):
        self.shotno = shotno
#        self.taPDataRaw = []
#        self.taPData = []
#        self.taRDataRaw = []
#        self.time = None
#        self.tStart = tStart
#        self.tStop = tStop
        
        # names of poloidal and radial sensors
        self.namesTAPol=['TA01_S1P', 'TA01_S2P', 'TA01_S3P', 'TA02_S1P', 'TA02_S2P', 'TA02_S3P', 'TA03_S1P', 'TA03_S2P', 'TA03_S3P', 'TA04_S1P', 'TA04_S2P', 'TA04_S3P', 'TA05_S1P', 'TA05_S2P', 'TA05_S3P', 'TA06_S1P', 'TA06_S2P', 'TA06_S3P', 'TA07_S1P', 'TA07_S2P', 'TA07_S3P', 'TA08_S1P', 'TA08_S2P', 'TA08_S3P', 'TA09_S1P', 'TA09_S2P', 'TA09_S3P', 'TA10_S1P', 'TA10_S2P', 'TA10_S3P'];
        self.namesTARad=['TA01_S2R', 'TA02_S2R', 'TA03_S2R', 'TA04_S2R', 'TA05_S2R', 'TA06_S2R', 'TA07_S2R', 'TA08_S2R', 'TA09_S2R', 'TA10_S2R']
#        self.phi=_np.zeros(30);
#        self.phiRc=_np.zeros(10);
        
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
                temp=_copy(_process.boxCar(data=self.taPolRaw[i][:],c=50))
                temp=_copy(_process.boxCar(data=temp,c=50))
                temp=_copy(_process.boxCar(data=temp,c=10))
                self.taPolRawFit.append(temp)
                self.taPolData.append(self.taPolRaw[i]-temp)
                
                if i < 10:
                    temp=_copy(_process.boxCar(data=self.taRadRaw[i][:],c=50))
                    temp=_copy(_process.boxCar(data=temp,c=50))
                    temp=_copy(_process.boxCar(data=temp,c=10))
                    self.taRadRawFit.append(temp)
                    self.taRadData.append(self.taRadRaw[i]-temp)
        else:
            _sys.exit("You must specify a correct smoothing algorithm.  Exiting code...")
                     
    

        # plot
        if plotSample==True:
            self.plotPol();
        if plotAll==True:
            self.plotAll();
        
    def plotPol(self, i=0, alsoPlotRawAndFit=True):
        """ Plot one of the PA1 plots.  based on index, i. """
        p1=_plot.plot();
#        p1.xLim=[self.tStart,self.tStop]
        p1.xLabel='time [ms]'
        p1.yLabel=r'Gauss'
        p1.title=str(self.shotno)+'. '+str(self.namesTAPol[i])+' data'
        
        # smoothed data
        p1.yData.append(self.taPolData[i]); #[self.taPDataRaw[i]]#
        p1.xData.append(self.taPolTime*1000);
        p1.yLegendLabel.append('smoothed')

        if alsoPlotRawAndFit==True:
            # raw data
            p1.yData.append(self.taPolRaw[i])
            p1.xData.append(self.taPolTime*1000);
            p1.yLegendLabel.append('raw')
            
            # fit data (which is subtracted from raw)
            p1.yData.append(self.taPolRawFit[i])
            p1.xData.append(self.taPolTime*1000);
            p1.yLegendLabel.append('fit')

        # plot
        p1.plot()
        
    def plotAll(self,alsoPlotRawAndFit=True):
        """
        Plots poloidal sensor data for all 64 sensors
        
        Warning, 64 plots is tough on memory.  
        """
        for k in range(0,30):
            self.plotPol(k,alsoPlotRawAndFit);
  

class jumperData:
    """
    Measures jumper data between toroidal sections    
    
    work in progress
    """

            
class gpuControlData:
    """
    reads control data from Caliban's control files
    """
    def __init__(self,shotno=96496, tStart=1, tStop=8, plot=False, password='', forceDownload=False, calibrateTime=False, plotEVERYTHING=False):
        # note, shootno=None loads the data generated from FAKE_INPUT        
        import _feedBackTools as fbt
        import os     

        self.shotno=shotno;
        self._tStart=tStart
        self._tStop=tStop
        
        if _ON_HBTEP_SERVER==False:
            # ssh data must be transfered if operating remotely.
            # this copies control data to your local computer.
        
            # change data directories because we are operating remotely
            fbt._CONTROL_CODE_PATH=_pref._LOCAL_DATA_DIR
            fbt._DATA_PATH=_pref._LOCAL_DATA_DIR
        
            # check to see if the file has previously been downloaded to local directory.  if not, it is downloaded
            if shotno!=None:
                filePath=_pref._LOCAL_DATA_DIR + "fbsettings_" +str(int(shotno))+'.py';
            else:
                filePath=_pref._LOCAL_DATA_DIR + 'fbsettings.py';
            print filePath
            if os.path.isfile(filePath)==False or forceDownload==True:
                self._downloadCDFromCaliban(password=password)
            
        ## Load time from control files.  
        # TODO(JOHN) recalibrate time offset
        timeOffset = 0.0340/1000
        self.time=fbt.get_ctrl_times(shotno)+timeOffset; ## NOTE:  For some reason, there is a shift in the time data on the GPU vs. the Tree.  Not sure why.  I "roughly" correct for it in the next line.  Calibrated with hbt.compareBP_Bn1_GPU(97088)

        # get data
        analogIn=fbt.get_ctrl_ai(shotno); # TODO(JOHN) remove _numCh 
#        try:
#            analogInRaw=fbt.get_ctrl_airaw(shotno);
#        except Exception:
#            print("raw file not present.  skipping...")     
#            pass
        analogOut=fbt.get_ctrl_ao(shotno);
        mAmp=fbt.get_ctrl_mamp(shotno);
        mPhase=fbt.get_ctrl_mphase(shotno);
        mFreq=fbt.get_ctrl_mfreq(shotno);
            
        # enforce that units are in seconds (NOT milliseconds)
        if tStop > 1:
            tStop=tStop*1e-3
            tStart=tStart*1e-3
            
        # trim time and data
        temp,analogOut=_trimTime(self.time,list(analogOut),tStart,tStop)
        temp,analogIn=_trimTime(self.time,list(analogIn),tStart,tStop)
        temp,mAmp=_trimTime(self.time,list(mAmp),tStart,tStop)
        temp,mPhase=_trimTime(self.time,list(mPhase),tStart,tStop)
        temp,mFreq=_trimTime(self.time,list(mFreq),tStart,tStop)
        self.time=temp
        
        # distribute trimmed data to class variables
        self.BPS9Voltage=analogOut[43]; #43
        self.mAmpCosSec1=mAmp[0];
        self.mAmpSinSec1=mAmp[1];
        self.mAmpSec1=_np.sqrt(self.mAmpCosSec1**2 + self.mAmpSinSec1**2)
        self.mAmpCosSec2=mAmp[2];
        self.mAmpSinSec2=mAmp[3];
        self.mAmpSec2=_np.sqrt(self.mAmpCosSec2**2 + self.mAmpSinSec2**2)
        self.mAmpCosSec3=mAmp[4];
        self.mAmpSinSec3=mAmp[5];
        self.mAmpSec3=_np.sqrt(self.mAmpCosSec3**2 + self.mAmpSinSec3**2)
        self.mAmpCosSec4=mAmp[6];
        self.mAmpSinSec4=mAmp[7];
        self.mAmpSec4=_np.sqrt(self.mAmpCosSec4**2 + self.mAmpSinSec4**2)
        self.mPhaseSec1=mPhase[0];
        self.mPhaseSec2=mPhase[2];
        self.mPhaseSec3=mPhase[4];
        self.mPhaseSec4=mPhase[6];
        self.mFreqSec1=mFreq[0];
        self.mFreqSec2=mFreq[2];
        self.mFreqSec3=mFreq[4];
        self.mFreqSec4=mFreq[6];

        # initialize BPS9 voltage plot
        self.plotOfBPS9Voltage=_plot.plot();
        self.plotOfBPS9Voltage.yLim=[-11,11]
        self.plotOfBPS9Voltage.yLabel='V'
        self.plotOfBPS9Voltage.xLabel='Time [ms]'
        self.plotOfBPS9Voltage.subtitle='BP Voltage, GPU Output'
        self.plotOfBPS9Voltage.title=str(self.shotno);
        self.plotOfBPS9Voltage.yData.append(self.BPS9Voltage)
        self.plotOfBPS9Voltage.xData.append(self.time*1000)
        self.plotOfBPS9Voltage.yLegendLabel.append('BPS9')
        
        # init mode amp plot
        self.plotOfAmplitudes=_plot.plot();
        self.plotOfAmplitudes.yLabel='G'
        self.plotOfAmplitudes.xLabel='Time [ms]'
        self.plotOfAmplitudes.subtitle='n=1 mode amplitude'
        self.plotOfAmplitudes.title=str(self.shotno);
        self.plotOfAmplitudes.yData.append(self.mAmpSec1*1e4)
        self.plotOfAmplitudes.xData.append(self.time*1000)
        self.plotOfAmplitudes.yLegendLabel.append('GPU-Sec1')
        if plotEVERYTHING==True:
            self.plotOfAmplitudes.yData.append(self.mAmpSec2*1e4)
            self.plotOfAmplitudes.yData.append(self.mAmpSec3*1e4)
            self.plotOfAmplitudes.yData.append(self.mAmpSec4*1e4)
            self.plotOfAmplitudes.xData.append(self.time*1000)
            self.plotOfAmplitudes.xData.append(self.time*1000)
            self.plotOfAmplitudes.xData.append(self.time*1000)
            self.plotOfAmplitudes.yLegendLabel.append('FBSensors-GPU Sec2')
            self.plotOfAmplitudes.yLegendLabel.append('FBSensors-GPU Sec3')
            self.plotOfAmplitudes.yLegendLabel.append('FBSensors-GPU Sec4')
         
        # init mode phase plot
        self.plotOfPhase=_plot.plot();
        self.plotOfPhase.yLabel='Radians'
        self.plotOfPhase.xLabel='Time [ms]'
        self.plotOfPhase.subtitle='n=1 mode phase'
        self.plotOfPhase.title=str(self.shotno);
        self.plotOfPhase.yLim=[-_np.pi, _np.pi]  
        self.plotOfPhase.yData.append(_process.wrapPhase(self.mPhaseSec1))
        self.plotOfPhase.linestyle.append('')
        self.plotOfPhase.marker.append('.')
        self.plotOfPhase.xData.append(self.time*1000)
        self.plotOfPhase.yLegendLabel.append('FBSensors-GPU Sec1')
        if plotEVERYTHING==True:
            self.plotOfPhase.yData.append(_process.wrapPhase(self.mPhaseSec2))
            self.plotOfPhase.yData.append(_process.wrapPhase(self.mPhaseSec3))
            self.plotOfPhase.yData.append(_process.wrapPhase(self.mPhaseSec4))
            self.plotOfPhase.xData.append(self.time*1000)
            self.plotOfPhase.xData.append(self.time*1000)
            self.plotOfPhase.xData.append(self.time*1000)
            self.plotOfPhase.yLegendLabel.append('FBSensors-GPU Sec2')
            self.plotOfPhase.yLegendLabel.append('FBSensors-GPU Sec3')
            self.plotOfPhase.yLegendLabel.append('FBSensors-GPU Sec4')
            self.plotOfPhase.linestyle.append('')
            self.plotOfPhase.marker.append('.')
            self.plotOfPhase.linestyle.append('')
            self.plotOfPhase.marker.append('.')
            self.plotOfPhase.linestyle.append('')
            self.plotOfPhase.marker.append('.')

        # init mode frequency plot
        self.plotOfFreq=_plot.plot();
        self.plotOfFreq.yLabel='kHz'
        self.plotOfFreq.xLabel='Time [ms]'
        self.plotOfFreq.subtitle='n=1 mode freq'
        self.plotOfFreq.title=str(self.shotno);
        self.plotOfFreq.yData.append(self.mFreqSec1*1e-3)
        self.plotOfFreq.xData.append(self.time*1000)
        self.plotOfFreq.yLegendLabel.append('FBSensors-GPU Sec1')
        if plotEVERYTHING==True:
            self.plotOfFreq.yData.append(self.mFreqSec2*1e-3)
            self.plotOfFreq.yData.append(self.mFreqSec3*1e-3)
            self.plotOfFreq.yData.append(self.mFreqSec4*1e-3)
            self.plotOfFreq.xData.append(self.time*1000)
            self.plotOfFreq.xData.append(self.time*1000)
            self.plotOfFreq.xData.append(self.time*1000)
            self.plotOfFreq.yLegendLabel.append('FBSensors-GPU Sec2')
            self.plotOfFreq.yLegendLabel.append('FBSensors-GPU Sec3')
            self.plotOfFreq.yLegendLabel.append('FBSensors-GPU Sec4')

        if plot == True:
            _plot.subPlot([self.plotOfAmplitudes,self.plotOfPhase,self.plotOfFreq,self.plotOfBPS9Voltage])

        if calibrateTime==True:
            self.calibrateCPCIAndGPUData()

    def calibrateCPCIAndGPUData(self):
        """
        This plots
        1)  Pre-amp Voltage out from the GPU going to BPS9 (CPCI)
        2)  BPS9 volltage (CPCI)
        3)  BPS9 requested voltage sent out to BPS9 (GPU)
        The goal is that by plotting all of these side by side, it will be easier 
        to make sure that their times all line up
        """
        # enforce that units are in seconds (NOT milliseconds)
        tStop=self._tStop
        tStart=self._tStart
        if tStop > 1:
            tStop=tStop*1e-3
            tStart=tStart*1e-3
        BP=bpData(self.shotno,tStart,tStop)
        _plot.subPlot([BP.plotOfGPUVoltageRequest,BP.plotOfBPS9Voltage,self.plotOfBPS9Voltage])
            
    def _downloadCDFromCaliban(self,password=''):
        """
        Downloads control data files from Caliban
        """
        
        if password=='' or password == None:
            #password = raw_input("Enter spitzer password:  ")
            password=_rwData.getPwd(system=_pref._HBT_SERVER_NAME,username=_pref._HBT_SERVER_USERNAME);

        # address where all real data is stored (i.e. each shot has a shot number)
#        _REMOTE_DATA_PATH='/opt/hbt/data/control'
        # address where code on caliban in located.  this is where the FAKE_INPUT shots are stored
        REMOTE_CODE_PATH='/home/brooks/TokaMac/control' # this will need to be adjusted for every user
        
        # open connection
        sshCon = _rwData.scpData(password=password,port=22,username=_pref._HBT_SERVER_USERNAME,address=_pref._HBT_SERVER_ADDRESS)

        # _copy data
        if self.shotno=='' or self.shotno==None:
            sshCon.downloadFile('%s/ao_store.dat' % (REMOTE_CODE_PATH), localFilePath=_FILEDIR+'ao_store.dat' )
            sshCon.downloadFile('%s/ai_store.dat' % (REMOTE_CODE_PATH), localFilePath=_FILEDIR+'ai_store.dat' )
            try:
                sshCon.downloadFile('%s/airaw_store.dat' % (REMOTE_CODE_PATH))
            except:
                print("raw file not present.  skipping...")
                pass
            sshCon.downloadFile('%s/mamp_store.dat' % (REMOTE_CODE_PATH), localFilePath=_FILEDIR+'mamp_store.dat' )
            sshCon.downloadFile('%s/mphase_store.dat' % (REMOTE_CODE_PATH), localFilePath=_FILEDIR+'mphase_store.dat' )
            sshCon.downloadFile('%s/mfreq_store.dat' % (REMOTE_CODE_PATH), localFilePath=_FILEDIR+'mfreq_store.dat' )
            sshCon.downloadFile('%s/fb_store.dat' % (REMOTE_CODE_PATH), localFilePath=_FILEDIR+'fb_store.dat' )
            sshCon.downloadFile('%s/fbsettings.py' % (REMOTE_CODE_PATH), localFilePath=_FILEDIR+'fbsettings.py' )
        else:
            sshCon.downloadFile('%s/ao_store_%d.dat' % (_REMOTE_DATA_PATH, self.shotno), localFilePath='%s/ao_store_%d.dat' % (_FILEDIR, self.shotno))
            sshCon.downloadFile('%s/ai_store_%d.dat' % (_REMOTE_DATA_PATH, self.shotno), localFilePath='%s/ai_store_%d.dat' % (_FILEDIR, self.shotno))
            sshCon.downloadFile('%s/mamp_store_%d.dat' % (_REMOTE_DATA_PATH, self.shotno), localFilePath='%s/mamp_store_%d.dat' % (_FILEDIR, self.shotno))
            sshCon.downloadFile('%s/mphase_store_%d.dat' % (_REMOTE_DATA_PATH, self.shotno), localFilePath='%s/mphase_store_%d.dat' % (_FILEDIR, self.shotno))
            sshCon.downloadFile('%s/mfreq_store_%d.dat' % (_REMOTE_DATA_PATH, self.shotno), localFilePath='%s/mfreq_store_%d.dat' % (_FILEDIR, self.shotno))
            sshCon.downloadFile('%s/fbsettings_%d.py' % (_REMOTE_DATA_PATH, self.shotno), localFilePath='%s/fbsettings_%d.py'  % (_FILEDIR, self.shotno))
            print '%s/fb_store_%d.dat' % (_REMOTE_DATA_PATH, self.shotno)
            try:
                sshCon.downloadFile('%s/fb_store_%d.dat' % (_REMOTE_DATA_PATH, self.shotno), localFilePath='%s/fb_store_%d.dat'  % (_FILEDIR, self.shotno))
            except Exception:
                print("fb_store file not present.  skipping...")     
                pass

        # close connection
        sshCon.closeConnection();
        
        
    def GPUBPData(self,cpciShotno=None,operatingMode='currentControl'):
        """
        Load fb_store_<shotno>.dat from GPU.  This is the active feedback data file
        
        Notes
        -----
        Sometimes this function loads a fb_store.dat file that was generated 
        from fake data (FAKE_INPUT) instead of from a real plasma shot.  In 
        those cases, the CPCI data used in the plots below should be pulled
        from the same shotno used to create the FAKE_INPUT.  This is what the
        variable, cpciShotno, does.  
        
        note: 97598 is one of first GPU/CPSI combined shots
        """
        class fb:
            """ container object for feedback data """
        self.fb=fb;
        
        import _feedBackTools as _fbt
        DT=0.000006; # should be 0.00006
        timeOffset=-1.1;# calibrated on shotno 97601.  sig gen on input with 200 Hz triangular wave.  
        
        GPUBP=_fbt.get_fb(self.shotno);
        
        # the time array is mostly zeros.  removing these from time AND from the data
        indices=_np.where(GPUBP[0,:]!=0)[0]  
        self.fb.time=GPUBP[0,indices]*DT+timeOffset*1e-3;  # multipling by DT converts sequential numbers into microseconds.  

        # set feedback variables from loaded data
        print _np.shape(GPUBP)
        self.fb.Vact=GPUBP[1,indices]
        self.fb.Vreq=GPUBP[2,indices]
        self.fb.Iact=GPUBP[3,indices]
        self.fb.Ireq=GPUBP[4,indices]
        self.fb.fb=GPUBP[5,indices]
        self.fb.ff=GPUBP[6,indices]
        self.fb.P=GPUBP[7,indices]
        self.fb.I=GPUBP[8,indices]
        self.fb.D=GPUBP[9,indices]      
        self.fb.error=GPUBP[10,indices] 
        self.fb.wfb=GPUBP[11,indices]
        self.fb.wreq=GPUBP[12,indices]
        self.fb.wact=GPUBP[13,indices]
        
        # load supplementary data
        if self.shotno!=None and cpciShotno==None:
            cpciShotno=self.shotno
        bp=bpData(cpciShotno,self._tStart,self._tStop)
        nMode=nModeData(cpciShotno,self._tStart,self._tStop)
        
        # initiate BP voltage plot
        self.plotOfVin=_plot.plot();
        self.plotOfVin.yLabel='V'
        self.plotOfVin.xLabel='Time [ms]'
        self.plotOfVin.subtitle='BPS9 Voltages'
        self.plotOfVin.title=str(self.shotno);
        self.plotOfVin.yLim=[-120, 120]  
        self.plotOfVin.yData.append(self.fb.Vact) 
        self.plotOfVin.xData.append(self.fb.time*1e3)
        self.plotOfVin.yLegendLabel.append(r'GPU-V$_{measured}$')
        self.plotOfVin.yData.append(self.fb.Vreq) 
        self.plotOfVin.xData.append(self.fb.time*1e3)
        self.plotOfVin.yLegendLabel.append(r'GPU-V$_{request}$')
        self.plotOfVin.yData.append(bp.bps9Voltage) 
        self.plotOfVin.xData.append(bp.time*1e3)
        self.plotOfVin.yLegendLabel.append(r'CPCI-V$_{measured}$')
        self.plotOfVin.yData.append(bp.bps9GPURequestVoltage) 
        self.plotOfVin.xData.append(bp.time*1e3)
        self.plotOfVin.yLegendLabel.append(r'CPCI-V$_{request}$')
        if operatingMode=='voltageControl':
            self.plotOfIin.yData.append(self.fb.error) 
            self.plotOfIin.xData.append(self.fb.time*1e3)
            self.plotOfIin.yLegendLabel.append('GPU-error')
        
        # initialize BP current plot
        self.plotOfIin=_plot.plot();
        self.plotOfIin.yLabel='A'
        self.plotOfIin.xLabel='Time [ms]'
        self.plotOfIin.subtitle='BPS9 Currents'
        self.plotOfIin.title=str(self.shotno);
        self.plotOfIin.yLim=[-20, 50]  
        self.plotOfIin.yData.append(self.fb.Iact) 
        self.plotOfIin.xData.append(self.fb.time*1e3)
        self.plotOfIin.yLegendLabel.append('GPU-I$_{measured}$')
        self.plotOfIin.yData.append(self.fb.Ireq) 
        self.plotOfIin.xData.append(self.fb.time*1e3)
        self.plotOfIin.yLegendLabel.append('GPU-I$_{request}$')
        self.plotOfIin.yData.append(bp.bps9Current) 
        self.plotOfIin.xData.append(bp.time*1e3)
        self.plotOfIin.yLegendLabel.append('CPCI-I$_{measured}$')
        if operatingMode=='currentControl':
            self.plotOfIin.yData.append(self.fb.error) 
            self.plotOfIin.xData.append(self.fb.time*1e3)
            self.plotOfIin.yLegendLabel.append('GPU-error')
        
        # initialize plot of gains
        self.plotOfGains=_plot.plot();
        self.plotOfGains.yLabel=''
        self.plotOfGains.xLabel='Time [ms]'
        self.plotOfGains.subtitle='Gains'
        self.plotOfGains.title=str(self.shotno);
        self.plotOfGains.yLim=[-50, 100]  
        self.plotOfGains.yData.append(self.fb.fb) 
        self.plotOfGains.xData.append(self.fb.time*1e3)
        self.plotOfGains.yLegendLabel.append('fb=P+I')
        self.plotOfGains.yData.append(self.fb.ff) 
        self.plotOfGains.xData.append(self.fb.time*1e3)
        self.plotOfGains.yLegendLabel.append('ff')
        self.plotOfGains.yData.append(self.fb.P) 
        self.plotOfGains.xData.append(self.fb.time*1e3)
        self.plotOfGains.yLegendLabel.append('P')
        self.plotOfGains.yData.append(self.fb.I) 
        self.plotOfGains.xData.append(self.fb.time*1e3)
        self.plotOfGains.yLegendLabel.append('I')
        
        # initialize plot of freq
        self.plotOfFBFreq=_plot.plot();
        self.plotOfFBFreq.yLabel='kHz'
        self.plotOfFBFreq.xLabel='Time [ms]'
        self.plotOfFBFreq.subtitle='n=1 mode Frequency'
        self.plotOfFBFreq.title=str(self.shotno);
        self.plotOfFBFreq.yLim=[-10, 20]  
        self.plotOfFBFreq.yData.append(self.fb.wact*1e-3) 
        self.plotOfFBFreq.xData.append(self.fb.time*1e3)
        self.plotOfFBFreq.yLegendLabel.append(r'GPU-$\omega_{measured}$')
        self.plotOfFBFreq.yData.append(self.fb.wfb*1e-3) 
        self.plotOfFBFreq.xData.append(self.fb.time*1e3)
        self.plotOfFBFreq.yLegendLabel.append(r'GPU-$\omega_{fb}$')
        self.plotOfFBFreq.yData.append(self.fb.wreq*1e-3) 
        self.plotOfFBFreq.xData.append(self.fb.time*1e3)
        self.plotOfFBFreq.yLegendLabel.append(r'GPU-$\omega_{request}$')
        self.plotOfFBFreq.yData.append(nMode.n1FreqStrongFilter*1e-3) 
        self.plotOfFBFreq.xData.append(nMode.time*1e3)
        self.plotOfFBFreq.yLegendLabel.append(r'CPCI-$\omega_{measured}$')
        if operatingMode=='frequencyControl':
            self.plotOfFBFreq.yData.append(self.fb.error*1e-3) 
            self.plotOfFBFreq.xData.append(self.fb.time*1e3)
            self.plotOfFBFreq.yLegendLabel.append('GPU-error')

        # add to mode amplitude plot
        self.plotOfAmplitudes.yData.append(nMode.n1Amp) #*1e4
        self.plotOfAmplitudes.xData.append(nMode.time*1000)
        self.plotOfAmplitudes.yLegendLabel.append('CPCI')
        
        # add to mode phase plot
        self.plotOfPhase.yData.append(nMode.n1Phase)
        self.plotOfPhase.xData.append(nMode.time*1000)
        self.plotOfPhase.yLegendLabel.append('CPCI')
        
        # construct fb subplot
        self.plotOfFeedback=_plot.subPlot([self.plotOfVin,
                                           self.plotOfIin,
                                           self.plotOfGains,
                                           self.plotOfFBFreq,
                                           self.plotOfAmplitudes,
                                           self.plotOfPhase])


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
    
    
# TODO:There are more and better ways to write blacklist functions.  Write more


###############################################################################
### HBTEP shot data format

def loadShotData(shotno=96635,tStart=0.0,tStop=8.0,forceDownload=False,gpu=False,BP=True,TP=False):
    
    """
    Returns shotData structure/class
    
    If the data has been previously downloaded, the data is loaded locally. 
    Otherwise, the data is downloaded (and stored for future access) from Spitzer.
    """
    import os
    
    # check if shotno is a single entry or an array of entries
    if type(shotno)==_np.ndarray:
        m=len(shotno);
    else:
        m=1;
        shotno=_np.array([shotno])
        
    # init output list
    data=[None]*m
    
    # iterate through 1 or more shotnumbers
    for i in range(0,m):
                
        # filepath to save pickle file
        filePath=_FILEDIR + str(int(shotno[i]))+'.pickle'
        
        # check if file is alreay downloaded.  download permanently otherwise.  also downloads file if specifically requested.
        if os.path.isfile(filePath)==False or forceDownload==True:
            print "downloading "+str(int(shotno[i]))+" shot data from Spitzer"
            tags=['nMode','r','q','Ip'] # all HBT data will have this information
            # not all HBT data will have the following.  Check first before loading.
            if gpu==True:
                tags.append('gpu')
            if BP==True:
                tags.append('BP')
            if TP==True:
                tags.append('TP')

            # download shot data to picke file
            _shotData(shotno[i],save=True,tags=tags)
            
        # load file from pickle file
        with open(filePath) as f:  
            data[i] = _pk.load(f)
            
        # trim time to desired range
        if True:
            iStart=_process.find_nearest(data[i].time,tStart/1000.);
            iStop=_process.find_nearest(data[i].time,tStop/1000.)
            
            try:
                data[i].time=data[i].time[iStart:iStop];
                data[i].n1Amp=data[i].n1Amp[iStart:iStop];
                data[i].n1Phase=data[i].n1Phase[iStart:iStop];
                data[i].n1PhaseVelRaw=data[i].n1PhaseVelRaw[iStart:iStop];
                data[i].n1PhaseVelStrongFilter=data[i].n1PhaseVelStrongFilter[iStart:iStop];
                
            except AttributeError:
                print "No nMode data. Skipping"
                
            try:
                data[i].minorRadius=data[i].minorRadius[iStart:iStop];
                data[i].majorRadius=data[i].majorRadius[iStart:iStop];
            except AttributeError:
                print "No radial data. Skipping"
                                
            try:
                data[i].q=data[i].q[iStart:iStop];
                data[i].qUpsample=data[i].qUpsample[iStart:iStop];
            except AttributeError:
                print "No q data. Skipping"   
                            
            try:
                data[i].BPV=data[i].BPV[iStart:iStop];
                data[i].BPI=data[i].BPI[iStart:iStop];
            except AttributeError:
                print "No BP data. Skipping"
                
            try:
                data[i].Ip=data[i].Ip[iStart:iStop];
            except AttributeError:
                print "No Ip data. Skipping"
                
    # return data
    if m==1:
        return data[0]
    else:
        return data
    # TODO(John) i would like to combine the loadShotData function and the shotData class, but returning a loaded class is hard to do from that same class.
                
       
class _shotData:
    """
    data structure for storing multiple data types for a single shot number
    """
#    def __init__(self,shotno=None,tStart=1.0,tStop=8.0,tags=['BP','nMode','r','q'],save=True):
    def __init__(self,shotno=None,tags=['nMode','r','q','Ip'],save=True):

        self.shotno=shotno;
        self.tags=tags;
        tStartFixed=0.0;  # all data downloaded starts at this time
        tStopFixed=10.0;  # all data downloaded stops at this time
        
        import HBTTools as hbt
#        import loadDataTools as ldt

        if "nMode" in tags:
            nMode=hbt.nModeAnalysis(shotno=shotno,tStart=tStartFixed,tStop=tStopFixed)
            
            self.time=nMode.time;
            self.n1Amp=nMode.n1Amp;
            self.n1Phase=nMode.n1Phase;
            self.n1PhaseVelRaw=nMode.n1PhaseVelRaw;
            self.n1PhaseVelStrongFilter=nMode.n1PhaseVelStrongFilter;  
            self.plotOfAmps=nMode.plotOfAmps
            self.plotOfN1Amp=nMode.plotOfN1Amp
            self.plotOfN1Phase=nMode.plotOfN1Phase
            self.plotOfN1Freq=nMode.plotOfN1Freq
            self.plotOfN1PhaseAndAmp=nMode.plotOfPhaseAmp
           
        if "BP" in tags:
#            self.BP=ldt.BPData(shotno,tStart=tStartFixed,tStop=tStopFixed,loadMethod='new')
            BP=BPData(shotno,tStart=tStartFixed,tStop=tStopFixed,loadMethod='new')
            
            self.BPV=BP.voltageBPS9;
            self.BPI=BP.currentBPS9;
            self.plotOfBPV=BP.plotOfVoltage
            self.plotOfBPI=BP.plotOfCurrent
            
        if "r" in tags:
            r=hbt.plasmaRadius(shotno,tStart=tStartFixed,tStop=tStopFixed)
            
            self.minorRadius=r.minorRadius;
            self.majorRadius=r.majorRadius;
            self.plotOfMajorRadius=r.plotOfMajorRadius;
            self.plotOfMinorRadius=r.plotOfMinorRadius;
            
        if "q" in tags:
            q=hbt.qStar(shotno,tStart=tStartFixed,tStop=tStopFixed)
            
            self.q=q.q;
            self.qUpsample=q.qUpsample;
            self.plotOfQ=q.plotOfQ
            
        if "gpu" in tags:
            gpu=gpuControlData(shotno,tStart=tStartFixed,tStop=tStopFixed, password='',download=True)

            self.gpuBPV=gpu.BPS9Voltage
            self.gpuAoBPV=gpu.aoBP
            self.gpuTime=gpu.time
            self.gpuN1Amp=gpu.mAmpSec4
            self.gpuN1Freq=gpu.mFreqSec4
            self.gpuN1Phase=gpu.mPhaseSec4
            self.plotOfGpuAmp=gpu.plotOfAmplitudes
            self.plotOfGpuBPV=gpu.plotOfBPS9Voltage
            self.plotOfGpuFreq=gpu.plotOfFreq
            self.plotOfGpuPhase=gpu.plotOfPhase

            self.gpuIndices=self.getTimeIndices()
            
        if "IP" in tags or "Ip" in tags:
            IPdata=IP(shotno,tStart=tStartFixed,tStop=tStopFixed) 
            self.Ip=IPdata.ip
            self.plotOfIp=IPdata.plotOfIP
            
        if "TP" in tags:
            TP=TPData(shotno,tStart=tStartFixed,tStop=tStopFixed)
            self.TPS2vf=TP.Vf2
            self.TPS8vf=TP.Vf8
            self.TPS2n=TP.ne2
            self.TPS8n=TP.ne8
            self.TPS2kTe=TP.kTe2
            self.TPS8kTe=TP.kTe8
        
        if save==True:
            self.save()
            
        
    def getTimeIndices(self):
        """
        provides indices of CPCI time that matches GPU time.  effectively downsamples CPCI data to GPU data
        """              
        gpuIndices=_np.zeros(len(self.gpuTime),dtype=int)
        for i in range(0,len(self.gpuTime)):
            gpuIndices[i]=_process.find_nearest(self.time,self.gpuTime[i])
            
        return gpuIndices
        
    def save(self):
#        import pickle
        filePath=_FILEDIR+str(int(self.shotno))+'.pickle'
#            pk.saveAsPickle(fileName)
        with open(filePath, 'w') as f:  # Python 3: open(..., 'wb')
            _pk.dump(self, f)
            
        
###############################################################################
### Concatenated shotno data format
"""
A different type of shot data structure.  This structure concatenates hbtep 
data so that it is easier to anayze data across multiple shot numbers. 
"""
class appendShotData:
    """
    Combines list of shotData() into a class structure of single arrays per data type.
    This makes processing lots of shot data easier as everything is an appended array.
    """
    
    def __init__(self,data):
        self.time=_np.zeros(0)
        self.n1Amp=_np.zeros(0)
        self.n1Phase=_np.zeros(0)
        self.n1PhaseVelRaw=_np.zeros(0)
        self.n1PhaseVelStrongFilter=_np.zeros(0)
        self.minorRadius=_np.zeros(0)
        self.majorRadius=_np.zeros(0)
        self.q=_np.zeros(0)
        self.qUpsample=_np.zeros(0)
        self.BPV=_np.zeros(0)
        self.BPI=_np.zeros(0)
        self.shotNos=_np.zeros(0,dtype=int)
        
        for i in range(0,len(data)):
            dataIn=data[i];
            
            try:
#                print "len(data.time) = %d" % len(dataIn.time)
#                print "len(data.q) = %d" % len(dataIn.q)
                self.time=_np.append(self.time,dataIn.time);
                self.n1Amp=_np.append(self.n1Amp,dataIn.n1Amp);
                self.n1Phase=_np.append(self.n1Phase,dataIn.n1Phase);
                self.n1PhaseVelRaw=_np.append(self.n1PhaseVelRaw,dataIn.n1PhaseVelRaw);
                self.n1PhaseVelStrongFilter=_np.append(self.n1PhaseVelStrongFilter,dataIn.n1PhaseVelStrongFilter);  
                self.shotNos=_np.append(self.shotNos,_np.array([int(dataIn.shotno)]*len(dataIn.time)))
            except AttributeError:
                print "shotno %d has no nMode data. Skipping" % (dataIn.shotno)
                
            try:
                self.minorRadius=_np.append(self.minorRadius,dataIn.minorRadius);
                self.majorRadius=_np.append(self.majorRadius,dataIn.majorRadius);
            except AttributeError:
                print "shotno %d has no radius data. Skipping" % (dataIn.shotno)
                                
            try:
#                self.q=_np.append(self.q,dataIn.q);
                # note that q had a different time base.  qUpsample is safe.
                self.qUpsample=_np.append(self.qUpsample,dataIn.qUpsample);
            except AttributeError:
                print "shotno %d has no q data. Skipping" % (dataIn.shotno)  
                            
            try:
                self.BPV=_np.append(self.BPV,dataIn.BPV);
                self.BPI=_np.append(self.BPI,dataIn.BPI);
            except AttributeError:
                print "shotno %d has no BP data. Skipping" % (dataIn.shotno)
                
            try:
                self.TPS2n=_np.append(self.TPS2n,dataIn.TPS2n);
                self.TPS2kTe=_np.append(self.TPS2kTe,dataIn.TPS2kTe);
                self.TPS2vf=_np.append(self.TPS2vf,dataIn.TPS2vf);
                self.TPS8n=_np.append(self.TPS8n,dataIn.TPS8n);
                self.TPS8kTe=_np.append(self.TPS8kTe,dataIn.TPS8kTe);
                self.TPS8vf=_np.append(self.TPS8vf,dataIn.TPS8vf);
            except AttributeError:
                print "shotno %d has no TP data. Skipping" % (dataIn.shotno)
                
    def filterData(self,indices): #e.g. indices=_np.where(self.n1Amp>2.)
        """
        Using a function call like indices=_np.where(self.n1Amp>2.), an array
        of indices can be generated.  This function removes all data for self 
        that does not match these indices
        """
        
        try:
            self.time=self.time[indices];
            self.n1Amp=self.n1Amp[indices];
            self.n1Phase=self.n1Phase[indices];
            self.n1PhaseVelRaw=self.n1PhaseVelRaw[indices];
            self.n1PhaseVelStrongFilter=self.n1PhaseVelStrongFilter[indices];
        except AttributeError:
            print "data has no nMode data. Skipping" 
            
        try:
            self.minorRadius=self.minorRadius[indices];
            self.majorRadius=self.majorRadius[indices];
        except AttributeError:
            print "data has no radius data. Skipping" 
                            
        try:
#            self.q=self.q[indices];  # note that q has a different time base
            self.qUpsample=self.qUpsample[indices];
        except AttributeError:
            print "data has no q data. Skipping" 
                        
        try:
            self.BPV=self.BPV[indices];
            self.BPI=self.BPI[indices];
        except AttributeError:
            print "data has no BP data. Skipping"
            
        try:
            self.TPS2n=self.TPS2n[indices];
            self.TPS2kTe=self.TPS2kTe[indices];
            self.TPS2vf=self.TPS2vf[indices];
            self.TPS8n=self.TPS8n[indices];
            self.TPS8kTe=self.TPS8kTe[indices];
            self.TPS8vf=self.TPS8vf[indices];
        except AttributeError:
            print "data has no TP data. Skipping"
            
            
            
###############################################################################
### Processed data from HBTEP
class nModeData:
    """
    This function performs a least squares fit to a toroidal array of sensors and analyzes n=1 and n=2 modes.  Mode amplitude, phase, and phase velocity. 
    In addtion, this code generates a perturbed B_pol(t) measurement as observed by a sensor at location, phi0
    Function uses either 30 poloidal TA sensors or 10 poloidal FB sensors. 
    
    Note:  I want to create a "mModeAnalysis" equivalent for measuring "m" mode numebrs.
    """    

    def Bn1(self,phi0=0):
        """
        Generates a pretend B_{n=1} signal at the toroidal location, phi0
        Not presently in use
        """
        return self.x[1,:]*_np.sin(self.phi0)+self.x[2,:]*_np.cos(self.phi0)
        
    def __init__(self,shotno=96530,tStart=None,tStop=None,plot=False,phi0=0,nModeSensor='FB',method='leastSquares'):
        
        self.shotno=shotno
        
        self.title = 'shotno = %d.  %s sensor.  n=1,2 mode analysis' % (shotno,nModeSensor)

        if nModeSensor=='TA':
            ## load TA data
            temp=taData(self.shotno,tStart,tStop);
            data=temp.taPolData
            self.time=temp.tbPolTime
            phi=temp.phi
            [n,m]=_np.shape(data)
        elif nModeSensor=='FB':
            ## load FB data
            temp=fbData(self.shotno,tStart=tStart,tStop=tStop);
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
            # self.n0Offset=_np.zeros(m)
            self.n1Amp=_np.zeros(m)
            self.n1Phase=_np.zeros(m)
            self.n2Amp=_np.zeros(m)
            self.n2Phase=_np.zeros(m)
            # TODO(John): remove for loop and convert into all matrix math 
            # should simplify code and make it run faster
            for j in range(0,m):
                y=_np.zeros(n);
                for i in range(0,n):
                    y[i]=data[i][j]*1e4
                x[:,j]=Ainv.dot(y)
                # self.n0Offset=self.x[0,j]
                self.n1Amp[j]=_np.sqrt(x[1,j]**2+x[2,j]**2)
                self.n2Amp[j]=_np.sqrt(x[3,j]**2+x[4,j]**2)
                self.n1Phase[j]=_np.arctan2(x[1,j],x[2,j])
                self.n2Phase[j]=_np.arctan2(x[3,j],x[4,j])
            self._x=x
            self.n1Phase*=-1  # for some reason, the slope of phase had the wrong sign.  this corrects that.
            self.n2Phase*=-1  # for some reason, the slope of phase had the wrong sign.  this corrects that.
            
            ## Calculate frequency (in Hz) using second order deriv 
            self.n1Freq=_np.gradient(_process.unwrapPhase(self.n1Phase))/_np.gradient(self.time)/(2*_np.pi)
            
            # boxcar filter of frequency
            self.n1FreqWeakFilter=_process.boxCar(self.n1Freq,10)
            self.n1FreqTimeWeakFilter = self.time        
            self.n1FreqStrongFilter=_process.boxCar(self.n1Freq,40)
            self.n1FreqTimeStrongFilter = self.time
            
            # boxcar filter of n1 amplitude
            self.n1AmpFiltered=_process.boxCar(self.n1Amp,30)
        else:
            _sys.exit("Invalid mode analysis method provided.")
        
        ## mode amplitude plots  
        self.plotOfAmps=_plot.plot()
        self.plotOfAmps.yData=[self.n1Amp,self.n2Amp,self.n1AmpFiltered]
        self.plotOfAmps.xData=[self.time*1000,self.time*1000,self.time*1000]
        if nModeSensor=='TA':
            self.plotOfAmps.yLegendLabel=['TA Sensors, n=1','TA Sensors, n=2','TA Sensors, n=2, filtered']
        elif nModeSensor=='FB':
            self.plotOfAmps.yLegendLabel=['FB Sensors, n=1','FB Sensors, n=2','FB Sensors, n=2, filtered']
        self.plotOfAmps.title=str(self.title)
        self.plotOfAmps.xLabel='ms'
        self.plotOfAmps.yLabel='T'    

        # n=1 mode amplitude
        self.plotOfN1Amp=_plot.plot()
        self.plotOfN1Amp.subtitle='n=1 mode amp'
        self.plotOfN1Amp.yLim=[0,10]
        self.plotOfN1Amp.yData=[self.n1Amp]
        self.plotOfN1Amp.xData=[self.time*1000]
        if nModeSensor=='TA':
            self.plotOfN1Amp.yLegendLabel=['TA Sensors']
        elif nModeSensor=='FB':
            self.plotOfN1Amp.yLegendLabel=['FB Lower Sensors']
        self.plotOfN1Amp.title=str(self.title)
        self.plotOfN1Amp.xLabel='ms'
        self.plotOfN1Amp.yLabel='G'                                      
                     
        # n=1 mode phase
        self.plotOfN1Phase=_plot.plot()
        self.plotOfN1Phase.subtitle='n=1 mode phase'
        self.plotOfN1Phase.yLim=[-_np.pi,_np.pi]
        self.plotOfN1Phase.yData=[self.n1Phase]
        self.plotOfN1Phase.xData=[self.time*1000]
        self.plotOfN1Phase.linestyle=['']
        self.plotOfN1Phase.marker=['.']
        if nModeSensor=='TA':
            self.plotOfN1Phase.yLegendLabel=['TA Sensors']
        elif nModeSensor=='FB':
            self.plotOfN1Phase.yLegendLabel=['FB Lower Sensors']
        self.plotOfN1Phase.title=str(self.title)
        self.plotOfN1Phase.xLabel='ms'
        self.plotOfN1Phase.yLabel='phi'
        
        # n=1 mode freq
        self.plotOfN1Freq=_plot.plot()        
        self.plotOfN1Freq.subtitle='n=1 mode frequency'
        self.plotOfN1Freq.yData=[self.n1Freq/1000.,self.n1FreqStrongFilter/1000.,self.n1FreqWeakFilter/1000.]
        self.plotOfN1Freq.xData=[self.time*1000,self.time*1000,self.time*1000]
        self.plotOfN1Freq.yLegendLabel=['raw','strong filter','weak filter']
        self.plotOfN1Freq.title=str(self.title)
        self.plotOfN1Freq.xLabel='ms'
        self.plotOfN1Freq.yLabel='kHz'
        self.plotOfN1Freq.yLim=[-20,20]
                   
        # hybrid plot of phase AND amplitude
        # TODO(John) implement in new plot function
        self.plotOfPhaseAmp=_plot.plot() 
        self.plotOfPhaseAmp.yData=[self.n1Phase]
        self.plotOfPhaseAmp.xData=[self.time*1000]
        self.plotOfPhaseAmp.colorData=[self.n1AmpFiltered]#[self.n1Amp]
        self.plotOfPhaseAmp.linestyle=['']
        self.plotOfPhaseAmp.marker=['.']
        self.plotOfPhaseAmp.subtitle='n=1 Phase and Filtered Amplitude'
        self.plotOfPhaseAmp.title=str(shotno)
        self.plotOfPhaseAmp.xLabel='ms'
        self.plotOfPhaseAmp.yLabel=r'$\phi$'
        self.plotOfPhaseAmp.zLabel='Gauss'
        self.plotOfPhaseAmp.yLegendLabel=['TA sensors']
        self.plotOfPhaseAmp.plotType='scatter'
        self.plotOfPhaseAmp.yLim=[-_np.pi,_np.pi]
        
#       TODO fix this code        
#        mx=_np.max(self.n1AmpFiltered)
#        lCutoff=2.5
#        uCutoff=8.
#        cm = _processt.singleColorMapWithLowerAndUpperCutoffs(lowerCutoff=lCutoff/mx,upperCutoff=uCutoff/mx)
#        self.plotOfPhaseAmp.cmap=cm
                                
        ## plot data
        if plot==True:
            self.plotSlice(index=int(m/4));
            self.plotSlice(index=int(m/2));
            _plot.subPlot([self.plotOfAmps,self.plotOfN1Phase,self.plotOfN1Freq])
#            self.plotOfPhaseAmp.plot();  #TODO fix this.  still having issues

    def plotSlice(self,index=0):
        """
        Plots fit data for a single time value
        """
        j=index;
        [n,m]=_np.shape(self._data)
        y=_np.zeros(n);
        for i in range(0,n):
                y[i]=self._data[i][j]*1e4
        self.p1=_plot.plot()
        phi=_np.linspace(self._phi[0],self._phi[-1],100)
        n1Fit=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)
        n2Fit=self._x[0,j]+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)
        fitTotal=self._x[0,j]+self._x[1,j]*_np.sin(phi)+self._x[2,j]*_np.cos(phi)+self._x[3,j]*_np.sin(2*phi)+self._x[4,j]*_np.cos(2*phi)

        # plot
        self.p1.yData=[y,n1Fit,n2Fit,fitTotal]
        self.p1.xData=[self._phi,phi,phi,phi]
        self.p1.linestyle=['','-','-','-']
        self.p1.marker=['.','','','']
        self.p1.yLegendLabel=['raw','n=1','n=2','n=1 + n=2']
        self.p1.color=['black','red','blue','green']
        self.p1.title='t='+str(self.time[j]*1000)+'ms.  '+str(self.shotno)
        self.p1.plot()
        

class mModeData:
    """
    This function performs a least squares fit to a poloidal array of sensors and analyzes m=2,3 and 4 modes.  Mode amplitude, phase, and phase velocity. 
    In addtion, this code generates a perturbed B_pol(t) measurement as observed by a sensor at location, theta0
    Function uses either 32 poloidal PA1 or PA2 sensors
    
    
    """  
    def __init__(self,shotno=96530,tStart=None,tStop=None,plot=False,theta0=0,sensor='PA1'):
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
        
        # TODO:  add frequency data, raw and smoothed
        
        # plot amplitudes
        self.plotOfAmplitudes=_plot.plot()
        self.plotOfAmplitudes.yData=[self.m1Amp,self.m2Amp,self.m3Amp,self.m4Amp,self.m5Amp]
        self.plotOfAmplitudes.xData=[self.time*1000,self.time*1000,self.time*1000,self.time*1000,self.time*1000]
        self.plotOfAmplitudes.yLegendLabel=[r'$|B_{pol, m=1}|$',r'$|B_{pol, m=2}|$',r'$|B_{pol, m=3}|$',r'$|B_{pol, m=4}|$',r'$|B_{pol, m=5}|$']
        self.plotOfAmplitudes.title=str(self.title)
        self.plotOfAmplitudes.xLabel='ms'
        self.plotOfAmplitudes.yLabel='G'
        
        if plot == True:
            self.plotSlice(index=int(m/4));
            self.plotSlice(index=int(m/2));
            self.plotOfAmplitudes.plot()
        
        
    def plotSlice(self,index=0):
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
        
        p1.yData=[y,m1Fit,m2Fit,m3Fit,m4Fit,m5Fit,fitTotal]
        p1.xData=[self._theta,theta,theta,theta,theta,theta,theta]
        p1.linestyle=['','-','-','-','-','-','-']
        p1.marker=['.','','','','','','']
        p1.yLegendLabel=['raw','m=1','m=2','m=3','m=4','m=5','m=1-5']
        p1.title='t=%.3f ms. %s ' % (self.time[j]*1000, self.title)
        p1.plot()
        
