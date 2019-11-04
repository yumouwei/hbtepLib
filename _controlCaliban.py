"""
Functions related to the caliban gpu feedback computer
"""

import numpy as _np
import hbtepLib as _hbt
#_plot=_hbt.plot
from socket import gethostname
import os
import sys

# Important for headless runn
#if 'DISPLAY' in os.environ and sys.version_info <= (3,0): #python3 check
if 'DISPLAY' in os.environ: # should be python3 compatable
    import matplotlib.pyplot as _plt
    _plot = _hbt.plot


if gethostname()=='spitzer':
    _LOCAL_DATA_PATH='/opt/hbt/data/control' 
elif os.environ['USER'] == 'rian':
    _LOCAL_DATA_PATH='/home/rian/Documents/control_data'

class calibanData:
    """
    reads control data from Caliban's control files
    under ACTIVE developement
    
    Parameters
    ----------
    shotno : int or None
        shot number of desired data
        if None, the code loads the simulated data from the /control directory
    tStart : float
        time (in seconds) to trim data before
        default is 0 ms
    tStop : float
        time (in seconds) to trim data after
        default is 10 ms
    plot : bool
        plots all relevant plots if true
        default is False
    password : str
        ssh password to gpu computer
    forceDownload : bool
        causes data to all be downloaded locally even if already there
    plotEVERYTHING : bool
        ...
    _LOCAL_DATA_DIR : str
        string to your remote directory to download the gpu data to
    cpciShotno : int
        If a cpci shotnumber if supplied, it will plot the CPCI mode data
        alongside the gpu data
    
    Attributes
    ----------
    
    Notes
    -----
    Presently, the function is setup to run remotely but could be modified to 
    operate locally without too much effort
    """
    def __init__(self,shotno=99591, tStart=0*1e-3, tStop=10*1e-3, plot=False, #shotno=96496
              password='', forceDownload=True, 
              plotEVERYTHING=False, remoteDataDir='',cpciShotno=None):

        import os     
        
        self.shotno=shotno;
        self.cpciShotno=cpciShotno
        
        # handles cases where shotno and cpciShotno are different or not defined
        if self.shotno!=None and self.cpciShotno==None:
            self.cpciShotno=self.shotno
            self.nModeData=_hbt.get.nModeData(shotno)
            #self.bpData=_hbt.get.bpData(shotno)
            self._cpci=True
        elif self.cpciShotno!=None:
            self.nModeData=_hbt.get.nModeData(cpciShotno)
            #self.bpData=_hbt.get.bpData(cpciShotno)
            self._cpci=True
        else:
            self._cpci=False
            
        # enforces a minimum time of 1ms
        if tStart<1.0e-3:
            tStart=1.0e-3
        
        self._tStart=tStart
        self._tStop=tStop
        self._plotEVERYTHING=plotEVERYTHING 
        
        # if operating remotely on another computer (not the HBT server)
        if remoteDataDir!='': #_hbt.get._ON_HBTEP_SERVER==False or 
            # ssh data must be transfered if operating remotely.
            # this copies control data to your local computer.
        
            # check to see if the file has previously been downloaded to local directory.  if not, it is downloaded
            if shotno!=None:
                filePath=remoteDataDir + "/fbsettings_" +str(int(shotno))+'.py';
            else:
                filePath=remoteDataDir + '/fbsettings.py';
            print(filePath)
            if os.path.isfile(filePath)==False or forceDownload==True:
#                print('Downloading %d data' % shotno)
                _downloadCDFromCaliban(shotno,password=password)
            dataDir=remoteDataDir
            
        # if operating locally on the HBT server
        else:
            dataDir=_LOCAL_DATA_PATH
            
        # load data and time from control files
        print("data directory = %s"%dataDir)
        self.totalNumSamples=get_totalNumSamples(shotno,dataPath=dataDir)
        analogIn=getAI(shotno,self.totalNumSamples,dataPath=dataDir)
        time=getTime(self.totalNumSamples)
        self.time=time
        analogOut=getAO(shotno,self.totalNumSamples,dataPath=dataDir);
        mAmp=getModeAmp(shotno,self.totalNumSamples,dataPath=dataDir);
        #mPhase=getModePhase(shotno,self.totalNumSamples,dataPath=dataDir);
        mPhase=getModePhase(shotno,dataPath=dataDir);
        mFreq=getModeFreq(shotno,self.totalNumSamples,dataPath=dataDir);
        
        #necessary step if running with EUV data
        analogIn = _np.transpose(analogIn)
        analogOut = _np.transpose(analogOut)
        mAmp = _np.transpose(mAmp)
        mPhase = _np.transpose(mPhase)
        mFreq = _np.transpose(mFreq)
        
        # trim time and data to specified time range
        temp,self.analogOut=_hbt.get._trimTime(time,list(_np.transpose(analogOut)),tStart,tStop)
        temp,self.analogIn=_hbt.get._trimTime(time,list(_np.transpose(analogIn)),tStart,tStop)
        temp,mAmp=_hbt.get._trimTime(time,list(_np.transpose(mAmp)),tStart,tStop)
        self.modeAmp=mAmp
        temp,mPhase=_hbt.get._trimTime(time,list(_np.transpose(mPhase)),tStart,tStop)
        self.modePhase=mPhase
        time,mFreq=_hbt.get._trimTime(time,list(_np.transpose(mFreq)),tStart,tStop)
        self.modeFreq=mFreq
        self.time=time
    
        # distribute trimmed data to class variables
        mAmpCosSec1=mAmp[0];
        mAmpSinSec1=mAmp[1];
        self.n1ModeAmpSec1=_np.sqrt(mAmpCosSec1**2 + mAmpSinSec1**2)
        mAmpCosSec2=mAmp[2];
        mAmpSinSec2=mAmp[3];
        self.n1ModeAmpSec2=_np.sqrt(mAmpCosSec2**2 + mAmpSinSec2**2)
        mAmpCosSec3=mAmp[4];
        mAmpSinSec3=mAmp[5];
        self.n1ModeAmpSec3=_np.sqrt(mAmpCosSec3**2 + mAmpSinSec3**2)
        mAmpCosSec4=mAmp[6];
        mAmpSinSec4=mAmp[7];
        self.n1ModeAmpSec4=_np.sqrt(mAmpCosSec4**2 + mAmpSinSec4**2)
        self.n1ModePhaseSec1=mPhase[0];
        self.n1ModePhaseSec2=mPhase[2];
        self.n1ModePhaseSec3=mPhase[4];
        self.n1ModePhaseSec4=mPhase[6];
        self.n1ModeFreqSec1=mFreq[0];
        self.n1ModeFreqSec2=mFreq[2];
        self.n1ModeFreqSec3=mFreq[4];
        self.n1ModeFreqSec4=mFreq[6];
        
        if plot == True:
            self.plot()
            
    def plot(self,plot=True):
        # Not necessarilly necessary to do voltage probe
        #return _plot.subPlot([self.plotOfAmplitudes(),self.plotOfPhase(),self.plotOfFreq()],plot=plot)
        return _plot.subPlot([self.plotOfAmplitudes(),self.plotOfPhase(),self.plotOfFreq()],plot=plot)
        
    def plotOfAmplitudes(self):
        # init mode amp plot
        
        p1=_hbt.plot.plot(title=str(self.shotno),
                    xLabel='Time (ms)',
                    yLabel='G',
                    yLim=[0,20],
                    subtitle='Mode amplitude, n=1')
        p1.addTrace(self.time*1000,self.n1ModeAmpSec4*1e4,yLegendLabel='GPU-FB_S4P')
        if self._cpci==True:
            p1.addTrace(self.nModeData.time*1e3,self.nModeData.n1Amp,yLegendLabel='CPCI-FB_S4P')
        if self._plotEVERYTHING==True:
            
            p1.addTrace(self.time*1000,self.n1ModeAmpSec1*1e4,yLegendLabel='GPU-Sec1')
            p1.addTrace(self.time*1000,self.n1ModeAmpSec2*1e4,yLegendLabel='GPU-Sec2')
            p1.addTrace(self.time*1000,self.n1ModeAmpSec3*1e4,yLegendLabel='GPU-Sec3')
            # p1.addTrace(self.time*1000,self.n1ModeAmpSec4*1e4,yLegendLabel='GPU-Sec4')
        return p1
         
    def plotOfPhase(self):
        # init mode phase plot
        p1=_hbt.plot.plot(title=str(self.shotno),
                    xLabel='Time (ms)',
                    yLabel='Radias',
                    yLim=[-_np.pi, _np.pi]  ,
                    subtitle='Mode phase, n=1')
        p1.addTrace(self.time*1000,_hbt.process.wrapPhase(self.n1ModePhaseSec4),
              yLegendLabel='GPU-FB_S4P',marker='.',linestyle='')
        if self._cpci==True:
            p1.addTrace(self.nModeData.time*1e3,self.nModeData.n1Phase,yLegendLabel='CPCI-FB_S4P',
                  marker='.',linestyle='',alpha=.4,markerSize=.1)#,markerSize=.1
        if self._plotEVERYTHING==True:
            p1.addTrace(self.time*1000,_hbt.process.wrapPhase(self.n1ModePhaseSec1),
                  yLegendLabel='GPU-Sec1',marker='.',linestyle='')
            p1.addTrace(self.time*1000,_hbt.process.wrapPhase(self.n1ModePhaseSec2),
                  yLegendLabel='GPU-Sec2',marker='.',linestyle='')
            p1.addTrace(self.time*1000,_hbt.process.wrapPhase(self.n1ModePhaseSec3),
                  yLegendLabel='GPU-Sec3',marker='.',linestyle='')
        return p1

    def plotOfFreq(self):
        # init mode frequency plot
        p1=_hbt.plot.plot(title=str(self.shotno),
                    xLabel='Time (ms)',
                    yLabel='kHz',
                    #yLim=[-_np.pi, _np.pi]  ,
                    subtitle='Mode frequency, n=1')
#        p1.addTrace(self.time*1000,self.n1ModeFreqSec4*1e-3, yLegendLabel='GPU-FB_S4P')
        if self._cpci==True:
            p1.addTrace(self.nModeData.time*1e3,self.nModeData.n1Freq*1e-3,yLegendLabel='CPCI-FB_S4P')
        if self._plotEVERYTHING==True:
            p1.addTrace(self.time*1000,self.n1ModeFreqSec1*1e-3, yLegendLabel='GPU-Sec1')
            p1.addTrace(self.time*1000,self.n1ModeFreqSec2*1e-3, yLegendLabel='GPU-Sec2')
            p1.addTrace(self.time*1000,self.n1ModeFreqSec3*1e-3, yLegendLabel='GPU-Sec3')
        return p1
  

def calcLeastSquaresMatrix(shotno):
    """
    Calculates the psuedo-inverse matrix used in the least squares 
    operation: x = A^-1 * b
    
    Parameters
    ----------
    shotno : int
        Looks at a recent shot number to 1) exclude the broken FB sensors
        from A and 2) get the theta and phi coordinates for those that ARE
        functional.  
        
    Returns
    -------
    invMtx : numpy.array (8 by number of working sensors)
        A^-1 matrix from the least squares operation: x = A^-1 * b
    outText : str
        String version of invMtx that can be copy and pasted into fbsettigs.h
    
    Example
    -------
    ::
        
        aInv,aInvString=hbt.get.fbData(101169)
    
    
    Code for sanity checking
    ------------------------    
    ### this commented code uses the old fbtools.py code to do roughly the 
    ### same thing.  I leave it here as a sanity check for those who follow
    import fbtools as fbt; 
        
    # 2D list of functioning sensor names
    sensors=[['FB01_S1P', 'FB02_S1P', 'FB03_S1P', 'FB04_S1P', 'FB05_S1P', 'FB06_S1P', 'FB07_S1P', 'FB08_S1P', 'FB09_S1P', 'FB10_S1P'], ['FB01_S2P', 'FB02_S2P', 'FB03_S2P', 'FB04_S2P', 'FB05_S2P', 'FB07_S2P', 'FB08_S2P', 'FB09_S2P', 'FB10_S2P'], ['FB01_S3P', 'FB02_S3P', 'FB03_S3P', 'FB05_S3P', 'FB06_S3P', 'FB07_S3P', 'FB08_S3P', 'FB09_S3P', 'FB10_S3P'], ['FB01_S4P', 'FB02_S4P', 'FB03_S4P', 'FB04_S4P', 'FB05_S4P', 'FB06_S4P', 'FB07_S4P', 'FB09_S4P', 'FB10_S4P']]
    #b2,b2a=fbt.get_IN_MODE_MATRIX(True)
    b,temp=fbt.get_IN_MODE_MATRIX_johnEdit(sensors) # not working?
    #a=hbt.get.fbData(101169)
    #b3,temp=fbt.get_IN_MODE_MATRIX_johnEdit(a.fbPolNames)
    
    
    """
    
    # get FB sensor data
    a=_hbt.get.fbData(101169)
    phi=a.phi # list of arrays of sensor phi coordinates
    theta=a.theta # list of arrays of sensor theta coordinates
    n=len(a.fbPolNames[0])+len(a.fbPolNames[1])+len(a.fbPolNames[2])+len(a.fbPolNames[3]) # total number of sensors
    
    # calculate least squares matrix
    outMtx=_np.zeros((n,8))
    n1=0
    for i in range(0,4):
        #print i
        mtx=_np.zeros((len(a.fbPolNames[i]),2))
        mtx[:,0]=1.0*_np.cos(3*_np.array(theta[i])-_np.array(phi[i]))
        mtx[:,1]=1.0*_np.sin(3*_np.array(theta[i])-_np.array(phi[i]))
        outMtx[n1:n1+len(a.fbPolNames[i]),i*2:i*2+2]=mtx
        n1+=len(a.fbPolNames[i])
        
    # pseudo-invert least squares matrix
    invMtx=_np.linalg.pinv(outMtx)
    
    # create text string to be placed in fbsettings.h
    outText=''
    outText+='#define IN_MODE_MATRIX {'
    for i in range(0,8):
        print(i)
        if i!=0:
            outText+="\t\t\t\t\t\t"
        outText+="{"
        for j in range(0,len(invMtx[i])):
            print(j)
            outText+="%.5f"%invMtx[i,j]
            if j!=len(invMtx[i])-1:
                outText+=",\t"
        if i!=7:
            outText+="},\\\n"
        else:
            outText+="}}"
        
    return invMtx,outText


def getAI(shotno,numColumns=37,dataPath=_LOCAL_DATA_PATH,plot=False):#,numCols=37):
    """ 
    Read data tfrom ai_store_<shotno>.dat. 
    
    Parameters
    ----------
    shotno : int
        Shot number
    numColumns : int
        Number of columns of the data in ai_store
    dataPath : str
        File path where ai_store is located
    plot : bool
        Plots all data
        
    Returns
    -------
    data : 2D np.ndarray
        2D array with dimension (numColumns x time).
        Contains all of the analog input data as recorded in the GPU.
    
    """
    
    if shotno==None:  #
        data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/ai_store.dat'%(dataPath),numColumns=numColumns,dataType=_np.float32)
    else:
        data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/ai_store_%d.dat'%(dataPath,shotno),numColumns=numColumns,dataType=_np.float32)#* 10. / INT16_MAX 

    if plot==True:
        
        time=getTime(len(data))
        fig,ax=_plt.subplots()
        ax.plot(time,data)
        ax.set_xlabel('time (ms)')
        _plt.show()

    return data


def getAO(shotno,numColumns=40,dataPath=_LOCAL_DATA_PATH,plot=False):#numCols=64):
    """ 
    Read data tfrom ao_store_<shotno>.dat. 
    
    Parameters
    ----------
    shotno : int
        Shot number
    numColumns : int
        Number of columns of the data in ao_store
    dataPath : str
        File path where ao_store is located
    plot : bool
        Plots all data
        
    Returns
    -------
    data : 2D np.ndarray
        2D array with dimension (numColumns x time).
        Contains all of the analog output data as recorded in the GPU.
    
    Notes
    -----
    
    #TODO(John) I'm not sure if this function works correctly...
    """
    INT16_MAX = _np.iinfo(_np.int16).max
    if shotno==None:
        data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/ao_store.dat'%(dataPath),numColumns=numColumns,dataType=_np.int16)* 10. / INT16_MAX 
    else:
        data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/ao_store_%d.dat'%(dataPath,shotno),numColumns=numColumns,dataType=_np.int16)* 10. / INT16_MAX 

    if plot==True:
        time=getTime(len(data))
        fig,ax=_plt.subplots()
        ax.plot(time,data)
        ax.set_xlabel('time (ms)')
        _plt.show()
                
    return data

def getFeedback(shotno,numCols=14,dataPath=_LOCAL_DATA_PATH,plot=False):
    """ 
    Read data tfrom fb_store_<shotno>.dat. 
    
    Notes
    -----
    This file must have been previously downloaded to your computer using the
    _downloadCDFromCaliban() function.
    """
    # for some reason, the data matrix isn't divisible by 14 (don't know why).  the code below coorects for that.
    if shotno==None:
        a=_hbt.readWrite.readBinaryFileInto2DMatrix('%s/fb_store.dat'%(dataPath),numColumns=1,dataType=_np.float32)
    else:
        a=_hbt.readWrite.readBinaryFileInto2DMatrix('%s/fb_store_%d.dat'%(dataPath,shotno),numColumns=1,dataType=_np.float32)
    b=_np.remainder(len(a),numCols)
    c= a[:-b].reshape((-1,numCols))
    
    # there is a lot of dead data on either end (don't konw why).  clipping
    if False:
        d=c[:,0]
        iStart=_np.where(d>0)[0][0]
        iStop=int(_np.where(d==d.max())[0][0])
        fb= c[iStart:iStop+1,:]
    else:
        iStop=int(_np.where(c[:,0]==c[:,0].max())[0][0])
        fb=c[:iStop+1,:]
    
    if plot == True:
        fig,ax=_plt.subplots(numCols-1,sharex=True)
        # assuming the first column is time
        for j in range(0,numCols-1):
            ax[j].plot(fb[:,0],fb[:,j+1])
        _plt.show()
    
    return fb
    

def getModeAmp(shotno,numColumns=8,dataPath=_LOCAL_DATA_PATH,plot=False):#numCols=8):
    """ 
    Read mode amplitude data from mamp_store_<shotno>.dat. 
    
    Parameters
    ----------
    shotno : int
        Shot number
    numColumns : int
        Number of columns of the data in mamp_store
    dataPath : str
        File path where mamp_store is located
    plot : bool
        Plots all data
        
    Returns
    -------
	data : 2D np.ndarray
		2D array with dimension (numColumns x time).
		Contains all of the mode amplitude data as recorded in the GPU.
	
	"""
    if shotno==None:
        data = _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mamp_store.dat'%(dataPath),numColumns=numColumns,dataType=_np.float32)
    else:
        data = _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mamp_store_%d.dat'%(dataPath,shotno),numColumns=numColumns,dataType=_np.float32)
    
    if plot==True:
        time=getTime(len(data))
        fig,ax=_plt.subplots(2,sharex=True)
        ax[0].plot(time*1e3,data)
        ax[1].set_xlabel('time (ms)')
        for i in range(numColumns/2):
            ax[1].plot(time*1e3,_np.sqrt(data[:,i*2]**2+data[:,i*2+1]**2))
        _plt.show()

    return data


def getModeFreq(shotno,numColumns=8,dataPath=_LOCAL_DATA_PATH,plot=False):
    """ 
    Read mode freq (Hz) data from mfreq_store_<shotno>.dat. 
    
    Parameters
    ----------
    shotno : int
        Shot number
    numColumns : int
        Number of columns of the data in mfreq_store
    dataPath : str
        File path where mfreq_store is located
    plot : bool
        Plots all data
        
    Returns
    -------
    data : 2D np.ndarray
        2D array with dimension (numColumns x time).
        Contains all of the mode freq data as recorded in the GPU.
    
    """
    if shotno==None:
        data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mfreq_store.dat'%(dataPath),numColumns=numColumns,dataType=_np.float32)
    else:
        data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mfreq_store_%d.dat'%(dataPath,shotno),numColumns=numColumns,dataType=_np.float32)
        
    if plot==True:
        time=getTime(len(data))
        fig,ax=_plt.subplots()    
        ax.plot(time,data)
        ax.set_xlabel('time (ms)')
        _plt.show()
        
    return data


def getModePhase(shotno,numColumns=8,dataPath=_LOCAL_DATA_PATH,plot=False):
    """ 
    Read mode phase (rad) data from mphase_store_<shotno>.dat. 
    
    Parameters
    ----------
    shotno : int
        Shot number
    numColumns : int
        Number of columns of the data in mphase_store
    dataPath : str
        File path where mphase_store is located
    plot : bool
        Plots all data
        
    Returns
    -------
    data : 2D np.ndarray
        2D array with dimension (numColumns x time).
        Contains all of the mode phase data as recorded in the GPU.
    """
    if shotno==None:
        data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mphase_store.dat'%(dataPath),numColumns=numColumns,dataType=_np.float32)
    else:
        data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mphase_store_%d.dat'%(dataPath,shotno),numColumns=numColumns,dataType=_np.float32)
        
    if plot==True:
        
        time=getTime(len(data))
        fig,ax=_plt.subplots()    
        ax.plot(time,data%(2*_np.pi))
        ax.set_xlabel('time (ms)')
        _plt.show()
        
    return data


def getTime(numSamples=1231,offsetSamples=-165,cycleTime=6e-6): # -165 sample offset appears correct for comparing input signals
    return _np.arange(offsetSamples,numSamples+offsetSamples)*cycleTime


def get_totalNumSamples(shotno,dataPath=_LOCAL_DATA_PATH):
    """ 
    Get the total number of samples associated with all GPU data (with the 
    exception of fb_data).
    
    Notes
    -----
    Files must have been previously downloaded to your computer using the
    _downloadCDFromCaliban() function.
    """
    print(shotno)
    if shotno==None:(totalNumSamples,n)=_np.shape(_hbt.readWrite.readBinaryFileInto2DMatrix('%s/mamp_store.dat'%(dataPath),numColumns=8,dataType=_np.float32))
    else:
        (totalNumSamples,n)=_np.shape(_hbt.readWrite.readBinaryFileInto2DMatrix('%s/mamp_store_%d.dat'%(dataPath,shotno),numColumns=8,dataType=_np.float32))
    return totalNumSamples

        
def _downloadCDFromCaliban(shotno,
                   password='',
                   localCodePath='/home/brooks/TokaMac/control',
                   remoteFileDir='/home/john/shotData',
                   hbtServerUsername="brooks",
                   localDataPath='/opt/hbt/data/control'  ):
    """
    Downloads control data files from control computer using ssh to remote computer
    """
    
    # get password
    if password=='' or password == None:
        #password = raw_input("Enter spitzer password:  ")
#            password=_hbt._rwDataTools.getPwd(systemName=_hbt._hbtPreferences._HBT_SERVER_NAME,userName=hbtServerUsername); #username=_pref.hbtServerUsername
        try:
            password=_hbt._rwDataTools.getPwd(systemName=_hbt._hbtPreferences._HBT_SERVER_NAME,userName=hbtServerUsername); #username=_pref.hbtServerUsername
        except:
            print("Your password is not stored within the keyring.  See hbtepLib.readWrite.getPwd() for details.")

    # open connection
#    print('Downloading %d' % shotno)
    sshCon = _hbt._rwDataTools.scpData(password=password,port=22,username=hbtServerUsername,address=_hbt._hbtPreferences._HBT_SERVER_ADDRESS) #username=_pref.hbtServerUsername

    # _copy data
    if shotno=='' or shotno==None:
        # copy the simulated feedback files
        sshCon.downloadFile('%s/ao_store.dat' % (localCodePath), localFilePath=remoteFileDir+'/ao_store.dat' )
        sshCon.downloadFile('%s/ai_store.dat' % (localCodePath), localFilePath=remoteFileDir+'/ai_store.dat' )
        try:
            sshCon.downloadFile('%s/airaw_store.dat' % (localCodePath))
        except:
            print("raw file not present.  skipping...")
            pass
        sshCon.downloadFile('%s/mamp_store.dat' % (localCodePath), localFilePath=remoteFileDir+'/mamp_store.dat' )
        sshCon.downloadFile('%s/mphase_store.dat' % (localCodePath), localFilePath=remoteFileDir+'/mphase_store.dat' )
        sshCon.downloadFile('%s/mfreq_store.dat' % (localCodePath), localFilePath=remoteFileDir+'/mfreq_store.dat' )
        sshCon.downloadFile('%s/fb_store.dat' % (localCodePath), localFilePath=remoteFileDir+'/fb_store.dat' )
        sshCon.downloadFile('%s/fbsettings.py' % (localCodePath), localFilePath=remoteFileDir+'/fbsettings.py' )
    else:
        # copy the feedback files associated with actual shot numbers
        print(localDataPath)
        sshCon.downloadFile('%s/ao_store_%d.dat' % (localDataPath, shotno), localFilePath='%s/ao_store_%d.dat' % (remoteFileDir, shotno))
        sshCon.downloadFile('%s/ai_store_%d.dat' % (localDataPath, shotno), localFilePath='%s/ai_store_%d.dat' % (remoteFileDir, shotno))
        sshCon.downloadFile('%s/mamp_store_%d.dat' % (localDataPath, shotno), localFilePath='%s/mamp_store_%d.dat' % (remoteFileDir, shotno))
        sshCon.downloadFile('%s/mphase_store_%d.dat' % (localDataPath, shotno), localFilePath='%s/mphase_store_%d.dat' % (remoteFileDir, shotno))
        sshCon.downloadFile('%s/mfreq_store_%d.dat' % (localDataPath, shotno), localFilePath='%s/mfreq_store_%d.dat' % (remoteFileDir, shotno))
        try:
            sshCon.downloadFile('%s/fbsettings_%d.py' % (localDataPath, shotno), localFilePath='%s/fbsettings_%d.py'  % (remoteFileDir, shotno))
        except:
            print("fb_settings.py file not present.  skipping...")
            pass
        
        try:
            # print '%s/fb_store_%d.dat' % (localDataPath, shotno)
            sshCon.downloadFile('%s/fb_store_%d.dat' % (localDataPath, shotno), localFilePath='%s/fb_store_%d.dat'  % (remoteFileDir, shotno))
        except Exception:
            print("fb_store file not present.  skipping...")     
            pass

    # close connection
    sshCon.closeConnection();


def prepAwg(waveform):
    """
    Creates feedforward signal to be piped to the GPU's analog output by the
    ./do_awg function in the /control directory.  The output file is 
    'awgdata.dat'
    
    waveform : 2D numpy.ndarray
        2D numpy array with dimensnions (NUMBER_OF_SAMPLES,AO_CHANNELS)
    
    Example
    -------
    ::
        
        import numpy as np
        import matplotlib.pyplot as plt
        
        CYCLE_TIME=6e-6
        AO_CHANNELS=60
        
        # create time
        SAMPLES=int(10e-3//CYCLE_TIME)+2 # 10 ms worth of samples
        t=np.arange(0,10e-3+CYCLE_TIME,CYCLE_TIME)-1e-3; 
        
        # create signals
        y1=np.zeros(len(t))
        y1[t>1.5e-3]=1.
        y1[t>2e-3]+=np.sin(2*np.pi*(t[t>2e-3]-2e-3)*4000)
        index1=41;
            
        y2=np.zeros(len(t))
        y2[t>1.25e-3]=-1.
        index2=43;
        
        # add signals to 2D waveform matrix
        waveform = np.zeros((SAMPLES,AO_CHANNELS))
        waveform[:,index1]=y1
        waveform[:,index2]=y2
        
        # plot all waveform signals as a sanity check
        plt.plot(t,waveform)
        
        prepAwg(waveform)
    """
    import struct
    INT16_MAX = _np.iinfo(_np.int16).max
    AO_CHANNELS=60 # number of channels (i.e. width of waveform matrix)
    
    fh = open('awgdata.dat', 'wb') 
    fh.write(struct.pack('=h', AO_CHANNELS))
    fh.write((waveform * INT16_MAX / 10).astype(_np.int16).tostring())
