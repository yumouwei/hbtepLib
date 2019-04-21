"""
Functions related to the caliban gpu feedback computer
"""

import numpy as _np
import matplotlib.pyplot as _plt
import hbtepLib as _hbt
_plot=_hbt.plot

class gpuControlData:
	"""
	reads control data from Caliban's control files
	under ACTIVE developement
	
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
	password : str
		ssh password to gpu computer
	forceDownload : bool
		causes data to all be downloaded locally even if already there
	calibrateTime : bool
		...
	plotEVERYTHING : bool
		...
	operatingMode : str
		'' or 'currentControl' or 'frequencyControl'
	
	Attributes
	----------
	
	Notes
	-----
	"""
	def __init__(self,shotno=99591, tStart=0*1e-3, tStop=10*1e-3, plot=False, #shotno=96496
			  password='', forceDownload=True, 
			  plotEVERYTHING=False, _LOCAL_DATA_DIR='/home/john/shotData',cpciShotno=None):

		import os	 
		
		self.shotno=shotno;
		self.cpciShotno=cpciShotno
		if self.shotno!=None and self.cpciShotno==None:
			self.cpciShotno=self.shotno
			self.nModeData=_hbt.get.nModeData(shotno)
			self.bpData=_hbt.get.bpData(shotno)
			self.cpci=True
		elif self.cpciShotno!=None:
			self.nModeData=_hbt.get.nModeData(cpciShotno)
			self.bpData=_hbt.get.bpData(cpciShotno)
			self.cpci=True
		else:
			self.cpci=False
		if tStart<1.0e-3:
			tStart=1.0e-3
		self._tStart=tStart
		self._tStop=tStop
		self._plotEVERYTHING=plotEVERYTHING 
		
		if _hbt.get._ON_HBTEP_SERVER==False:
			# ssh data must be transfered if operating remotely.
			# this copies control data to your local computer.
		
			# check to see if the file has previously been downloaded to local directory.  if not, it is downloaded
			if shotno!=None:
				filePath=_LOCAL_DATA_DIR + "/fbsettings_" +str(int(shotno))+'.py';
			else:
				filePath=_LOCAL_DATA_DIR + '/fbsettings.py';
			print filePath
			if os.path.isfile(filePath)==False or forceDownload==True:
#				print('Downloading %d data' % shotno)
				_downloadCDFromCaliban(shotno,password=password)
			
		# load data and time from control files
		time=get_ctrl_times(shotno)
		self.time=time
		self.TOTAL_SAMPLES=get_total_samples(shotno)
		analogIn=get_ctrl_ai(shotno,self.TOTAL_SAMPLES)
		analogOut=get_ctrl_ao(shotno,self.TOTAL_SAMPLES);
		mAmp=get_ctrl_mamp(shotno,self.TOTAL_SAMPLES);
		mPhase=get_ctrl_mphase(shotno,self.TOTAL_SAMPLES);
		mFreq=get_ctrl_mfreq(shotno,self.TOTAL_SAMPLES);
	
		# trim time and data
		temp,self.analogOut=_hbt.get._trimTime(time,list(_np.transpose(analogOut)),tStart,tStop)
		temp,analogIn=_hbt.get._trimTime(time,list(_np.transpose(analogIn)),tStart,tStop)
		temp,mAmp=_hbt.get._trimTime(time,list(_np.transpose(mAmp)),tStart,tStop)
		temp,mPhase=_hbt.get._trimTime(time,list(_np.transpose(mPhase)),tStart,tStop)
		time,mFreq=_hbt.get._trimTime(time,list(_np.transpose(mFreq)),tStart,tStop)
		self.time=time

		# distribute trimmed data to class variables
		self.BP1VoltageReq=self.analogOut[41]; #41 or 43 (SOUTH_CPCI_10 channels 41 or 43)
		self.BP2VoltageReq=self.analogOut[43]; #41 or 43 (SOUTH_CPCI_10 channels 41 or 43)
		mAmpCosSec1=mAmp[0];
		mAmpSinSec1=mAmp[1];
		self.mAmpSec1=_np.sqrt(mAmpCosSec1**2 + mAmpSinSec1**2)
		mAmpCosSec2=mAmp[2];
		mAmpSinSec2=mAmp[3];
		self.mAmpSec2=_np.sqrt(mAmpCosSec2**2 + mAmpSinSec2**2)
		mAmpCosSec3=mAmp[4];
		mAmpSinSec3=mAmp[5];
		self.mAmpSec3=_np.sqrt(mAmpCosSec3**2 + mAmpSinSec3**2)
		mAmpCosSec4=mAmp[6];
		mAmpSinSec4=mAmp[7];
		self.mAmpSec4=_np.sqrt(mAmpCosSec4**2 + mAmpSinSec4**2)
		self.mPhaseSec1=mPhase[0];
		self.mPhaseSec2=mPhase[2];
		self.mPhaseSec3=mPhase[4];
		self.mPhaseSec4=mPhase[6];
		self.mFreqSec1=mFreq[0];
		self.mFreqSec2=mFreq[2];
		self.mFreqSec3=mFreq[4];
		self.mFreqSec4=mFreq[6];
		
		if plot == True:
			self.plot()
		
	def plot(self,plot=True):
		return _plot.subPlot([self.plotOfAmplitudes(),self.plotOfPhase(),self.plotOfFreq(),self.plotOfProbeReqVoltage()],plot=plot)

	def plotOfProbeReqVoltage(self):
		# initialize BPS9 voltage plot
		p1=_hbt.plot.plot(title=str(self.shotno),
					xLabel='Time [ms]',
					yLabel='V',
#					yLim=[-11,11],
					subtitle='Probe Voltage Request')
		p1.addTrace(self.time*1000,self.BP1VoltageReq,yLegendLabel='GPU-BP1_Req.')
		p1.addTrace(self.time*1000,self.BP2VoltageReq,yLegendLabel='GOU-BP2_Req.')
#		if self.cpci==True:
#			p1.addTrace(self.bpData.time*1e3,self.bpData.bps9Voltage,yLegendLabel='CPCI')
		
#		p1.yLegendLabel.append('GPU-BPS9')
		return p1
		
	def plotOfAmplitudes(self):
		# init mode amp plot
		
		p1=_hbt.plot.plot(title=str(self.shotno),
					xLabel='Time (ms)',
					yLabel='G',
					yLim=[0,20],
					subtitle='Mode amplitude, n=1')
		p1.addTrace(self.time*1000,self.mAmpSec4*1e4,yLegendLabel='GPU-FB_S4P')
		if self.cpci==True:
			p1.addTrace(self.nModeData.time*1e3,self.nModeData.n1Amp,yLegendLabel='CPCI-FB_S4P')
		if self._plotEVERYTHING==True:
			
			p1.addTrace(self.time*1000,self.mAmpSec1*1e4,yLegendLabel='GPU-Sec1')
			p1.addTrace(self.time*1000,self.mAmpSec2*1e4,yLegendLabel='GPU-Sec2')
			p1.addTrace(self.time*1000,self.mAmpSec3*1e4,yLegendLabel='GPU-Sec3')
			# p1.addTrace(self.time*1000,self.mAmpSec4*1e4,yLegendLabel='GPU-Sec4')
		return p1
		 
	def plotOfPhase(self):
		# init mode phase plot
		p1=_hbt.plot.plot(title=str(self.shotno),
					xLabel='Time (ms)',
					yLabel='Radias',
					yLim=[-_np.pi, _np.pi]  ,
					subtitle='Mode phase, n=1')
		p1.addTrace(self.time*1000,_hbt.process.wrapPhase(self.mPhaseSec4),
			  yLegendLabel='GPU-FB_S4P',marker='.',linestyle='')
		if self.cpci==True:
			p1.addTrace(self.nModeData.time*1e3,self.nModeData.n1Phase,yLegendLabel='CPCI-FB_S4P',
				  marker='.',linestyle='',alpha=.4,markerSize=.1)#,markerSize=.1
		if self._plotEVERYTHING==True:
			p1.addTrace(self.time*1000,_hbt.process.wrapPhase(self.mPhaseSec1),
				  yLegendLabel='GPU-Sec1',marker='.',linestyle='')
			p1.addTrace(self.time*1000,_hbt.process.wrapPhase(self.mPhaseSec2),
				  yLegendLabel='GPU-Sec2',marker='.',linestyle='')
			p1.addTrace(self.time*1000,_hbt.process.wrapPhase(self.mPhaseSec3),
				  yLegendLabel='GPU-Sec3',marker='.',linestyle='')
		return p1

	def plotOfFreq(self):
		# init mode frequency plot
		p1=_hbt.plot.plot(title=str(self.shotno),
					xLabel='Time (ms)',
					yLabel='kHz',
					#yLim=[-_np.pi, _np.pi]  ,
					subtitle='Mode frequency, n=1')
#		p1.addTrace(self.time*1000,self.mFreqSec4*1e-3, yLegendLabel='GPU-FB_S4P')
		if self.cpci==True:
			p1.addTrace(self.nModeData.time*1e3,self.nModeData.n1Freq*1e-3,yLegendLabel='CPCI-FB_S4P')
		if self._plotEVERYTHING==True:
			p1.addTrace(self.time*1000,self.mFreqSec1*1e-3, yLegendLabel='GPU-Sec1')
			p1.addTrace(self.time*1000,self.mFreqSec2*1e-3, yLegendLabel='GPU-Sec2')
			p1.addTrace(self.time*1000,self.mFreqSec3*1e-3, yLegendLabel='GPU-Sec3')
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
        print i
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


def get_ctrl_ai(shotno,TOTAL_SAMPLES=None,LOCAL_FILE_DIR='/home/john/shotData',plot=False):#,numCols=37):
	""" 
	Read data tfrom ai_store_<shotno>.dat. 
	
	Notes
	-----
	This file must have been previously downloaded to your computer using the
	_downloadCDFromCaliban() function.
	"""
	
	#INT16_MAX = _np.iinfo(_np.int16).max
	#TODO(John) I'm not sure if this function works correctly...
	if type(TOTAL_SAMPLES)==type(None):
		TOTAL_SAMPLES=get_total_samples(shotno)
	if shotno==None:  #
		data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/ai_store.dat'%(LOCAL_FILE_DIR),numColumns=37,dataType=_np.float32)
	else:
#		return _hbt.readWrite.readBinaryFileInto2DMatrix('%s/ai_store_%d.dat'%(LOCAL_FILE_DIR,shotno),numRows=TOTAL_SAMPLES,dataType=_np.float32)
		data= _hbt.readWrite.readBinaryFileInto2DMatrix('%s/ai_store_%d.dat'%(LOCAL_FILE_DIR,shotno),numColumns=37,dataType=_np.float32)#* 10. / INT16_MAX 

	if plot==True:
		fig,ax=_plt.subplots()
		ax.plot(data*1e4)
		ax.set_ylabel('Gauss')

	return data


def get_ctrl_ao(shotno,TOTAL_SAMPLES=None,LOCAL_FILE_DIR='/home/john/shotData'):#numCols=64):
	""" 
	Read data tfrom ao_store_<shotno>.dat. 
	
	Notes
	-----
	This file must have been previously downloaded to your computer using the
	_downloadCDFromCaliban() function.
	"""
	if type(TOTAL_SAMPLES)==type(None):
		TOTAL_SAMPLES=get_total_samples(shotno)
	
	INT16_MAX = _np.iinfo(_np.int16).max
	if shotno==None:
		return _hbt.readWrite.readBinaryFileInto2DMatrix('%s/ao_store.dat'%(LOCAL_FILE_DIR),numRows=TOTAL_SAMPLES,dataType=_np.int16)* 10. / INT16_MAX 
	else:
		return _hbt.readWrite.readBinaryFileInto2DMatrix('%s/ao_store_%d.dat'%(LOCAL_FILE_DIR,shotno),numRows=TOTAL_SAMPLES,dataType=_np.int16)* 10. / INT16_MAX 
		  

def get_fb(shotno,numCols=14,LOCAL_FILE_DIR='/home/john/shotData'):
	""" 
	Read data tfrom fb_store_<shotno>.dat. 
	
	Notes
	-----
	This file must have been previously downloaded to your computer using the
	_downloadCDFromCaliban() function.
	"""
	# for some reason, the data matrix isn't divisible by 14 (don't know why).  the code below coorects for that.
	if shotno==None:
		a=_hbt.readWrite.readBinaryFileInto2DMatrix('%s/fb_store.dat'%(LOCAL_FILE_DIR),numColumns=1,dataType=_np.float32)
	else:
		a=_hbt.readWrite.readBinaryFileInto2DMatrix('%s/fb_store_%d.dat'%(LOCAL_FILE_DIR,shotno),numColumns=1,dataType=_np.float32)
	b=_np.remainder(len(a),numCols)
	c= a[:-b].reshape((-1,numCols))
	
	# there is a lot of dead data on either end (don't konw why).  clipping
	d=c[:,0]
	iStart=_np.where(d>0)[0][0]
	iStop=int(_np.where(d==d.max())[0][0])
	return c[iStart:iStop+1,:]
	

def get_ctrl_mamp(shotno,TOTAL_SAMPLES=None,LOCAL_FILE_DIR='/home/john/shotData',plot=False):#numCols=8):
	""" 
	Read mode amplitude data from mamp_store_<shotno>.dat. 
	
	Notes
	-----
	This file must have been previously downloaded to your computer using the
	_downloadCDFromCaliban() function.
	"""
	if type(TOTAL_SAMPLES)==type(None):
		TOTAL_SAMPLES=get_total_samples(shotno)
	if shotno==None:
		data = _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mamp_store.dat'%(LOCAL_FILE_DIR),numRows=TOTAL_SAMPLES,dataType=_np.float32)
	else:
		data = _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mamp_store_%d.dat'%(LOCAL_FILE_DIR,shotno),numRows=TOTAL_SAMPLES,dataType=_np.float32)
	
	if plot==True:
		fig,ax=_plt.subplots(2,sharex=True)
		b=_np.zeros((len(data),4))
		b[:,0]=1e4*_np.sqrt(data[:,0]**2+data[:,1]**2)
		b[:,1]=1e4*_np.sqrt(data[:,2]**2+data[:,3]**2)
		b[:,2]=1e4*_np.sqrt(data[:,4]**2+data[:,5]**2)
		b[:,3]=1e4*_np.sqrt(data[:,6]**2+data[:,7]**2)
		for i in range(4):
			ax[0].plot(data[:,2*i],label='FB_S%iP cos'%(i+1))
			ax[0].plot(data[:,2*i+1],label='FB_S%iP sin'%(i+1))
			ax[1].plot(b[:,i],label='FB_S%iP'%(i+1))
#			plt.plot(b[:,i],label='%d'%i)
		ax[0].legend()
		ax[1].legend()

	return data


def get_ctrl_mfreq(shotno,TOTAL_SAMPLES=None,LOCAL_FILE_DIR='/home/john/shotData'):
	""" 
	Read mode freq (Hz) data from mamp_store_<shotno>.dat. 
	
	Notes
	-----
	This file must have been previously downloaded to your computer using the
	_downloadCDFromCaliban() function.
	"""
	if type(TOTAL_SAMPLES)==type(None):
		TOTAL_SAMPLES=get_total_samples(shotno)
		
	if shotno==None:
		return _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mfreq_store.dat'%(LOCAL_FILE_DIR),numRows=TOTAL_SAMPLES,dataType=_np.float32)
	else:
		return _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mfreq_store_%d.dat'%(LOCAL_FILE_DIR,shotno),numRows=TOTAL_SAMPLES,dataType=_np.float32)


def get_ctrl_mphase(shotno,TOTAL_SAMPLES=None,LOCAL_FILE_DIR='/home/john/shotData'):
	""" 
	Read mode phase (rad) data from mamp_store_<shotno>.dat. 
	
	Notes
	-----
	This file must have been previously downloaded to your computer using the
	_downloadCDFromCaliban() function.
	"""
	if type(TOTAL_SAMPLES)==type(None):
		TOTAL_SAMPLES=get_total_samples(shotno)
		
	if shotno==None:
		return _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mphase_store.dat'%(LOCAL_FILE_DIR),numRows=TOTAL_SAMPLES,dataType=_np.float32)
	else:
		return _hbt.readWrite.readBinaryFileInto2DMatrix('%s/mphase_store_%d.dat'%(LOCAL_FILE_DIR,shotno),numRows=TOTAL_SAMPLES,dataType=_np.float32)


def get_ctrl_times(shotno,time_offset=-166*6e-6,CYCLE_TIME=6e-6,LOCAL_FILE_DIR='/home/john/shotData'):
	""" 
	Get the GPU time data associated with each shot number
	
	Notes
	-----
	Files must have been previously downloaded to your computer using the
	_downloadCDFromCaliban() function.
	"""
	TOTAL_SAMPLES=get_total_samples(shotno)
	return _np.arange(0, TOTAL_SAMPLES) * CYCLE_TIME+time_offset


def get_total_samples(shotno,LOCAL_FILE_DIR='/home/john/shotData'):
	""" 
	Get the total number of samples associated with all GPU data (with the 
	exception of fb_data).
	
	Notes
	-----
	Files must have been previously downloaded to your computer using the
	_downloadCDFromCaliban() function.
	"""
	print(shotno)
	if shotno==None:(TOTAL_SAMPLES,n)=_np.shape(_hbt.readWrite.readBinaryFileInto2DMatrix('%s/mamp_store.dat'%(LOCAL_FILE_DIR),numColumns=8,dataType=_np.float32))
	else:
		(TOTAL_SAMPLES,n)=_np.shape(_hbt.readWrite.readBinaryFileInto2DMatrix('%s/mamp_store_%d.dat'%(LOCAL_FILE_DIR,shotno),numColumns=8,dataType=_np.float32))
	return TOTAL_SAMPLES

		
def _downloadCDFromCaliban(shotno,
				   password='',
				   REMOTE_CODE_PATH='/home/brooks/TokaMac/control',
				   LOCAL_FILE_DIR='/home/john/shotData',
				   _HBT_SERVER_USERNAME="brooks",
				   _REMOTE_DATA_PATH='/opt/hbt/data/control'  ):
	"""
	Downloads control data files from control computer using ssh
	"""
	
	# get password
	if password=='' or password == None:
		#password = raw_input("Enter spitzer password:  ")
		password=_hbt._rwDataTools.getPwd(systemName=_hbt._hbtPreferences._HBT_SERVER_NAME,userName=_HBT_SERVER_USERNAME); #username=_pref._HBT_SERVER_USERNAME
		try:
			password=_hbt._rwDataTools.getPwd(systemName=_hbt._hbtPreferences._HBT_SERVER_NAME,userName=_HBT_SERVER_USERNAME); #username=_pref._HBT_SERVER_USERNAME
		except:
			print("Your password is not stored within the keyring.  See hbtepLib.readWrite.getPwd() for details.")

	# open connection
#	print('Downloading %d' % shotno)
	sshCon = _hbt._rwDataTools.scpData(password=password,port=22,username=_HBT_SERVER_USERNAME,address=_hbt._hbtPreferences._HBT_SERVER_ADDRESS) #username=_pref._HBT_SERVER_USERNAME

	# _copy data
	if shotno=='' or shotno==None:
		# copy the simulated feedback files
		sshCon.downloadFile('%s/ao_store.dat' % (REMOTE_CODE_PATH), localFilePath=LOCAL_FILE_DIR+'/ao_store.dat' )
		sshCon.downloadFile('%s/ai_store.dat' % (REMOTE_CODE_PATH), localFilePath=LOCAL_FILE_DIR+'/ai_store.dat' )
		try:
			sshCon.downloadFile('%s/airaw_store.dat' % (REMOTE_CODE_PATH))
		except:
			print("raw file not present.  skipping...")
			pass
		sshCon.downloadFile('%s/mamp_store.dat' % (REMOTE_CODE_PATH), localFilePath=LOCAL_FILE_DIR+'/mamp_store.dat' )
		sshCon.downloadFile('%s/mphase_store.dat' % (REMOTE_CODE_PATH), localFilePath=LOCAL_FILE_DIR+'/mphase_store.dat' )
		sshCon.downloadFile('%s/mfreq_store.dat' % (REMOTE_CODE_PATH), localFilePath=LOCAL_FILE_DIR+'/mfreq_store.dat' )
		sshCon.downloadFile('%s/fb_store.dat' % (REMOTE_CODE_PATH), localFilePath=LOCAL_FILE_DIR+'/fb_store.dat' )
		sshCon.downloadFile('%s/fbsettings.py' % (REMOTE_CODE_PATH), localFilePath=LOCAL_FILE_DIR+'/fbsettings.py' )
	else:
		# copy the feedback files associated with actual shot numbers
		print(_REMOTE_DATA_PATH)
		sshCon.downloadFile('%s/ao_store_%d.dat' % (_REMOTE_DATA_PATH, shotno), localFilePath='%s/ao_store_%d.dat' % (LOCAL_FILE_DIR, shotno))
		sshCon.downloadFile('%s/ai_store_%d.dat' % (_REMOTE_DATA_PATH, shotno), localFilePath='%s/ai_store_%d.dat' % (LOCAL_FILE_DIR, shotno))
		sshCon.downloadFile('%s/mamp_store_%d.dat' % (_REMOTE_DATA_PATH, shotno), localFilePath='%s/mamp_store_%d.dat' % (LOCAL_FILE_DIR, shotno))
		sshCon.downloadFile('%s/mphase_store_%d.dat' % (_REMOTE_DATA_PATH, shotno), localFilePath='%s/mphase_store_%d.dat' % (LOCAL_FILE_DIR, shotno))
		sshCon.downloadFile('%s/mfreq_store_%d.dat' % (_REMOTE_DATA_PATH, shotno), localFilePath='%s/mfreq_store_%d.dat' % (LOCAL_FILE_DIR, shotno))
		sshCon.downloadFile('%s/fbsettings_%d.py' % (_REMOTE_DATA_PATH, shotno), localFilePath='%s/fbsettings_%d.py'  % (LOCAL_FILE_DIR, shotno))
		
		try:
			# print '%s/fb_store_%d.dat' % (_REMOTE_DATA_PATH, shotno)
			sshCon.downloadFile('%s/fb_store_%d.dat' % (_REMOTE_DATA_PATH, shotno), localFilePath='%s/fb_store_%d.dat'  % (LOCAL_FILE_DIR, shotno))
		except Exception:
			print("fb_store file not present.  skipping...")	 
			pass

	# close connection
	sshCon.closeConnection();
