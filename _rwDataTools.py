"""
_rwDataTools.py - read write generic data functions and related tools (not HBT 
specific)

"""


import numpy as _np
import _hbtPreferences as _pref
#reload(_pref)
import pickle as _pk


_HBT_SERVER_USERNAME='brooks'

###############################################################################
### Read/write to pickle
# TODO implement these functions

      
###############################################################################
### Text file data storage
      
def writeDataToFile(data,filename='out.csv',delimiter=',',
                    writeNumRowsAtHeader=False):
    """
    writes data to text file
    
    Parameters
    ----------
    data : numpy.ndarray (1D or 2D)
        data to be written to file
    delimeter : str
        delimeter between entries
    writeNumRowsAtHeader : bool

    Notes
    -----
    input data must be a 1 or 2D numpy array only.  
    
    References
    ----------
    https://stackoverflow.com/questions/37289951/python-write-to-csv-line-by-line
    """
#    import csv
    n=data.ndim
    out = open(filename, 'w');
    if writeNumRowsAtHeader==True:
        out.write(str(len(data))+'\n')
    if n==1:
        for i in range(len(data)):
            out.write(str(data[i])+delimiter+'\n')
    if n==2:
        for i in range(len(data)):
            for j in range(len(data[0])):
                out.write(str(data[i,j])+delimiter)
                if j==len(data[0])-1:
                    out.write('\n')
    out.close()
    
# TODO need readDataFromFile code
    
###############################################################################
### binary files
def readBinaryFileInto2DMatrix(fileName,numColumns=None,numRows=None,dataType=_np.float32):
	"""
	Reads a binary data file into a 1D or 2D numpy array, depending on the 
	number of columns specified.
	
	# fileName='/home/john/Downloads/ai_store.dat'
	
	Example #1
	----------
	::
		
		# write and then read 1D binary file
		
		fileName='out.dat'
		a=np.arange(0,21,dtype=np.float32)
		write2DNumpyArrayToSeqBinaryFile(a,fileName)
		data=readSeqBinaryFileInto2DMatrix(fileName)
			
	Example #2
	----------
	::
		
		# write and then read 2D binary file
		
		fileName='out.dat'
		numColumns=3
		a=np.arange(0,21,dtype=np.float32).reshape(-1,numColumns)
		write2DNumpyArrayToSeqBinaryFile(a,fileName)
		data=readSeqBinaryFileInto2DMatrix(fileName,numColumns=numColumns)
	"""
	if type(numColumns)==type(None) and type(numRows)==type(None):
		numColumns=1
	
	if type(numRows)==type(None):
		data= _np.fromfile(fileName, dataType).reshape((-1,numColumns))
	elif type(numColumns)==type(None):
		data= _np.fromfile(fileName, dataType).reshape((numRows,-1))
	(m,n)=_np.shape(data)
	if n==1:
		return data[:,0]
	else:
		return data	
    
###############################################################################
### keyring password management for python
       
def setPwd(password,system=_pref._HBT_SERVER_NAME,username=_HBT_SERVER_USERNAME):
    """ 
    Encrypts password using keyring, a password management tool.  
    
    Parameters
    ----------
    password : str
        ssh password for hbtep server
    system : str
        name of ssh hbtep server.  
    username : str
        name of login on ssh hbtep server.  
    
    NOTES
    -----
    this function also ONLY works on a work station where the OS-based 
    function keyring is installed and the password has already been set for 
    that user.  i'm using this on ubuntu.  not sure if it'll work on windows
    """
    import keyring
    keyring.set_password(system,username,password)
    
    
def getPwd(systemName=_pref._HBT_SERVER_NAME,userName=_HBT_SERVER_USERNAME):
    """ 
    Returns unencrypted password
    
    Parameters
    ----------
    password : str
        ssh password for hbtep server
    system : str
        name of ssh hbtep server.  
    username : str
        name of login on ssh hbtep server.  
        
    Returns
    -------
    : str
        ssh password 
    
    NOTES
    -----
    this function also ONLY works on a work station where the OS-based 
    function keyring is installed and the password has already been set for 
    that user.  i'm using this on ubuntu.  not sure if it'll work on windows
    """
    import keyring as kr
    return str(kr.get_password(systemName, userName))

    
###############################################################################
### ssh and scp code

class scpData:
    """
    ssh and scp combination toolkit for file transfers
    
    Parameters
    ----------
    password : str
        password for ssh server
    address : str
        address for ssh server
    port : int
        port number for ssh server
    username : str
        username for ssh server
 
    Attributes
    ----------
    ssh : 'paramiko.SSHClient'
    scp : 'scp.SCPClient'
    
    Notes
    -----    
    
    References
    ----------
    https://gist.github.com/stonefury/06ab3531a1c30c3b998a
    https://github.com/jbardin/scp.py
    
    Examples
    --------
    
    """        
#    def __init__(self,password=getPwd(),username=_pref._HBT_SERVER_USERNAME,
#                 address=_pref._HBT_SERVER_ADDRESS,port=22):
                     
    def __init__(self,password,username=_HBT_SERVER_USERNAME,
                 address=_pref._HBT_SERVER_ADDRESS,port=22):
        from paramiko import SSHClient
        from scp import SCPClient
        self.ssh = SSHClient()
        self.ssh.load_system_host_keys()
        self.ssh.connect(address, port=port, username=username, 
                         password=password)
        self.scp = SCPClient(self.ssh.get_transport())
        
    def downloadFile(self,remoteFilePath,localFilePath=''):
#        try:
        self.scp.get(remoteFilePath, localFilePath)
#        except:
#            print("%s not present.  skipping..." % remoteFilePath)
#            pass
        
    def uploadFile(self,localFilePath,remoteFilePath):
        self.scp.put(localFilePath,remoteFilePath)
        
    def uploadFolder(self,localDirPath,remoteDirPath):
        """ upload folder and all contents """
        self.scp.put(localDirPath,remoteDirPath,recursive=True)
        
    def closeConnection(self):
        self.scp.close();
        self.ssh.close();
        
    def __del__(self):
        """ upon deletion of object, close ssh/scp connections """
        self.closeConnection();


    
###############################################################################
### save/read pickle files
    
def saveToPickle(data,fileName):
    """
    Save data (any object) to file
    
    Parameters
    ----------
    data : any object or class instance
        data to be saved
    fileName : str
        the COMPLETE file and directory name.  Also please include an extension
        e.g. fileName='/home/john/shotData/98030.pickle'
        
    Notes
    -----
    Pickles are very large files and shouldn't be used excessively

    """
    fileHandler = open(fileName,'w')
    _pk.dump(data,fileHandler)
    fileHandler.close()
    
    
def loadFromPickle(fileName):
    """
    Load and return data (any object) from a pickle file
    
    Parameters
    ----------
    fileName : str
        the COMPLETE file and directory name.  Also please include an extension
        e.g. fileName='/home/john/shotData/98030.pickle'
        
    Returns
    -------
    data : any object or class instance
        loaded data
        
    Notes
    -----
    Pickles are very large files and shouldn't be used excessively

    """
    fileHandler = open(fileName,'r')
    data=_pk.load(fileHandler)
    return data
    
  
###############################################################################
### dictionaries
  
class dictionary:
    """
    Wrapper for python's dictionary tool
    
    Notes
    -----
    It's probably easier to just use the dictionary tool without this wrapper.
    However, I plan on keeping this wrapper around because it's nice to remind
    myself of the dictionaries functionality
    
    References
    ----------
    https://www.tutorialspoint.com/python/python_dictionary.htm
    """
    def __init__(self,listOfKeys=[],listOfEntries=[]):
        m=len(listOfKeys)
        
        # initialize empty dictionary
        self.dict={}
        
        # add elements
        if listOfEntries==[]:
#            for i in range(0,m):
            self.initializeKeysWithNumpyArrays(listOfKeys)
        else:
            for i in range(0,m):
                self.addEntry(listOfKeys[i],listOfEntries[i])
            
    def initializeKeysWithNumpyArrays(self,listOfKeys):
        for i in range(0,len(listOfKeys)):
            self.addEntry(listOfKeys[i],_np.zeros(0));
            
    def indexEntry(self,key):
        """ pulls entry from dictionary """
        return self.dict[key]
        
    def addEntry(self, key, entry):
        """ add entry """
        self.dict[key]=entry
        
    def updateEntry(self, key, entry):
        """ same as add entry but the key should already exist """
        self.addEntry(key,entry) 
        
    def delEntry(self,key):
        """ delete entry """
        del self.dict[key]
        
    def returnListOfKeys(self):
        """ returns a list of the keys """
        return self.dict.keys()
        
    def appendNumpyArrayToEntry(self,key,entry):
        self.dict[key]=_np.append(self.dict[key],entry)