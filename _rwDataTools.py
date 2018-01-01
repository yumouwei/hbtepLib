"""
_rwDataTools.py - read write generic data functions and related tools (not HBT 
specific)

"""


#import numpy as _npd
import _hbtPreferences as _pref
reload(_pref)
import pickle as _pk


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
    
# TODO do the same as above for binary as well (for better file storage)
    
    
###############################################################################
### keyring password management for python
       
def setPwd(password,system=_pref._HBT_SERVER_NAME,username=_pref._HBT_SERVER_USERNAME):
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
    
    
def getPwd(system=_pref._HBT_SERVER_NAME,username=_pref._HBT_SERVER_USERNAME):
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
    import keyring
    return str(keyring.get_password(system, username))

    
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
    def __init__(self,password=getPwd(),username=_pref._HBT_SERVER_USERNAME,
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