
from __future__ import print_function, absolute_import, division
from hbtepold.find_shots import get_calm_period
from hbtepold.mdsplus import signals_from_shot
from hbtepold.misc import poly_fit_weights 
#import MDSplus
import logging
import math
import numpy as np
import scipy.integrate
from scipy.linalg import block_diag
import os
import types
import matplotlib.pyplot as plt

_fileDir='/home/john/shotData'


#DATA_PATH='/opt/hbt/data/control'
#REF_PATH='/home/john'
REF_PATH=_fileDir
DATA_PATH=_fileDir
#REF_PATH='/home/pq/Repository/fbsystem/control'
CFG_CACHE=dict()
INT16_MAX = np.iinfo(np.int16).max
AO_CHANNELS = 64
AI_CHANNELS = 96

log = logging.getLogger('fbtools')

class FBConfiguration(object):
    '''A container object to hold FB configurations'''
    pass

def get_calibration_matrix(cfg):
    '''Return sensor calibration factors
    
    The returned matrix converts the raw sensors readings from
    volt-seconds to Tesla, taking into account amplifier gain,
    number of turns and polarity.
    '''
    
    tree = MDSplus.Tree('hbtep2', -1)
    cal = np.identity(len(cfg.SENSORS))
    
    for (i, name) in enumerate(cfg.SENSORS):
        na = tree.getNode('.sensors.magnetic:%s:na' % name).data()
        gain = tree.getNode('.sensors.magnetic:%s:gain' % name).data()
        polarity = tree.getNode('.sensors.magnetic:%s:polarity' % name).data()
        
        cal[i,i] /= na * gain * polarity
                 
    return cal

def get_acq_channel_offset_map():
    '''Mapping from ACQ channels to memory offsets in LL mode
    
    Channels are indexed starting with 0.
    '''
    
    fh = open('/opt/hbt/data/acq196_channel_map')
    
    map_ = np.zeros((96), dtype=np.int16)
    for line in fh:
        (ch, offset) = [ int(x) for x in line.split() ]
        map_[ch-1] = offset
    
    return map_
    
def slope_fit(pts_cnt):
    '''Return weights to determine slope from *pts* points'''
    
    A = np.empty((pts_cnt, 2))
    for i in np.arange(2):
        A[:,i] = np.arange(pts_cnt)**i
        
    A_i = np.linalg.pinv(A)
    
    return A_i[1,::-1]
 

def get_ctrl_times(shotno):
    '''Return control system sampling times'''

    cfg = get_fb_config(shotno)
    return np.arange(-cfg.OFFSET_SAMPLES, 
                      cfg.TOTAL_SAMPLES-cfg.OFFSET_SAMPLES) * cfg.CYCLE_TIME


def get_ctrl_ai(shotno):
    '''Return control system input in shot *shotno*'''
    
    cfg = get_fb_config(shotno)
    if shotno is None:
        path = '%s/ai_store.dat' % DATA_PATH
    else:
        path = '%s/ai_store_%d.dat' % (DATA_PATH, shotno)
#     return np.fromfile(path, np.float32).reshape((-1, len(cfg.SENSORS))).T
    tmp = np.fromfile(path, np.float32).reshape((-1, len(cfg.SENSORS)+0 + 0)).T #+1 is added
    #print(len(cfg.SENSORS))
    ### the +1 is because sensors count was not sync with the sensors blacklist.  This is fixed after shot 87384
    return tmp[:-1,:]
    
"""(JOHN) I added this to look at all of AI data """
def get_ctrl_ai_numCh(shotno,n=96):
    '''Return control system input in shot *shotno*'''
    
#    cfg = get_fb_config(shotno)
    if shotno is None:
        path = '%s/ai_store.dat' % DATA_PATH
    else:
        path = '%s/ai_store_%d.dat' % (DATA_PATH, shotno)
#     return np.fromfile(path, np.float32).reshape((-1, len(cfg.SENSORS))).T
#    print path
    tmp2 = np.fromfile(path,np.float32);  
    print ("np.size(data)=%d" % tmp2.shape)
    tmp = tmp2.reshape((-1, n + 0)).T 
    #print(len(cfg.SENSORS))
    ### the +1 is because sensors count was not sync with the sensors blacklist.  This is fixed after shot 87384
    return tmp[:,:]


def get_ctrl_airaw(shotno):
    '''Return raw control system input in shot *shotno*''' # added by John 10/08/17
    
    cfg = get_fb_config(shotno)
    if shotno is None:
        path = 'airaw_store.dat' 
    else:
        path = '%s/airaw_store_%d.dat' % (DATA_PATH, shotno)
#     return np.fromfile(path, np.float32).reshape((-1, len(cfg.SENSORS))).T
    tmp = np.fromfile(path, np.float32).reshape((-1, len(cfg.SENSORS) + 0)).T 
    #print(len(cfg.SENSORS))
    ### the +1 is because sensors count was not sync with the sensors blacklist.  This is fixed after shot 87384
    return tmp[:-1,:]



def get_ctrl_mamp(shotno):
    '''Return mode amplitudes in shot *shotno*'''
    
    cfg = get_fb_config(shotno)
    if shotno is None:
        path = '%s/mamp_store.dat' % REF_PATH
    else:
        path = '%s/mamp_store_%d.dat' % (DATA_PATH, shotno)
    return np.fromfile(path, np.float32).reshape((-1, len(cfg.IN_MODE_MATRIX))).T
    
def get_fb(shotno, n=14):
    '''Return mode amplitudes in shot *shotno*'''
    
#    cfg = get_fb_config(shotno)
    if shotno is None or shotno is '':
        path = '%s/fb_store.dat' % REF_PATH
    else:
        path = '%s/fb_store_%d.dat' % (DATA_PATH, shotno)
        
    # this next bit of code corrects the fact that the data returned isn't always divisible by n, the num of channels.  not sure why.  
    a=np.fromfile(path, np.float32)
    b=np.remainder(len(a),n)
    if b!= 0:   
        b=n-b;
    a=np.append(a,np.zeros(b))
    
    # return data, formated into a matrix
    return a.reshape((-1,n)).T


def get_ctrl_mphase(shotno):
    '''Return control system mode phases in shot *shotno*'''
    
    cfg = get_fb_config(shotno)
    if shotno is None:
        path = '%s/mphase_store.dat' % REF_PATH
    else:
        path = '%s/mphase_store_%d.dat' % (DATA_PATH, shotno)
    return np.fromfile(path, np.float32).reshape((-1, len(cfg.IN_MODE_MATRIX))).T


def get_ctrl_mfreq(shotno):
    '''Return control system mode frequencies in shot *shotno*'''
    
    cfg = get_fb_config(shotno)
    if shotno is None:
        path = '%s/mfreq_store.dat' % REF_PATH
    else:
        path = '%s/mfreq_store_%d.dat' % (DATA_PATH, shotno)
    return np.fromfile(path, np.float32).reshape((-1, len(cfg.IN_MODE_MATRIX))).T

            
def get_ctrl_ao(shotno):
    '''Return control system output in shot *shotno*'''
    
    if shotno is None:
        path = '%s/ao_store.dat'% REF_PATH
    else:
        path = '%s/ao_store_%d.dat' % (DATA_PATH, shotno)
        
    raw = np.fromfile(path, np.int16).reshape((-1, AO_CHANNELS)).T
    return raw.astype(np.float) * 10 / INT16_MAX                        

class NoFBConfiguration(Exception):
    '''Raised if there is no feedback configuration for a shot'''
    
    pass


def get_fb_config(shotno):
    '''Return feedback configuration in shot *shotno*'''
    
    if shotno in CFG_CACHE:
        return CFG_CACHE[shotno]
    
    cfg = FBConfiguration()
    
    if shotno is None:
        execfile('%s/fbsettings.py'% REF_PATH , cfg.__dict__)
        return cfg
     
    fname =  '%s/fbsettings_%d.py' % (DATA_PATH, shotno)
    print(fname)
    if not os.path.exists(fname):
        raise NoFBConfiguration()
        
    execfile(fname, cfg.__dict__)
    CFG_CACHE[shotno] = cfg
    
    return cfg

def same_fb_config(cfg1, cfg2):
    '''Return true if configurations are equal'''
    
    if cfg1 is cfg2:
        return True

    if cfg1 is None or cfg2 is None:
        return False

    return bool(cmp_fb_settings(cfg1, cfg2))

def cmp_fb_settings(cfg1, cfg2):
    '''Return settings that have different values in *cfg1* and *cfg2*'''
    
    names1 = set(x for (x,y) in cfg1.__dict__.iteritems() 
                 if not x.startswith('__') and x.isupper()
                 and not isinstance(y, (types.ModuleType, types.FunctionType)))
    names2 = set(x for (x,y) in cfg2.__dict__.iteritems() 
                 if not x.startswith('__') and x.isupper()
                 and not isinstance(y, (types.ModuleType, types.FunctionType)))
    
    # Settings only present in one shot
    diff = names1.symmetric_difference(names2)
        
    # Settings with different values
    for n in names1.intersection(names2):
        at1 = getattr(cfg1, n)
        at2 = getattr(cfg2, n)
        if isinstance(at1, np.ndarray):
            if not np.allclose(at1, at2):
                diff.add(n)
        elif at1 != at2:
            diff.add(n)
        
    return diff

def calc_freq(times, phase, freq_dt):
    '''Calculate frequency from given phase array
    
    Frequency is calculated by least-squares fitting to a linear
    phase shift over the last *phase_dt* seconds.
    '''
    
    dt = times[1] - times[0]
    assert np.allclose(times[1:] - times[:-1], dt)
    freq_samples = int(round(freq_dt / dt)) 
    freq_weights = slope_fit(freq_samples)
    freq = np.empty_like(phase)
    
    if len(freq.shape) == 1:
        freq = np.convolve(phase, freq_weights, mode='full')[:freq.shape[1]]
    else:
        for j in xrange(freq.shape[0]):
            freq[j] = np.convolve(phase[j], freq_weights, mode='full')[:freq.shape[1]]
            
    return freq / (dt * 2 * math.pi)


def get_perturbed_signals(tree, fbcfg, t1, t2):
    '''Get perturbed signals for given feedback configuration'''
    
    # We need to retrieve before target time for equilibrium fitting to work properly.
    (times, signals) = signals_from_shot(tree, [ '.sensors.magnetic:%s' % x for x in fbcfg.SENSORS ],
                                         t_start=t1 - fbcfg.FIT_SAMPLES*fbcfg.CYCLE_TIME,
                                         t_end=t2)

    # Subtract equilibrium
    weights = poly_fit_weights(1, fbcfg.FIT_SAMPLES).astype(np.float32)
    for data in signals:
        data -= np.convolve(data, weights, mode='full')[:len(data)]
               
    # Restrict to target time
    idx = np.abs(times - t1).argsort()[0]    
    return (times[idx:], signals[:, idx:])

def get_SensorInfo(romer=False):
#    romer=True
    '''Read HBT-EP sensor information'''          
    if romer:
        suffix = 'ROMER'
    else:
        suffix = 'Nominal'            
    data_info = dict()
    filedir=os.path.dirname(__file__)  
    #########
    radImportance=49.0
    polImportance=1.0
    #  use relative dir because it might change in home directory while debugging
#        for filename in ('%s/sensor_info_PA_%s.txt' % (filedir,suffix), 
#                         '%s/sensor_info_FB_%s.txt' % (filedir,suffix)): 
    #filename='%s/../data/sensor_info_FB_%s.txt' % (filedir,suffix)
    filename='/home/john/Dropbox/Columbia/Research/sensor_info_FB_%s.txt' % (suffix)
    with open(filename,'r') as fh:                
        for line in fh:# Skip over comments
            if (line.split()[0] =='Columns'):
                break    
        # Store column names and positions
        line = fh.next()
        column_map = [ (pos, name[1:-1]) for (pos, name) 
                      in enumerate(line.split()) ]   
        fh.next()   # Skip one line     
        # Read actual data
        for line in fh:
            fields = line.split()
            data = dict( (column_name, float(fields[column_index]))
                         for (column_index, column_name) in column_map[1:] )
            sensor = fields[0]
            R = np.sqrt(data['loc_x']**2 + data['loc_y']**2)
            phi = np.arctan2(data['loc_y'], data['loc_x'])
            theta = np.arctan2(data['loc_z'], R-0.92)
            data['loc_r']= np.sqrt((R-0.92)**2+data['loc_z']**2)
            data['loc_phi']=phi
            data['loc_theta']=theta
            data['loc_MR']=R
            # the direction in the file has sign error
            data['n_x']=-data['n_x'];data['n_y']=-data['n_y'];data['n_z']=-data['n_z'];
            if sensor[-1]=='R':
                data['importance']=radImportance
            else:
                data['importance']=polImportance
            data_info[sensor]=data   
    return data_info  

def get_IN_MODE_MATRIX_select_sensors(sensor_list): #### get the in mode matrix with selected sensors
    sensor_info = get_SensorInfo()
    used_sensors_list = []
    p_in_mtx_inv_list = []
    p_in_mtx_list = []
    bad_sensor_list=get_Bad_Sensorlist()
    for i in xrange(4): # For each toroidal array
        tmp = [ 'FB%02d_S%dP' % (j,i+1) for j in np.arange(10)+1 ]
        sensors = []
        for name in tmp:
            if name in bad_sensor_list:
                continue
            if name in sensor_list:
                sensors.append(name)
                used_sensors_list.append(name)
        theta=[]
        phi=[]
        for name in sensors:
            theta.append(sensor_info[name]['loc_theta'])
            phi.append(sensor_info[name]['loc_phi'])
        theta=np.array(theta)
        phi=np.array(phi)
        ############ calculate the mode projection
        p_in_mtx=np.zeros((len(sensors),2))
        p_in_mtx[:,0]=1.0*np.cos(3*theta-phi)
        p_in_mtx[:,1]=1.0*np.sin(3*theta-phi)
  #      p_in_mtx[:,0]+=0.8*np.cos(2*theta-phi)
  #      p_in_mtx[:,1]+=0.8*np.sin(2*theta-phi) 
        p_in_mtx_list.append(p_in_mtx)       
        p_in_mtx_inv_list.append(np.linalg.pinv(p_in_mtx))
    p_in_mtx=block_diag(*p_in_mtx_list)
    p_in_mtx_i=block_diag(*p_in_mtx_inv_list)
    ################# radial sensors
#     r_in_mtx_inv_list=[]
#     r_in_mtx_list=[]
#     for i in xrange(4): # For each toroidal array
#         tmp = [ 'FB%02d_S%dR' % (j,i+1) for j in np.arange(10)+1 ]
#         sensors = []
#         for name in tmp:
#             if name in sensor_list:
#                 sensors.append(name)
#         theta=[]
#         phi=[]
#         for name in sensors:
#             theta.append(sensor_info[name]['loc_theta'])
#             phi.append(sensor_info[name]['loc_phi'])
#         theta=np.array(theta)
#         phi=np.array(phi)
#         ############ calculate the mode projection
#         r_in_mtx=np.zeros((len(sensors),2))
#         r_in_mtx[:,0]=1.0*np.cos(3*theta-phi)
#         r_in_mtx[:,1]=1.0*np.sin(3*theta-phi)
#    #     r_in_mtx[:,0]+=0.8*np.cos(2*theta-phi)
#    #     r_in_mtx[:,1]+=0.8*np.sin(2*theta-phi)
#         r_in_mtx_list.append(r_in_mtx)
#         r_in_mtx_inv_list.append(np.linalg.pinv(r_in_mtx))
#     r_in_mtx=block_diag(*r_in_mtx_list)
#     r_in_mtx_i=block_diag(*r_in_mtx_inv_list)
    ##########################
    #in_mtx_i=block_diag(p_in_mtx_i,r_in_mtx_i)
    in_mtx_i=p_in_mtx_i  ### only use the poloidal sensors
    return (in_mtx_i, used_sensors_list)

def get_IN_MODE_MATRIX(IfOriginal=False):
    sensor_info=get_SensorInfo()    
    ######### poloidal sensors
    sensors=[]
    Badlist=get_Bad_Sensorlist()
    p_in_mtx_inv_list=[]
    p_in_mtx_list=[]
    for i in xrange(4): # For each toroidal array
#         sensors.extend([ 'FB%02d_S%dP' % (j,i+1) for j in np.arange(10)+1 ])
        sensors=[ 'FB%02d_S%dP' % (j,i+1) for j in np.arange(10)+1 ]
        for name in Badlist:
            if name not in sensors:
                continue
            idx = sensors.index(name)
            del sensors[idx]
        theta=[]
        phi=[]
        for name in sensors:
            theta.append(sensor_info[name]['loc_theta'])
            phi.append(sensor_info[name]['loc_phi'])
        theta=np.array(theta)
        phi=np.array(phi)
        ############ calculate the mode projection
        p_in_mtx=np.zeros((len(sensors),2))
        p_in_mtx[:,0]=1.0*np.cos(3*theta-phi)
        p_in_mtx[:,1]=1.0*np.sin(3*theta-phi)
  #      p_in_mtx[:,0]+=0.8*np.cos(2*theta-phi)
  #      p_in_mtx[:,1]+=0.8*np.sin(2*theta-phi) 
        p_in_mtx_list.append(p_in_mtx)       
        p_in_mtx_inv_list.append(np.linalg.pinv(p_in_mtx))
    p_in_mtx=block_diag(*p_in_mtx_list)
    p_in_mtx_i=block_diag(*p_in_mtx_inv_list)
    ################# radial sensors
    sensors=[]
    r_in_mtx_inv_list=[]
    r_in_mtx_list=[]
    for i in xrange(4): # For each toroidal array
#         sensors.extend([ 'FB%02d_S%dR' % (j,i+1) for j in np.arange(10)+1 ])
        sensors=[ 'FB%02d_S%dR' % (j,i+1) for j in np.arange(10)+1 ]
        for name in Badlist:
            if name not in sensors:
                continue
            idx = sensors.index(name)
            del sensors[idx]
        theta=[]
        phi=[]
        for name in sensors:
            theta.append(sensor_info[name]['loc_theta'])
            phi.append(sensor_info[name]['loc_phi'])
        theta=np.array(theta)
        phi=np.array(phi)
        ############ calculate the mode projection
        r_in_mtx=np.zeros((len(sensors),2))
        r_in_mtx[:,0]=1.0*np.cos(3*theta-phi)
        r_in_mtx[:,1]=1.0*np.sin(3*theta-phi)
   #     r_in_mtx[:,0]+=0.8*np.cos(2*theta-phi)
   #     r_in_mtx[:,1]+=0.8*np.sin(2*theta-phi)
        r_in_mtx_list.append(r_in_mtx)
        r_in_mtx_inv_list.append(np.linalg.pinv(r_in_mtx))
    r_in_mtx=block_diag(*r_in_mtx_list)
    r_in_mtx_i=block_diag(*r_in_mtx_inv_list)
    ##########################
    #in_mtx_i=block_diag(p_in_mtx_i,r_in_mtx_i)
    in_mtx_i=p_in_mtx_i
    if IfOriginal: ### also return the original mapping  matrix
        #in_mtx=block_diag(p_in_mtx,r_in_mtx)
        in_mtx=p_in_mtx
        return (in_mtx_i,in_mtx)
    else:
        return in_mtx_i

def get_CCInfo():
    # get ref_phi angle
    sensor_info=get_SensorInfo()
    phi0=sensor_info['FB01_S1R']['loc_phi']+18*np.pi/180 ### offset to the right
    CC_info=dict()
    for sect in range(1,11):        
        for posi in range(1,5):
            CC=dict()
            name='FB%02d_C%d' % (sect,posi)
            phi=phi0+(sect-1)*2*np.pi/10
            dphi=31*np.pi/360
            if (posi==1):
                theta1=-111.07
                theta2=-58.6
            elif(posi==2):
                theta1=-54.0
                theta2=-0.5
            elif(posi==3):
                theta1=0.5
                theta2=54.0
            else:
                theta1=58.6
                theta2=111.07
            CC['loc_phi_1']=phi-dphi/2
            CC['loc_phi_2']=phi+dphi/2
            CC['loc_theta_1']=theta1*np.pi/180.0
            CC['loc_theta_2']=theta2*np.pi/180.0
            CC_info[name]=CC
    return CC_info    

def get_OUT_MODE_MATRIX_select_coils(coil_list): ### get the output matrix based on selected sensors
    def effective_fn(theta1,theta2,phi1,phi2,m,n,R,r):
        theta=(theta1+theta2)/2.0
        phi=(phi1+phi2)/2.0
        dtheta=theta2-theta1
        dphi=phi2-phi1
        rlt=r*R*2*np.sin(m*dtheta/2)*np.sin(m*theta-n*phi)/m
        rlt+=(r**2)*np.sin((m+1)*dtheta/2)*np.sin((m+1)*theta-phi)/(m+1)
        rlt+=(r**2)*np.sin((m-1)*dtheta/2)*np.sin((m-1)*theta-phi)/(m-1)
        rlt*=2*np.sin(dphi/2)
        return rlt
    def norm_fn(theta1,theta2,phi1,phi2,R,r):
        dphi=phi2-phi1
        dtheta=theta2-theta1
        rlt=r*R*dphi*dtheta
        rlt+=(r**2)*dphi*(np.sin(theta2)-np.sin(theta1))
        return rlt        
    ### define the local function to calculate the factor modification
    def get_cc_value(theta1,theta2,phi1,phi2):
        R=0.92
        r=0.16
        eff=1.0*effective_fn(theta1,theta2,phi1,phi2,3,1,R,r)
#        eff+=0.8*effective_fn(theta1,theta2,phi1,phi2,2,1,R,r)
#         eff+=0.8*effective_fn(theta1,theta2,phi1-30.0/36*np.pi,phi2-30.0/36*np.pi,4,1,R,r)
#         eff+=0.8*effective_fn(theta1,theta2,phi1-8.0/36*np.pi,phi2-8.0/36*np.pi,5,1,R,r)
        base=norm_fn(theta1,theta2,phi1,phi2,R,r)
        return eff/base
    #############
    cc_info=get_CCInfo()
    ccs=[]
    total = len(coil_list)
    out_mtx=np.zeros([total,4*2])
    cnt=0;
    coils_used_list = []
    for i in xrange(4): # For each toroidal array
#         ccs.extend([ 'FB%02d_C%d' % (j,i+1) for j in np.arange(10)+1 ])
        ccs=[ 'FB%02d_C%d' % (j,i+1) for j in np.arange(10)+1 ]        
        for name in ccs:
            if name in coil_list:
                coils_used_list.append(name)
                cc=cc_info[name]
                theta1=cc['loc_theta_1']
                theta2=cc['loc_theta_2']
                phi1=cc['loc_phi_1']
                phi2=cc['loc_phi_2']
                #### radial and poloidal sensors will apply the same mode shape for now.
                ### the 0,2 is cos,  1,3 is sin,   consistent with the MODE_IN_MTX
                ### get_cc_value calculate the sin value
                out_mtx[cnt,2*i+1]=get_cc_value(theta1,theta2,phi1,phi2)
      #          out_mtx[cnt,2*i+9]=get_cc_value(theta1,theta2,phi1,phi2)
                ### cosx = sin(x +pi/2)  equivalent to  (phi-pi/2)
                out_mtx[cnt,2*i]=get_cc_value(theta1,theta2,phi1-np.pi/2,phi2-np.pi/2)
      #          out_mtx[cnt,2*i+8]=get_cc_value(theta1,theta2,phi1-np.pi/2,phi2-np.pi/2)
                cnt+=1
    return (out_mtx,coils_used_list)


def get_OUT_MODE_MATRIX():
    def effective_fn(theta1,theta2,phi1,phi2,m,n,R,r):
        theta=(theta1+theta2)/2.0
        phi=(phi1+phi2)/2.0
        dtheta=theta2-theta1
        dphi=phi2-phi1
        rlt=r*R*2*np.sin(m*dtheta/2)*np.sin(m*theta-n*phi)/m
        rlt+=(r**2)*np.sin((m+1)*dtheta/2)*np.sin((m+1)*theta-phi)/(m+1)
        rlt+=(r**2)*np.sin((m-1)*dtheta/2)*np.sin((m-1)*theta-phi)/(m-1)
        rlt*=2*np.sin(dphi/2)
        return rlt
    def norm_fn(theta1,theta2,phi1,phi2,R,r):
        dphi=phi2-phi1
        dtheta=theta2-theta1
        rlt=r*R*dphi*dtheta
        rlt+=(r**2)*dphi*(np.sin(theta2)-np.sin(theta1))
        return rlt        
    ### define the local function to calculate the factor modification
    def get_cc_value(theta1,theta2,phi1,phi2):
        R=0.92
        r=0.16
        eff=1.0*effective_fn(theta1,theta2,phi1,phi2,3,1,R,r)
#        eff+=0.8*effective_fn(theta1,theta2,phi1,phi2,2,1,R,r)
#         eff+=0.8*effective_fn(theta1,theta2,phi1-30.0/36*np.pi,phi2-30.0/36*np.pi,4,1,R,r)
#         eff+=0.8*effective_fn(theta1,theta2,phi1-8.0/36*np.pi,phi2-8.0/36*np.pi,5,1,R,r)
        base=norm_fn(theta1,theta2,phi1,phi2,R,r)
        return eff/base
    #############
    cc_info=get_CCInfo()
    ccs=[]
    out_mtx=np.zeros([40,4*2])
    cnt=0;
    for i in xrange(4): # For each toroidal array
#         ccs.extend([ 'FB%02d_C%d' % (j,i+1) for j in np.arange(10)+1 ])
        ccs=[ 'FB%02d_C%d' % (j,i+1) for j in np.arange(10)+1 ]        
        for name in ccs:
            cc=cc_info[name]
            theta1=cc['loc_theta_1']
            theta2=cc['loc_theta_2']
            phi1=cc['loc_phi_1']
            phi2=cc['loc_phi_2']
            #### radial and poloidal sensors will apply the same mode shape for now.
            ### the 0,2 is cos,  1,3 is sin,   consistent with the MODE_IN_MTX
            ### get_cc_value calculate the sin value
            out_mtx[cnt,2*i+1]=get_cc_value(theta1,theta2,phi1,phi2)
  #          out_mtx[cnt,2*i+9]=get_cc_value(theta1,theta2,phi1,phi2)
            ### cosx = sin(x +pi/2)  equivalent to  (phi-pi/2)
            out_mtx[cnt,2*i]=get_cc_value(theta1,theta2,phi1-np.pi/2,phi2-np.pi/2)
  #          out_mtx[cnt,2*i+8]=get_cc_value(theta1,theta2,phi1-np.pi/2,phi2-np.pi/2)
            cnt+=1
    return out_mtx

def get_Bad_Sensorlist():
    filedir=os.path.dirname(__file__ )
    BlacklistFileName="sensor_blacklistQian.txt"
    #filename="%s/../data/%s" % (filedir,BlacklistFileName)
    filename="/home/john/Dropbox/Columbia/Research/sensor_blacklistQian.txt"
    with open(filename,'r') as fh:
        NameBlacklist=[]
        #skip comments
        for line in fh:
            if (line.split()[0]=='Blacklist'):
                break
        for line in fh:
            if line[0]=='#':
                continue
            NameBlacklist.append(line.split()[0])
    return NameBlacklist

if __name__ == "__main__":
    inmtx=get_IN_MODE_MATRIX()
    outmtx=get_OUT_MODE_MATRIX()
    cc_info=get_CCInfo()
    ccs=[]
    for i in xrange(4): # For each toroidal array
        ccs.extend([ 'FB%02d_C%d' % (j,i+1) for j in np.arange(10)+1 ])
    theta=np.zeros(len(ccs))
    phi=np.zeros(len(ccs))
    cnt=0;
#    for name in ccs:
#        cc=cc_info[name]
#        theta[cnt]=(cc['loc_theta_1']+cc['loc_theta_2'])/2.0
#        phi[cnt]=(cc['loc_phi_1']+cc['loc_phi_2'])/2.0
#        cnt+=1;
    data=np.reshape(outmtx[:,1],[4,10])
    plt.contour(data)
    plt.show()
    
        
          
        
        
        
    
        
