"""
this library contains a number of useful functions that are plasma physics 
related functions, models, and analysis.
"""
            
###############################################################################
### import libraries
            
# common libraries
import numpy as _np
import os # X11 PRotection
#if os.environ.has_key('DISPLAY'):  # for headless operation
import matplotlib.pyplot as _plt
import _plotTools as _plot
import scipy.sparse as sp

# hbtepLib library
#import _plotTools as _plot
            
###############################################################################
### misc functions
            
def _findNearestForWeighting(array,value):
    """
    search through array and returns the index of the cell closest to the value.
    in addition, returns if the found value is larger (+1) or smaller (-1) than the provided value
    ref: http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    idx = (_np.abs(array-value)).argmin()
    if value>array[idx]:
        sign=+1;
    else:
        sign=-1;
    
    return idx,sign #array[idx] #index instead of 
    
###############################################################################
### Langmuir probe calculation
	

def langmuirProbeSimulation(V=_np.arange(-150,151,1),A_probe=0.00032258,V_plasma=50,T_elec=30,T_ion=30,density=1e18,ionMassNumber=2.014102,plot=True):
	"""
	Produces a langmuir probe I-V plot based on several known values
	
	Parameters
	----------
	V : numpy.ndarray
		voltage array from some negative value to some positive value
	A_probe : float
		probe area in m^2
	V_plasma : float
		plasma voltage in volts
	T_elec : float
		electron temperature in eV
	T_ion : float
		ion temperature in eV
	density : float
		plasma density
	ionMassNumber : float
		atomic mass number of the ion.  deuterium = 2.014102
		
	Returns
	-------
	I : numpy.ndarray
		probe current from both electrons and ions
		
	References
	----------
	https://doi.org/10.1119/1.2772282
	
	Example
	-------
	:
		dV=1
		V=np.arange(-150,150+dV,dV)
		a=langmuirProbeSimulation(V)
	"""
	import numpy as np

	# constants
	eV=1.60218e-19;	# eV
	q=1.6e-19;		# fundamental charge
	amu=1.66054e-27; 	# 1 amu to kg
	m_elec=9.109e-31; 	# mass of an electron
	
	# convert temperatures from eV to Joules
	T_elec*=eV
	T_ion*=eV
	
	# ions
	if True:
		m_ion=ionMassNumber*amu
		v_ion_thermal=np.sqrt(8*T_ion/(np.pi*m_ion))
		if T_elec > T_ion*5:
			I_ion_sat=0.6*q*density*np.sqrt(T_elec/m_ion)*A_probe
		else:
			I_ion_sat=0.25*q*density*v_ion_thermal*A_probe
		I_ion=np.zeros(len(V))
		I_ion[V<V_plasma]=-I_ion_sat
		I_ion[V>=V_plasma]=-I_ion_sat*np.exp(q*(V_plasma-V[V>=V_plasma])/(T_ion))
		
	# electrons
	if True:
		v_elec_thermal=np.sqrt(8*T_elec/(np.pi*m_elec))
		I_elec_sat=0.25*q*density*v_elec_thermal*A_probe
		I_elec=np.zeros(len(V))
		I_elec[V<V_plasma]=I_elec_sat*np.exp(-q*(V_plasma-V[V<V_plasma])/(T_elec))
		I_elec[V>=V_plasma]=I_elec_sat
		
	# total current
	I=I_elec+I_ion
		
	if plot==True:
		import matplotlib.pyplot as plt
		fig,ax=plt.subplots()
		ax.plot(V,I_ion,label="Ion current")
		ax.plot(V,I_elec,label="Elec. current")
		ax.plot(V,I,label="Total current")
		_plot.finalizeSubplot(ax,xlabel='Bias voltage (V)',ylabel='Probe current (A)')
		
	return I
    

class langmuirProbe:
    """
    Langmuir probe analysis.  Returns temperature and density of a plasma 
    with a provided I-V profile. 
    The code fits an exponential fit between the floating potential and the 
    plasma potential in order provide the temperature and density.  
    
    Parameters
    ----------
    V : np.array() of floats
        Voltage of probe.  Units in volts.
    I : np.array() of floats
        Current of probe.  Units in amps.
    expRegionMinVoltage : float or NoneType
        Minimum potential to be used in the exponential fit
        if None, you will given a plot of the I-V profile and asked for the value
    expRegionMaxVoltage : float or NoneType
        Maximum potential to be used in the exponential fit
        if None, you will given a plot of the I-V profile and asked for the value
    ionSatRegionMaxVoltage : float or NoneType
        Maximum voltage value to consider for the ion saturation region  
        if None, you will given a plot of the I-V profile and asked for the value
    area : float
        probe area in cubic meters.  
        HBT-EP's BP: 0.75 inch diameter, half sphere ->  
                     4*pi*(0.0254*.75/2)^2/2 = 0.000580644 m^3
    plot : bool
        True = plots final fit 
        False = does not
    expFitGuess : tuple of three floats
        provide a guess for the exponential fit
        expFitGuess = (a,b,c) 
        where the function is a*_np.exp(b*x)+c
        a = the amplitude
        b = exp const
        c = y-offset
        expFitGuess = (6, 0.05, -5) by default 
        
    Attributes
    ----------
    #TODO(John) add attributes
    
    Notes
    -----
    Using the exponential offset from the current data gives an uncomfortably 
    large density value.  This is probably not the best way to do this.  
        
    """
    
    
    def __init__(self, V, I, expRegionMinVoltage=None, expRegionMaxVoltage=None, ionSatRegionMaxVoltage=None,
                 area=0.000580644,
                 plot=False, expFitGuess=(6, 20, -5)):
        ## physical constants
#        eV=1.60218e-19;
#        mi=1.6737236 * 10**(-27) * 2
#        q=1.6e-19
        
        # parameters
        self.probeArea=area
        self.expRegionMinVoltage=expRegionMinVoltage
        self.expRegionMaxVoltage=expRegionMaxVoltage
        self.ionSatRegionMaxVoltage=ionSatRegionMaxVoltage
        
        # ensure V and I are arrays
        self.V=_np.array(V); 
        self.I=_np.array(I)
        
        # sort V in ascending order
        i=_np.argsort(self.V)
        self.V=self.V[i]
        self.I=self.I[i]
        
        # initialize plot
        p1=_plot.plot()
        p1.addTrace(xData=self.V,yData=self.I,marker='.',linestyle='',
                         yLegendLabel='raw data')
        
        # if expRegionMinVoltage or expRegionMaxVoltage were not specified, the code will plot the I-V
        # profile and ask that you provide the floating and/or plasma 
        # potential.  These values are used as the lower and upper limits for
        # the exponential fit
        if expRegionMinVoltage==None or expRegionMaxVoltage==None or ionSatRegionMaxVoltage==None:
            p1.plot()
            _plt.show()
            _plt.pause(1.0) # a pause is required for the plot to show correctly
            if expRegionMinVoltage==None:
                self.expRegionMinVoltage=float(raw_input("Please provide the approximate lower voltage (units in volts) limit to be used in the exp fit by looking at the I-V plot:  "))
            if expRegionMaxVoltage==None:
                self.expRegionMaxVoltage=float(raw_input("Please provide the approximate upper voltage (units in volts) limit to be used in the exp fit by looking at the I-V plot:  "))
            if ionSatRegionMaxVoltage==None:
                self.ionSatRegionMaxVoltage=float(raw_input("Please provide the approximate maximum voltage (units in volts) limit for the ion saturation region by looking at the I-V plot:  "))
        
        # exp. curve fit setup
        from scipy.optimize import curve_fit

        # perform exp curve fit
        popt, pcov = curve_fit(self._exponenial_func, V, I, p0=expFitGuess)
        self.expFitParameters=popt
#        print popt
        
        # temperature
        self.temperatureInEV=self.calcTempInEV(popt[1])#q*popt[1]/eV
        print("Temperature = " + str(self.temperatureInEV) + ' eV')
        
        # density calculation from the exp fit offset
        self.densityFromExpFitOffset=self.calcDensity(popt[2],
                                                      probeArea=self.probeArea,
                                                      temperatureInEV=self.temperatureInEV)
        print("Density from the exp. fit offset current = " + str(self.densityFromExpFitOffset) + ' m^3')
        
        # density calculation from averaging the values in the ion sat region
        i=self.V<self.ionSatRegionMaxVoltage
        aveIonSatCurrent=_np.average(self.I[i])
        print("Average Current in the ion sat. region = " + str(aveIonSatCurrent))
        self.densityFromAveIonSatRegion=self.calcDensity(aveIonSatCurrent,
                                                      probeArea=self.probeArea,
                                                      temperatureInEV=self.temperatureInEV)
        print("Density from the average current in the ion sat. region = " + str(self.densityFromAveIonSatRegion) + ' m^3')
        
        # optional plot
        if plot==True:
            self.plot()
            
            
    def calcTempInEV(self, expFitCoeffWithVoltUnits):
        """
        Calulates temperature from langmuir exp fit
        
        Parameters
        ----------
        expFitCoeffWithVoltUnits : float
            
        """
        # constants
        eV=1.60218e-19;
        q=1.6e-19
        
        # temperature in eV
        return q*expFitCoeffWithVoltUnits/eV
        
        
    def calcDensity(self, ionSatCurrent, probeArea, temperatureInEV):
        """
        Calulates density
        
        Parameters
        ----------
        expFitCoeffWithVoltUnits : float
            
        """
        # constants
        eV=1.60218e-19;
        q=1.6e-19
        mi=1.6737236 * 10**(-27) * 2
        
        # thermal velocity
        vth=_np.sqrt(2*temperatureInEV*eV/mi)
        
        # density
        return 4*_np.abs(ionSatCurrent)/q/probeArea/vth
        
        
    def plot(self):
        """ plot raw data and exp fit """
        
        # exp fit
        xFit=_np.arange(self.expRegionMinVoltage,self.expRegionMaxVoltage,0.1)
        yFit=self._exponenial_func(xFit,self.expFitParameters[0],
                                  self.expFitParameters[1],
                                  self.expFitParameters[2])
        
        # extrapolated exp fit 
        xFitExtrap=_np.arange(_np.min(self.V),_np.max(self.V),0.1)
        yFitExtrap=self._exponenial_func(xFitExtrap,self.expFitParameters[0],
                                  self.expFitParameters[1],
                                  self.expFitParameters[2])
        
        # generate plot
        p1=_plot.plot(title='I-V Profile',xLabel='Probe Voltage [V]',
                      yLabel='Probe Current [A]')
        p1.addTrace(xData=self.V,yData=self.I,marker='.',linestyle='',
                         yLegendLabel='Raw data',alpha=0.5)
        p1.addTrace(xData=xFit,yData=yFit, yLegendLabel='Exp fit')
        p1.addTrace(xData=xFitExtrap,yData=yFitExtrap, 
                    yLegendLabel='Extrapolated exp fit',linestyle=':')
        p1.plot()
        
        return p1
        
        
    def _exponenial_func(self,x, a, b, c):
        """ exponential function """
        return a*_np.exp(x/b)+c
            
    
    
###############################################################################
### PIC code

class picCode:
    """
    PIC code solver for plasma

    Parameters
    ----------
    N : int
        number of particles to be created.  default = 128
    dt : float
        time step.  default = 0.1
    writeSteps : int
        writes data every dtWrite iterations.  default = 10
    tEnd : float
        time to end calculations.  default = 8*pi
    M : int
        number of grid points for the various fields.  default = 128
    L : float
        x-axis spatial domain length.  default = 2*pi
    vxBounds : list of two floats
        bounds for v_x for plotting purposes.  does not impose any limit on
        calculations.  default = [-5,5]
    Bz : float
        z-axis magnetic static magnetic field.  
    xInit : str or np.array(N) of floats
        'uniform' - positions particles uniformly across the x-axis spatial
                    domain (default)
        'random' - positions particles randomly across the x-axis spatial
                    domain (default)
    vxInit : str or np.array(N) of floats
        'norm' - normal distribution with v0 as the std. dev. (default)
        'sin' - n=1 sine wave distribution with v0 as the amplitude
        '2stream' - two stream distribution with v0 as the amplitude
        'ring' - ring distribution with v0 as the radius
    vyInit : str or np.array(N) of floats
        'norm' - normal distribution with v0 as the std. dev. (default)
        'cos' - n=1 cosine wave distribution with v0 as the amplitude
        '2stream' - two stream distribution with v0 as the amplitude
        'ring' - ring distribution with v0 as the radius
    charge : str or np.array(N) of floats
        'allelectrons' - all negative charge
        'allions' - all positive charge
        'mix' - even/odd mix of particles
    plot : bool
        causes the results to be plotted
    titleAdendum : str
        a string that is added to the figure title
    order : 
        order of charge weighting and force weighting
    qOverM : float
        ratio of charge (q) to mass (M) for the particles.  default = 1
    v0 : float
        amplitude of initial velocity if using a default distribution
    numBins : int
        number of bins for distribution function histogram.  default = 20
    
    Attributes
    ----------
    #TODO(John) fill in attributes
    """
    
    def chargeWeighting(self,x,q,Xj,order=1):
        """
        charge weighting methods, 1st and 0th order.
        converts charge distribution to scalar rho.
        
        Parameters
        ----------
        x : np.array
            array of particle locations
        q : np.array
            the charge of each particle
        Xj : np.array
            x-coordinate for the center of each grid
        order : int
            order of field weighting algorithm
            order = 0 is picks the x-coordinate that is nearest
            order = 1 (default) uses linear inerpolation
        """
        M=len(Xj);
        N=len(x);
        rho=_np.zeros(M);
        
        if order == 1:
            dx=Xj[1]-Xj[0];
            for i in range(0,N):
                idx,sign=_findNearestForWeighting(Xj,x[i])
                if idx+sign >= M:
                    rho[idx]+=q[i]*_np.abs(Xj[idx]+dx-x[i])/dx;
                    rho[0]+=q[i]*_np.abs(Xj[idx]-x[i])/dx;
                elif idx+sign < 0:
                    rho[0]+=q[i]*_np.abs(Xj[idx]-dx-x[i])/dx;
                    rho[-1]+=q[i]*_np.abs(Xj[idx]-x[i])/dx;
                else:
                    rho[idx]+=q[i]*_np.abs(Xj[idx+sign]-x[i])/dx;
                    rho[idx+sign]+=q[i]*_np.abs(Xj[idx]-x[i])/dx;
                        
        if order == 0:
            for i in range(0,N):
                [idx,sign]=_findNearestForWeighting(Xj,x[i])
                rho[idx]+=q[i]
                
        # enforce quasi-neutrality even if n_e != n_i
        rho-=_np.average(rho)  
        return rho/dx
    #    return rho
        
        
    def fieldWeighting(self,x, Xj, E, order=1):
        """
        1st and 0th order force (field) weighting methods.
        Converts electric field at each grid location to the electric field
        at each particle location.
        
        Parameters
        ----------
        x : np.array
            array of particle locations
        Xj : np.array
            x-coordinate for the center of each grid
        E : np.array
            electric field at each grid location
        order : int
            order of field weighting algorithm
            order = 0 is picks the x-coordinate that is nearest
            order = 1 (default) uses linear inerpolation
        """
        ## interpolate Ej to get Fi
        N=len(x);
        if order == 1:
            Ei=_np.interp(x,Xj,E)
            
        elif order == 0:
            Ei=_np.zeros(N);
            for i in range(0,N):
                [idx,sign]=_findNearestForWeighting(Xj,x[i]);
    #            Ei[i]=q[i]*E[idx];
                Ei[i]=E[idx];
        return Ei
    
    
    def __init__(self,N=128,dt=0.1,writeSteps=10,tEnd=8*_np.pi,M=128,
                 L=2*_np.pi,vxBounds=[-5,5],Bz=0,xInit='uniform',vxInit='norm',
                 vyInit=1,charge='allions',plot=False,titleAdendum='',order=1,
                 qOverM=1, v0=None, numBins=20):
        
        ## initialize
        self.N=N;
        self.dt=dt;

        self.tEnd=tEnd;
        if M==None:
            self.M=N;
        else:
            self.M=M;
        self.vxBounds=vxBounds;
        self.qOverM=1;
        self.L=L 
        self.xBounds=[-0.5*L,0.5*L];
        numBins=numBins;

        
        ## variables
        # time
        self.time=_np.arange(0,tEnd+dt,dt);
        
        # title
        self.title=titleAdendum+'. N='+str(N)+'. M='+str(M)+". dt=" + "%.3f" % dt+'. Order=' + str(order) +'. L=' + "%.3f" % L +'.' 
        
        # distribution function 
        self.fv=_np.zeros([numBins,len(self.time)]);
        
        # distribution function x-axis
        self.fvX=_np.zeros(numBins+1);
        
        # x-axis spatial domain for each particle
        self.x=_np.zeros([self.N,len(self.time)])
        
        # x-axis velocity for each particle
        self.vx=_np.zeros([self.N,len(self.time)])
        
        # y-axis velocity for each partile
        self.vy=_np.zeros([self.N,len(self.time)])
        
        # charge distribution  at each grid location
        self.rho=_np.zeros([self.M,len(self.time)])
        
        # electric field at each grid location
        self.Ex=_np.zeros([self.M,len(self.time)])
        
        # electric field at each particle
        self.Ei=_np.zeros([self.N,len(self.time)])
        
        # z-axis magnetic field
        self.Bz=Bz;
        
        # gyroradius
        self.wc=Bz*self.qOverM;
        
        # electric potential at each grid location
        self.phi=_np.zeros([self.M,len(self.time)]) 
        
        ## measured quantities
        # electric-field potential energy
        self.EE=_np.zeros(len(self.time))
        
        # potential energy        
        self.KE=_np.zeros(len(self.time)) 
    
        # x coordinate initialization
        if xInit=='uniform':
            self.x[:,0]=_np.random.uniform(high=self.xBounds[1],low=self.xBounds[0],size=N);
        elif isinstance(xInit,_np.ndarray) or isinstance(xInit,list):
            self.x[:,0]=xInit;
        elif xInit == 'even':
            dxN=self.L/N;
            self.x[:,0]=_np.linspace(self.xBounds[0]+dxN/2.,self.xBounds[1]-dxN/2,N) # center x coordinate of each grid cell
            
        # vx coordinate initialization
        if vxInit=='norm':
            self.vx[:,0]=_np.random.normal(loc=0,scale=2.0/2.35482,size=N);  #vFWHM=2.0
        elif isinstance(vxInit,_np.ndarray) or isinstance(vxInit,list):
            self.vx[:,0]=vxInit;
        elif vxInit=='sin':
            self.vx[:,0]=v0*_np.sin(self.x[:,0]*2*_np.pi/self.L)
        elif vxInit=='2stream':
            self.vx[::2,0]=v0[0]+v0[1]*_np.sin(self.x[::2,0]*2*_np.pi/self.L)
            self.vx[1::2,0]=-v0[0]-v0[1]*_np.sin(self.x[1::2,0]*2*_np.pi/self.L)
        elif vxInit=='ring':
            theta=_np.random.uniform(high=2*_np.pi,low=0,size=N)
            self.vx[:,0]=v0*_np.cos(theta);
            
        # vy coordinate initialization
        if vyInit=='cos':
            self.vy[:,0]=v0*_np.cos(self.x[:,0]*2*_np.pi/self.L)
        elif isinstance(vyInit,_np.ndarray) or isinstance(vyInit,list):
            self.vy[:,0]=vyInit;
        elif vyInit=='ring':
            self.vy[:,0]=v0*_np.sin(theta);
            
        if N >=100:
            a=_np.histogram(self.vx[:,0],numBins,range=vxBounds)
            self.fv[:,0]=a[0];
            self.fvX=a[1];
            
        # charge initialization
        q=_np.zeros((N,))
        if charge == 'mix':
            q[::2]=1
            q[1::2]=-1
        elif charge == 'allions':
            q+=1;
        elif charge == 'allelectrons':
            q-=1;
        elif isinstance(q,_np.ndarray) or isinstance(q,list):
            q=charge

        # grid cell spacing and dimensioning
        self.dx=self.L/M;
        dx2=self.dx**2;
        oneOver2Dx=1./(2.*self.dx)
        self.Xj=_np.linspace(self.xBounds[0]+self.dx/2.,self.xBounds[1]-self.dx/2,M) # center x coordinate of each grid cell

        ## initial calculations
        # calculate charge
        self.rho[:,0]=self.chargeWeighting(self.x[:,0],q,self.Xj,order=order);
        
        # calculate Phi, electric potential
        A=sp.diags([1, -2, 1], [0,1,2], shape=(M-1, M)).toarray()
        A[M-2,0]=1;
        Ainv=_np.linalg.pinv(A);
        Ainv=sp.csc_matrix(Ainv)
        self.phi[:,0]=-dx2*Ainv.dot(self.rho[1:M,0])
        
        # calculate E, electric field
        ExFDMatrix=sp.diags([-1, 0, 1], [-1, 0, 1], shape=(M, M)).toarray()
        ExFDMatrix[0,M-1]=-1
        ExFDMatrix[M-1,0]=1
        ExFDMatrix=sp.csc_matrix(ExFDMatrix)
        self.Ex[:,0]=-1/(2*self.dx)*ExFDMatrix.dot(self.phi[:,0])
        
        # calculate energy
        self.EE[0]=_np.sum(self.Ex[:,0]**2)*self.dx;
        
        # kinetic energy
        self.KE[0]=0.5*_np.sum(_np.square(self.vx[:,0])+_np.square(self.vy[:,0]))*self.dx
            
        ## evolve in time
        for i in range(1,len(self.time)):
            
            # print time every integer step
            if _np.remainder(self.time[i],1.)==0:
                print("t = " + str(self.time[i]))
            
            ## Solve force
            self.Ei[:,i-1]=self.fieldWeighting(self.x[:,i-1], self.Xj, self.Ex[:,i-1], order=order);
            
            ## evolve v
            # first half v step
            vxprime=self.vx[:,i-1]+dt/2.*self.Ei[:,i-1]
            vyprime=self.vy[:,i-1]
            # Full rotation (above z)
            vxdprime=vxprime*_np.cos(self.wc*dt)+vyprime*_np.sin(self.wc*dt);
            vydprime=-vxprime*_np.sin(self.wc*dt)+vyprime*_np.cos(self.wc*dt);
            # second half v step 
            self.vx[:,i]=vxdprime+dt/2.*self.Ei[:,i-1]
            self.vy[:,i]=vydprime
            
            ## evolve x
            self.x[:,i]=self.vx[:,i]*dt+self.x[:,i-1];        
    
            ## periodic BCs, impose
            self.x[_np.where(self.x[:,i]>self.xBounds[1]),i]=self.x[_np.where(self.x[:,i]>self.xBounds[1]),i]-self.L;  #while np.sum(x>xBounds[1])>0:
            self.x[_np.where(self.x[:,i]<self.xBounds[0]),i]=self.x[_np.where(self.x[:,i]<self.xBounds[0]),i]+self.L;  #while np.sum(x<xBounds[0])>0:
                
            ## solve for grid-dependent fields and scalars
            self.rho[:,i]=self.chargeWeighting(self.x[:,i],q,self.Xj);
            self.phi[:,i]=-dx2*Ainv.dot(self.rho[1:M,i])
            self.Ex[:,i]=-oneOver2Dx*ExFDMatrix.dot(self.phi[:,i])
            
            ## calculate energies
            self.KE[i]=0.5*_np.sum(_np.square(self.vx[:,i])+_np.square(self.vy[:,i]))*self.dx
            self.EE[i]=_np.sum(self.Ex[:,i]**2)*self.dx;
            
            ## velocity distribution
            a=_np.histogram(self.vx[:,i],numBins,range=vxBounds)
            self.fv[:,i]=a[0];
            
        if plot==True:
            self.animateTraj()
            
            
    def plotKE(self):
        """
        plots total kinetic energy with time
        
        Parameters
        ----------
        """
        p1=_plot.plot.plot(title=self.title,xLabel='time',yLabel='KE')
        p1.addTrace(xData=self.time,yData=self.KE)
        p1.plot()
        
    def plotEE(self):
        """
        plots total electric potential energy with time
        
        Parameters
        ----------
        """
        p1=_plot.plot.plot(title=self.title,xLabel='time',yLabel='EE')
        p1.addTrace(xData=self.time,yData=self.EE)
        p1.plot()
        
    def plotSingleX(self,index=0):
        """
        plots a single particle's x-location with time
        
        Parameters
        ----------
        index : int
            index of particle to be plotted
        """
        p1=_plot.plot.plot(title=self.title,xLabel='time',yLabel='x position')
        p1.addTrace(xData=self.time,yData=self.x[index,:])
        p1.plot()
        
    def plotSingleVX(self,index=0):
        """
        plots a single particle's x-velocity with time
        
        Parameters
        ----------
        index : int
            index of particle to be plotted
        """
        p1=_plot.plot.plot(title=self.title,xLabel='time',yLabel='v_x')
        p1.addTrace(xData=self.time,yData=self.vx[index,:])
        p1.plot()
        
    def plotSingleVY(self,index=0):
        """
        plots a single particle's y-velocity with time
        
        Parameters
        ----------
        index : int
            index of particle to be plotted
        """
        p1=_plot.plot.plot(title=self.title,xLabel='time',yLabel='v_y')
        p1.addTrace(xData=self.time,yData=self.vy[index,:])
        p1.plot()
        
    def plotPhase(self,timeIndex=0, v='vx'):
        """
        plots phase (x vs v_x) for a single time
        
        Parameters
        ----------
        timeIndex : int
            time index of particle to be plotted
        v : str
            'vx' - plots x-velocity
            'vy' - plots y-velocity
        """
        p1=_plot.plot.plot(title=self.title+' t='+str(self.time[timeIndex]),xLabel='x',yLabel='v_x')
        p1.addTrace(xData=self.x[:,timeIndex],yData=self.vx[:,timeIndex],marker='.',linestyle='')
        p1.plot()
        
    def animateTraj(self,stepPause=0.01):
        """
        animates all particle trajectories in phase space
        
        Parameters
        ----------
        stepPause : float
            time delay between frames
        """
        fig, ax = _plt.subplots()
        ax.set_xlim(self.xBounds) 
        ax.set_ylim(self.vxBounds) 
        ax.set_xlabel('x')
        ax.set_ylabel(r'$v_{x}$')
        x=self.x[:,0];
        y=self.vx[:,0];
        points, = ax.plot(x, y, marker='o', linestyle='None')
        for i in range(0,len(self.time)): #
            print(self.time[i])
            points.set_data(self.x[:,i], self.vx[:,i])
            _plt.title('t='+str(self.time[i]))
            _plt.pause(stepPause);
            
    def animateFv(self,stepPause=0.01,plotFV0=False):
        """
        animates distribution function in time
        
        Parameters
        ----------
        stepPause : float
            time delay between frames
        plotFV0 : bool
            True - Plots initial FV behind the animation
        """
        _plt.figure()
        _plt.ylabel(r'$f(v_x)$')
        _plt.xlabel(r'$v_{x}$')
        for i in range(0,len(self.time)): #
            _plt.xlim(self.vxBounds)
            _plt.pause(stepPause);
            _plt.cla()
            if plotFV0==True:
                _plt.hist(self.vx[:,0],bins=30,range=self.vxBounds)
            _plt.hist(self.vx[:,i],bins=30,range=self.vxBounds)
            print(self.time[i])
            _plt.title('t='+str(self.time[i]))
    
    
###############################################################################
### Misc. plasma models	
			
def thetaCorrection(	shotno,
					theta=_np.linspace(-_np.pi,_np.pi,100),
					tStart=2e-3,
					tStop=5e-3,
					plot=False):
	
	"""
	This function corrects the theta coordinate (theta) for the non-cylindrical
	nature of the HFS and LFS magnetic fields.  This is based on equation 4.4 
	in Jeff's thesis
	
	Work in progress
	"""
	
	# constants
	mu0=4e-7*_np.pi
	
	# libraries
	from _getHBTData import ipData
	from _getHBTData import capBankData
	import _plotTools as _plot
	import pandas as pd
	
	
	def lambdaCalc(Bv, Ip):
		# Jeff's thesis, equations 4.5 and 4.7
		a=0.15
		R0=0.92
		term=4*_np.pi*R0*Bv/mu0/Ip-_np.log(8*R0/a)+1.5
		L=(term+1)*a/R0
		L[L>1.0]=1.0
		L[L<-1.0]=-1.0
		return L
	
	def thetaStarCalc(theta, L):
		return theta-L*_np.sin(theta)
		
	# get plasma current
	Ip=ipData(shotno,tStart=tStart,tStop=tStop).ip
#	Ip=ipData.ip
#	timeIp=ipData.time
	
	# get cap bank data
	capData=capBankData(shotno,tStart=tStart,tStop=tStop)
	vfCurrent=capData.vfBankCurrent
	ohCurrent=capData.	ohBankCurrent
	time=capData.vfTime
	
	# calculate B fields at R_0 = 0.92m.  (Using static values from Jeff's code)
	Bv_vf=vfCurrent*(-2.6602839e-06)
	Bv_oh=ohCurrent*(1.9580808e-08)*(-1) # -1 so that the field subtracts 
	
	# net field
	Bv=Bv_vf+Bv_oh
	
	# calculate lambda
	L=lambdaCalc(Bv,Ip)
	
	# theta correction
	try:
		m=len(theta)
	except:
		m=1
	thetaStar=_np.zeros((len(L),m))
	for i in range(len(L)):
		thetaStar[i]=thetaStarCalc(theta,L[i]) #theta-L[i]*_np.sin(theta)
		
	# init dataframe
	dfData=pd.DataFrame(	data=thetaStar,
						index=time,
						columns=theta)
	
	if plot==True:
		
		if m!=1:
			fig,ax=_plt.subplots()
			_plot.contourPlot(		ax,
								x=dfData.index*1e3,
								y=dfData.columns,
								z=dfData.transpose(),
								levels=_np.arange(-3,3+0.5,0.5),
								xlabel='Time (ms)',
								ylabel='Theta (rad)',
								zlabel='Theta corrected (rad)',
								zticklabels=_np.arange(-3,3+0.5,0.5),
								ztickLabels=_np.arange(-3,3+0.5,0.5),
								fill=False
								)
		if True:
			fig,ax=_plt.subplots()
			
			ax.plot(time*1e3,L,label='L(t)')
			_plot.finalizeSubplot(	ax,
								xlabel=r'Time (ms)',
								ylabel=r'L(t)',
								)
			
			
		if True:
			fig,ax=_plt.subplots()
			theta=_np.arange(-_np.pi,_np.pi,0.1)
			
			for L in _np.linspace(-1,1,11):
				ts=thetaStarCalc(theta,L)
				
				ax.plot(theta,ts,label=r'$\lambda$=%0.2f'%L)
			_plot.finalizeSubplot(	ax,
								xlabel=r'$\theta$',
								ylabel=r'$\theta^*$',
								)
					

	
			
def currentDensityModel(iP,q_limiter,r,r_limiter=0.15,q_offset=0.9,plot=False,verbose=False):
	"""
	Calculates a tokamak's current density using Wesson's model
	
	Parameters
	----------
	ip : float
		plasma current
	q_limiter : float
		safety factor at r_limiter (minor radius at the limiter)
	r : numpy.ndarray
		radial coordinate in meters.  should range from 0 to r_wall
	r_wall : float
		minor radius at the wall in meters.  r_wall=0.16 in HBT-EP
	r_limiter : float
		minor radius at the limiter in meters.  r_wall=0.15 in HBT-EP
	q_offset : float
		safety factor at r=0
	plot : bool
		plot results
	verbose : bool
		print misc output
		
	Returns
	-------
	j : numpy.ndarray
		current density as a function of minor radius in amps per meter squared
	r : numpy.ndarray
		radial coordinate in meters
		
	Example
	-------
	
	::
		
		r_wall=0.16
		r=np.linspace(0,r_wall,1001)
		currentDensityModel(10e3,3,r_wall=r_wall)
	
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	
	l=q_limiter/q_offset-1
	
	def firstOrderIntegration(x,y):
		""" numerical intergration """
		dx=x[1]-x[0]
		return np.sum(dx*y)
	
	def wessonCurrentModel(r,params):
		"""
		Wesson's current model, page 114 in his 2004 book
		"""
		j0=params[0] # j(r=0)=j0
		r0=params[1] # plasma edge (last closed flux surface)
		l=params[2]
		j=j0*(1-(r/r0)**2)**(l)
		j[np.where(r>r0)]=0
		return j
	
	def calcCurrentProfileFromIP(r,r_limiter,radialFunction,params,iP,j0GuessLeft=1e5,j0GuessRight=1e7,j0Guess=1e6,errorTol=1e-6,plot=True,verbose=True):
		"""
		The references only provide I_P and do not provide j(r=0).  This 
		subfunction makes a guess at j(0) and calculates j(r) with the provided
		q-profile function.  It then iterates until the integral is equal to IP.  
		
		Parameters
		----------
		r : numpy.array
			radial coordinate array
		r_limiter : float
			radial location of the limiter
		radialFunction : function(r,params)
			returns radial density current distribution
		params : list
			list of parameters to pass to radialFunction
		iP : float
			plasma current [amps]
		j0GuessLeft : float
			lower bound of j0 guess value
		j0GuessRight : float
			upper bound of j0 guess value
		errorTol : float
			error tolerance to end iteration
			
		Return
		------
		j : np.array
			radial current density where it's intergal = iP
		
		References
		----------
		http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
		"""
		
		j=radialFunction(r,[j0Guess,params[1],params[2]])
		ITotal=firstOrderIntegration(r,j*r)*2*np.pi
		error=(ITotal-iP)/iP
		
		count=0
		if verbose==True:
			print('Starting iterative solver to calculated current density given the plasma current')
		while(np.abs(error)>errorTol):
			count+=1
		
			if error<0:
				j0GuessLeft=j0Guess
			else:
				j0GuessRight=j0Guess
				
			j0Guess=np.mean([j0GuessLeft,j0GuessRight])
			
			j=radialFunction(r,[j0Guess,params[1],params[2]])
			
			ITotal=firstOrderIntegration(r,j*r)*2*np.pi
			error=(ITotal-iP)/iP
			if verbose==True:
				print('count: %d, \t error: %.6f \t guess: %.3f, \t I: %.1f' % (count,error,j0Guess,ITotal))
	
		if plot==True:
			fig,ax=plt.subplots()
			ax.plot(r*100.,j/10000.)
			_plot.finalizeSubplot(ax,xlabel='minor radius (cm)',ylabel=r'current density ($A/cm^2$)')
			plt.show()
		return j
	
	j=calcCurrentProfileFromIP(r,r_limiter=r_limiter,iP=iP,
						   radialFunction=wessonCurrentModel, 
						   params=[1,r_limiter,l],j0Guess=2783578.873,plot=plot,verbose=verbose)
	
	return j

def qProfile_cylindricalApproximation(r,j,iP,r_limiter=0.15,R=0.92,BT=0.35,q_limiter=3.0,q_offset=0.9,plot=False):
	"""
	Recommended q-profile model in Ivanov's 2014 paper. 
	
	Parameters
	----------
	r : numpy.ndarray
		minor radial coordinate in meters.  should range from 0 to r_wall
	j : numpy.ndarray
		current density (amps per meter squared) as a function of radius
	iP : float
		plasma current (in amps)
	r_limiter : float
		radial location of the limiter
	R : float
		major radius (in meters)
	BT : float
		toroidal magnetic field strength (in Tesla)
	q_limiter : float
		safety factor at the limiter
	q_offset : float
		safety factor at r=0
	
	Example
	-------
		import numpy as np
		r_wall=0.16
		iP=14e3
		r=np.linspace(0,r_wall,1001)
		q_offset=0.9
		r_limiter=0.15
		q_limiter=3.0
		j=currentDensityModel(iP=iP,r=r,r_wall=r_wall,q_offset=q_offset,r_limiter=r_limiter,q_limiter=q_limiter)
		q=qProfile_cylindricalApproximation(r,j,iP=iP,r_limiter=r_limiter,q_offset=q_offset,q_limiter=q_limiter,plot=True)
		
	Notes
	-----
	The original source for q(r) is only valid for r<=a.  
	To correct for this, I solved \int B_{\theta} dl = \mu I_p and 
	q=\frac{rB_z}{RB_{\theta}} to provide q all the way out to r=b.
	
	References
	----------
	https://doi.org/10.1063/1.4897174
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	
	## physical constants
	mu0=4*np.pi*1e-7
	
	l=q_limiter/q_offset-1
	q=  2*(l+1)*BT/(mu0*j[0]*R)*(r/r_limiter)**2/(1-(1-(r/r_limiter)**2)**(l+1))
	q[0]=q[1]
	i=np.where(q>0)[0]
	for k in range(i[-1]+1,len(q)):
		q[k]=2*np.pi*r[k]**2*BT/(R*mu0*iP)
		
	if plot==True:
		fig,ax=plt.subplots(2)
		ax[0].plot(r*100.,j/10000.)
		_plot.finalizeSubplot(ax[0],ylabel=r'$A/cm^2$',subtitle='Current density')
		ax[1].plot(r*100,q)
		_plot.finalizeSubplot(ax[1],xlabel='minor radius (cm)',ylabel=r'q',subtitle='Safety factor')
		_plot.finalizeFigure(fig,figSize=[6,6/1.6])
		plt.show()
	return q,j,r
	