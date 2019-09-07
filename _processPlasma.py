"""
this library contains a number of useful functions that are general purpose 
data processing functions.
"""
            
###############################################################################
### import libraries
            
# common libraries
import numpy as _np
import os # X11 PRotection
if os.environ.has_key('DISPLAY'): 
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
        vth=np.sqrt(2*temperatureInEV*eV/mi)
        
        # density
        return 4*np.abs(ionSatCurrent)/q/probeArea/vth
        
        
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
        
