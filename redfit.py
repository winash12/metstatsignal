import pandas as pd
import numpy as np
import math
from scipy.optimize import brent,minimize_scalar
from scipy import signal,special
import sys
import traceback
from scipy.stats  import chisquare,chi2,ttest_ind,sem
from pymicra.signal import test_reverse_arrangement
import matplotlib.pyplot as plt

def main():
    
    global freq,varx,nsim
    data = np.loadtxt("x.dat",skiprows=0)
    data1 = np.loadtxt("y.dat",skiprows=0)
    t = data[:,0] 
    x = data[:,1]

    y = data1[:,1]

    #tt,xx= testAWSData()
    gxx = createLombScargleSpectrum(t,x)
    #fig,ax = plt.subplots()

    #ax.plot(freq[1:],gxx)
    
    #plt.savefig('normalized_periodogram.png')
    
    #plt.show()

    varx = freq[1]*np.sum(gxx)

    tau = gettau(t,x)

    grxx = performMCSimulation(t,tau)

    grxxavg = np.mean(grxx,axis=0)

    varrx = freq[1]*np.sum(grxxavg)
    facx = varx/varrx

    grxxavg = facx*grxxavg

    avgdt = np.sum(np.diff(t))/t.size
    rho = np.exp(-avgdt/tau)
    rhosq = rho*rho
    freq_ang=2.*np.pi*freq
    gredthx = (1.0 -rhosq)/(1.0-2.0*rho*np.cos(np.pi*freq_ang[1:]/freq_ang[-1])+rhosq)
    varrx = freq[1]*np.sum(gredthx)
    facx = varx/varrx

    gredthx = facx*gredthx
    corrx = grxxavg/gredthx

    gxxc = gxx/corrx


    np.sort(grxx,axis=1)

    idx90 = int(0.90* nsim)
    idx95 = int(0.95* nsim)
    idx99 = int(0.99* nsim)

    ci90 = grxx[idx90,:]/corrx
    ci95 = grxx[idx95,:]/corrx
    ci99 = grxx[idx99,:]/corrx


    statistic,p_value=ttest_ind(gxxc,grxx)

    dof = getdof(1,7)
    print(dof)

    try:
        chi_square_test_statistic,p_value = chisquare(gxxc,gredthx)
    except:
        pass
    print(chi2.ppf(1-0.10,dof))
    print(chi2.ppf(1-0.05,dof))
    print(chi2.ppf(1-0.01,dof))
    alphacritx= 1.0/75
    faccritx = chi2.isf(alphacritx,dof)
    print(faccritx)
    test_reverse_arrangement(gxxc,points_number=gxxc.size,alpha=0.1)

def testAWSData():

    data = pd.read_csv('N 2021-2022-ICHENN7-WestMambalamPWS.csv',parse_dates=[['Date','Time']],dayfirst=True)
    
    times = pd.to_datetime(data['Date_Time'],unit='ms')

    
    minutes = times.view(np.int64)/6E10


    df1 = (data[['Pressure_hPa']])
    xx = df1.to_numpy()
    xx = xx.ravel()
    times = minutes.to_numpy()
    times = minutes.ravel()
    return times,xx

    
def performMCSimulation(t,tau):
    global nsim
    nsim = 1000
    grxx = []
    for i in range(0,nsim):
        red = makear1(t,tau)
        psd = createLombScargleSpectrum(t,red)
        varrx = freq[1]*np.sum(psd)
        facx = varx/varrx
        psd = psd*facx
        grxx.append(psd)
    grxx = np.asarray(grxx,dtype=np.float64)
    return grxx
    
def getdof(window,noverlap):

    c50 = [0.5,0.3,0.167,0.250,0.096]
    
    c2 = 2.0 *c50[window]*c50[window]
    denom = 1.0+c2 -c2/noverlap
    neff = noverlap/denom
    dof = 2.0 *neff
    return dof
    
def createLombScargleSpectrum(t,x):
    global freq
    noverlap = 7
    ofac = 4.
    tChunked,xChunked = chunkDataIntoSegments(noverlap,t,x)
    psdSeg=[]
    freq = makeFrequencyVector1(tChunked[0,:],4,1)
    freq_angular = 2*np.pi*freq
    psd = []
    for tseg,xseg in zip(tChunked,xChunked):
        detrendedSegment = signal.detrend(xseg,type='linear')
        ww = createWelchWindow(tseg)
        detrendedWindowedX = detrendedSegment*ww
        psdseg = signal.lombscargle(tseg,detrendedWindowedX,freq_angular[1:])
        psd.append(psdseg)

    psd = np.asarray(psd,dtype=float)
    psd = np.sum(psd,axis=0)
    tpx = np.mean(np.diff(t)*tChunked.shape[1])
    dfx = 1/(4.*tpx)
    scaleFactorForAutoSpectrum = 2.0/(7*dfx*4*np.sum(ww*ww))
    psd = psd*scaleFactorForAutoSpectrum
    return psd

def makear1(t,tau):

    n = np.size(t)

    red = np.zeros(n)
    red[0] = 0.

    dt = np.diff(t)

    sigma = 1.0 - np.exp(-2.0 *dt/tau)
    sigma = np.sqrt(sigma)
    gasdev = np.random.normal(0,sigma,n-1)
    red[1:] = red[0:-1]*np.exp(-dt/tau)+sigma * gasdev
    return red
    
def gettau(t,x):
    avgdt = np.sum(np.diff(t))/(t.size-1)
    noverlap = 7
    tChunked,xChunked = chunkDataIntoSegments(noverlap,t,x)
    rhosum = 0.
    for tseg,xseg in zip(tChunked,xChunked):
        detrendedSegment = rmtrend(tseg,xseg)
        tauA,rhoA = tauest(tseg,detrendedSegment)
        rho = (rhoA*(float(tseg.size)-1.0)+1.0)/(float(tseg.size)-4.0)
        rhosum += rho

    rho = rhosum/noverlap
    tau = -avgdt/np.log(rho)
    return tau

def chunkDataIntoSegments(noverlap,t,x):
    nseg = int(2*t.size/(noverlap+1))

    tChunked = np.zeros(shape = (noverlap,nseg))
    xChunked = np.zeros(shape = (noverlap,nseg))
    istart = np.linspace(0, math.ceil(noverlap * nseg / 2), num=noverlap, endpoint=False, dtype=np.int32)
    jstart = np.linspace(0, nseg, num=nseg, endpoint=False, dtype=np.int32)

    k = istart[:, np.newaxis] + jstart[np.newaxis,:]
    tChunked = t[k]
    xChunked = x[k]
    return tChunked,xChunked


def minls(t,x):

    a = 1/math.e
    t1=t
    x1 = x
    def leastSquares(a):
        dt = np.diff(t1)
        ls = np.sum((x1[1:] - x1[:-1]*np.sign(a)*a**dt)**2)
        return ls
    try:
        amin = brent(leastSquares,brack=(-1,1),tol=3.0e-8,maxiter=10)
        #amin = minimize_scalar(leastSquares, bounds=[0, 1], method='bounded').x
    except Warning:
        print(traceback.format_exc())
    return amin
    
def tauest(t,x):

    variance = np.var(x)

    fac = np.sqrt(variance)

    xscaled = x/fac

    dt = (t[-1]-t[0])/np.size(t)

    rho = rhoest(xscaled)

    scalt = -np.log(rho)/dt

    tscaled = t*scalt
    amin = minls(tscaled,xscaled)

    tau = -1.0/(scalt*np.log(amin))
    rhoavg = np.exp(-dt/tau)
    return tau, rhoavg

def rhoest(x):

    sum1 = np.zeros(np.size(x))
    sum1 = np.zeros(np.size(x))
    sum1a = 0
    sum2a = 0

    sum1 = x[1:]*x[0:-1]
    sum2 = x[1:]**2
    cumulsum1 = np.sum(sum1)
    cumulsum2 = np.sum(sum2)

    rho=cumulsum1/cumulsum2
    return rho

def rmtrend(t,x):

    sx = 0.
    sy = 0.
    st2 = 0.
    b = 0.
    nseg = t.size
    sx = np.sum(t)
    sy = np.sum(x)
    z = 0.
    b = 0.
    detrendedSegment = np.zeros(nseg)
    sxoss = sx/nseg
    for i in range(nseg):
        z = t[i] - sxoss
        st2 += z*z
        b += z * x[i]
    bsum = b/st2

    a = (sy - sx *bsum)/nseg
    detrendedSegment = x - (a + bsum*t)
    return detrendedSegment


def createWelchWindow(t):

    nseg = float(t.size)

    fac1 = (nseg/2.0) - 0.5

    fac2 = 1.0/((nseg/2.0) + 0.5)

    fac3 = nseg - 1.0

    fac4 = (2.0 * np.pi)/(nseg -1.0)

    tlen = t[-1] - t[0]

    ww = np.zeros(t.size)

    jeff = nseg * (t[0:]-t[0])/tlen

    ww = 1.0 - ((jeff-fac1)* fac2)**2.0

    wwsq2 = (ww*ww)
    sumw2 = np.sum(wwsq2)
    scal = np.sqrt(nseg/sumw2)
    ww = ww*scal
    return ww

def makeFrequencyVector1(t,ofac=4.,hifac=1.):

    dt = np.mean(np.diff(t))
    tpx = dt * t.size
    flo =1.0 / (tpx*ofac)

    fhi = hifac / (2*dt)

    df = flo
    nf = fhi /df + 1
    nf = int(nf)

    freq = np.linspace(flo, fhi, nf)
    #flow = (2*dt)**(-1)/(len(t)*ofac)
    #fhigh = hifac/(2*dt)+(1/(2*dt))/ (len(t)*ofac)
    #nfi = ((2*dt)**(-1))/ (len(t)*ofac)
    #f = np.arange(flow,fhigh,nfi)
    return freq


def makeFrequencyVector(t,ofac=4.,hifac=1.):

    dt = np.mean(np.diff(t))

    flo = (1/(2*dt*300))/ (np.size(t)*ofac)

    fhi = hifac / (2*dt*300)
    df = flo
    nf = int((fhi - flo) / df + 1)

    freq = np.linspace(flo, fhi, nf)
    flow = (2*dt*300)**(-1)/(len(t)*ofac)
    fhigh = hifac/(2*dt*300)+(1/(2*300*dt))/ (len(t)*ofac)
    nfi = ((2*dt*300)**(-1))/ (len(t)*ofac)
    f = np.arange(flow,fhigh,nfi)
    return f

main()
