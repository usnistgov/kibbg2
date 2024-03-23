import numpy as np

def airDensity(temp,press,relHumid):     # in degC, hPa, %
    Ma= 28.96546e-3
    #Z = 0.999603
    Z = compressibility(temp,press,relHumid)
    xv =relwaterPress(temp,press,relHumid)
    R=  8.314510
    T=273.15 + temp
    p = press* 100
    ratio = 0.3780
    rho_air = p * Ma/(Z*R*T)*(1.0-xv*ratio )
    return rho_air
    
    
def buoyCorr(temp,press,relHumid, DensityOfMass):
    md = DensityOfMass
    rho_air = airDensity(temp,press,relHumid)
    buoycorr=1.0-rho_air/md
    return buoycorr

def satWaterPress(temp): # temp ind degC
    T=273.15 + temp
    A= 1.2378847e-5
    B=-1.9121316e-2
    C=33.93711047
    D=-6.3431645e3
    psv = np.exp(A*T*T+B*T+C+D/T)
    return psv
    
def waterPress(temp,press,relHumid):     # in degC, hPa, %
    p = press*100
    rH = relHumid/100
    alpha = 1.00062
    beta = 3.14e-8
    gamma=5.6e-7
    f = alpha+beta*p+gamma*temp*temp
    #f=1.0
    sp = satWaterPress(temp)
    return sp*f*rH
    

def relwaterPress(temp,press,relHumid):     # in degC, hPa, %
    wP = waterPress(temp,press,relHumid)
    p = press*100
    return wP/p


def compressibility(temp,press,relHumid):     # in degC, hPa, %
    xv =  relwaterPress(temp,press,relHumid)
    p = press*100
    T=273.15 + temp
    a0 =1.58123e-6
    a1=-2.9331e-8
    a2=1.1043e-10
    b0 = 5.707e-6
    b1=-2.051e-8
    c0=1.9898e-4
    c1=-2.376e-6
    d=1.83e-11
    e=-0.765e-8
    
    Z = 1-p/T*(a0+a1*temp+a2*temp*temp+(b0+b1*temp)*xv+(c0+c1*temp)*xv*xv)+p*p/(T*T)*(d+e*xv*xv)
    return Z
    
