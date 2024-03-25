import numpy as np
import scipy.special

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
    


def FitLikeACanadianOrthoMaster(t,xin,v,V,S,xmin=-5,xmax=5,order=4,x0=0):
     #
    #Fitting like a Canadian. 
    #
    # V = V0 + b0 *v + b1*v*x + b2*v*x**2 + b2*v*x**3

    #X  =np.matrix((np.ones(len(t)),v,v*x,v*x**2)).T

    mypol = scipy.special.eval_chebyt
    x     =  (xin-xmin)/(0.5*(xmax-xmin))-1
    x0s   =  (x0-xmin)/(0.5*(xmax-xmin))-1
    
    X = np.zeros((len(t),len(set(S))+1+order))
    X[:,0]=t-min(t)
    for i in range(1,order+1):
        #X = np.c_[X,v*scipy.special.eval_legendre(i,x)]
        X[:,i] = v*mypol(i,x)
        
    i=i+1
    oss=S[0]
    for n,s in enumerate(S):
        if s!=oss:
            i=i+1
        X[n,i] = v[n]
        
            
        oss=s
    X  =np.matrix(X)
    C =((X.T*X).I)
    fit_pars=C*X.T*np.matrix(V).T
    fit_vals =np.array( X*fit_pars)[:,0]  
    C2 = np.dot(V-fit_vals,V-fit_vals)
    NDF = len(V)
    uncert =np.sqrt( C2/NDF)
    
    vals = np.array(fit_pars)[order+1:]
    cor=0
    for i in range(1,order+1):
        cor += fit_pars[i,0]*mypol(i,x0s)
        
    ts =[]
    for s in set(S):
        ts.append(np.mean([i for i,j in zip(t,S) if j==s]))
        
    return np.array(ts),np.array(vals+cor)[:,0],C2,NDF




