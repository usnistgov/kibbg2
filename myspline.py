import numpy as np
import scipy.integrate as integrate

class BSpline:
    def __init__(self,nrknots,order,coeff=[],lam=0.0):
        self.nrknots = nrknots
        self.coeff   = coeff
        self.order   = order
        self.lam     = lam
        self.knots = [self.ki(i,self.order,self.nrknots) 
                      for i in range(0,self.nrknots+self.order+2)  ]
    
    def val(self,x):
        """ for a single value x """
        su=0
        if x<self.knots[0]:
            x=self.knots[0]
        if x>self.knots[-1]:
            x=self.knots[-1]
        for i in range(len(self.knots)-1-self.order):
            su+= self.coeff[i]*self.Ni(i,self.order,x)
        return su        
    
    def val2(self,x):
        """ for a single value x """
        su=0
        for i in range(len(self.knots)-1-self.order):
            su+= self.coeff[i]*self.Ni(i,self.order,x)
        return su*su        

    def calcX(self,x):
        X =[]
        for i in range(len(self.knots)-1-self.order):
             X.append([self.Ni(i,self.order,xx) for xx in x])
        return np.matrix(X).T
    
    def fit(self,x,y):
        """ 
        will generate new coefficients such that the spline
        goes best through x,yvectors
        """
        self.X = self.calcX(x)
        Y = np.matrix(y).T
        self.C = self.X.T*self.X
        CI = np.linalg.inv(self.C)
        self.coeff=CI*self.X.T*Y
        res= self.calc(x)-y
        C2r = np.dot(res,res)
        if self.lam!=0:
            Pen = integrate.quad(self.der().der().val2,0,1-1e-6,limit=5000)
            ret = C2r+Pen[0]
        else:
            ret = C2r
        return ret
 
    def calc(self,x):
        """ for a vector x
            Note only use if you can make sure that knots[0]<=x<=knots[-1]
        """
        X = self.calcX(x)
        return np.array(X*self.coeff)[:,0]
    
    def der(self):
        newsol =[]
        for i in range(len(self.knots)-1-self.order):
            if self.knots[i+self.order]!=self.knots[i]:
                newsol.append(self.order*(self.coeff[i,0]-self.coeff[i-1,0])/(self.knots[i+self.order]-self.knots[i]))
            else:
                newsol.append(0)
        ret = BSpline(self.nrknots,self.order-1,np.matrix(newsol).T)
        return ret
        
    def Ni(self,i,d,t):
        if d==0:
            if self.knots[i]<=t and t<self.knots[i+1]:
                return 1
            return 0
        else:
            kp = self.Ni(i,d-1,t)
            kn = self.Ni(i+1,d-1,t)
            if (self.knots[i+d]-self.knots[i])!=0:
                wp = (t-self.knots[i])/(self.knots[i+d]-self.knots[i]) 
            else:
                wp=0
            if (self.knots[i+d+1]-self.knots[i+1])!=0:
                wn = (self.knots[i+d+1]-t)/(self.knots[i+d+1]-self.knots[i+1])
            else:
                wn=0
            return kp*wp+\
                   kn*wn
                   
    def ki(self,i,d,n):
        if i<d:
            return 0
        else: 
            if i<=n:
                return (i-d)/(n+1-d)
            else:
                return 1
            
            
def xscale(x,r=1.001):
    """
    scales t such it runs from 0. to 1/r.
    return xs,xmin,xmax
    xs is scaled from 0 to t
    xmin and xmax are used to unscale

    """
    xmin = np.min(x)
    xmax = np.max(x)
    xs   =  (x-xmin)/(r*(xmax-xmin))
    return xs,xmin,xmax


def xfscale(x,xmin,xmax,r=1.001):
    """
    scales t such it runs from 0 to 1/r.
    but forces xmin and xmax and does not take it from the dataser
    return xs
    xs is scaled from 0 to t

    """
    xs   =  (x-xmin)/(r*(xmax-xmin))
    xs = np.clip(xs,0,1.0/r)
    return xs

    
def xunscale(xs,xmin,xmax,r=1.001):
    """
    unscales x
    """
    x   =  r*xs*(xmax-xmin)+xmin
    return x