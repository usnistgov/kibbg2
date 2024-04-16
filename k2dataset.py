# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 11:34:34 2024

@author: schlammi
"""
import numpy as np
import os
import k2tools
from pathlib import Path
import configparser
import datetime
from statsmodels.robust.scale import huber


class MyFiles():
    def __init__(self,c):
        self.c=c
        self.myfn={}
        self.maxGrpMem=-1
    def addFn(self,fn,grp):
        if grp not in self.myfn:
            self.myfn[grp]=[]
        self.myfn[grp].append(str(fn))
    def readGrp(self,grp,Vmul=1):
        if grp in self.myfn:
            for fi in self.myfn[grp]:
                self.readFn(grp, fi,Vmul)
                self.maxGrpMem=grp
                
    def readFn(self,grp,fi,Vmul=1):
        pass
    def totGrps(self):
        if not self.myfn:
            return 0
        return max(list(self.myfn))
    
    def clear(self):
        self.myfn={}
        self.maxGrpMem=-1
        
    
    

class MyVelos(MyFiles):
    def __init__(self,c):
        super(MyVelos, self).__init__(c)
        self.sco=0
        self.adata=[]
        self.sincdata=[]
        self.blfit=[]
        self.maxgrp=-1
        self.fit_pars=[]
        self.T0=4.34
        
    
        
    def readFn(self,grp,fi,Vmul=1000):
        data = np.loadtxt(fi)
        self.sincdata=[]
        t1 = data[:,0]
        t2 = data[:,1]
        z1 = data[:,2]
        z2 = data[:,3]
        v = data[:,5]
        V = data[:,6]*Vmul*self.c.Vcal
        S = np.ones(len(data[:,5]))*self.sco
        G = np.ones(len(data[:,5]))*grp
        self.sco+=1
        t =0.5*(t1+t2)
        z =0.5*(z1+z2)
        dset = np.c_[t,z,v,V,S,G] # 0-t, 1-z,2-v,3-V 4-S
        if len(self.adata)==0:
            self.adata  =np.array(dset)
        else:
            self.adata=np.r_[self.adata,dset]
    
    def  makesinc(self):
        hars=[1,2,3,4,5,6]
        T0= 4.240118026626671
        ix = np.where(self.adata[:,4]==0)[0]
        u=k2tools.findT0(self.adata[ix,0],self.adata[ix,2],\
                         hars,self.c.T0)
        tau =0.1
        self.sincdata = np.array(self.adata)
        smax = int(np.max(self.adata[:,4]))
        for s in range(smax+1):
            ix = np.where(self.adata[:,4]==s)[0]
            ov,av,pv,c2v,fv=k2tools.fit_sine(self.adata[ix,0],self.adata[ix,2],T0,hars)
            oV,aV,pV,c2V,fV=k2tools.fit_sine(self.adata[ix,0],self.adata[ix,3],T0,hars)
            resv=self.adata[ix,2]-fv
            resV=self.adata[ix,3]-fV
            
            ava,pva = k2tools.myatten(T0,tau,hars,av,pv)
            aVa,pVa = k2tools.myatten(T0,tau,hars,aV,pV)
            
            self.sincdata[ix,2] =k2tools.calcAP(self.adata[ix,0],T0,hars,ov,ava,pva)+resv
            self.sincdata[ix,3] =k2tools.calcAP(self.adata[ix,0],T0,hars,oV,aVa,pVa)+resV
            
            
        
    
    def fitMe(self,order=4,usesinc=False):
        if self.maxGrpMem<0:
            return
        self.z0=0
        self.order=order
        self.zmin = np.min(self.adata[:,1])
        self.zmax = np.max(self.adata[:,1])        
        self.tmin = np.min(self.adata[:,0])
        self.tmax = np.max(self.adata[:,0])    
        
        if usesinc == False:
        
            self.blfit,self.C2,self.NDF,self.fit_pars = \
                k2tools.FitLikeACanadianOrthoMaster(self.adata,\
                            zmin=self.zmin,zmax=self.zmax,\
                                order=self.order,z0=self.z0)
        else:
            if np.shape(self.sincdata)!=np.shape(self.adata):
                self.makesinc()
            self.blfit,self.C2,self.NDF,self.fit_pars = \
                k2tools.FitLikeACanadianOrthoMaster(self.sincdata,\
                            zmin=self.zmin,zmax=self.zmax,\
                                order=self.order,z0=self.z0)
                
        self.maxgrp =int(max(self.blfit[:,2]))
        self.cov=[]
        self.piecewise=[]
        for k in range(0,self.maxgrp):
            ix=np.where(np.logical_or(self.blfit[:,2]==k,self.blfit[:,2]==k+1))[0]
            if len(self.blfit[ix,0])>2:
                pfs,cov=np.polyfit(self.blfit[ix,0],self.blfit[ix,1],1,cov=True)
            elif len(self.blfit[ix,0])==2:
                pfs=np.polyfit(self.blfit[ix,0],self.blfit[ix,1],1)
                cov = np.zeros((2,2))
                cov[0,0] = 7.67111872e-15
                cov[1,1] = 2.68464606e-08
            else:
                pfs = np.zeros(2)
                pfs[1] =self.blfit[ix,1]
                cov = np.zeros((2,2))
                cov[0,0] = 7.67111872e-15
                cov[1,1] = 2.68464606e-08
                
            self.cov.append(cov)
            self.piecewise.append(pfs)
      
    def getgrp(self,t):
        o = self.blfit[0,:]
        if t<o[0]:
            return int(o[2])
        for i in self.blfit[:-1]:
            if t<i[0] and t>=o[0]:
                break
            o =i
        nr = o[2]
        if nr>=self.maxgrp-1:
            nr =self.maxgrp-1
        return int(nr)
    
    def getgrpA(self,t):
        grps=[]    
        for tt in t:
            grps.append(self.getgrp(tt))
        return np.array(grps)
        
    def getBlAndUnc(self,t):
        avals=np.array([])
        uvals=np.array([])
        for tt in t:
            grp = self.getgrp(tt)
            pf =self.piecewise[grp]
            cov=self.cov[grp]
            X =np.c_[[tt],np.ones(1)]
            unc=np.sqrt(np.diag(X@cov@X.T))
            y=X@pf
            avals=np.r_[avals,y]
            uvals=np.r_[uvals,unc]
        return avals,uvals
    
                
        
class MyForces(MyFiles):
    def __init__(self,c):
        super(MyForces, self).__init__(c)
        self.sco=0
        self.data=[]
        self.maxgrp=-1
        self.maxS=int(0)
        self.adatalen =0

    def readFn(self,grp,fi,Vmul=1):
        data = np.loadtxt(fi)
        t = data[:,0]
        I = Vmul*data[:,1]*self.c.Vcal/self.c.R*1e6
        z = data[:,2]
        S = np.ones(len(data[:,2]))*self.sco
        G = np.ones(len(data[:,2]))*grp
        self.sco+=1
        dset = np.c_[t,z,I,S,G]
        if len(self.data)==0:
            self.data  =np.array(dset)
        else:
            self.data=np.r_[self.data,dset]
        self.maxgrp=int(np.max(self.data[:,4]))
        self.maxS=int(np.max(self.data[:,3]))

    def aveForce(self):
        self.adata =[]
        for s in range(int(self.maxS)+1):
            ix = np.where(self.data[:,3]==s)[0]
            means = np.mean( self.data[ix,:],axis=0)        #order 0-time, 1-z, 2-Current, 3-S, 4-Grp
            stds  = np.std( self.data[ix,:],axis=0,ddof =1) #order 5-time, 6-z, 7-Current, 8-S, 9-Grp
            N = len(self.data[ix,:])
            if N>1:
                rtN=np.sqrt(N)
            else:
                rtN=1
            stds=stds/rtN
            if len(self.adata)==0:
                self.adata =np.array(np.r_[means,stds])
            else:
                self.adata=np.c_[self.adata,np.r_[means,stds]]
        if len(np.shape(self.adata))==2:
            self.adata =self.adata.T
        else:
            self.adata =np.expand_dims(self.adata , axis=0)
        self.adatalen =len(self.adata)
        
class Mass:
    def __init__(self,Velos,myOns,myOffs,Env,usebl=True,useg=True,\
                 usedens=True,covk=1,excl3=False):
        self.covk = covk
        self.myEnv = Env
        self.myVelos = Velos
        self.myOns = myOns
        self.myOffs =myOffs
        self.c = self.myVelos.c

        ix_on = np.where(self.myOns.adata[:,4]<self.myVelos.maxGrpMem)[0]

        
        on_t = self.myOns.adata[ix_on,0]
        on_z = self.myOns.adata[ix_on,1]
        bl,blUnc = self.myVelos.getBlAndUnc(on_t)
        bl_cor1 = k2tools.calcProfile(
            self.myVelos.fit_pars,self.myVelos.order,\
            on_z,self.myVelos.zmin,self.myVelos.zmax)
        bl_cor2 = k2tools.calcProfile(
            self.myVelos.fit_pars,self.myVelos.order,\
            self.myVelos.z0*np.ones(len(on_z)),self.myVelos.zmin,self.myVelos.zmax)
        if usebl:
            bl = bl+bl_cor1-bl_cor2
        else:
            bl=1
        if useg:
            g=self.c.g
        else:
            g=1
        on_F = self.myOns.adata[ix_on,2]*bl/g
        on_Func = self.myOns.adata[ix_on,7]*bl/g
        on_grp = self.myOns.adata[ix_on,4]

        self.on_d = np.c_[on_t,on_z,on_F,on_Func,on_grp]
        
        ix_of = np.where(self.myOffs.adata[:,4]<self.myVelos.maxGrpMem)[0]

        
        of_t = self.myOffs.adata[ix_of,0]
        of_z = self.myOffs.adata[ix_of,1]
        bl,blUnc = self.myVelos.getBlAndUnc(of_t)
        bl_cor1 = k2tools.calcProfile(
            self.myVelos.fit_pars,self.myVelos.order,\
            of_z,self.myVelos.zmin,self.myVelos.zmax)
        bl_cor2 = k2tools.calcProfile(
            self.myVelos.fit_pars,self.myVelos.order,\
            self.myVelos.z0*np.ones(len(of_z)),self.myVelos.zmin,self.myVelos.zmax)
        if usebl:
            bl = bl+bl_cor1-bl_cor2
        else:
            bl=1
        if useg:
            g=self.c.g
        else:
            g=1
        
        of_F = self.myOffs.adata[ix_of,2]*bl/g
        of_Func = self.myOffs.adata[ix_of,7]*bl/g
        of_grp = self.myOffs.adata[ix_of,4]
        self.of_d = np.c_[of_t,of_z,of_F,of_Func,of_grp]

        arow=[]
        for g in range(1+int(np.max(of_grp))):
            ofix = np.where(self.of_d[:,4]==g)[0]
            onix = np.where(self.on_d[:,4]==g)[0]
            of_d =self.of_d[ofix,:]
            on_d =self.on_d[onix,:]

            for a1,b,a2 in zip(of_d[:-1,:],on_d,of_d[1:,:]):
                ta1=a1[0]
                tb=b[0]
                ta2=a2[0]
                f = (tb-ta1)/(ta2-ta1)
                val =a1[2]*(1-f)+a2[2]*f
                unc =np.sqrt(a1[3]**2*(1-f)**2+a2[3]**2*f**2)
                grp =b[4]
                z = (1-f)*a1[1]+f*a2[1]
                row= np.r_[tb,z,val,unc,grp]
                arow.append(row)
        self.ofa_d= np.array(arow)
        diffs=[]
        for g in range(1+int(np.max(of_grp))):
            ofix = np.where(self.ofa_d[:,4]==g)[0]
            onix = np.where(self.on_d[:,4]==g)[0]
            of_d =self.ofa_d[ofix,:]
            on_d =self.on_d[onix,:]

            for a,b in zip(of_d,on_d):
                t = 0.5*(b[0]+a[0])
                airdens = np.interp(t,self.myEnv.edata[:,0],\
                                    self.myEnv.edata[:,4])
                if usedens:
                    denscorr=1+airdens/(self.c.dens)
                else:
                    denscorr=False
                di= (b[2]-a[2])*denscorr
                su =b[2]+a[2]
                z = 0.5*(a[1]+b[1])
                unc =np.sqrt(a[3]**2+b[3]**2)
                grp =b[4]
                diffs.append(np.r_[t,z,di,unc,grp])

        if len(diffs)>6 and excl3==True:
            diffs = np.array(diffs)
            mean,sig= huber(diffs[:,2])
            ix = np.where(np.abs(diffs[:,2]-mean)<5*sig)[0]
            self.dif_d=diffs[ix,:]
        else:                                
            self.dif_d= np.array(diffs)      
        self.avemass = sum(self.dif_d[:,2]/self.dif_d[:,3]**2)/sum(1/self.dif_d[:,3]**2)
        self.uncmass = 1/np.sqrt(sum(1/self.dif_d[:,3]**2))*self.covk
        self.c2=sum((self.dif_d[:,2]-self.avemass)**2/self.dif_d[:,3]**2)
        self.ndf =  len(self.dif_d[:,2])-1
        if self.ndf>1:
            self.br = np.sqrt(self.c2/self.ndf)
            if self.br>2:
                self.uncmass=self.uncmass*self.br
                self.dif_d[:,3] =self.dif_d[:,3] *self.br
            

        
        

        
class MyConfig:
    def __init__(self):
        self.clear()
        
    def clear(self):       
        self.hacconfig=False
        self.t0= datetime.datetime(2024,1,1)
        self.bd0=''
        self.mydict={}
        
    def setbd0(self,bd0):
        self.bd0=bd0
        self.config = configparser.ConfigParser()
        self.config.read(os.path.join(bd0,'config.ini'))
        verstr= self.config['Measurement']['SoftwareVersion']
        self.ver=1
        if verstr.startswith('"'):
            verstr = self.trim(verstr)
        self.mydict['Version'] = float(verstr)
        self.ver = self.mydict['Version']
        self.mydict['Title']           = self.parse(self.config['Measurement']['MeasDesc'],'string')       
        if 'SerialNo' in self.config['Measurement']:
            self.mydict['SerialNo']        = self.parse(self.config['Measurement']['SerialNo'],'string')       
        else:
            self.mydict['SerialNo'] = 'not provided'
        self.mydict['StartTime']       = self.parse(self.config['Measurement']['StartTime'],'time')
        self.mydict['Location']        = self.parse(self.config['Measurement']['Location'],'string')
        self.mydict['Nominal']        = self.parse(self.config['Measurement']['NominalMassGrams'],'float')
        self.mydict['R']               = self.parse(self.config['Calibration']['ResOhm'],'float')
        self.mydict['Runc']            = self.parse(self.config['Calibration']['ResUncPpm'],'float')*1e-6*self.mydict['R']
        self.mydict['ResCalDate']      = self.parse(self.config['Calibration']['ResCalDate'],'date')
        self.mydict['DVMREAD']         = self.parse(self.config['Calibration']['DvmCalValV'],'float')
        self.mydict['DvmCalDate']      = self.parse(self.config['Calibration']['DvmCalDate'],'date')
        self.mydict['ZENER']           = self.parse(self.config['Calibration']['ZenerRefValV'],'float')
        self.mydict['Zunc']            = self.parse(self.config['Calibration']['ZenerUncNv'],'float')*1e-9
        self.mydict['g']               = self.parse(self.config['Calibration']['LocGravity'],'float')
        self.mydict['gunc']            = self.parse(self.config['Calibration']['LocGravityUncPpm'],'float')*self.mydict['g']*1e-6
        self.mydict['Dens']            = self.parse(self.config['Calibration']['MassDensity'],'float')        
        self.mydict['VerticalityDate'] = self.parse(self.config['VelocityMode']['VerticalityDate'],'date')
        self.mydict['SinePeriodSec']   = self.parse(self.config['VelocityMode']['SinePeriodSec'],'float')

        self.title     = self.mydict['Title']       
        self.R         = self.mydict['R'] 
        self.dens      = self.mydict['Dens']
        self.Runc      = self.mydict['Runc'] 
        self.g         = self.mydict['g'] 
        self.gunc      = self.mydict['gunc'] 
        self.Vcal      = self.mydict['ZENER']/self.mydict['DVMREAD']     #multiply all V readings with that
        self.T0        = self.mydict['SinePeriodSec']
        self.hacconfig = True
    
    def trim(self,inp):
        if self.ver==1:
            return inp.replace('"', '')
        else:
            return inp

    def scaletime(self,tim):
        """
        internally in the k2 viewer time will be stored as days since 1/1/2024
        """
        return (tim-self.t0).total_seconds()/24/3600
    
    def parsetime(self,inp):
        if self.ver==1:
            start=datetime.datetime.strptime(self.trim(inp),'%m/%d/%Y %I:%M:%S %p')
        else:
            start=datetime.datetime.strptime(inp,'%m/%d/%Y %H:%M:%S')
        return self.scaletime(start)
            
    def parsedate(self,inp):
        start=datetime.datetime.strptime(self.trim(inp),'%b %d %Y')
        return self.scaletime(start)
    
    def parse(self,inp,ty):
        if ty.upper()=='STRING':
            return self.trim(inp)
        if ty.upper()=='TIME':
            return self.parsetime(inp)
        if ty.upper()=='DATE':
            return self.parsedate(inp)
        if ty.upper()=='FLOAT':
            return float(inp)
      
        
class MyEnv():
    def __init_(self):
        self.clear()
    
    def clear(self):
        self.hasEnv=0
        self.edata = []
        
    
    def setbd0(self,bd0,guesslen):
        self.bd0=bd0
        self.guesslen = guesslen

        if self.bd0=='':     
            self.makeFakeEnv()
            return            
        if 'PRTData.dat' not in os.listdir(self.bd0):
            self.makeFakeEnv()
            return
        fn =os.path.join(self.bd0,'PRTData.dat')
        da = np.loadtxt(fn,skiprows=1,usecols=[1,2,3])
        da =np.vstack((np.arange(len(da[:,0])).T*10,da.T)).T
        ix = np.where(np.logical_and(\
            np.logical_and(da[:,1]!=0,(da[:,2]!=0)),(da[:,3]!=0)))[0]
        da = da[ix,:]
        if np.shape(da)[0]==0:
            self.makeFakeEnv()
            return
        dens = k2tools.airDensity(da[:,3],da[:,2],da[:,1])
        self.edata = np.vstack((da.T,dens)).T
        if np.shape(self.edata)[0]>0:
            self.hasEnv=1
        #edata #0 = time 1 = humid 2=press 3 =temp
            
    def makeFakeEnv(self):
        rows = int((self.guesslen*1.2)/10)
        da = np.zeros((rows,4))
        for n in range(rows):
            da[n,0] = n*10
            da[n,1] = 40
            da[n,2] = 999.915
            da[n,3] = 20.0
        dens = k2tools.airDensity(da[:,3],da[:,2],da[:,1])
        self.edata = np.vstack((da.T,dens)).T
        if np.shape(self.edata)[0]>0:
            self.hasEnv=2
            
        
        
        
    
        
class k2Set():
    def __init__(self,mutex):
        """
        k2Set(mutex,ver)
        version =1.0
        """
        self.covk=1
        self.bd0=''
        self.mutex =mutex
        self.c = MyConfig()
        self.myVelos=MyVelos(self.c)
        self.myOns=MyForces(self.c)
        self.myOffs=MyForces(self.c)
        self.myEnv=MyEnv()
        self.Mass=0
        self.clear()
    
    def setcoverage(self,k):
        self.covk=2
        
    def clear(self):
        self.clearRefMass()
        self.myEnv.clear()
        self.myVelos.clear()
        self.myOns.clear()
        self.myOffs.clear()
        self.Mass=0
        
    def clearRefMass(self):
        self.mutex.lock()
        self.refMass=0
        self.hasRefMass = False
        self.mutex.unlock()
        
    def setRefMass(self,refmass):
        self.mutex.lock()
        self.refMass=refmass
        self.hasRefMass = True
        self.mutex.unlock()


    def calcMass(self,excl3=False,usebl=True,useg=True,usedens=True):
        self.Mass= Mass(self.myVelos,self.myOns,self.myOffs,self.myEnv,\
                        usebl,useg,usedens,covk=self.covk,excl3=excl3)
    
        
    def readEnv(self):
        self.myEnv.setbd0(self.bd0,self.guesslen)
        

    def setbd0(self,bd0):
        self.bd0= bd0
        self.c.setbd0(bd0)
        self.ver = self.c.ver
        self.allfiles = self.readNSort()
        self.assignFiles()
        self.totGrps  = max(self.myVelos.totGrps(),self.myOns.totGrps(),\
                      self.myOffs.totGrps())
               
    def readNSort(self):
        if self.bd0=='':
            return
        bd1 = os.path.join(self.bd0,'Force mode')
        bd2 = os.path.join(self.bd0,'Velocity mode')
        if self.ver>=1.1:
            l1 = [i for i in os.listdir(bd1)]
            l2 = [i for i in os.listdir(bd2)]
            la = sorted(l1+l2)
            u =[]
            for i in la:
                if i in l1:
                    u.append(os.path.join(bd1,i))
                elif i in l2:
                    u.append(os.path.join(bd2,i))
            tbeg=datetime.datetime.strptime(la[0].split('_')[0],'%Y%m%d%H%M%S')
            tend=datetime.datetime.strptime(la[-1].split('_')[0],'%Y%m%d%H%M%S')
            self.guesslen = (tbeg-tend).total_seconds()
            return [str(i) for i in u]
        elif self.ver>=1:
            all=[os.path.join(bd1,i) for i in os.listdir(bd1)] + [os.path.join(bd2,i) for i in os.listdir(bd2)]
            ap = [Path(i) for i in all]
            u = sorted(ap,key=os.path.getmtime)
            start = os.path.getmtime(u[0])
            stop = os.path.getmtime(u[-1])
            self.guesslen = stop-start
        return [str(i) for i in u]
    
    def getftype(self,wpath):
        """
        returns
        1 for velo
        2 for mass on
        3 for mass off
        0 else
        """
        fn = os.path.split(str(wpath))[-1]
        parts=fn.split('_')
        if parts[-1].upper()=='VMDATA.DAT':
            return 1
        if parts[-1].upper().startswith('ON') and parts[-1].upper().endswith('FMDATA.DAT'):
            return 2
        if parts[-1].upper().startswith('OFF') and parts[-1].upper().endswith('FMDATA.DAT'):
            return 3
        return 0

    
    def assignFiles(self):
        if self.bd0=='':
            return
        if len(self.allfiles)<2:
            return
        vgrp=0
        fgrp=-1
        omode=1
        self.myVelos=MyVelos(self.c)
        self.myOns=MyForces(self.c)
        self.myOffs=MyForces(self.c)
        for k in self.allfiles:
            nmode = self.getftype(k) 
            if nmode==1:
                if omode==1:
                    self.myVelos.addFn(k,vgrp)
                elif omode==2 or omode==3:
                    vgrp+=1
                    self.myVelos.addFn(k,vgrp)
            if (nmode==2 or nmode==3) and omode==1:
                fgrp+=1
            if nmode==2:
                self.myOns.addFn(k,fgrp)
            if nmode==3:
                self.myOffs.addFn(k,fgrp)
            omode=nmode
        
    def tmul(self):
        if len(self.myVelos.blfit)>1:
            maxt = self.myVelos.tmax
        else:
            maxt=1
        if maxt <300:
            mul=1 
            label='t/s'
        elif maxt<7200:
            mul=60
            label='t/min'
        elif maxt<48*3600:
            mul=3600
            label='t/hrs'
        else:            
            mul=3600*24
            label='t/days'
        return 1/mul,label
            

        
