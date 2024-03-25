# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 11:34:34 2024

@author: schlammi
"""
import numpy as np
import os
import k2tools
from pathlib import Path
import myspline
import myfilters

class k2DataSet():
    def __init__(self,mutex):
        self.bd0 = ''
        self.hasresult=False
        self.mutex =mutex
        self.clearForce()
        self.clearVelo()
        self.clearEnv()
        
        
    def setbd0(self,bd0):
        self.bd0= bd0
    def readAll(self):
        self.readForce()
        self.readEnv()

    def clearForce(self):
        self.mutex.lock()
        self.ton=[]
        self.Uon=[]
        self.tof=[]
        self.Uof=[]
        self.taon=[]
        self.Iaon=[]
        self.taof=[]
        self.Iaof=[]
        self.mutex.unlock()

    def clearEnv(self):
        self.hasEnv=False
        self.mutex.lock()
        self.edata = []
        self.mutex.unlock()

    def countForceFiles(self):
        co=0
        if 'Force mode' not in os.listdir(self.bd0):
            return 0
        bd = os.path.join(self.bd0,'Force mode')
        for files in os.listdir(bd):
            if files[7:9].upper()=='OF':
                co+=1
            elif files[7:9].upper()=='ON':
                co+=1
        return co
                
            
        
    def readForce(self):
        co=0
        self.clearForce()
        if self.bd0=='':
            return            
        if 'Force mode' not in os.listdir(self.bd0):
            return
        bd = os.path.join(self.bd0,'Force mode')
        ton=[]
        Uon=[]
        tof=[]
        Uof=[]
        Gon=[]
        Gof=[]
        taon=[]
        Uaon=[]
        taof=[]
        Uaof=[]
        cugrp=0
        lastread=-1 #0 off,1 off
        allfiles = [os.path.split(str(i))[-1] for i in 
                    sorted(Path(bd).iterdir(), key=os.path.getmtime)]
        for files in allfiles:
            if files[7:9].upper()=='OF':
                currentread =0
                if currentread==lastread:
                    cugrp+=1
                lastread=currentread
                da  =np.loadtxt(os.path.join(bd,files))
                tof.append( da[:,0])
                Uof.append( da[:,1])
                Gof.append(cugrp)
            elif files[7:9].upper()=='ON':
                currentread =1
                if currentread==lastread:
                    cugrp+=1 
                lastread=currentread
                da  =np.loadtxt(os.path.join(bd,files))
                ton.append( da[:,0])
                Uon.append( da[:,1])
                Gon.append(cugrp)
            if co==0:
                self.fohdr,self.ti = self.readForHdr(bd,files)
                if 'Resistor (ohm)' in self.fohdr:
                    self.R = self.fohdr['Resistor (ohm)']
                else:
                    self.R =10000.032685
                if 'local gravity (m/s^2)'  in self.fohdr:
                    self.g = self.fohdr['local gravity (m/s^2)']
                else:
                    self.g =9.8010
                    
                self.dens=8.03
                    
            taon=[]
            Uaon=[]
            taof=[]
            Uaof=[]                              
            for a,b in zip(ton,Uon):
                taon.append(np.mean(a))
                Uaon.append(np.mean(b))
            for a,b in zip(tof,Uof):
                taof.append(np.mean(a))
                Uaof.append(np.mean(b))                
            self.mutex.lock()
            self.Uon = Uon
            self.Uof = Uof
            self.tof = tof
            self.ton = ton
            self.Iaon = np.array(Uaon)/self.R*1e6
            self.Iaof = np.array(Uaof)/self.R*1e6
            self.taof = np.array(taof)
            self.taon = np.array(taon)
            if len(self.ton)>1:
                self.cton = np.concatenate(self.ton)
                self.cIon = np.concatenate(self.Uon)/self.R*1e6
            else:
                self.cton = np.array(self.ton).flatten()
                self.cIon = np.array(self.Uon).flatten()/self.R*1e6
            if len(self.tof)>1:   
                self.ctof = np.concatenate(self.tof)
                self.cIof = np.concatenate(self.Uof)/self.R*1e6
            else:
                self.ctof = np.array(self.tof).flatten()
                self.cIof = np.array(self.Uof).flatten()/self.R*1e6
            Fdata=np.vstack((\
                            np.hstack((self.taof,self.taon)),\
                            np.hstack((self.Iaof,self.Iaon)),\
                            np.hstack((Gof,Gon)),\
                            0*np.hstack((Gof,Gon)),\
                        )).T
            ix=Fdata[:,0].argsort()
            self.Fdata = Fdata[ix,:]
            self.mutex.unlock()
            yield co
            co=co+1
            
    def hasOn(self):
        return len(self.Uon)>0
        
    def hasOff(self):
        return len(self.Uof)>0


    def readEnv(self):
        self.clearEnv()
        if self.bd0=='':          
            return            
        if 'PRTData.dat' not in os.listdir(self.bd0):
            return
        fn =os.path.join(self.bd0,'PRTData.dat')
        da = np.loadtxt(fn,skiprows=1,usecols=[1,2,3])
        da =np.vstack((np.arange(len(da[:,0])).T*10,da.T)).T
        dens = k2tools.airDensity(da[:,3],da[:,2],da[:,1])
        self.mutex.lock()
        self.edata = np.vstack((da.T,dens)).T
        if np.shape(self.edata)[0]>0:
            self.hasEnv=True
        self.mutex.unlock()
        
    def readForHdr(self,bd,fn):
        with open(os.path.join(bd,fn)) as input_file:
            head = [next(input_file) for _ in range(6)]
        if len(head[1])>1:
            self.title=head[1][1:].strip()
        else:
            self.title='No run title set'
        fields=[o.strip() for o in head[3].split('|')]
        values=[float(f) for f  in head[4].split()[1:]]
        return dict(zip(fields,values)),head[1][1:]

    def clearVelo(self):
        self.mutex.lock()
        self.vt1=[]
        self.vt2=[]
        self.vt=[]
        self.vz1=[]
        self.vz2=[]
        self.vz=[]
        self.vv=[]
        self.vV=[]
        self.vS=[]
        self.ft = []
        self.fv = []
        self.vblt = []
        self.vblv = []

 
        self.mutex.unlock()
        
    def countVeloFiles(self):
        co=0
        if 'Velocity mode' not in os.listdir(self.bd0):
            return co
        bd = os.path.join(self.bd0,'Velocity mode')
        files=os.listdir(bd)
        for f in files:
            if f.endswith('VMData.dat'):
                co+=1
        return co
    

    def readVelo(self):
        self.clearVelo()
        if self.bd0=='':
            return            
        if 'Velocity mode' not in os.listdir(self.bd0):
            return
      
        bd = os.path.join(self.bd0,'Velocity mode')
        files=os.listdir(bd)
        Sco=0
        files = [str(i) for i in sorted(Path(bd).iterdir(), key=os.path.getmtime)]
        for f in files:
            if f.endswith('VMData.dat'):
                data = np.loadtxt(os.path.join(bd,f))
                self.mutex.lock()
                self.vt1 = np.r_[self.vt1,data[:,0]]
                self.vt2 = np.r_[self.vt2,data[:,1]]
                self.vz1 = np.r_[self.vz1,data[:,2]]
                self.vz2 = np.r_[self.vz2,data[:,3]]
                self.vv = np.r_[self.vv,data[:,5]]
                self.vV = np.r_[self.vV,data[:,6]]
                self.vS = np.r_[self.vS,np.ones(len(data[:,4]))*Sco]
                self.vt =0.5*(self.vt1+self.vt2)
                self.vz =0.5*(self.vz1+self.vz2)
                self.mutex.unlock()
                #print(f,min(data[:,0]))
                Sco+=1
                yield Sco
    
    def fitVelo(self,order=4,nrknots=5): 
        if len(self.vt1)<20:
            return
        zmin = min(self.vz)
        zmax = max(self.vz)
        vblt,vblv,fC2,fNDF = \
        k2tools.FitLikeACanadianOrthoMaster(
            self.vt,self.vz,self.vv,self.vV,self.vS,zmin,zmax,order)
        self.maxt = self.getmaxt()
        xs = myspline.xfscale(vblt,0,self.maxt)
        self.myspl = myspline.BSpline(nrknots,2)
        self.myspl.fit(xs,vblv)
        sblt_ =np.linspace(0,0.99,100)
        sblv =self.myspl.calc(sblt_)
        sblt= myspline.xunscale(sblt_,0,self.maxt)        
        self.mutex.lock()
        self.vblt = vblt
        self.vblv = vblv
        self.fC2 = fC2
        self.fNDF = fNDF        
        self.sblt = sblt
        self.sblv = sblv        
        self.mutex.unlock()
        self.calcforce()
        
    def calcforce(self): 
        ts = myspline.xfscale(self.Fdata[:,0],0,self.maxt)
        bls = self.myspl.calc(ts)
        self.Fdata[:,3]=self.Fdata[:,1]*bls # in mN
        time=[]
        vals=[]
        unc=[]
        for grp in range(int(max(self.Fdata[:,2]))):
            data0 = self.Fdata[np.where(self.Fdata[:,2]==grp)[0],:]
            val,var=myfilters.SwanFilterWithVar(data0[:,3],2)
            val =val*(2)
            var = var*4
            tm = np.mean(data0[:,0])
            airdens = np.interp(tm, self.edata[:,0] , self.edata[:,4])
            denscorr=1+airdens/(self.dens*1000)
            time.append(tm)
            vals.append(val/self.g*denscorr)
            unc.append(np.sqrt(var)/self.g*denscorr)
        Fresult=np.vstack((time,vals,unc)).T
        self.Fresult=Fresult
        self.hasresult=True
        
        
        
        
    def getmaxt(self):
        mv,mn,mf =0,0,0,
        if len(self.vt1)>0:
            mv = np.max(self.vt1)
        if len(self.ton)>0:
            mn = np.max(self.cton)
        if len(self.tof)>0:
            mf = np.max(self.ctof)
        maxt = max((mv,mn,mf))
        return maxt        
        
    def tmul(self):
        maxt = self.getmaxt()
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
            
                   
        
