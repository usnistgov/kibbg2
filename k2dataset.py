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
import configparser
import datetime



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
    
    
class MyVelos(MyFiles):
    def __init__(self,c):
        super(MyVelos, self).__init__(c)
        self.sco=0
        self.adata=[]
        self.blfit=[]
        self.maxgrp=0
        
    def readFn(self,grp,fi,Vmul=1000):
        data = np.loadtxt(fi)
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
        dset = np.c_[t,z,v,V,S,G]
        if len(self.adata)==0:
            self.adata  =np.array(dset)
        else:
            self.adata=np.r_[self.adata,dset]
            
    def fitMe(self,order=4):
        self.z0=0
        self.order=order
        self.zmin = np.min(self.adata[:,1])
        self.zmax = np.max(self.adata[:,1])        
        self.tmin = np.min(self.adata[:,0])
        self.tmax = np.max(self.adata[:,0])        
        
        self.blfit,self.C2,self.NDF,self.fit_pars = \
            k2tools.FitLikeACanadianOrthoMaster(self.adata,\
                        zmin=self.zmin,zmax=self.zmax,\
                            order=self.order,z0=self.z0)
        self.maxgrp =int(max(self.blfit[:,2]))
        self.cov=[]
        self.piecewise=[]
        for k in range(0,self.maxgrp):
            ix=np.where(np.logical_or(self.blfit[:,2]==k,self.blfit[:,2]==k+1))[0]
            pfs,cov=np.polyfit(self.blfit[ix,0],self.blfit[ix,1],1,cov=True)
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
        self.maxgrp=0
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
        for s in range(int(self.maxS)):
            ix = np.where(self.data[:,3]==s)[0]
            means = np.mean( self.data[ix,:],axis=0)        #order 0-time, 1-z, 2-Current, 3-S, 4-Grp
            stds  = np.std( self.data[ix,:],axis=0,ddof =1) #order 5-time, 6-z, 7-Current, 8-S, 9-Grp
            N = len(self.data[ix,:])
            stds=stds/np.sqrt(N)
            if len(self.adata)==0:
                self.adata  =np.array(np.r_[means,stds])
            else:
                self.adata=np.c_[self.adata,np.r_[means,stds]]
        if len(np.shape(self.adata))==2:
            self.adata =self.adata.T
        else:
            self.adata =np.expand_dims(self.adata , axis=0)
        self.adatalen =len(self.adata)
        
class Mass:
    def __init__(self,Velos,myOns,myOffs):
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
        bl = bl+bl_cor1-bl_cor2
        on_F = self.myOns.adata[ix_on,2]*bl/self.c.g
        on_Func = self.myOns.adata[ix_on,7]*bl/self.c.g
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
        bl = bl+bl_cor1-bl_cor2
        of_F = self.myOffs.adata[ix_of,2]*bl/self.c.g
        of_Func = self.myOffs.adata[ix_of,7]*bl/self.c.g
        of_grp = self.myOffs.adata[ix_of,4]
        self.of_d = np.c_[of_t,of_z,of_F,of_Func,of_grp]

        arow=[]
        for g in range(1+int(np.max(of_grp))):
            ofix = np.where(self.of_d[:,4]==g)[0]
            onix = np.where(self.on_d[:,4]==g)[0]
            #print(oifx)
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
            #print(oifx)
            of_d =self.ofa_d[ofix,:]
            on_d =self.on_d[onix,:]

            for a,b in zip(of_d,on_d):
                t = 0.5*(b[0]+a[0])
                di= b[2]-a[2]
                su =b[2]+a[2]
                z = 0.5*(a[1]+b[1])
                unc =np.sqrt(a[3]**2+b[3]**2)
                grp =b[4]
                diffs.append(np.r_[t,z,di,unc,grp])
        self.dif_d= np.array(diffs)        
        

        
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
        self.mydict['Title']      = self.parse(self.config['Measurement']['MeasDesc'],'string')       
        self.mydict['StartTime']  = self.parse(self.config['Measurement']['StartTime'],'time')
        self.mydict['Location']   = self.parse(self.config['Measurement']['Location'],'string')
        self.mydict['R']          = self.parse(self.config['Calibration']['ResOhm'],'float')
        self.mydict['Runc']       = self.parse(self.config['Calibration']['ResUncPpm'],'float')*1e-6*self.mydict['R']
        self.mydict['ResCalDate'] = self.parse(self.config['Calibration']['ResCalDate'],'date')
        self.mydict['DVMREAD']    = self.parse(self.config['Calibration']['DvmCalValV'],'float')
        self.mydict['DvmCalDate'] = self.parse(self.config['Calibration']['DvmCalDate'],'date')
        self.mydict['ZENER']      = self.parse(self.config['Calibration']['ZenerRefValV'],'float')
        self.mydict['Zunc']       = self.parse(self.config['Calibration']['ZenerUncNv'],'float')*1e-9
        self.mydict['g']          = self.parse(self.config['Calibration']['LocGravity'],'float')
        self.mydict['gunc']       = self.parse(self.config['Calibration']['LocGravityUncPpm'],'float')*self.mydict['g']*1e-6
        self.mydict['Dens']       = self.parse(self.config['Calibration']['MassDensity'],'float')        
        self.mydict['VerticalityDate'] = self.parsedate(self.config['VelocityMode']['VerticalityDate'])

        self.R = self.mydict['R'] 
        self.Runc = self.mydict['Runc'] 
        self.g = self.mydict['g'] 
        self.gunc = self.mydict['gunc'] 
        self.Vcal = self.mydict['ZENER']/self.mydict['DVMREAD']     #multiply all V readings with that
        self.hacconfig=True
    
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
            
class k2Set():
    def __init__(self,mutex):
        """
        k2Set(mutex,ver)
        version =1.0
        """
        self.bd0=''
        self.mutex =mutex
        self.c = MyConfig()
        self.myVelos=MyVelos(self.c)
        self.myOns=MyForces(self.c)
        self.myOffs=MyForces(self.c)
        self.clearEnv()
        self.Mass=0
        
    def calcMass(self):
        self.Mass= Mass(self.myVelos,self.myOns,self.myOffs)
    
    def clearEnv(self):
        self.hasEnv=False
        self.mutex.lock()
        self.edata = []
        self.mutex.unlock()
        
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
        if self.ver>=1:
            bd1 = os.path.join(self.bd0,'Force mode')
            bd2 = os.path.join(self.bd0,'Velocity mode')
            all=[os.path.join(bd1,i) for i in os.listdir(bd1)] + [os.path.join(bd2,i) for i in os.listdir(bd2)]
            ap = [Path(i) for i in all]
            u = sorted(ap,key=os.path.getmtime)
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
            

        


class k2DataSet():
    def __init__(self,mutex):
        self.bd0 = ''
        self.relUncTypeB = 3
        self.relUncTypeA = 9e99
        self.relUncTot = 9e99
        self.hasresult=False
        self.hasvelofit=False
        self.mutex =mutex
        self.clearForce()
        self.clearVelo()
        self.clearEnv()
        
    def setRelUncs(self):
        self.relDict={}
        self.relDict['Resistor']=1.4
        self.relDict['Voltmeter']=1.4
        self.relDict['mass position']=1 
        self.relDict['local acceleration']=2 
        self.relDict['verticality']=0.5 
        
        
        
    def setRefMass(self,refmass):
        self.mutex.lock()
        self.refMass=refmass
        self.hasRefMass = True
        self.mutex.unlock()
        
    def setbd0(self,bd0):
        self.bd0= bd0
    def readAll(self):
        self.readForce()
        self.readEnv()

    def clearForce(self):
        self.mutex.lock()
        self.title='no title set'
        self.relUnTypeA = 9e99
        self.relUncTot = 9e99
        self.ton=[]
        self.Uon=[]
        self.tof=[]
        self.Uof=[]
        self.taon=[]
        self.Iaon=[]
        self.taof=[]
        self.Iaof=[]
        self.hasRefMass = False
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
        self.velogrp=0
        self.hasvelofit=False
        self.hasspline=False
        self.nrknots=2
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
                data = np.loadtxt(f)
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

                Sco+=1
                self.velogrp=Sco
                yield Sco
    
    def fitVelo(self,order=4,nrknots=-1): 
        self.mutex.lock()
        self.hasvelofit=False
        self.xs = np.array([])
        self.vblt = np.array([])
        self.vblv = np.array([])
        self.fC2 = 0
        self.fNDF = 0        
        self.sblt = np.array([])
        self.sblv = np.array([])        
        self.mutex.unlock()

        if len(self.vt1)<20:
            return
        zmin = min(self.vz)
        zmax = max(self.vz)
        vblt,vblv,fC2,fNDF = \
        k2tools.FitLikeACanadianOrthoMaster(
            self.vt,self.vz,self.vv,self.vV,self.vS,zmin,zmax,order)
        self.maxt = self.getmaxt()
        xs = myspline.xfscale(vblt,0,self.maxt)
        if len(vblv)>6 and self.velogrp>6:
            if nrknots==-1:
                mknr = min(12,len(vblv))-2
                knr =range(2,mknr)       
                c2=[]
                for kn in knr:

                    myspl = myspline.BSpline(kn,3,lam=1)
                    c2.append(myspl.fit(xs,vblv))
                nrknots =knr[np.argmin(c2)]
            self.nrknots=nrknots
            self.myspl = myspline.BSpline(nrknots,3,lam=0)
            self.myspl.fit(xs,vblv)
            sblt_ =np.linspace(0,0.99,100)
            sblv =self.myspl.calc(sblt_)
            sblt= myspline.xunscale(sblt_,0,self.maxt)   
            self.hasspline=True
        else:
            self.nrknots=2
            sblt_ =np.linspace(0,0.99,100)
            pf = np.polyfit(xs , vblv,1)
            sblv = np.poly1d(pf)(sblt_)
            sblt= myspline.xunscale(sblt_,0,self.maxt) 
            
        self.mutex.lock()
        self.xs = xs
        self.vblt = vblt
        self.vblv = vblv
        self.fC2 = fC2
        self.fNDF = fNDF        
        self.sblt = sblt
        self.sblv = sblv  
        self.hasvelofit=True
        self.mutex.unlock()
        
    def calcvelo(self,ts):
        if self.hasspline:
            return self.myspl.calc(ts)
        pf = np.polyfit(self.xs , self.vblv,1)
        return np.poly1d(pf)(ts)
      
    def calcforce(self):   
        self.mass=float("nan")
        self.massunc=float("nan")
        self.hasresult=False
        if self.hasOn==False:
            return
        if self.hasOff==False:
            return
        if self.hasvelofit==False:
            return
            
        ts = myspline.xfscale(self.Fdata[:,0],0,self.maxt)
        bls = self.calcvelo(ts)
        self.Fdata[:,3]=self.Fdata[:,1]*bls*1000*(1-6e-6) # in mN
        time=[]
        vals=[]
        unc=[]
        for grp in range(1+int(max(self.Fdata[:,2]))):
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
        denom =np.sum(1/self.Fresult[:,2]**2)
        
        if denom !=0 and not np.isnan(denom):
            self.mass = np.sum(self.Fresult[:,1]/self.Fresult[:,2]**2)/denom  
            absunc = 3e-6*self.mass
            self.massunc = np.sqrt(1/denom+absunc**2)
            self.relUncTypeA = np.sqrt(1/denom)/self.mass*1e6 
            self.relUncTot = np.sqrt(self.relUncTypeA**2+self.relUncTypeB**2)
        else:
            self.mass =np.mean(self.Fresult[:,1])
            absunc = 3e-6*self.mass
            self.massunc=float("nan")
            self.relUncTypeA =9e99
            self.relUncTot = 9e99

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
            
                   
        
