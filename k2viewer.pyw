import os,sys
import traceback
import ctypes
import xlwt,xlrd,xlutils.copy
import datetime
from shutil import copyfile

##https://www.pythonguis.com/tutorials/pyqt-basic-widgets/
##
##https://www.geeksforgeeks.org/pyqt5-qtabwidget/
#https://realpython.com/python-pyqt-qthread/
#https://stackoverflow.com/questions/6783194/background-thread-with-qthread-in-pyqt
#
from PyQt5.QtCore import (
    Qt,
    QSize,
    QMutex, 
    QObject, 
    QThread, 
    pyqtSignal, 
    pyqtSlot )

from PyQt5.QtGui import QFont,QIcon
from PyQt5.QtWidgets import (
    QApplication,
    QCheckBox,
    QLabel,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QWidget,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QSpacerItem,
    QSizePolicy,
    QStatusBar,
    QSpinBox,
    QDoubleSpinBox,
    QAbstractSpinBox,
    QProgressBar
)
import mplwidget
import numpy as np
import time
import sqlite3
import k2dataset
import k2tools
try:
       import pyi_splash
       pyi_splash.update_text('UI Loaded ...')
       pyi_splash.close()
except:
    pass

mutex = QMutex()
kda =   k2dataset.k2Set(mutex)    
kda.setcoverage(2) 
kda.clear()



class Worker(QObject):
    finished = pyqtSignal()
    intReady = pyqtSignal(int,int,int)

    def __init__(self,excl3,order,usesinc):
        self.order =order
        self.excl3=excl3
        self.usesinc = usesinc
        super(QObject, self).__init__()

    @pyqtSlot()
    def procCounter(self): # A slot takes no params
        self.intReady.emit(0,0,0)
        maxgrp = kda.totGrps
        kda.readEnv()
        self.intReady.emit(2,0,0) 
        Npl = int((maxgrp+1)//20)
        if Npl==0: Npl=1
        for k in range(maxgrp+1):
            kda.myVelos.readGrp(k,Vmul=1000)
            kda.myOns.readGrp(k)
            kda.myOffs.readGrp(k)
            if k>=1 and k%Npl==0:
                kda.myVelos.fitMe(order=self.order,usesinc=self.usesinc)
                kda.myOns.aveForce()
                kda.myOffs.aveForce()
            if k>=1 and k%Npl==0:
                kda.calcMass()
            if k>Npl:
                self.intReady.emit(1,k+1,maxgrp+1)          
        kda.myVelos.fitMe(order=self.order,usesinc=self.usesinc)
        kda.myOns.aveForce()
        kda.myOffs.aveForce()
        kda.calcMass(excl3=self.excl3)
            
        self.intReady.emit(99,0,0) 
        
        self.finished.emit()
        
# Subclass QMainWindow to customize your application's main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        #if os.getcwd().startswith('Z:\\BB'):    
            #self.bd = r"K:\TableTopWattBalance\KIBB-g2\DATA"
        #else:
            
            
        self.foBig     = QFont('Arial',10)
        self.foBigBold = QFont('Arial',10)
        self.foBigBold.setBold(True)

        

        
        self.setWindowIcon(QIcon('k2viewer.png'))
        self.bd='..\DATA'
        self.idle=True
        self.statust=time.time()

        self.thread = QThread()
        self.setWindowTitle("Kibb-g2 Viewer")
        self.setFixedSize(QSize(1300, 600))
        
        self.mytable = QTableWidget(2,4)
        self.mytable.setHorizontalHeaderItem(0,\
                QTableWidgetItem("run"))
        self.mytable.setHorizontalHeaderItem(1,\
                QTableWidgetItem("title"))
        self.mytable.setHorizontalHeaderItem(2,\
                QTableWidgetItem("mass/mg"))
        self.mytable.setHorizontalHeaderItem(3,\
                QTableWidgetItem("unc/mg"))

        self.mytable.setColumnWidth(0, 100)
        self.mytable.setColumnWidth(1, 200)
        self.mytable.setColumnWidth(2, 120)
        self.mytable.setColumnWidth(3, 70)
        

        self.mytable.setFixedWidth(400)
        self.Brefresh = QPushButton("Reload")
        
        ### The plot windows
        
        self.mplfor      = mplwidget.MplWidget2()
        self.mplenv      = mplwidget.MplWidget4()
        self.mplvel      = mplwidget.MplWidget(rightax=True)
        self.mplmass     = mplwidget.MplWidget(rightax=True)
        self.mplprofile  = mplwidget.MplWidget(rightax=True)
        
        ### check and spin boxes
        
        self.cbShowVolt = QCheckBox()
        self.cbUseSync  = QCheckBox()
        self.cbMvsZ     = QCheckBox()
        self.cbExc3sig  = QCheckBox()
        self.sbOrder    = QSpinBox()
        self.sbMass     = QDoubleSpinBox()
        self.sbOrder.setValue(6)
        self.sbOrder.setMinimum(1)
        self.sbOrder.setMaximum(10)
        self.cbUseSync.setChecked(True)
        self.cbExc3sig.setChecked(True)
      
        self.sbMass.setMinimumWidth(100)
        self.sbMass.setMinimum(0)
        self.sbMass.setMaximum(99999)
        self.sbMass.setDecimals(4)
        self.sbMass.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.sbMass.setKeyboardTracking(False)
        
        ### global labels
        
        self.sblabel  = QLabel("Click on a run") 
        self.lares    = QLabel("") 
        self.laUncTot = QLabel("n/a")
        self.laMass   = QLabel("")
        self.laUnc    = QLabel("")
        self.laUncB    = QLabel("")
        self.laMass2   = QLabel("")
        self.laTotUnc  = QLabel("")
        self.lacov     = QLabel("(k={0})".format(kda.covk))
        self.laNoEnv   = QLabel("")

        
        self.laUa= []
        self.laUaMaxRows =8
        self.laUaMaxCols =3
        for i in range(self.laUaMaxRows):
            row=[]
            for j in range(self.laUaMaxCols):
                row.append( QLabel(""))
            self.laUa.append(row)
        for i in range(self.laUaMaxRows):
            for j in range(self.laUaMaxCols):
                self.laUa[i][j].setText('')

                
        self.Uncdict ={'Resistance':0.21,'Voltage': 1.0,'Mass position': 1.0,
                       'g': 2.0, 'Verticality': 0.5,
                       'Type A': -1, 'Total': -1, 'Balance mechanics': -1}
        
        
        self.resultLabels=[
            'Serial Number',
            'Weight designation',
            'True mass',
            'Assumed density',
            'Conventional mass',
            'Deviation from nominal',
            'Total uncertainy',
            'Tolerance (Class 3)',
            'Temperatue',
            'Barometric pressure',
            'Humidity'\
            ]
        self.laResult=[]

        for i in self.resultLabels:
            row=[]
            for j in range(3):
                if j==0:
                    la =QLabel(i)
                else:
                    la=QLabel('')
                    la.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
                la.setFont(self.foBig)
                row.append(la )
            self.laResult.append(row)
                
        
        self.laMass.setFont(self.foBigBold)
        self.laUnc.setFont(self.foBigBold)
        self.laUncB.setFont(self.foBigBold)
        
        ### Status bar

        self.progressBar = QProgressBar()
        self.progressBar.setMaximumWidth(400)
        self.tabWidget = MyTabWidget(self) 
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        widget = QWidget(self)
        widget.setLayout(QHBoxLayout())
        widget.layout().addWidget(self.sblabel)
        widget.layout().addWidget(self.progressBar)
   
        #self.statusBar.setStyleSheet("border 
        
        self.statusBar.addPermanentWidget(widget) 
        
        
        
        self.statusBar.showMessage('Welcome to k2viewer',5000)

        layout   = QVBoxLayout()
        hlayout  = QHBoxLayout()
        hlayout2  = QHBoxLayout()
        vlayout1 = QVBoxLayout()
        vlayout2 = QVBoxLayout()
        
    
        
        widget = QWidget()
        widget.setLayout(layout)
        layout.addLayout(hlayout)
        hlayout.addLayout(vlayout1)
        hlayout.addLayout(vlayout2)
        
        vlayout2.addLayout(hlayout2)
        vlayout2.addWidget(self.tabWidget)

        vlayout1.addWidget(self.mytable)
        vlayout1.addWidget(self.Brefresh)
        
        l2 = QLabel() 
        l2.setText("order")
        hlayout2.addWidget(l2)
        hlayout2.addWidget(self.sbOrder)
        hlayout2.addWidget(self.cbUseSync)
        hlayout2.addWidget(QLabel('use sinc'))
        
        hSpacer = QSpacerItem(20, 2,QSizePolicy.Expanding,
                                     QSizePolicy.Minimum)
        
        hlayout2.addItem(hSpacer)
     
        self.setCentralWidget(widget)
        self.loadTable()

        self.mytable.clicked.connect(self.on_table_clicked)
        self.Brefresh.clicked.connect(self.loadTable)
        self.cbShowVolt.clicked.connect(self.plotForce)
        self.cbUseSync.clicked.connect(self.recalcvelo)
        self.cbMvsZ.clicked.connect(self.plotMass)
        self.cbExc3sig.clicked.connect(self.plotMass)
        self.tabWidget.tabs.currentChanged.connect(self.replot)
        self.sbOrder.valueChanged.connect(self.recalcvelo)
        self.sbMass.valueChanged.connect(self.gotmassval)
        

    def createdb(self):
        connection = sqlite3.connect('k2viewer.db')
        cursor = connection.cursor()
        cursor.execute("""
                       CREATE TABLE k2data (run TEXT UNIQUE, value FLOAT,
                    uncertainty FLOAT, title TEXT, time TIMESTAMP, airdens FLOAT)""")
                       
        connection.close()



    def loadTable(self):
         c = k2dataset.MyConfig()
         self.mytable.clearContents()
         self.mytable.setRowCount(0)
         if os.path.isfile('k2viewer.db')==False:
             self.createdb()
             
         connection = sqlite3.connect('k2viewer.db')
         cursor = connection.cursor()

  
         for s1 in sorted([ f.path for f in os.scandir(self.bd) \
                           if f.is_dir() ],reverse=True):
            yymm=os.path.split(s1)[-1]
            if len(yymm)==4:
                for s2 in sorted([ f.path for f in os.scandir(s1)\
                                  if f.is_dir() ],reverse=True):
                    day=os.path.split(s2)[-1]
                    if len(day)==2:
                        for s3 in sorted([ f.path for f in os.scandir(s2)\
                                          if f.is_dir() ],reverse=True):
                            if os.path.isfile(os.path.join(s3,'config.ini'))==False:
                                continue
                            c.setbd0(s3)
                            letter=os.path.split(s3)[-1]
                            if len(letter)==1:
                                run =yymm+day+letter
                                row_number = self.mytable.rowCount()
                                self.mytable.insertRow(row_number)
                                self.mytable.setItem(row_number,0,\
                                 QTableWidgetItem(str(run)))
                                self.mytable.setItem(row_number,1,QTableWidgetItem(str(c.title)))
                                dbentry = cursor.execute("""SELECT run,value,
                                uncertainty,title FROM k2data 
                                WHERE run="{0}";""".format(run)).fetchall()
                                if len(dbentry)==0:
                                    continue
                                dbentry=dbentry[0]
                                if dbentry[1]>-9e96:
                                    nitem =QTableWidgetItem('{0:,.4f}'.format(dbentry[1]) )
                                else:
                                    nitem =QTableWidgetItem('n/a')
                                nitem.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
                                self.mytable.setItem(row_number,2,nitem)
                                if dbentry[2]>-9e96:
                                    nitem =QTableWidgetItem('{0:6.4f}'.format(dbentry[2]))
                                else:
                                    nitem =QTableWidgetItem('n/a')                                
                                nitem.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
                                self.mytable.setItem(row_number,3,nitem)
                                #self.mytable.setItem(row_number,4,QTableWidgetItem(dbentry[3]))

         connection.close()
         
    def plotForce(self):
        if kda.myOffs.maxGrpMem<0:
            return
        self.mplfor.canvas.ax1.clear()
        self.mplfor.canvas.ax2.clear()
        if self.cbShowVolt.isChecked():
            scale = kda.myOffs.c.R/1000
            yla1='U(on)/mV'
            yla2='U(off)/mV'          
        else:
            scale =1
            yla1='I(on)/uA'
            yla2='I(off)/uA'
        mutex.lock()
        tmul,tla = kda.tmul()
        # if kda.myOns.maxS>=1:
        #     p1,=self.mplfor.canvas.ax1.plot(
        #         kda.myOns.data[:,0]*tmul,\
        #         kda.myOns.data[:,2]-mon,'r.')
        # if kda.myOffs.maxS>=1:
        #     p2,=self.mplfor.canvas.bx1.plot(\
        #         kda.myOffs.data[:,0]*tmul,\
        #         kda.myOffs.data[:,2]-mof,'b.')
        if kda.myOns.adatalen>0:
            _=self.mplfor.canvas.ax1.errorbar(
                kda.myOns.adata[:,0]*tmul,\
                kda.myOns.adata[:,2]*scale,\
                kda.myOns.adata[:,7]*scale,\
                    fmt='ro')           
        if kda.myOffs.adatalen>0:
            _=self.mplfor.canvas.ax2.errorbar(
                kda.myOffs.adata[:,0]*tmul,\
                kda.myOffs.adata[:,2]*scale,\
                kda.myOffs.adata[:,7]*scale,\
                    fmt='bs')
            
        self.mplfor.canvas.setsamexscale()
        mutex.unlock()
        self.mplfor.canvas.ax1.set_ylabel(yla1)
        self.mplfor.canvas.ax2.set_ylabel(yla2)
        self.mplfor.canvas.ax2.set_xlabel(tla)
        self.mplfor.canvas.draw() 
            
    def plotEnv(self):
        self.mplenv.canvas.ax1.clear()
        self.mplenv.canvas.ax2.clear()
        self.mplenv.canvas.ax3.clear()
        self.mplenv.canvas.ax4.clear()
        if kda.myEnv.hasEnv==False:
            return
        mutex.lock()
        tmul,tla = kda.tmul()
        self.mplenv.canvas.ax1.plot(kda.myEnv.edata [:,0]*tmul,\
                                    kda.myEnv.edata[:,1],'r.')
        self.mplenv.canvas.ax2.plot(kda.myEnv.edata [:,0]*tmul,\
                                    kda.myEnv.edata[:,2],'b.')
        self.mplenv.canvas.ax3.plot(kda.myEnv.edata [:,0]*tmul,\
                                    kda.myEnv.edata[:,3],'g.')
        self.mplenv.canvas.ax4.plot(kda.myEnv.edata [:,0]*tmul,\
                                    kda.myEnv.edata[:,4],'m.')
        mutex.unlock()
        self.mplenv.canvas.ax1.ticklabel_format(useOffset=False)
        self.mplenv.canvas.ax2.ticklabel_format(useOffset=False)
        self.mplenv.canvas.ax3.ticklabel_format(useOffset=False)
        self.mplenv.canvas.ax4.ticklabel_format(useOffset=False)
        self.mplenv.canvas.ax1.xaxis.set_ticklabels([])
        self.mplenv.canvas.ax3.xaxis.set_ticklabels([])
        self.mplenv.canvas.ax1.tick_params(axis='x',direction='inout')
        self.mplenv.canvas.ax2.tick_params(axis='x',direction='inout')
        self.mplenv.canvas.ax3.tick_params(axis='x',direction='inout')
        self.mplenv.canvas.ax4.tick_params(axis='x',direction='inout')
        self.mplenv.canvas.ax1.set_ylabel('rel. humid (%)')
        self.mplenv.canvas.ax2.set_ylabel('press. (hPa)')
        self.mplenv.canvas.ax3.set_ylabel('temp (degC)')
        self.mplenv.canvas.ax4.set_ylabel('air dens. (kg/m^3)')
        self.mplenv.canvas.ax2.set_xlabel(tla)  
        self.mplenv.canvas.ax4.set_xlabel(tla)  
        self.mplenv.canvas.draw() 


    def plotVelocity(self):
        self.mplvel.canvas.ax1.clear()
        mutex.lock()
        tmul,tla = kda.tmul()
        if len(kda.myVelos.blfit)>0:
            p1,=self.mplvel.canvas.ax1.plot(
                kda.myVelos.blfit[:,0]*tmul,\
                kda.myVelos.blfit[:,1],'b.')
            if kda.myVelos.maxgrp>=1:
                tt = np.linspace(kda.myVelos.tmin,kda.myVelos.tmax,400)
                val,unc =  kda.myVelos.getBlAndUnc(tt) 
                p1,=self.mplvel.canvas.ax1.plot(
                    tt*tmul,\
                   val,'k-')
                self.mplvel.canvas.ax1.fill_between(
                    tt*tmul,val-unc,val+unc,fc='r',alpha=0.1)
            
        
        #tmul,tla = kda.tmul()
        
#        if len(kda.vblv)>=2:
 #           p1,=self.mplvel.canvas.ax1.plot(
  #              kda.vblt*tmul,\
   #             np.abs(kda.vblv)*1e3,'b.')
    #        p2,=self.mplvel.canvas.ax1.plot(
     #           kda.sblt*tmul,\
      #          np.abs(kda.sblv)*1e3,'r-')
        mutex.unlock()
        self.mplvel.canvas.bx1.ticklabel_format(useOffset=False)
        self.mplvel.canvas.ax1.ticklabel_format(useOffset=False)        
        be,en =self.mplvel.canvas.ax1.get_ylim()
        self.mplvel.canvas.draw() 
        me =0.5*(be+en)
        self.mplvel.canvas.bx1.set_ylim((be/me-1)*1e6,(en/me-1)*1e6)
        self.mplvel.canvas.ax1.set_xlabel(tla)    
        self.mplvel.canvas.ax1.set_ylabel('Bl/Tm')
        self.mplvel.canvas.bx1.set_ylabel('rel. change /ppm')
        self.mplvel.canvas.draw() 

    def plotMass(self):
        if kda.Mass==0:
            return
        kda.calcMass(excl3=self.cbExc3sig.isChecked())
        self.populateUnc()
        self.mplmass.canvas.ax1.clear()
        mutex.lock()
        tmul,tla = kda.tmul()
        if self.cbMvsZ.isChecked()==False:
            self.mplmass.canvas.ax1.errorbar(
                kda.Mass.dif_d[:,0]*tmul,kda.Mass.dif_d[:,2],\
                    kda.Mass.dif_d[:,3],fmt='bo')
        else:
            self.mplmass.canvas.ax1.errorbar(
                kda.Mass.dif_d[:,1],kda.Mass.dif_d[:,2],\
                    kda.Mass.dif_d[:,3],fmt='bo')
            tla='z/mm'
        mutex.unlock()
        be,en =self.mplmass.canvas.ax1.get_xlim()
        mean = kda.Mass.avemass
        sig  =  kda.Mass.uncmass
        sigB = self.totuncB/1000
        sigAll = self.totuncabs/1000
        self.laMass.setText('{0:8.4f} mg'.format(mean))
        self.laUnc.setText('\u00B1 {0:4.1f} \u00B5g (Type A)'.format(sig*1000))
        self.laUncB.setText('\u00B1 {0:4.1f} \u00B5g (Type B)'.format(sigB*1000))
        if kda.myEnv.hasEnv==2:
            self.laNoEnv.setText('No Env data available')
        else:
            self.laNoEnv.setText('')
            
        self.mplmass.canvas.ax1.plot((be,en),\
                        (mean,mean),c='k',linestyle='dashed',lw=2)
        self.mplmass.canvas.ax1.plot((be,en),(mean+sig,mean+sig),\
                                      c='r',linestyle='dashdot')
        self.mplmass.canvas.ax1.plot((be,en),(mean-sig,mean-sig),\
                                      c='r',linestyle='dashdot')
      
        self.mplmass.canvas.ax1.plot((be,en),(mean-sigAll,mean-sigAll),\
                                      c='b',linestyle='dashdot')
      
        self.mplmass.canvas.ax1.plot((be,en),(mean+sigAll,mean+sigAll),\
                                      c='b',linestyle='dashdot')
        self.mplmass.canvas.ax1.fill_between((be,en), (mean-sigAll,mean-sigAll),\
                                 (mean+sigAll,mean+sigAll),color='b', alpha=0.2)
      
        self.mplmass.canvas.ax1.fill_between((be,en), (mean-sig,mean-sig),\
                                 (mean+sig,mean+sig),color='r', alpha=0.2)
        self.mplmass.canvas.ax1.set_xlim(be,en)
        if kda.hasRefMass==True:
             be,en =self.mplmass.canvas.ax1.get_xlim()
             self.mplmass.canvas.ax1.plot((be,en),
                                          (kda.refMass,kda.refMass),c='m',
                                          linestyle='dotted',lw=4)
             self.mplmass.canvas.ax1.set_xlim(be,en)
            
        be,en =self.mplmass.canvas.ax1.get_ylim()
            
        #me =0.5*(be+en)
        me = mean
        if kda.hasRefMass:
            me =kda.refMass
        self.mplmass.canvas.bx1.set_ylim((be-me)*1e3,(en-me)*1e3)
        self.mplmass.canvas.ax1.set_xlabel(tla)    
        self.mplmass.canvas.ax1.set_ylabel('mass /mg')
        self.mplmass.canvas.bx1.set_ylabel('deviaton  /\u00B5g')
        self.mplmass.canvas.bx1.ticklabel_format(useOffset=False)
        self.mplmass.canvas.ax1.ticklabel_format(useOffset=False)
        self.mplmass.canvas.draw() 
        
        
    def plotProfile(self):
        if kda.myVelos.maxGrpMem<0:
            return
        if len(kda.myVelos.fit_pars)<2:
            return
        self.mplprofile.canvas.ax1.clear()
        z=np.linspace(-1.5,1.5,200)
        bl=k2tools.calcProfile(kda.myVelos.fit_pars,kda.myVelos.order,\
            z,kda.myVelos.zmin,kda.myVelos.zmax,withOffset=True)
        self.mplprofile.canvas.ax1.plot(z,bl,'r-')
        be,en =self.mplprofile.canvas.ax1.get_ylim()
        me =0.5*(be+en)
        self.mplprofile.canvas.bx1.set_ylim((be/me-1)*1e6,(en/me-1)*1e6)
        self.mplprofile.canvas.ax1.set_ylabel('Bl (Tm)')
        self.mplprofile.canvas.bx1.set_ylabel('rel. change /ppm')
        self.mplprofile.canvas.ax1.set_xlabel('z (mm)')
        self.mplprofile.canvas.draw() 
        
    def populateUnc(self):
        mean = kda.Mass.avemass
        sig  =  kda.Mass.uncmass
        self.Uncdict['Balance mechanics']=1e-3/mean*1e6
        self.Uncdict['Type A']=sig/mean*1e6/kda.covk
        cumAll=0
        cumB=0
        for k,v in self.Uncdict.items():
            if k=='Total':
                continue
            if k!='Type A':
                cumB += v*v
            cumAll+=v*v
        self.Uncdict['Total'] = np.sqrt(cumAll)
        self.totuncrel = np.sqrt(cumAll)
        self.totuncB = np.sqrt(cumB)*1e-6*mean*1000*kda.covk
        self.totuncabs = np.sqrt(cumAll)*1e-6*mean*1000*kda.covk
        
        
    def plotUnc(self):
        if kda.Mass==0:
            return
        self.populateUnc()
        row=0
        mean = kda.Mass.avemass
        for cat,rel in sorted(self.Uncdict.items(),\
                              key=lambda item: item[1],reverse=True):
            if cat=='Total':
                prow =self.laUaMaxRows-1
            else:
                prow=row
                row=row+1
            self.laUa[prow][0].setText('{0:14}'.format(cat))
            self.laUa[prow][1].setText('{0:4.1f} ppm'.format(rel*kda.covk))
            self.laUa[prow][2].setText('{0:6.1f} \u00B5g'.format(rel*1e-3*mean*kda.covk))
            if cat=='Total':
                myFont=self.foBigBold
            else:
                myFont=self.foBig
            self.laUa[prow][0].setFont(myFont)
            self.laUa[prow][1].setFont(myFont)
            self.laUa[prow][2].setFont(myFont)
            
        self.laMass2.setText('{0:8.4f} mg'.format(mean))
        self.laTotUnc.setText('\u00B1  {0:6.1f} \u00B5g (k={1})'\
                              .format( self.totuncabs,kda.covk))

            
        


    def gotmassval(self,val):
        if val!=0:
            kda.setRefMass(val)
            self.plotMass()
        


    def recalcvelo(self):
        order = int(self.sbOrder.value() )
        if kda.myVelos.maxGrpMem>0:
            kda.myVelos.fitMe(order,usesinc=self.cbUseSync.isChecked())
            kda.calcMass(excl3=self.cbExc3sig.isChecked())
            self.replot()

    def convmass(self,truemass,density):
        aird = 1.2
        steeld =8000
        conv = truemass*(1-aird/density+aird/steeld)
        return conv
    
    def getol(self,nom):
        """ reports toleance in mg, nom is in g"""
        if nom==20:
            return 0.35
        elif nom==10:
            return 0.25
        elif nom==5:
            return 0.18
        elif nom==3:
            return 0.15        
        elif nom==2:
            return 0.13
        elif nom==1:
            return 0.1
        elif nom==0.5:
            return 0.08
        elif nom==0.3:
            return 0.07
        elif nom==0.2:
            return 0.06
        elif nom==0.1:
            return 0.05
        elif nom==0.05:
            return 0.042
        elif nom==0.02:
            return 0.035
        elif nom==0.01:
            return 0.030
        return 0
        
    
    def WriteExcelHdr(self,fifn):
        if os.path.isfile('templ.xls'):
            copyfile('templ.xls',fifn)
            return
        book = xlwt.Workbook(encoding="utf-8")

        sheet1 = book.add_sheet("Sheet 1")

        sheet1.write(0,  0, "Serial Number")
        sheet1.write(0,  1, "Designated Weight")
        sheet1.write(0,  2, "True Mass")
        sheet1.write(0,  3, "Assumed Dens.")
        sheet1.write(0,  4, "Conv. Mass")
        sheet1.write(0,  5, "Dev. from Nom.")
        sheet1.write(0,  6, "Total Unc.")
        sheet1.write(0,  7, "Tolerance")
        sheet1.write(0,  8, "Temp.")
        sheet1.write(0 , 9, "Press.")
        sheet1.write(0, 10, "Humid.")

        sheet1.write(1,  0, "")
        sheet1.write(1,  1, "")
        sheet1.write(1,  2, "g")
        sheet1.write(1,  3, "g/cm3")
        sheet1.write(1,  4, "g")
        sheet1.write(1,  5, "mg")
        sheet1.write(1,  6, "mg")
        sheet1.write(1,  7, "mg")
        sheet1.write(1,  8, "degC")
        sheet1.write(1,  9, "mmHg")
        sheet1.write(1, 10, "%rel")
        
        book.save(fifn)
            
    def WriteExcel(self):
        if kda.Mass==0:
            return
        now = datetime.datetime.now()
        fn = 'kibbg2_{0}.xls'.format(now.strftime('%Y%m'))
        fifn=os.path.join(r'..\results',fn)
        if os.path.isfile(fifn)==False:
            self.WriteExcelHdr(fifn)
        rb = xlrd.open_workbook(fifn,formatting_info=True)
        r_sheet = rb.sheet_by_index(0) 
        r = r_sheet.nrows
        wb = xlutils.copy.copy(rb) 
        sheet = wb.get_sheet(0) 
        ser = kda.c.mydict['SerialNo']
        sheet.write(r,0,now.strftime('%m/%d/%Y'))
        sheet.write(r,1,now.strftime('%H:%M:%S'))
        sheet.write(r,2, ser)
        nom =kda.c.mydict['Nominal']
        if nom>=0.995:
            sheet.write(r,3, '{0}'.format(nom))
        else:
            sheet.write(r,3, '{0}'.format(nom))            
        m =kda.Mass.avemass
        if m>=995:
            sheet.write(r,4, '{0:10.7f}'.format(m/1000))
        else:
            sheet.write(r,4, '{0:10.9f}'.format(m/1000))
            
        conv = self.convmass(m,kda.c.dens)
        sheet.write(r,5, '{0:10.7f}'.format(kda.c.dens/1000))
        if conv>=995:
            sheet.write(r,6, '{0:10.7f}'.format(m/1000))
        else:
            sheet.write(r,6, '{0:10.9f}'.format(m/1000))
        sheet.write(r,7, '{0:8.4f}'.format(conv-nom*1000))
        unc =self.totuncabs
        sheet.write(r,8, '{0:6.4f}'.format(unc/1000))
        tol = self.getol(nom)            
        sheet.write(r,9, '{0:6.4f}'.format(tol/1000))
        if kda.myEnv.hasEnv==2:
            sheet.write(r,10, 'nominal')
            sheet.write(r,11, 'nominal')
            sheet.write(r,12,'nominal')
        else:
            temp = np.mean(kda.myEnv.edata[:,3])
            sheet.write(r,10, '{0:6.3f}'.format(temp))          
            press = np.mean(kda.myEnv.edata[:,2])/1.33322
            sheet.write(r,11, '{0:6.3f}'.format(press))
            hum = np.mean(kda.myEnv.edata[:,1])
            sheet.write(r,12, '{0:6.2f}'.format(hum))
        wb.save(fifn)
    
    def plotReport(self):    
        if kda.Mass==0:
            for i in self.laResult:
                i[1].setText('')
                i[2].setText('')
            return
        else:
            self.populateUnc()
            ser = kda.c.mydict['SerialNo']
            self.laResult[0][1].setText(ser)
            nom =kda.c.mydict['Nominal']
            if nom>=0.995:
                self.laResult[1][1].setText('{0}'.format(nom))
                self.laResult[1][2].setText('g')
            else:                
                self.laResult[1][1].setText('{0}'.format(nom*1000))
                self.laResult[1][2].setText('mg')
            m =kda.Mass.avemass
            if m>=995:
                self.laResult[2][1].setText('{0:10.7f}'.format(m/1000))
                self.laResult[2][2].setText('g')
            else:
                self.laResult[2][1].setText('{0:10.4f}'.format(m))
                self.laResult[2][2].setText('mg')
            conv = self.convmass(m,kda.c.dens)
            self.laResult[3][1].setText('{0:4.1f}'.format(kda.c.dens/1000))
            self.laResult[3][2].setText('g/cm\u00B3')
            if conv>=995:
                self.laResult[4][1].setText('{0:10.7f}'.format(conv/1000))
                self.laResult[4][2].setText('g')
            else:
                self.laResult[4][1].setText('{0:10.4f}'.format(conv))
                self.laResult[4][2].setText('mg')
            
            self.laResult[5][1].setText('{0:8.4f}'.format(conv-nom*1000))
            self.laResult[5][2].setText('mg')
               
            unc =self.totuncabs
            self.laResult[6][1].setText('{0:6.4f}'.format(unc/1000))
            self.laResult[6][2].setText('mg')
            tol = self.getol(nom)            
            self.laResult[7][1].setText('{0:6.4f}'.format(tol))
            self.laResult[7][2].setText('mg')

            if kda.myEnv.hasEnv==2:
                suf=' (n/a)'
            else:
                suf=''
            
            temp = np.mean(kda.myEnv.edata[:,3])
            self.laResult[8][1].setText('{0:6.3f}'.format(temp))
            self.laResult[8][2].setText('\u00b0C'+suf)
                
            press = np.mean(kda.myEnv.edata[:,2])/1.33322
            self.laResult[9][1].setText('{0:6.3f}'.format(press))
            self.laResult[9][2].setText('mm Hg'+suf)
            hum = np.mean(kda.myEnv.edata[:,1])
            self.laResult[10][1].setText('{0:6.2f}'.format(hum))
            self.laResult[10][2].setText('% rel'+suf)
            
                
        
        


    def replot(self):
        currentIndex=self.tabWidget.tabs.currentIndex()
        tat = self.tabWidget.tabs.tabText(currentIndex)
        if tat=='Force':
           self.plotForce()
        if tat=='Environmentals':
            self.plotEnv()
        elif tat=='Velocity':
            self.plotVelocity()
        elif tat=='Mass':
            self.plotMass()
        elif tat=='Profile':
            self.plotProfile()
        elif tat=='Uncertainty':
            self.plotUnc()    
        elif tat=='Report':
            self.plotReport()
   
            
    def updateTable(self):
        if kda.Mass==0:
            return
        self.populateUnc()
        mass = kda.Mass.avemass
        massunc = self.totuncabs
        title =  kda.c.title
        title = title.replace('"', '\'')
        nitem =QTableWidgetItem('{0:,.4f}'.format(mass) )
        nitem.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
        self.mytable.setItem(self.calcrow,2,nitem)
        nitem =QTableWidgetItem('{0:6.4f}'.format(massunc))
        nitem.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
        self.mytable.setItem(self.calcrow,3,nitem)
        self.mytable.setItem(self.calcrow,1,QTableWidgetItem(title))

        mycmd ="""
        replace into k2data (run,value, uncertainty,title)
        values  ("{0}",{1},{2},"{3}");""".\
           format(self.runid,mass,massunc,title)
        connection = sqlite3.connect('k2viewer.db')
        cursor = connection.cursor()
        #print(mycmd)
        cursor.execute(mycmd)
        connection.commit()
        connection.close()    
   
    
   
    def readStatus(self,ix,cur,tot):
        if ix==0: 
            self.start=time.time()
            return
        elif ix==1:
            self.progressBar.setValue(int(100*cur/tot))
            self.sblabel.setText('reading {0} sets'.format(tot))
        elif ix==2:
            self.progressBar.setValue(0)
            self.sblabel.setText('reading Environmentals')
        elif ix==99:
            self.sblabel.setText('reading of {0} done'.format(kda.bd0))
            self.statusBar.showMessage('Data available',5000)
            self.progressBar.setValue(0)
            self.populateUnc()
            self.updateTable()
            try:
                self.WriteExcel()
            except:
                pass

        self.replot()
        self.idle=True
#           
            
        

    def on_table_clicked(self,item):
        if self.idle:
            self.idle=False
            self.sbMass.clear()
            self.calcrow = item.row()
            self.runid =self.mytable.item(item.row(), 0).text()
            filePath = os.path.join(self.bd,self.runid[0:4],self.runid[4:6],self.runid[6])
            
            kda.clear()
            kda.setbd0(filePath)
            order = int(self.sbOrder.value() )
            usesinc=self.cbUseSync.isChecked()
            excl3=self.cbExc3sig.isChecked()
    
            self.obj = Worker(excl3,order,usesinc)  # no parent!
            self.thread = QThread()  # no parent!
            self.obj.intReady.connect(self.readStatus)
            self.obj.moveToThread(self.thread)
            self.obj.finished.connect(self.thread.quit)
            self.thread.started.connect(self.obj.procCounter)
            self.thread.start()
        else:
            self.statusBar.showMessage('Wait till last run is processed',5000)
         
class MyTabWidget(QWidget): 
    def __init__(self, parent): 
        super(QWidget, self).__init__(parent)  
        self.layout = QVBoxLayout(self) 
        # Initialize tab screen 
        self.tabs     = QTabWidget() 
        self.tabForce = QWidget() 
        self.tabEnv   = QWidget() 
        self.tabVelo  = QWidget() 
        self.tabMass  = QWidget() 
        self.tabUnc     = QWidget() 
        self.tabProfile = QWidget() 
        self.tabReport = QWidget() 
        self.tabs.resize(300, 200) 
        # Add tabs 
  
        # Create Force tab 
        self.tabForce.layout = QHBoxLayout()
        tab1ctrl = QVBoxLayout()
        l1 = QLabel() 
        l1.setText("Forcemode") 
        h1 = QHBoxLayout()
        l2 = QLabel() 
        l2.setText("show voltage")
        h1.addWidget(l2)
        h1.addWidget(parent.cbShowVolt)   
        verticalSpacer = QSpacerItem(20, 40, 
                                     QSizePolicy.Minimum, 
                                     QSizePolicy.Expanding)

        tab1ctrl.addWidget(l1)
        tab1ctrl.addLayout(h1)
        tab1ctrl.addItem(verticalSpacer)
        
        
        self.tabForce.layout.addLayout(tab1ctrl)
        self.tabForce.layout.addWidget(parent.mplfor) 
        self.tabForce.setLayout(self.tabForce.layout) 
        
        # Create Env tab 
        self.tabEnv.layout = QHBoxLayout()
        self.tabEnv.layout.addWidget(parent.mplenv)
        self.tabEnv.setLayout(self.tabEnv.layout) 
    
        # Create tVelo tab 
        self.tabVelo.layout = QHBoxLayout()
        tabVeloctrl = QVBoxLayout()
        l1 = QLabel() 
        l1.setText("Velocitymode") 
        h1 = QHBoxLayout()      
        verticalSpacer = QSpacerItem(20, 40, 
                                     QSizePolicy.Minimum, 
                                     QSizePolicy.Expanding)
        tabVeloctrl.addWidget(l1)
        tabVeloctrl.addLayout(h1)
        tabVeloctrl.addItem(verticalSpacer)
        self.tabVelo.layout.addLayout(tabVeloctrl)
        self.tabVelo.layout.addWidget(parent.mplvel) 
        self.tabVelo.setLayout(self.tabVelo.layout) 

        # Create Mass tab 
        self.tabMass.layout = QHBoxLayout()
        tabMassctrl = QVBoxLayout()
        l4b = QLabel() 
        l4b.setText("ref. mass (mg):")    
        verticalSpacer = QSpacerItem(20, 40, 
                                     QSizePolicy.Minimum, 
                                     QSizePolicy.Expanding)
        l4a = QLabel() 
        l4a.setText("mass")
        h1 = QHBoxLayout()

        h1.addWidget(parent.cbExc3sig) #parent.cbMvsZ)
        h1.addWidget(QLabel('Excl. outlier'))
        tabMassctrl.addWidget(l4a)
        tabMassctrl.addLayout(h1) 

        tabMassctrl.addWidget(l4b)
        tabMassctrl.addWidget(parent.sbMass)
        tabMassctrl.addWidget(QLabel("Measured Mass:")) 
        tabMassctrl.addWidget(parent.laMass)
        tabMassctrl.addWidget(parent.laUnc)
        tabMassctrl.addWidget(parent.laUncB)
        tabMassctrl.addWidget(parent.lacov)
        tabMassctrl.addWidget(parent.laNoEnv)
        tabMassctrl.addItem(verticalSpacer)
   
        self.tabMass.layout.addLayout(tabMassctrl)
        self.tabMass.layout.addWidget(parent.mplmass)
        self.tabMass.setLayout(self.tabMass.layout) 
    
        # Create Uncertainty tab 
        self.tabUnc.layout =  QGridLayout()
        
        self.tabUnc.layout.setColumnStretch(parent.laUaMaxRows+1, parent.laUaMaxCols)
  
        la10 = QLabel('Item')
        la11 = QLabel('rel. unc.')
        la12 = QLabel('uncertainty')
        la10.setFont(parent.foBig)
        la11.setFont(parent.foBig)
        la12.setFont(parent.foBig)
 
        la00 = QLabel('measured mass')
        la00.setFont(parent.foBigBold)
        parent.laMass2.setFont(parent.foBigBold)
        parent.laTotUnc.setFont(parent.foBigBold)

        self.tabUnc.layout.addWidget(la00,0,0)
        self.tabUnc.layout.addWidget(parent.laMass2,0,1)
        self.tabUnc.layout.addWidget(parent.laTotUnc,0,2)

        self.tabUnc.layout.addWidget(la10,1,0)
        self.tabUnc.layout.addWidget(la11,1,1)
        self.tabUnc.layout.addWidget(la12,1,2)
        for i in range(parent.laUaMaxRows):
            for j in range(parent.laUaMaxCols):
                self.tabUnc.layout.addWidget(parent.laUa[i][j],i+2,j)

        self.tabUnc.setLayout(self.tabUnc.layout)
        # Create tab6 tab 
        self.tabProfile.layout = QHBoxLayout()
        self.tabProfile.layout.addWidget(parent.mplprofile)
        self.tabProfile.setLayout(self.tabProfile.layout) 
        
        # Create Report tab 
        self.tabReport.layout =  QGridLayout()
        
        self.tabReport.layout.setColumnStretch(len(parent.resultLabels), 2)
        for n,i in enumerate(parent.laResult):
            for m,j in enumerate(i):
                self.tabReport.layout.addWidget(j,n,m)
        
        self.tabReport.setLayout(self.tabReport.layout) 
  
        # Add tabs to widget 
        self.layout.addWidget(self.tabs) 
        self.setLayout(self.layout) 
        self.tabs.addTab(self.tabMass, "Mass") 
        self.tabs.addTab(self.tabForce, "Force") 
        self.tabs.addTab(self.tabEnv, "Environmentals") 
        self.tabs.addTab(self.tabVelo, "Velocity") 
        self.tabs.addTab(self.tabProfile, "Profile")
        self.tabs.addTab(self.tabUnc, "Uncertainty") 
        self.tabs.addTab(self.tabReport, "Report") 



def excepthook(exc_type, exc_value, exc_tb):
    tb = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
    print("error catched!:")
    print("error message:\n", tb)
    fi=open('k2viewer.dbg','w')
    fi.write("error catched!:\n")
    fi.write("error message:\n"+ tb)
    fi.close()
    QApplication.quit()
    # or QtWidgets.QApplication.exit(0)

myappid = 'mycompany.myproduct.subproduct.version' # arbitrary string
ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
app = QApplication(sys.argv)
sys.excepthook = excepthook

window = MainWindow()
window.show()
app.exec()
sys.exit()
