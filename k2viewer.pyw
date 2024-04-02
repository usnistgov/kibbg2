import os,sys
import traceback

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

from PyQt5.QtGui import QFont
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
    QAbstractSpinBox
)
import mplwidget
import numpy as np
import time
import sqlite3
import k2dataset
try:
       import pyi_splash
       pyi_splash.update_text('UI Loaded ...')
       pyi_splash.close()
except:
    pass

mutex = QMutex()
kda =   k2dataset.k2DataSet(mutex)     


class Worker(QObject):
    finished = pyqtSignal()
    intReady = pyqtSignal(int,int,int)


    @pyqtSlot()
    def procCounter(self): # A slot takes no params
        self.intReady.emit(0,0,0)
        tot = kda.countForceFiles()
        for k in kda.readForce():
            self.intReady.emit(1,k,tot)
        kda.readEnv()
        self.intReady.emit(2,0,0)
        tot = kda.countVeloFiles()
        for k in kda.readVelo():
            self.intReady.emit(3,k,tot) 
        
        self.intReady.emit(99,0,0)
        self.finished.emit()
        
# Subclass QMainWindow to customize your application's main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        #if os.getcwd().startswith('Z:\\BB'):    
            #self.bd = r"K:\TableTopWattBalance\KIBB-g2\DATA"
        #else:
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
                QTableWidgetItem("mass/mg"))
        self.mytable.setHorizontalHeaderItem(2,\
                QTableWidgetItem("unc/mg"))
        self.mytable.setHorizontalHeaderItem(3,\
                QTableWidgetItem("title"))

        self.mytable.setColumnWidth(0, 100)
        self.mytable.setColumnWidth(1, 120)
        self.mytable.setColumnWidth(2, 70)
        self.mytable.setColumnWidth(3, 200)

        self.mytable.setFixedWidth(400)
        self.Brefresh = QPushButton("Reload")
        
        ### The plot windows
        
        self.mplfor  = mplwidget.MplWidget(rightax=True)
        self.mplenv  = mplwidget.MplWidget4()
        self.mplvel  = mplwidget.MplWidget(rightax=True)
        self.mplmass = mplwidget.MplWidget(rightax=True)
        
        ### check and spin boxes
        
        self.cbForMean  = QCheckBox()
        self.sbOrder    = QSpinBox()
        self.sbSupports = QSpinBox()
        self.sbMass     = QDoubleSpinBox()
        self.sbOrder.setValue(4)
        self.sbOrder.setMinimum(1)
        self.sbOrder.setMaximum(10)
        
        self.sbSupports.setValue(5)
        self.sbSupports.setMinimum(2)
        self.sbSupports.setMaximum(20)
   
        self.sbMass.setMinimumWidth(100)
        self.sbMass.setMinimum(0)
        self.sbMass.setMaximum(99999)
        self.sbMass.setDecimals(4)
        self.sbMass.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.sbMass.setKeyboardTracking(False)
        
        ### global labels
        
        self.sblabel  = QLabel("Click on a run") 
        self.lares    = QLabel("") 
        self.laUncA   = QLabel("n/a")
        self.laUncTot = QLabel("n/a")
        
        
        ### Status bar

        self.tabWidget = MyTabWidget(self) 
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        #self.statusBar.setStyleSheet("border :3px solid black;") 
        
        #self.sblabel.setStyleSheet("border :2px solid blue;") 
        self.statusBar.addPermanentWidget(self.sblabel) 
        
        
        
        self.statusBar.showMessage('Welcome to k2viewer',5000)

        layout   = QVBoxLayout()
        hlayout  = QHBoxLayout()
        vlayout1 = QVBoxLayout()
        vlayout2 = QVBoxLayout()
        
    
        
        widget = QWidget()
        widget.setLayout(layout)
        layout.addLayout(hlayout)
        hlayout.addLayout(vlayout1)
        hlayout.addLayout(vlayout2)
        vlayout2.addWidget(self.tabWidget)
        vlayout1.addWidget(self.mytable)
        vlayout1.addWidget(self.Brefresh)
     
        self.setCentralWidget(widget)
        self.loadTable()

        self.mytable.clicked.connect(self.on_table_clicked)
        self.Brefresh.clicked.connect(self.loadTable)
        self.cbForMean.clicked.connect(self.plotForce)
        self.tabWidget.tabs.currentChanged.connect(self.replot)
        self.sbOrder.valueChanged.connect(self.recalcvelo)
        self.sbSupports.valueChanged.connect(self.recalcvelo)
        self.sbMass.valueChanged.connect(self.gotmassval)
        

    

    def loadTable(self):
         self.mytable.clearContents()
         self.mytable.setRowCount(0)
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
                            letter=os.path.split(s3)[-1]
                            if len(letter)==1:
                                run =yymm+day+letter
                                row_number = self.mytable.rowCount()
                                self.mytable.insertRow(row_number)
                                self.mytable.setItem(row_number,0,\
                                 QTableWidgetItem(str(run)))
                                dbentry = cursor.execute("""SELECT run,value,
                                uncertainty,title FROM k2data 
                                WHERE run="{0}";""".format(run)).fetchall()
                                if len(dbentry)==0:
                                    continue
                                dbentry=dbentry[0]
                                #print(dbentry)
                                if dbentry[1]>-9e96:
                                    nitem =QTableWidgetItem('{0:,.4f}'.format(dbentry[1]*1e3) )
                                else:
                                    nitem =QTableWidgetItem('n/a')
                                nitem.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
                                self.mytable.setItem(row_number,1,nitem)
                                if dbentry[2]>-9e96:
                                    nitem =QTableWidgetItem('{0:6.4f}'.format(dbentry[2]*1e3))
                                else:
                                    nitem =QTableWidgetItem('n/a')                                
                                nitem.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
                                self.mytable.setItem(row_number,2,nitem)
                                self.mytable.setItem(row_number,3,QTableWidgetItem(dbentry[3]))

         connection.close()
         
    def plotForce(self):
        p1=0
        p2=0
        self.mplfor.canvas.ax1.clear()
        if self.cbForMean.isChecked():
            mof = np.mean(kda.cIof)
            mon = np.mean(kda.cIon)
        else:
            mof=0
            mon=0            
        mutex.lock()
        tmul,tla = kda.tmul()
        if kda.hasOff():
            p1,=self.mplfor.canvas.ax1.plot(
                kda.ctof*tmul,\
                kda.cIof-mof,'b.')
        if kda.hasOn():
            p2,=self.mplfor.canvas.ax1.plot(
                kda.cton*tmul,\
                    kda.cIon-mon,'r.')
        if len(kda.Iaon)>1:
            p3,=self.mplfor.canvas.ax1.plot(
                kda.taon*tmul,\
                    kda.Iaon-mon,'m.')
        if len(kda.Iaof)>1:
            p4,=self.mplfor.canvas.ax1.plot(
                kda.taof*tmul,\
                    kda.Iaof-mof,'c.')
            
        mutex.unlock()
        if self.cbForMean.isChecked():
            if p1!=0 and p2!=0:
                self.mplfor.canvas.ax1.legend((p1,p2),\
                        ('mass off-{0:12.6f} uA (left)'.format(mof),\
                          'mass on-{0:12.6f} uA (right)'.format(mon)))
            elif p1!=0 and p2==0:
                self.mplfor.canvas.ax1.legend((p1),\
                        ('mass off-{0:12.6f} uA (left)'.format(mof)))
            elif p1==0 and p2!=0:
                self.mplfor.canvas.ax1.legend((p2),\
                        ('mass on-{0:12.6f} uA (right)'.format(mon)))
        else:
            if p1!=0 and p2!=0:
                self.mplfor.canvas.ax1.legend((p1,p2),\
                        ('mass off (left)','mass on (right)'))
            elif p1!=0 and p2==0:
                self.mplfor.canvas.ax1.legend((p1),\
                        ('mass off (left)'))
            elif p1==0 and p2!=0:
                self.mplfor.canvas.ax1.legend((p2),\
                        ('mass on (right)'))

        self.mplmass.canvas.bx1.ticklabel_format(useOffset=False)
        self.mplfor.canvas.ax1.ticklabel_format(useOffset=False)
        self.mplfor.canvas.ax1.set_xlabel(tla)    
        self.mplfor.canvas.ax1.set_ylabel('I(off)/uA')
        self.mplfor.canvas.bx1.set_ylabel('I(on)/uA')

        self.mplfor.canvas.draw() 
            
    def plotEnv(self):
        self.mplenv.canvas.ax1.clear()
        self.mplenv.canvas.ax2.clear()
        self.mplenv.canvas.ax3.clear()
        self.mplenv.canvas.ax4.clear()
        if kda.hasEnv==False:
            return
        mutex.lock()
        tmul,tla = kda.tmul()
        self.mplenv.canvas.ax1.plot(kda.edata [:,0]*tmul,kda.edata[:,1],'r.')
        self.mplenv.canvas.ax2.plot(kda.edata [:,0]*tmul,kda.edata[:,2],'b.')
        self.mplenv.canvas.ax3.plot(kda.edata [:,0]*tmul,kda.edata[:,3],'g.')
        self.mplenv.canvas.ax4.plot(kda.edata [:,0]*tmul,kda.edata[:,4],'m.')
        mutex.unlock()
        self.mplenv.canvas.ax1.ticklabel_format(useOffset=False)
        self.mplenv.canvas.ax2.ticklabel_format(useOffset=False)
        self.mplenv.canvas.ax3.ticklabel_format(useOffset=False)
        self.mplenv.canvas.ax4.ticklabel_format(useOffset=False)
        self.mplenv.canvas.ax1.xaxis.set_ticklabels([])
        self.mplenv.canvas.ax2.xaxis.set_ticklabels([])
        self.mplenv.canvas.ax3.xaxis.set_ticklabels([])
        self.mplenv.canvas.ax1.tick_params(axis='x',direction='inout')
        self.mplenv.canvas.ax2.tick_params(axis='x',direction='inout')
        self.mplenv.canvas.ax3.tick_params(axis='x',direction='inout')
        self.mplenv.canvas.ax1.set_ylabel('rel. humid (%)')
        self.mplenv.canvas.ax2.set_ylabel('press. (hPa)')
        self.mplenv.canvas.ax3.set_ylabel('temp (degC)')
        self.mplenv.canvas.ax4.set_ylabel('air dens. (kg/m^3)')
        self.mplenv.canvas.ax4.set_xlabel(tla)  
        self.mplenv.canvas.draw() 


    def plotVelocity(self):
        self.mplvel.canvas.ax1.clear()
        mutex.lock()
        tmul,tla = kda.tmul()
        
        if len(kda.vblv)>=2:
            p1,=self.mplvel.canvas.ax1.plot(
                kda.vblt*tmul,\
                np.abs(kda.vblv)*1e3,'b.')
            p2,=self.mplvel.canvas.ax1.plot(
                kda.sblt*tmul,\
                np.abs(kda.sblv)*1e3,'r-')
        mutex.unlock()
        be,en =self.mplvel.canvas.ax1.get_ylim()
        self.mplvel.canvas.draw() 
        me =0.5*(be+en)
        self.mplvel.canvas.bx1.set_ylim((be/me-1)*1e6,(en/me-1)*1e6)
        self.mplvel.canvas.ax1.set_xlabel(tla)    
        self.mplvel.canvas.ax1.set_ylabel('Bl/Tm')
        self.mplvel.canvas.bx1.set_ylabel('rel. change /ppm')
        self.mplvel.canvas.bx1.ticklabel_format(useOffset=False)
        self.mplvel.canvas.ax1.ticklabel_format(useOffset=False)
        self.mplvel.canvas.draw() 

    def plotMass(self):
        if kda.hasresult==False:
            return
        self.mplmass.canvas.ax1.clear()
        mutex.lock()
        tmul,tla = kda.tmul()
        lines,cols = np.shape(kda.Fresult)
        if lines==0: 
            return
        mean = kda.mass
        sig = kda.massunc
        #mean = np.sum(kda.Fresult[:,1]/kda.Fresult[:,2]**2)/np.sum(1/kda.Fresult[:,2]**2)
        #sig = np.sqrt(1/np.sum(1/kda.Fresult[:,2]**2))
        self.mplmass.canvas.ax1.errorbar(
            kda.Fresult[:,0]*tmul,\
            kda.Fresult[:,1],kda.Fresult[:,2],fmt='bo')
        #mint=np.min(kda.Fresult[:,0]*tmul)
        #maxt=np.max(kda.Fresult[:,0]*tmul)
        mutex.unlock()
        be,en =self.mplmass.canvas.ax1.get_xlim()
        self.mplmass.canvas.ax1.plot((be,en),(mean,mean),c='k',linestyle='dashed',lw=2)
        self.mplmass.canvas.ax1.plot((be,en),(mean+sig,mean+sig),\
                                     c='r',linestyle='dashdot')
        self.mplmass.canvas.ax1.plot((be,en),(mean-sig,mean-sig),\
                                     c='r',linestyle='dashdot')
        self.mplmass.canvas.ax1.fill_between((be,en), (mean-sig,mean-sig),\
                                (mean+sig,mean+sig),color='r', alpha=0.2)
        self.mplmass.canvas.ax1.set_xlim(be,en)
        if kda.hasRefMass:
            be,en =self.mplmass.canvas.ax1.get_xlim()
            self.mplmass.canvas.ax1.plot((be,en),
                                         (kda.refMass,kda.refMass),c='m',
                                         linestyle='dotted',lw=4)
            self.mplmass.canvas.ax1.set_xlim(be,en)
            
        be,en =self.mplmass.canvas.ax1.get_ylim()
            
        #self.mplmass.canvas.draw() 
        me =0.5*(be+en)
        if kda.hasRefMass:
            me =kda.refMass
        self.mplmass.canvas.bx1.set_ylim((be/me-1)*1e6,(en/me-1)*1e6)
        self.mplmass.canvas.ax1.set_xlabel(tla)    
        self.mplmass.canvas.ax1.set_ylabel('mass /mg')
        self.mplmass.canvas.bx1.set_ylabel('rel. change /ppm')
        self.mplmass.canvas.bx1.ticklabel_format(useOffset=False)
        self.mplmass.canvas.ax1.ticklabel_format(useOffset=False)
        self.mplmass.canvas.draw() 
        
        
    def plotUnc(self):
        #print('uu: ',kda.relUncTypeA ,kda.relUncTot )
        if kda.relUncTypeA <9e99:
            self.laUncA.setText('{0:3.1f}'.format(kda.relUncTypeA))
        else:
            self.laUncA.setText("n/a")
        if kda.relUncTot<9e99:
           self.laUncTot.setText('{0:3.1f}'.format(kda.relUncTot))
        else:
            self.laUncTot.setText("n/a")


    def gotmassval(self,val):
        kda.setRefMass(val)
        self.plotMass()
        


    def recalcvelo(self):
        order = int(self.sbOrder.value() )
        supports = int(self.sbSupports.value() )
        kda.fitVelo(order,supports)
        #self.plotVelocity()
        kda.calcforce()
        self.updateTable()
        self.replot()




    def replot(self):
        currentIndex=self.tabWidget.tabs.currentIndex()
        tat = self.tabWidget.tabs.tabText(currentIndex)
        if tat=='Force':
            self.plotForce()
        elif tat=='Environmentals':
            self.plotEnv()
        elif tat=='Velocity':
            self.plotVelocity()
        elif tat=='Mass':
            self.plotMass()
        elif tat=='Uncertainty':
            self.plotUnc()
    
    def upStatus(self,msg,dstime):
        now =time.time()
        if now-self.statust>1:
            self.statusBar.showMessage(msg,1000)
            self.statust=time.time()
   
            
    def updateTable(self):
        #if kda.hasresult:# and np.isnan(kda.mass)==False:
        if np.isnan(kda.mass):
            nitem =QTableWidgetItem('n/a') 
        else:
            nitem =QTableWidgetItem('{0:,.4f}'.format(kda.mass) )
        nitem.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
        self.mytable.setItem(self.calcrow,1,nitem)
        if np.isnan(kda.massunc):
            nitem =QTableWidgetItem('n/a') 
        else:
            nitem =QTableWidgetItem('{0:6.4f}'.format(kda.massunc))
        nitem.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
        self.mytable.setItem(self.calcrow,2,nitem)
        self.mytable.setItem(self.calcrow,3,QTableWidgetItem(kda.title))
        mass = kda.mass
        massunc = kda.massunc
        if np.isnan(mass):
            mass=-9e99
        if np.isnan(massunc):
            massunc=-9e99

        mycmd ="""
        replace into k2data (run,value, uncertainty,title)
        values  ("{0}",{1},{2},"{3}");""".\
           format(self.runid,mass/1000,massunc/1000,kda.title)
        connection = sqlite3.connect('k2viewer.db')
        cursor = connection.cursor()
        cursor.execute(mycmd)
        connection.commit()
        connection.close()    
   
    
   
    def readStatus(self,ix,cur,tot):
        if ix==0: 
            self.start=time.time()
            return
        if ix==1:
            self.upStatus('Reading Force Files {0}/{1}'.\
                                       format(cur,tot),100)
            if cur>2:
                self.replot()
        elif ix==2:
            self.upStatus('Reading Environmentals',1000)
        elif ix==3:
            self.upStatus('Reading Velo Files {0}/{1}'.\
                                       format(cur,tot),100)
        elif ix==99:
            dt = time.time()- self.start 
            kda.fitVelo( int(self.sbOrder.value() ),-1)
            kda.calcforce()

            self.upStatus('Reading completed in {0:3.1f} seconds'.\
                                       format(dt),2000)
            self.sblabel.setText(kda.title)
            self.sbSupports.setValue(kda.nrknots)
            self.idle=True
            self.updateTable()
            self.replot()
            
            
            
            
        

    def on_table_clicked(self,item):
        if self.idle:
            self.idle=False
            self.sbMass.clear()
            self.calcrow = item.row()
            self.runid =self.mytable.item(item.row(), 0).text()
            filePath = os.path.join(self.bd,self.runid[0:4],self.runid[4:6],self.runid[6])
            kda.setbd0(filePath)
            self.obj = Worker()  # no parent!
            self.thread = QThread()  # no parent!
            self.obj.intReady.connect(self.readStatus)
            self.obj.moveToThread(self.thread)
            self.obj.finished.connect(self.thread.quit)
            self.thread.started.connect(self.obj.procCounter)
            self.thread.start()
        else:
            self.upStatus('Wait till last run is processed',4000)
         
class MyTabWidget(QWidget): 
    def __init__(self, parent): 
        super(QWidget, self).__init__(parent) 
  
        self.layout = QVBoxLayout(self) 

        # Initialize tab screen 
        self.tabs = QTabWidget() 
        self.tab1 = QWidget() 
        self.tab2 = QWidget() 
        self.tab3 = QWidget() 
        self.tab4 = QWidget() 
        self.tab5 = QWidget() 
        self.tabs.resize(300, 200) 
  
        # Add tabs 
  
        # Create first tab 
        self.tab1.layout = QHBoxLayout()
        tab1ctrl = QVBoxLayout()
        l1 = QLabel() 
        l1.setText("Forcemode") 
        h1 = QHBoxLayout()
        l2 = QLabel() 
        l2.setText("subtract mean")
        h1.addWidget(l2)
        h1.addWidget(parent.cbForMean)
        
        verticalSpacer = QSpacerItem(20, 40, 
                                     QSizePolicy.Minimum, 
                                     QSizePolicy.Expanding)

        tab1ctrl.addWidget(l1)
        tab1ctrl.addLayout(h1)
        tab1ctrl.addItem(verticalSpacer)
        
        
        self.tab1.layout.addLayout(tab1ctrl)
        self.tab1.layout.addWidget(parent.mplfor) 
        self.tab1.setLayout(self.tab1.layout) 
        
        # Create sceond tab 
        self.tab2.layout = QHBoxLayout()
        self.tab2.layout.addWidget(parent.mplenv)
        self.tab2.setLayout(self.tab2.layout) 
    
        # Create third tab 
        self.tab3.layout = QHBoxLayout()
        tab3ctrl = QVBoxLayout()
        l1 = QLabel() 
        l1.setText("Velocitymode") 
        h1 = QHBoxLayout()
        l2 = QLabel() 
        l2.setText("order")
        h1.addWidget(l2)
        h1.addWidget(parent.sbOrder)
        
        h2 = QHBoxLayout()
        l3 = QLabel() 
        l3.setText("supports")
        h2.addWidget(l3)
        h2.addWidget(parent.sbSupports)

        
        verticalSpacer = QSpacerItem(20, 40, 
                                     QSizePolicy.Minimum, 
                                     QSizePolicy.Expanding)



        tab3ctrl.addWidget(l1)
        tab3ctrl.addLayout(h1)
        tab3ctrl.addLayout(h2)
        tab3ctrl.addItem(verticalSpacer)
        
        self.tab3.layout.addLayout(tab3ctrl)
        self.tab3.layout.addWidget(parent.mplvel) 
        self.tab3.setLayout(self.tab3.layout) 
        

        # Create forth tab 
        self.tab4.layout = QHBoxLayout()
        tab4ctrl = QVBoxLayout()
        
        #h4 = QHBoxLayout()
        l4b = QLabel() 
        l4b.setText("ref. mass (mg):")
        #h4.addWidget(l4b)
        #h4.addWidget(parent.sbMass)
        
        verticalSpacer = QSpacerItem(20, 40, 
                                     QSizePolicy.Minimum, 
                                     QSizePolicy.Expanding)
        l4a = QLabel() 
        l4a.setText("mass")

        tab4ctrl.addWidget(l4a)
        tab4ctrl.addWidget(l4b)
        tab4ctrl.addWidget(parent.sbMass)
        tab4ctrl.addItem(verticalSpacer)
        
        self.tab4.layout.addLayout(tab4ctrl)
        self.tab4.layout.addWidget(parent.mplmass)
        self.tab4.setLayout(self.tab4.layout) 
    
        # Create fifth tab 
        self.tab5.layout =  QGridLayout()
        self.tab5.layout.setColumnStretch(12, 3)
        self.tab5.layout.addWidget(QLabel('Item'),0,0)
        self.tab5.layout.addWidget(QLabel('rel. unc/ppm'),0,1)
        self.tab5.layout.addWidget(QLabel('unc/ug'),0,2)
        
        self.tab5.layout.addWidget(QLabel('Resistor'),1,0)
        self.tab5.layout.addWidget(QLabel('1.4'),1,1)
        
        self.tab5.layout.addWidget(QLabel('Voltmeter'),2,0)
        self.tab5.layout.addWidget(QLabel('1.0'),2,1)
 
        self.tab5.layout.addWidget(QLabel('mass position'),3,0)
        self.tab5.layout.addWidget(QLabel('1.0'),3,1)

        self.tab5.layout.addWidget(QLabel('g'),4,0)
        self.tab5.layout.addWidget(QLabel('2.0'),4,1)
        
        self.tab5.layout.addWidget(QLabel('verticality'),5,0)
        self.tab5.layout.addWidget(QLabel('0.5'),5,1)

        self.tab5.layout.addWidget(QLabel('statistical'),6,0)
        self.tab5.layout.addWidget(parent.laUncA,6,1)

        ltot = QLabel('total')
        self.tab5.layout.addWidget(ltot,7,0)
        self.tab5.layout.addWidget(parent.laUncTot,7,1)
        
        myFont=QFont()
        myFont.setBold(True)
        parent.laUncTot.setFont(myFont)
        ltot.setFont(myFont)


        self.tab5.setLayout(self.tab5.layout)
            
 
  
        # Add tabs to widget 
        self.layout.addWidget(self.tabs) 
        self.setLayout(self.layout) 
        self.tabs.addTab(self.tab1, "Force") 
        self.tabs.addTab(self.tab2, "Environmentals") 
        self.tabs.addTab(self.tab3, "Velocity") 
        self.tabs.addTab(self.tab4, "Mass") 
        self.tabs.addTab(self.tab5, "Uncertainty") 



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


app = QApplication(sys.argv)
sys.excepthook = excepthook

window = MainWindow()
window.show()
app.exec()
sys.exit()
