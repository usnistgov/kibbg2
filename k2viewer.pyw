import os,sys
##https://www.pythonguis.com/tutorials/pyqt-basic-widgets/
##
##https://www.geeksforgeeks.org/pyqt5-qtabwidget/
#https://realpython.com/python-pyqt-qthread/
#https://stackoverflow.com/questions/6783194/background-thread-with-qthread-in-pyqt
#
from PyQt5.QtCore import QSize, QMutex, QObject, QThread, pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDateEdit,
    QDateTimeEdit,
    QDial,
    QDoubleSpinBox,
    QFontComboBox,
    QLabel,
    QLCDNumber,
    QLineEdit,
    QMainWindow,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QSlider,
    QSpinBox,
    QTimeEdit,
    QVBoxLayout,
    QHBoxLayout,
    QWidget,
    QTreeView,
    QFileSystemModel,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QSpacerItem,
    QSizePolicy,
    QStatusBar
)
from pathlib import Path
import k2tools
import mplwidget
import numpy as np
import time
from threading import *
try:
       import pyi_splash
       pyi_splash.update_text('UI Loaded ...')
       pyi_splash.close()
except:
    pass

mutex = QMutex()

class k2DataSet():
    def __init__(self):
        self.bd0 = ''
        self.clearForce()
        self.clearEnv()
        
    def setbd0(self,bd0):
        self.bd0= bd0
    def readAll(self):
        self.readForce()
        self.readEnv()

    def clearForce(self):
        mutex.lock()
        self.ton=[]
        self.Uon=[]
        self.tof=[]
        self.Uof=[]
        mutex.unlock()

    def clearEnv(self):
        mutex.lock()
        self.edata = []
        mutex.unlock()

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
        if self.bd0=='':
            self.clearForce()
            return            
        if 'Force mode' not in os.listdir(self.bd0):
            self.clearForce()
            return
        bd = os.path.join(self.bd0,'Force mode')
        ton=[]
        Uon=[]
        tof=[]
        Uof=[]
        for files in os.listdir(bd):
            if files[7:9].upper()=='OF':
                da  =np.loadtxt(os.path.join(bd,files))
                tof.append( da[:,0])
                Uof.append( da[:,1])
             
            elif files[7:9].upper()=='ON':
                da  =np.loadtxt(os.path.join(bd,files))
                ton.append( da[:,0])
                Uon.append( da[:,1])
            if co==0:
                self.fohdr,self.ti = self.readForHdr(bd,files)
                if 'Resistor (ohm)' in self.fohdr:
                    self.R = self.fohdr['Resistor (ohm)']
                else:
                    self.R =10000.032685
            mutex.lock()
            self.Uon = Uon
            self.Uof = Uof
            self.tof = tof
            self.ton = ton
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
            mutex.unlock()
            yield co
            co=co+1
            
    def hasOn(self):
        return len(self.Uon)>0
        
    def hasOff(self):
        return len(self.Uof)>0

    def hasEnv(self):
        return len(self.edata)>0

    def readEnv(self):
        if self.bd0=='':
            self.cleaEnv()
            return            
        if 'PRTData.dat' not in os.listdir(self.bd0):
            self.cleaEnv()
            return
        fn =os.path.join(self.bd0,'PRTData.dat')
        da = np.loadtxt(fn,skiprows=1,usecols=[1,2,3])
        da =np.vstack((np.arange(len(da[:,0])).T*10,da.T)).T
        dens = k2tools.airDensity(da[:,3],da[:,2],da[:,1])
        mutex.lock()
        self.edata = np.vstack((da.T,dens)).T
        mutex.unlock()
        
    def readForHdr(self,bd,fn):
        with open(os.path.join(bd,fn)) as input_file:
            head = [next(input_file) for _ in range(6)]
        fields=[o.strip() for o in head[3].split('|')]
        values=[float(f) for f  in head[4].split()[1:]]
        return dict(zip(fields,values)),head[1][1:]

    def clearVelo(self):
        mutex.lock()
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
 
        mutex.unlock()
        
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
                mutex.lock()
                self.vt1 = np.r_[self.vt1,data[:,0]]
                self.vt2 = np.r_[self.vt2,data[:,1]]
                self.vz1 = np.r_[self.vz1,data[:,2]]
                self.vz2 = np.r_[self.vz2,data[:,3]]
                self.vv = np.r_[self.vv,data[:,5]]
                self.vV = np.r_[self.vV,data[:,6]]
                self.vS = np.r_[self.vS,np.ones(len(data[:,4]))*Sco]
                self.vt =0.5*(self.vt1+self.vt2)
                self.vz =0.5*(self.vz1+self.vz2)
                mutex.unlock()
                #print(f,min(data[:,0]))
                Sco+=1
                yield Sco
    
    def fitVelo(self,order): 
        if len(self.vt1)<20:
            return
        zmin = min(self.vz)
        zmax = max(self.vz)
        ft,fv,fC2,fNDF = \
        k2tools.FitLikeACanadianOrthoMaster(
            self.vt,self.vz,self.vv,self.vV,self.vS,zmin,zmax,order)
        mutex.lock()
        self.ft = ft
        self.fv = fv
        self.fC2 = fC2
        self.fNDF = fNDF        
        mutex.unlock()
        
        


kda =   k2DataSet()     


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
        
        kda.fitvelo()
        self.intReady.emit(99,0,0)
        self.finished.emit()
        
# Subclass QMainWindow to customize your application's main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        if os.getcwd().startswith('Z:\\BB'):    
            self.bd = r"K:\TableTopWattBalance\KIBB-g2\DATA"
        else:
            self.bd='..\DATA'

        self.thread = QThread()
        self.setWindowTitle("Kibb-g2 Viewer")
        self.setFixedSize(QSize(1200, 600))
        
        self.mytable = QTableWidget(2,1) 
        self.mytable.setFixedWidth(200)
        self.Brefresh = QPushButton("Reload")
        
        ### The plot windows
        
        self.mplfor =mplwidget.MplWidget()
        self.mplenv =mplwidget.MplWidget4()
        self.mplvel =mplwidget.MplWidget()
        
        self.cbForMean = QCheckBox()

        self.tabWidget = MyTabWidget(self) 
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
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
    

    def loadTable(self):
         self.mytable.clearContents()
         self.mytable.setRowCount(0)
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
                                row_number = self.mytable.rowCount()
                                self.mytable.insertRow(row_number)
                                self.mytable.setItem(row_number,0,\
                                 QTableWidgetItem(str(yymm+day+letter)))
         
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
        if kda.hasOff():
            p1,=self.mplfor.canvas.ax1.plot(
                kda.ctof,\
                kda.cIof-mof,'b.')
        if kda.hasOn():
            p2,=self.mplfor.canvas.ax1.plot(
                kda.cton,\
                    kda.cIon-mon,'r.')
        mutex.unlock()
        self.mplfor.canvas.ax1.set_xlabel('t/s')    
        self.mplfor.canvas.ax1.set_ylabel('I/uA')
        if self.cbForMean.isChecked():
            if p1!=0 and p2!=0:
                self.mplfor.canvas.ax1.legend((p1,p2),\
                        ('mass off-{0:12.6f} uA'.format(mof),\
                          'mass on-{0:12.6f} uA'.format(mon)))
            elif p1!=0 and p2==0:
                self.mplfor.canvas.ax1.legend((p1),\
                        ('mass off-{0:12.6f} uA'.format(mof)))
            elif p1==0 and p2!=0:
                self.mplfor.canvas.ax1.legend((p2),\
                        ('mass on-{0:12.6f} uA'.format(mon)))
        else:
            if p1!=0 and p2!=0:
                self.mplfor.canvas.ax1.legend((p1,p2),\
                        ('mass off','mass on'))
            elif p1!=0 and p2==0:
                self.mplfor.canvas.ax1.legend((p1),\
                        ('mass off'))
            elif p1==0 and p2!=0:
                self.mplfor.canvas.ax1.legend((p2),\
                        ('mass on'))

        self.mplfor.canvas.draw() 
            
    def plotEnv(self):
        self.mplenv.canvas.ax1.clear()
        self.mplenv.canvas.ax2.clear()
        self.mplenv.canvas.ax3.clear()
        self.mplenv.canvas.ax4.clear()
        if kda.hasEnv()==False:
            return
        mutex.lock()
        self.mplenv.canvas.ax1.plot(kda.edata [:,0],kda.edata[:,1],'r.')
        self.mplenv.canvas.ax2.plot(kda.edata [:,0],kda.edata[:,2],'b.')
        self.mplenv.canvas.ax3.plot(kda.edata [:,0],kda.edata[:,3],'g.')
        self.mplenv.canvas.ax4.plot(kda.edata [:,0],kda.edata[:,4],'m.')
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


    def plotVelocity(self):
        self.mplvel.canvas.ax1.clear()
        mutex.lock()
        if len(kda.ft)>2:
            p1,=self.mplvel.canvas.ax1.plot(
                kda.ft,\
                np.abs(kda.fv)*1e3,'b.')
        mutex.unlock()
        self.mplfor.canvas.ax1.set_xlabel('t/s')    
        self.mplfor.canvas.ax1.set_ylabel('Bl/Tm')
        self.mplfor.canvas.draw() 




    def replot(self):
        currentIndex=self.tabWidget.tabs.currentIndex()
        tat = self.tabWidget.tabs.tabText(currentIndex)
        if tat=='Force':
            self.plotForce()
        elif tat=='Environmentals':
            self.plotEnv()
        elif tat=='Velocity':
            self.plotVelocity()
            

    def readStatus(self,ix,cur,tot):
        if ix==0: 
            self.start=time.time()
            return
        if ix==1:
            self.statusBar.showMessage('Reading Force Files {0}/{1}'.\
                                       format(cur,tot),100)
            if cur>2:
                self.replot()
        elif ix==2:
            self.statusBar.showMessage('Reading Environmentals',1000)
        elif ix==3:
            self.statusBar.showMessage('Reading Velo Files {0}/{1}'.\
                                       format(cur,tot),100)
        elif ix==99:
            dt = time.time()- self.start 
            self.statusBar.showMessage('Reading completed in {0:3.1f} seconds'.\
                                       format(dt),2000)
            
        #self.replot()
            
        

    def on_table_clicked(self,item):
        oo=self.mytable.item(item.row(), item.column()).text()
        filePath = os.path.join(self.bd,oo[0:4],oo[4:6],oo[6])
        kda.setbd0(filePath)
        self.obj = Worker()  # no parent!
        self.thread = QThread()  # no parent!
        # 2 - Connect Worker`s Signals to Form method slots to post data.
        self.obj.intReady.connect(self.readStatus)
        # 3 - Move the Worker object to the Thread object
        self.obj.moveToThread(self.thread)
        # 4 - Connect Worker Signals to the Thread slots
        self.obj.finished.connect(self.thread.quit)
        # 5 - Connect Thread started signal to Worker operational slot method
        self.thread.started.connect(self.obj.procCounter)

       # * - Thread finished signal will close the app if you want!
       #self.thread.finished.connect(app.exit)

       # 6 - Start the thread
        self.thread.start()
      
class MyTabWidget(QWidget): 
    def __init__(self, parent): 
        super(QWidget, self).__init__(parent) 
        self.layout = QVBoxLayout(self) 
  
        # Initialize tab screen 
        self.tabs = QTabWidget() 
        self.tab1 = QWidget() 
        self.tab2 = QWidget() 
        self.tab3 = QWidget() 
        self.tabs.resize(300, 200) 
  
        # Add tabs 
        self.tabs.addTab(self.tab1, "Force") 
        self.tabs.addTab(self.tab2, "Environmentals") 
        self.tabs.addTab(self.tab3, "Velocity") 
  
        # Create first tab 
        self.tab1.layout = QHBoxLayout(self)
        tab1ctrl = QVBoxLayout(self)
        
        l1 = QLabel() 
        l1.setText("Forcemode") 
        h1 = QHBoxLayout(self)
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
        self.tab2.layout = QHBoxLayout(self)
        self.tab2.layout.addWidget(parent.mplenv)
        self.tab2.setLayout(self.tab2.layout) 
    
        # Create third tab 
        self.tab3.layout = QHBoxLayout(self)
        tab3ctrl = QVBoxLayout(self)
        
        l1 = QLabel() 
        l1.setText("Velocitymode") 
        h1 = QHBoxLayout(self)
        l2 = QLabel() 
        l2.setText("order")
        h1.addWidget(l2)
        #h1.addWidget(parent.cbForMean)
        
        verticalSpacer = QSpacerItem(20, 40, 
                                     QSizePolicy.Minimum, 
                                     QSizePolicy.Expanding)

        tab3ctrl.addWidget(l1)
        tab3ctrl.addLayout(h1)
        tab3ctrl.addItem(verticalSpacer)
        
        
        
        
        self.tab3.layout.addLayout(tab3ctrl)
        self.tab3.layout.addWidget(parent.mplvel) 
        self.tab3.setLayout(self.tab3.layout) 
        

  
        # Add tabs to widget 
        self.layout.addWidget(self.tabs) 
        self.setLayout(self.layout) 

app = QApplication(sys.argv)
window = MainWindow()
window.show()
app.exec()
