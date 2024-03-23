import os,sys
##https://www.pythonguis.com/tutorials/pyqt-basic-widgets/
##
##https://www.geeksforgeeks.org/pyqt5-qtabwidget/
from PyQt5.QtCore import QSize, Qt
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
    QSizePolicy
)


import mplwidget
import numpy as np
try:
       import pyi_splash
       pyi_splash.update_text('UI Loaded ...')
       pyi_splash.close()
except:
    pass

# Subclass QMainWindow to customize your application's main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        if os.getcwd().startswith('Z:\\BB'):    
            self.bd = r"K:\TableTopWattBalance\KIBB-g2\DATA"
        else:
            self.bd='..\DATA'

        self.clearVariables()
        self.setWindowTitle("Kibb-g2 Viewer")
        self.setFixedSize(QSize(1200, 600))
        self.mytable = QTableWidget(2,1) 
        self.mytable.setFixedWidth(200)
        self.Brefresh = QPushButton("Reload")
        self.mplw1 =mplwidget.MplWidget()
        self.cbForMean = QCheckBox()

        self.tabWidget = MyTabWidget(self) 
 

        layout   = QVBoxLayout()
        hlayout  = QHBoxLayout()
        vlayout1 = QVBoxLayout()
        
        vlayout1.addWidget(self.mytable)
        vlayout1.addWidget(self.Brefresh)
        hlayout.addLayout(vlayout1)
    
        hlayout.addWidget(self.tabWidget)
        layout.addLayout(hlayout)
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)
        self.loadTable()

        self.mytable.clicked.connect(self.on_table_clicked)
        self.Brefresh.clicked.connect(self.loadTable)
        self.cbForMean.clicked.connect(self.plotForce)

    def clearVariables(self):
        self.ton=[]
        self.Uon=[]
        self.tof=[]
        self.Uof=[]


    def loadTable(self):
         self.mytable.clearContents()
         self.mytable.setRowCount(0)
         for s1 in sorted([ f.path for f in os.scandir(self.bd) \
                           if f.is_dir() ],reverse=True):
            yymm=os.path.split(s1)[-1]
            #print(yymm)
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
        self.mplw1.canvas.ax1.clear()
        if self.cbForMean.isChecked():
            if len(self.tof)>=1:
                Uf =np.concatenate(self.Uof)*1e3
                p1,=self.mplw1.canvas.ax1.plot(
                    np.concatenate(self.tof),\
                    Uf-np.mean(Uf),'b.')
            if len(self.ton)>=1:
                Un =np.concatenate(self.Uon)*1e3
                p2,=self.mplw1.canvas.ax1.plot(
                    np.concatenate(self.ton),\
                    Un-np.mean(Un),'r.')
        else:
        
            if len(self.tof)>=1:
                p1,=self.mplw1.canvas.ax1.plot(
                    np.concatenate(self.tof),\
                    np.concatenate(self.Uof)*1e3,'b.')
            if len(self.ton)>=1:
                p2,=self.mplw1.canvas.ax1.plot(
                    np.concatenate(self.ton),\
                    np.concatenate(self.Uon)*1e3,'r.')

        self.mplw1.canvas.ax1.set_xlabel('t/s')    
        self.mplw1.canvas.ax1.set_ylabel('U/mV')
        if self.cbForMean.isChecked():
            if p1!=0 and p2!=0:
                self.mplw1.canvas.ax1.legend((p1,p2),\
                        ('mass off-{0:10.3f} mV'.format(np.mean(Uf)),\
                         'mass on-{0:10.3f} mV'.format(np.mean(Un))))
            elif p1!=0 and p2==0:
                self.mplw1.canvas.ax1.legend((p1),\
                        ('mass off-{0:10.3f} mV'.format(np.mean(Uf))))
            elif p1==0 and p2!=0:
                self.mplw1.canvas.ax1.legend((p2),\
                        ('mass on-{0:10.3f} mV'.format(np.mean(Un))))
        else:
            if p1!=0 and p2!=0:
                self.mplw1.canvas.ax1.legend((p1,p2),\
                        ('mass off','mass on'))
            elif p1!=0 and p2==0:
                self.mplw1.canvas.ax1.legend((p1),\
                        ('mass off'))
            elif p1==0 and p2!=0:
                self.mplw1.canvas.ax1.legend((p2),\
                        ('mass on'))

        self.mplw1.canvas.draw() 
        
        
        
    def readForce(self,bd):
        self.ton=[]
        self.Uon=[]
        self.tof=[]
        self.Uof=[]
        self.mplw1.canvas.ax1.clear() 
        for files in os.listdir(bd):
            if files[7:9].upper()=='OF':
                da  =np.loadtxt(os.path.join(bd,files))
                self.tof.append( da[:,0])
                self.Uof.append( da[:,1])
             
            elif files[7:9].upper()=='ON':
                da  =np.loadtxt(os.path.join(bd,files))
                self.ton.append( da[:,0])
                self.Uon.append( da[:,1])
        currentIndex=self.tabWidget.tabs.currentIndex()
        tat = self.tabWidget.tabs.tabText(currentIndex)
        if tat=='Force':
            self.plotForce()

    def on_table_clicked(self,item):
        oo=self.mytable.item(item.row(), item.column()).text()
        
        filePath = os.path.join(self.bd,oo[0:4],oo[4:6],oo[6])
        if 'Force mode' in os.listdir(filePath):
            self.readForce(os.path.join(filePath,'Force mode'))
   
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
        self.tabs.addTab(self.tab2, "Future") 
        self.tabs.addTab(self.tab3, "Distant Future") 
  
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
        self.tab1.layout.addWidget(parent.mplw1) 
        self.tab1.setLayout(self.tab1.layout) 
  
        # Add tabs to widget 
        self.layout.addWidget(self.tabs) 
        self.setLayout(self.layout) 

app = QApplication(sys.argv)
window = MainWindow()
window.show()
app.exec()
