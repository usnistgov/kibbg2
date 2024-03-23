import sys,os
#https://www.pythonguis.com/tutorials/pyqt-basic-widgets/
#
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
)
from PyQt5.QtCore import QRect

import mplwidget
from numpy import loadtxt

# Subclass QMainWindow to customize your application's main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        #self.bd = r"K:\TableTopWattBalance\KIBB-g2\DATA"
        self.bd='..\DATA'

        self.setWindowTitle("Kibb-g2 Viewer")
        self.setFixedSize(QSize(1200, 600))
        self.mytree = QTreeView() 
        self.mytree.setFixedWidth(200)
        self.Brefresh = QPushButton("Reload")
        self.mplw1 =mplwidget.MplWidget()
 

        layout   = QVBoxLayout()
        hlayout  = QHBoxLayout()
        vlayout1 = QVBoxLayout()
        
        vlayout1.addWidget(self.mytree)
        vlayout1.addWidget(self.Brefresh)
        hlayout.addLayout(vlayout1)
    
        hlayout.addWidget(self.mplw1)
        layout.addLayout(hlayout)
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)
        self.loadTree()

        self.mytree.clicked.connect(self.on_treeView_clicked)
        self.Brefresh.clicked.connect(self.loadTree)



    def loadTree(self):
         self.model = QFileSystemModel()
         self.model.setRootPath(self.bd)
         self.model.sort(0,1)
         
         self.mytree.setModel(self.model)
         self.mytree.setRootIndex(self.model.index(self.bd))
         self.mytree.header().resizeSection(0,200)
         self.mytree.hideColumn(1)
         self.mytree.hideColumn(2)
         self.mytree.hideColumn(3)
         
    def plotForce(self,bd):
        self.mplw1.canvas.ax1.clear() 
        for files in os.listdir(bd):
            if files[7:9].upper()=='OF':
                da  =loadtxt(os.path.join(bd,files))
                p1,=self.mplw1.canvas.ax1.plot(da[:,0],da[:,1]*1e3,'b.')
            elif files[7:9].upper()=='ON':
                da  =loadtxt(os.path.join(bd,files))
                p2,=self.mplw1.canvas.ax1.plot(da[:,0],da[:,1]*1e3,'r.')
                
        self.mplw1.canvas.ax1.set_xlabel('t/s')    
        self.mplw1.canvas.ax1.set_ylabel('U/mV')    
        
        self.mplw1.canvas.ax1.legend((p1,p2),\
        ('mass off','mass on'))
        self.mplw1.canvas.draw() 
        

    def on_treeView_clicked(self,index):
        indexItem = self.model.index(index.row(), 0, index.parent())
        filePath = self.model.filePath(indexItem)
        #print(filePath)
        if 'Force mode' in os.listdir(filePath):
            self.plotForce(os.path.join(filePath,'Force mode'))
            
        #if  re.match(".+[0-9]{8}[A-Z]$",str(filePath))!=None:
        #    self.msg('Loading File Names in {0}'.format(str(filePath)))
        #    self.update_list(str(filePath))



app = QApplication(sys.argv)
window = MainWindow()
window.show()
app.exec()
