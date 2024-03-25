
from PyQt5.QtWidgets import QWidget,QSizePolicy,QVBoxLayout,\
                        QHBoxLayout,QLabel
                        
from matplotlib.backends.backend_qt5agg import FigureCanvas,\
    NavigationToolbar2QT


from matplotlib.figure import Figure
#from matplotlib.pyplot import subplots

class MplCanvas(FigureCanvas):
    def __init__(self,rightax=False):
        self.fig = Figure()
        self.ax1 = self.fig.add_subplot(111)
        if rightax:
            self.bx1 = self.ax1.twinx()
        self.fig.subplots_adjust(left=0.2,bottom=0.15,right=0.8)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                    QSizePolicy.Expanding,
                                    QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

class MplWidget(QWidget):
    def __init__(self, parent = None,rightax=False):
        QWidget.__init__(self, parent)
        self.canvas = MplCanvas(rightax)
        self.vbl = QVBoxLayout()
        self.hbl = QHBoxLayout()
        self.ntb = NavigationToolbar2QT(self.canvas,parent)
        #self.label = QLabel('k2 viewer v0.1')
        self.hbl.addWidget(self.ntb)
        #self.hbl.addWidget(self.label)
        self.vbl.addWidget(self.canvas)
        self.vbl.addLayout(self.hbl)
        self.setLayout(self.vbl)

class MplCanvas4(FigureCanvas):
    def __init__(self):
        self.fig = Figure()
        ll=0.15
        bb=0.15
        h=0.2
        w=0.8
        self.ax1 = self.fig.add_axes((ll,bb+3*h,w,h))
        self.ax2 = self.fig.add_axes((ll,bb+2*h,w,h))
        self.ax3 = self.fig.add_axes((ll,bb+1*h,w,h))
        self.ax4 = self.fig.add_axes((ll,bb+0*h,w,h))
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                    QSizePolicy.Expanding,
                                    QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

class MplWidget4(QWidget):
    def __init__(self, parent = None):
        QWidget.__init__(self, parent)
        self.canvas = MplCanvas4()
        self.vbl = QVBoxLayout()
        self.hbl = QHBoxLayout()
        self.ntb = NavigationToolbar2QT(self.canvas,parent)
        self.label = QLabel('k2 viewer v0.1')
        self.hbl.addWidget(self.ntb)
        self.hbl.addWidget(self.label)
        self.vbl.addWidget(self.canvas)
        self.vbl.addLayout(self.hbl)
        self.setLayout(self.vbl)
