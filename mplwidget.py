
from PyQt5.QtWidgets import QWidget,QSizePolicy,QVBoxLayout,\
                        QHBoxLayout,QLabel
                        
from matplotlib.backends.backend_qt5agg import FigureCanvas,\
    NavigationToolbar2QT


from matplotlib.figure import Figure

class MplCanvas(FigureCanvas):
    def __init__(self):
        self.fig = Figure()
        self.ax1 = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.2,bottom=0.15,right=0.8)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                    QSizePolicy.Expanding,
                                    QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

class MplWidget(QWidget):
    def __init__(self, parent = None):
        QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.vbl = QVBoxLayout()
        self.hbl = QHBoxLayout()
        self.ntb = NavigationToolbar2QT(self.canvas,parent)
        self.label = QLabel('k2 viewer v0.1')
        self.hbl.addWidget(self.ntb)
        self.hbl.addWidget(self.label)
        self.vbl.addWidget(self.canvas)
        self.vbl.addLayout(self.hbl)
        self.setLayout(self.vbl)
