import os
import sys
import logging
import re
import glob
import gzip
import pickle
import numpy as np
import pyqtgraph as pg
from scipy import ndimage as ndi
from collections import defaultdict
from PyQt6 import QtCore, QtWidgets, QtGui
from p3fc.lib.gui import Ui_MainWindow
from p3fc.lib.utility import read_pilatus_cbf, read_pilatus_tif, read_pilatus_tif_gz, get_run_info, pilatus_pad,\
                             convert_frame_APS_Bruker, convert_frame_SP8_Bruker, convert_frame_SP8_Bruker_gz,\
                             convert_frame_DLS_Bruker, write_bruker_frame, bruker_header
# todo
# use tth to calculate beamcenter offset on rotation
# clear patches when loading new folder -> or setup a dict structure to keep them in order!

class FillRectROI(pg.RectROI):
    def __init__(self, pos, size, brush=None, **args):
        pg.RectROI.__init__(self, pos, size, **args)

        if brush is None:
            brush = (255, 255, 255, 75)
        self.setBrush(brush)
    
    def set_handles(self, size=5, width=0):
        for handle in self.getHandles():
            handle.radius = size
            handle.pen.setWidth(width)

    def setBrush(self, *args, **kwargs):
        """
        Set the brush to use when drawing the ROI shape.
        For arguments, see :func:`mkBrush <pyqtgraph.mkBrush>`.
        """
        self.brush = pg.mkBrush(*args, **kwargs)
        self.currentBrush = self.brush
        self.update()
        
    def paint(self, p, opt, widget):
        # Note: don't use self.boundingRect here, because subclasses may need to redefine it.
        r = QtCore.QRectF(0, 0, self.state['size'][0], self.state['size'][1]).normalized()
        
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        p.setPen(self.currentPen)
        p.setBrush(self.currentBrush)
        p.translate(r.left(), r.top())
        p.scale(r.width(), r.height())
        p.drawRect(0, 0, 1, 1)

class FillCircleROI(pg.CircleROI):
    def __init__(self, pos, size=None, radius=None, brush=None, **args):
        if size is None:
            if radius is None:
                raise TypeError("Must provide either size or radius.")
            size = (radius*2, radius*2)
        pg.EllipseROI.__init__(self, pos, size, aspectLocked=True, **args)
        
        if brush is None:
            brush = (255, 255, 255, 75)
        self.setBrush(brush)

    def set_handles(self, size=5, width=0):
        for handle in self.getHandles():
            handle.radius = size
            handle.pen.setWidth(width)
    
    def setBrush(self, *args, **kwargs):
        """
        Set the brush to use when drawing the ROI shape.
        For arguments, see :func:`mkBrush <pyqtgraph.mkBrush>`.
        """
        self.brush = pg.mkBrush(*args, **kwargs)
        self.currentBrush = self.brush
        self.update()

    def paint(self, p, opt, widget):
        r = self.boundingRect()
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        p.setPen(self.currentPen)
        p.setBrush(self.currentBrush)
        p.scale(r.width(), r.height())## workaround for GL bug
        r = QtCore.QRectF(r.x()/r.width(), r.y()/r.height(), 1,1)
        p.drawEllipse(r)

class Main_GUI(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        logging.debug(self.__class__.__name__)
        super(QtWidgets.QMainWindow, self).__init__()
        self.setupUi(self)
        self.status = QtWidgets.QLabel()
        self.status.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.statusBar.addWidget(self.status, 1)
        
        # essentials
        self.init_icons()
        self.init_styles()
        self.init_vars()
        self.set_tooltips()
        
        # window icon
        self.setWindowIcon(self.windowIcon)
        
        # installEventFilters allow us to use the 'line_edit' widget as a button 
        # and catch the FocusIn event to switch between data-path and sfrm-path
        # to be altered by the Filesystem browser
        self.paths_active = self.le_input
        self.le_input.installEventFilter(self)
        self.le_output.installEventFilter(self)
        
        # initiate the filebrowser
        self.init_file_browser()
        
        # init FrameView class
        pg.setConfigOptions(imageAxisOrder='row-major', background='k', leftButtonPan=True)
        self.glwidget.viewport().setAttribute(QtCore.Qt.WidgetAttribute.WA_AcceptTouchEvents, False) 

        self.glwidget.setAspectLocked(True)
        self.plt = self.glwidget.addPlot()
        self.plt.setAspectLocked(True)
        self.plt.hideAxis('bottom')
        self.plt.hideAxis('left')

        self.scatter_bc = pg.ScatterPlotItem(symbol = 'x',
                                             size = 16,
                                             brush = pg.mkBrush((255, 255, 255)),
                                             pen = pg.mkPen((0, 0, 0)))
        self.scatter_bc.setZValue(1000)
        self.plt.addItem(self.scatter_bc)

        # Monkey-patch the image to use our custom hover function.
        # This is generally discouraged (you should subclass ImageItem instead),
        # but it works for a very simple use like this.
        self.img = pg.ImageItem()
        self.plt.addItem(self.img)
        self.colormaps = sorted(pg.colormap.listMaps())
        
        self.cmap = pg.colormap.get(self.colormap, skipCache=True)
        self.img.setLookupTable( self.cmap.getLookupTable(nPts=256) )
        #self.img.setColorMap(self.colormap)
        self.img.setZValue(-2)
        self.img.hoverEvent = self.imageHoverEvent
        self.patches = defaultdict(list)
        self.patches_reset_size()
        
        # link GUI to functions
        self.tb_convert.clicked.connect(self.start_conversion)
        self.cb_link.stateChanged.connect(self.check_path_link)
        self.tb_mask_save.clicked.connect(self.mask_prepare_writing)
        self.hs_mask_int.valueChanged.connect(self.mask_change_frame_max_int)
        self.tb_mask_next_img.clicked.connect(lambda: self.mask_change_image_rel(inc =  1))
        self.tb_mask_prev_img.clicked.connect(lambda: self.mask_change_image_rel(inc = -1))
        self.cb_mask_fname.currentIndexChanged.connect(self.mask_change_image_abs)
        #self.tb_mask_reset.clicked.connect(self.FVObj.reset_patches)
        self.tabWidget.currentChanged.connect(self.on_tab_change)
        
        self.action_add_circle.triggered.connect(self.patches_circs_add)
        self.action_rem_circle.triggered.connect(self.patches_circs_rem)
        self.action_flip_image.triggered.connect(self.change_image)
        self.action_set_wavelength.triggered.connect(self.set_wavelength)
        self.action_set_twotheta.triggered.connect(self.set_twotheta)
        
        # disable the draw-mask tabWidget
        # enable if valid images are loaded
        self.tabWidget.setTabEnabled(1, False)
        self.menu_mask.setEnabled(False)
        
        # hide progress and status-bar on startup
        self.pb_convert.hide()
        self.statusBar.hide()
    
    def set_wavelength(self):
        if not self.action_set_wavelength.isChecked():
            return
        val, ok = QtWidgets.QInputDialog.getDouble(self, 'Set Experimental Parameter', 'Wavelength [\u212B]', value=0.2486, min=0.0, max=1.0, decimals=4, step=0.0001)
        if ok:
            self.exp_wavelength = val

    def set_twotheta(self):
        val, ok = QtWidgets.QInputDialog.getDouble(self, 'Set Experimental Parameter', '2-Theta offset [%]', value=self.SP8_tth_corr*100, min=-10.0, max=10.0, decimals=1, step=0.1)
        if ok:
            self.SP8_tth_corr = round(val / 100, 3)
            self.change_image()

    def set_tooltips(self):
        logging.debug(self.__class__.__name__)
        # add tooltips
        self.tb_convert.setToolTip('Start the conversion')
        self.le_output.setToolTip('Current output directory.\nIf unlinked: Select to specify the target output directory using the file-browser.\nManual editing is allowed, non-existing paths will be created recursively.')
        self.le_input.setToolTip('Current input directory.\nIf unlinked: Select to specify the target input directory using the file-browser.')
        self.cb_link.setToolTip('Link/Unlink output directory and input directory.\nIf linked: the output directory follows the input directory (plus added suffix).\nIf unlinked: Select either to specify the target directory using the filebrowser.')
        self.cb_overwrite.setToolTip('Overwrite existing files in the output directory?')
        
        self.action_set_wavelength.setToolTip('Check and manually set the wavelength, uncheck to use the .inf information.')
        self.action_set_twotheta.setToolTip('Check and manually set an 2-Theta offset, uncheck to use the .inf information.')

        self.action_add_circle.setToolTip('Add a pair of circles. Use the green circle to unmask regions.')
        self.action_rem_circle.setToolTip('Remove the last Circle pair.')
        self.action_write_bruker_sfrm.setToolTip('Write a bruker .sfrm file?')
        self.action_write_numpy_npy.setToolTip('Write numpy .npy file?')
        self.action_show_matplotlib.setToolTip('Check to plot and show the final mask using matplotlib.')
        self.action_use_padding.setToolTip('Check to pad the mask to a multiple of 8 (SAINT).')
    
    def init_file_browser(self):
        logging.debug(self.__class__.__name__)
        # use the QFileSystemModel
        self.model = QtGui.QFileSystemModel()
        self.model.setReadOnly(True)
        self.model.setRootPath('')
        # currently only shows directories
        # use:  | QtCore.QDir.AllEntries
        # to show files
        self.model.setFilter(QtCore.QDir.Filter.AllDirs | QtCore.QDir.Filter.NoDotAndDotDot)# | QtCore.QDir.Filter.AllEntries)
        
        # set treeView to use the QFileSystemModel
        self.treeView.setAnimated(False)
        self.treeView.setUniformRowHeights(True)
        self.treeView.setModel(self.model)
        self.treeView.sortByColumn(0, QtCore.Qt.SortOrder.AscendingOrder)
        # don't strech last column
        self.treeView.header().setStretchLastSection(True)
        # stretch first
        self.treeView.header().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        # hide 'size' and 'type' header columns
        self.treeView.setColumnHidden(1, True)
        self.treeView.setColumnHidden(2, True)
        self.treeView.setColumnHidden(3, True)
        
        # scroll the treeview to start_dir
        # apply start_dir as current dir by calling the 'on_click' function
        index = self.model.index(os.getcwd())
        self.on_treeView_clicked(index)
        self.treeView.setCurrentIndex(index)
        self.treeView.setExpanded(index, True)
        # delay the scrollTo call to allow the treeview
        # to expand all sections leading to the directory
        QtCore.QTimer.singleShot(500, self.current_row_changed)
        
    def current_row_changed(self):
        index = self.model.index(os.getcwd())
        self.treeView.scrollTo(index, QtWidgets.QAbstractItemView.ScrollHint.PositionAtCenter)
        self.treeView.setFocus()

    def init_icons(self):
        logging.debug(self.__class__.__name__)
        # icons
        self.windowIcon = self.style().standardIcon(getattr(QtWidgets.QStyle.StandardPixmap, 'SP_BrowserReload'))
        self.ic_MessageBoxWarning = self.style().standardIcon(getattr(QtWidgets.QStyle.StandardPixmap, 'SP_MessageBoxWarning'))
        self.ic_MessageBoxCritical = self.style().standardIcon(getattr(QtWidgets.QStyle.StandardPixmap, 'SP_MessageBoxCritical'))
        self.ic_MessageBoxInformation = self.style().standardIcon(getattr(QtWidgets.QStyle.StandardPixmap, 'SP_MessageBoxInformation'))
        self.ic_MessageBoxQuestion = self.style().standardIcon(getattr(QtWidgets.QStyle.StandardPixmap, 'SP_MessageBoxQuestion'))
        
    def init_styles(self):
        logging.debug(self.__class__.__name__)
        self.pb_style = ('QProgressBar        {text-align: center; border: 1px solid grey; border-radius: 2px}'
                         'QProgressBar:chunk  {background: qlineargradient(x1: 0, y1: 0.5, x2: 1, y2: 0.5, stop: 0 rgb(  0, 171, 164), stop: 1 rgb( 55, 160, 203));}')
        
        #self.tb_style = ('QToolButton          {background-color: rgb(240, 240, 240); color: rgb(  0,   0,   0); border: 1px solid rgb( 75,  75,  75); border-radius: 2px}'
        #                 'QToolButton:hover    {background-color: rgb(255, 255, 255); color: rgb(  0,   0,   0); border: 1px solid rgb( 75,  75,  75)}'
        #                 'QToolButton:pressed  {background-color: rgb(255, 255, 255); color: rgb(  0,   0,   0); border: 1px solid rgb( 75,  75,  75)}'
        #                 'QToolButton:checked  {background-color: rgb(200, 200, 200); color: rgb(  0,   0,   0); border: 1px solid rgb( 75,  75,  75)}'
        #                 'QToolButton:disabled {background-color: rgb(220, 200, 200); color: rgb(  0,   0,   0); border: 1px solid rgb( 75,  75,  75)}')
                                 
        self.le_style_coupled = self.le_input.styleSheet()#('QLineEdit      {background-color: rgb(240, 240, 240)}')
        
        self.le_style_single = ('QLineEdit       {border-width: 2px; border-style: solid; border-color: rgb(100, 255, 100)}'
                                'QLineEdit:hover {border-width: 2px; border-style: solid; border-color: rgb(100, 255, 100)}'
                                'QLineEdit:focus {border-width: 2px; border-style: solid; border-color: rgb(100, 255, 100)}')
        
        # apply style sheets
        self.le_input.setStyleSheet(self.le_style_coupled)
        self.le_output.setStyleSheet(self.le_style_coupled)
        #self.tb_convert.setStyleSheet(self.tb_style)
        #self.tb_mask_save.setStyleSheet(self.tb_style)
        #self.tb_mask_reset.setStyleSheet(self.tb_style)
        self.pb_convert.setStyleSheet(self.pb_style)
        #self.tabWidget.setStyleSheet('QTabBar::tab { height: 25px; width: 200px; }')
        #self.tb_mask_prev_img.setStyleSheet(self.tb_style)
        #self.tb_mask_next_img.setStyleSheet(self.tb_style)
        self.tb_mask_prev_img.setArrowType(QtCore.Qt.ArrowType.LeftArrow)
        self.tb_mask_next_img.setArrowType(QtCore.Qt.ArrowType.RightArrow)

    def init_vars(self):
        logging.debug(self.__class__.__name__)
        '''
         
        '''
        self.exp_pixelsize = 172e-6  # m
        self.exp_beamcenter_x = None # m
        self.exp_beamcenter_y = None # m CCD_SPATIAL_BEAM_POSITION=501.05 528.57;
        self.exp_distance = None     # m SCAN_DET_RELZERO=2.000 0.000   130.00;
        self.exp_wavelength = None   # Ang SCAN_WAVELENGTH=0.2482;
        self.SP8_tth_corr = 0.0      # Correction factor for 2-theta offset [%]
        self.current_tth = 0.0       # Indicator to change the patches to new positions
        self.reset_patches = True
        self.mask_negative = True
        self.flag_reset_view = False
        self.fRnum = None
        self.fStem = None
        self.runList = []
        self.framesList = []
        self.patches_base = []
        self.patches_circs = []
        self.currentIndex = 0
        self.currentFrame = None
        self.path_mask = None
        self.path_inf = None
        self.suffix = '_sfrm'
        self.colormap = 'viridis'
        #self.colormap = 'CET-C5s'
        self.patch_size_default = 100
        self.patch_size_increment = 50
        self.handle_size = 12
        self.handle_width = 5

        _cmap = pg.colormap.get(self.colormap)
        _color_05 = _cmap.map(0.5, mode='qcolor')
        _color_03 = _cmap.map(0.3, mode='qcolor')
        _color_00 = _cmap.map(0.0, mode='qcolor')
        #self.pen_unmask = pg.mkPen((51, 255, 153, 255), width=3)
        self.pen_unmask = pg.mkPen(_color_00, width=3)
        self.pen_mask = pg.mkPen(_color_00, width=3)
        self.patch_parameter = {'scaleSnap':False,'translateSnap':False,'movable':True,'rotateSnap':False,
                                'hoverPen':pg.mkPen(_color_03, width=5),
                                'handlePen':pg.mkPen(_color_05),
                                'handleHoverPen':pg.mkPen(_color_03, width=8)}
        #_color_05.setAlphaF(0.5)
        #_color_03.setAlphaF(0.5)
        #_color_00.setAlphaF(0.5)
        self.brush_mask = pg.mkBrush((255, 255, 255, 100))
        self.brush_unmask = pg.mkBrush((83, 255, 69, 100))
        self.brush_saved = pg.mkBrush((0, 255, 231, 100))
        self.brush_unsaved = pg.mkBrush((255, 255, 255, 100))
        
        # some hardcoded limits that might make sense
        self.hs_mask_int.setMinimum(1)
        self.hs_mask_int.setMaximum(100)
        
        #########################################
        ##  Add new format identifiers here!   ##
        #########################################
        self.exts = ('*_*.tif', '*_*.cbf', '*_*.tif.gz')
        self.availableFormats = [self.format_SP8,
                                 self.format_SP8_gz,
                                 self.format_APS,
                                 self.format_DLS]
    
    ##############################################
    ##         Frame Format definitions         ##
    ##############################################
    def check_format(self):
        logging.debug(self.__class__.__name__)
        for aCheckFunc in self.availableFormats:
            if aCheckFunc():
                return True
            else:
                continue
        return False
    
    def format_DLS(self):
        logging.debug(self.__class__.__name__)
        '''
        Check the first file if reformatting to Bruker name format is possible
        any_name_#run_#frame.tif -> any_name_rr_ffff.sfrm
        '''
        try:
            fhead, fname = os.path.split(self.currentFrame)
            bname, ext = os.path.splitext(fname)
            if not ext == '.cbf':
                return False
            # open file and check: _diffrn.id DLS_I19-1
            with open(self.currentFrame, 'rb') as oFrame:
                try:
                    id = re.search(rb'_diffrn.id\s+(?P<id>.+)', oFrame.read(2048)).group('id').decode().strip()
                except AttributeError:
                    return False
            if not id == 'DLS_I19-1':
                return False
            fstm, rnum, fnum, flen = get_run_info(bname)
            self.fRnum = rnum                         # Run number
            self.fStem = fstm                         # Frame name up to the run number
            self.fStar = '{:>0{w}}.'.format(1, w=flen)# Number indicating start of a run
            self.fInfo = (1679, 1475, 0, np.int32)    # Frame info (rows, cols, offset)
            self.fSite = 'DLS'                        # Facility identifier
            self.fFunc = read_pilatus_cbf             # Frame read function (from _Utility)
            self.fRota = False                        # rotate the frame upon conversion?
            self.detector_type = 'PILATUS'            # detector type for SAINT
            return True
        except (ValueError, IndexError):
            return False
    
    def format_APS(self):
        logging.debug(self.__class__.__name__)
        '''
        Check the first file if reformatting to Bruker name format is possible
        any_name_#run_#frame.tif -> any_name_rr_ffff.sfrm
        '''
        try:
            fhead, fname = os.path.split(self.currentFrame)
            bname, ext = os.path.splitext(fname)
            if not ext == '.tif':
                return False
            # open file and check S/N: 10-0147
            with open(self.currentFrame, 'rb') as oFrame:
                SN = re.search(rb'S/N\s+(?P<SN>\d+\-\d+)', oFrame.read(128)).group('SN').decode()
            if not SN == '10-0147':
                return False
            fstm, rnum, fnum, flen = get_run_info(bname)
            self.fRnum = rnum                         # Run number
            self.fStem = fstm                         # Frame name up to the run number
            self.fStar = '{:>0{w}}.'.format(1, w=flen)# Number indicating start of a run
            self.fInfo = (1043, 981, 4096, np.int32)  # Frame info (rows, cols, offset)
            self.fSite = 'APS'                        # Facility identifier
            self.fFunc = read_pilatus_tif             # Frame read function (from _Utility)
            self.fRota = True                         # rotate the frame upon conversion?
            self.detector_type = 'PILATUS'            # detector type for SAINT
            return True
        except (ValueError, IndexError):
            return False
    
    def format_SP8(self):
        logging.debug(self.__class__.__name__)
        '''
        Check the first file if name is compatible with SPring-8 convention
        e.g. any_name_rrfff.tif, where rr is the 2 digit run numer: 00 - 99
        fff is the 3 digit frame number: 001 - 999
        '''
        try:
            fhead, fname = os.path.split(self.currentFrame)
            bname, ext = os.path.splitext(fname)
            if not ext == '.tif':
                return False
            # open file and check S/N: 10-0163
            with open(self.currentFrame, 'rb') as oFrame:
                SN = re.search(rb'S/N\s+(?P<SN>\d+\-\d+)', oFrame.read(128)).group('SN').decode()
            if not SN == '10-0163':
                return False
            fstm, rnum, fnum, flen = get_run_info(bname)
            ########################################
            ## USE FNUM TO DEFINE START OF RUN!!! ##
            ########################################
            self.fRnum = rnum                         # Run number
            self.fStem = fstm                         # Frame name up to the run number
            self.fStar = '{:>0{w}}.'.format(1, w=flen)# Number indicating start of a run
            self.fInfo = (1043, 981, 4096, np.int32)  # Frame info (rows, cols, offset)
            self.fSite = 'SP8'                        # Facility identifier
            self.fFunc = read_pilatus_tif             # Frame read function (from _Utility)
            self.fRota = True                         # rotate the frame upon conversion?
            self.detector_type = 'PILATUS'            # detector type for SAINT
            return True
        except ValueError:
            return False
    
    def format_SP8_gz(self):
        logging.info(self.__class__.__name__)
        '''
        Check the first file if name is compatible with SPring-8 convention
        e.g. any_name_rrfff.tif, where rr is the 2 digit run numer: 00 - 99
        fff is the 3 digit frame number: 001 - 999
        '''
        try:
            fhead, fname = os.path.split(self.currentFrame)
            zname, ext = os.path.splitext(fname)
            if not ext == '.gz':
                return False
            # open file and check S/N: 10-0163
            with gzip.open(self.currentFrame, 'rb') as oFrame:
                SN = re.search(rb'S/N\s+(?P<SN>\d+\-\d+)', oFrame.read(128)).group('SN').decode()
            if not SN == '10-0163':
                return False
            bname, ext = os.path.splitext(zname)
            fstm, rnum, fnum, flen = get_run_info(bname)
            ########################################
            ## USE FNUM TO DEFINE START OF RUN!!! ##
            ########################################
            self.fRnum = rnum                         # Run number
            self.fStem = fstm                         # Frame name up to the run number
            self.fStar = '{:>0{w}}.'.format(1, w=flen)# Number indicating start of a run
            self.fInfo = (1043, 981, 4096, np.int32)  # Frame info (rows, cols, offset)
            self.fSite = 'SP8_gz'                     # Facility identifier
            self.fFunc = read_pilatus_tif_gz          # Frame read function (from _Utility)
            self.fRota = True                         # rotate the frame upon conversion?
            self.detector_type = 'PILATUS'            # detector type for SAINT
            return True
        except ValueError:
            return False
    ##############################################
    ##       END Frame Format definitions       ##
    ##############################################
    
    def read_inf(self):
        # check if info file exists
        if os.path.isfile(self.path_inf):
            # extract header information
            with open(self.path_inf) as rFile:
                infoFile = rFile.read()
                self.exp_beamcenter_y, self.exp_beamcenter_x = [float(i) for i in re.search(r'CCD_SPATIAL_BEAM_POSITION\s*=\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*;', infoFile).groups()]
                self.exp_wavelength = float(re.search(r'SCAN_WAVELENGTH\s*=\s*(\d+\.\d+)\s*;', infoFile).groups()[0])
                self.exp_tth, inf_distance = [float(i) for i in re.search(r'SCAN_DET_RELZERO\s*=\s*-*\d+\.\d+\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*;', infoFile).groups()]
                self.exp_distance = inf_distance * 1e-3
                # apply SP8 2-theta correction factor
                self.exp_tth += self.exp_tth * self.SP8_tth_corr
                if self.exp_tth == self.current_tth:
                    self.reset_patches = False
                else:
                    self.current_tth = self.exp_tth
                    self.reset_patches = True
                # Add POBI 2-theta offset
                offset_tth = np.tan(np.deg2rad(self.exp_tth)) * self.exp_distance / self.exp_pixelsize
                self.exp_beamcenter_x += offset_tth
        else:
            logging.warning(f'WARNING: Info file {os.path.basename(self.path_inf)} is missing')
            self.exp_beamcenter_x = self.img_dim_x/2
            self.exp_beamcenter_y = self.img_dim_y/2

    def mask_prepare_writing(self):
        logging.debug(self.__class__.__name__)
        oPath = os.path.abspath(self.le_output.text())
        self.create_output_directory(oPath)
        self.mask_write()
    
    def mask_change_image_abs(self, idx):
        logging.debug(self.__class__.__name__)
        self.currentFrame = os.path.abspath(self.runList[idx])
        self.check_format()
        self.change_image()
    
    def mask_change_image_rel(self, inc):
        logging.debug(self.__class__.__name__)
        check_idx = self.currentIndex + int(inc)
        if check_idx < 0 or check_idx >= len(self.runList):
            return
        self.currentIndex += int(inc)
        # setCurrentIndex calls self.mask_change_image_abs
        self.cb_mask_fname.setCurrentIndex(self.currentIndex)
    
    def mask_change_frame_max_int(self):
        #logging.debug(self.__class__.__name__)
        self.img.setLevels(levels=[-2, self.hs_mask_int.value()])
        
    def eventFilter(self, obj, event):
        #logging.debug(self.__class__.__name__)
        '''
         
        '''
        if event.type() == QtCore.QEvent.Type.FocusIn:
            if (obj == self.le_input or obj == self.le_output) and self.cb_link.isChecked():
                self.le_input.setStyleSheet(self.le_style_coupled)
                self.le_output.setStyleSheet(self.le_style_coupled)
                self.paths_active = self.le_input
                self.treeView.scrollTo(self.model.index(self.le_input.text()), QtWidgets.QAbstractItemView.ScrollHint.PositionAtCenter)
            elif obj == self.le_output and not self.cb_link.isChecked():
                self.le_input.setStyleSheet(self.le_style_coupled)
                self.le_output.setStyleSheet(self.le_style_single)
                self.paths_active = self.le_output
                self.treeView.scrollTo(self.model.index(self.le_output.text()), QtWidgets.QAbstractItemView.ScrollHint.PositionAtCenter)
            elif obj == self.le_input and not self.cb_link.isChecked():
                self.le_input.setStyleSheet(self.le_style_single)
                self.le_output.setStyleSheet(self.le_style_coupled)
                self.paths_active = self.le_input
                self.treeView.scrollTo(self.model.index(self.le_input.text()), QtWidgets.QAbstractItemView.ScrollHint.PositionAtCenter)
        return super(Main_GUI, self).eventFilter(obj, event)

    def popup_window(self, _title, _text, _info):
        logging.debug(self.__class__.__name__)
        '''
         _icon:
            QtWidgets.QMessageBox.Icon.NoIcon      0 the message box does not have any icon.
            QtWidgets.QMessageBox.Icon.Information 1 an icon indicating that the message is nothing out of the ordinary.
            QtWidgets.QMessageBox.Icon.Warning     2 an icon indicating that the message is a warning, but can be dealt with.
            QtWidgets.QMessageBox.Icon.Critical    3 an icon indicating that the message represents a critical problem.
            QtWidgets.QMessageBox.Icon.Question    4 an icon indicating that the message is asking a question.
        '''
        if _title.upper() == 'INFORMATION':
            _wicon = self.windowIcon#self.ic_MessageBoxInformation
            _icon = QtWidgets.QMessageBox.Icon.Information
        elif _title.upper() == 'WARNING':
            _wicon = self.windowIcon#self.ic_MessageBoxWarning
            _icon = QtWidgets.QMessageBox.Icon.Warning
        elif _title.upper() == 'CRITICAL':
            _wicon = self.windowIcon#self.ic_MessageBoxCritical
            _icon = QtWidgets.QMessageBox.Icon.Critical
        elif _title.upper() == 'QUESTION':
            _wicon = self.windowIcon#self.ic_MessageBoxQuestion
            _icon = QtWidgets.QMessageBox.Icon.Question
        else:
            _wicon = self.windowIcon
            _icon = QtWidgets.QMessageBox.Icon.NoIcon
        msgBox = QtWidgets.QMessageBox()
        msgBox.setWindowIcon(_wicon)
        msgBox.setIcon(_icon)
        msgBox.setWindowTitle(_title)
        msgBox.setText(_text)
        msgBox.setInformativeText(_info)
        msgBox.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
        msgBox.exec()
    
    def check_path_link(self):
        logging.debug(self.__class__.__name__)
        '''
         switch between input / output path manipulation
         this function ONLY takes care about the stylesheet
         and determines which one is active after toggling
         the link checkbox 'self.cb_link.isChecked()'
          - if checked, both paths reference to 'self.le_input'
          - 'self.paths_active' links to the currently active object
          - checked: 'self.le_input' is active (always)
          - unchecked: 'self.le_output' is activated by default,
            can be switched
        '''
        self.le_input.setStyleSheet(self.le_style_coupled)
        self.le_output.setStyleSheet(self.le_style_coupled)
        if self.cb_link.isChecked():
            self.paths_active = self.le_input
            self.le_output.setText(self.le_input.text() + self.suffix)
        else:
            self.paths_active = self.le_output
            self.paths_active.setStyleSheet(self.le_style_single)
    
    def reset_view(self):
        self.plt.setXRange(0, self.img_dim_x, padding=0)
        self.plt.setYRange(0, self.img_dim_y, padding=0)

    def change_image(self):
        _, data = self.fFunc(self.currentFrame, *self.fInfo)
        if self.fRota:
            data = np.rot90(data, k=1, axes=(1, 0))
        if self.action_flip_image.isChecked():
            data = np.flipud(data)
        self.img.setImage(data, rotate=self.fRota)
        self.img_dim_y, self.img_dim_x = data.shape
        if self.flag_reset_view:
            self.reset_view()
        self.mask_change_frame_max_int()
        self.add_resolution_label()
        
        #iPath = os.path.abspath(self.le_input.text())
        oPath = os.path.abspath(self.le_output.text())
        #self.path_inf = os.path.join(iPath, '{}_{:>02}_{}inf'.format(self.fStem, int(self.fRnum), self.fStar))
        self.path_inf = f'{os.path.splitext(os.path.splitext(self.currentFrame)[0])[0]}.inf'
        self.path_mask = os.path.join(oPath, '{}_xa_{:>02}_0001.sfrm'.format(self.fStem, int(self.fRnum)))
        self.path_patches = os.path.join(oPath, '{}_xa_{:>02}_0001.msk'.format(self.fStem, int(self.fRnum)))
        self.read_inf()
        self.patches_clear()
        self.patches_reset_size()
        self.patches_load()
        self.patches_add()
        self.add_beamcenter()

    def add_beamcenter(self):
        if self.exp_beamcenter_x is None:
            return
        if self.exp_beamcenter_y is None:
            return
        
        self.scatter_bc.setData([self.exp_beamcenter_x], [self.exp_beamcenter_y])

    def on_tab_change(self, idx):
        '''
         update frame only if Frameviewer tab is opened
        '''
        if idx == 1:
            self.mask_change_image_abs(self.currentIndex)
            self.menu_mask.setEnabled(True)
            self.glwidget.setFocus()
        else:
            self.menu_mask.setEnabled(False)
            return
    
    def mask_add_obj(self, obj, val):
        '''
        Circles and Ellipses
        Note: returnMappedCoords is not yet supported for this ROI type.
        
        Workaround taken from:
        https://groups.google.com/g/pyqtgraph/c/fcysRvIcJi8
        https://groups.google.com/g/pyqtgraph/c/-kNPXxDeERs
        
        Still produces erroneously unmasked regions so we are not going to use ellipses
        The tilting of the rectangle should be only minute, showing only single unmasked pixels
        Application of scipy.ndimage.binary_erosion() before writing the mask should make it smooth and clean.
        - This is bad!
        - Still not supported in pyqtgraph version 0.13.3
        '''
        cols, rows = self.img.image.shape
        m = np.mgrid[:cols,:rows]
        possx = m[0,:,:]
        possy = m[1,:,:]
        possx.shape = cols, rows
        possy.shape = cols, rows
        mpossx = obj.getArrayRegion(possx, self.img).astype(int)
        mpossy = obj.getArrayRegion(possy, self.img).astype(int)
        self.msk[mpossx, mpossy] = val

    def mask_write(self):
        ##############################
        # DEMANDING FOR MANY PATCHES #
        # !!!!!!! TWEAK HERE !!!!!!! #
        #      -> mask_add_obj()     #
        ##############################
        self.msk = np.ones(self.img.image.shape)
        self.patches['circles'] = []
        newlist = sorted(self.patches_circs, key=lambda x: x[0].size().manhattanLength(), reverse=True)
        for obj, val in newlist:
            self.mask_add_obj(obj, val)
            self.patches['circles'].append(['circ', obj.pos(), obj.size(), obj.angle(), val])
        ##########################

        self.patches['base'] = []
        for obj, val in self.patches_base:
            self.mask_add_obj(obj, val)
            if isinstance(obj, pg.graphicsItems.ROI.RectROI):
                self.patches['base'].append(['rect', obj.pos(), obj.size(), obj.angle(), val])
            elif isinstance(obj, pg.graphicsItems.ROI.CircleROI):
                self.patches['base'].append(['circ', obj.pos(), obj.size(), obj.angle(), val])
        # interpolation fails -> erode the mask
        self.msk = ndi.binary_erosion(self.msk)
        
        # mask negatives?
        if self.mask_negative:
            self.msk[self.img.image < 0] = 0

        # get the frame saint ready
        # - pad with zeros
        if self.action_use_padding.isChecked():
            self.msk, offset_rows, offset_cols = pilatus_pad(self.msk, fill=0)

        header = bruker_header()
        # fill known header entries
        header['NCOLS']       = [self.msk.shape[1]]                 # Number of pixels per row; number of mosaic tiles in X; dZ/dX
        header['NROWS']       = [self.msk.shape[0]]                 # Number of rows in frame; number of mosaic tiles in Y; dZ/dY value
        #header['CCDPARM'][:] = [1.47398, 36.60, 359.8295, 0.0, 163810.0] # readnoise, electronsperadu, electronsperphoton, bruker_bias, bruker_fullscale
        #header['DETTYPE'][:] = ['CMOS-PHOTONII', 37.037037, 1.004, 0, 0.425, 0.035, 1]
        header['DETTYPE'][:]  = [self.detector_type, 10.0, 1.0, 0, 0.0, 0.0, 1] # dettype pix512percm cmtogrid circular brassspacing windowthickness accuratetime
        #header['SITE']       = ['Aarhus Huber Diffractometer']           # Site name
        #header['MODEL']      = ['Microfocus X-ray Source']               # Diffractometer model
        #header['TARGET']     = ['Ag Ka']                                 # X-ray target material)
        #header['SOURCEK']    = [50.0]                                    # X-ray source kV
        #header['SOURCEM']    = [0.880]                                   # Source milliamps
        #header['WAVELEN'][:] = [0.560860, 0.559420, 0.563810]            # Wavelengths (average, a1, a2)
        if self.action_set_wavelength.isChecked():
            header['WAVELEN'][:] = [self.exp_wavelength, self.exp_wavelength, self.exp_wavelength] # Wavelengths (average, a1, a2)
        else:
            header['WAVELEN'][:] = [1.0, 1.0, 1.0]                            # Wavelengths (average, a1, a2)
        #header['CORRECT']    = ['INTERNAL, s/n: A110247']                # Flood correction filename
        #header['DARK']       = ['INTERNAL, s/n: A110247']                # Dark current frame name
        #header['WARPFIL']    = ['LINEAR']                                # Spatial correction filename
        #header['LINEAR'][:]  = [1.00, 0.00]                              # bruker_linearscale, bruker_linearoffset
        #header['PHD'][:]     = [0.68, 0.051]                             # Phosphor efficiency, phosphor thickness
        #header['OCTMASK'][:] = [0, 0, 0, 767, 767, 1791, 1023, 1023]

        # write the frame
        if self.action_write_bruker_sfrm.isChecked():
            write_bruker_frame(self.path_mask, header, np.flipud(self.msk))
        
        # save mask as numpy npy file
        if self.action_write_numpy_npy.isChecked():
            np.save(os.path.splitext(self.path_mask)[0], np.flipud(self.msk))
        
        # dump patches dict
        self.patches_save()

        self.cb_mask_stored.setChecked(True)

        # reload image to set saved mask colors
        self.change_image()

        # DEBUG
        # show mask in matplotlib
        if self.action_show_matplotlib.isChecked():
            import matplotlib.pyplot as plt
            plt.axis('off')
            plt.subplots_adjust(0,0,1,1,0,0)
            plt.imshow(np.flipud(self.msk))
            plt.show()

    def patches_reset_size(self):
        self.patch_size_current = self.patch_size_default
    
    def patches_load(self):
        if os.path.exists(self.path_patches):
            self.cb_mask_stored.setChecked(True)
            with open(self.path_patches, 'rb') as rf:
                self.patches = pickle.load(rf)
        else:
            self.cb_mask_stored.setChecked(False)
            if self.reset_patches or len(self.patches['base']) == 0:
                x = self.exp_beamcenter_x
                y = self.exp_beamcenter_y
                self.patches['base'] = []
                self.patches['base'].append(['rect', (x-10, -10), (20, y), 0.0, 0])
                self.patches['base'].append(['circ', (x-20, y-20), (40, 40), 0.0, 0])
    
    def patches_add(self):
        # set brush for saved masks
        if self.cb_mask_stored.isChecked():
            active_brush = self.brush_saved
        else:
            active_brush = self.brush_unsaved
        # add patches
        if 'base' in self.patches:
            for name, pos, size, angle, msk in self.patches['base']:
                if name == 'rect':
                    r_roi = FillRectROI(pos=pos, size=size, angle=angle, pen=self.pen_mask, brush=active_brush, sideScalers=True, **self.patch_parameter)
                    r_roi.addRotateHandle((0.0,1.0), center=(0.5,0.0))
                    r_roi.addRotateHandle((1.0,0.0), center=(0.5,1.0))
                    r_roi.addScaleHandle((0.5,0.0), center=(0.5,1.0))
                    r_roi.addScaleHandle((0.0,0.5), center=(1.0,0.5))
                    r_roi.addScaleHandle((0.0,0.0), center=(1.0,1.0))
                    r_roi.setZValue(100)
                    r_roi.set_handles(size=self.handle_size, width=self.handle_width)
                    self.plt.addItem(r_roi)
                    self.patches_base.append((r_roi, msk))
                elif name == 'circ':
                    c_roi = FillCircleROI(pos=pos, size=size, angle=angle, pen=self.pen_mask, brush=active_brush, **self.patch_parameter)
                    c_roi.setZValue(101)
                    c_roi.set_handles(size=self.handle_size, width=self.handle_width)
                    c_roi.sigRegionChangeFinished.connect(self.patches_adjust_size)
                    self.plt.addItem(c_roi)
                    self.patches_base.append((c_roi, msk))
        if 'circles' in self.patches:
            for idx, (name, pos, size, angle, msk) in enumerate(self.patches['circles']):
                if msk:
                    c_roi = FillCircleROI(pos=pos, size=size, angle=angle, pen=self.pen_unmask, brush=self.brush_unmask, **self.patch_parameter)
                else:
                    c_roi = FillCircleROI(pos=pos, size=size, angle=angle, pen=self.pen_mask, brush=self.brush_mask, **self.patch_parameter)
                c_roi.set_handles(size=self.handle_size, width=self.handle_width)
                c_roi.setZValue(idx)
                c_roi.sigRegionChangeFinished.connect(self.patches_circs_sort)
                self.plt.addItem(c_roi)
                self.patches_circs.append((c_roi, msk))
    
    def patches_adjust_size(self):
        for patch, mask in self.patches_base:
            if isinstance(patch, pg.graphicsItems.ROI.CircleROI):
                if patch.size()[0] > self.patch_size_current:
                    self.patch_size_current = patch.size()[0]

    def patches_save(self):
        with open(self.path_patches, 'wb') as wf:
            pickle.dump(self.patches, wf)
    
    def patches_circs_add(self):
        #x = self.img_dim_x/2 - self.patch_size_current/2 - self.patch_size_increment/2
        #y = self.img_dim_y/2 - self.patch_size_current/2 - self.patch_size_increment/2
        for patch, mask in self.patches_base:
            if isinstance(patch, pg.graphicsItems.ROI.CircleROI):
                x = patch.pos().x() + patch.size()[0]/2 - self.patch_size_current/2 - self.patch_size_increment/2
                y = patch.pos().y() + patch.size()[0]/2 - self.patch_size_current/2 - self.patch_size_increment/2
        self.patch_size_current += self.patch_size_increment
        patch_add = FillCircleROI((x,y), (self.patch_size_current,self.patch_size_current), pen=self.pen_unmask, brush=self.brush_unmask, **self.patch_parameter)
        patch_add.set_handles(size=self.handle_size, width=self.handle_width)
        patch_add.sigRegionChangeFinished.connect(self.patches_circs_sort)
        self.plt.addItem(patch_add)
        self.patches_circs.append((patch_add, 1))
        
        x = x - self.patch_size_increment/2
        y = y - self.patch_size_increment/2
        self.patch_size_current += self.patch_size_increment
        patch_sub = FillCircleROI((x,y), (self.patch_size_current,self.patch_size_current), pen=self.pen_mask, brush=self.brush_mask, **self.patch_parameter)
        patch_sub.set_handles(size=self.handle_size, width=self.handle_width)
        patch_sub.sigRegionChangeFinished.connect(self.patches_circs_sort)
        self.plt.addItem(patch_sub)
        self.patches_circs.append((patch_sub, 0))
        
        self.patches_circs_sort()
    
    def patches_circs_rem(self):
        if self.patches_circs:
            p,_ = self.patches_circs.pop()
            self.plt.removeItem(p)
            p,_ = self.patches_circs.pop()
            self.plt.removeItem(p)
            self.patch_size_current -= 2 * self.patch_size_increment
            self.patches_circs_sort()
    
    def patches_circs_sort(self):
        if self.patches_circs:
            newlist = sorted(self.patches_circs, key=lambda x: x[0].size().manhattanLength(), reverse=True)
            for idx, (obj,_) in enumerate(newlist):
                obj.setZValue(idx)
                size = obj.size().manhattanLength()/2
                if size >= self.patch_size_current:
                    self.patch_size_current = size
        else:
            self.patch_size_current = self.patch_size_default
    
    def patches_clear(self):
        if self.reset_patches:
            self.patches['base'] = []
            self.patches['circles'] = []
        while len(self.patches_base) > 0:
            p,_ = self.patches_base.pop()
            self.plt.removeItem(p)
        while len(self.patches_circs) > 0:
            p,_ = self.patches_circs.pop()
            self.plt.removeItem(p)
    
    def keyPressEvent(self, event):
        k = event.key()

        if not self.tabWidget.currentIndex() == 1:
            return

        # bind 'c' to cycle colormaps
        if k == QtCore.Qt.Key.Key_C:
            if event.modifiers() == QtCore.Qt.KeyboardModifier.ShiftModifier:
                inc = -1
            else:
                inc = +1
            idx = self.colormaps.index(self.colormap) + inc
            if idx >= len(self.colormaps):
                idx = 0
            elif idx < 0:
                idx = len(self.colormaps) - 1
            self.colormap = self.colormaps[idx]
            
            self.cmap = pg.colormap.get(self.colormap, skipCache=True)
            self.img.setLookupTable( self.cmap.getLookupTable(nPts=256) )
            logging.info(self.colormap)
            #self.img.setColorMap(self.colormaps[idx])
        elif k == QtCore.Qt.Key.Key_A:
            if event.modifiers() == QtCore.Qt.KeyboardModifier.ShiftModifier:
                self.patches_circs_rem()
            else:
                self.patches_circs_add()
        elif k == QtCore.Qt.Key.Key_N or k == QtCore.Qt.Key.Key_Right:
            self.mask_change_image_rel(inc = 1)
        elif k == QtCore.Qt.Key.Key_P or k == QtCore.Qt.Key.Key_Left:
            self.mask_change_image_rel(inc = -1)
        elif k == QtCore.Qt.Key.Key_R:
            self.reset_view()

    def add_resolution_label(self):
        font = QtGui.QFont()
        font.setPixelSize(12)
        self.res_label = pg.TextItem(anchor=(1.0,1.0), color=(255,255,255,255), fill=(0,0,0,255))
        self.res_label.setText('')
        self.res_label.setFont(font)
        self.plt.addItem(self.res_label)
        self.res_label.hide()
        self.res_label.setToolTip('d [\u212B]')
        self.res_label.setPos(self.img_dim_x, 0)
    
    def imageHoverEvent(self, event):
        """Hover event linked to cormap
        and should only be active while
        either or both maps are displayed """
        if event.isExit():
            self.res_label.hide()
            return
        
        if self.exp_wavelength is None:
            return
        if self.exp_distance is None:
            return
        if self.exp_beamcenter_x is None:
            return
        if self.exp_beamcenter_y is None:
            return

        # hoverEvent is only called if
        # either or both maps are active
        # -> always show on isEnter event
        if event.isEnter():
            self.res_label.show()

        # cormap displays the product of both corrections
        # but the individual values can be retrieves from their arrays
        # -> _polcor and _solang
        cu_y,cu_x = map(int, event.pos())
        dx = (cu_x - self.exp_beamcenter_y) * self.exp_pixelsize
        dy = (cu_y - self.exp_beamcenter_x) * self.exp_pixelsize
        tth = np.arctan(np.sqrt(dx**2 + dy**2)/self.exp_distance)
        # Conversion factor keV to Angstrom: 12.398
        # sin(t)/l: np.sin(Theta) / lambda -> (12.398/geo_energy)
        stl = np.sin(tth/2)/self.exp_wavelength

        # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
        dsp = 1/(2*stl)

        _text = f'd[\u212B]: {dsp:.2f}'
        self.res_label.setText(_text)

    def on_treeView_clicked(self, index):
        logging.debug(self.__class__.__name__)
        '''
         Current find run_list implementation is flawed yet simple:
         - iterate over all frames found and find the 001. (SP8) or 
           0001. (APS) entries that mark the beginning of a new run
         - if the first frame is missing THIS WILL CRASH!
           hence, to be fixed later, if ever.
        '''
        self.indexItem = self.model.index(index.row(), 0, index.parent())
        self.curPath = os.path.abspath(self.model.filePath(self.indexItem))
        
        if self.cb_link.isChecked():
            self.le_input.setText(self.curPath)
            self.le_output.setText(self.curPath + self.suffix)
        
        elif self.paths_active == self.le_input:
            self.le_input.setText(self.curPath)
        
        elif self.paths_active == self.le_output:
            self.le_output.setText(self.curPath)
            return
        
        self.tb_convert.setText('Convert Images')
        self.tb_convert.setEnabled(False)

        # find files
        fDir = QtCore.QDir()
        fDir.setPath(self.curPath)
        fDir.setNameFilters(self.exts)
        fDir.setFilter(QtCore.QDir.Filter.Files | QtCore.QDir.Filter.NoDotAndDotDot)
        fDir.setSorting(QtCore.QDir.SortFlag.Name)
        nFrames = fDir.count()
        
        if nFrames > 0:
            self.framesList = [i.absoluteFilePath() for i in fDir.entryInfoList()]
            
            # Check frame format
            self.currentFrame = self.framesList[0]
            self.currentIndex = 0
            if not self.check_format():
                self.currentFrame = None
                return
            
            # Incorrect/Incomplete runs may end in empty self.runList
            # - e.g. if first frame is missing it's not considered a run!
            # - self.fStar is updated by self.check_format()
            self.runList = sorted([os.path.abspath(f) for f in self.framesList if self.fStar in f])
            
            # generate the mask list here would save calling check_format a lot!
            # - getting the run name however is non-trivial due to different naming conventions!
            # - here: simple counting solution - bad idea!
            #self.mList = [os.path.join(self.le_output.text(), '{}_xa_{:>02}_0001.sfrm'.format(self.fStem, i)) for i in range(len((self.runList)))]
            if len(self.runList) == 0:
                return
            
            # clearing and adding to combobox triggers it's .currentIndexChanged()
            # block signals to not call self.mask_change_image_abs
            self.cb_mask_fname.blockSignals(True)
            
            # For some reason removing the last item from the qcombobox (e.g. using clear) and filling it afterwards (addItem/s) crashes the program:
            # *** Terminating app due to uncaught exception 'NSRangeException', reason: '*** -[__NSArrayM objectAtIndexedSubscript:]: index 0 beyond bounds for empty array'
            # *** First throw call stack:
            # libc++abi: terminating due to uncaught exception of type NSException
            #
            # CURRENT WORKAROUND:
            # set max count to 1 -> truncate down to 1
            # set max count to number of items
            # if first item exists -> rename
            # else add item in loop
            #
            self.cb_mask_fname.setMaxCount(1)
            self.cb_mask_fname.setMaxCount(len(self.runList))
            for idx, txt in enumerate([os.path.basename(i) for i in self.runList]):
                if idx == 0 and  self.cb_mask_fname.count() > 0:
                    self.cb_mask_fname.setItemText(idx, txt)
                else:
                    self.cb_mask_fname.addItem(txt)

            # clear combobox
            #self.cb_mask_fname.clear()
            
            # add runs to combobox
            #self.cb_mask_fname.addItems([os.path.basename(i) for i in self.runList])
                    
            self.cb_mask_fname.blockSignals(False)
            
            # if we are here we may allow conversion
            # - the check for the .inf files (SP8 data) is done
            #   by the actual conversion function!
            self.tb_convert.setText('Convert {} Images'.format(nFrames))
            self.tabWidget.setTabEnabled(1, True)
            self.tb_convert.setEnabled(True)
            
        elif self.paths_active == self.le_input:
            self.tabWidget.setTabEnabled(1, False)
        else:
            logging.warning('You should not be able to read this message!')
                    
    def create_output_directory(self, aPath):
        logging.debug(self.__class__.__name__)
        # create output file path
        if not os.path.exists(aPath):
            os.makedirs(aPath)
    
    def disable_user_input(self, toggle):
        logging.debug(self.__class__.__name__)
        self.cb_link.setDisabled(toggle)
        self.cb_overwrite.setDisabled(toggle)
        self.le_input.setDisabled(toggle)
        self.le_output.setDisabled(toggle)
        self.treeView.setDisabled(toggle)
        
    def start_conversion(self):
        logging.debug(self.__class__.__name__)
        '''
          - assign data/sfrm paths
          - get files (again)
          - check for files
          - check if facility is set and proceed accordingly
             - check image name format
             - create sfrm directory
             - start conversion for facility
        '''
        #####################################
        ##     THIS MIGHT BE REDUNDANT     ##
        ## CHECK IF REASSIGNMENT IS NEEDED ##
        #####################################
        # QLineEdit current text to path
        path_input = os.path.abspath(self.le_input.text())
        path_output = os.path.abspath(self.le_output.text())
        #####################################
        
        # check if there are any files
        if not self.framesList:
            self.popup_window('Information', 'No suitable image files found.', 'Please check path.')
            return
        
        # Make directories recursively
        self.create_output_directory(path_output)
        
        # disable main window elements
        # re-enabled after conversion finished
        # -> at end of 'start_conversion'
        self.disable_user_input(True)
        
        # check if overwrite flag is set and
        # pass it on to the conversion function
        overwrite_flag = self.cb_overwrite.isChecked()
        
        # create a pool of workers
        #  - pool.apply_async, map doesn't work since we need to specify the output directory!
        #  - the list 'results' together with 'callback=results.append' is used to track the conversion progress
        #     - 'while len(results) != _todo' constantly checks results to update the progressbar
        #  - pool.close(): close the pool when there are no more files left in 'self.framesList'
        #  - pool.join(): wait for remaining processes to finish
        #########################################
        ##  Add new format identifiers here!   ##
        #########################################
        # fork here according to specified facility
        #  - conversion: what _Utility.py function to call
        #  - parameters: parameters for the conversion function
        #     - path_output, dimension1, dimension2, overwrite_flag
        #     - more if needed, e.g. SP8 2-th correction value
        if self.fSite == 'APS':
            rows, cols, offset, dtype = self.fInfo
            conversion = convert_frame_APS_Bruker
            beamflux = {}
            for f in glob.glob(os.path.join(path_input,'*_flux.txt')):
                with open(f) as ofile:
                    beamflux[int(f.split('_')[-2])] = [int(float(x)) for x in ofile.read().split()[1::2]]
            args = [path_output]
            kwargs = {'rows':rows, 'cols':cols, 'offset':offset, 'overwrite':overwrite_flag, 'beamflux':beamflux}
        elif self.fSite == 'SP8':
            rows, cols, offset, dtype = self.fInfo
            # check data collection timestamp
            with open(self.currentFrame, 'rb') as ofile:
                year = int(re.search(rb'(\d{4}):\d{2}:\d{2}\s+\d{2}:\d{2}:\d{2}', ofile.read(64)).group(1).decode())
            if year < 2019:
                self.SP8_tth_corr = 0.048
            conversion = convert_frame_SP8_Bruker
            args = [path_output]
            # change wavelength
            source_w = None
            if self.action_set_wavelength.isChecked():
                source_w = self.exp_wavelength
            kwargs = {'tth_corr':self.SP8_tth_corr, 'rows':rows, 'cols':cols, 'offset':offset, 'overwrite':overwrite_flag, 'source_w':source_w}
        elif self.fSite == 'SP8_gz':
            rows, cols, offset, dtype = self.fInfo
            # check data collection timestamp
            with gzip.open(self.currentFrame, 'rb') as ofile:
                year = int(re.search(rb'(\d{4}):\d{2}:\d{2}\s+\d{2}:\d{2}:\d{2}', ofile.read(64)).group(1).decode())
            if year < 2019:
                self.SP8_tth_corr = 0.048
            conversion = convert_frame_SP8_Bruker_gz
            args = [path_output]
            # change wavelength
            source_w = None
            if self.action_set_wavelength.isChecked():
                source_w = self.exp_wavelength
            kwargs = {'tth_corr':self.SP8_tth_corr, 'rows':rows, 'cols':cols, 'offset':offset, 'overwrite':overwrite_flag, 'source_w':source_w}
        elif self.fSite == 'DLS':
            rows, cols, offset, dtype = self.fInfo
            conversion = convert_frame_DLS_Bruker
            args = [path_output]
            kwargs = {'rows':rows, 'cols':cols, 'offset':offset, 'overwrite':overwrite_flag}
        else:
            self.popup_window('Information', 'Unknown facility!', '')
            return
        
        self.tb_convert.hide()
        self.pb_convert.show()
        self.statusBar.show()
        
        # Now uses QRunnable and QThreadPool instead of multiprocessing.pool()
        self.num_to_convert = len(self.framesList)
        self.converted = []
        self.pool = QtCore.QThreadPool()
        for fname in self.framesList:
            worker = self.__class__.Threading(conversion, fname, args, kwargs)
            worker.signals.finished.connect(self.conversion_process)
            self.pool.start(worker)
        
        # switch view to mask drawing
        self.tabWidget.setCurrentIndex(1)

    class Threading(QtCore.QRunnable):
        class Signals(QtCore.QObject):
            '''
             Custom signals can only be defined on objects derived from QObject
            '''
            finished = QtCore.pyqtSignal(bool)
    
        def __init__(self, fn_conversion, file_name, fn_args, fn_kwargs):
            '''
             fn_conversion: Conversion function
             file_name:     File name to convert
             fn_args:       Arguments to pass to the function
             fn_kwargs:     Keywords to pass to the function
            '''
            super(self.__class__, self).__init__()
            self.conversion = fn_conversion
            self.name = file_name
            self.args = fn_args
            self.kwargs = fn_kwargs
            self.signals = self.__class__.Signals()
        
        def run(self):
            # conversion: returns True/False
            # signal to conversion_process to track the process
            self.signals.finished.emit(self.conversion(self.name, *self.args, **self.kwargs))
    
    def conversion_process(self, finished):
        self.converted.append(finished)
        num_converted = len(self.converted)
        progress = float(num_converted) / float(self.num_to_convert) * 100.0
        self.pb_convert.setValue(int(round(progress,0)))
        self.status.setText('{}'.format(os.path.basename(self.framesList[num_converted-1])))
        # conversion finished
        if num_converted == self.num_to_convert:
            self.popup_window('Information', 'Successfully converted {} images!'.format(np.count_nonzero(self.converted)), '')
            self.statusBar.hide()
            self.pb_convert.hide()
            self.tb_convert.show()
            # enable main window elements
            self.disable_user_input(False)
        
    def closeEvent(self, event):
        logging.debug(self.__class__.__name__)
        '''
        User clicks the 'x' mark in window
        '''
        self.exitApp()

    def exitApp(self):
        logging.debug(self.__class__.__name__)
        sys.exit()
