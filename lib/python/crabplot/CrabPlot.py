#!/usr/bin/python
# -*- coding: utf-8 -*-
# 

################################
# 
# class CrabPlot()
# 
#   Example: 
#            TODO
# 
#   Last-update: 
#                2017-12-05
#                2018-01-02 plot_line, plot_text
#                2018-05-25 plot_xy errorbarlinestyle
#                2018-10-04 set_xcharsize set_ycharsize
# 
#   Notes:
#     -- to use TeX, we need to install 'texlive-latex-extra' and 'texlive-fonts-recommended'
# 
################################

try:
    import pkg_resources
except ImportError:
    raise SystemExit("Error! Failed to import pkg_resources!")

pkg_resources.require("numpy")
pkg_resources.require("astropy>=1.3")
pkg_resources.require("matplotlib")
#pkg_resources.require("wcsaxes") # http://wcsaxes.readthedocs.io/en/latest/getting_started.html #20180611 commented out because it is merged into astropy

# 
# pip-2.7 install --user --upgrade matplotlib==2.0.1
# ---- a bug in get_grid_positions()
# open /Users/dzliu/Library/Python/2.7/lib/python/site-packages/matplotlib/gridspec.py
#      cellHeights = [cellH] * nrows
# 

import os
import sys
import re
import glob
import inspect
import math
import numpy
import scipy
import scipy.interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import astropy
import astropy.io.ascii as asciitable
#from astropy import units
#from astropy.io import fits
#from astropy.wcs import WCS
#import wcsaxes
#from scipy import spatial
#from scipy.spatial import KDTree
from pprint import pprint
from copy import copy


try: 
    import matplotlib
except ImportError:
    raise SystemExit("Error! Failed to import matplotlib!")

import platform
if sys.version_info < (3, 0):
    if platform.system() == 'Darwin':
        matplotlib.use('Qt5Agg') # sudo port install py27-pyqt5
    else:
        matplotlib.use('TkAgg')

try: 
    from matplotlib import pyplot
except ImportError:
    raise SystemExit("Error! Failed to import pyplot from matplotlib!")

try: 
    from matplotlib.colors import hex2color, rgb2hex
except ImportError:
    raise SystemExit("Error! Failed to import hex2color, rgb2hex from matplotlib.colors!")

try:
    from matplotlib.patches import Ellipse, Circle, Rectangle
except ImportError:
    raise SystemExit("Error! Failed to import Ellipse, Circle, Rectangle from matplotlib.patches!")

try:
    from astropy.visualization import MinMaxInterval, PercentileInterval, AsymmetricPercentileInterval, SqrtStretch, PowerStretch, ImageNormalize
except ImportError:
    raise SystemExit("Error! Failed to import MinMaxInterval, PercentileInterval, AsymmetricPercentileInterval, SqrtStretch, PowerStretch, ImageNormalize from astropy.visualization!")
# ImageNormalize requires astropy>=1.3

import matplotlib.gridspec as gridspec

import matplotlib.ticker as ticker

from matplotlib.ticker import FormatStrFormatter, NullFormatter, AutoMinorLocator, LogLocator



import warnings

warnings.filterwarnings("ignore",".*GUI is implemented.*")

import matplotlib as mpl # https://matplotlib.org/users/customizing.html
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
#mpl.rcParams['axes.grid'] = True
#mpl.rcParams['grid.color'] = 'b0b0b0'
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['grid.linewidth'] = 0.25
mpl.rcParams['grid.alpha'] = 0.8













# 
class CrabPlot(object):
    # 
    def __init__(self, x = None, y = None, xerr = None, yerr = None, xlog = False, ylog = False, xtitle = None, ytitle = None, xrange = [], yrange = [], 
                       image_data = None, image_wcs = None, figure_size = None, figure_dpi = 90, position = None, label = '', 
                       N_panel_per_row = None, N_panel_per_col = None):
        # 
        # open plot
        if figure_size is None:
            self.Figure_size = self.default_figure_size()
        else:
            self.Figure_size = figure_size
        if figure_dpi is None:
            self.Figure_dpi = self.default_figure_dpi()
        else:
            self.Figure_dpi = figure_dpi
        self.Plot_device = pyplot.figure(figsize=self.Figure_size, dpi=self.Figure_dpi) # set figure size 9.0 x 8.0 inches, 90 pixels per inch. 
        self.Plot_panels = [] # each item is a dict {'label': '', 'panel': axis, 'position': position}
        self.Plot_grids = None
        self.Plot_data = {}
        self.Annotation = []
        self.Image_data = image_data
        self.Image_wcs = image_wcs
        self.N_panel_per_row = N_panel_per_row
        self.N_panel_per_col = N_panel_per_col
        # 
        # <TODO> position = [0.10, 0.10, 0.85, 0.85]
        # 
        # plot image
        if image_data:
            self.plot_image(position = position, image_data = image_data, image_wcs = image_wcs)
        elif x is not None and y is not None:
            if len(x)>0 and len(y)==len(x):
                self.plot_xy(position = position, x = x, y = y, xerr = xerr, yerr = yerr, xlog = xlog, ylog = ylog, xtitle = xtitle, ytitle = ytitle, xrange = xrange, yrange = yrange)
        # 
        # world
        self.World = {}
        self.World['My Type'] = 'CrabPlot'
        # 
        # get variable name 
        # -- see http://stackoverflow.com/questions/1690400/getting-an-instance-name-inside-class-init
        self.World['My Name'] = ""
        self.World['My Names'] = []
        tmp_frame = inspect.currentframe().f_back
        tmp_variables = tmp_frame.f_globals.copy()
        tmp_variables.update(tmp_frame.f_locals)
        #tmp_variables = dict(tmp_frame.f_globals.items() + tmp_frame.f_locals.items())
        for tmp_name, tmp_variable in tmp_variables.items():
            if isinstance(tmp_variable, self.__class__):
                if hash(self) == hash(tmp_variable):
                    self.World['My Names'].append(tmp_name)
        if len(self.World['My Names']) > 0:
            self.World['My Name'] = self.World['My Names'][0]
    # 
    def check_array(self, arr):
        # check array type numpy.ndarray
        if type(arr) is not list and type(arr) is not numpy.ndarray:
            arr = [arr]
        if type(arr) is list:
            arr = numpy.array(arr)
        return arr
    # 
    def check_scalar(self, arr):
        # check scalar
        if type(arr) is list or type(arr) is numpy.ndarray:
            if len(arr) > 0:
                arr = arr[0]
        if type(arr) is list or type(arr) is numpy.ndarray:
            if len(arr) > 0:
                arr = arr[0]
        return arr
    # 
    def check_label(self, label):
        # check label, set to the number of panels if it is None
        if label is None:
            label = '%d'%(len(self.Plot_panels)+1)
            if len(self.Plot_panels) > 0:
                while label in [t['label'] for t in self.Plot_panels]:
                    label = '%s_'%(label)
        return label
    # 
    def get_grid_position_in_gridspec(self):
        # 
        # override gridspec.py get_grid_position()
        # see https://github.com/Goldcap/Constellation/blob/master/trunk/testing/matplotlib-1.0.1/build/lib.linux-x86_64-2.6/matplotlib/gridspec.py
        # 
        debug = 0
        # 
        nrows, ncols = self.Plot_grids.get_geometry()
        #nrows = long(nrows)
        #ncols = long(ncols)
        nrows = int(nrows) # Python 3 has no 'long' type anymore
        ncols = int(ncols) # Python 3 has no 'long' type anymore
        # 
        subplot_params = self.Plot_grids.get_subplot_params(self.Plot_device)
        left = subplot_params.left # the full image left,right,top,bottom margin
        right = subplot_params.right # the full image left,right,top,bottom margin
        top = subplot_params.top # the full image left,right,top,bottom margin
        bottom = subplot_params.bottom # the full image left,right,top,bottom margin
        wspace = subplot_params.wspace # horizontal cell spacing
        hspace = subplot_params.hspace # vertical cell spacing
        totWidth = right-left # all-panel width
        totHeight = top-bottom # all-panel height
        # 
        # calculate accumulated heights of columns
        cellH = totHeight/(nrows+hspace*(nrows-1)) # assuming panels are uniformly distributed
        sepH = hspace*cellH # hspace is relative to each cell height
        # 
        if self.Plot_grids._row_height_ratios is not None:
            netHeight = cellH * nrows
            tr = float(sum(self.Plot_grids._row_height_ratios))
            cellHeights = [netHeight*r/tr for r in self.Plot_grids._row_height_ratios]
        else:
            cellHeights = [cellH] * nrows
        
        sepHeights = [0] + ([sepH] * (nrows-1))
        
        cellHs = numpy.add.accumulate(numpy.ravel(list(zip(sepHeights, cellHeights))))
        # 20180112 dzliu - adding margins (the margins are also relative to the size of each cell, similar to hspace and wspace)
        margin_top_of_current_row = 0.0
        margin_bottom_of_current_row= 0.0
        for j in range(nrows):
            margin_top_of_cells_in_a_row = [0.0]
            margin_bottom_of_cells_in_a_row = [0.0]
            for k in range(ncols):
                i = j * ncols + k
                if i < len(self.Plot_panels):
                    margin_top_of_cells_in_a_row.append(self.Plot_panels[i]['margin-top'])
                    margin_bottom_of_cells_in_a_row.append(self.Plot_panels[i]['margin-bottom'])
            margin_top_of_current_row = numpy.nanmax(numpy.array(margin_top_of_cells_in_a_row))
            margin_bottom_of_current_row = numpy.nanmax(numpy.array(margin_bottom_of_cells_in_a_row))
            if debug>=1:
                print('CrabPlot::get_grid_position_in_gridspec() adding vertical margins for cell at col %d row %d: top = %s, bottom = %s'%(k+1, j+1, margin_top_of_current_row, margin_bottom_of_current_row))
            # now we consider the margin-top/bottom of each cell
            cellHs[2*j+0] = cellHs[2*j+0] + margin_top_of_current_row*cellH
            cellHs[2*j+1] = cellHs[2*j+1] - margin_bottom_of_current_row*cellH
        # -- dzliu -- now added margins
        
        # calculate accumulated widths of rows
        cellW = totWidth/(ncols + wspace*(ncols-1))
        sepW = wspace*cellW
        
        if self.Plot_grids._col_width_ratios is not None:
            netWidth = cellW * ncols
            tr = float(sum(self.Plot_grids._col_width_ratios))
            cellWidths = [netWidth*r/tr for r in self.Plot_grids._col_width_ratios]
        else:
            cellWidths = [cellW] * ncols
        
        sepWidths = [0] + ([sepW] * (ncols-1))
        
        cellWs = numpy.add.accumulate(numpy.ravel(list(zip(sepWidths, cellWidths))))
        # 20180112 dzliu - adding margins (the margins are also relative to the size of each cell, similar to hspace and wspace)
        margin_left_of_current_col = 0.0
        margin_right_of_current_col= 0.0
        for k in range(ncols):
            margin_left_of_cells_in_a_col = [0.0]
            margin_right_of_cells_in_a_col = [0.0]
            for j in range(nrows):
                i = j * ncols + k
                if i < len(self.Plot_panels):
                    margin_left_of_cells_in_a_col.append(self.Plot_panels[i]['margin-left'])
                    margin_right_of_cells_in_a_col.append(self.Plot_panels[i]['margin-right'])
            margin_left_of_current_col = numpy.nanmax(numpy.array(margin_left_of_cells_in_a_col))
            margin_right_of_current_col = numpy.nanmax(numpy.array(margin_right_of_cells_in_a_col))
            if debug>=1:
                print('CrabPlot::get_grid_position_in_gridspec() adding horizontal margins for cell at col %d row %d: left = %s, right = %s'%(k+1, j+1, margin_left_of_current_col, margin_right_of_current_col))
            # now we consider the margin-left/right of each cell
            cellWs[2*k+0] = cellWs[2*k+0] + margin_left_of_current_col*cellW
            cellWs[2*k+1] = cellWs[2*k+1] - margin_right_of_current_col*cellW
        # -- dzliu -- now added margins
        
        # 
        # debug
        if debug>=1:
            print('CrabPlot::get_grid_position_in_gridspec() cell position x coordinates (lower, upper) pairs: ', cellHs)
            print('CrabPlot::get_grid_position_in_gridspec() cell position y coordinates (lower, upper) pairs: ', cellWs)
        
        figTops = [top - cellHs[2*rowNum] for rowNum in range(nrows)]
        figBottoms = [top - cellHs[2*rowNum+1] for rowNum in range(nrows)]
        figLefts = [left + cellWs[2*colNum] for colNum in range(ncols)]
        figRights = [left + cellWs[2*colNum+1] for colNum in range(ncols)]
        
        # 
        # debug
        if debug>=1:
            print('CrabPlot::get_grid_position_in_gridspec() figBottoms: ', figBottoms)
            print('CrabPlot::get_grid_position_in_gridspec() figTops: ', figTops)
            print('CrabPlot::get_grid_position_in_gridspec() figLefts: ', figLefts)
            print('CrabPlot::get_grid_position_in_gridspec() figRights: ', figRights)
        
        return figBottoms, figTops, figLefts, figRights
    # 
    def get_panel_position_in_gridspec(self, num1, num2=None, return_all=False):
        # 
        # override gridspec.py get_position()
        # num1 is the index, if num2 is given, then select [num1:num2]
        # 
        nrows, ncols = self.Plot_grids.get_geometry()
        #nrows = long(nrows)
        #ncols = long(ncols)
        nrows = int(nrows) # Python 3 has no 'long' type anymore
        ncols = int(ncols) # Python 3 has no 'long' type anymore
        
        #print('CrabPlot::get_panel_position_in_gridspec() left, bottom, top, right', left, bottom, top, right)
        figBottoms, figTops, figLefts, figRights = self.get_grid_position_in_gridspec()
        rowNum, colNum =  divmod(num1, ncols)
        figBottom = figBottoms[rowNum]
        figTop = figTops[rowNum]
        figLeft = figLefts[colNum]
        figRight = figRights[colNum]
        
        if num2 is not None:
            
            rowNum2, colNum2 =  divmod(num2, ncols)
            figBottom2 = figBottoms[rowNum2]
            figTop2 = figTops[rowNum2]
            figLeft2 = figLefts[colNum2]
            figRight2 = figRights[colNum2]
            
            figBottom = min(figBottom, figBottom2)
            figLeft = min(figLeft, figLeft2)
            figTop = max(figTop, figTop2)
            figRight = max(figRight, figRight2)
            
        figbox = matplotlib.transforms.Bbox.from_extents(figLeft, figBottom, figRight, figTop)
        
        if return_all:
            return figbox, rowNum, colNum, nrows, ncols
        else:
            return figbox
    # 
    #def get_left_right_top_bottom_margins(self, i):
    #    # 20180112 load left, right, top, bottom margins for the panel i.
    #    if i < len(self.Plot_panels):
    #        if left is not None:
    #            if not numpy.isnan(left):
    #                self.Plot_panels[i]['margin-left'] = left
    #                print('CrabPlot::get_left_right_top_bottom_margins() Setting left margin to %s for panel index %d'%(left, i))
    #            else:
    #                del self.Plot_panels[i]['margin-left']
    #                left = None
    #        else:
    #            if 'margin-left' in self.Plot_panels[i]:
    #                left = self.Plot_panels[i]['margin-left']
    #            else:
    #                left = None
    #        # 
    #        if right is not None:
    #            if not numpy.isnan(right):
    #                self.Plot_panels[i]['margin-right'] = right
    #                print('CrabPlot::get_left_right_top_bottom_margins() Setting right margin to %s for panel index %d'%(right, i))
    #            else:
    #                del self.Plot_panels[i]['margin-right']
    #                right = None
    #        else:
    #            if 'margin-right' in self.Plot_panels[i]:
    #                right = self.Plot_panels[i]['margin-right']
    #            else:
    #                right = None
    #        # 
    #        if top is not None:
    #            if not numpy.isnan(top):
    #                self.Plot_panels[i]['margin-top'] = top
    #                print('CrabPlot::get_left_right_top_bottom_margins() Setting top margin to %s for panel index %d'%(top, i))
    #            else:
    #                del self.Plot_panels[i]['margin-top']
    #                top = None
    #        else:
    #            if 'margin-top' in self.Plot_panels[i]:
    #                top = self.Plot_panels[i]['margin-top']
    #            else:
    #                top = None
    #        # 
    #        if bottom is not None:
    #            if not numpy.isnan(bottom):
    #                self.Plot_panels[i]['margin-bottom'] = bottom
    #                print('CrabPlot::get_left_right_top_bottom_margins() Setting bottom margin to %s for panel index %d'%(bottom, i))
    #            else:
    #                del self.Plot_panels[i]['margin-bottom']
    #                bottom = None
    #        else:
    #            if 'margin-bottom' in self.Plot_panels[i]:
    #                bottom = self.Plot_panels[i]['margin-bottom']
    #            else:
    #                bottom = None
    #    return left, right, top, bottom
    # 
    def add_panel(self, position = None, label = None, projection = None, xrange = [], yrange = [], debug = False):
        # aim: 
        #     if position is None, add one panel and adjust panels to new grid, ignoring panels with given position. 
        #     if position is not None, then add one panel according to the position without adjust other panels. 
        if position is None:
            # check previous panels without given position
            n_panel = 0
            for i in range(len(self.Plot_panels)):
                if self.Plot_panels[i]['position'] is None:
                    n_panel = n_panel + 1
            # add one panel number
            n_panel = n_panel + 1
            # make grid number
            grid_ny = int(numpy.sqrt(float(n_panel)))
            grid_nx = int(numpy.round((float(n_panel)/grid_ny)))
            while (n_panel > grid_nx*grid_ny):
                grid_nx = grid_nx + 1
            # but if user provided self.N_panel_per_row, then fix to it
            if self.N_panel_per_row is not None:
                grid_nx = self.N_panel_per_row
                grid_ny = numpy.round((float(n_panel)/grid_nx))
                while (n_panel > grid_nx*grid_ny):
                    grid_ny = grid_ny + 1
            # but also if user provided self.N_panel_per_row, then fix to it (in priori)
            if self.N_panel_per_col is not None:
                grid_ny = self.N_panel_per_col
                grid_nx = numpy.round((float(n_panel)/grid_ny))
                while (n_panel > grid_nx*grid_ny):
                    grid_nx = grid_nx + 1
            # make gridspec
            self.Plot_grids = gridspec.GridSpec(grid_nx, grid_ny)
            # debug
            if debug:
                print('CrabPlot::add_panel() DEBUG: add_subplot(%d, %d, %d)'%(grid_nx, grid_ny, n_panel))
                print('CrabPlot::add_panel() DEBUG: self.Plot_grids.get_geometry() = %s %s'%(self.Plot_grids.get_geometry()))
            # adjust previous panels
            for i in range(len(self.Plot_panels)):
                if self.Plot_panels[i]['position'] is None:
                    if debug:
                        print('CrabPlot::add_panel() DEBUG: Plot_panels[%d][\'panel\'].get_geometry() = %s'%(i, self.Plot_panels[i]['panel'].get_geometry()))
                    #<1># self.Plot_panels[i]['panel'].set_position(self.Plot_grids[i].get_position(self.Plot_device)) # see -- http://stackoverflow.com/questions/22881301/changing-matplotlib-subplot-size-position-after-axes-creation
                    #<1># self.Plot_panels[i]['panel'].set_subplotspec(self.Plot_grids[i])
                    #<2># self.Plot_panels[i]['panel'].change_geometry(grid_nx, grid_ny, i+1) # see -- http://stackoverflow.com/questions/22881301/changing-matplotlib-subplot-size-position-after-axes-creation
                    # 
                    # dzliu overriding gridspec.py get_position()
                    prev_position_in_gridspec = self.get_panel_position_in_gridspec(i)
                    self.Plot_panels[i]['panel'].set_position(prev_position_in_gridspec) # see -- http://stackoverflow.com/questions/22881301/changing-matplotlib-subplot-size-position-after-axes-creation
                    self.Plot_panels[i]['panel'].set_subplotspec(self.Plot_grids[i])
                    
            # add new panel
            panel = self.Plot_device.add_subplot(grid_nx, grid_ny, n_panel)
        else:
            panel = self.Plot_device.add_axes(position, projection = projection)
        # tick label format
        #panel.ticklabel_format(axis='both', style='sci', scilimits=(-2,+2)) # not working for log!
        # determine label
        label = self.check_label(label)
        # append to self.Plot_panels
        self.Plot_panels.append( { 'label': label, 'panel': panel, 'position': position, 
                                   'x': None, 'y': None, 'xerr': None, 'yerr': None, 'xlog': False, 'ylog': False, 'xrange': xrange, 'yrange': yrange, 
                                   'image_data': None, 'image_wcs': None, 
                                   'margin-left':0.0, 'margin-bottom':0.0, 'margin-right':0.0, 'margin-top':0.0 } )
        # set ticks
        self.set_xticksize(panel=len(self.Plot_panels))
        self.set_yticksize(panel=len(self.Plot_panels))
        # 
        # dzliu overriding gridspec.py get_position()
        i = len(self.Plot_panels)-1
        current_position_in_gridspec = self.get_panel_position_in_gridspec(i)
        self.Plot_panels[i]['panel'].set_position(current_position_in_gridspec)
        self.Plot_panels[i]['panel'].set_subplotspec(self.Plot_grids[i])
        #self.Plot_panels[i]['panel'].xaxis.tick_top()
        #self.Plot_panels[i]['panel'].yaxis.tick_right()
        # return
        return self.Plot_panels[-1]
    # 
    def default_figsize(self):
        return (9.0,8.0)
    # 
    def default_figure_size(self):
        return self.default_figsize()
    # 
    def default_figure_dpi(self):
        return 90
    # 
    def default_fontname_for_axis_title(self):
        return ['NGC','Helvetica','sans-serif']
    # 
    def default_fontsize_for_axis_title(self):
        return 20 - (9.0-numpy.max(self.Figure_size)) / (9.0-5.0) * (18-14)
        # if figure size 9.0inch, then font size 18
        # if figure size 5.0inch, then font size 14
    # 
    def default_fontname_for_axis_ticks(self):
        return ['NGC','Helvetica','sans-serif']
    # 
    def default_fontsize_for_axis_ticks(self):
        return 18.5 - (9.0-numpy.max(self.Figure_size)) / (9.0-5.0) * (16.5-12.5)
        # if figure size 9.0inch, then font size 16.5
        # if figure size 5.0inch, then font size 12.5
    # 
    # def set_xrange
    def set_xrange(self, xrange, ax=None, panel=None):
        if len(xrange)<2:
            return
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    # set_xrange for the last panel
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['xrange'] = xrange
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_xrange for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['xrange'] = xrange
                else:
                    return
            if ax is not None:
                ax.set_xlim(xrange)
    # 
    # def set_yrange
    def set_yrange(self, yrange, ax=None, panel=None):
        if len(yrange)<2:
            return
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    # set_yrange for the last panel
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['yrange'] = yrange
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_yrange for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['yrange'] = yrange
                else:
                    return
            if ax is not None:
                ax.set_ylim(yrange)
    # 
    # def set_xtitle
    def set_xtitle(self, xtitle, ax=None, panel=None, fontsize=None, fontname=None, **kwargs):
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = -1 # if panel is None, we set for the last active panel
                if panel == 0:
                    # set_xtitle for all panels
                    for panel in range(1,len(self.Plot_panels)+1,1):
                        self.set_xtitle(xtitle, panel=panel, fontsize=fontsize, fontname=fontname, **kwargs)
                elif panel < 0:
                    # set_xtitle for the last panel
                    ax = self.Plot_panels[len(self.Plot_panels)+panel]['panel']
                    self.Plot_panels[len(self.Plot_panels)+panel]['xtitle'] = xtitle
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_xtitle for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['xtitle'] = xtitle
                else:
                    return
            reset_ticksfont = False
            if fontsize is None:
                fontsize = self.default_fontsize_for_axis_title()
                reset_ticksfont = True
            if fontname is None:
                fontname = self.default_fontname_for_axis_title()
            if ax is not None:
                ax.set_xlabel(xtitle, fontsize=fontsize, fontname=fontname, **kwargs)
                if reset_ticksfont == True:
                    self.set_xticksfont(fontsize=fontsize-1) # when setting title fontsize, we will also reset ticks labelsize.
    # 
    # def set_ytitle
    def set_ytitle(self, ytitle, ax=None, panel=None, fontsize=None, fontname=None, **kwargs):
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = -1 # if panel is None, we set for the last active panel
                if panel == 0:
                    # set_ytitle for all panels
                    for panel in range(1,len(self.Plot_panels)+1,1):
                        self.set_ytitle(ytitle, panel=panel, fontsize=fontsize, fontname=fontname, **kwargs)
                elif panel < 0:
                    # set_ytitle for the last panel
                    ax = self.Plot_panels[len(self.Plot_panels)+panel]['panel']
                    self.Plot_panels[len(self.Plot_panels)+panel]['ytitle'] = ytitle
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_ytitle for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['ytitle'] = ytitle
                else:
                    return
            reset_ticksfont = False
            if fontsize is None:
                fontsize = self.default_fontsize_for_axis_title()
                reset_ticksfont = True
            if fontname is None:
                fontname = self.default_fontname_for_axis_title()
            if ax is not None:
                ax.set_ylabel(ytitle, fontsize=fontsize, fontname=fontname, **kwargs)
                if reset_ticksfont == True:
                    self.set_yticksfont(fontsize=fontsize-1) # when setting title fontsize, we will also reset ticks labelsize.
    # 
    # def set_xtickfont
    def set_xtickfont(self, ax=None, panel=None, fontsize=None, fontname=None):
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = -1 # if panel is None, we set for the last active panel
                if panel == 0:
                    # set_xtickfont for all panels
                    for panel in range(1,len(self.Plot_panels)+1,1):
                        self.set_xtickfont(panel=panel, fontsize=fontsize, fontname=fontname)
                elif panel < 0:
                    # set_xtickfont for the last panel
                    ax = self.Plot_panels[len(self.Plot_panels)+panel]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_xtickfont for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                else:
                    return
            if fontsize is None:
                fontsize = self.default_fontsize_for_axis_ticks()
            if fontname is None:
                fontname = self.default_fontname_for_axis_ticks()
            if ax is not None:
                if fontname != '' or fontsize > 0:
                    for label in ax.get_xticklabels():
                        if fontname != '': label.set_family(fontname)
                        if fontsize > 0: label.set_fontsize(fontsize)
        return
    def set_xticksfont(self, ax=None, panel=None, fontsize=None, fontname=None):
        self.set_xtickfont(panel=panel, fontsize=fontsize, fontname=fontname)
        return
    # 
    # def set_ytickfont
    def set_ytickfont(self, ax=None, panel=None, fontsize=None, fontname=None):
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = -1 # if panel is None, we set for the last active panel
                if panel == 0:
                    # set_ytickfont for all panels
                    for panel in range(1,len(self.Plot_panels)+1,1):
                        self.set_ytickfont(panel=panel, fontsize=fontsize, fontname=fontname)
                elif panel < 0:
                    # set_ytickfont for the last panel
                    ax = self.Plot_panels[len(self.Plot_panels)+panel]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_ytickfont for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                else:
                    return
            if fontsize is None:
                fontsize = self.default_fontsize_for_axis_ticks()
            if fontname is None:
                fontname = self.default_fontname_for_axis_ticks()
            if ax is not None:
                if fontname != '' or fontsize > 0:
                    for label in ax.get_yticklabels():
                        if fontname != '': label.set_family(fontname)
                        if fontsize > 0: label.set_fontsize(fontsize)
        return
    def set_yticksfont(self, ax=None, panel=None, fontsize=None, fontname=None):
        self.set_ytickfont(panel=panel, fontsize=fontsize, fontname=fontname)
        return
    # 
    # def set_xticksize
    def set_xticksize(self, ax=None, panel=None, majorsize=None, minorsize=None, majorwidth=None, minorwidth=None, majordirection='in', minordirection='in'):
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = -1 # if panel is None, we set for the last active panel
                if panel == 0:
                    # set_xticksize for all panels
                    for panel in range(1,len(self.Plot_panels)+1,1):
                        self.set_xticksize(panel=panel, majorsize=majorsize, minorsize=minorsize, majorwidth=majorwidth, minorwidth=minorwidth, majordirection=majordirection, minordirection=minordirection)
                elif panel < 0:
                    # set_xticksize for the last panel
                    ax = self.Plot_panels[len(self.Plot_panels)+panel]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_xticksize for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                else:
                    return
        if majorsize is None:
            majorsize = 10
        if minorsize is None:
            minorsize = 5
        if ax is not None:
            ax.tick_params(axis='x', which='major', direction=majordirection, length=majorsize, width=majorwidth)
            ax.tick_params(axis='x', which='minor', direction=minordirection, length=minorsize, width=minorwidth)
    # 
    # def set_yticksize
    def set_yticksize(self, ax=None, panel=None, majorsize=None, minorsize=None, majorwidth=None, minorwidth=None, majordirection='in', minordirection='in'):
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = -1 # if panel is None, we set for the last active panel
                if panel == 0:
                    # set_yticksize for all panels
                    for panel in range(1,len(self.Plot_panels)+1,1):
                        self.set_yticksize(panel=panel, majorsize=majorsize, minorsize=minorsize, majorwidth=majorwidth, minorwidth=minorwidth, majordirection=majordirection, minordirection=minordirection)
                elif panel < 0:
                    # set_yticksize for the last panel
                    ax = self.Plot_panels[len(self.Plot_panels)+panel]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_yticksize for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                else:
                    return
        if majorsize is None:
            majorsize = 10
        if minorsize is None:
            minorsize = 5
        if ax is not None:
            ax.tick_params(axis='y', which='major', direction=majordirection, length=majorsize, width=majorwidth)
            ax.tick_params(axis='y', which='minor', direction=minordirection, length=minorsize, width=minorwidth)
    # 
    # def set_xcharsize
    def set_xcharsize(self, ax=None, panel=None, charsize=None, minortickscharsize=None, axislabelcharsize=None):
        # charsize is major ticks char size
        # minortickscharsize is minor ticks char size
        # axislabelcharsize is axis label (i.e. title) char size
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = -1 # if panel is None, we set for the last active panel
                if panel == 0:
                    # set_xcharsize for all panels
                    for panel in range(1,len(self.Plot_panels)+1,1):
                        self.set_xcharsize(panel=panel, charsize=charsize, minortickscharsize=minortickscharsize, axislabelcharsize=axislabelcharsize)
                elif panel < 0:
                    # set_xcharsize for the last panel
                    ax = self.Plot_panels[len(self.Plot_panels)+panel]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_xcharsize for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                else:
                    return
        if ax is not None:
            if charsize is not None:
                ax.tick_params(axis='x', which='major', labelsize=charsize)
            if minortickscharsize is not None:
                ax.tick_params(axis='x', which='minor', labelsize=minortickscharsize)
            if axislabelcharsize is not None:
                if ax.get_xlabel() is not None:
                    ax.set_xlabel(ax.get_xlabel(), fontsize=axislabelcharsize)
    # 
    # def set_ycharsize
    def set_ycharsize(self, ax=None, panel=None, charsize=None, minortickscharsize=None, axislabelcharsize=None):
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = -1 # if panel is None, we set for the last active panel
                if panel == 0:
                    # set_ycharsize for all panels
                    for panel in range(1,len(self.Plot_panels)+1,1):
                        self.set_ycharsize(panel=panel, charsize=charsize, minortickscharsize=minortickscharsize, axislabelcharsize=axislabelcharsize)
                elif panel < 0:
                    # set_ycharsize for the last panel
                    ax = self.Plot_panels[len(self.Plot_panels)+panel]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
                    # set_ycharsize for the specified panel
                    ax = self.Plot_panels[panel-1]['panel']
                else:
                    return
        if ax is not None:
            if charsize is not None:
                ax.tick_params(axis='y', which='major', labelsize=charsize)
            if minortickscharsize is not None:
                ax.tick_params(axis='y', which='minor', labelsize=minortickscharsize)
            if axislabelcharsize is not None:
                if ax.get_ylabel() is not None:
                    ax.set_ylabel(ax.get_ylabel(), fontsize=axislabelcharsize)
    # 
    # def set_grid_hspace
    def set_grid_hspace(self, hspace = 0.05):
        if len(self.Plot_panels)>0:
            self.Plot_grids.update(hspace=hspace)
    # 
    # def set_grid_vspace
    def set_grid_wspace(self, wspace = 0.05):
        if len(self.Plot_panels)>0:
            self.Plot_grids.update(wspace=wspace)
    # 
    # def set_figure_margin
    def set_figure_margin(self, left=None, bottom=None, right=None, top=None):
        # set the margin to the full image
        self.Plot_device.subplots_adjust(left=left, bottom=bottom, right=right, top=top)
        print('CrabPlot::set_figure_margin() left=%s, bottom=%s, right=%s, top=%s'%(left, bottom, right, top))
    # 
    # def set_panel_margin
    def set_panel_margin(self, i, left=None, bottom=None, right=None, top=None):
        if i >= 0 and i < len(self.Plot_panels):
            if left is not None:
                self.Plot_panels[i]['margin-left'] = float(left)
            if bottom is not None:
                self.Plot_panels[i]['margin-bottom'] = float(bottom)
            if right is not None:
                self.Plot_panels[i]['margin-right'] = float(right)
            if top is not None:
                self.Plot_panels[i]['margin-top'] = float(top)
        # then the figure will be updated only when called get_panel_position_in_gridspec
        ## dzliu overriding gridspec.py get_position()
        #current_position_in_gridspec = self.get_panel_position_in_gridspec(i)
        #self.Plot_panels[i]['panel'].set_position(current_position_in_gridspec)
        #self.Plot_panels[i]['panel'].set_subplotspec(self.Plot_grids[i])
    # 
    # def set_margin
    def set_margin(self, ax=None, panel=None, left=None, bottom=None, right=None, top=None):
        # if ax is not None, then get panel number and set margin for it
        # in this case, we must need len(self.Plot_panels)>0
        if ax is not None and len(self.Plot_panels)>0:
            panel = None
            if panel is None:
                for i in range(len(self.Plot_panels)):
                    if ax == self.Plot_panels[i]:
                        panel = i+1
            if panel is None:
                raise Exception('Error! CrabPlot::set_margin() failed to set margin for ax = %s! It is not in self.Plot_panels?'%(ax))
                return
        # else if panel number is given
        elif panel is not None and len(self.Plot_panels)>0:
            # and if panel>0, then set panel margin
            if panel>0:
                self.set_panel_margin(panel-1, left=left, bottom=bottom, right=right, top=top)
            # otherwise if panel<=0, set figure margin
            else:
                self.set_figure_margin(left=left, bottom=bottom, right=right, top=top)
        # else if nothing is given, set margins to the last panel
        # in this case, we must also need len(self.Plot_panels)>0
        elif len(self.Plot_panels)>0:
            self.set_panel_margin(len(self.Plot_panels)-1, left=left, bottom=bottom, right=right, top=top)
        # otherwise set figure margin
        else:
            self.set_figure_margin(left=left, bottom=bottom, right=right, top=top)
    # 
    # def totuple
    def convert_array_to_tuple(self,a):
        try:
            return tuple(self.convert_array_to_tuple(i) for i in a)
        except TypeError:
            return a
    # 
    def get_color_by_value(self, data_array, min_value = None, max_value = None, input_value = None, log = False, cmap = matplotlib.cm.cool):
        # add color bar by data_array
        data_arr = numpy.array(data_array)
        if min_value is None:
            min_value = numpy.min(data_arr)
        if max_value is None:
            max_value = numpy.max(data_arr)
        if type(log) is bool:
            log = int(log)
        if log > 0:
            data_norm = matplotlib.colors.LogNorm(vmin=min_value, vmax=max_value, clip=True)
        else:
            data_norm = matplotlib.colors.Normalize(vmin=min_value, vmax=max_value, clip=True)
        #print('vmin=', min_value, 'vmax=', max_value, 'log=', log)
        color_mapper = matplotlib.cm.ScalarMappable(norm=data_norm, cmap=cmap)
        if input_value is None:
            color_by_data_array = self.convert_array_to_tuple(color_mapper.to_rgba(data_arr))
            return color_by_data_array
        else:
            if type(input_value) is not list or type(input_value) is not numpy.array or type(input_value) is not numpy.ndarray:
                color_by_input_value = self.convert_array_to_tuple(color_mapper.to_rgba(numpy.array([input_value])))
                color_by_input_value = color_by_input_value[0]
            else:
                color_by_input_value = self.convert_array_to_tuple(color_mapper.to_rgba(input_value))
            return color_by_input_value
    # 
    def get_panel_ax(self, ax, current, overplot):
        # this function solves the combination of the three input arguments: 
        #   ax (a matplotlib axis variable), 
        #   current (a number), 
        #   and overplot (a boolean value)
        # 
        # check inputs
        if current is None:
            current = 0
        if overplot is None:
            overplot = False
        # 
        # check if a direct matplotlib axis var is given or not
        if ax is not None:
            plot_panel_ax = ax
            plot_panel_xy = None
            current = -1
            overplot = True
            for i in range(len(self.Plot_panels)):
                if self.Plot_panels[i]['panel'] == ax:
                    current = i+1
                    break
        else:
            # get current panel or add new panel
            if len(self.Plot_panels) == 0:
                plot_panel_xy = self.add_panel()
                plot_panel_ax = plot_panel_xy['panel']
                current = len(self.Plot_panels)
                overplot = False
            # get current panel
            elif current > 0: 
                plot_panel_xy = self.Plot_panels[current-1]
                plot_panel_ax = plot_panel_xy['panel']
            # get last panel for overplotting
            elif overplot: 
                plot_panel_xy = self.Plot_panels[len(self.Plot_panels)-1]
                plot_panel_ax = plot_panel_xy['panel']
                current = len(self.Plot_panels)
            else:
                plot_panel_xy = self.add_panel()
                plot_panel_ax = plot_panel_xy['panel']
                current = len(self.Plot_panels)
                overplot = False
        # 
        # tweak ax
        #plot_panel_ax.ticklabel_format(useOffset=False) # only works for ScalarFormatter
        # 
        # return
        return plot_panel_ax, current, overplot
    # 
    def get_upper_limit_marker(self):
        # see -- https://matplotlib.org/users/path_tutorial.html
        vertices = []
        codes = []
        vertices.append((-1.0,+0.0)) ; codes.append(matplotlib.path.Path.MOVETO)
        vertices.append((+1.0,+0.0)) ; codes.append(matplotlib.path.Path.LINETO)
        vertices.append((+0.0,+0.0)) ; codes.append(matplotlib.path.Path.MOVETO)
        vertices.append((+0.0,-3.0)) ; codes.append(matplotlib.path.Path.LINETO)
        vertices.append((-0.8,-1.9)) ; codes.append(matplotlib.path.Path.LINETO)
        vertices.append((+0.0,-3.0)) ; codes.append(matplotlib.path.Path.LINETO)
        vertices.append((+0.8,-1.9)) ; codes.append(matplotlib.path.Path.LINETO)
        vertices.append((+0.0,-3.0)) ; codes.append(matplotlib.path.Path.MOVETO)
        vertices.append((+0.0,+0.0)) ; codes.append(matplotlib.path.Path.STOP)
        return matplotlib.path.Path(vertices, codes)
    # 
    def plot_xy(self, x, y, xerr = None, yerr = None, xlog = None, ylog = None, xrange = None, yrange = None, 
                ax = None, NormalizedCoordinate = False, overplot = False, current = 0, 
                position = None, label = None, 
                dataname = '', 
                xtitle = None, ytitle = None, 
                fmt = None, 
                symbol = 'o', symsize = 3, thick = 2, 
                marker = None, size = None, color = '#1b81e5', fillstyle = None, 
                linestyle = 'None', linewidth = None, drawstyle = None, 
                facecolor = None, edgecolor = None, edgewidth = None, alpha = 1.0, zorder = 5, 
                errorbarlinestyle = '', 
                uplims = None, lolims = None, 
                margin = None, 
                grid = None, 
                verbose = 1, 
                **kwargs):
        # inputs
        # -- ax is a direct 
        # 
        # check x y type
        x = self.check_array(x)
        y = self.check_array(y)
        # store x y
        if dataname != '':
            self.Plot_data[dataname] = numpy.column_stack((x,y))
            if verbose >= 1:
                print('CrabPlot::plot_xy() Stored the input x and y into self.Plot_data[%s]!'%(dataname))
        # check x y dimension
        if x.shape != y.shape:
            print('CrabPlot::plot_xy() Error! The input x and y do not have the same shape!')
            return
        # get panel ax
        plot_panel_ax, current, overplot = self.get_panel_ax(ax, current, overplot)
        # set grid
        if grid is not None:
            plot_panel_ax.grid(grid)
        # set margin if x y title are given
        if xtitle is not None:
            self.set_margin(bottom=0.23)
        if ytitle is not None:
            self.set_margin(left=0.18)
        # check xlog ylog
        if xlog is not None:
            if type(xlog) is not int:
                xlog = int(xlog)
        else:
            xlog = int(plot_panel_ax.get_xscale()=='log')
        if ylog is not None:
            if type(ylog) is not int:
                ylog = int(ylog)
        else:
            ylog = int(plot_panel_ax.get_yscale()=='log')
        # parse symbol (does not override marker)
        if symbol is not None and marker is None:
            if symbol == 'square':
                marker = 's'
            elif symbol == 'cross':
                marker = 'x'
            elif symbol == 'Cross':
                marker = 'X'
            elif symbol == 'diamond':
                marker = 'd'
            elif symbol == 'Diamond':
                marker = 'D'
            elif symbol == 'open square' or symbol == 'open squares' or symbol == 'empty square' or symbol == 'empty squares':
                marker = 's'
                facecolor = 'none'
                if color is not None:
                    edgecolor = color
            elif symbol == 'open circle' or symbol == 'open circles' or symbol == 'empty circle' or symbol == 'empty circles':
                marker = 'o'
                facecolor = 'none'
                if color is not None:
                    edgecolor = color
            elif symbol == 'filled square' or symbol == 'filled squares' or symbol == 'solid square' or symbol == 'solid squares':
                marker = 's'
            elif symbol == 'filled circle' or symbol == 'filled circles' or symbol == 'solid circle' or symbol == 'solid circles':
                marker = 'o'
            elif symbol == 'upper limit' or symbol == 'upper limits':
                #marker = u'$\u2193$'
                marker = self.get_upper_limit_marker()
                if color is not None:
                    edgecolor = color
                facecolor = 'none'
                if symsize is not None:
                    symsize = symsize * 3
                #uplims = True
                if yerr is not None:
                    yerr = None
            else:
                marker = symbol
        # parse symsize (does not override size)
        if symsize is not None and size is None:
            size = symsize*3
        # parse thick
        if thick is not None:
            if edgewidth is None:
                edgewidth = thick
            if linewidth is None:
                linewidth = thick
        # plot data
        if uplims is None:
            if xlog>0 and ylog>0:
                if fmt:
                    plot_panel_ax.loglog(x, y, fmt, label=label, zorder=zorder)
                else:
                    plot_panel_ax.loglog(x, y, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha, label=label, zorder=zorder)
            elif xlog>0 and ylog<=0:
                plot_panel_ax.semilogx(x, y, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha, label=label, zorder=zorder)
            elif xlog<=0 and ylog>0:
                plot_panel_ax.semilogy(x, y, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha, label=label, zorder=zorder)
            else:
                plot_panel_ax.plot(x, y, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha, label=label, zorder=zorder)
            #plot_panel_xy['panel'].scatter(x, y, marker=marker, s=size**2, c=color, edgecolors=edgecolor, linewidths=linewidth, alpha=alpha)
        # plot errorbars
        if xerr is not None or yerr is not None:
            plot_var_errorbar = plot_panel_ax.errorbar(x, y, xerr=xerr, yerr=yerr, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, drawstyle=drawstyle, linewidth=linewidth, alpha=alpha, label=label, zorder=zorder, **kwargs)
            # set errorbar line style if given by the user
            if errorbarlinestyle != '':
                plot_var_errorbar[-1][0].set_linestyle(errorbarlinestyle) # [-1][0] is the LineCollection objects of the errorbar lines -- https://stackoverflow.com/questions/22995797/can-matplotlib-errorbars-have-a-linestyle-set?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
        #if uplims is not None:
        #    plot_panel_ax.errorbar(x, y, xerr=xerr, yerr=0.8*y, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, drawstyle=drawstyle, linewidth=linewidth, alpha=alpha, zorder=zorder, uplims=uplims, lolims=lolims, **kwargs)
        # axis ticks format
        # -- http://matplotlib.org/examples/ticks_and_spines/tick-locators.html
        # -- http://stackoverflow.com/questions/33126126/matplotlib-minor-ticks
        if xlog>0:
            #plot_panel_xy['panel'].xaxis.set_tick_params(which='major', length=8.0)
            #plot_panel_xy['panel'].xaxis.set_tick_params(which='minor', length=4.0)
            plot_panel_ax.xaxis.set_major_locator(LogLocator(base=10,numticks=30))
            plot_panel_ax.xaxis.set_minor_locator(LogLocator(base=10,numticks=30,subs=numpy.arange(2.0,10.0,1.0)))
            plot_panel_ax.xaxis.set_major_formatter(ticker.FuncFormatter(self.major_formatter_function_for_logarithmic_axis))
            #plot_panel_ax.xaxis.set_minor_formatter(NullFormatter())
            #plot_panel_ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            #plot_panel_ax.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
        else:
            plot_panel_ax.ticklabel_format(style='plain',axis='x',useOffset=False) # 20180122
        if ylog>0:
            #plot_panel_xy['panel'].yaxis.set_tick_params(which='major', length=8.0)
            #plot_panel_xy['panel'].yaxis.set_tick_params(which='minor', length=4.0)
            plot_panel_ax.yaxis.set_major_locator(LogLocator(base=10,numticks=30))
            plot_panel_ax.yaxis.set_minor_locator(LogLocator(base=10,numticks=30,subs=numpy.arange(2.0,10.0,1.0)))
            plot_panel_ax.yaxis.set_major_formatter(ticker.FuncFormatter(self.major_formatter_function_for_logarithmic_axis))
            #plot_panel_ax.yaxis.set_minor_formatter(NullFormatter())
            #plot_panel_ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            #plot_panel_ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
            ##plot_panel_xy['panel'].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ##plot_panel_xy['panel'].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        else:
            plot_panel_ax.ticklabel_format(style='plain',axis='y',useOffset=False) # 20180122
        # 
        #print('CrabPlot::plot_xy() debug xlog=%s ylog=%s'%(xlog, ylog)) # 20180122, axis useOffset problem.
        #print('CrabPlot::plot_xy() debug xlog=%s ylog=%s'%(xlog, ylog)) # 20180122, axis useOffset problem.
        #print('CrabPlot::plot_xy() debug xlog=%s ylog=%s'%(xlog, ylog)) # 20180122, axis useOffset problem.
        # 
        # set tick label format, scitific notation
        #plot_panel_ax.ticklabel_format(axis='both', style='sci', scilimits=(-2,+2)) # not working for log!
        # 
        # set title, range, etc.
        if xtitle is not None:
            if xtitle != '':
                #plot_panel_xy['panel'].set_xlabel(xtitle, fontsize=14)
                self.set_xtitle(xtitle, ax=plot_panel_ax)
                if current>0: self.Plot_panels[current-1]['xtitle'] = xtitle
        if ytitle is not None:
            if ytitle != '':
                #plot_panel_xy['panel'].set_ylabel(ytitle, fontsize=14)
                self.set_ytitle(ytitle, ax=plot_panel_ax)
                if current>0: self.Plot_panels[current-1]['ytitle'] = ytitle
        if xrange is not None:
            if len(xrange) >= 2:
                plot_panel_ax.set_xlim(xrange)
                if current>0: self.Plot_panels[current-1]['xrange'] = xrange
        if yrange is not None:
            if len(yrange) >= 2:
                plot_panel_ax.set_ylim(yrange)
                if current>0: self.Plot_panels[current-1]['yrange'] = yrange
        # store data
        if current>0:
            self.Plot_panels[current-1]['x'] = x
            self.Plot_panels[current-1]['y'] = y
            self.Plot_panels[current-1]['xerr'] = xerr
            self.Plot_panels[current-1]['yerr'] = yerr
            self.Plot_panels[current-1]['xlog'] = xlog
            self.Plot_panels[current-1]['ylog'] = ylog
    # 
    @ticker.FuncFormatter
    def major_formatter_function_for_logarithmic_axis(x, pos):
        if x > 1000 or x < 0.001:
            return '10$^{%g}$'%(numpy.log10(x))
        else:
            return '%g'%(x)
    # 
    def plot_line(self, x0, y0, x1 = None, y1 = None, xlog = None, ylog = None, xrange = [], yrange = [], 
                    ax = None, NormalizedCoordinate = False, overplot = True, current = 0, 
                    position = None, label = None, 
                    dataname = '', 
                    xtitle = None, ytitle = None, 
                    xtitlefontsize = 14, ytitlefontsize = 14, 
                    linestyle = 'solid', 
                    margin = None, 
                    thick = 0, linewidth = 0, 
                    **kwargs):
        # get panel ax
        plot_panel_ax, current, overplot = self.get_panel_ax(ax, current, overplot)
        # set linewidth and thick
        if thick > 0:
            if linewidth <= 0:
                linewidth = thick
        if linewidth <= 0:
            linewidth = 1.0 # default linewidth
        # check x1 y1
        # we can input either (x0,y0), or (x0,y0,x1,y1).
        if x1 is not None and y1 is not None:
            # check x y type must be scalar in this case
            x0 = self.check_scalar(x0)
            y0 = self.check_scalar(y0)
            x1 = self.check_scalar(x1)
            y1 = self.check_scalar(y1)
            # plot line
            if plot_panel_ax:
                if NormalizedCoordinate is True:
                    #xrange_for_plot = plot_panel_ax.get_xlim()
                    #yrange_for_plot = plot_panel_ax.get_ylim()
                    #x0_for_plot = x0 * (xrange_for_plot[1]-xrange_for_plot[0]) + xrange_for_plot[0]
                    #y0_for_plot = y0 * (yrange_for_plot[1]-yrange_for_plot[0]) + yrange_for_plot[0]
                    plot_one_line = matplotlib.lines.Line2D([x0,x1], [y0,y1], linestyle=linestyle, linewidth=linewidth, transform=plot_panel_ax.transAxes, label=label, **kwargs)
                else:
                    plot_one_line = matplotlib.lines.Line2D([x0,x1], [y0,y1], linestyle=linestyle, linewidth=linewidth, label=label, **kwargs)
                # 
                plot_panel_ax.add_line(plot_one_line)
        else:
            # when x1 and y1 are None, then we assume x0 and y0 are arrays
            if x0 is not list and x0 is not numpy.array and x0 is not numpy.ndarray:
                x0 = numpy.array(x0)
            if y0 is not list and y0 is not numpy.array and y0 is not numpy.ndarray:
                y0 = numpy.array(y0)
            if x0.shape != y0.shape:
                print('Error! plot_line x0 y0 do not have the same dimension!')
                return
            # check non-positive values if ylog
            if ylog is not None:
                ylog = int(ylog)
                if ylog > 0:
                    if numpy.count_nonzero(y0<=0.0):
                        y0[y0<=0.0] = numpy.nan
            # store x y
            if dataname != '':
                self.Plot_data[dataname] = numpy.column_stack((x0,y0))
                #print('CrabPlot::plot_line() Stored the input x and y into self.Plot_data[%s]!'%(dataname))
            # sort
            x_sorted_index = x0.argsort()
            x0 = x0[x_sorted_index]
            y0 = y0[x_sorted_index]
            # plot line
            if plot_panel_ax:
                if NormalizedCoordinate is True:
                    #xrange_for_plot = plot_panel_ax.get_xlim()
                    #yrange_for_plot = plot_panel_ax.get_ylim()
                    #x0_for_plot = x0 * (xrange_for_plot[1]-xrange_for_plot[0]) + xrange_for_plot[0]
                    #y0_for_plot = y0 * (yrange_for_plot[1]-yrange_for_plot[0]) + yrange_for_plot[0]
                    plot_one_line = matplotlib.lines.Line2D(x0, y0, linestyle=linestyle, linewidth=linewidth, transform=plot_panel_ax.transAxes, label=label, **kwargs)
                else:
                    plot_one_line = matplotlib.lines.Line2D(x0, y0, linestyle=linestyle, linewidth=linewidth, label=label, **kwargs)
                # 
                plot_panel_ax.add_line(plot_one_line)
        # 
        # set log, title, range, etc.
        if xlog is not None:
            xlog = int(xlog)
            if xlog > 0:
                plot_panel_ax.set_xscale('log')
                plot_panel_ax.xaxis.set_major_formatter(ticker.FuncFormatter(self.major_formatter_function_for_logarithmic_axis))
                if current>0: self.Plot_panels[current-1]['xlog'] = xlog
        if ylog is not None:
            ylog = int(ylog)
            if ylog > 0:
                plot_panel_ax.set_yscale('log')
                plot_panel_ax.yaxis.set_major_formatter(ticker.FuncFormatter(self.major_formatter_function_for_logarithmic_axis))
                if current>0: self.Plot_panels[current-1]['ylog'] = ylog
        if xtitle is not None:
            if xtitle != '':
                #plot_panel_ax.set_xlabel(xtitle, fontsize=xtitlefontsize)
                self.set_xtitle(xtitle, ax=plot_panel_ax)
                if current>0: self.Plot_panels[current-1]['xtitle'] = xtitle
        if ytitle is not None:
            if ytitle != '':
                #plot_panel_ax.set_ylabel(ytitle, fontsize=ytitlefontsize)
                self.set_ytitle(ytitle, ax=plot_panel_ax)
                if current>0: self.Plot_panels[current-1]['ytitle'] = ytitle
        if xrange is not None:
            if len(xrange) >= 2:
                plot_panel_ax.set_xlim(xrange)
                if current>0: self.Plot_panels[current-1]['xrange'] = xrange
        if yrange is not None:
            if len(yrange) >= 2:
                plot_panel_ax.set_ylim(yrange)
                if current>0: self.Plot_panels[current-1]['yrange'] = yrange
    # 
    def plot_text(self, x0, y0, text_input, ax = None, current = None, overplot = None, NormalizedCoordinate = False, **kwargs):
        # check x y type
        x0 = self.check_scalar(x0)
        y0 = self.check_scalar(y0)
        # check if a direct matplotlib axis var is given or not
        #if ax is not None:
        #    current = -1
        #    overplot = True
        #    plot_panel_ax = ax
        #    plot_panel_xy = None
        #else:
        #    # get current panel or add new panel
        #    if len(self.Plot_panels) == 0:
        #        plot_panel_xy = self.add_panel(position = position, label = label)
        #        plot_panel_ax = plot_panel_xy['panel']
        #        current = len(self.Plot_panels)
        #        current = 1
        #        overplot = False
        #    # get current panel
        #    elif current>0:
        #        plot_panel_xy = self.Plot_panels[current-1]
        #        plot_panel_ax = plot_panel_xy['panel']
        #    # get last panel for overplotting
        #    elif overplot: 
        #        plot_panel_xy = self.Plot_panels[-1]
        #        plot_panel_ax = plot_panel_xy['panel']
        #        current = len(self.Plot_panels)
        ## get plot panel 'ax' variable
        #if ax is not None:
        #    plot_panel_ax = ax
        #else:
        #    plot_panel_xy = self.Plot_panels[-1]
        #    plot_panel_ax = plot_panel_xy['panel']
        ax, current, overplot = self.get_panel_ax(ax, current, overplot)
        # plot line
        if ax:
            if NormalizedCoordinate is True:
                ax.text(x0, y0, text_input, transform=ax.transAxes, **kwargs) # verticalalignment='center', horizontalalignment='left'
            else:
                ax.text(x0, y0, text_input, **kwargs) # verticalalignment='center', horizontalalignment='left'
    # 
    def xyouts(self, x0, y0, text_input, ax = None, current = None, overplot = True, Data = False, NormalizedCoordinate = True, **kwargs):
        if Data is True:
            NormalizedCoordinate = False
        ax, current, overplot = self.get_panel_ax(ax, current, overplot)
        self.plot_text(x0, y0, text_input, ax = ax, NormalizedCoordinate = NormalizedCoordinate, **kwargs)
    # 
    def plot_hist(self, x, y, ax = None, current = None, overplot = True, xtitle = None, ytitle = None, xlog = None, ylog = None, useTex = None, 
                    margin = None, xrange = None, yrange = None, 
                    **kwargs):
        # get panel ax
        plot_panel_ax, current, overplot = self.get_panel_ax(ax, current, overplot)
        # set margin if x y title are given
        if xtitle is not None:
            self.set_margin(bottom=0.2)
        if ytitle is not None:
            self.set_margin(left=0.15)
        # set log
        log = False
        if xlog is not None:
            xlog = int(xlog)
            if xlog > 0:
                plot_panel_ax.set_xscale('log')
                plot_panel_ax.xaxis.set_major_formatter(ticker.FuncFormatter(self.major_formatter_function_for_logarithmic_axis))
                if current>0: self.Plot_panels[current-1]['xlog'] = xlog
                #<TODO># how to apply xlog?
            else:
                plot_panel_ax.ticklabel_format(style='plain',axis='x',useOffset=False) # 20180122
        else:
            plot_panel_ax.ticklabel_format(style='plain',axis='x',useOffset=False) # 20180122
        if ylog is not None:
            ylog = int(ylog)
            if ylog > 0:
                plot_panel_ax.set_yscale('log')
                plot_panel_ax.yaxis.set_major_formatter(ticker.FuncFormatter(self.major_formatter_function_for_logarithmic_axis))
                if current>0: self.Plot_panels[current-1]['ylog'] = ylog
                log = True
            else:
                plot_panel_ax.ticklabel_format(style='plain',axis='x',useOffset=False) # 20180122
        else:
            plot_panel_ax.ticklabel_format(style='plain',axis='x',useOffset=False) # 20180122
        # plot histogram using the matplotlib bar() function
        plot_panel_ax.bar(x, y, log=log, **kwargs)
        # set titles
        if xtitle is not None:
            if xtitle != '':
                #plot_panel_ax.set_xlabel(xtitle, fontsize=xtitlefontsize)
                self.set_xtitle(xtitle, ax=plot_panel_ax)
                if current>0: self.Plot_panels[current-1]['xtitle'] = xtitle
        if ytitle is not None:
            if ytitle != '':
                #plot_panel_ax.set_ylabel(ytitle, fontsize=ytitlefontsize)
                self.set_ytitle(ytitle, ax=plot_panel_ax)
                if current>0: self.Plot_panels[current-1]['ytitle'] = ytitle
        # set ranges
        if xrange is not None:
            if len(xrange) >= 2:
                plot_panel_ax.set_xlim(xrange)
                if current>0: self.Plot_panels[current-1]['xrange'] = xrange
        if yrange is not None:
            if len(yrange) >= 2:
                plot_panel_ax.set_ylim(yrange)
                if current>0: self.Plot_panels[current-1]['yrange'] = yrange
    # 
    def plot_image(self, image_data, image_wcs = None, position = None, label = None):
        # check image_data type
        image_data = self.check_array(image_data)
        # check image_data dimension
        if len(image_data.shape) != 2:
            print('CrabPlot::plot_image() Error! The input image_data does not have a shape of 2D!')
            return
        # add panel
        plot_panel_im = self.add_panel(position = position, label = label, projection = image_wcs)
        # plot image
        background_mu = numpy.nanmean(self.Image_data)
        background_sigma = numpy.nanstd(self.Image_data)
        plot_panel_im['panel'].imshow(self.Image_data, origin = 'lower', cmap = 'binary', aspect = 'equal', 
                                      norm = ImageNormalize(self.Image_data, vmin=background_mu-0.5*background_sigma, vmax=background_mu+3.0*background_sigma))
        # store data
        plot_panel_im['x'] = x
        plot_panel_im['y'] = y
        plot_panel_im['xerr'] = xerr
        plot_panel_im['yerr'] = yerr
        plot_panel_im['xlog'] = xlog
        plot_panel_im['ylog'] = ylog
    # 
    def plot_data_file(self, input_data_file, xcol = 1, ycol = 2, redshift = None, linestyle = 'none', 
                            xclip = None, yclip = None, **kwargs):
        # we also allow 'xclip', which is a list of tuples, and we clip the data by filtering out where x in each xclip tuple/range.
        #print(input_data_file)
        input_data_table = asciitable.read(input_data_file, Reader=asciitable.NoHeader, delimiter=' ', guess=False)
        if input_data_table:
            if len(input_data_table.colnames) >= numpy.max(numpy.array([xcol,ycol])):
                input_x = input_data_table.field(input_data_table.colnames[xcol-1])
                input_y = input_data_table.field(input_data_table.colnames[ycol-1])
                if redshift is not None:
                    input_x = input_x * (1.0+float(redshift))
                if xclip is not None:
                    if len(xclip) > 0:
                        for clip_x in xclip:
                            if len(clip_x) >= 2:
                                #print('CrabPlot::plot_data_file() Clipping data from %s to %s'%(clip_x[0],clip_x[1]))
                                clip_mask = (input_x >= clip_x[0]) & (input_x < clip_x[1])
                                clip_args = numpy.argwhere(clip_mask)
                                if len(clip_args) > 0:
                                    input_x = numpy.delete(input_x, clip_args)
                                    input_y = numpy.delete(input_y, clip_args)
                if linestyle == 'none':
                    #print(input_x, input_y)
                    self.plot_xy(input_x, input_y, **kwargs)
                else:
                    self.plot_line(input_x, input_y, linestyle = linestyle, **kwargs)
    # 
    def spline(self, input_x, input_y, output_x, xlog=0, ylog=0, outputxlog=None, outputylog=None, **kwargs):
        # spline
        # note that we can input xlog, ylog, outputxlog, outputylog to deal with logarithm input/output needs
        # xlog>0 means the input_x array is in linear space and will be converted to log space when doing the spline
        # ylog>0 means the input_y array is in linear space and will be converted to log space when doing the spline
        # outputxlog>0 means the output_x array is in linear space and will be converted to log space when doing the spline
        # outputylog>0 means the output_y array is in linear space but is in log space when doing the spline, so we need to convert it back to linear space after the spline.
        input_x_coord = self.check_array(input_x)
        input_y_value = self.check_array(input_y)
        output_x_coord = self.check_array(output_x)
        if xlog>0:
            input_x_mask = (input_x_coord<=0.0)
            input_x_coord[input_x_mask] = numpy.nan
            input_x_coord = numpy.log10(input_x)
        if ylog>0:
            input_y_mask = (input_y_value<=0.0)
            input_y_value[input_y_mask] = numpy.nan
            input_y_value = numpy.log10(input_y_value)
        # 
        if outputxlog is None:
            outputxlog = xlog
        if outputylog is None:
            outputylog = ylog
        # 
        if outputxlog>0:
            output_x_mask = (output_x_coord<=0.0)
            output_x_coord[output_x_mask] = numpy.nan
            output_x_coord = numpy.log10(output_x_coord)
        # 
        #print(numpy.column_stack((input_x_coord, input_y_value)))
        input_mask = numpy.isnan(input_y_value)
        ###spl = UnivariateSpline(input_x_coord, input_y_value, w=~input_mask)
        ###output_y_value = spl(output_x_coord)
        ##output_y_value = scipy.interpolate.spline(input_x_coord, input_y_value, output_x_coord, **kwargs) # order=3, kind='smoothest', conds=None
        spl = scipy.interpolate.CubicSpline(input_x_coord, input_y_value, **kwargs) # axis=0, bc_type='not-a-knot', extrapolate=None
        output_y_value = spl(output_x_coord)
        # 
        # deal with data out of X range
        output_mask = (output_x_coord<numpy.nanmin(input_x_coord)) | (output_x_coord>numpy.nanmax(input_x_coord))
        output_y_value[output_mask] = numpy.nan
        # 
        if outputylog>0:
            output_y = numpy.power(10,output_y_value)
            output_y[output_mask] = 0.0 #<TODO># fill with 0.0 instead of nan
        else:
            output_y = output_y_value
        # 
        return output_y
    # 
    def ax(self, ax = None, current = True, overplot = True):
        plot_panel_ax, current, overplot = self.get_panel_ax(ax, current, overplot)
        return plot_panel_ax
    # 
    def show(self, block=True):
        self.Plot_device.canvas.draw()
        self.Plot_device.show()
        pyplot.show(block=block)
        #self.Plot_device.waitforbuttonpress()
    # 
    def savefig(self, filename):
        self.Plot_device.savefig(filename)
        print('Output to "%s"!'%(filename))
    # 
    def savepdf(self, filename):
        if not filename.endswith('.pdf') and not filename.endswith('.PDF'):
            filename = filename + '.pdf'
        self.Plot_device.savefig(filename)
        print('Output to "%s"!'%(filename))
    # 
    def save(self, filename):
        if filename.endswith('.png') or filename.endswith('.PNG'):
            self.savefig(filename)
        elif filename.endswith('.jpg') or filename.endswith('.JPG'):
            self.savefig(filename)
        elif filename.endswith('.jpeg') or filename.endswith('.JPEG'):
            self.savefig(filename)
        else:
            self.savepdf(filename)
        #self.Plot_device.savefig(filename)
    # 
    def close(self):
        self.clear()
    # 
    def clear(self):
        if self.Plot_device is not None:
            self.Plot_device.clf() # fig.clear() is a synonym for fig.clf()
            pyplot.close(self.Plot_device)










# 
def plot_line(ax, x0, y0, x1 = None, y1 = None, NormalizedCoordinate = False, **kwargs):
    if x1 is not None and y1 is not None:
        if NormalizedCoordinate is True:
            plot_one_line = matplotlib.lines.Line2D([x0,x1], [y0,y1], transform=ax.transAxes, **kwargs)
        else:
            plot_one_line = matplotlib.lines.Line2D([x0,x1], [y0,y1], **kwargs)
        # 
        ax.add_line(plot_one_line)
    else:
        # when x1 and y1 are None, then we assume x0 and y0 are arrays
        if x0 is not list and x0 is not numpy.array and x0 is not numpy.ndarray:
            x0 = numpy.array(x0)
        if y0 is not list and y0 is not numpy.array and y0 is not numpy.ndarray:
            y0 = numpy.array(y0)
        if x0.shape != y0.shape:
            print('Error! plot_line x0 y0 do not have the same dimension!')
            return
        # sort
        x_sorted_index = x0.argsort()
        x0 = x0[x_sorted_index]
        y0 = y0[x_sorted_index]
        # 
        if NormalizedCoordinate is True:
            plot_one_line = matplotlib.lines.Line2D(x0, y0, transform=ax.transAxes, **kwargs)
        else:
            plot_one_line = matplotlib.lines.Line2D(x0, y0, **kwargs)
        # 
        ax.add_line(plot_one_line)
        return plot_one_line


# 
def plot_text(ax, x0, y0, text_input, NormalizedCoordinate = False, **kwargs):
    if NormalizedCoordinate is True:
        ax.text(x0, y0, text_input, transform=ax.transAxes, **kwargs) # verticalalignment='center', horizontalalignment='left'
    else:
        ax.text(x0, y0, text_input, **kwargs) # verticalalignment='center', horizontalalignment='left'







def adjust_points_not_too_close(x0, y0, mask = None, too_close_distance = 0.0, max_adjust_distance = 1.0): 
    x0_adjusted = copy(x0)
    y0_adjusted = copy(y0)
    i0 = numpy.array(range(len(x0)))
    if mask is None:
        mask = (i0>=0)
    for i in i0:
        # check nan and mask
        if (x0[i]==x0[i]) and (y0[i]==y0[i]) and (mask[i]==True):
            # choose a 'too_close_distance'
            if too_close_distance <= 0.0:
                imask = (i0[mask]>=0)
            else:
                imask = (numpy.abs(x0[mask]-x0[i])<=too_close_distance) & (numpy.abs(y0[mask]-y0[i])<=too_close_distance)
            # 
            iargs = numpy.argwhere(imask)
            if(len(iargs)>0):
                x1dis = (x0[mask]-x0[i])
                y1dis = (y0[mask]-y0[i])
                x1wei = numpy.exp(-(x1dis**2)/(2.0*(too_close_distance/(2*numpy.sqrt(2*numpy.log(2))))**2)) # 1/(x1dis[imask])**2 # 
                y1wei = numpy.exp(-(y1dis**2)/(2.0*(too_close_distance/(2*numpy.sqrt(2*numpy.log(2))))**2)) # 1/(y1dis[imask])**2 # 
                x1cen = numpy.nansum(x1wei[imask]*x0[mask][imask])/numpy.nansum(x1wei[imask])
                y1cen = numpy.nansum(y1wei[imask]*y0[mask][imask])/numpy.nansum(y1wei[imask])
                #print(x0[i], y0[i], x1cen, y1cen)
                x0_adjusted[i] = x0[i]+(x0[i]-x1cen)*(numpy.nansum(x1wei[imask])-1)*max_adjust_distance
                y0_adjusted[i] = y0[i]+(y0[i]-y1cen)*(numpy.nansum(x1wei[imask])-1)*max_adjust_distance
                #print(x0[i], y0[i], x1cen, y1cen, x0_adjusted[i], y0_adjusted[i])
                #if (i+1)==102:
                #    print(numpy.column_stack((i0[mask][imask]+1, x0[mask][imask], y0[mask][imask], x1wei[imask], y1wei[imask])))
                #    print(i+1, x0[i], y0[i], x1cen, y1cen, x0_adjusted[i], y0_adjusted[i])
    return x0_adjusted, y0_adjusted















