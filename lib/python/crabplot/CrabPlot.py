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
pkg_resources.require("wcsaxes") # http://wcsaxes.readthedocs.io/en/latest/getting_started.html

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
        matplotlib.use('Qt5Agg')
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

from matplotlib.ticker import FormatStrFormatter, NullFormatter, AutoMinorLocator, LogLocator



import warnings

warnings.filterwarnings("ignore",".*GUI is implemented.*")













# 
class CrabPlot(object):
    # 
    def __init__(self, x = None, y = None, xerr = None, yerr = None, xlog = False, ylog = False, xtitle = None, ytitle = None, xrange = [], yrange = [], 
                       image_data = None, image_wcs = None, figure_size = (9.0,8.0), figure_dpi = 90, position = None, label = ''):
        # 
        # open plot
        self.Plot_device = pyplot.figure(figsize=figure_size, dpi=figure_dpi) # set figure size 9.0 x 8.0 inches, 90 pixels per inch. 
        self.Plot_panels = [] # each item is a dict {'label': '', 'panel': axis, 'position': position}
        self.Plot_grids = None
        self.Plot_data = {}
        self.Annotation = []
        self.Image_data = image_data
        self.Image_wcs = image_wcs
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
        # 
        nrows, ncols = self.Plot_grids.get_geometry()
        #nrows = long(nrows)
        #ncols = long(ncols)
        nrows = int(nrows) # Python 3 has no 'long' type anymore
        ncols = int(ncols) # Python 3 has no 'long' type anymore
        
        subplot_params = self.Plot_grids.get_subplot_params(self.Plot_device)
        left = subplot_params.left
        right = subplot_params.right
        bottom = subplot_params.bottom
        top = subplot_params.top
        wspace = subplot_params.wspace
        hspace = subplot_params.hspace
        totWidth = right-left
        totHeight = top-bottom
        
        # calculate accumulated heights of columns
        cellH = totHeight/(nrows + hspace*(nrows-1))
        sepH = hspace*cellH
        
        if self.Plot_grids._row_height_ratios is not None:
            netHeight = cellH * nrows
            tr = float(sum(self.Plot_grids._row_height_ratios))
            cellHeights = [netHeight*r/tr for r in self.Plot_grids._row_height_ratios]
        else:
            cellHeights = [cellH] * nrows
            
        sepHeights = [0] + ([sepH] * (nrows-1))
        cellHs = numpy.add.accumulate(numpy.ravel(list(zip(sepHeights, cellHeights))))
        
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
        
        figTops = [top - cellHs[2*rowNum] for rowNum in range(nrows)]
        figBottoms = [top - cellHs[2*rowNum+1] for rowNum in range(nrows)]
        figLefts = [left + cellWs[2*colNum] for colNum in range(ncols)]
        figRights = [left + cellWs[2*colNum+1] for colNum in range(ncols)]
        
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
    def add_panel(self, position = None, label = None, projection = None, debug = False):
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
                    self.Plot_panels[i]['panel'].set_position(self.get_panel_position_in_gridspec(i)) # see -- http://stackoverflow.com/questions/22881301/changing-matplotlib-subplot-size-position-after-axes-creation
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
                                   'x': None, 'y': None, 'xerr': None, 'yerr': None, 'xlog': False, 'ylog': False, 'xrange': [], 'yrange': [], 
                                   'image_data': None, 'image_wcs': None } )
        # set ticks
        self.set_xticksize(panel=len(self.Plot_panels))
        self.set_yticksize(panel=len(self.Plot_panels))
        # return
        return self.Plot_panels[-1]
    # 
    def default_fontname_for_axis_title(self):
        return ['NGC','Helvetica','sans-serif']
    # 
    def default_fontsize_for_axis_title(self):
        return 14
    # 
    def default_fontname_for_axis_ticks(self):
        return ['NGC','Helvetica','sans-serif']
    # 
    def default_fontsize_for_axis_ticks(self):
        return 13
    # 
    # def set_xrange
    def set_xrange(self, xrange, ax=None, panel=None):
        if len(xrange)<2:
            return
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['xrange'] = xrange
                elif panel > 0 and panel <= len(self.Plot_panels):
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
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['yrange'] = yrange
                elif panel > 0 and panel <= len(self.Plot_panels):
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
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['xtitle'] = xtitle
                elif panel > 0 and panel <= len(self.Plot_panels):
                    ax = self.Plot_panels[panel-1]['panel']
                    self.Plot_panels[panel-1]['xtitle'] = xtitle
                else:
                    return
            if fontsize is None:
                fontsize = self.default_fontsize_for_axis_title()
            if fontname is None:
                fontname = self.default_fontname_for_axis_title()
            if ax is not None:
                ax.set_xlabel(xtitle, fontsize=fontsize, fontname=fontname, **kwargs)
                self.set_xticksfont(fontsize=fontsize-1)
    # 
    # def set_ytitle
    def set_ytitle(self, ytitle, ax=None, panel=None, fontsize=None, fontname=None, **kwargs):
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
                    ax = self.Plot_panels[panel-1]['panel']
                else:
                    return
            if fontsize is None:
                fontsize = self.default_fontsize_for_axis_title()
            if fontname is None:
                fontname = self.default_fontname_for_axis_title()
            if ax is not None:
                ax.set_ylabel(ytitle, fontsize=fontsize, fontname=fontname, **kwargs)
                self.set_yticksfont(fontsize=fontsize-1)
    # 
    # def set_xtickfont
    def set_xtickfont(self, ax=None, panel=None, fontsize=None, fontname=None):
        if len(self.Plot_panels)>0:
            if ax is None:
                if panel is None:
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
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
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
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
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
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
                    panel = len(self.Plot_panels)
                    ax = self.Plot_panels[panel-1]['panel']
                elif panel > 0 and panel <= len(self.Plot_panels):
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
    def plot_xy(self, x, y, xerr = None, yerr = None, xlog = None, ylog = None, xrange = [], yrange = [], 
                ax = None, NormalizedCoordinate = False, overplot = False, current = 0, 
                position = None, label = None, 
                dataname = '', 
                xtitle = None, ytitle = None, 
                fmt = None, 
                symbol = 'o', symsize = 3, thick = 2, 
                marker = None, size = None, color = 'blue', fillstyle = None, 
                linestyle = 'None', linewidth = None, drawstyle = None, 
                facecolor = None, edgecolor = None, edgewidth = None, alpha = 1.0, zorder = None, 
                uplims = None, lolims = None, 
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
            print('CrabPlot::plot_xy() Stored the input x and y into self.Plot_data[%s]!'%(dataname))
        # check x y dimension
        if x.shape != y.shape:
            print('CrabPlot::plot_xy() Error! The input x and y do not have the same shape!')
            return
        # get panel ax
        plot_panel_ax, current, overplot = self.get_panel_ax(ax, current, overplot)
        # set grid
        plot_panel_ax.grid(False)
        # check xlog ylog
        if xlog is not None:
            if type(xlog) is bool:
                xlog = int(xlog)
        else:
            xlog = (plot_panel_ax.get_xscale()=='log')
        if ylog is not None:
            if type(ylog) is bool:
                ylog = int(ylog)
        else:
            ylog = (plot_panel_ax.get_yscale()=='log')
        # parse symbol (does not override marker)
        if symbol is not None and marker is None:
            if symbol == 'square':
                marker = 's'
            elif symbol == 'open square' or symbol == 'open squares':
                marker = 's'
                facecolor = 'none'
                if color is not None:
                    edgecolor = color
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
                    plot_panel_ax.loglog(x, y, fmt, zorder=zorder)
                else:
                    plot_panel_ax.loglog(x, y, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha, zorder=zorder)
            elif xlog>0 and ylog<=0:
                plot_panel_ax.semilogx(x, y, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha, zorder=zorder)
            elif xlog<=0 and ylog>0:
                plot_panel_ax.semilogy(x, y, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha, zorder=zorder)
            else:
                plot_panel_ax.plot(x, y, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha, zorder=zorder)
            #plot_panel_xy['panel'].scatter(x, y, marker=marker, s=size**2, c=color, edgecolors=edgecolor, linewidths=linewidth, alpha=alpha)
        # plot errorbars
        if xerr is not None or yerr is not None:
            plot_panel_ax.errorbar(x, y, xerr=xerr, yerr=yerr, marker=marker, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, drawstyle=drawstyle, linewidth=linewidth, alpha=alpha, zorder=zorder, **kwargs)
        #if uplims is not None:
        #    plot_panel_ax.errorbar(x, y, xerr=xerr, yerr=0.8*y, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, drawstyle=drawstyle, linewidth=linewidth, alpha=alpha, zorder=zorder, uplims=uplims, lolims=lolims, **kwargs)
        # axis ticks format
        # -- http://matplotlib.org/examples/ticks_and_spines/tick-locators.html
        # -- http://stackoverflow.com/questions/33126126/matplotlib-minor-ticks
        #if xlog>0:
        #    plot_panel_xy['panel'].xaxis.set_tick_params(which='major', length=8.0)
        #    plot_panel_xy['panel'].xaxis.set_tick_params(which='minor', length=4.0)
        if ylog>0:
            #plot_panel_xy['panel'].yaxis.set_tick_params(which='major', length=8.0)
            #plot_panel_xy['panel'].yaxis.set_tick_params(which='minor', length=4.0)
            plot_panel_ax.yaxis.set_major_locator(LogLocator(base=10,numticks=30))
            plot_panel_ax.yaxis.set_minor_locator(LogLocator(base=10,numticks=30,subs=numpy.arange(2.0,10.0,1.0)))
            plot_panel_ax.yaxis.set_minor_formatter(NullFormatter())
            ##print(plot_panel_xy['panel'].yaxis.get_major_formatter())
            ##plot_panel_xy['panel'].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ##plot_panel_xy['panel'].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
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
    def plot_line(self, x0, y0, x1 = None, y1 = None, xlog = None, ylog = None, xrange = [], yrange = [], 
                    ax = None, NormalizedCoordinate = False, overplot = True, current = 0, 
                    position = None, label = None, 
                    dataname = '', 
                    xtitle = None, ytitle = None, 
                    xtitlefontsize = 14, ytitlefontsize = 14, 
                    linestyle = 'solid', 
                    **kwargs):
        # get panel ax
        plot_panel_ax, current, overplot = self.get_panel_ax(ax, current, overplot)
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
                    plot_one_line = matplotlib.lines.Line2D([x0,x1], [y0,y1], linestyle=linestyle, transform=plot_panel_ax.transAxes, **kwargs)
                else:
                    plot_one_line = matplotlib.lines.Line2D([x0,x1], [y0,y1], linestyle=linestyle, **kwargs)
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
                    plot_one_line = matplotlib.lines.Line2D(x0, y0, linestyle=linestyle, transform=plot_panel_ax.transAxes, **kwargs)
                else:
                    plot_one_line = matplotlib.lines.Line2D(x0, y0, linestyle=linestyle, **kwargs)
                # 
                plot_panel_ax.add_line(plot_one_line)
        # 
        # set log, title, range, etc.
        if xlog is not None:
            if xlog > 0:
                plot_panel_ax.set_xscale('log')
                if current>0: self.Plot_panels[current-1]['xlog'] = xlog
        if ylog is not None:
            if ylog > 0:
                plot_panel_ax.set_yscale('log')
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
        if ax is not None:
            current = -1
            overplot = True
            plot_panel_ax = ax
            plot_panel_xy = None
        else:
            # get current panel or add new panel
            if len(self.Plot_panels) == 0:
                plot_panel_xy = self.add_panel(position = position, label = label)
                plot_panel_ax = plot_panel_xy['panel']
                current = len(self.Plot_panels)
                current = 1
                overplot = False
            # get current panel
            elif current>0:
                plot_panel_xy = self.Plot_panels[current-1]
                plot_panel_ax = plot_panel_xy['panel']
            # get last panel for overplotting
            elif overplot: 
                plot_panel_xy = self.Plot_panels[-1]
                plot_panel_ax = plot_panel_xy['panel']
                current = len(self.Plot_panels)
        # get plot panel 'ax' variable
        if ax is not None:
            plot_panel_ax = ax
        else:
            plot_panel_xy = self.Plot_panels[-1]
            plot_panel_ax = plot_panel_xy['panel']
        # plot line
        if plot_panel_ax:
            if NormalizedCoordinate is True:
                plot_panel_ax.text(x0, y0, text_input, transform=plot_panel_ax.transAxes, **kwargs) # verticalalignment='center', horizontalalignment='left'
            else:
                plot_panel_ax.text(x0, y0, text_input, **kwargs) # verticalalignment='center', horizontalalignment='left'
    # 
    def xyouts(self, x0, y0, text_input, ax = None, current = None, overplot = True, NormalizedCoordinate = False, **kwargs):
        ax, current, overplot = self.get_panel_ax(ax, current, overplot)
        self.plot_text(x0, y0, text_input, ax = ax, NormalizedCoordinate = NormalizedCoordinate, **kwargs)
    # 
    def plot_hist(self, x, y, ax = None, current = None, overplot = True, xtitle = None, ytitle = None, useTex = None, **kwargs):
        plot_panel_ax, current, overplot = self.get_panel_ax(ax, current, overplot)
        plot_panel_ax.bar(x, y, **kwargs)
        if xtitle is not None:
            if xtitle != '':
                #plot_panel_ax.set_xlabel(xtitle, fontsize=xtitlefontsize)
                self.set_xtitle(xtitle, ax=plot_panel_ax, useTex=useTex)
                if current>0: self.Plot_panels[current-1]['xtitle'] = xtitle
        if ytitle is not None:
            if ytitle != '':
                #plot_panel_ax.set_ylabel(ytitle, fontsize=ytitlefontsize)
                self.set_ytitle(ytitle, ax=plot_panel_ax, useTex=useTex)
                if current>0: self.Plot_panels[current-1]['ytitle'] = ytitle
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
                                print('plot_data_file: clipping data from %s to %s'%(clip_x[0],clip_x[1]))
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
    def show(self, block=True):
        self.Plot_device.canvas.draw()
        self.Plot_device.show()
        pyplot.show(block=block)
        #self.Plot_device.waitforbuttonpress()
    # 
    def savefig(self, filename):
        self.Plot_device.savefig(filename)
    # 
    def savepdf(self, filename):
        if not filename.endswith('.pdf') and not filename.endswith('.PDF'):
            filename = filename + '.pdf'
        self.Plot_device.savefig(filename)
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















