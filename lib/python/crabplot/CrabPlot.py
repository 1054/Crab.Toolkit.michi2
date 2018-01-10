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
import astropy
import astropy.io.ascii as asciitable
#from astropy import units
#from astropy.io import fits
#from astropy.wcs import WCS
#import wcsaxes
#from scipy import spatial
#from scipy.spatial import KDTree
from pprint import pprint


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
        nrows = long(nrows)
        ncols = long(ncols)
        
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
        nrows = long(nrows)
        ncols = long(ncols)
        
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
        # determine label
        label = self.check_label(label)
        # append to self.Plot_panels
        self.Plot_panels.append( { 'label': label, 'panel': panel, 'position': position, 
                                   'x': None, 'y': None, 'xerr': None, 'yerr': None, 'xlog': False, 'ylog': False, 'xrange': [], 'yrange': [], 
                                   'image_data': None, 'image_wcs': None } )
        return self.Plot_panels[-1]
    # def totuple
    def convert_array_to_tuple(self,a):
        try:
            return tuple(self.convert_array_to_tuple(i) for i in a)
        except TypeError:
            return a
    # 
    def get_color_by_value(self, data_array, min_value = None, max_value = None, log = False, cmap = matplotlib.cm.cool):
        # add color bar by data_array
        data_arr = numpy.array(data_array)
        if min_value is None:
            min_value = numpy.min(data_arr)
        if max_value is None:
            max_value = numpy.max(data_arr)
        if log is True:
            data_norm = matplotlib.colors.LogNorm(vmin=min_value, vmax=max_value, clip=True)
        else:
            data_norm = matplotlib.colors.Normalize(vmin=min_value, vmax=max_value, clip=True)
        color_mapper = matplotlib.cm.ScalarMappable(norm=data_norm, cmap=cmap)
        color_by_data_array = self.convert_array_to_tuple(color_mapper.to_rgba(data_arr))
        return color_by_data_array
    # 
    def plot_xy(self, x, y, xerr = None, yerr = None, xlog = False, ylog = False, xrange = [], yrange = [], 
                xtitle = None, ytitle = None, 
                fmt = None, 
                symbol = 'o', size = 3, color = 'blue', fillstyle = None, 
                linestyle = 'None', linewidth = None, drawstyle = None, 
                facecolor = None, edgecolor = None, edgewidth = None, alpha = 1.0, 
                position = None, label = None, current = False, overplot = False):
        # check x y type
        x = self.check_array(x)
        y = self.check_array(y)
        # check x y dimension
        if x.shape != y.shape:
            print('CrabPlot::plot_xy() Error! The input x and y do not have the same shape!')
            return
        # get current panel or add new panel
        if len(self.Plot_panels) == 0:
            current = False
            overplot = False
        if current or overplot:
            # get current panel
            plot_panel_xy = self.Plot_panels[-1]
        else:
            # add new panel
            plot_panel_xy = self.add_panel(position = position, label = label)
        # set range
        if len(xrange) >= 2:
            plot_panel_xy['panel'].set_xlim(xrange[0], xrange[1])
        if len(yrange) >= 2:
            plot_panel_xy['panel'].set_ylim(yrange[0], yrange[1])
        # set grid
        plot_panel_xy['panel'].grid(False)
        # plot data
        if xlog and ylog:
            if fmt:
                plot_panel_xy['panel'].loglog(x, y, fmt)
            else:
                plot_panel_xy['panel'].loglog(x, y, marker=symbol, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha)
        elif xlog and not ylog:
            plot_panel_xy['panel'].semilogx(x, y, marker=symbol, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha)
        elif not xlog and ylog:
            plot_panel_xy['panel'].semilogy(x, y, marker=symbol, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha)
        else:
            plot_panel_xy['panel'].plot(x, y, marker=symbol, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, linewidth=linewidth, drawstyle=drawstyle, alpha=alpha)
            #plot_panel_xy['panel'].scatter(x, y, marker=symbol, s=size**2, c=color, edgecolors=edgecolor, linewidths=linewidth, alpha=alpha)
            if xerr or yerr:
                plot_panel_xy['panel'].errorbar(x, y, xerr = xerr, yerr = yerr, marker=symbol, markersize=size, color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markeredgewidth=edgewidth, fillstyle=fillstyle, linestyle=linestyle, drawstyle=drawstyle, linewidth=linewidth, alpha=alpha)
        # axis ticks format
        # -- http://matplotlib.org/examples/ticks_and_spines/tick-locators.html
        # -- http://stackoverflow.com/questions/33126126/matplotlib-minor-ticks
        if xlog:
            plot_panel_xy['panel'].xaxis.set_tick_params(which='major', length=8.0)
            plot_panel_xy['panel'].xaxis.set_tick_params(which='minor', length=4.0)
        if ylog:
            plot_panel_xy['panel'].yaxis.set_tick_params(which='major', length=8.0)
            plot_panel_xy['panel'].yaxis.set_tick_params(which='minor', length=4.0)
            plot_panel_xy['panel'].yaxis.set_major_locator(LogLocator(base=10,numticks=30))
            plot_panel_xy['panel'].yaxis.set_minor_locator(LogLocator(base=10,numticks=30,subs=numpy.arange(2.0,10.0,1.0)))
            plot_panel_xy['panel'].yaxis.set_minor_formatter(NullFormatter())
            #print(plot_panel_xy['panel'].yaxis.get_major_formatter())
            #plot_panel_xy['panel'].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            #plot_panel_xy['panel'].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        # plot title
        if xtitle:
            plot_panel_xy['panel'].set_xlabel(xtitle, fontsize=14)
        if ytitle:
            plot_panel_xy['panel'].set_ylabel(ytitle, fontsize=14)
        # store data
        plot_panel_xy['x'] = x
        plot_panel_xy['y'] = y
        plot_panel_xy['xerr'] = xerr
        plot_panel_xy['yerr'] = yerr
        plot_panel_xy['xlog'] = xlog
        plot_panel_xy['ylog'] = ylog
        if len(xrange) >= 2:
            plot_panel_xy['xrange'] = xrange
        if len(yrange) >= 2:
            plot_panel_xy['yrange'] = yrange
    # 
    def plot_line(self, x0, y0, x1 = None, y1 = None, ax = None, NormalizedCoordinate = False, overplot = True, **kwargs):
        # check x1 y1
        # we can input either (x0,y0), or (x0,y0,x1,y1).
        if x1 is not None and y1 is not None:
            # check x y type must be scalar in this case
            x0 = self.check_scalar(x0)
            y0 = self.check_scalar(y0)
            x1 = self.check_scalar(x1)
            y1 = self.check_scalar(y1)
            # get plot panel 'ax' variable
            if ax is not None:
                plot_panel_ax = ax
            else:
                if overplot: 
                    plot_panel_xy = self.Plot_panels[-1]
                else:
                    plot_panel_xy = self.add_panel(position = position, label = label)
                plot_panel_ax = plot_panel_xy['panel']
            # plot line
            if plot_panel_ax:
                if NormalizedCoordinate is True:
                    #xrange_for_plot = plot_panel_ax.get_xlim()
                    #yrange_for_plot = plot_panel_ax.get_ylim()
                    #x0_for_plot = x0 * (xrange_for_plot[1]-xrange_for_plot[0]) + xrange_for_plot[0]
                    #y0_for_plot = y0 * (yrange_for_plot[1]-yrange_for_plot[0]) + yrange_for_plot[0]
                    plot_one_line = matplotlib.lines.Line2D([x0,x1], [y0,y1], transform=plot_panel_ax.transAxes, **kwargs)
                else:
                    plot_one_line = matplotlib.lines.Line2D([x0,x1], [y0,y1], **kwargs)
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
            # sort
            x_sorted_index = x0.argsort()
            x0 = x0[x_sorted_index]
            y0 = y0[x_sorted_index]
            # get plot panel 'ax' variable
            if ax is not None:
                plot_panel_ax = ax
            else:
                if overplot: 
                    plot_panel_xy = self.Plot_panels[-1]
                else:
                    plot_panel_xy = self.add_panel(position = position, label = label)
                plot_panel_ax = plot_panel_xy['panel']
            # plot line
            if plot_panel_ax:
                if NormalizedCoordinate is True:
                    #xrange_for_plot = plot_panel_ax.get_xlim()
                    #yrange_for_plot = plot_panel_ax.get_ylim()
                    #x0_for_plot = x0 * (xrange_for_plot[1]-xrange_for_plot[0]) + xrange_for_plot[0]
                    #y0_for_plot = y0 * (yrange_for_plot[1]-yrange_for_plot[0]) + yrange_for_plot[0]
                    plot_one_line = matplotlib.lines.Line2D(x0, y0, transform=plot_panel_ax.transAxes, **kwargs)
                else:
                    plot_one_line = matplotlib.lines.Line2D(x0, y0, **kwargs)
                # 
                plot_panel_ax.add_line(plot_one_line)
    # 
    def plot_text(self, x0, y0, text_input, ax = None, NormalizedCoordinate = False, **kwargs):
        # check x y type
        x0 = self.check_scalar(x0)
        y0 = self.check_scalar(y0)
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
    def xyouts(self, x0, y0, text_input, ax = None, NormalizedCoordinate = False, **kwargs):
        self.plot_text(x0, y0, text_input, ax = ax, NormalizedCoordinate = NormalizedCoordinate, **kwargs)
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
    def show(self, block=True):
        self.Plot_device.canvas.draw()
        self.Plot_device.show()
        pyplot.show(block=block)
        #self.Plot_device.waitforbuttonpress()
    # 










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


















