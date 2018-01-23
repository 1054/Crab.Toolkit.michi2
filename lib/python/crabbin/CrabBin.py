#!/usr/bin/env python
# 

#####################################
# 
# Functions for data binning
# 
#   Last update: 
#            20180122, initialized
# 
#####################################

try:
    import pkg_resources
except ImportError:
    raise SystemExit("Error! Failed to import pkg_resources!")

pkg_resources.require("numpy")
pkg_resources.require("scipy")
pkg_resources.require("astropy")

import numpy
import scipy
from scipy import optimize, interpolate
from pprint import pprint
from copy import copy


def crab_bin_compute_param_chisq_histogram(chisq_array, param_array, min = None, max = None, nbin = None, step = None, log = False, delta_chisq = 2.3, verbose = 1):
    # 
    # compute chisq min
    chisq_array_copy = numpy.array(copy(chisq_array))
    chisq_min = numpy.nanmin(chisq_array_copy)
    # 
    # get the best min chisq param value
    param_mask_min_chisq = (numpy.abs(chisq_array_copy-chisq_min)<1e-6)
    param_best_min_chisq = param_array[param_mask_min_chisq]
    if len(param_best_min_chisq)>1:
        param_best_min_chisq = numpy.mean(param_best_min_chisq)
    # 
    # remove nan
    param_array_copy = numpy.array(copy(param_array))
    param_array_nonan = numpy.array(copy(param_array))
    param_array_nonan[numpy.isnan(param_array_nonan)] = -99
    # 
    # compute param log
    if log == True: 
        param_log_mask = (param_array_copy>0)
        param_log_mask2 = (param_array_copy<=0)
        param_array_copy[param_log_mask] = numpy.log10(param_array_copy[param_log_mask])
        param_array_copy[param_log_mask2] = numpy.nan
        param_array_nonan[param_log_mask] = numpy.log10(param_array_nonan[param_log_mask]) # NoNaN array
        if param_best_min_chisq>0:
            param_best_min_chisq = numpy.log10(param_best_min_chisq)
    # 
    # compute param min max
    param_min = numpy.nanmin(param_array_copy)
    param_max = numpy.nanmax(param_array_copy)
    # 
    # apply user param min max
    if min is not None:
        param_min = min
    if max is not None:
        param_max = max
    # 
    # apply log
    if log == True: 
        param_min = numpy.log10(param_min)
    if log == True: 
        param_max = numpy.log10(param_max)
    # 
    # set bin numb
    param_bin_numb = 50 #<TODO># 
    if nbin is not None:
        param_bin_numb = nbin
    # 
    # compute bin step
    param_bin_step = (param_max-param_min)/param_bin_numb
    # 
    # get xrange, the parameter range where chisq < chisq_min+delta_chisq
    xrange = numpy.array([param_min, param_max])
    valid = False
    # 
    # optimize xrange and bin step
    xclip_mask = (chisq_array_copy <= chisq_min+delta_chisq)
    if len(numpy.argwhere(xclip_mask)) > 0:
        xclip_min = numpy.nanmin(param_array_copy[xclip_mask])
        xclip_max = numpy.nanmax(param_array_copy[xclip_mask])
        if xclip_max > xclip_min:
            xrange = [xclip_min, xclip_max]
            valid = True
            param_bin_numb = 15 #<TODO># 
            if nbin is not None:
                param_bin_numb = nbin
            param_bin_step = float(xclip_max - xclip_min) / param_bin_numb #<TODO># bin step, do 15. 
            param_bin_numb = int(float(param_max-param_min)/param_bin_step)+1
    # 
    # apply user input bin step
    if step is not None:
        param_bin_step = step
        param_bin_numb = int((param_max-param_min)/param_bin_step)+1 #<TODO># int()
    # 
    # optimize yrange
    yrange = numpy.array([chisq_min, chisq_min+delta_chisq]) # [1/(chisq_min+delta_chisq), 1/(chisq_min)]
    # 
    # prepare binning
    param_bin_x = []
    param_bin_y = []
    # 
    if verbose>=1:
        print('Binning param, min max %s %s, valid range for (chi2<chi2min+%s) %s, step %s, nbin %s'%(param_min, param_max, delta_chisq, xrange, param_bin_step, param_bin_numb))
    # 
    for i in range(param_bin_numb):
        param_bin_lower = numpy.nan
        param_bin_upper = numpy.nan
        if i == param_bin_numb-1:
            param_bin_lower = param_min+float(i)*param_bin_step
            param_bin_upper = param_min+float(i+1)*param_bin_step # param_max
            param_bin_mask = (param_array_nonan>=param_bin_lower) & (param_array_nonan<=param_bin_upper)
        else:
            param_bin_lower = param_min+float(i)*param_bin_step
            param_bin_upper = param_min+float(i+1)*param_bin_step
            param_bin_mask = (param_array_nonan>=param_bin_lower) & (param_array_nonan<param_bin_upper)
        # check if there are data points in a parameter bin
        param_bin_args = numpy.argwhere(param_bin_mask)
        if verbose>=2:
            print('Binning from param value %s to %s, step %s, count %d'%(param_bin_lower, param_bin_upper, param_bin_step, len(param_bin_args)))
        if len(param_bin_args) > 0:
            chisq_bin_array = chisq_array_copy[param_bin_mask]
            chisq_bin_min = numpy.nanmin(chisq_bin_array)
            chisq_bin_max = numpy.nanmax(chisq_bin_array)
            param_bin_x.append(param_min+i*param_bin_step)
            param_bin_y.append(chisq_bin_min)
    # 
    # get best value
    param_bin_x = numpy.array(param_bin_x)
    param_bin_y = numpy.array(param_bin_y)
    param_best_x = numpy.nan
    if len(param_bin_y) > 0:
        param_bin_mask = (numpy.abs(param_bin_y-numpy.nanmin(param_bin_y))<1e-6)
        param_best_x = param_bin_x[param_bin_mask]
        if len(param_best_x) > 1:
            param_best_x = numpy.mean(param_best_x)
    else:
        print('*******************************************************')
        print('Error! param_bin_x param_bin_y could not be determined!')
        print('*******************************************************')
    # 
    # prepare output
    param_stats = {}
    param_stats['min'] = param_min
    param_stats['max'] = param_max
    param_stats['valid'] = valid
    param_stats['median'] = numpy.nan
    param_stats['sigma'] = numpy.nan
    param_stats['L68'] = numpy.nan
    param_stats['H68'] = numpy.nan
    if valid is True:
        param_stats['median'] = (xrange[0]+xrange[1])/2.0
        param_stats['sigma'] = (xrange[1]-xrange[0])/2.0
        param_stats['L68'] = xrange[0]
        param_stats['H68'] = xrange[1]
    param_stats['best_in_range'] = param_best_x
    param_stats['best_min_chisq'] = param_best_min_chisq
    param_stats['best'] = param_best_min_chisq
    param_stats['hist_x'] = param_bin_x
    param_stats['hist_y'] = param_bin_y
    param_stats['bin_step'] = param_bin_step
    param_stats['bin_numb'] = param_bin_numb
    param_stats['xrange'] = xrange # for plotting purpose
    param_stats['yrange'] = yrange # for plotting purpose
    param_stats['delta_chisq'] = delta_chisq
    param_stats['minimum_chisq'] = chisq_min
    # 
    # return
    return param_stats








