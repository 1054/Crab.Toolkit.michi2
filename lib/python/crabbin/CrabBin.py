#!/usr/bin/env python
# 

#####################################
# 
# Functions for data binning
# 
#   Last update: 
#            20180122, initialized
#            20200519, dynamical binning, optimized for the range of interest
#            20200710, fixed param_global_best = numpy.log10(param_global_best)
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
from copy import copy, deepcopy


def check_array(arr):
    # check array type numpy.ndarray
    if type(arr) is not list and type(arr) is not numpy.ndarray:
        arr = [arr]
    if type(arr) is list:
        arr = numpy.array(arr)
    return arr



def spline(input_x, input_y, output_x, xlog=False, ylog=False, outputxlog=False, outputylog=False, \
        fill_value=0.0, **kwargs):
    # spline
    # note that we can input xlog, ylog, outputxlog, outputylog to deal with logarithm input/output needs
    # xlog>0 means the input_x array is in linear space and will be converted to log space when doing the spline
    # ylog>0 means the input_y array is in linear space and will be converted to log space when doing the spline
    # outputxlog>0 means the output_x array is in linear space and will be converted to log space when doing the spline
    # outputylog>0 means the output_y array is in linear space but is in log space when doing the spline, so we need to convert it back to linear space after the spline.
    input_x_coord = check_array(input_x)
    input_y_value = check_array(input_y)
    output_x_coord = check_array(output_x)
    if xlog:
        input_x_mask = (input_x_coord<=0.0)
        input_x_coord[input_x_mask] = numpy.nan
        input_x_coord = numpy.log10(input_x)
    if ylog:
        input_y_mask = (input_y_value<=0.0)
        input_y_value[input_y_mask] = numpy.nan
        input_y_value = numpy.log10(input_y_value)
    # 
    #if outputxlog is None:
    #    outputxlog = xlog
    #if outputylog is None:
    #    outputylog = ylog
    # 
    if xlog:
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
    output_mask = numpy.logical_or( (output_x_coord<numpy.nanmin(input_x_coord)), (output_x_coord>numpy.nanmax(input_x_coord)) )
    output_y_value[output_mask] = numpy.nan
    # 
    if ylog:
        if outputylog:
            output_y = output_y_value
        else:
            output_y = numpy.power(10, output_y_value)
            output_y[output_mask] = fill_value #<TODO># fill with 0.0 instead of nan
    else:
        if outputylog:
            output_y = numpy.log10(output_y_value)
        else:
            output_y = output_y_value
    # 
    return output_y



def interp(input_x, input_y, output_x, xlog=False, ylog=False, outputxlog=False, outputylog=False, \
        kind='linear', bounds_error=False, fill_value=numpy.nan, **kwargs):
    # interp
    # note that we can input xlog, ylog, outputxlog, outputylog to deal with logarithm input/output needs
    # xlog>0 means the input_x array is in linear space and will be converted to log space when doing the interp
    # ylog>0 means the input_y array is in linear space and will be converted to log space when doing the interp
    # outputxlog>0 means the output_x array is in linear space and will be converted to log space when doing the interp
    # outputylog>0 means the output_y array is in linear space but is in log space when doing the interp, so we need to convert it back to linear space after the spline.
    input_x_coord = check_array(input_x)
    input_y_value = check_array(input_y)
    output_x_coord = check_array(output_x)
    if xlog>0:
        input_x_mask = (input_x_coord<=0.0)
        input_x_coord[input_x_mask] = numpy.nan
        input_x_coord = numpy.log10(input_x)
    if ylog>0:
        input_y_mask = (input_y_value<=0.0)
        input_y_value[input_y_mask] = numpy.nan
        input_y_value = numpy.log10(input_y_value)
    # 
    #if outputxlog is None:
    #    outputxlog = xlog
    #if outputylog is None:
    #    outputylog = ylog
    # 
    if xlog:
        output_x_mask = (output_x_coord<=0.0)
        output_x_coord[output_x_mask] = numpy.nan
        output_x_coord = numpy.log10(output_x_coord)
    # 
    #print(numpy.column_stack((input_x_coord, input_y_value)))
    input_mask = numpy.isnan(input_y_value)
    f = scipy.interpolate.interp1d(input_x_coord, input_y_value, kind=kind, bounds_error=bounds_error, fill_value=fill_value, **kwargs) # axis=0, bc_type='not-a-knot', extrapolate=None
    output_y_value = f(output_x_coord)
    # 
    # deal with data out of X range
    output_mask = numpy.logical_or( (output_x_coord<numpy.nanmin(input_x_coord)), (output_x_coord>numpy.nanmax(input_x_coord)) )
    output_y_value[output_mask] = numpy.nan
    # 
    if ylog:
        if outputylog:
            output_y = output_y_value
        else:
            output_y = numpy.power(10, output_y_value)
            output_y[output_mask] = fill_value #<TODO># fill with 0.0 instead of nan
    else:
        if outputylog:
            output_y = numpy.log10(output_y_value)
        else:
            output_y = output_y_value
    # 
    return output_y



def crab_bin_compute_param_chisq_histogram(chisq_array, param_array, \
        min = None, max = None, nbin = None, nbinmax = 10000, step = None, smooth = None, 
        log = False, delta_chisq = 2.3, verbose = 1):
    # 
    # compute chisq min
    chisq_array_copy = numpy.array(deepcopy(chisq_array))
    chisq_global_min = numpy.nanmin(chisq_array_copy)
    # 
    # compute param min
    param_array_copy = numpy.array(deepcopy(param_array))
    # 
    # appply log if set
    if log:
        with numpy.errstate(invalid='ignore'):
            param_log_mask = (param_array_copy>0)
            param_array_copy[param_log_mask] = numpy.log10(param_array_copy[param_log_mask])
            param_array_copy[~param_log_mask] = numpy.nan
    # 
    param_global_min = numpy.nanmin(param_array_copy)
    param_global_max = numpy.nanmax(param_array_copy)
    # 
    # get the best min chisq param value
    mask_global_min_chisq = (numpy.abs(chisq_array_copy-chisq_global_min)<1e-6)
    param_global_best = param_array[mask_global_min_chisq]
    if len(param_global_best)>1:
        param_global_best = numpy.mean(param_global_best)
    if log:
        param_global_best = numpy.log10(param_global_best) # fixed 20200710
    # 
    # compute param min max
    param_min = numpy.nanmin(param_array_copy)
    param_max = numpy.nanmax(param_array_copy)
    # 
    # remove nan
    param_array_nonan = deepcopy(param_array_copy)
    param_array_nonan[numpy.isnan(param_array_copy)] = param_global_min-99
    # 
    # apply user param min max
    if min is not None:
        if log: 
            param_min = numpy.log10(min)
        else:
            param_min = min
    if max is not None:
        if log: 
            param_max = numpy.log10(max)
        else:
            param_max = max
    # 
    # select param array in range and param value not to be zero if log
    mask_in_range = numpy.logical_and.reduce( ( ~numpy.isnan(param_array_copy), (param_array_nonan>=param_min), (param_array_nonan<=param_max) ) )
    chisq_array_copy[~mask_in_range] = numpy.nan
    # 
    # recompute chisq array min
    chisq_min = numpy.nanmin(chisq_array_copy)
    chisq_max = numpy.nanmax(chisq_array_copy)
    chisq_in_range_min = chisq_min
    chisq_in_range_max = chisq_max
    # 
    # filter out chisq array nan
    chisq_array_nonan = deepcopy(chisq_array_copy)
    chisq_array_nonan[numpy.isnan(chisq_array_nonan)] = chisq_max+99
    # 
    # get the best min chisq param value
    param_in_range_best = param_array_copy[(numpy.abs(chisq_array_nonan-chisq_min)<1e-6)]
    if len(param_in_range_best)>1:
        param_in_range_best = numpy.mean(param_in_range_best)
    # 
    # mask valid data within the (chisq <= chisq_min+delta_chisq) range of interest
    valid_mask = (chisq_array_nonan <= chisq_min+delta_chisq)
    # 
    # prepare the 'valid' flag for any valid data with chisq <= chisq_min+delta_chisq or not
    valid = (numpy.count_nonzero(valid_mask) > 0)
    # 
    # define xrange yrange for the (chisq <= chisq_min+delta_chisq) range of interest
    xrange = numpy.array([param_min, param_max])
    yrange = numpy.array([chisq_min, chisq_max])
    valid_xrange = numpy.array([numpy.nan, numpy.nan])
    valid_yrange = numpy.array([numpy.nan, numpy.nan])
    smooth_xrange = numpy.array([numpy.nan, numpy.nan])
    smooth_yrange = numpy.array([numpy.nan, numpy.nan])
    if valid:
        valid_xrange = numpy.array([param_min, param_max])
        valid_yrange = numpy.array([chisq_min, chisq_min+delta_chisq])
    # 
    # set a flag to redefine bin edges later
    redefine_bin_edges = True
    # 
    # prepare to compute bin step and edges
    if step is not None:
        param_bin_step = step
        param_bin_numb = int(numpy.ceil((param_max-param_min)/param_bin_step))
        param_bin_edge = numpy.linspace(param_min, param_max, num=param_bin_numb+1, endpoint=True)
    elif nbin is not None:
        param_bin_numb = nbin
        param_bin_step = (param_max-param_min)/param_bin_numb
        param_bin_edge = numpy.linspace(param_min, param_max, num=param_bin_numb+1, endpoint=True)
        redefine_bin_edges = False
    else:
        param_bin_numb = 100 #<TODO># 
        param_bin_step = (param_max-param_min)/param_bin_numb
        param_bin_edge = numpy.linspace(param_min, param_max, num=param_bin_numb+1, endpoint=True)
    # 
    # redefine bin step (optimize for valid data within chisq <= chisq_min+delta_chisq range)
    # if user has input a step use it
    # if user has input a delta_chisq, we recompute the bin step within that chisq range.
    if valid and redefine_bin_edges:
        xclip_min = numpy.nanmin(param_array_copy[valid_mask])
        xclip_max = numpy.nanmax(param_array_copy[valid_mask])
        if xclip_min == xclip_max:
            # if there is only one single param solution, how to compute the width?? 
            #<TODO># for now, apply a width of 0.1 * value
            valid_xrange = [xclip_min-0.05*xclip_min, xclip_max+0.05*xclip_max]
        else:
            valid_xrange = [xclip_min, xclip_max]
            # 
            #param_bin_numb = 15 #<TODO># 
            #param_bin_step2 = float(xclip_max - xclip_min) / param_bin_numb # do a bin step of 15 in the delta_chisq range. 
            #if param_bin_step2 < param_bin_step:
            #    param_bin_step = param_bin_step2 # if the bin step is not better, use the old bin step from the whole param range
            #    param_bin_numb = int(numpy.ceil(float(param_max-param_min)/param_bin_step))
            # 
            xclip3_min = numpy.nanmin(param_array_copy[1./chisq_array_nonan >= 0.362*numpy.max(1./chisq_array_nonan)]) # expand 2 times delta_chisq
            xclip3_max = numpy.nanmax(param_array_copy[1./chisq_array_nonan >= 0.362*numpy.max(1./chisq_array_nonan)]) # expand 2 times delta_chisq
            xclip2_min = numpy.nanmin(param_array_copy[1./chisq_array_nonan >= 0.618*numpy.max(1./chisq_array_nonan)]) # expand 2 times delta_chisq
            xclip2_max = numpy.nanmax(param_array_copy[1./chisq_array_nonan >= 0.618*numpy.max(1./chisq_array_nonan)]) # expand 2 times delta_chisq
            param_bin_edge = []
            # print(
            #     param_min, xclip3_min, xclip2_min, xclip_min, '\n', 
            #     xclip_max, xclip2_max, xclip3_max, param_max, '\n',
            #     ) # debugging 20221101
            xcursor = param_min
            if xcursor < xclip3_min:
                param_bin_edge.extend(numpy.arange(xcursor, xclip3_min, param_bin_step).tolist()) # 
                xcursor = xclip3_min
            if xcursor < xclip2_min:
                param_bin_edge.extend(numpy.arange(xcursor, xclip2_min, param_bin_step*0.8).tolist()) # 
                xcursor = xclip2_min
            if xcursor < xclip_min:
                param_bin_edge.extend(numpy.arange(xcursor, xclip_min, param_bin_step*0.7).tolist()) # 
                xcursor = xclip_min
            if xcursor < xclip_max:
                param_bin_edge.extend(numpy.arange(xcursor, xclip_max, param_bin_step*0.6).tolist()) # 
                xcursor = xclip_max
            if xcursor < xclip2_max:
                param_bin_edge.extend(numpy.arange(xcursor, xclip2_max, param_bin_step*0.7).tolist()) # using finer grid for the range of interest
                xcursor = xclip2_max
            if xcursor < xclip3_max:
                param_bin_edge.extend(numpy.arange(xcursor, xclip3_max, param_bin_step*0.8).tolist()) # 
                xcursor = xclip3_max
            if xcursor < param_max:
                param_bin_edge.extend(numpy.arange(xcursor, param_max, param_bin_step).tolist()) # 
                xcursor = param_max
            param_bin_edge = numpy.array(param_bin_edge)
            param_bin_numb = len(param_bin_edge)-1
            # 
    #print('param_bin_edge:', param_bin_edge)
    # 
    # apply user input nbinmax
    #if nbinmax is not None:
    #    if param_bin_numb > nbinmax:
    #        param_bin_numb = nbinmax
    #        param_bin_step = float(param_max - param_min) / param_bin_numb
    # 
    # 
    # prepare binning
    param_bin_x = [] # histogram x, i.e., param value
    param_bin_y = [] # histogram y, i.e., minchisq value
    param_bin_dx = [] # histogram x bin width
    # 
    if verbose>=1:
        #print('Binning param, min max %s %s, valid range for (chi2<chi2min+%s) %s, step %s, nbin %s'%(param_min, param_max, delta_chisq, xrange, param_bin_step, param_bin_numb))
        #print('Binning param, min max %s %s, delta_chisq %s, chisq_min+delta_chisq %s, valid range %s %s, step %s, nbin %s'%(param_min, param_max, delta_chisq, chisq_min+delta_chisq, valid_xrange[0], valid_xrange[1], param_bin_step, param_bin_numb))
        print('Binning param, min max %s %s, chisq_min %s, delta_chisq %s, chisq_min+delta_chisq %s, valid range %s %s, nbin %s'%(param_min, param_max, chisq_min, delta_chisq, chisq_min+delta_chisq, valid_xrange[0], valid_xrange[1], param_bin_numb))
    # 
    param_bin_lower = numpy.nan
    param_bin_upper = numpy.nan
    for i in range(len(param_bin_edge)-1):
        param_bin_lower = param_bin_edge[i] # param_min+float(i)*param_bin_step
        param_bin_upper = param_bin_edge[i+1] # param_min+float(i+1)*param_bin_step
        if i == len(param_bin_edge)-1-1:
            param_bin_mask = numpy.logical_and( (param_array_nonan>=param_bin_lower), (param_array_nonan<=param_bin_upper) )
        else:
            param_bin_mask = numpy.logical_and( (param_array_nonan>=param_bin_lower), (param_array_nonan<param_bin_upper) )
        param_bin_mask_delta_chisq = numpy.logical_and( param_bin_mask, (chisq_array_nonan <= chisq_min+delta_chisq) )
        # check if there are data points in a parameter bin
        if numpy.count_nonzero(param_bin_mask) > 0:
            chisq_bin_array = chisq_array_copy[param_bin_mask]
            chisq_bin_min = numpy.nanmin(chisq_bin_array)
            chisq_bin_max = numpy.nanmax(chisq_bin_array)
            param_bin_x.append(param_bin_lower)
            param_bin_y.append(chisq_bin_min)
            param_bin_dx.append(param_bin_upper-param_bin_lower)
            if verbose>=2:
                print('Binning from param value %s to %s,  bin minchisq %s, count %d, within delta_chisq count %s'%(\
                        param_bin_lower, param_bin_upper, chisq_bin_min, \
                        numpy.count_nonzero(param_bin_mask), numpy.count_nonzero(param_bin_mask_delta_chisq) ) )
        else:
            if verbose>=2:
                print('Binning from param value %s to %s, count %d, within delta_chisq count %s'%(\
                        param_bin_lower, param_bin_upper, \
                        numpy.count_nonzero(param_bin_mask), numpy.count_nonzero(param_bin_mask_delta_chisq) ) )
    param_bin_x = numpy.array(param_bin_x)
    param_bin_y = numpy.array(param_bin_y)
    param_bin_dx = numpy.array(param_bin_dx)
    
    # smooth if needed
    if smooth is not None:
        try:
            smooth = float(smooth)
        except:
            smooth = 0
        if smooth > 1.0:
            #param_bin_y = numpy.convolve(param_bin_y, numpy.ones(int(smooth))/float(smooth), mode='same')
            gausshalflen = int(numpy.ceil(5.*smooth)/2)
            gausslen = gausshalflen*2+1
            gaussx = numpy.arange(gausslen)-float(gausshalflen)
            gaussy = numpy.exp(-0.5 * (gaussx/(smooth/2.35482))**2)
            param_bin_y = numpy.convolve(param_bin_y, gaussy/numpy.sum(gaussy), mode='same')
            param_bin_y = param_bin_y - numpy.min(param_bin_y) + chisq_min
    
    # 
    # spline the histogram to a finer grid
    #spline_x = numpy.arange(param_min, param_max+0.5*param_bin_step, param_bin_step)
    #spline_y = spline(param_bin_x, param_bin_y, spline_x) # here we do not need xlog ylog because param_bin_x and param_bin_y have already been converted to log
    smooth_x = copy(param_bin_x)
    smooth_y = copy(param_bin_y)
    smooth_dx = copy(param_bin_dx)
    # 
    # and get smooth xrange yrange
    if valid and len(param_bin_x)>1:
        #xclip_min = numpy.nanmin(param_array_copy[chisq_array_nonan <= chisq_min+2.0*delta_chisq]) # expand 2 times delta_chisq
        #xclip_max = numpy.nanmax(param_array_copy[chisq_array_nonan <= chisq_min+2.0*delta_chisq]) # expand 2 times delta_chisq
        #smooth_x = []
        #if param_min < xclip_min:
        #    smooth_x.extend(numpy.arange(param_min, xclip_min, 0.5*param_bin_step).tolist()) # using 2 times finer grid
        #if xclip_min < xclip_max:
        #    smooth_x.extend(numpy.arange(xclip_min, xclip_max, 0.2*param_bin_step).tolist()) # using 5 times finer grid for the range of interest
        #if xclip_max < param_max:
        #    smooth_x.extend(numpy.arange(xclip_max, param_max, 0.5*param_bin_step).tolist()) # using 2 times finer grid
        #smooth_x = numpy.array(smooth_x)
        #smooth_x = numpy.arange(param_min, param_max+0.25*param_bin_step, 0.25*param_bin_step) # using 4 times finer grid
        #smooth_y = interp(param_bin_x, param_bin_y, smooth_x, ylog=True) # here we do not need xlog because param_bin_x have already been converted to log, but we need ylog
        #print(list(zip(param_bin_x, param_bin_y)))
        #print(list(zip(smooth_x, smooth_y)))
        # 
        smooth_x = param_bin_edge[0:-1]
        #print('dzliu debugging 20230822 param_bin_x, param_bin_y:', list(zip(numpy.round(param_bin_x, 3), numpy.round(param_bin_y, 3))))
        smooth_y = interp(param_bin_x, param_bin_y, smooth_x, ylog=True) # here we do not need xlog because param_bin_x have already been converted to log, but we need ylog
        #print('dzliu debugging 20230822 smooth_x, smooth_y:', list(zip(numpy.round(smooth_x, 3), numpy.round(smooth_y, 3))))
        with numpy.errstate(invalid='ignore'):
            valid_indicies = numpy.argwhere(smooth_y <= chisq_min+delta_chisq).flatten()
            left_index = valid_indicies[0]
            right_index = valid_indicies[-1]+1 # use the right edge of the bin, if the bin is not the last bin
            smooth_xrange[0] = param_bin_edge[left_index]
            smooth_xrange[1] = param_bin_edge[right_index]
            smooth_dx = numpy.concatenate((numpy.diff(smooth_x), [param_bin_edge[-1]-param_bin_edge[-2]]), axis=None)
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
    # 
    # outputlog
    #if log:
    #    param_min = 10**param_min
    #    param_max = 10**param_max
    # 
    # 
    # 
    # prepare output
    param_stats = {}
    param_stats['valid'] = valid
    param_stats['min'] = param_min
    param_stats['max'] = param_max
    param_stats['median'] = numpy.nan
    param_stats['sigma'] = numpy.nan
    param_stats['L68'] = numpy.nan
    param_stats['H68'] = numpy.nan
    param_stats['global_min_chisq'] = chisq_global_min
    param_stats['global_best'] = param_global_best
    param_stats['in_range_min_chisq'] = chisq_in_range_min
    param_stats['in_range_best'] = param_in_range_best
    param_stats['in_range_min'] = param_min
    param_stats['in_range_max'] = param_max
    param_stats['best'] = param_global_best
    param_stats['bin_step'] = param_bin_step
    param_stats['bin_numb'] = param_bin_numb
    param_stats['bin_edge'] = param_bin_edge
    param_stats['hist_dx'] = param_bin_dx
    param_stats['hist_x'] = param_bin_x
    param_stats['hist_y'] = param_bin_y
    #param_stats['spline_x'] = spline_x
    #param_stats['spline_y'] = spline_y
    param_stats['smooth_dx'] = smooth_dx
    param_stats['smooth_x'] = smooth_x
    param_stats['smooth_y'] = smooth_y
    param_stats['xrange'] = xrange # best param range, = [param_min, param_max] if no valid data.
    param_stats['yrange'] = yrange # best chisq range, = [chisq_min, chisq_max] if no valid data.
    param_stats['valid_xrange'] = valid_xrange # (chisq <= chisq_min+delta_chisq) param range
    param_stats['valid_yrange'] = valid_yrange # (chisq <= chisq_min+delta_chisq) chisq range
    param_stats['smooth_xrange'] = smooth_xrange # splined (chisq <= chisq_min+delta_chisq) param range
    param_stats['smooth_yrange'] = smooth_yrange # splined (chisq <= chisq_min+delta_chisq) chisq range
    #param_stats['plot_xrange'] = xrange # param range for plotting, = [param_min, param_max] if no valid data.
    #param_stats['plot_yrange'] = yrange # chisq range for plotting, = [chisq_min, chisq_max] if no valid data.
    param_stats['delta_chisq'] = delta_chisq
    param_stats['min_chisq'] = chisq_min
    param_stats['max_chisq'] = chisq_max
    param_stats['minimum_chisq'] = chisq_min
    param_stats['threshold_chisq'] = chisq_min + delta_chisq
    # 
    #if valid is True:
    #    param_stats['median'] = (xrange[0]+xrange[1])/2.0
    #    param_stats['sigma'] = (xrange[1]-xrange[0])/2.0
    #    param_stats['L68'] = xrange[0]
    #    param_stats['H68'] = xrange[1]
    # 
    if valid is True:
        if not numpy.isnan(smooth_xrange[0]) and not numpy.isnan(smooth_xrange[1]):
            param_stats['median'] = (smooth_xrange[0]+smooth_xrange[1])/2.0
            param_stats['sigma'] = (smooth_xrange[1]-smooth_xrange[0])/2.0
            param_stats['L68'] = smooth_xrange[0]
            param_stats['H68'] = smooth_xrange[1]
            param_stats['xrange'][0] = smooth_xrange[0]
            param_stats['xrange'][1] = smooth_xrange[1]
        else:
            param_stats['median'] = (valid_xrange[0]+valid_xrange[1])/2.0
            param_stats['sigma'] = (valid_xrange[1]-valid_xrange[0])/2.0
            param_stats['L68'] = valid_xrange[0]
            param_stats['H68'] = valid_xrange[1]
            param_stats['xrange'][0] = valid_xrange[0]
            param_stats['xrange'][1] = valid_xrange[1]
    # 
    # return
    return param_stats








