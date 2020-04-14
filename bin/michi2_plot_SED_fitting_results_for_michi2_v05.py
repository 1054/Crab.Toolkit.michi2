#!/usr/bin/env python
# 

import os
import sys

#os.system('cp /Users/dzliu/Softwares/Python/lib/crab/crabplot/CrabPlot.py /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/lib/python/crabplot/CrabPlot.py')

sys.path.insert(1, os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabtable')
sys.path.insert(1, os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabplot')
sys.path.insert(1, os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabbin')
#sys.path.insert(1, '/Users/dzliu/Softwares/Python/lib/crab/crabtable')
#sys.path.insert(1, '/Users/dzliu/Softwares/Python/lib/crab/crabplot')
#sys.path.insert(1, '/Users/dzliu/Softwares/Python/lib/crab/crabbin')

from CrabTable import *
from CrabPlot import *
from CrabBin import *

import glob
import math
import numpy
import scipy
import astropy
from astropy import units
from astropy.io import fits
import astropy.io.ascii as asciitable
import re
import json
from copy import copy

#from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)

import matplotlib as mpl # https://matplotlib.org/users/customizing.html
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['axes.grid'] = True
#mpl.rcParams['grid.color'] = 'b0b0b0'
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['grid.linewidth'] = 0.25
mpl.rcParams['grid.alpha'] = 0.8






#########################################
#               Constants               #
#########################################

Delta_chisq_of_interest = 5.89 # for 5 interested parameters, 68.3% confidence (1-sigma).
# Delta_chisq_of_interest is the allowed chi-square range for computing errors in interested parameters.
# -- See http://www.astro.sunysb.edu/metchev/PHY515/astrostatistics.html
# -- See also Press (1992) "" Chap. 15 P. 6 (press1992chap15p6.pdf)
# -- The value can also be calculated by the C++ code "michi2_compute_delta_chisq_1sigma", which is written based on the code in Press (1992). 





#########################################
#               Functions               #
#########################################

def analyze_chisq_distribution(param_dict, verbose = 0, Plot_engine = None, Output_dir = ''):
    # Plot_engine must be the CrabPlot class
    if 'Lib_file' in param_dict and \
        'Lib_name' in param_dict and \
        'Lib_numb' in param_dict and \
        'Par_name' in param_dict and \
        'Col_numb' in param_dict and \
        'value' in param_dict and \
        'chisq' in param_dict :
        # 
        if verbose >= 1:
            print('Analyzing the chi-square distribution for parameter "%s" in library %s from file "%s"'%(param_dict['Par_name'], param_dict['Lib_name'], param_dict['Lib_file']))
            #print(param_dict)
        # 
        if Output_dir != '':
            if not os.path.isdir(Output_dir):
                os.makedirs(Output_dir)
            if not Output_dir.endswith(os.sep):
                Output_dir = Output_dir + os.sep
            
        # 
        #print(big_chisq_data_table.colnames)
        #chisq_array = big_chisq_data_table[big_chisq_data_table.colnames[2-1]] # the second column of the input 'big_chisq_data_table' is the chisq column. -- TODO: MUST MAKE SURE USER DO NOT CHANGE THE FORMAT!
        #coeff_array = big_chisq_data_table[big_chisq_data_table.colnames[2+2*param_dict['Lib_numb']-1]] # the coefficient/normalization column, e.g., for LIB(2), the 2+2*(2)th column is a2, which is the normalization of LIB(2).
        #param_array = big_chisq_data_table[big_chisq_data_table.colnames[param_dict['Col_numb']-1]] # the parameter column.
        #print('Taking column %s as the chisq_array'%(2))
        #print('Taking column %s as the coeff_array'%(2+2*param_dict['Lib_numb']))
        #print('Taking column %s as the param_array'%(param_dict['Col_numb']))
        # 
        # check arrays
        param_array = copy(param_dict['value'])
        chisq_array = copy(param_dict['chisq'])
        #if 'Degree_of_freedom' in param_dict:
        #    chisq_array = chisq_array / param_dict['Degree_of_freedom'] # make it reduced-chi-square <20191120>
        chisq_min = numpy.nanmin(chisq_array)
        # 
        # copy
        #param_array_nonan = copy(param_array)
        #param_array_nonan[numpy.isnan(param_array_nonan)] = -99
        # 
        # apply range -- no, do not cut the array, but just adjust the plotting range.
        #if 'range' in param_dict:
        #    if len(param_dict['range'])>=2:
        #        param_range_mask = (param_array >= param_dict['range'][0]) & (param_array <= param_dict['range'][1])
        #        chisq_array = chisq_array[param_range_mask]
        #        param_array = param_array[param_range_mask]
        #        param_array_nonan = param_array_nonan[param_range_mask]
        # 
        # set plot_xrange to the user-specified values
        plot_xrange = None
        param_min = None
        param_max = None
        if 'range' in param_dict:
            if len(param_dict['range'])>=2:
                param_min = param_dict['range'][0]
                param_max = param_dict['range'][1]
        # 
        # bin param
        param_log = False
        if 'Log_calc' in param_dict:
            if param_dict['Log_calc'] == True:
                param_log = True
        # 
        # crab_bin_compute_param_chisq_histogram for delta_chisq = 2.3 (2p)
        #verbose = 1
        param_stats_2p = crab_bin_compute_param_chisq_histogram(chisq_array, param_array, min = param_min, max = param_max, 
                            delta_chisq = 2.3, log = param_log, verbose = verbose)
        #print('param_stats_2p.valid', param_stats_2p['valid'])
        #print('param_stats_2p.threshold_chisq', param_stats_2p['threshold_chisq'], '1/x', 1/param_stats_2p['threshold_chisq'])
        #print('param_stats_2p.xrange', param_stats_2p['xrange']) # L68
        #print('param_stats_2p.yrange', param_stats_2p['yrange']) # H68
        #print('param_stats_2p.global_min_chisq', param_stats_2p['global_min_chisq'], '1/x', 1/param_stats_2p['global_min_chisq'])
        #print('param_stats_2p.global_best', param_stats_2p['global_best'])
        #print('param_stats_2p.in_range_min_chisq', param_stats_2p['in_range_min_chisq'])
        #print('param_stats_2p.in_range_best', param_stats_2p['in_range_best'])
        #print('param_stats_2p.in_range_min', param_stats_2p['in_range_min'])
        #print('param_stats_2p.in_range_max', param_stats_2p['in_range_max'])
        #print('param_stats_2p.min', param_stats_2p['min'])
        #print('param_stats_2p.max', param_stats_2p['max'])
        #print('param_stats_2p.median', param_stats_2p['median'])
        #print('param_stats_2p.best', param_stats_2p['best'])
        #print('param_stats_2p.sigma', param_stats_2p['sigma'])
        #print('param_stats_2p.L68', param_stats_2p['L68'])
        #print('param_stats_2p.H68', param_stats_2p['H68'])
        #sys.exit()
        if 'Par_file' in param_dict:
            # remove previous file
            if os.path.isfile(Output_dir+'best-fit_param_'+param_dict['Par_file']+'.txt'):
                os.system('mv %s %s.backup'%(Output_dir+'best-fit_param_'+param_dict['Par_file']+'.txt', 
                                             Output_dir+'best-fit_param_'+param_dict['Par_file']+'.txt'))
            if param_stats_2p['valid'] is True:
                param_median = param_stats_2p['median']
                param_best = param_stats_2p['best']
                param_sigma = param_stats_2p['sigma']
                param_L68 = param_stats_2p['L68']
                param_H68 = param_stats_2p['H68']
            else:
                param_median = 0.0
                param_best = 0.0
                param_sigma = 0.0
                param_L68 = 0.0
                param_H68 = 0.0
            #asciitable.write(numpy.column_stack((param_median, param_best, param_sigma, param_L68, param_H68)), 
            #                        Output_dir+'best-fit_param_'+param_dict['Par_file']+'.txt', Writer=asciitable.Ipac, 
            #                                names=['param_median', 'param_best', 'param_sigma', 'param_L68', 'param_H68'], 
            #                                formats={'param_median': '%20.10g', 'param_best': '%20.10g', 'param_sigma': '%20.10g', 'param_L68': '%20.10g', 'param_H68': '%20.10g'}, 
            #                                    delimiter='    ', overwrite = True)
            asciitable.write(numpy.column_stack((param_median, param_best, param_sigma, param_L68, param_H68)), 
                                    Output_dir+'best-fit_param_'+param_dict['Par_file']+'.txt', Writer=asciitable.Ipac, 
                                            names=['param_median', 'param_best', 'param_sigma', 'param_L68', 'param_H68'], 
                                            formats={'param_median': '%20.10g', 'param_best': '%20.10g', 'param_sigma': '%20.10g', 'param_L68': '%20.10g', 'param_H68': '%20.10g'}, 
                                             overwrite = True)
        # 
        # crab_bin_compute_param_chisq_histogram for plotting
        param_stats = crab_bin_compute_param_chisq_histogram(chisq_array, param_array, min = param_min, max = param_max, delta_chisq = Delta_chisq_of_interest, log = param_log, verbose = verbose)
        # 
        param_bin_x = param_stats['hist_x']
        param_bin_y = param_stats['hist_y']
        param_bin_step = param_stats['bin_step']
        # 
        plot_xrange = param_stats['xrange'] # xrange is the param range where chi-sq < min-chi-sq + 2.3
        plot_yrange = param_stats['yrange']
        #plot_xrange = [plot_xrange[0]-(plot_xrange[1]-plot_xrange[0])*0.50, plot_xrange[1]+(plot_xrange[1]-plot_xrange[0])*0.50] # extend the range for plotting.
        if param_stats['valid']:
            # zoomed in plot_xrange (right panel)
            plot_xrange = [ param_stats['L68'] - 10*param_stats['bin_step'], 
                            param_stats['H68'] + 10*param_stats['bin_step'] ]
        ##if plot_xrange[0] < param_stats['min']: plot_xrange[0] = param_stats['min']
        ##if plot_xrange[1] > param_stats['max']: plot_xrange[1] = param_stats['max']
        # invert y
        plot_yrange = [1.0/plot_yrange[1], 1.0/plot_yrange[0]]
        plot_yrange = numpy.log10(plot_yrange)
        plot_yrange = [plot_yrange[0]-(plot_yrange[1]-plot_yrange[0])*0.50, plot_yrange[1]+(plot_yrange[1]-plot_yrange[0])*0.05] # extend the range for plotting.
        plot_yrange = numpy.power(10,plot_yrange)
        # 
        xlog = None
        #if 'Log_plot' in param_dict:
        #    if param_dict['Log_plot'] == True:
        #        xlog = 1 # not working for matplotlib bar plot (i.e., CrabPlot plot_hist)!
        # 
        ylog = None
        # 
        # log
        if param_log is True:
            param_array_mask = (param_array>0)
            param_array_mask2 = (param_array<=0)
            param_array[param_array_mask] = numpy.log10(param_array[param_array_mask])
            param_array[param_array_mask2] = numpy.nan
        # 
        # verbose
        if verbose >= 2:
            #pprint(numpy.column_stack((param_bin_x, param_bin_y, 1/param_bin_y)))
            #print('------ xrange', plot_xrange)
            #print('------ yrange', plot_yrange, [1/plot_yrange[1],1/plot_yrange[0]])
            #print('------ param_stats.xrange', param_stats['xrange'])
            #print('------ param_stats.yrange', param_stats['yrange'], [1/param_stats['yrange'][1],1/param_stats['yrange'][0]])
            #print('------ param_stats.minimum_chisq', param_stats['minimum_chisq'])
            #print('------ param_stats.best_min_chisq', param_stats['best_min_chisq'])
            #print('------ param_stats.best', param_stats['best'])
            print('param_stats.min', param_stats['min'])
            print('param_stats.max', param_stats['max'])
            print('param_stats.xrange', param_stats['xrange'])
            print('param_stats.yrange', param_stats['yrange'], '(1/chisq ',1/param_stats['yrange'][1],' ',1/param_stats['yrange'][0],')')
            print('plotting xrange', plot_xrange)
            print('plotting yrange', plot_yrange)
            print('param_stats_2p.valid', param_stats_2p['valid'])
            print('param_stats_2p.threshold_chisq', param_stats_2p['threshold_chisq'], '1/x', 1/param_stats_2p['threshold_chisq'])
            print('param_stats_2p.xrange', param_stats_2p['xrange']) # L68
            print('param_stats_2p.yrange', param_stats_2p['yrange']) # H68
            print('param_stats_2p.global_min_chisq', param_stats_2p['global_min_chisq'], '1/x', 1/param_stats_2p['global_min_chisq'])
            print('param_stats_2p.global_best', param_stats_2p['global_best'])
            print('param_stats_2p.in_range_min_chisq', param_stats_2p['in_range_min_chisq'])
            print('param_stats_2p.in_range_best', param_stats_2p['in_range_best'])
            print('param_stats_2p.in_range_min', param_stats_2p['in_range_min'])
            print('param_stats_2p.in_range_max', param_stats_2p['in_range_max'])
            print('param_stats_2p.min', param_stats_2p['min'])
            print('param_stats_2p.max', param_stats_2p['max'])
            print('param_stats_2p.median', param_stats_2p['median'])
            print('param_stats_2p.best', param_stats_2p['best'])
            print('param_stats_2p.sigma', param_stats_2p['sigma'])
            print('param_stats_2p.L68', param_stats_2p['L68'])
            print('param_stats_2p.H68', param_stats_2p['H68'])
        #--
        #--TODO--20180123-10h44m-- when param_log is True, param_min can be zero!
        #--
        # 
        # 20180319: prevent param_stats_2p['xrange'] from being too small (1/20.0 of the plotting xrange in right panel) <TODO>
        if param_stats_2p['xrange'][1] - param_stats_2p['xrange'][0] < numpy.abs(plot_xrange[1]-plot_xrange[0])/20.0:
            param_stats_2p['xrange'][0] = param_stats_2p['median'] - numpy.abs(plot_xrange[1]-plot_xrange[0])/20.0
            param_stats_2p['xrange'][1] = param_stats_2p['median'] + numpy.abs(plot_xrange[1]-plot_xrange[0])/20.0
            if verbose >= 2:
                print('param_stats_2p.xrange', param_stats_2p['xrange']) # optimized xrange for L68-H68
                print('param_stats_2p.yrange', param_stats_2p['yrange']) # optimized xrange for L68-H68
        # 
        # Initialize a plot
        if Plot_engine is None:
            Plot_engine = CrabPlot(figure_size=(9.0,5.0))
            Plot_engine.set_margin(panel=0, top=0.96, bottom=0.04)
        # 
        # Plot xy (left panel)
        Plot_engine.plot_xy(param_array, 1/numpy.array(chisq_array), overplot = False, 
                                xtitle = param_dict['Par_name'], ytitle = r'$1/\chi^2$', useTex = True, 
                                size = 2.2, color='#1873cc', symbol = 'o')
        # 
        # Plot Cut_chi2 horizontal line
        Plot_engine.plot_line(param_stats['xrange'][0], 1/(param_stats['threshold_chisq']), 
                                param_stats['xrange'][1], 1/(param_stats['threshold_chisq']), 
                                overplot = True, color='gold', linestyle = 'dashed', linewidth = 4.0, alpha = 0.8, zorder=9)
        # 
        # Plot Cut_chi2 horizontal line (2p = 2.3)
        if param_stats_2p['valid']:
            Plot_engine.plot_line(param_stats_2p['xrange'][0], 1/(param_stats_2p['threshold_chisq']), 
                                    param_stats_2p['xrange'][1], 1/(param_stats_2p['threshold_chisq']), 
                                    overplot = True, color='orangered', linestyle = '--', dashes=(1.25,0.75), linewidth = 8.0, alpha = 0.5, zorder=9) # dashes=(0.5,0.1) means length of 0.5, space of 0.25
                                    # color: http://www.color-hex.com/color/1e90ff
                                    # https://stackoverflow.com/questions/35099130/change-spacing-of-dashes-in-dashed-line-in-matplotlib
        # 
        # Plot histogram (right panel)
        if False:
            if param_stats_2p['valid']:
                Plot_engine.plot_hist(param_bin_x, 1/numpy.array(param_bin_y), width = param_bin_step*1.5, align = 'edge', overplot = False, 
                                        xtitle = param_dict['Par_name'], ytitle = r'$1/\chi^2$', useTex = True, 
                                        xrange = plot_xrange, yrange = plot_yrange, xlog = xlog, ylog = ylog)
                Plot_engine.set_xcharsize(8.0)
                # 
                # Plot Cut_chi2 line
                Plot_engine.plot_line(plot_xrange[0], 1/(param_stats['threshold_chisq']), plot_xrange[1], 1/(param_stats['threshold_chisq']), overplot = True, linestyle = 'dashed')
                Plot_engine.plot_text(plot_xrange[1], plot_yrange[1]-0.02*(plot_yrange[1]-plot_yrange[0]), ' (zoomed) ', NormalizedCoordinate=False, overplot=True, horizontalalignment='right', verticalalignment='top')
                # 
                # Plot Cut_chi2 line (2p = 2.3)
                Plot_engine.plot_line(param_stats_2p['xrange'][0], 1/(param_stats_2p['threshold_chisq']), 
                                        param_stats_2p['xrange'][1], 1/(param_stats_2p['threshold_chisq']), 
                                        overplot = True, color='#1e90ff', linestyle = '--', dashes=(0.5,0.25), linewidth = 4.0, alpha = 0.8) # dashes=(0.5,0.1) means length of 0.5, space of 0.25
                                        # color: http://www.color-hex.com/color/1e90ff
                                        # https://stackoverflow.com/questions/35099130/change-spacing-of-dashes-in-dashed-line-in-matplotlib
            else:
                Plot_engine.plot_text(1.0-0.02, 1.00-0.02, ' (zoomed) ', NormalizedCoordinate=True, overplot=False, horizontalalignment='right', verticalalignment='top')
                Plot_engine.plot_text(0.5, 0.5, ' (No valid data) ', NormalizedCoordinate=True, overplot=True, horizontalalignment='center', verticalalignment='center')
        # 
    else:
        print('Error! analyze_chisq_distribution() got unaccepted inputs!')
        sys.exit()






def constrain_by_upper_limits(chisq_file, chisq_array, lib_dict):
    #if os.path.isfile('flagged_chi2_solution.txt'):
    #    os.system('mv flagged_chi2_solution.txt flagged_chi2_solution.txt.backup')
    #if os.path.isfile('flagged_chi2_solution_sorted_index.txt'):
    #    os.system('mv flagged_chi2_solution_sorted_index.txt flagged_chi2_solution_sorted_index.txt.backup')
    #os.system('bash -c \"rm -rf obj_* 2>/dev/null\"')
    #os.system('bash -c \"rm -rf dump_* 2>/dev/null\"')
    if os.path.isfile('extracted_flux.txt'):
        obs_data_table = asciitable.read('extracted_flux.txt')
        obs_wave = numpy.array(obs_data_table.field(obs_data_table.colnames[0]))
        obs_flux = numpy.array(obs_data_table.field(obs_data_table.colnames[1]))
        obs_error = numpy.array(obs_data_table.field(obs_data_table.colnames[2]))
        obs_undetection = (obs_flux < 3.0*obs_error) & (obs_error>0)
        obs_where_undetected = numpy.argwhere(obs_undetection)
        # 
        # check if obs data contains upper limits
        if len(obs_where_undetected) > 0:
            obs_wave_undetected = obs_wave[obs_undetection] / (1.0 + float(lib_dict['REDSHIFT']))
            obs_flux_undetected = obs_error[obs_undetection] * 5.0 # 5-sigma upper limit <TODO>
            # 
            # loop each input chi2 solution
            chisq_indices_sorted = numpy.argsort(chisq_array)
            i_constrain = 0
            i_constrain = 1436 #<TODO><DEBUG># 
            while i_constrain < len(chisq_array):
                print('constrain_by_upper_limits: dump_LIB_SEDs_to_files(%d) (%d)'%(chisq_indices_sorted[i_constrain], i_constrain))
                #Read_SED_LIB(chisq_file, chisq_array, lib_dict, chisq_indices_sorted[i_constrain])
                dump_LIB_SEDs_to_files(chisq_file = chisq_file, chisq_array = chisq_array, lib_dict = lib_dict, 
                                        dump_indices = chisq_indices_sorted[i_constrain], 
                                        output_numbers = 1, 
                                        output_prefix = 'dump_')
                SED_data_table = asciitable.read('dump_1/SED_SUM') # always read the minimum chi2 solution
                SED_x = numpy.array(SED_data_table.field(SED_data_table.colnames[0]))
                SED_y = numpy.array(SED_data_table.field(SED_data_table.colnames[1]))
                #SED_flux_to_constrain = scipy.interpolate.spline(SED_x, SED_y, obs_wave_undetected, order='1') # order=3, kind='smoothest', conds=None
                SED_flux_to_constrain = scipy.interpolate.interp1d(SED_x, SED_y, kind='nearest')(obs_wave_undetected)
                # 
                constrained = (SED_flux_to_constrain > obs_flux_undetected)
                asciitable.write(numpy.column_stack((obs_wave_undetected, obs_flux_undetected, SED_flux_to_constrain, constrained)), 
                                    sys.stdout, 
                                    Writer=asciitable.FixedWidthTwoLine, 
                                    delimiter='|', delimiter_pad=' ', bookend=True, 
                                    names=['obs_wave_undetected', 'obs_flux_undetected', 'SED_flux_to_constrain', 'constrained'])
                print('')
                # 
                where_to_constrain = numpy.argwhere(constrained)
                if len(where_to_constrain) > 0:
                    # this chi2 solution is not allowed by the upper limit
                    chisq_array[chisq_indices_sorted[i_constrain]] = 1e+99
                    os.system('echo %d >> flagged_chi2_solution.txt'%(chisq_indices_sorted[i_constrain]))
                    os.system('echo %d >> flagged_chi2_solution_sorted_index.txt'%(i_constrain))
                os.system('rm -rf dump_1') # always read the minimum chi2 solution
                if len(where_to_constrain) <= 0:
                    # ok, nothing to constrain, break
                    break
                i_constrain = i_constrain + 1
    return chisq_array
                





def random_sorted_chi2_index_dict(Cut_chi2_array, max = 50):
    # 
    # This function returns a list of index in the 'Cut_chi2_array', including the first one, the last one, and at most 50 random elements in between.
    # 
    Plot_chi2_max_number = max # Max chi2 solutions to plot, we plot the first Plot_chi2_max_number/2 and the last Plot_chi2_max_number/2 solutions, skip solutions in the middle.
    # 
    #Cut_chi2_array = numpy.random.random(30)
    Cut_chi2_array_size = len(Cut_chi2_array)
    # 
    # make sure Plot_chi2_max_number <= Cut_chi2_array_size
    if Plot_chi2_max_number > Cut_chi2_array_size:
        Plot_chi2_max_number = Cut_chi2_array_size
    # 
    # first randomly select Plot_chi2_max_number chisq solutions
    # then append minimum chisq solution to the Plot_chi2_indices
    Plot_chi2_indices = numpy.sort(numpy.argsort(numpy.random.random(Cut_chi2_array_size))[0:Plot_chi2_max_number])[::-1] # non-repeated index array from 0 to Plot_chi2_max_number
    Plot_chi2_indices[0] = Cut_chi2_array_size-1 if Plot_chi2_indices[0] != Cut_chi2_array_size-1 else Cut_chi2_array_size-1 # make sure the first element is always 'Cut_chi2_array_size-1', i.e., the worst chi-square solution
    Plot_chi2_indices[-1] = 0 if Plot_chi2_indices[-1] != 0 else 0 # make sure the last element is always '0', i.e., the mininum chi-square solution
    # 
    Plot_chi2_index_dict = {}
    for i in range(len(Plot_chi2_indices)): 
        Plot_chi2_index_dict['%d'%(Plot_chi2_indices[i])] = Cut_chi2_array[Plot_chi2_indices[i]] # Cut_chi2_array[i] #<BUG><FIXED><20181005><DZLIU># 
    # 
    return Plot_chi2_index_dict, Plot_chi2_indices





def dump_LIB_SEDs_to_files(chisq_file = '', chisq_array = [], lib_dict = {}, 
                            dump_indices = [], output_numbers = [], output_prefix = 'dump_', 
                            redshift = numpy.nan):
    # 
    #def Read_SED_LIB(DataFile, DataArray, InfoDict, chisq_indices_sorted, Cut_chi2_array_size = 1, Plot_chi2_index_dict = [])
    # 
    if numpy.isscalar(dump_indices):
        dump_indices = [dump_indices]
    # 
    if len(dump_indices) == 0:
        dump_indices = numpy.arange(len(chisq_array))
    # 
    if numpy.isscalar(output_numbers):
        output_numbers = [output_numbers]
    # 
    if len(output_numbers) == 0:
        output_numbers = numpy.arange(len(dump_indices)) + 1
    # 
    #output_prefix = 'dump_'
    # 
    for i in range(len(output_numbers)):
        # 
        # skip solutions between 11th to last 11th.
        #if i > Plot_chi2_max_number/2 and i<(dump_number-1-Plot_chi2_max_number/2):
        #    continue
        #if len(Plot_chi2_index_dict) > 0:
        #    if not ('%d'%i) in Plot_chi2_index_dict:
        #        continue
        # 
        # 
        if i >= len(dump_indices):
            continue
        # 
        # 
        # create output directory
        if not os.path.isdir('%s%d'%(output_prefix, output_numbers[i])):
            os.makedirs('%s%d'%(output_prefix, output_numbers[i]))
        # 
        # loop each SED LIB and dump LIB file (we do not overwrite existing files)
        for j in range(int(InfoDict['NLIB'])):
            # 
            if not os.path.isfile('%s%d/SED_LIB%d'%(output_prefix, output_numbers[i], j+1)):
                #BashCommand = 'cd dump/%d/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_SED ../%s %d %s SED_LIB%d'%\
                #                    (i+1, \
                #                        InfoDict['LIB%d'%(j+1)], \
                #                            DataArray['i%d'%(j+1)][dump_indices[i]], \
                #                                DataArray['a%d'%(j+1)][dump_indices[i]], \
                #                                    j+1)
                #print(BashCommand)
                #os.system(BashCommand)
                # 
                # do python way 20180113
                BashCommand = '%s/michi2_read_lib_SEDs.py %s %d %s%d > %s%d/log.txt'%\
                                ( os.path.dirname(os.path.abspath(__file__)), \
                                    chisq_file, \
                                        dump_indices[i]+1, \
                                            output_prefix, \
                                                output_numbers[i], \
                                                    output_prefix, \
                                                        output_numbers[i])
                print(BashCommand)
                os.system(BashCommand)
                # 
                # check
                #cd dump/8/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_SED ../lib.DL07.LoExCom.SED 140140 16.1932 SED_LIB4
                #/Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_SEDs.py fit_5.out 2451 c_2451
                #topcat -f ASCII c_2451/SED_LIB4 dump/8/SED_LIB4 &
                #checked that the two code give exactly the same result!
                #
                # how about the integrated IR luminosity?
                #cat dump/8/SED_LIB4.vLv_8_1000 # 8.8442327616e+03
                #cd dump/8
                #sm
                #load astroSfig.sm
                #data SED_LIB4 read {x 1 y 2}
                #calc_ltir x y # 3536.147921
                #calc 10**2.339198 * 16.1932 # 3536.150006 -- 2.339198 is the PAR3 in lib file, agreed with our manual integration! 
                # 
        # 
        # also dump chi2, line_number (in the chisq_file)
        if not os.path.isfile('%s%d/chi2.txt'%(output_prefix, output_numbers[i])):
            BashCommand = 'echo "%s" > %s%d/chi2.txt'%\
                            (chisq_array[dump_indices[i]], \
                                output_prefix, \
                                    output_numbers[i])
            print(BashCommand)
            os.system(BashCommand)
        # 
        if not os.path.isfile('%s%d/line_number.txt'%(output_prefix, output_numbers[i])):
            BashCommand = 'echo "%s" > %s%d/line_number.txt'%\
                            (dump_indices[i]+1, \
                                output_prefix, \
                                    output_numbers[i])
            print(BashCommand)
            os.system(BashCommand)
        # 
        # also dump redshift if possible
        if not numpy.isnan(redshift) and not os.path.isfile('%s%d/redshift.txt'%(output_prefix, output_numbers[i])):
            BashCommand = 'echo "%s" > %s%d/redshift.txt'%\
                            (redshift, \
                                output_prefix, \
                                    output_numbers[i])
            print(BashCommand)
            os.system(BashCommand)
        # 
        #return SED_x, SED_y




















######################################
#               MAIN PROGRAM             #
##########################################

if len(sys.argv) <= 1:
    
    print('Usage: michi2_plot_SED_fitting_results.py fit_5.out')
    sys.exit()

else:
    # 
    # Read user input
    SetOnlyPlotBestSED = False
    ConstrainByUpperLimits = False
    SourceName = ''
    PlotYRange = []
    PlotMaxSEDNumber = 15 #<TODO># 50
    UserInputFluxFile = ''
    UserOutputName = ''
    UserInputText = []
    UserInputColorForAGN = ''
    UserThickColorForAGN = 0.0
    UserInputThickForAGN = 1.5
    UserInputVerbose = 0
    iarg = 1
    while iarg < len(sys.argv):
        TempCmd = sys.argv[iarg].replace('--','-').lower()
        if TempCmd=='-only-plot-best-sed' or TempCmd=='-only-best':
            SetOnlyPlotBestSED = True
            print('Setting only plot best-fit!')
        elif TempCmd=='-constrain-by-upper-limits' or TempCmd=='-constrain':
            ConstrainByUpperLimits = True
            print('Setting constrain by upper limits!')
        elif TempCmd=='-source-name' or TempCmd=='-source' or TempCmd=='-s':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                SourceName = sys.argv[iarg]
                print('Setting SourceName = %s'%(SourceName))
        elif TempCmd=='-yrange':
            if iarg+2 < len(sys.argv):
                iarg = iarg + 1
                PlotYRange.append(float(sys.argv[iarg]))
                iarg = iarg + 1
                PlotYRange.append(float(sys.argv[iarg]))
                print('Setting PlotYRange = %s'%(PlotYRange))
        elif TempCmd=='-max-sed-number' or TempCmd=='-max':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                PlotMaxSEDNumber = int(sys.argv[iarg])
                print('Setting PlotMaxSEDNumber = %s'%(PlotMaxSEDNumber))
        elif TempCmd=='-flux-file' or TempCmd=='-flux' or TempCmd=='-f':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                UserInputFluxFile = sys.argv[iarg]
                print('Setting UserInputFluxFile = %s'%(UserInputFluxFile))
        elif TempCmd=='-output-file' or TempCmd=='-output-name' or TempCmd=='-output' or TempCmd=='-out' or TempCmd=='-o':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                UserOutputName = sys.argv[iarg]
                print('Setting UserOutputName = %s'%(UserOutputName))
        elif TempCmd=='-text':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                UserInputText.append(sys.argv[iarg])
                print('Setting UserInputText += \'%s\''%(sys.argv[iarg]))
        elif TempCmd=='-color-for-agn':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                UserInputColorForAGN = sys.argv[iarg]
                print('Setting UserInputColorForAGN = %s'%(sys.argv[iarg]))
        elif TempCmd=='-thick-for-agn':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                UserInputThickForAGN = float(sys.argv[iarg])
                print('Setting UserInputThickForAGN = %s'%(sys.argv[iarg]))
        elif TempCmd=='-verbose':
            UserInputVerbose += 1
            print('Setting UserInputVerbose = %s'%(UserInputVerbose))
        else:
            DataFile = sys.argv[iarg]
        iarg = iarg + 1
    # 
    # check if output figure already exists or not, see if we overwrite or not.
    Output_dir = ""
    Output_name, Output_extension = os.path.splitext(DataFile)
    if UserOutputName != '':
        if UserOutputName.find(os.sep) >= 0:
            Output_dir = os.path.dirname(UserOutputName) + os.sep
            Output_name = os.path.basename(UserOutputName)
            if not os.path.isdir(Output_dir):
                os.makedirs(Output_dir)
        else:
            Output_name = UserOutputName
        if Output_name.endswith('.pdf') or Output_name.endswith('.PDF') or \
           Output_name.endswith('.eps') or Output_name.endswith('.EPS'):
            Output_name, Output_extension = os.path.splitext(Output_name)
        else:
            Output_name = UserOutputName
            Output_extension = 'pdf'
    # 
    # Read chi2 table
    #DataFile = sys.argv[1]
    if DataFile == '' or not os.path.isfile(DataFile):
        print('Error! The input fitted chi2 data file "%s" was not found!'%(DataFile))
        sys.exit()
    print('# Reading "%s"'%(DataFile))
    DataTable = CrabTable(DataFile, verbose=1)
    # 
    # Read fit info
    InfoFile = DataFile + '.info'
    if not os.path.isfile(InfoFile):
        print('Error! The input fitted chi2 info file "%s" was not found!'%(InfoFile))
        sys.exit()
    print('# Reading "%s"'%(InfoFile))
    InfoDict = CrabTableReadInfo(InfoFile, verbose=0)
    CheckInfoDictOK = True
    for InfoKey in ['OBS', 'NLIB', 'OUT']:
        if InfoKey not in InfoDict:
            print('Error! Key "%s" is not in the InfoFile "%s"!'%(InfoKey, InfoFile))
            CheckInfoDictOK = False
    if CheckInfoDictOK is False:
        sys.exit()
    # 
    print(InfoDict)
    # 
    # Get Redshift
    Redshift = float(InfoDict['REDSHIFT'])
    # 
    # Get Redshift
    if 'DISTANCE' in InfoDict:
        Distance = float(InfoDict['DISTANCE'])
    else:
        Distance = -99
    # 
    # Set default OBS data file name <TODO><DEBUG>
    if UserInputFluxFile == '' and os.path.isfile('extracted_flux.txt'):
        UserInputFluxFile = 'extracted_flux.txt'
    # 
    # Read OBS data file
    if UserInputFluxFile == '':
        DataTable_obs = asciitable.read(InfoDict['OBS'])
    else:
        DataTable_obs = asciitable.read(UserInputFluxFile)
    try:
        Wavelength_obs = DataTable_obs[DataTable_obs.colnames[0]].data
        Flux_obs = DataTable_obs[DataTable_obs.colnames[1]].data
        FluxErr_obs = DataTable_obs[DataTable_obs.colnames[2]].data
        Detection_mask = (Flux_obs>=2.0*FluxErr_obs)
        UpperLimits_mask = (Flux_obs<2.0*FluxErr_obs)
        #print(Wavelength_obs)
    except Exception as err:
        print(err)
    # 
    # Fix data table header problem
    DataHeaders = []
    with open(DataFile,'r') as fp:
        while True:
            DataLine = fp.readline()
            if not DataLine:
                break
            if DataLine.startswith('#'):
                DataLineSplit = DataLine.replace('#','').strip().split()
                if len(DataLineSplit) == len(DataTable.TableHeaders):
                    DataHeaders = DataLineSplit
            else:
                break
        fp.close()
    # 
    # Read big chi2 data table columns
    DataArray = {}
    for i in range(len(DataHeaders)):
        #if DataHeaders[i] == 'chi2':
        #    Data_chi2 = DataTable.getColumn(i+1)
        for j in range(int(InfoDict['NLIB'])):
            if DataHeaders[i] == 'i%d'%(j+1):
                DataArray['i%d'%(j+1)] = DataTable.getColumn(i+1)
            elif DataHeaders[i] == 'a%d'%(j+1):
                DataArray['a%d'%(j+1)] = DataTable.getColumn(i+1)
            elif DataHeaders[i] == 'i0':
                DataArray['i0'] = DataTable.getColumn(i+1)
            elif DataHeaders[i] == 'a0':
                DataArray['a0'] = DataTable.getColumn(i+1)
            elif DataHeaders[i] == 'chi2':
                DataArray['chi2'] = DataTable.getColumn(i+1)
    # 
    # Determine degree of freedom
    DegreeOfFreedom = numpy.count_nonzero(Detection_mask) - int(InfoDict['NLIB'])
    # 
    # 
    # 
    # 
    # Constrain DataArray by upper limits <20180202>
    # <20190104> Oops, I forgot what it is.
    if ConstrainByUpperLimits:
        #DataArray['chi2'] = constrain_by_upper_limits(DataFile, DataArray['chi2'], InfoDict)
        constrained = (DataArray['a2']>3.2)
        where_to_constrain = numpy.argwhere(constrained)
        DataArray['chi2'][where_to_constrain] = 1e99
    # 
    # 
    # 
    # 
    # Sort chi2 table
    All_chi2_array = DataArray['chi2']
    All_chi2_indices_sorted = numpy.argsort(All_chi2_array)
    All_chi2_array_size = len(All_chi2_indices_sorted)
    # 
    All_chi2_minimum_value = numpy.nanmin(All_chi2_array)
    Cut_chi2_threshold = All_chi2_minimum_value + Delta_chisq_of_interest # 5 parameter, 68% confidence
    # 
    Cut_chi2_array_size = len(numpy.argwhere(All_chi2_array<=Cut_chi2_threshold)) # select how many fitting solution with chi2 <= 'Cut_chi2_threshold'
    Cut_chi2_array_size = All_chi2_array_size if All_chi2_array_size<Cut_chi2_array_size else Cut_chi2_array_size # do not exceed the total number of chi2
    Cut_chi2_array = All_chi2_array[All_chi2_indices_sorted[0:Cut_chi2_array_size]] # the cut-and-sorted chi2 array, note that when selecting subscript/index with :, the upper index is not included.
    # 
    Min_chi2 = numpy.nanmin(Cut_chi2_array)
    Max_chi2 = numpy.nanmax(Cut_chi2_array)
    # 
    Plot_chi2_linewidth = numpy.sqrt(1.44/float(Cut_chi2_array_size)) #<TODO># tune line width
    Plot_SED_linewidth = 1.0
    #print('')
    print('Selecting %d chi2 solutions with chi2 <= min(chi2)+%s'%(Cut_chi2_array_size, Delta_chisq_of_interest))
    # 
    if not SetOnlyPlotBestSED:
        if not os.path.isfile(Output_dir+'Plot_chi2_index_dict.json') or not os.path.isfile(Output_dir+'Plot_chi2_indices.json'):
            Plot_chi2_index_dict, Plot_chi2_indices = random_sorted_chi2_index_dict(Cut_chi2_array, max = PlotMaxSEDNumber) # we plot 50 chi2 solution curves
            with open(Output_dir+'Plot_chi2_index_dict.json', 'w') as fp:
                json.dump(Plot_chi2_index_dict, fp, indent=4)
                fp.close()
            with open(Output_dir+'Plot_chi2_indices.json', 'w') as fp:
                json.dump(Plot_chi2_indices.tolist(), fp, indent=4)
                fp.close()
        else:
            with open(Output_dir+'Plot_chi2_index_dict.json', 'r') as fp:
                Plot_chi2_index_dict = json.load(fp)
                fp.close()
            with open(Output_dir+'Plot_chi2_indices.json', 'r') as fp:
                Plot_chi2_indices = numpy.array(json.load(fp))
                fp.close()
    else:
        Plot_chi2_index_dict = {}
        Plot_chi2_index_dict['0'] = Cut_chi2_array[0]
        Plot_chi2_indices = [0] # note -- this is the index of Cut_chi2_array, and 'Cut_chi2_array' is the cut-and-sorted array of 'All_chi2_array'. 
        # 
        # also tune plotting linewidth
        Plot_chi2_linewidth = 1.5
        Plot_SED_linewidth = 2.0
    #print('Will plot chi2 solution indices %s'%(Plot_chi2_index_dict))
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    #########################################
    # Now prepare to plot the best-fit SEDs #
    #########################################
    # 
    # Get SED (dump to subdirectories and files)
    #Read_SED_LIB(DataFile, DataArray, InfoDict, All_chi2_indices_sorted, Cut_chi2_array_size, Plot_chi2_index_dict)
    #print('')
    print('Dumping library SEDs')
    dump_LIB_SEDs_to_files(chisq_file = DataFile, chisq_array = All_chi2_array, 
                            lib_dict = InfoDict, 
                            dump_indices = All_chi2_indices_sorted[Plot_chi2_indices], 
                            output_numbers = numpy.array(Plot_chi2_indices)+1, 
                            output_prefix = Output_dir+'dump'+os.sep, 
                            redshift = Redshift)
    # 
    # Wait for a long time
    # 
    # set color styles
    Color_list = ['cyan', 'gold', 'red', 'blue', 'purple', 'blue']
    Color_preset = {'DL07.HiExCom': 'red', 'DL07.LoExCom': 'blue', 'Radio': 'purple', 'MullaneyAGN': 'gold'}
    Plot_engine = CrabPlot(figure_size=(8.0,5.0))
    Plot_engine.set_margin(top=0.92, bottom=0.16, left=0.12, right=0.96)
    Count_label_chi2 = 0 # to count the chi-square label printed on the figure, make sure there are not too many labels.
    Count_label_rchi2 = 0 # to count the reduced-chi-square label printed on the figure, make sure there are not too many labels.
    Count_plot_chi2 = 0
    #Count_plot_rchi2 = 0 # using Count_plot_chi2 instead
    # 
    # 
    # Then plot SEDs
    #print('')
    for i in Plot_chi2_indices:
        # 
        # alpha by chi2
        if i >= len(Cut_chi2_array): continue
        print('Plotting chi2=%s, reduced-chi2=%s, dump id %d'%(Cut_chi2_array[i], Cut_chi2_array[i]/DegreeOfFreedom, i+1))
        Min_chi2_log = numpy.log10(Min_chi2)
        Max_chi2_log = numpy.log10(Max_chi2)
        Min_chi2_for_plot = numpy.power(10, Min_chi2_log-(Max_chi2_log-Min_chi2_log)*0.8)
        Max_chi2_for_plot = numpy.power(10, Max_chi2_log+(Max_chi2_log-Min_chi2_log)*0.3)
        Plot_chi2_alpha = Plot_engine.get_color_by_value([Min_chi2_for_plot, Max_chi2_for_plot], 
                                                         input_value=Cut_chi2_array[i], 
                                                         log=1, 
                                                         cmap=matplotlib.cm.get_cmap('gray_r'))[0]
        #print('Plot_chi2_alpha: ', Plot_chi2_alpha)
        # 
        # 
        # plot each single SED component
        for j in range(int(InfoDict['NLIB'])):
            
            xclip = None
            if InfoDict['LIB%d'%(j+1)].find('FSPS')>=0: 
                xclip = [(50,numpy.inf)]
            elif InfoDict['LIB%d'%(j+1)].find('BC03')>=0: 
                xclip = [(50,numpy.inf)]
            elif InfoDict['LIB%d'%(j+1)].find('Radio')>=0: 
                xclip = [(-numpy.inf,2e3)]
            
            Plot_lib_color = Color_list[j]
            Plot_lib_thick = Plot_chi2_linewidth
            Plot_lib_alpha = Plot_chi2_alpha
            for Preset_lib_color in Color_preset:
                if InfoDict['LIB%d'%(j+1)].find(Preset_lib_color)>=0: 
                    Plot_lib_color = Color_preset[Preset_lib_color]
                if InfoDict['LIB%d'%(j+1)].find('AGN')>=0 and UserInputColorForAGN != '': 
                    Plot_lib_color = UserInputColorForAGN # <20180515> allow user to set AGN component color
                if InfoDict['LIB%d'%(j+1)].find('AGN')>=0 and UserInputThickForAGN > 0.0: 
                    Plot_lib_thick = UserInputThickForAGN # <20180515> allow user to set AGN component thickness
            
            Plot_engine.plot_data_file(Output_dir+'dump'+os.sep+'%d/SED_LIB%d'%(i+1,j+1), xlog=1, ylog=1, xclip=xclip, current=1, \
                                dataname=Output_dir+'dump'+os.sep+'%d_SED_LIB%d'%(i+1,j+1), 
                                redshift = Redshift, 
                                linestyle='dashed', linewidth=Plot_lib_thick, color=Plot_lib_color, alpha=Plot_lib_alpha)
        # 
        Min_chi2_for_plot = numpy.power(10, Min_chi2_log-(Max_chi2_log-Min_chi2_log)*0.05)
        Max_chi2_for_plot = numpy.power(10, Max_chi2_log+(Max_chi2_log-Min_chi2_log)*0.85)
        Plot_chi2_color = Plot_engine.get_color_by_value([Min_chi2_for_plot, Max_chi2_for_plot], 
                                                    input_value=Cut_chi2_array[i], 
                                                    log=1, 
                                                    cmap=matplotlib.cm.get_cmap('gray'))
        #print('Plot_chi2_color: ', Plot_chi2_color)
        # 
        # 
        # plot total SED
        Plot_engine.plot_data_file(Output_dir+'dump'+os.sep+'%d/SED_SUM'%(i+1), xlog=1, ylog=1, current=1, \
                            dataname=Output_dir+'dump'+os.sep+'%d_SED_SUM'%(i+1), 
                            redshift = Redshift, 
                            linestyle='solid', linewidth=Plot_SED_linewidth, 
                            color=Plot_chi2_color, alpha=1.0, zorder=8) # alpha=Plot_chi2_alpha
        # 
        # 
        # count++
        Count_plot_chi2 = Count_plot_chi2 + 1
        # 
        # 
        # show chi2 text on the figure
        if not SetOnlyPlotBestSED:
            if i == Cut_chi2_array_size-1:
                Plot_engine.xyouts(0.05, 0.95, r'$\chi^2:$', NormalizedCoordinate=True, useTex=True)
            if i == 0:
                Plot_engine.xyouts(0.09, 0.95-0.03*(Count_label_chi2), '......', NormalizedCoordinate=True, color=Plot_chi2_color) # i == 0 is the minimum chisq
                Count_label_chi2 = Count_label_chi2 + 1
            if Count_plot_chi2 % int((Cut_chi2_array_size/7)+1) == 0 or i == 0 or i == Cut_chi2_array_size-1:
                #print('Plotting label at', 0.09, 0.95-0.03*(Cut_chi2_array_size-1-i), 'chi2 = %.1f'%(Cut_chi2_array[i]))
                Plot_engine.xyouts(0.09, 0.95-0.03*(Count_label_chi2), '%.1f'%(Cut_chi2_array[i]), NormalizedCoordinate=True, useTex=True, color=Plot_chi2_color)
                Count_label_chi2 = Count_label_chi2 + 1
        # 
        # 
        # show reduced-chi2 text on the figure
        if not SetOnlyPlotBestSED:
            if i == Cut_chi2_array_size-1:
                Plot_engine.xyouts(0.05+0.10, 0.95, r'$\chi_{r.}^2:$', NormalizedCoordinate=True, useTex=True)
            if i == 0:
                Plot_engine.xyouts(0.09+0.10, 0.95-0.03*(Count_label_rchi2), '......', NormalizedCoordinate=True, color=Plot_chi2_color) # i == 0 is the minimum chisq
                Count_label_rchi2 = Count_label_rchi2 + 1
            if Count_plot_chi2 % int((Cut_chi2_array_size/7)+1) == 0 or i == 0 or i == Cut_chi2_array_size-1:
                #print('Plotting label at', 0.09, 0.95-0.03*(Cut_chi2_array_size-1-i), 'reduced-chi2 = %.12f'%(Cut_chi2_array[i]/DegreeOfFreedom))
                Plot_engine.xyouts(0.09+0.10, 0.95-0.03*(Count_label_rchi2), '%.2f'%(Cut_chi2_array[i]/DegreeOfFreedom), NormalizedCoordinate=True, useTex=True, color=Plot_chi2_color)
                Count_label_rchi2 = Count_label_rchi2 + 1
        # 
        # 
        # show redshift (z) and source name on the figure
        if not SetOnlyPlotBestSED:
            # if plot a range of solutions with chi-square lower than then minimum_chisq + delta_chisq
            if i == 0:
                # 
                PlotTextPosY = 0.95
                # 
                Plot_engine.xyouts(0.05+0.10+0.10, 0.95, r'$z=%s$'%(Redshift), NormalizedCoordinate=True, useTex=True)
                PlotTextPosY = 0.95
                # 
                if SourceName != '':
                    Plot_engine.xyouts(0.97, 0.90, SourceName, NormalizedCoordinate=True, fontsize=16, horizontalalignment='right')
                    PlotTextPosY = 0.90 # line spacing 0.05
                # 
                #<20180216># allow user input text with the "-text" argument
                if len(UserInputText) > 0:
                    for UserInputTextIndex in range(len(UserInputText)):
                        UserInputTextUseTeX = (UserInputText[UserInputTextIndex].find('$')>=0)
                        Plot_engine.xyouts(0.05, PlotTextPosY-0.05*(UserInputTextIndex+1), UserInputText[UserInputTextIndex], NormalizedCoordinate=True, useTex=UserInputTextUseTeX, fontsize=15, horizontalalignment='right')
        else:
            # if only plot the best solution, then we plot no chi2 on the figure, but only source name
            if i == 0:
                # 
                PlotTextPosY = 0.90
                # 
                if SourceName != '':
                    Plot_engine.xyouts(0.05, 0.90, SourceName, NormalizedCoordinate=True, fontsize=15)
                    PlotTextPosY = 0.825 # line spacing 0.075
                # 
                Plot_engine.xyouts(0.05, PlotTextPosY, '$z=%s$'%(Redshift), NormalizedCoordinate=True, useTex=True, fontsize=15)
                #<20180216># allow user input text with the "-text" argument
                if len(UserInputText) > 0:
                    for UserInputTextIndex in range(len(UserInputText)):
                        UserInputTextUseTeX = (UserInputText[UserInputTextIndex].find('$')>=0)
                        Plot_engine.xyouts(0.05, PlotTextPosY-0.075*(UserInputTextIndex+1), UserInputText[UserInputTextIndex], NormalizedCoordinate=True, useTex=UserInputTextUseTeX, fontsize=15)
            
            # when only plotting the best solution, we also read the "dump/1/SED_SUM" (wavelength is restframe) and convert it to obsframe (only wavelength) and save it as an output. 
            if i == 0:
                best_SED_data_table = asciitable.read(Output_dir+'dump'+os.sep+'%d/SED_SUM'%(i+1))
                with open(Output_dir+'dump'+os.sep+'%d/redshift.txt'%(i+1), 'r') as ifp:
                    best_SED_redshift = float(ifp.read())
                with open(Output_dir+'best-fit_SED_%s.txt'%(SourceName), 'w') as ofp:
                    ofp.write('# %25s %25s\n'%('wavelength_obsframe_um', 'flux_obsframe_mJy'))
                    for isedrow in range(len(best_SED_data_table)):
                        ofp.write('%27.10f %25.10e\n'%(best_SED_data_table.field(0)[isedrow]*(1.0+best_SED_redshift), best_SED_data_table.field(1)[isedrow]))
                    print('Output to "%s"!'%(Output_dir+'best-fit_SED_%s.txt'%(SourceName)))
            
        #break
    # 
    # 
    # 
    # 
    # 
    # 
    # Plot OBS data points
    try:
        Plot_engine.plot_xy(Wavelength_obs[Detection_mask], Flux_obs[Detection_mask], yerr=FluxErr_obs[Detection_mask], dataname='obs', overplot=True, symbol='o', symsize=2, thick=1.5, capsize=4, alpha=0.6, zorder=10, verbose=UserInputVerbose)
        Plot_engine.plot_xy(Wavelength_obs[UpperLimits_mask], 3.0*FluxErr_obs[UpperLimits_mask], dataname='upper limits', overplot=True, symbol='upper limits', symsize=3, thick=1.25, alpha=0.5, zorder=9, verbose=UserInputVerbose)
    except Exception as err:
        print(err)
    # 
    # 
    Plot_engine.set_xrange([0.1, 1e6])
    Plot_engine.set_yrange([1e-6, 1e4])
    if len(PlotYRange) == 2:
        Plot_engine.set_yrange(PlotYRange)
    else:
        if Redshift < 0.1:
            Plot_engine.set_yrange([0.8e-3, 1e6])
        if Redshift < 0.05:
            Plot_engine.set_yrange([0.8e-2, 1e7])
        if Redshift < 0.01:
            Plot_engine.set_yrange([0.8e-1, 1e8])
        #if Redshift < 0.003:
        #    Plot_engine.set_yrange([1e3,1e13])
    Plot_engine.set_xtitle(r'Observed-frame wavelength [$\mu$m]')
    Plot_engine.set_ytitle(r'Flux density [mJy]')
    Plot_engine.set_xcharsize(charsize=12, axislabelcharsize=16)
    Plot_engine.set_ycharsize(charsize=12, axislabelcharsize=16)
    Plot_engine.savepdf(Output_dir+Output_name+'.pdf')
    #Plot_engine.show()
    Plot_engine.close()
    print('Output to "%s"!'%(Output_dir+Output_name+'.pdf'))
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    ###############################################################
    # Now plot another figure of chi-square distribution analysis #
    ###############################################################
    # 
    # Now we analyze the following quantities in the big chi-square data table 'DataTable'
    # -- stellar mass
    Stellar_mass_dict = {}
    Stellar_age_dict = {}
    Stellar_EBV_dict = {}
    LTIR_warm_dust_dict = {}
    LTIR_cold_dust_dict = {}
    LFIR_warm_dust_dict = {}
    LFIR_cold_dust_dict = {}
    LFIR122_warm_dust_dict = {}
    LFIR122_cold_dust_dict = {}
    Umin_warm_dust_dict = {}
    Umin_cold_dust_dict = {}
    Mass_warm_dust_dict = {}
    Mass_cold_dust_dict = {}
    Lumin_AGN_dict = {}
    LTIR_total_dust_dict = {}
    LFIR_total_dust_dict = {}
    LFIR122_total_dust_dict = {}
    Mass_total_dust_dict = {}
    fPDR_total_dust_dict = {}
    Umean_total_dust_dict = {}
    dust_emissivity_beta_warm_dust_dict = {}
    dust_emissivity_beta_cold_dust_dict = {}
    dust_emissivity_beta_dict = {}
    dust_temperature_warm_dust_dict = {}
    dust_temperature_cold_dust_dict = {}
    dust_temperature_dict = {}
    # 
    # define constants
    pi = numpy.pi
    dL = 1.0/(4*pi) # if Redshift is not given, then we do not apply 4*pi*dL**2 to the quantities below. 
    #Redshift = float(InfoDict['REDSHIFT'])
    print('')
    print('')
    print('z = %s'%(Redshift))
    if Redshift > 0.0:
        # 
        # compute lumdist
        dL = cosmo.luminosity_distance(Redshift).value
        # 
        # compute lumdist
        #lumdist_command = '%s/lumdist'%(os.path.dirname(os.path.realpath(__file__)))
        #if not os.path.isfile(lumdist_command):
        #    print('Error! "lumdist" command was not found at %s! Maybe you have not fully downloaded the code from "https://github.com/1054/Crab.Toolkit.michi2"?')
        #    print('We need that to compute the luminosity distance!')
        #    sys.exit()
        #dL_command = '%s -simple %s'%(lumdist_command, Redshift)
        #dL_str = os.popen(dL_command).read()
        ##print(dL_str)
        #if dL_str == '':
        #    print('Error! Failed to run "lumdist" in the Terminal! Sorry! TODO: You can set "dL = ..." in "%s".'%(InfoFile))
        #    print('We need that to compute the luminosity distance!')
        #    sys.exit()
        #dL = float(dL_str)
        print('dL = %s [Mpc] (lumdist)'%(dL))
    if Distance > 0.0:
        dL = Distance
        print('dL = %s [Mpc] (set by user)'%(dL))
    # 
    # get parameter list for each lib
    Lib_number = int(InfoDict['NLIB'])
    Lib_params = {}
    Num_params = {}
    Col_number = 2+2*Lib_number+1 # this quantity indicates the column number in the big chi-square data table. The first two columns are always i0 and chi2, so the LIB params start from column 3.
    for j in range(Lib_number):
        Lib_name = 'LIB%d'%(j+1)
        Lib_dict = CrabTableReadInfo(InfoDict[Lib_name], verbose=0)
        #print(Lib_dict)
        # 
        # read the number of parameters from the SED LIB file
        Key_NPAR = '# NPAR'
        if Key_NPAR in Lib_dict:
            Num_params[Lib_name] = int(Lib_dict[Key_NPAR])
        else:
            print("Error! \"%s\" was not found in \"%s\"!"%(Key_NPAR, InfoDict[Lib_name]))
            sys.exit()
        # 
        # read the title of each parameter from the SED LIB file
        Lib_params[Lib_name] = []
        for k in range(Num_params[Lib_name]):
            Key_TPAR = '# TPAR%d'%(k+1)
            if Key_TPAR in Lib_dict:
                Lib_params[Lib_name].append(Lib_dict[Key_TPAR])
            else:
                print("Error! \"%s\" was not found in \"%s\"!"%(Key_TPAR, InfoDict[Lib_name]))
                sys.exit()
            # 
            # check stellar mass
            if InfoDict[Lib_name].find('BC03.Padova1994') >= 0 or \
                InfoDict[Lib_name].find('FSPS.CSP') >= 0:
                if 'Mass' == Lib_dict[Key_TPAR]:
                    # Stellar mass SED library Y column unit is solar luminosity per Hertz, 
                    # flux unit is mJy, 
                    # 
                    Stellar_mass_dict['Lib_file'] = InfoDict[Lib_name]
                    Stellar_mass_dict['Lib_name'] = Lib_name
                    Stellar_mass_dict['Lib_numb'] = j+1
                    Stellar_mass_dict['Par_name'] = r'$\log_{10} \ M_{*} \ \mathrm{[\mathrm{M}_{\odot}]}$' # Lib_dict[Key_TPAR]
                    Stellar_mass_dict['Par_file'] = 'Mstar'
                    Stellar_mass_dict['Col_numb'] = Col_number
                    Stellar_mass_dict['Log_calc'] = True
                    Stellar_mass_dict['range'] = numpy.power(10,[7.0,13.5])
                    Stellar_mass_dict['value'] = DataArray['a%d'%(j+1)] / (3.839e33*1e26/(4*pi*dL**2*9.52140e48)) * DataTable.getColumn(Col_number) / (1+Redshift)
                    Stellar_mass_dict['chisq'] = DataArray['chi2']
                    Stellar_mass_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                elif 'Age' == Lib_dict[Key_TPAR]:
                    Stellar_age_dict['Lib_file'] = InfoDict[Lib_name]
                    Stellar_age_dict['Lib_name'] = Lib_name
                    Stellar_age_dict['Lib_numb'] = j+1
                    Stellar_age_dict['Par_name'] = r'$\log_{10} \ \mathrm{Age} \ \mathrm{[Gyr]}$' # Lib_dict[Key_TPAR]
                    Stellar_age_dict['Par_file'] = 'Age'
                    Stellar_age_dict['Col_numb'] = Col_number
                    Stellar_age_dict['Log_calc'] = True
                    Stellar_age_dict['range'] = [0.05, 15.0]
                    Stellar_age_dict['value'] = DataTable.getColumn(Col_number)
                    Stellar_age_dict['chisq'] = DataArray['chi2']
                    Stellar_age_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                elif 'EBV' == Lib_dict[Key_TPAR]:
                    Stellar_EBV_dict['Lib_file'] = InfoDict[Lib_name]
                    Stellar_EBV_dict['Lib_name'] = Lib_name
                    Stellar_EBV_dict['Lib_numb'] = j+1
                    Stellar_EBV_dict['Par_name'] = 'E(B-V)' # Lib_dict[Key_TPAR]
                    Stellar_EBV_dict['Par_file'] = 'EBV'
                    Stellar_EBV_dict['Col_numb'] = Col_number
                    Stellar_EBV_dict['Log_calc'] = False
                    Stellar_EBV_dict['range'] = [0.0, 1.2]
                    Stellar_EBV_dict['value'] = DataTable.getColumn(Col_number)
                    Stellar_EBV_dict['chisq'] = DataArray['chi2']
                    Stellar_EBV_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
            # 
            # check dust properties
            elif InfoDict[Lib_name].find('DL07.') >= 0 and InfoDict[Lib_name].find('.HiExCom') > 0:
                if 'lgLTIR' == Lib_dict[Key_TPAR]:
                    LTIR_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LTIR_warm_dust_dict['Lib_name'] = Lib_name
                    LTIR_warm_dust_dict['Lib_numb'] = j+1
                    LTIR_warm_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{IR}}$ (warm) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LTIR_warm_dust_dict['Par_file'] = 'LIR_warm'
                    LTIR_warm_dust_dict['Col_numb'] = Col_number
                    LTIR_warm_dust_dict['Log_calc'] = True
                    LTIR_warm_dust_dict['range'] = numpy.power(10,[6.0,14.5])
                    LTIR_warm_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn(Col_number)) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LTIR_warm_dust_dict['chisq'] = DataArray['chi2']
                    LTIR_warm_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                elif 'lgLFIR' == Lib_dict[Key_TPAR]:
                    LFIR_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LFIR_warm_dust_dict['Lib_name'] = Lib_name
                    LFIR_warm_dust_dict['Lib_numb'] = j+1
                    LFIR_warm_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{FIR,\,40-400{\mu}\mathrm{m}}}$ (warm) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LFIR_warm_dust_dict['Par_file'] = 'LFIR_warm'
                    LFIR_warm_dust_dict['Col_numb'] = Col_number
                    LFIR_warm_dust_dict['Log_calc'] = True
                    LFIR_warm_dust_dict['range'] = numpy.power(10,[6.0,14.5])
                    LFIR_warm_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn(Col_number)) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LFIR_warm_dust_dict['chisq'] = DataArray['chi2']
                    LFIR_warm_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                elif 'lgLFIR122' == Lib_dict[Key_TPAR]:
                    LFIR122_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LFIR122_warm_dust_dict['Lib_name'] = Lib_name
                    LFIR122_warm_dust_dict['Lib_numb'] = j+1
                    LFIR122_warm_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{FIR,\,40-122{\mu}\mathrm{m}}}$ (warm) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LFIR122_warm_dust_dict['Par_file'] = 'LFIR122_warm'
                    LFIR122_warm_dust_dict['Col_numb'] = Col_number
                    LFIR122_warm_dust_dict['Log_calc'] = True
                    LFIR122_warm_dust_dict['range'] = numpy.power(10,[6.0,14.5])
                    LFIR122_warm_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn(Col_number)) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LFIR122_warm_dust_dict['chisq'] = DataArray['chi2']
                    LFIR122_warm_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                elif 'Umin' == Lib_dict[Key_TPAR]:
                    Umin_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Umin_warm_dust_dict['Lib_name'] = Lib_name
                    Umin_warm_dust_dict['Lib_numb'] = j+1
                    Umin_warm_dust_dict['Par_name'] = '$U_{\mathrm{min}}$ (warm)' # Lib_dict[Key_TPAR]
                    Umin_warm_dust_dict['Par_file'] = 'Umin_warm'
                    Umin_warm_dust_dict['Col_numb'] = Col_number
                    Umin_warm_dust_dict['Log_plot'] = True # 'Log_plot', plot X axis in log scale
                    Umin_warm_dust_dict['range'] = [0.08,50.0]
                    Umin_warm_dust_dict['value'] = DataTable.getColumn(Col_number)
                    Umin_warm_dust_dict['chisq'] = DataArray['chi2']
                    Umin_warm_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                    Mass_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Mass_warm_dust_dict['Lib_name'] = Lib_name
                    Mass_warm_dust_dict['Lib_numb'] = j+1
                    Mass_warm_dust_dict['Par_name'] = r'$\log_{10} \ M_{\mathrm{dust}}$ (warm) [$\mathrm{M}_{\odot}$]'
                    Mass_warm_dust_dict['Par_file'] = 'Mdust_warm'
                    Mass_warm_dust_dict['Col_numb'] = 2+2*(j+1)
                    Mass_warm_dust_dict['Log_calc'] = True
                    Mass_warm_dust_dict['range'] = numpy.power(10,[5.0,12.0])
                    Mass_warm_dust_dict['value'] = DataArray['a%d'%(j+1)] * dL**2 / (1+Redshift) # Mdust #NOTE# no need to multiply a '4*pi'!
                    Mass_warm_dust_dict['chisq'] = DataArray['chi2']
                    Mass_warm_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                    if True:
                        if InfoDict[Lib_name].find('DL07.2010.03.18') > 0: 
                            Umin_warm_dust_dict['range'][1] = 90.0 #<20180319>#
                    # 
            elif InfoDict[Lib_name].find('DL07.') >= 0 and InfoDict[Lib_name].find('.LoExCom') > 0:
                if 'lgLTIR' == Lib_dict[Key_TPAR]:
                    LTIR_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LTIR_cold_dust_dict['Lib_name'] = Lib_name
                    LTIR_cold_dust_dict['Lib_numb'] = j+1
                    LTIR_cold_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{IR}}$ (cold) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LTIR_cold_dust_dict['Par_file'] = 'LIR_cold'
                    LTIR_cold_dust_dict['Col_numb'] = Col_number
                    LTIR_cold_dust_dict['Log_calc'] = True
                    LTIR_cold_dust_dict['range'] = numpy.power(10,[6.0,14.5])
                    LTIR_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn(Col_number)) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LTIR_cold_dust_dict['chisq'] = DataArray['chi2']
                    LTIR_cold_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                elif 'lgLFIR' == Lib_dict[Key_TPAR]:
                    LFIR_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LFIR_cold_dust_dict['Lib_name'] = Lib_name
                    LFIR_cold_dust_dict['Lib_numb'] = j+1
                    LFIR_cold_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{FIR,\,40-400{\mu}\mathrm{m}}}$ (cold) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LFIR_cold_dust_dict['Par_file'] = 'LFIR_cold'
                    LFIR_cold_dust_dict['Col_numb'] = Col_number
                    LFIR_cold_dust_dict['Log_calc'] = True
                    LFIR_cold_dust_dict['range'] = numpy.power(10,[6.0,14.5])
                    LFIR_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn(Col_number)) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LFIR_cold_dust_dict['chisq'] = DataArray['chi2']
                    LFIR_cold_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                elif 'lgLFIR122' == Lib_dict[Key_TPAR]:
                    LFIR122_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LFIR122_cold_dust_dict['Lib_name'] = Lib_name
                    LFIR122_cold_dust_dict['Lib_numb'] = j+1
                    LFIR122_cold_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{FIR,\,40-122{\mu}\mathrm{m}}}$ (cold) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LFIR122_cold_dust_dict['Par_file'] = 'LFIR122_cold'
                    LFIR122_cold_dust_dict['Col_numb'] = Col_number
                    LFIR122_cold_dust_dict['Log_calc'] = True
                    LFIR122_cold_dust_dict['range'] = numpy.power(10,[6.0,14.5])
                    LFIR122_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn(Col_number)) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LFIR122_cold_dust_dict['chisq'] = DataArray['chi2']
                    LFIR122_cold_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                elif 'Umin' == Lib_dict[Key_TPAR]:
                    Umin_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Umin_cold_dust_dict['Lib_name'] = Lib_name
                    Umin_cold_dust_dict['Lib_numb'] = j+1
                    Umin_cold_dust_dict['Par_name'] = r'$U_{\mathrm{min}}$ (cold)' # Lib_dict[Key_TPAR]
                    Umin_cold_dust_dict['Par_file'] = 'Umin_cold'
                    Umin_cold_dust_dict['Col_numb'] = Col_number
                    Umin_cold_dust_dict['Log_plot'] = True # 'Log_plot', plot X axis in log scale
                    Umin_cold_dust_dict['range'] = [0.08,50.0]
                    Umin_cold_dust_dict['value'] = DataTable.getColumn(Col_number)
                    Umin_cold_dust_dict['chisq'] = DataArray['chi2']
                    Umin_cold_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                    Mass_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Mass_cold_dust_dict['Lib_name'] = Lib_name
                    Mass_cold_dust_dict['Lib_numb'] = j+1
                    Mass_cold_dust_dict['Par_name'] = r'$\log_{10} \ M_{\mathrm{dust}}$ (cold) [$\mathrm{M}_{\odot}$]'
                    Mass_cold_dust_dict['Par_file'] = 'Mdust_cold'
                    Mass_cold_dust_dict['Col_numb'] = 2+2*(j+1)
                    Mass_cold_dust_dict['Log_calc'] = True
                    Mass_cold_dust_dict['range'] = numpy.power(10,[5.0,12.0])
                    Mass_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * dL**2 / (1+Redshift) # Mdust # Mdust #NOTE# no need to multiply a '4*pi'!
                    Mass_cold_dust_dict['chisq'] = DataArray['chi2']
                    Mass_cold_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                    if True:
                        if InfoDict[Lib_name].find('DL07.2010.03.18') > 0: 
                            Umin_cold_dust_dict['range'][1] = 90.0 #<20180319>#
            # 
            # check AGN properties
            elif InfoDict[Lib_name].find('MullaneyAGN') >= 0:
                # Mullaney AGN
                #   by integrating the AGN template (AGN_TYPE=2) from 1um to 1000um, we get an integration of 5133.913101
                if 'AGN_TYPE' == Lib_dict[Key_TPAR].upper():
                    Lumin_AGN_dict['Lib_file'] = InfoDict[Lib_name]
                    Lumin_AGN_dict['Lib_name'] = Lib_name
                    Lumin_AGN_dict['Lib_numb'] = j+1
                    Lumin_AGN_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{AGN}}$ [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    Lumin_AGN_dict['Par_file'] = 'LAGN'
                    Lumin_AGN_dict['Col_numb'] = Col_number
                    Lumin_AGN_dict['Log_calc'] = True
                    Lumin_AGN_dict['range'] = numpy.power(10,[0.0,14.5])
                    Lumin_AGN_dict['value'] = DataArray['a%d'%(j+1)] * 5133.913101 * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    Lumin_AGN_dict['chisq'] = DataArray['chi2']
                    Lumin_AGN_dict['Degree_of_freedom'] = DegreeOfFreedom
            
            elif InfoDict[Lib_name].find('SiebenmorgenAGN') >= 0:
                # Siebenmorgen AGN
                #   by integrating the AGN template from 1um to 1000um, we get an integration of Lbol
                if 'AGN_TYPE' == Lib_dict[Key_TPAR].upper():
                    Lumin_AGN_dict['Lib_file'] = InfoDict[Lib_name]
                    Lumin_AGN_dict['Lib_name'] = Lib_name
                    Lumin_AGN_dict['Lib_numb'] = j+1
                    Lumin_AGN_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{AGN}}$ [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    Lumin_AGN_dict['Par_file'] = 'LAGN'
                    Lumin_AGN_dict['Col_numb'] = Col_number
                    Lumin_AGN_dict['Log_calc'] = True
                    Lumin_AGN_dict['range'] = numpy.power(10,[0.0,14.5])
                    #Lumin_AGN_dict['value'] = DataArray['a%d'%(j+1)] * 5133.913101 * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    Lumin_AGN_dict['value'] = DataArray['a%d'%(j+1)] * DataTable.getColumn(Col_number+6) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    Lumin_AGN_dict['chisq'] = DataArray['chi2']
                    Lumin_AGN_dict['Degree_of_freedom'] = DegreeOfFreedom
            
            # 
            # check dust modified blackbody properties
            elif InfoDict[Lib_name].find('.MBB.') >= 0:
                # MBB dust
                if 'BETA' == Lib_dict[Key_TPAR].upper():
                    has_two_same_lib = False
                    if 'LIB%d'%(j+1+1) in InfoDict:
                        # check if there are contiguous two MBB LIBs, if yes, then check whether current one is in the first place
                        if InfoDict['LIB%d'%(j+1+1)] == InfoDict[Lib_name]:
                            dust_emissivity_beta_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                            dust_emissivity_beta_warm_dust_dict['Lib_name'] = Lib_name
                            dust_emissivity_beta_warm_dust_dict['Lib_numb'] = j+1
                            dust_emissivity_beta_warm_dust_dict['Par_name'] = r'$\beta$ (warm)'
                            dust_emissivity_beta_warm_dust_dict['Par_file'] = 'beta_warm'
                            dust_emissivity_beta_warm_dust_dict['Col_numb'] = Col_number
                            dust_emissivity_beta_warm_dust_dict['range'] = [1.0,3.5]
                            dust_emissivity_beta_warm_dust_dict['value'] = numpy.array(DataTable.getColumn(Col_number)).astype(float) # 
                            dust_emissivity_beta_warm_dust_dict['chisq'] = DataArray['chi2']
                            dust_emissivity_beta_warm_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                            has_two_same_lib = True
                    if 'LIB%d'%(j+1-1) in InfoDict:
                        # check if there are contiguous two MBB LIBs, if yes, then check whether current one is in the second place
                        if InfoDict['LIB%d'%(j+1-1)] == InfoDict[Lib_name]:
                            # if previous LIB is the same dust.MBB.SED
                            dust_emissivity_beta_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                            dust_emissivity_beta_cold_dust_dict['Lib_name'] = Lib_name
                            dust_emissivity_beta_cold_dust_dict['Lib_numb'] = j+1
                            dust_emissivity_beta_cold_dust_dict['Par_name'] = r'$\beta$ (cold)'
                            dust_emissivity_beta_cold_dust_dict['Par_file'] = 'beta_cold'
                            dust_emissivity_beta_cold_dust_dict['Col_numb'] = Col_number
                            dust_emissivity_beta_cold_dust_dict['range'] = [1.0,3.5]
                            dust_emissivity_beta_cold_dust_dict['value'] = numpy.array(DataTable.getColumn(Col_number)).astype(float) # 
                            dust_emissivity_beta_cold_dust_dict['chisq'] = DataArray['chi2']
                            dust_emissivity_beta_cold_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                            has_two_same_lib = True
                    # 
                    if has_two_same_lib == False:
                        dust_emissivity_beta_dict['Lib_file'] = InfoDict[Lib_name]
                        dust_emissivity_beta_dict['Lib_name'] = Lib_name
                        dust_emissivity_beta_dict['Lib_numb'] = j+1
                        dust_emissivity_beta_dict['Par_name'] = r'$\beta$'
                        dust_emissivity_beta_dict['Par_file'] = 'beta'
                        dust_emissivity_beta_dict['Col_numb'] = Col_number
                        dust_emissivity_beta_dict['range'] = [1.0,3.5]
                        dust_emissivity_beta_dict['value'] = numpy.array(DataTable.getColumn(Col_number)).astype(float) # 
                        dust_emissivity_beta_dict['chisq'] = DataArray['chi2']
                        dust_emissivity_beta_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                if 'T_DUST' == Lib_dict[Key_TPAR].upper():
                    has_two_same_lib = False
                    if 'LIB%d'%(j+1+1) in InfoDict:
                        # check if there are contiguous two MBB LIBs, if yes, then check whether current one is in the first place
                        if InfoDict['LIB%d'%(j+1+1)] == InfoDict[Lib_name]:
                            dust_temperature_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                            dust_temperature_warm_dust_dict['Lib_name'] = Lib_name
                            dust_temperature_warm_dust_dict['Lib_numb'] = j+1
                            dust_temperature_warm_dust_dict['Par_name'] = r'$T_{dust}$ (warm)'
                            dust_temperature_warm_dust_dict['Par_file'] = 'T_dust_warm'
                            dust_temperature_warm_dust_dict['Col_numb'] = Col_number
                            dust_temperature_warm_dust_dict['range'] = [10.0,100.0]
                            dust_temperature_warm_dust_dict['value'] = numpy.array(DataTable.getColumn(Col_number)).astype(float) # 
                            dust_temperature_warm_dust_dict['chisq'] = DataArray['chi2']
                            dust_temperature_warm_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                            has_two_same_lib = True
                    if 'LIB%d'%(j+1-1) in InfoDict:
                        # check if there are contiguous two MBB LIBs, if yes, then check whether current one is in the second place
                        if InfoDict['LIB%d'%(j+1-1)] == InfoDict[Lib_name]:
                            dust_temperature_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                            dust_temperature_cold_dust_dict['Lib_name'] = Lib_name
                            dust_temperature_cold_dust_dict['Lib_numb'] = j+1
                            dust_temperature_cold_dust_dict['Par_name'] = r'$T_{dust}$ (cold)'
                            dust_temperature_cold_dust_dict['Par_file'] = 'T_dust_cold'
                            dust_temperature_cold_dust_dict['Col_numb'] = Col_number
                            dust_temperature_cold_dust_dict['range'] = [10.0,100.0]
                            dust_temperature_cold_dust_dict['value'] = numpy.array(DataTable.getColumn(Col_number)).astype(float) # 
                            dust_temperature_cold_dust_dict['chisq'] = DataArray['chi2']
                            dust_temperature_cold_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                            has_two_same_lib = True
                    # 
                    if has_two_same_lib == False:
                        dust_temperature_dict['Lib_file'] = InfoDict[Lib_name]
                        dust_temperature_dict['Lib_name'] = Lib_name
                        dust_temperature_dict['Lib_numb'] = j+1
                        dust_temperature_dict['Par_name'] = r'$T_{dust}$'
                        dust_temperature_dict['Par_file'] = 'T_dust'
                        dust_temperature_dict['Col_numb'] = Col_number
                        dust_temperature_dict['range'] = [10.0,100.0]
                        dust_temperature_dict['value'] = numpy.array(DataTable.getColumn(Col_number)).astype(float) # 
                        dust_temperature_dict['chisq'] = DataArray['chi2']
                        dust_temperature_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
                if 'L_DUST' == Lib_dict[Key_TPAR].upper():
                    has_two_same_lib = False
                    if 'LIB%d'%(j+1+1) in InfoDict:
                        # check if there are contiguous two MBB LIBs, if yes, then check whether current one is in the first place
                        if InfoDict['LIB%d'%(j+1+1)] == InfoDict[Lib_name]:
                            # if next LIB is the same dust.MBB.SED
                            LTIR_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                            LTIR_warm_dust_dict['Lib_name'] = Lib_name
                            LTIR_warm_dust_dict['Lib_numb'] = j+1
                            LTIR_warm_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{IR}}$ (warm) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                            LTIR_warm_dust_dict['Par_file'] = 'LIR_warm'
                            LTIR_warm_dust_dict['Col_numb'] = Col_number
                            LTIR_warm_dust_dict['Log_calc'] = True
                            LTIR_warm_dust_dict['range'] = numpy.power(10,[6.0,14.5])
                            LTIR_warm_dust_dict['value'] = DataArray['a%d'%(j+1)] * (DataTable.getColumn(Col_number)) * dL**2 / (1+Redshift) # Note: no 4*pi, see LIB.SED
                            LTIR_warm_dust_dict['chisq'] = DataArray['chi2']
                            LTIR_warm_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                            has_two_same_lib = True
                    # 
                    if 'LIB%d'%(j+1-1) in InfoDict:
                        # check if there are contiguous two MBB LIBs, if yes, then check whether current one is in the second place
                        if InfoDict['LIB%d'%(j+1-1)] == InfoDict[Lib_name]:
                            # if previous LIB is the same dust.MBB.SED
                            LTIR_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                            LTIR_cold_dust_dict['Lib_name'] = Lib_name
                            LTIR_cold_dust_dict['Lib_numb'] = j+1
                            LTIR_cold_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{IR}}$ (cold) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                            LTIR_cold_dust_dict['Par_file'] = 'LIR_cold'
                            LTIR_cold_dust_dict['Col_numb'] = Col_number
                            LTIR_cold_dust_dict['Log_calc'] = True
                            LTIR_cold_dust_dict['range'] = numpy.power(10,[6.0,14.5])
                            LTIR_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * (DataTable.getColumn(Col_number)) * dL**2 / (1+Redshift) # Note: no 4*pi, see LIB.SED
                            LTIR_cold_dust_dict['chisq'] = DataArray['chi2']
                            LTIR_cold_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                            has_two_same_lib = True
                    # 
                    if has_two_same_lib == False:
                        LTIR_total_dust_dict['Lib_file'] = InfoDict[Lib_name]
                        LTIR_total_dust_dict['Lib_name'] = Lib_name
                        LTIR_total_dust_dict['Lib_numb'] = j+1
                        LTIR_total_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{IR}}$ [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                        LTIR_total_dust_dict['Par_file'] = 'LIR_total'
                        LTIR_total_dust_dict['Col_numb'] = Col_number
                        LTIR_total_dust_dict['Log_calc'] = True
                        LTIR_total_dust_dict['range'] = numpy.power(10,[6.0,14.5])
                        LTIR_total_dust_dict['value'] = DataArray['a%d'%(j+1)] * (DataTable.getColumn(Col_number)) * dL**2 / (1+Redshift) # Note: no 4*pi, see LIB.SED
                        LTIR_total_dust_dict['chisq'] = DataArray['chi2']
                        LTIR_total_dust_dict['Degree_of_freedom'] = DegreeOfFreedom
                    # 
            # 
            # finished checking library properties
            # 
            # count the column number in the big chi-square data table. All params in each LIB are listed in the big chi-square data table.
            Col_number = Col_number + 1
    # 
    # Total LIR
    if 'value' in LTIR_warm_dust_dict and 'value' in LTIR_cold_dust_dict:
        LTIR_total_dust_dict = copy(LTIR_warm_dust_dict)
        LTIR_total_dust_dict['value'] = LTIR_warm_dust_dict['value'] + LTIR_cold_dust_dict['value']
        LTIR_total_dust_dict['Lib_file'] = [LTIR_warm_dust_dict['Lib_file'], LTIR_cold_dust_dict['Lib_file']]
        LTIR_total_dust_dict['Lib_name'] = [LTIR_warm_dust_dict['Lib_name'], LTIR_cold_dust_dict['Lib_name']]
        LTIR_total_dust_dict['Col_numb'] = [LTIR_warm_dust_dict['Col_numb'], LTIR_cold_dust_dict['Col_numb']]
        LTIR_total_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{IR}}$ (total) [$\mathrm{L}_{\odot}$]'
        LTIR_total_dust_dict['Par_file'] = 'LIR_total'
    # 
    # Total LIR 40-400um
    if 'value' in LFIR_warm_dust_dict and 'value' in LFIR_cold_dust_dict:
        LFIR_total_dust_dict = copy(LFIR_warm_dust_dict)
        LFIR_total_dust_dict['value'] = LFIR_warm_dust_dict['value'] + LFIR_cold_dust_dict['value']
        LFIR_total_dust_dict['Lib_file'] = [LFIR_warm_dust_dict['Lib_file'], LFIR_cold_dust_dict['Lib_file']]
        LFIR_total_dust_dict['Lib_name'] = [LFIR_warm_dust_dict['Lib_name'], LFIR_cold_dust_dict['Lib_name']]
        LFIR_total_dust_dict['Col_numb'] = [LFIR_warm_dust_dict['Col_numb'], LFIR_cold_dust_dict['Col_numb']]
        LFIR_total_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{FIR},\,40-400{\mu}\mathrm{m}}$ (total) [$\mathrm{L}_{\odot}$]'
        LFIR_total_dust_dict['Par_file'] = 'LIR_total_40_400'
    # 
    # Total LIR 40-122um
    if 'value' in LFIR122_warm_dust_dict and 'value' in LFIR122_cold_dust_dict:
        LFIR122_total_dust_dict = copy(LFIR122_warm_dust_dict)
        LFIR122_total_dust_dict['value'] = LFIR122_warm_dust_dict['value'] + LFIR122_cold_dust_dict['value']
        LFIR122_total_dust_dict['Lib_file'] = [LFIR122_warm_dust_dict['Lib_file'], LFIR122_cold_dust_dict['Lib_file']]
        LFIR122_total_dust_dict['Lib_name'] = [LFIR122_warm_dust_dict['Lib_name'], LFIR122_cold_dust_dict['Lib_name']]
        LFIR122_total_dust_dict['Col_numb'] = [LFIR122_warm_dust_dict['Col_numb'], LFIR122_cold_dust_dict['Col_numb']]
        LFIR122_total_dust_dict['Par_name'] = r'$\log_{10} \ L_{\mathrm{FIR},\,40-122{\mu}\mathrm{m}}$ (total) [$\mathrm{L}_{\odot}$]'
        LFIR122_total_dust_dict['Par_file'] = 'LIR_total_40_122'
    # 
    # Total Mdust
    if 'value' in Mass_warm_dust_dict and 'value' in Mass_cold_dust_dict:
        Mass_total_dust_dict = copy(Mass_warm_dust_dict)
        Mass_total_dust_dict['value'] = Mass_warm_dust_dict['value'] + Mass_cold_dust_dict['value']
        Mass_total_dust_dict['Lib_file'] = [Mass_warm_dust_dict['Lib_file'], Mass_cold_dust_dict['Lib_file']]
        Mass_total_dust_dict['Lib_name'] = [Mass_warm_dust_dict['Lib_name'], Mass_cold_dust_dict['Lib_name']]
        Mass_total_dust_dict['Col_numb'] = [Mass_warm_dust_dict['Col_numb'], Mass_cold_dust_dict['Col_numb']]
        Mass_total_dust_dict['Par_name'] = r'$\log_{10} \ M_{\mathrm{dust}}$ (total) [$\mathrm{M}_{\odot}$]'
        Mass_total_dust_dict['Par_file'] = 'Mdust_total'
    # 
    # Total fPDR Umean
    #   the Draine & Li 2007 IRSF (interstellar radiation field) model is like Dale 2005, 
    #   LoExComponent has a ISRF (or U) of Umin * G0, while
    #   HiExComponent has ISRF (or U) described by a exponential function within range of Umin and Umax, 
    #   usually Umax is fixed to 1e6 * G0, and the slope of exponential function is fixed to -2. 
    #   So the mean U value can be computed by 
    #       
    #       (1-gamma) * Umin + gamma * Umin * (ln(1e6/Umin)/(1-Umin/1e6))
    #       
    #       in which, gamma = aoe_Hi / (aoe_Lo+aoe_Hi) = a2/(a1+a2) = Mdust2/(Mdust1+Mdust2)
    #       
    #       e.g. Aniano et al. 2012
    # 
    if 'value' in Mass_warm_dust_dict and 'value' in Mass_cold_dust_dict:
        fPDR_total_dust_dict = copy(Umin_warm_dust_dict)
        fPDR_total_dust_dict['value'] = Mass_warm_dust_dict['value'] / (Mass_warm_dust_dict['value'] + Mass_cold_dust_dict['value'])
        fPDR_total_dust_dict['Lib_file'] = [Mass_warm_dust_dict['Lib_file'], Mass_cold_dust_dict['Lib_file']]
        fPDR_total_dust_dict['Lib_name'] = [Mass_warm_dust_dict['Lib_name'], Mass_cold_dust_dict['Lib_name']]
        fPDR_total_dust_dict['Col_numb'] = [Mass_warm_dust_dict['Col_numb'], Mass_cold_dust_dict['Col_numb']]
        fPDR_total_dust_dict['Par_name'] = r'$\log_{10} \ f_{\mathrm{PDR}}$ (total)'
        fPDR_total_dust_dict['Par_file'] = 'fPDR_total'
        fPDR_total_dust_dict['Log_calc'] = True
        fPDR_total_dust_dict['range'] = [1e-6,1.0]
        #fPDR_total_dust_dict['Log_calc'] = False
        #fPDR_total_dust_dict['range'] = [0.0,1.0]
        # 
        # Note that this is not exactly 'gamma'
    # 
    if 'value' in Mass_warm_dust_dict and 'value' in Mass_cold_dust_dict:
        Umean_total_dust_dict = copy(Umin_warm_dust_dict)
        Umean_total_dust_dict['value'] = (1-fPDR_total_dust_dict['value']) * Umin_cold_dust_dict['value'] + fPDR_total_dust_dict['value'] * Umin_warm_dust_dict['value'] * (numpy.log(1e6/Umin_warm_dust_dict['value'])/(1-Umin_warm_dust_dict['value']/1e6))
        Umean_total_dust_dict['Lib_file'] = [Mass_warm_dust_dict['Lib_file'], Mass_cold_dust_dict['Lib_file']]
        Umean_total_dust_dict['Lib_name'] = [Mass_warm_dust_dict['Lib_name'], Mass_cold_dust_dict['Lib_name']]
        Umean_total_dust_dict['Col_numb'] = [Mass_warm_dust_dict['Col_numb'], Mass_cold_dust_dict['Col_numb']]
        Umean_total_dust_dict['Par_name'] = r'$\left<U\right>$ (total)'
        Umean_total_dust_dict['Par_file'] = 'Umean_total'
        Umean_total_dust_dict['range'] = [0.08,50.0]
        if InfoDict[Mass_warm_dust_dict['Lib_name']].find('DL07.2010.03.18') > 0: Umean_total_dust_dict['range'][1] = 90.0 #<20180319>#
        # 
        # Note: calc_Umean
        #   set Umin = 1.0
        #   set fPDR = 0.01
        #   calc ( lg(1e6/Umin) / (1-Umin/1e6) ) * Umin * fPDR + Umin * (1-fPDR)
    # 
    # analyze 
    print('Num_params', Num_params)
    print('Lib_params', Lib_params)
    Plot_engine = CrabPlot(figure_size=(14.0,10.0))
    Plot_engine.set_margin(panel=0, top=0.96, bottom=0.04, left=0.06, right=0.96)
    if 'value' in Stellar_mass_dict:
        analyze_chisq_distribution(Stellar_mass_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Stellar_age_dict:
        analyze_chisq_distribution(Stellar_age_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Stellar_EBV_dict:
        analyze_chisq_distribution(Stellar_EBV_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Lumin_AGN_dict:
        analyze_chisq_distribution(Lumin_AGN_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Umin_warm_dust_dict:
        analyze_chisq_distribution(Umin_warm_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Umin_cold_dust_dict:
        analyze_chisq_distribution(Umin_cold_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in LTIR_warm_dust_dict:
        analyze_chisq_distribution(LTIR_warm_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in LTIR_cold_dust_dict:
        analyze_chisq_distribution(LTIR_cold_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Mass_warm_dust_dict:
        analyze_chisq_distribution(Mass_warm_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Mass_cold_dust_dict:
        analyze_chisq_distribution(Mass_cold_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in LTIR_total_dust_dict:
        analyze_chisq_distribution(LTIR_total_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in LFIR_total_dust_dict:
        analyze_chisq_distribution(LFIR_total_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in LFIR122_total_dust_dict:
        analyze_chisq_distribution(LFIR122_total_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Mass_total_dust_dict:
        analyze_chisq_distribution(Mass_total_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in fPDR_total_dust_dict:
        analyze_chisq_distribution(fPDR_total_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Umean_total_dust_dict:
        analyze_chisq_distribution(Umean_total_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in dust_emissivity_beta_warm_dust_dict:
        analyze_chisq_distribution(dust_emissivity_beta_warm_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in dust_emissivity_beta_cold_dust_dict:
        analyze_chisq_distribution(dust_emissivity_beta_cold_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in dust_emissivity_beta_dict:
        analyze_chisq_distribution(dust_emissivity_beta_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in dust_temperature_warm_dust_dict:
        analyze_chisq_distribution(dust_temperature_warm_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in dust_temperature_cold_dust_dict:
        analyze_chisq_distribution(dust_temperature_cold_dust_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in dust_temperature_dict:
        analyze_chisq_distribution(dust_temperature_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    Plot_engine.set_xcharsize(panel=0, charsize=11, axislabelcharsize=16) # all panels
    Plot_engine.set_ycharsize(panel=0, charsize=11, axislabelcharsize=16) # all panels
    Plot_engine.savepdf(Output_dir+Output_name+'.chisq.pdf')
    #Plot_engine.show()
    Plot_engine.close()
    print('Output to "%s"!'%(Output_dir+Output_name+'.chisq.pdf'))
    # 
    






