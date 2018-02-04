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
import re
import json
from copy import copy






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

def analyze_chisq_distribution(param_dict, verbose = 1, Plot_engine = None):
    # Plot_engine must be the CrabPlot class
    if 'Lib_file' in param_dict and \
        'Lib_name' in param_dict and \
        'Lib_numb' in param_dict and \
        'Par_name' in param_dict and \
        'Col_numb' in param_dict and \
        'value' in param_dict and \
        'chisq' in param_dict :
        if verbose>=1:
            print('Analyzing the chi-square distribution for parameter "%s" in library %s from file "%s"'%(param_dict['Par_name'], param_dict['Lib_name'], param_dict['Lib_file']))
            #print(param_dict)
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
        # set xrange to the user-specified values
        xrange = None
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
        param_stats_2p = crab_bin_compute_param_chisq_histogram(chisq_array, param_array, min = param_min, max = param_max, delta_chisq = 2.3, log = param_log)
        if 'Par_file' in param_dict:
            # remove previous file
            if os.path.isfile('best-fit_param_'+param_dict['Par_file']+'.txt'):
                os.system('mv %s %s.backup'%('best-fit_param_'+param_dict['Par_file']+'.txt', 
                                            'best-fit_param_'+param_dict['Par_file']+'.txt'))
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
            asciitable.write(numpy.column_stack((param_median, param_best, param_sigma, param_L68, param_H68)), 
                                    'best-fit_param_'+param_dict['Par_file']+'.txt', Writer=asciitable.Ipac, 
                                            names=['param_median', 'param_best', 'param_sigma', 'param_L68', 'param_H68'], 
                                            formats={'param_median': '%20.10g', 'param_best': '%20.10g', 'param_sigma': '%20.10g', 'param_L68': '%20.10g', 'param_H68': '%20.10g'}, 
                                                delimiter='    ', overwrite = True)
        # 
        # crab_bin_compute_param_chisq_histogram for plotting
        param_stats = crab_bin_compute_param_chisq_histogram(chisq_array, param_array, min = param_min, max = param_max, delta_chisq = Delta_chisq_of_interest, log = param_log)
        # 
        param_bin_x = param_stats['hist_x']
        param_bin_y = param_stats['hist_y']
        param_bin_step = param_stats['bin_step']
        # 
        xrange = param_stats['xrange'] # xrange is the param range where chi-sq < min-chi-sq + 2.3
        yrange = param_stats['yrange']
        #xrange = [xrange[0]-(xrange[1]-xrange[0])*0.50, xrange[1]+(xrange[1]-xrange[0])*0.50] # extend the range for plotting.
        if param_stats['valid']:
            xrange = [ param_stats['L68'] - 10*param_stats['bin_step'], 
                       param_stats['H68'] + 10*param_stats['bin_step'] ]
        ##if xrange[0] < param_stats['min']: xrange[0] = param_stats['min']
        ##if xrange[1] > param_stats['max']: xrange[1] = param_stats['max']
        # invert y
        yrange = [1.0/yrange[1], 1.0/yrange[0]]
        yrange = numpy.log10(yrange)
        yrange = [yrange[0]-(yrange[1]-yrange[0])*0.50, yrange[1]+(yrange[1]-yrange[0])*0.05] # extend the range for plotting.
        yrange = numpy.power(10,yrange)
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
        #pprint(numpy.column_stack((param_bin_x, param_bin_y, 1/param_bin_y)))
        #print('------ xrange', xrange)
        #print('------ yrange', yrange, [1/yrange[1],1/yrange[0]])
        #print('------ param_stats.xrange', param_stats['xrange'])
        #print('------ param_stats.yrange', param_stats['yrange'], [1/param_stats['yrange'][1],1/param_stats['yrange'][0]])
        #print('------ param_stats.minimum_chisq', param_stats['minimum_chisq'])
        #print('------ param_stats.best_min_chisq', param_stats['best_min_chisq'])
        #print('------ param_stats.best', param_stats['best'])
        print('param_stats.min', param_stats['min'])
        print('param_stats.max', param_stats['max'])
        print('param_stats.xrange', param_stats['xrange'])
        print('param_stats.yrange', param_stats['yrange'], [1/param_stats['yrange'][1],1/param_stats['yrange'][0]])
        print('plotting xrange', xrange)
        print('plotting yrange', yrange)
        #--
        #--TODO--20180123-10h44m-- when param_log is True, param_min can be zero!
        #--
        # 
        # Initialize a plot
        if Plot_engine is None:
            Plot_engine = CrabPlot(figure_size=(9.0,5.0))
            Plot_engine.set_margin(panel=0, top=0.96, bottom=0.04)
        # 
        # Plot xy (left panel)
        Plot_engine.plot_xy(param_array, 1/numpy.array(chisq_array), overplot = False, 
                                xtitle = param_dict['Par_name'], ytitle = '$1/\chi^2$', useTex = True, 
                                size = 2.2, color='#1873cc', symbol = 'o')
        # 
        # Plot Cut_chi2 line
        Plot_engine.plot_line(param_stats['xrange'][0], 1/(chisq_min+Delta_chisq_of_interest), 
                                param_stats['xrange'][0], 1/(chisq_min+Delta_chisq_of_interest), 
                                overplot = True, linestyle = 'dashed')
        # 
        # Plot histogram (right panel)
        Plot_engine.plot_hist(param_bin_x, 1/numpy.array(param_bin_y), width = param_bin_step*1.5, align = 'edge', overplot = False, 
                                xtitle = param_dict['Par_name'], ytitle = '$1/\chi^2$', useTex = True, 
                                xrange = xrange, yrange = yrange, xlog = xlog, ylog = ylog)
        # 
        # Plot Cut_chi2 line
        Plot_engine.plot_line(xrange[0], 1/(chisq_min+Delta_chisq_of_interest), xrange[1], 1/(chisq_min+Delta_chisq_of_interest), overplot = True, linestyle = 'dashed')
        Plot_engine.plot_text(xrange[1], yrange[1]-0.02*(yrange[1]-yrange[0]), ' (zoomed) ', NormalizedCoordinate=False, overplot=True, horizontalalignment='right', verticalalignment='top')
        # 
        # Plot Cut_chi2 line (2p = 2.3)
        Plot_engine.plot_line(param_stats_2p['xrange'][0], 1/(chisq_min+2.3), 
                                param_stats_2p['xrange'][1], 1/(chisq_min+2.3), 
                                overplot = True, color='#1e90ff', linestyle = 'dotted')
                                # color: http://www.color-hex.com/color/1e90ff
        # 
    else:
        print('Error! analyze_chisq_distribution() got unaccepted inputs!')
        sys.exit()






def constrain_by_upper_limits(chisq_file, chisq_array, lib_dict):
    if os.path.isfile('flagged_chi2_solution.txt'):
        os.system('mv flagged_chi2_solution.txt flagged_chi2_solution.txt.backup')
    if os.path.isfile('flagged_chi2_solution_sorted_index.txt'):
        os.system('mv flagged_chi2_solution_sorted_index.txt flagged_chi2_solution_sorted_index.txt.backup')
    os.system('bash -c \"rm -rf obj_* 2>/dev/null\"')
    os.system('bash -c \"rm -rf dump_* 2>/dev/null\"')
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
            obs_wave_undetected = obs_wave[obs_undetection] / (1.0 + float(Input_info_dict['REDSHIFT']))
            obs_flux_undetected = obs_error[obs_undetection] * 5.0 # 5-sigma upper limit <TODO>
            # 
            # loop each input chi2 solution
            chisq_indices_sorted = numpy.argsort(chisq_array)
            i_constrain = 0
            while i_constrain < len(chisq_array):
                print('constrain_by_upper_limits: dump_LIB_SEDs_to_files(Input_info_dict, %d)'%(chisq_indices_sorted[i_constrain]))
                #Read_SED_LIB(DataFile, DataArray, Input_info_dict, chisq_indices_sorted[i_constrain])
                dump_LIB_SEDs_to_files(chisq_file = chisq_file, chisq_array = chisq_array, lib_dict = lib_dict, 
                                        dump_indices = chisq_indices_sorted[i_constrain], 
                                        output_numbers = 1, 
                                        output_prefix = 'dump')
                SED_data_table = asciitable.read('dump_1/SED_SUM') # always read the minimum chi2 solution
                SED_x = numpy.array(SED_data_table.field(SED_data_table.colnames[0]))
                SED_y = numpy.array(SED_data_table.field(SED_data_table.colnames[1]))
                #SED_flux_to_constrain = scipy.interpolate.spline(SED_x, SED_y, obs_wave_undetected, order='1') # order=3, kind='smoothest', conds=None
                SED_flux_to_constrain = scipy.interpolate.interp1d(SED_x, SED_y, kind='nearest')(obs_wave_undetected)
                # 
                where_constraint = (SED_flux_to_constrain > obs_flux_undetected)
                which_to_constrain = numpy.argwhere(where_constraint)
                if len(which_to_constrain) > 0:
                    # this chi2 solution is not allowed by the upper limit
                    chisq_array[chisq_indices_sorted[i_constrain]] = 1e+99
                    os.system('echo %d >> flagged_chi2_solution.txt'%(chisq_indices_sorted[i_constrain]))
                    os.system('echo %d >> flagged_chi2_solution_sorted_index.txt'%(i_constrain))
                os.system('rm -rf dump_1') # always read the minimum chi2 solution
                if len(which_to_constrain) <= 0:
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
    Plot_chi2_indices = numpy.sort(numpy.argsort(numpy.random.random(Cut_chi2_array_size))[0:Plot_chi2_max_number])[::-1] # non-repeated index array from 0 to Plot_chi2_max_number
    Plot_chi2_indices[0] = Cut_chi2_array_size-1 if Plot_chi2_indices[0] != Cut_chi2_array_size-1 else Cut_chi2_array_size-1 # make sure the first element is always 'Cut_chi2_array_size-1', i.e., the worst chi-square solution
    Plot_chi2_indices[-1] = 0 if Plot_chi2_indices[-1] != 0 else 0 # make sure the last element is always '0', i.e., the mininum chi-square solution
    # 
    Plot_chi2_index_dict = {}
    for i in range(len(Plot_chi2_indices)): 
        Plot_chi2_index_dict['%d'%(Plot_chi2_indices[i])] = Cut_chi2_array[i]
    # 
    return Plot_chi2_index_dict, Plot_chi2_indices





def dump_LIB_SEDs_to_files(chisq_file = '', chisq_array = [], lib_dict = {}, 
                            dump_indices = [], output_numbers = [], output_prefix = 'dump', 
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
    #output_prefix = 'dump'
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
        if not os.path.isdir('%s_%d'%(output_prefix, output_numbers[i])):
            os.mkdir('%s_%d'%(output_prefix, output_numbers[i]))
        # 
        # loop each SED LIB and dump LIB file (we do not overwrite existing files)
        for j in range(int(InfoDict['NLIB'])):
            # 
            if not os.path.isfile('%s_%d/SED_LIB%d'%(output_prefix, output_numbers[i], j+1)):
                #BashCommand = 'cd obj_%d/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_SED ../%s %d %s SED_LIB%d'%\
                #                    (i+1, \
                #                        InfoDict['LIB%d'%(j+1)], \
                #                            DataArray['i%d'%(j+1)][dump_indices[i]], \
                #                                DataArray['a%d'%(j+1)][dump_indices[i]], \
                #                                    j+1)
                #print(BashCommand)
                #os.system(BashCommand)
                # 
                # do python way 20180113
                BashCommand = 'michi2_read_lib_SEDs.py %s %d %s_%d > %s_%d/log.txt'%\
                                ( \
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
                #cd obj_8/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_SED ../lib.DL07.LoExCom.SED 140140 16.1932 SED_LIB4
                #/Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_SEDs.py fit_5.out 2451 c_2451
                #topcat -f ASCII c_2451/SED_LIB4 obj_8/SED_LIB4 &
                #checked that the two code give exactly the same result!
                #
                # how about the integrated IR luminosity?
                #cat obj_8/SED_LIB4.vLv_8_1000 # 8.8442327616e+03
                #cd obj_8
                #sm
                #load astroSfig.sm
                #data SED_LIB4 read {x 1 y 2}
                #calc_ltir x y # 3536.147921
                #calc 10**2.339198 * 16.1932 # 3536.150006 -- 2.339198 is the PAR3 in lib file, agreed with our manual integration! 
                # 
        # 
        # also dump chi2, line_number (in the chisq_file)
        if not os.path.isfile('%s_%d/chi2.txt'%(output_prefix, output_numbers[i])):
            BashCommand = 'echo "%s" > %s_%d/chi2.txt'%\
                            (chisq_array[dump_indices[i]], \
                                output_prefix, \
                                    output_numbers[i])
            print(BashCommand)
            os.system(BashCommand)
        # 
        if not os.path.isfile('%s_%d/line_number.txt'%(output_prefix, output_numbers[i])):
            BashCommand = 'echo "%s" > %s_%d/line_number.txt'%\
                            (dump_indices[i]+1, \
                                output_prefix, \
                                    output_numbers[i])
            print(BashCommand)
            os.system(BashCommand)
        # 
        # also dump redshift if possible
        if not numpy.isnan(redshift) and not os.path.isfile('%s_%d/redshift.txt'%(output_prefix, output_numbers[i])):
            BashCommand = 'echo "%s" > %s_%d/redshift.txt'%\
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
    SourceName = ''
    PlotYRange = []
    iarg = 1
    while iarg < len(sys.argv):
        TempCmd = sys.argv[iarg].replace('--','-').lower()
        if sys.argv[iarg]=='-only-plot-best-sed' or sys.argv[iarg]=='-only-best':
            SetOnlyPlotBestSED = True
            print('Setting only plot best-fit!')
        elif sys.argv[iarg]=='-source-name' or sys.argv[iarg]=='-source':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                SourceName = sys.argv[iarg]
                print('Setting SourceName = %s'%(SourceName))
        elif sys.argv[iarg]=='-yrange':
            if iarg+2 < len(sys.argv):
                iarg = iarg + 1
                PlotYRange.append(float(sys.argv[iarg]))
                iarg = iarg + 1
                PlotYRange.append(float(sys.argv[iarg]))
                print('Setting PlotYRange = %s'%(PlotYRange))
        else:
            DataFile = sys.argv[iarg]
        iarg = iarg + 1
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
    # 
    # 
    # 
    # Constrain DataArray by upper limits <20180202>
    if True == False:
        DataArray['chi2'] = constrain_by_upper_limits(DataFile, DataArray['chi2'], InfoDict)
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
    print('Selecting %d chi2 solutions with chi2 <= min(chi2)+%s'%(Cut_chi2_array_size, Delta_chisq_of_interest))
    # 
    if not SetOnlyPlotBestSED:
        if not os.path.isfile('Plot_chi2_index_dict.json') or not os.path.isfile('Plot_chi2_indices.json'):
            Plot_chi2_index_dict, Plot_chi2_indices = random_sorted_chi2_index_dict(Cut_chi2_array) # we plot 50 chi2 solution curves
            with open('Plot_chi2_index_dict.json', 'w') as fp:
                json.dump(Plot_chi2_index_dict, fp, sort_keys=True, indent=4)
                fp.close()
            with open('Plot_chi2_indices.json', 'w') as fp:
                json.dump(Plot_chi2_indices.tolist(), fp)
                fp.close()
        else:
            with open('Plot_chi2_index_dict.json', 'r') as fp:
                Plot_chi2_index_dict = json.load(fp)
                fp.close()
            with open('Plot_chi2_indices.json', 'r') as fp:
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
    # check if output figure already exists or not, see if we overwrite or not.
    Output_name, Output_extension = os.path.splitext(DataFile)
    # 
    # Get Redshift
    Redshift = float(InfoDict['REDSHIFT'])
    # 
    # Get SED (dump to subdirectories and files)
    #Read_SED_LIB(DataFile, DataArray, InfoDict, All_chi2_indices_sorted, Cut_chi2_array_size, Plot_chi2_index_dict)
    dump_LIB_SEDs_to_files(chisq_file = DataFile, chisq_array = All_chi2_array, 
                            lib_dict = InfoDict, 
                            dump_indices = All_chi2_indices_sorted[Plot_chi2_indices], 
                            output_numbers = Plot_chi2_indices+1, 
                            output_prefix = 'obj', 
                            redshift = Redshift)
    # 
    # Wait for a long time
    # 
    # set color styles
    Color_list = ['cyan', 'gold', 'red', 'blue', 'purple']
    Plot_engine = CrabPlot(figure_size=(8.0,5.0))
    Plot_engine.set_margin(top=0.92, bottom=0.16, left=0.12, right=0.96)
    Count_label_chi2 = 0 # to count the chi-square label printed on the figure, make sure there are not too many labels.
    Count_plot_chi2 = 0
    # 
    # Then plot SEDs
    for i in Plot_chi2_indices[::-1]:
        # 
        # alpha by chi2
        print('Plotting chi2=%s obj_%d'%(Cut_chi2_array[i], i+1))
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
            if j == 0: xclip = [(50,numpy.inf)]
            elif j == 4: xclip = [(-numpy.inf,2e3)]
            Plot_engine.plot_data_file('obj_%d/SED_LIB%d'%(i+1,j+1), xlog=1, ylog=1, xclip=xclip, current=1, \
                                dataname='obj_%d_SED_LIB%d'%(i+1,j+1), 
                                redshift = Redshift, 
                                linestyle='dashed', linewidth=Plot_chi2_linewidth, color=Color_list[j], alpha=Plot_chi2_alpha)
        # 
        Min_chi2_for_plot = numpy.power(10, Min_chi2_log-(Max_chi2_log-Min_chi2_log)*0.05)
        Max_chi2_for_plot = numpy.power(10, Max_chi2_log+(Max_chi2_log-Min_chi2_log)*0.85)
        Color_chi2 = Plot_engine.get_color_by_value([Min_chi2_for_plot, Max_chi2_for_plot], 
                                                    input_value=Cut_chi2_array[i], 
                                                    log=1, 
                                                    cmap=matplotlib.cm.get_cmap('gray'))
        #print('Color_chi2: ', Color_chi2)
        # 
        # 
        # plot total SED
        Plot_engine.plot_data_file('obj_%d/SED_SUM'%(i+1), xlog=1, ylog=1, current=1, \
                            dataname='obj_%d_SED_SUM'%(i+1), 
                            redshift = Redshift, 
                            linestyle='solid', linewidth=Plot_SED_linewidth, color=Color_chi2, alpha=1.0, zorder=8) # alpha=Plot_chi2_alpha
        # 
        # 
        # count++
        Count_plot_chi2 = Count_plot_chi2 + 1
        # 
        # 
        # show chi2 text on the figure
        if not SetOnlyPlotBestSED:
            if i == Cut_chi2_array_size-1:
                Plot_engine.xyouts(0.05, 0.95, '$\chi^2:$', NormalizedCoordinate=True, useTex=True)
            if i == 0:
                Plot_engine.xyouts(0.09, 0.95-0.03*(Count_label_chi2), '......', NormalizedCoordinate=True, color=Color_chi2)
                Count_label_chi2 = Count_label_chi2 + 1
            if Count_plot_chi2 % int((Cut_chi2_array_size/7)+1) == 0 or i == 0 or i == Cut_chi2_array_size-1:
                #print('Plotting label at', 0.09, 0.95-0.03*(Cut_chi2_array_size-1-i), 'chi2 = %.1f'%(Cut_chi2_array[i]))
                Plot_engine.xyouts(0.09, 0.95-0.03*(Count_label_chi2), '%.1f'%(Cut_chi2_array[i]), NormalizedCoordinate=True, useTex=True, color=Color_chi2)
                Count_label_chi2 = Count_label_chi2 + 1
        # 
        # 
        # show redshift (z) and source name on the figure
        if not SetOnlyPlotBestSED:
            if i == 0:
                Plot_engine.xyouts(0.15, 0.95, '$z=%s$'%(Redshift), NormalizedCoordinate=True, useTex=True)
            if i == 0 and SourceName != '':
                Plot_engine.xyouts(0.97, 0.90, SourceName, NormalizedCoordinate=True, fontsize=16, horizontalalignment='right')
        else:
            if i == 0 and SourceName != '':
                Plot_engine.xyouts(0.05, 0.90, SourceName, NormalizedCoordinate=True, fontsize=15)
                Plot_engine.xyouts(0.20, 0.90, '$z=%s$'%(Redshift), NormalizedCoordinate=True, useTex=True, fontsize=15)
        #break
    # 
    # 
    # Then plot OBS data points
    DataTable_obs = asciitable.read(InfoDict['OBS'])
    #print(type(DataTable_obs))
    try:
        Wavelength_obs = DataTable_obs[DataTable_obs.colnames[0]].data
        Flux_obs = DataTable_obs[DataTable_obs.colnames[1]].data
        FluxErr_obs = DataTable_obs[DataTable_obs.colnames[2]].data
        Detection_mask = (Flux_obs>=2.0*FluxErr_obs)
        UpperLimits_mask = (Flux_obs<2.0*FluxErr_obs)
        #print(Wavelength_obs)
        Plot_engine.plot_xy(Wavelength_obs[Detection_mask], Flux_obs[Detection_mask], yerr=FluxErr_obs[Detection_mask], dataname='obs', overplot=True, symbol='open square', symsize=3, thick=1.5, capsize=4, zorder=10)
        Plot_engine.plot_xy(Wavelength_obs[UpperLimits_mask], 3.0*FluxErr_obs[UpperLimits_mask], dataname='upper limits', overplot=True, symbol='upper limits', symsize=3, thick=1.25, alpha=0.5, zorder=9)
    except Exception as err:
        print(err)
    # 
    # 
    Plot_engine.set_xrange([0.1,1e6])
    Plot_engine.set_yrange([1e-6,1e4])
    if len(PlotYRange) == 2:
        Plot_engine.set_yrange(PlotYRange)
    Plot_engine.set_xtitle('Observing-frame wavelength [um]')
    Plot_engine.set_ytitle('Flux density [mJy]')
    Plot_engine.savepdf(Output_name+'.pdf')
    #Plot_engine.show()
    Plot_engine.close()
    print('Output to "%s"!'%(Output_name+'.pdf'))
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
    LTIR_warm_dust_dict = {}
    LTIR_cold_dust_dict = {}
    Umin_warm_dust_dict = {}
    Umin_cold_dust_dict = {}
    Mass_warm_dust_dict = {}
    Mass_cold_dust_dict = {}
    Lumin_AGN_dict = {}
    LTIR_total_dust_dict = {}
    Mass_total_dust_dict = {}
    fPDR_total_dust_dict = {}
    Umean_total_dust_dict = {}
    # 
    # define constants
    pi = numpy.pi
    dL = 1.0/(4*pi) # if Redshift is not given, then we do not apply 4*pi*dL**2 to the quantities below. 
    Redshift = float(InfoDict['REDSHIFT'])
    print('z = %s'%(Redshift))
    if Redshift > 0.0:
        lumdist_command = '%s/lumdist'%(os.path.dirname(os.path.realpath(__file__)))
        if not os.path.isfile(lumdist_command):
            print('Error! "lumdist" command was not found at %s! Maybe you have not fully downloaded the code from "https://github.com/1054/Crab.Toolkit.michi2"?')
            print('We need that to compute the luminosity distance!')
            sys.exit()
        dL_command = '%s -simple %s'%(lumdist_command, Redshift)
        dL_str = os.popen(dL_command).read()
        #print(dL_str)
        if dL_str == '':
            print('Error! Failed to run "lumdist" in the Terminal! Sorry! TODO: You can set "dL = ..." in "%s".'%(InfoFile))
            print('We need that to compute the luminosity distance!')
            sys.exit()
        dL = float(dL_str)
        print('dL = %s [Mpc]'%(dL))
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
            if InfoDict[Lib_name].find('BC03.Padova1994') >= 0:
                if 'Mass' == Lib_dict[Key_TPAR]:
                    # Stellar mass SED library Y column unit is solar luminosity per Hertz, 
                    # flux unit is mJy, 
                    # 
                    Stellar_mass_dict['Lib_file'] = InfoDict[Lib_name]
                    Stellar_mass_dict['Lib_name'] = Lib_name
                    Stellar_mass_dict['Lib_numb'] = j+1
                    Stellar_mass_dict['Par_name'] = '$\log \ M_{*}$ [$\mathrm{M}_{\odot}$]' # Lib_dict[Key_TPAR]
                    Stellar_mass_dict['Par_file'] = 'Mstar'
                    Stellar_mass_dict['Col_numb'] = Col_number
                    Stellar_mass_dict['Log_calc'] = True
                    Stellar_mass_dict['range'] = numpy.power(10,[8.0,13.5])
                    Stellar_mass_dict['value'] = DataArray['a%d'%(j+1)] / (3.839e33*1e26/(4*pi*dL**2*9.52140e48)) * DataTable.getColumn(Col_number) / (1+Redshift)
                    Stellar_mass_dict['chisq'] = DataArray['chi2']
            elif InfoDict[Lib_name].find('DL07.HiExCom') >= 0:
                if 'lgLTIR' == Lib_dict[Key_TPAR]:
                    LTIR_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LTIR_warm_dust_dict['Lib_name'] = Lib_name
                    LTIR_warm_dust_dict['Lib_numb'] = j+1
                    LTIR_warm_dust_dict['Par_name'] = '$\log \ L_{\mathrm{IR}}$ (warm) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LTIR_warm_dust_dict['Par_file'] = 'LIR_warm'
                    LTIR_warm_dust_dict['Col_numb'] = Col_number
                    LTIR_warm_dust_dict['Log_calc'] = True
                    LTIR_warm_dust_dict['range'] = numpy.power(10,[9.0,14.5])
                    LTIR_warm_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn(Col_number)) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LTIR_warm_dust_dict['chisq'] = DataArray['chi2']
                elif 'Umin' == Lib_dict[Key_TPAR]:
                    Umin_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Umin_warm_dust_dict['Lib_name'] = Lib_name
                    Umin_warm_dust_dict['Lib_numb'] = j+1
                    Umin_warm_dust_dict['Par_name'] = '$U_{\mathrm{min}}$ (warm)' # Lib_dict[Key_TPAR]
                    Umin_warm_dust_dict['Par_file'] = 'Umin_warm'
                    Umin_warm_dust_dict['Col_numb'] = Col_number
                    Umin_warm_dust_dict['Log_plot'] = True # 'Log_plot', plot X axis in log scale
                    Umin_warm_dust_dict['range'] = [0.08,30.0]
                    Umin_warm_dust_dict['value'] = DataTable.getColumn(Col_number)
                    Umin_warm_dust_dict['chisq'] = DataArray['chi2']
                    # 
                    Mass_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Mass_warm_dust_dict['Lib_name'] = Lib_name
                    Mass_warm_dust_dict['Lib_numb'] = j+1
                    Mass_warm_dust_dict['Par_name'] = '$\log \ M_{\mathrm{dust}}$ (warm) [$\mathrm{M}_{\odot}$]'
                    Mass_warm_dust_dict['Par_file'] = 'Mdust_warm'
                    Mass_warm_dust_dict['Col_numb'] = 2+2*(j+1)
                    Mass_warm_dust_dict['Log_calc'] = True
                    Mass_warm_dust_dict['range'] = numpy.power(10,[7.0,12.0])
                    Mass_warm_dust_dict['value'] = DataArray['a%d'%(j+1)] * dL**2 / (1+Redshift) # Mdust #NOTE# no need to multiply a '4*pi'!
                    Mass_warm_dust_dict['chisq'] = DataArray['chi2']
                    # 
            elif InfoDict[Lib_name].find('DL07.LoExCom') >= 0:
                if 'lgLTIR' == Lib_dict[Key_TPAR]:
                    LTIR_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LTIR_cold_dust_dict['Lib_name'] = Lib_name
                    LTIR_cold_dust_dict['Lib_numb'] = j+1
                    LTIR_cold_dust_dict['Par_name'] = '$\log \ L_{\mathrm{IR}}$ (cold) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LTIR_cold_dust_dict['Par_file'] = 'LIR_cold'
                    LTIR_cold_dust_dict['Col_numb'] = Col_number
                    LTIR_cold_dust_dict['Log_calc'] = True
                    LTIR_cold_dust_dict['range'] = numpy.power(10,[9.0,14.5])
                    LTIR_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn(Col_number)) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LTIR_cold_dust_dict['chisq'] = DataArray['chi2']
                elif 'Umin' == Lib_dict[Key_TPAR]:
                    Umin_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Umin_cold_dust_dict['Lib_name'] = Lib_name
                    Umin_cold_dust_dict['Lib_numb'] = j+1
                    Umin_cold_dust_dict['Par_name'] = '$U_{\mathrm{min}}$ (cold)' # Lib_dict[Key_TPAR]
                    Umin_cold_dust_dict['Par_file'] = 'Umin_cold'
                    Umin_cold_dust_dict['Col_numb'] = Col_number
                    Umin_cold_dust_dict['Log_plot'] = True # 'Log_plot', plot X axis in log scale
                    Umin_cold_dust_dict['range'] = [0.08,30.0]
                    Umin_cold_dust_dict['value'] = DataTable.getColumn(Col_number)
                    Umin_cold_dust_dict['chisq'] = DataArray['chi2']
                    # 
                    Mass_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Mass_cold_dust_dict['Lib_name'] = Lib_name
                    Mass_cold_dust_dict['Lib_numb'] = j+1
                    Mass_cold_dust_dict['Par_name'] = '$\log \ M_{\mathrm{dust}}$ (cold) [$\mathrm{M}_{\odot}$]'
                    Mass_cold_dust_dict['Par_file'] = 'Mdust_cold'
                    Mass_cold_dust_dict['Col_numb'] = 2+2*(j+1)
                    Mass_cold_dust_dict['Log_calc'] = True
                    Mass_cold_dust_dict['range'] = numpy.power(10,[7.0,12.0])
                    Mass_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * dL**2 / (1+Redshift) # Mdust # Mdust #NOTE# no need to multiply a '4*pi'!
                    Mass_cold_dust_dict['chisq'] = DataArray['chi2']
            elif InfoDict[Lib_name].find('MullaneyAGN') >= 0 or \
                InfoDict[Lib_name].find('SiebenmorgenAGN') >= 0:
                # Mullaney AGN
                #   by integrating the AGN template (AGN_TYPE=2) from 1um to 1000um, we get an integration of 5133.913101
                if 'AGN_TYPE' == Lib_dict[Key_TPAR].upper():
                    Lumin_AGN_dict['Lib_file'] = InfoDict[Lib_name]
                    Lumin_AGN_dict['Lib_name'] = Lib_name
                    Lumin_AGN_dict['Lib_numb'] = j+1
                    Lumin_AGN_dict['Par_name'] = '$\log \ L_{\mathrm{AGN}}$ [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    Lumin_AGN_dict['Par_file'] = 'LAGN'
                    Lumin_AGN_dict['Col_numb'] = Col_number
                    Lumin_AGN_dict['Log_calc'] = True
                    Lumin_AGN_dict['range'] = numpy.power(10,[0.0,14.5])
                    Lumin_AGN_dict['value'] = DataArray['a%d'%(j+1)] * 5133.913101 * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    Lumin_AGN_dict['chisq'] = DataArray['chi2']
                    # 
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
        LTIR_total_dust_dict['Par_name'] = '$\log \ L_{\mathrm{IR}}$ (total) [$\mathrm{L}_{\odot}$]'
        LTIR_total_dust_dict['Par_file'] = 'LIR_total'
    # 
    # Total Mdust
    if 'value' in Mass_warm_dust_dict and 'value' in Mass_cold_dust_dict:
        Mass_total_dust_dict = copy(Mass_warm_dust_dict)
        Mass_total_dust_dict['value'] = Mass_warm_dust_dict['value'] + Mass_cold_dust_dict['value']
        Mass_total_dust_dict['Lib_file'] = [Mass_warm_dust_dict['Lib_file'], Mass_cold_dust_dict['Lib_file']]
        Mass_total_dust_dict['Lib_name'] = [Mass_warm_dust_dict['Lib_name'], Mass_cold_dust_dict['Lib_name']]
        Mass_total_dust_dict['Col_numb'] = [Mass_warm_dust_dict['Col_numb'], Mass_cold_dust_dict['Col_numb']]
        Mass_total_dust_dict['Par_name'] = '$\log \ M_{\mathrm{dust}}$ (total) [$\mathrm{M}_{\odot}$]'
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
        fPDR_total_dust_dict['Par_name'] = '$\log \ \delta_{\mathrm{PDR}}$ (total)'
        fPDR_total_dust_dict['Par_file'] = 'fPDR_total'
        fPDR_total_dust_dict['Log_calc'] = True
        fPDR_total_dust_dict['range'] = [1e-4,1.0]
    # 
    if 'value' in Mass_warm_dust_dict and 'value' in Mass_cold_dust_dict:
        Umean_total_dust_dict = copy(Umin_warm_dust_dict)
        Umean_total_dust_dict['value'] = (1-fPDR_total_dust_dict['value']) * Umin_cold_dust_dict['value'] + fPDR_total_dust_dict['value'] * Umin_warm_dust_dict['value'] * (numpy.log(1e6/Umin_warm_dust_dict['value'])/(1-Umin_warm_dust_dict['value']/1e6))
        Umean_total_dust_dict['Lib_file'] = [Mass_warm_dust_dict['Lib_file'], Mass_cold_dust_dict['Lib_file']]
        Umean_total_dust_dict['Lib_name'] = [Mass_warm_dust_dict['Lib_name'], Mass_cold_dust_dict['Lib_name']]
        Umean_total_dust_dict['Col_numb'] = [Mass_warm_dust_dict['Col_numb'], Mass_cold_dust_dict['Col_numb']]
        Umean_total_dust_dict['Par_name'] = '$\\left<U\\right>$ (total)'
        Umean_total_dust_dict['Par_file'] = 'Umean_total'
        Umean_total_dust_dict['range'] = [0.08,50.0]
    # 
    # analyze 
    print('Num_params', Num_params)
    print('Lib_params', Lib_params)
    Plot_engine = CrabPlot(figure_size=(14.0,10.0))
    Plot_engine.set_margin(panel=0, top=0.96, bottom=0.04, left=0.06, right=0.96)
    if 'value' in Stellar_mass_dict:
        analyze_chisq_distribution(Stellar_mass_dict, Plot_engine = Plot_engine)
    if 'value' in Lumin_AGN_dict:
        analyze_chisq_distribution(Lumin_AGN_dict, Plot_engine = Plot_engine)
    if 'value' in Umin_warm_dust_dict:
        analyze_chisq_distribution(Umin_warm_dust_dict, Plot_engine = Plot_engine)
    if 'value' in Umin_cold_dust_dict:
        analyze_chisq_distribution(Umin_cold_dust_dict, Plot_engine = Plot_engine)
    if 'value' in LTIR_warm_dust_dict:
        print(LTIR_warm_dust_dict)
        analyze_chisq_distribution(LTIR_warm_dust_dict, Plot_engine = Plot_engine)
    if 'value' in LTIR_cold_dust_dict:
        analyze_chisq_distribution(LTIR_cold_dust_dict, Plot_engine = Plot_engine)
    if 'value' in Mass_warm_dust_dict:
        analyze_chisq_distribution(Mass_warm_dust_dict, Plot_engine = Plot_engine)
    if 'value' in Mass_cold_dust_dict:
        analyze_chisq_distribution(Mass_cold_dust_dict, Plot_engine = Plot_engine)
    if 'value' in LTIR_total_dust_dict:
        analyze_chisq_distribution(LTIR_total_dust_dict, Plot_engine = Plot_engine)
    if 'value' in Mass_total_dust_dict:
        analyze_chisq_distribution(Mass_total_dust_dict, Plot_engine = Plot_engine)
    if 'value' in fPDR_total_dust_dict:
        analyze_chisq_distribution(fPDR_total_dust_dict, Plot_engine = Plot_engine)
    if 'value' in Umean_total_dust_dict:
        analyze_chisq_distribution(Umean_total_dust_dict, Plot_engine = Plot_engine)
    Plot_engine.savepdf(Output_name+'.chisq.pdf')
    #Plot_engine.show()
    Plot_engine.close()
    print('Output to "%s"!'%(Output_name+'.chisq.pdf'))
    # 
    






