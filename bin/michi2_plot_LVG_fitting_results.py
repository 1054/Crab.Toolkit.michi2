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
from astropy.table import Table
import re
import json
from copy import copy, deepcopy

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

import matplotlib.pyplot as plt





#########################################
#               Constants               #
#########################################

#Delta_chisq_of_interest = 5.89 # for 5 interested parameters, 68.3% confidence (1-sigma).
Delta_chisq_of_interest = 3.505883 # for 3 interested parameters, 68.3% confidence (1-sigma).
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
        param_array = deepcopy(param_dict['value'])
        chisq_array = deepcopy(param_dict['chisq'])
        chisq_min = numpy.nanmin(chisq_array)
        # 
        # copy
        #param_array_nonan = deepcopy(param_array)
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
        # 20210906
        if numpy.all(numpy.isclose(numpy.diff(param_array), 0.0, atol=1e-10*numpy.nanmin(param_array))):
            return
        # 
        # crab_bin_compute_param_chisq_histogram
        param_stats = crab_bin_compute_param_chisq_histogram(chisq_array, param_array, \
                            min = param_min, max = param_max, 
                            delta_chisq = Delta_chisq_of_interest, log = param_log, verbose = verbose)
        # 
        # crab_bin_compute_param_chisq_histogram for delta_chisq = 2.3 (2p)
        param_stats_2p = crab_bin_compute_param_chisq_histogram(chisq_array, param_array, 
                            min = param_min, max = param_max, 
                            delta_chisq = 2.3, log = param_log, verbose = verbose)
        # 
        # write to disk
        if 'Par_file' in param_dict:
            # remove previous file
            if os.path.isfile(Output_dir+'best-fit_param_'+param_dict['Par_file']+'.txt'):
                os.system('mv %s %s.backup'%(Output_dir+'best-fit_param_'+param_dict['Par_file']+'.txt', 
                                             Output_dir+'best-fit_param_'+param_dict['Par_file']+'.txt'))
            if os.path.isfile(Output_dir+'chi-square_table_'+param_dict['Par_file']+'.txt'):
                os.system('mv %s %s.backup'%(Output_dir+'chi-square_table_'+param_dict['Par_file']+'.txt', 
                                             Output_dir+'chi-square_table_'+param_dict['Par_file']+'.txt'))
            if param_stats['valid'] is True:
                param_median = param_stats['median']
                param_best = param_stats['best']
                param_sigma = param_stats['sigma']
                param_L68 = param_stats['L68']
                param_H68 = param_stats['H68']
            else:
                param_median = 0.0
                param_best = 0.0
                param_sigma = 0.0
                param_L68 = 0.0
                param_H68 = 0.0
            asciitable.write(numpy.column_stack((param_median, param_best, param_sigma, param_L68, param_H68)), 
                                    Output_dir+'best-fit_param_'+param_dict['Par_file']+'.txt', Writer=asciitable.Ipac, 
                                            names=['param_median', 'param_best', 'param_sigma', 'param_L68', 'param_H68'], 
                                            formats={'param_median': '%20.10g', 'param_best': '%20.10g', 'param_sigma': '%20.10g', 'param_L68': '%20.10g', 'param_H68': '%20.10g'}, 
                                            delimiter='    ', overwrite = True)
            asciitable.write(numpy.column_stack((chisq_array, param_array)), 
                                    Output_dir+'chi-square_table_'+param_dict['Par_file']+'.txt', 
                                            Writer=asciitable.FixedWidth, 
                                            bookend=True, delimiter=' ',
                                            names=['chisq_array', 'param_array'], 
                                            overwrite = True)
            with open(Output_dir+'chi-square_table_'+param_dict['Par_file']+'.txt', 'r+') as fp:
                fp.seek(0)
                fp.write('#')
        # 
        xlog = None
        #if 'Log_plot' in param_dict:
        #    if param_dict['Log_plot'] == True:
        #        xlog = 1 # not working for matplotlib bar plot (i.e., CrabPlot plot_hist)!
        # 
        ylog = None
        # 
        # verbose
        if verbose >= 2:
            print('param_stats.min', param_stats['min'])
            print('param_stats.max', param_stats['max'])
            print('param_stats.xrange', param_stats['xrange'])
            print('param_stats.yrange', param_stats['yrange'], '(1/chisq ',1/param_stats['yrange'][1],' ',1/param_stats['yrange'][0],')')
        # 
        # calc log for plotting
        if param_log is True:
            with numpy.errstate(invalid='ignore'):
                param_array_mask = (param_array>0)
                param_array_mask2 = (param_array<=0)
            param_array[param_array_mask] = numpy.log10(param_array[param_array_mask])
            param_array[param_array_mask2] = numpy.nan
        # 
        # range for plotting
        plot_xrange = [0.0, 0.0]
        plot_yrange = [0.0, 0.0]
        plot_xrange[0] = max(param_stats['min'], param_stats['L68'] - 15.0*param_stats['sigma'])
        plot_xrange[1] = min(param_stats['max'], param_stats['H68'] + 15.0*param_stats['sigma'])
        plot_yrange[0] = 0.0
        plot_yrange[1] = 1.0/param_stats['min_chisq'] * 1.1 # y-axis is 1/chisq rather than chisq
        # 
        # 
        # Initialize a plot
        if Plot_engine is None:
            Plot_engine = CrabPlot(figure_size=(9.0*0.8,5.0*0.8)) # figure_size=(9.0,5.0) # 20221103
            Plot_engine.set_margin(panel=0, top=0.96, bottom=0.08)
        # 
        # Plot xy (left panel)
        Plot_engine.plot_hist(param_stats['hist_x'], 
                              1./param_stats['hist_y'], 
                              width=param_stats['hist_dx']*0.9, 
                              overplot = False, 
                              #alpha=0.80, color='C0', # 20221103
                              alpha=0.5, color='#1e90ff', # 20221103
                              xtitle = param_dict['Par_name'], ytitle = r'$1/\chi^2$', useTex = True, 
                              #xrange = plot_xrange, yrange = plot_yrange, 
                              xlog = None, ylog = None)
        # 
        # print("param_stats['smooth_x']", param_stats['smooth_x']) # debugging 20221101
        Plot_engine.plot_hist(param_stats['smooth_x'], 
                              1./param_stats['smooth_y'], 
                              #alpha=0.25, color='C0', # 20221103
                              alpha=0.5, color='#1e90ff', # 20221103
                              overplot = True)
        # 
        # Plot Cut_chi2 line
        if param_stats['valid']:
            Plot_engine.fill_between(param_stats['xrange'], 
                                     #numpy.array([0.0,0.0])+1./(param_stats['yrange'][1]), 
                                     numpy.array([0.0,0.0]), 
                                     numpy.array([0.0,0.0])+1./(param_stats['yrange'][0]), 
                                     overplot = True, color='gold', linestyle = 'solid', capstyle = 'butt', alpha = 0.5, zorder=9)
                                     # color: http://www.color-hex.com/color/1e90ff
                                     # https://stackoverflow.com/questions/35099130/change-spacing-of-dashes-in-dashed-line-in-matplotlib
        # 
    else:
        print('Error! analyze_chisq_distribution() got unaccepted inputs!')
        sys.exit()





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





def dump_LIB_LVGs_to_files(chisq_file = '', chisq_array = [], lib_dict = {}, 
                            dump_indices = [], output_numbers = [], output_prefix = 'dump_', 
                            redshift = numpy.nan):
    # 
    #def Read_LVG_LIB(DataFile, DataArray, InfoDict, chisq_indices_sorted, Cut_chi2_array_size = 1, Plot_chi2_index_dict = [])
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
        # loop each LVG LIB and dump LIB file (we do not overwrite existing files)
        for j in range(int(InfoDict['NLIB'])):
            # 
            if not os.path.isfile('%s%d/LVG_LIB%d'%(output_prefix, output_numbers[i], j+1)):
                #BashCommand = 'cd dump/%d/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_LVG ../%s %d %s LVG_LIB%d'%\
                #                    (i+1, \
                #                        InfoDict['LIB%d'%(j+1)], \
                #                            DataArray['i%d'%(j+1)][dump_indices[i]], \
                #                                DataArray['a%d'%(j+1)][dump_indices[i]], \
                #                                    j+1)
                #print(BashCommand)
                #os.system(BashCommand)
                # 
                # do python way 20180113
                BashCommand = '%s/michi2_read_lib_LVGs.py %s %d %s%d > %s%d/log.txt'%\
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
                #cd dump/8/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_LVG ../lib.DL07.LoExCom.LVG 140140 16.1932 LVG_LIB4
                #/Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_LVGs.py fit_5.out 2451 c_2451
                #topcat -f ASCII c_2451/LVG_LIB4 dump/8/LVG_LIB4 &
                #checked that the two code give exactly the same result!
                #
                # how about the integrated IR luminosity?
                #cat dump/8/LVG_LIB4.vLv_8_1000 # 8.8442327616e+03
                #cd dump/8
                #sm
                #load astroSfig.sm
                #data LVG_LIB4 read {x 1 y 2}
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
        #return LVG_x, LVG_y





def convert_X_species_to_ID_Name_J_upper_J_lower(input_X_species):
    temp_value = input_X_species
    ID_species = int(temp_value/1e6); temp_value -= ID_species*1e6
    J_upper = int(temp_value/1e3); temp_value -= J_upper*1e3
    J_lower = int(temp_value)
    Name_dict = {}
    Name_dict['101'] = 'CO'
    Name_dict['102'] = 'C_atom' # check '/Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/data/make_lib_LVG/02_Make_Lib_LVG/makeliblvg_z_with_CO_and_C_atom.sm'
    if '%d'%(ID_species) in Name_dict:
        Name_species = Name_dict['%d'%(ID_species)]
    else:
        Name_species = 'Unknown'
    return ID_species, Name_species, J_upper, J_lower





















######################################
#               MAIN PROGRAM             #
##########################################

if len(sys.argv) <= 1:
    
    print('Usage: michi2_plot_LVG_fitting_results.py fit_5.out')
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
    UserInputDistance = 0.0 # in Mpc
    UserInputArea = 0.0 # in kpc^2
    UserInputNormalization = 0.0 # CO10 flux scale normalization, will multiply this to flux and mass
    UserInputChisqScale = 1.0
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
                UserInputText = UserInputText.append(sys.argv[iarg])
                print('Setting UserInputText += %s'%(sys.argv[iarg]))
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
        elif TempCmd=='-distance':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                UserInputDistance = float(sys.argv[iarg])
                print('Setting UserInputDistance = %s'%(sys.argv[iarg]))
        elif TempCmd=='-area':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                UserInputArea = float(sys.argv[iarg])
                print('Setting UserInputArea = %s'%(sys.argv[iarg]))
        elif TempCmd=='-norm' or TempCmd=='-normalization':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                UserInputNormalization = float(sys.argv[iarg])
                print('Setting UserInputNormalization = %s'%(sys.argv[iarg]))
        elif TempCmd=='-chisq-scale' or TempCmd=='-chisqscale':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                UserInputChisqScale = float(sys.argv[iarg])
                print('Setting UserInputChisqScale = %s'%(sys.argv[iarg]))
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
    else:
        if Output_name.find(os.sep) >= 0:
            Output_dir = os.path.dirname(Output_name) + os.sep
            Output_name = os.path.basename(Output_name)
            if not os.path.isdir(Output_dir):
                os.makedirs(Output_dir)
    #print('Output_dir', Output_dir)
    #sys.exit()
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
    Redshift = 0.0 # float(InfoDict['REDSHIFT'])
    # 
    # Read OBS data file
    if UserInputFluxFile == '':
        DataTable_obs = asciitable.read(InfoDict['OBS'])
    else:
        DataTable_obs = asciitable.read(UserInputFluxFile)
    try:
        Line_transition_obs = DataTable_obs[DataTable_obs.colnames[0]].data
        Flux_obs = DataTable_obs[DataTable_obs.colnames[1]].data
        FluxErr_obs = DataTable_obs[DataTable_obs.colnames[2]].data
        Detection_mask = (Flux_obs>=2.0*FluxErr_obs)
        UpperLimits_mask = (Flux_obs<2.0*FluxErr_obs)
        #print(Line_transition_obs)
    except Exception as err:
        print(err)
        raise
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
    # Rescale chisq if needed <20221101>
    if UserInputChisqScale > 0.0:
        DataArray['chi2'] *= UserInputChisqScale
        FluxErr_obs /= numpy.sqrt(UserInputChisqScale)
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
    print('')
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
    # Get LVG (dump to subdirectories and files)
    #Read_LVG_LIB(DataFile, DataArray, InfoDict, All_chi2_indices_sorted, Cut_chi2_array_size, Plot_chi2_index_dict)
    print('')
    print('Dumping library LVGs')
    dump_LIB_LVGs_to_files(chisq_file = DataFile, chisq_array = All_chi2_array, 
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
    Plot_engine = CrabPlot(figure_size=(8.0*0.8,5.0*0.8)) # figure_size=(8.0,5.0) # 20221103
    Plot_engine.set_margin(top=0.92, bottom=0.22, left=0.12, right=0.96)
    Count_label_chi2 = 0 # to count the chi-square label printed on the figure, make sure there are not too many labels.
    Count_label_rchi2 = 0 # to count the reduced-chi-square label printed on the figure, make sure there are not too many labels.
    Count_plot_chi2 = 0
    #Count_plot_rchi2 = 0 # using Count_plot_chi2 instead
    # 
    # 
    # Then plot LVGs
    print('')
    for i in Plot_chi2_indices:
        # 
        # alpha by chi2
        if i >= len(Cut_chi2_array): continue
        print('Plotting chi2=%s, reduced-chi2=%s, dump id %d'%(Cut_chi2_array[i], Cut_chi2_array[i]/numpy.count_nonzero(Detection_mask), i+1))
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
        # plot each single LVG component
        Sum_X = None
        Sum_Y = None
        for j in range(int(InfoDict['NLIB'])):
            
            xclip = None
            
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
            
            tb = Table.read(Output_dir+'dump'+os.sep+'%d/LVG_LIB%d'%(i+1,j+1), format='ascii')
            if j == 0:
                Sum_X = tb['X']
                Sum_Y = tb['Y']
            else:
                Sum_Y += tb['Y']
            
            #for k in range(len(tb)):
            #    Plot_engine.plot_xy([k+1],tb['Y'][k], color=Plot_lib_color, alpha=Plot_lib_alpha, overplot=True)
         
        Min_chi2_for_plot = numpy.power(10, Min_chi2_log-(Max_chi2_log-Min_chi2_log)*0.05)
        Max_chi2_for_plot = numpy.power(10, Max_chi2_log+(Max_chi2_log-Min_chi2_log)*0.85)
        Plot_chi2_color = Plot_engine.get_color_by_value([Min_chi2_for_plot, Max_chi2_for_plot], 
                                                    input_value=Cut_chi2_array[i], 
                                                    log=1, 
                                                    cmap=matplotlib.cm.get_cmap('gray'))
        #print('Plot_chi2_color: ', Plot_chi2_color)
        # 
        # 
        # plot total LVG
        #Plot_engine.plot_xy([k+1],tb['Y'][k], color=Plot_chi2_color, alpha=Plot_chi2_alpha, overplot=True)
        for k in range(len(Line_transition_obs)):
            match_obs_index = numpy.argwhere(numpy.abs(Sum_X-Line_transition_obs[k])<1e-6).flatten()
            if len(match_obs_index) > 0:
                match_obs_index = match_obs_index[0]
                Plot_engine.plot_xy([k+1], Sum_Y[match_obs_index], color=Plot_chi2_color, alpha=Plot_chi2_alpha, overplot=True)
        # 
        # 
        # count++
        Count_plot_chi2 = Count_plot_chi2 + 1
        text_starting_pos = 0.92 # 0.95 # 20221103
        text_line_spacing = 0.05 # 0.03 # 20221103
        # 
        # 
        # show chi2 text on the figure
        if not SetOnlyPlotBestSED:
            if i == Cut_chi2_array_size-1:
                Plot_engine.xyouts(0.05, text_starting_pos, r'$\chi^2:$', NormalizedCoordinate=True, useTex=True)
            if i == 0:
                Plot_engine.xyouts(0.10, text_starting_pos-text_line_spacing*(Count_label_chi2), '......', NormalizedCoordinate=True, color=Plot_chi2_color) # i == 0 is the minimum chisq
                Count_label_chi2 = Count_label_chi2 + 1
            if Count_plot_chi2 % int((Cut_chi2_array_size/7)+1) == 0 or i == 0 or i == Cut_chi2_array_size-1:
                Plot_engine.xyouts(0.10, text_starting_pos-text_line_spacing*(Count_label_chi2), '%.1f'%(Cut_chi2_array[i]), NormalizedCoordinate=True, useTex=True, color=Plot_chi2_color)
                Count_label_chi2 = Count_label_chi2 + 1
        # 
        # 
        # show reduced-chi2 text on the figure
        if not SetOnlyPlotBestSED:
            if i == Cut_chi2_array_size-1:
                Plot_engine.xyouts(0.05+0.13, text_starting_pos, r'$\chi_{r.}^2:$', NormalizedCoordinate=True, useTex=True)
            if i == 0:
                Plot_engine.xyouts(0.10+0.13, text_starting_pos-text_line_spacing*(Count_label_rchi2), '......', NormalizedCoordinate=True, color=Plot_chi2_color) # i == 0 is the minimum chisq
                Count_label_rchi2 = Count_label_rchi2 + 1
            if Count_plot_chi2 % int((Cut_chi2_array_size/7)+1) == 0 or i == 0 or i == Cut_chi2_array_size-1:
                Plot_engine.xyouts(0.10+0.13, text_starting_pos-text_line_spacing*(Count_label_rchi2), '%.2f'%(Cut_chi2_array[i]/numpy.count_nonzero(Detection_mask)), NormalizedCoordinate=True, useTex=True, color=Plot_chi2_color)
                Count_label_rchi2 = Count_label_rchi2 + 1
        # 
        # 
        # show redshift (z) and source name on the figure
        if not SetOnlyPlotBestSED:
            # if plot a range of solutions with chi-square lower than then minimum_chisq + delta_chisq
            if i == 0:
                # 
                Plot_engine.xyouts(0.05+0.13+0.13, text_starting_pos, r'$z=%s$'%(Redshift), NormalizedCoordinate=True, useTex=True)
                PlotTextPosY = text_starting_pos
                # 
                if SourceName != '':
                    #Plot_engine.xyouts(0.97, 0.90, SourceName, NormalizedCoordinate=True, fontsize=16, horizontalalignment='right') # 20221103
                    #PlotTextPosY = 0.90 # 20221103
                    Plot_engine.xyouts(0.94, 0.88, SourceName, NormalizedCoordinate=True, fontsize=16, horizontalalignment='right') # 20221103
                    PlotTextPosY = 0.88 # 20221103
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
                PlotTextPosY = 0.95
                # 
                if SourceName != '':
                    #Plot_engine.xyouts(0.05, 0.90, SourceName, NormalizedCoordinate=True, fontsize=15) # 20221103
                    #PlotTextPosY = 0.90 # 20221103
                    Plot_engine.xyouts(0.05, 0.88, SourceName, NormalizedCoordinate=True, fontsize=15) # 20221103
                    PlotTextPosY = 0.88 # 20221103
                # 
                #Plot_engine.xyouts(0.20, 0.90, '$z=%s$'%(Redshift), NormalizedCoordinate=True, useTex=True, fontsize=15)
                #<20180216># allow user input text with the "-text" argument
                if len(UserInputText) > 0:
                    for UserInputTextIndex in range(len(UserInputText)):
                        UserInputTextUseTeX = (UserInputText[UserInputTextIndex].find('$')>=0)
                        Plot_engine.xyouts(0.05, PlotTextPosY-0.05*(UserInputTextIndex+1), UserInputText[UserInputTextIndex], NormalizedCoordinate=True, useTex=UserInputTextUseTeX, fontsize=15)
            
            # when only plotting the best solution, we also read the "dump/1/LVG_SUM" (wavelength is restframe) and convert it to obsframe (only wavelength) and save it as an output. 
            if i == 0:
                best_LVG_data_table = asciitable.read(Output_dir+'dump'+os.sep+'%d/LVG_SUM'%(i+1))
                with open(Output_dir+'dump'+os.sep+'%d/redshift.txt'%(i+1), 'r') as ifp:
                    best_LVG_redshift = float(ifp.read())
                with open(Output_dir+'best-fit_LVG_%s.txt'%(SourceName), 'w') as ofp:
                    ofp.write('# %25s %25s\n'%('line_transition', 'flux_obsframe_Jykms'))
                    for isedrow in range(len(best_LVG_data_table)):
                        ofp.write('%27.10f %25.10e\n'%(best_LVG_data_table.field(0)[isedrow]*(1.0+best_LVG_redshift), best_LVG_data_table.field(1)[isedrow]))
                    print('Output to "%s"!'%(Output_dir+'best-fit_LVG_%s.txt'%(SourceName)))
            
        #break
    # 
    # 
    # 
    # 
    # 
    # 
    # Plot OBS data points
    try:
        X_tick_list = []
        X_ticklabel_list = []
        for k in range(len(Line_transition_obs)):
            ID_species, Name_species, J_upper, J_lower = convert_X_species_to_ID_Name_J_upper_J_lower(Line_transition_obs[k])
            if Detection_mask[k] == True:
                Plot_engine.plot_xy([k+1], Flux_obs[k], yerr=FluxErr_obs[k], dataname='obs', overplot=True, symbol='open square', symsize=3, thick=1.5, capsize=4, zorder=10, verbose=UserInputVerbose)
            else:
                Plot_engine.plot_xy([k+1], 3.0*FluxErr_obs[k], dataname='upper limits', overplot=True, symbol='upper limits', symsize=3, thick=1.25, alpha=0.5, zorder=9, verbose=UserInputVerbose)
            X_tick_list.append(k+1)
            X_ticklabel_list.append('%s(%d-%d)'%(Name_species, J_upper, J_lower))
        plt.xticks(X_tick_list, X_ticklabel_list, rotation=45) # rotation=75 # 20221103
    except Exception as err:
        print(err)
    # 
    # 
    #Plot_engine.set_xrange([0.1,1e6])
    #Plot_engine.set_yrange([1e-6,1e4])
    if len(PlotYRange) == 2:
        Plot_engine.set_yrange(PlotYRange)
    else:
        Plot_engine.expand_yrange(0.2)
    #Plot_engine.set_xtitle('X-species [ID,J_u,J_l]')
    Plot_engine.set_ytitle('Flux density [Jy km s$^{-1}$]')
    Plot_engine.set_xcharsize(charsize=12, axislabelcharsize=16)
    Plot_engine.set_ycharsize(charsize=12, axislabelcharsize=16)
    Plot_engine.savepdf(Output_dir+Output_name+'.pdf')
    Plot_engine.savefig(Output_dir+Output_name+'.png')
    #Plot_engine.show()
    Plot_engine.close()
    print('')
    #print('Output to "%s"!'%(Output_dir+Output_name+'.pdf'))
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
    Tkin_dict_list = []
    nH2_dict_list = []
    NH2_dict_list = []
    MH2_dict_list = []
    XCICO_dict_list = []
    # 
    # define constants
    pi = numpy.pi
    #dL = 1.0/(4*pi) # if Redshift is not given, then we do not apply 4*pi*dL**2 to the quantities below. 
    #Redshift = float(InfoDict['REDSHIFT'])
    AreaInKpcSquare = 1.0 # Libs are made with M_H2 from N_H2 within 1.0 kpc^2 area.
    NormalizationFactor = 1.0
    print('')
    print('')
    print('z = %s'%(Redshift))
    if UserInputDistance > 0.0:
        dL = UserInputDistance
        print('dL = %s [Mpc] (user input)'%(dL))
    elif Redshift > 0.0:
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
        print('dL = %s [Mpc]'%(dL))
    else:
        dL = 10.0 # Mpc # 20220923 use 10 Mpc for reference
        print('dL = %s [Mpc] (default value, use with caution!)'%(dL))
    # 
    if UserInputArea > 0.0:
        AreaInKpcSquare = UserInputArea
        print('area = %s [kpc^2] (user input)'%(AreaInKpcSquare))
    else:
        print('area = %s [kpc^2] (default value, use with caution!)'%(AreaInKpcSquare))
    # 
    if UserInputNormalization > 0.0:
        NormalizationFactor = UserInputNormalization
        print('normalization = %s (user input)'%(NormalizationFactor))
    # 
    # get parameter list for each lib
    Lib_number = int(InfoDict['NLIB'])
    Lib_params = {}
    Num_params = {}
    Col_number = 2+2*Lib_number+1 # this quantity indicates the column number in the big chi-square data table. The first two columns are always i0 and chi2, so the LIB params start from column 3.
    # 
    for j in range(Lib_number):
        # 
        Tkin_dict = {}
        nH2_dict = {}
        NH2_dict = {}
        MH2_dict = {}
        XCICO_dict = {}
        # 
        Lib_name = 'LIB%d'%(j+1)
        Lib_dict = CrabTableReadInfo(InfoDict[Lib_name], verbose=0)
        #print(Lib_dict)
        # 
        # read the number of parameters from the LVG LIB file
        Key_NPAR = '# NPAR'
        if Key_NPAR in Lib_dict:
            Num_params[Lib_name] = int(Lib_dict[Key_NPAR])
        else:
            print("Error! \"%s\" was not found in \"%s\"!"%(Key_NPAR, InfoDict[Lib_name]))
            sys.exit()
        # 
        # read the title of each parameter from the LVG LIB file
        Lib_params[Lib_name] = []
        for k in range(Num_params[Lib_name]):
            Key_TPAR = '# TPAR%d'%(k+1)
            if Key_TPAR in Lib_dict:
                Lib_params[Lib_name].append(Lib_dict[Key_TPAR])
            else:
                print("Error! \"%s\" was not found in \"%s\"!"%(Key_TPAR, InfoDict[Lib_name]))
                sys.exit()
            # 
            # check multi-lib suffix
            if Lib_number > 1:
                Par_suffix = '_%d'%(j+1)
            else:
                Par_suffix = ''
            # 
            # check Tkin nH2
            if InfoDict[Lib_name].endswith('.lvg'):
                if 'T_{kin}' == Lib_dict[Key_TPAR]:
                    Tkin_dict['Lib_file'] = InfoDict[Lib_name]
                    Tkin_dict['Lib_name'] = Lib_name
                    Tkin_dict['Lib_numb'] = j+1
                    Tkin_dict['Par_name'] = r'$T_{\mathrm{kin}}$ [$\mathrm{K}$]' + Par_suffix # Lib_dict[Key_TPAR]
                    Tkin_dict['Par_file'] = 'Tkin' + Par_suffix
                    Tkin_dict['Col_numb'] = Col_number
                    Tkin_dict['Log_calc'] = False
                    Tkin_dict['range'] = [5.0,300]
                    Tkin_dict['value'] = DataTable.getColumn(Col_number)
                    Tkin_dict['chisq'] = DataArray['chi2']
                    # 
                elif 'n_{H_2}' == Lib_dict[Key_TPAR]:
                    nH2_dict['Lib_file'] = InfoDict[Lib_name]
                    nH2_dict['Lib_name'] = Lib_name
                    nH2_dict['Lib_numb'] = j+1
                    nH2_dict['Par_name'] = r'$\log \ n_{\mathrm{H_2}}$ [$\mathrm{cm}^{-3}$]' + Par_suffix # Lib_dict[Key_TPAR]
                    nH2_dict['Par_file'] = 'nH2' + Par_suffix
                    nH2_dict['Col_numb'] = Col_number
                    nH2_dict['Log_calc'] = True
                    nH2_dict['range'] = numpy.power(10,[2.0, 6.5])
                    nH2_dict['value'] = DataTable.getColumn(Col_number)
                    nH2_dict['chisq'] = DataArray['chi2']
                    # 
                elif 'N_{H_2}' == Lib_dict[Key_TPAR]:
                    NH2_dict['Lib_file'] = InfoDict[Lib_name]
                    NH2_dict['Lib_name'] = Lib_name
                    NH2_dict['Lib_numb'] = j+1
                    NH2_dict['Par_name'] = r'$\log \ N_{\mathrm{H_2}}$ [$\mathrm{cm^{-2}}$]' + Par_suffix # Lib_dict[Key_TPAR]
                    NH2_dict['Par_file'] = 'N_H2' + Par_suffix
                    NH2_dict['Col_numb'] = Col_number
                    NH2_dict['Log_calc'] = True
                    NH2_dict['range'] = numpy.power(10,[19.0, 25.5])
                    NH2_dict['value'] = DataTable.getColumn(Col_number)
                    NH2_dict['chisq'] = DataArray['chi2']
                    # 
                    # MH2_dict['Lib_file'] = InfoDict[Lib_name]
                    # MH2_dict['Lib_name'] = Lib_name
                    # MH2_dict['Lib_numb'] = j+1
                    # MH2_dict['Par_name'] = r'$\log \ M_{\mathrm{H_2}}$ [$\mathrm{M_{\odot}}$]' + Par_suffix # Lib_dict[Key_TPAR]
                    # MH2_dict['Par_file'] = 'M_H2' + Par_suffix
                    # MH2_dict['Col_numb'] = Col_number
                    # MH2_dict['Log_calc'] = True
                    # MH2_dict['range'] = numpy.power(10,[4.0, 12.5])
                    # MH2_dict['value'] = DataArray['a%d'%(j+1)] * AreaInKpcSquare * 9.52140614e+42 * NH2_dict['value'] * 1.6828295e-57 / (dL/10)**2 * NormalizationFactor
                    # MH2_dict['chisq'] = DataArray['chi2']
                    # Notes: 
                    #   9.52140614e+42 is 1 kpc^2 in units of cm^2
                    #   1.6828295e-57 is the H2 mass (m_p*2+m_e) in solar mass unit
                    #   dL/10 because 10 Mpc is the fiducial distance used in library
                    #   1 kpc^2 is the fiducial area used in library
                    # 
                # elif 'M_{H_2}' == Lib_dict[Key_TPAR]: # bug in lib file, mass is incorrect
                #     MH2_dict['Lib_file'] = InfoDict[Lib_name]
                #     MH2_dict['Lib_name'] = Lib_name
                #     MH2_dict['Lib_numb'] = j+1
                #     MH2_dict['Par_name'] = r'$\log \ M_{\mathrm{H_2}}$ [$\mathrm{M_{\odot}}$]' + Par_suffix # Lib_dict[Key_TPAR]
                #     MH2_dict['Par_file'] = 'MH2' + Par_suffix
                #     MH2_dict['Col_numb'] = Col_number
                #     MH2_dict['Log_calc'] = True
                #     MH2_dict['range'] = numpy.power(10,[6.0, 12.5])
                #     MH2_dict['value'] = DataArray['a%d'%(j+1)] * DataTable.getColumn(Col_number) * AreaInKpcSquare * NormalizationFactor / (dL/10)**2
                #     MH2_dict['chisq'] = DataArray['chi2']
                elif 'X_{CICO}' == Lib_dict[Key_TPAR]: 
                    XCICO_dict['Lib_file'] = InfoDict[Lib_name]
                    XCICO_dict['Lib_name'] = Lib_name
                    XCICO_dict['Lib_numb'] = j+1
                    XCICO_dict['Par_name'] = r'[$\mathrm{C{\tt{I}}/CO}$]' + Par_suffix # Lib_dict[Key_TPAR]
                    XCICO_dict['Par_file'] = 'XCICO' + Par_suffix
                    XCICO_dict['Col_numb'] = Col_number
                    XCICO_dict['Log_calc'] = False
                    XCICO_dict['range'] = [0.0, 3.0]
                    XCICO_dict['value'] = DataTable.getColumn(Col_number)
                    XCICO_dict['chisq'] = DataArray['chi2']
            # 
            # finished checking library properties
            # 
            # count the column number in the big chi-square data table. All params in each LIB are listed in the big chi-square data table.
            Col_number = Col_number + 1
        # 
        Tkin_dict_list.append(Tkin_dict)
        nH2_dict_list.append(nH2_dict)
        NH2_dict_list.append(NH2_dict)
        MH2_dict_list.append(MH2_dict)
        XCICO_dict_list.append(XCICO_dict)
    # 
    # Total MH2
    Sum_MH2_dict = {}
    if Lib_number >= 2:
        for j in range(Lib_number):
            if len(MH2_dict_list[j]) > 0:
                if len(Sum_MH2_dict) == 0:
                    #print("Sum_MH2_dict['value'] = MH2_dict_list[%d]['value']"%(j))
                    #print(MH2_dict_list[j]['value'])
                    Sum_MH2_dict = deepcopy(MH2_dict_list[j])
                    Sum_MH2_dict['Lib_file'] = [Sum_MH2_dict['Lib_file']]
                    Sum_MH2_dict['Lib_name'] = [Sum_MH2_dict['Lib_name']]
                    Sum_MH2_dict['Col_numb'] = [Sum_MH2_dict['Col_numb']]
                else:
                    #print("Sum_MH2_dict['value'] += MH2_dict_list[%d]['value']"%(j))
                    #print(MH2_dict_list[j]['value'])
                    Sum_MH2_dict['value'] += MH2_dict_list[j]['value']
                    Sum_MH2_dict['Lib_file'].append(MH2_dict_list[j]['Lib_file'])
                    Sum_MH2_dict['Lib_name'].append(MH2_dict_list[j]['Lib_name'])
                    Sum_MH2_dict['Col_numb'].append(MH2_dict_list[j]['Col_numb'])
                    Sum_MH2_dict['Par_name'] = r'$\log \ M_{\mathrm{H_2}}$ (total) [$\mathrm{M}_{\odot}$]'
                    Sum_MH2_dict['Par_file'] = 'MH2_total'
    # 
    # analyze 
    print('Num_params', Num_params)
    print('Lib_params', Lib_params)
    Plot_engine = CrabPlot(figure_size=(8.0,5.5)) # figure_size=(14.0,10.0) # 20221103
    Plot_engine.set_margin(panel=0, top=0.99, bottom=0.11, left=0.10, right=0.99) # top=0.96, bottom=0.08, left=0.06, right=0.96 # 20221103
    Plot_engine.Plot_device.subplots_adjust(wspace=0.25, hspace=0.3) # 20221103
    for j in range(Lib_number):
        if 'value' in Tkin_dict_list[j]:
            analyze_chisq_distribution(Tkin_dict_list[j], Plot_engine = Plot_engine, Output_dir = Output_dir)
        if 'value' in nH2_dict_list[j]:
            analyze_chisq_distribution(nH2_dict_list[j], Plot_engine = Plot_engine, Output_dir = Output_dir)
        if 'value' in NH2_dict_list[j]:
            analyze_chisq_distribution(NH2_dict_list[j], Plot_engine = Plot_engine, Output_dir = Output_dir)
        if 'value' in MH2_dict_list[j]:
            analyze_chisq_distribution(MH2_dict_list[j], Plot_engine = Plot_engine, Output_dir = Output_dir)
        if 'value' in XCICO_dict_list[j]:
            analyze_chisq_distribution(XCICO_dict_list[j], Plot_engine = Plot_engine, Output_dir = Output_dir)
    if 'value' in Sum_MH2_dict:
        analyze_chisq_distribution(Sum_MH2_dict, Plot_engine = Plot_engine, Output_dir = Output_dir)
    Plot_engine.set_xcharsize(panel=0, charsize=11, axislabelcharsize=16) # all panels
    Plot_engine.set_ycharsize(panel=0, charsize=11, axislabelcharsize=16) # all panels
    Plot_engine.savepdf(Output_dir+Output_name+'.chisq.pdf')
    Plot_engine.savefig(Output_dir+Output_name+'.chisq.png')
    #Plot_engine.show()
    Plot_engine.close()
    #print('Output to "%s"!'%(Output_dir+Output_name+'.chisq.pdf'))
    # 
    






