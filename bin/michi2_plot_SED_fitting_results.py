#!/usr/bin/env python
# 

import os
import sys

#os.system('cp /Users/dzliu/Softwares/Python/lib/crab/crabplot/CrabPlot.py /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/lib/python/crabplot/CrabPlot.py')

sys.path.insert(1, os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabtable')
sys.path.insert(1, os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabplot')
#sys.path.insert(1, '/Users/dzliu/Softwares/Python/lib/crab/crabtable')
#sys.path.insert(1, '/Users/dzliu/Softwares/Python/lib/crab/crabplot')

from CrabTable import *
from CrabPlot import *

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

Delta_chisq_of_interest = 5.89





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
            print(param_dict)
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
        param_array_nonan = copy(param_array)
        param_array_nonan[numpy.isnan(param_array_nonan)] = -99
        # 
        # apply range -- no, do not cut the array, but just adjust the plotting range.
        #if 'range' in param_dict:
        #    if len(param_dict['range'])>=2:
        #        param_range_mask = (param_array >= param_dict['range'][0]) & (param_array <= param_dict['range'][1])
        #        chisq_array = chisq_array[param_range_mask]
        #        param_array = param_array[param_range_mask]
        #        param_array_nonan = param_array_nonan[param_range_mask]
        # 
        # bin param
        param_log = False
        if 'Log_calc' in param_dict:
            if param_dict['Log_calc'] == True:
                param_log = True
        if param_log == True: 
            param_log_mask = (param_array>0)
            param_log_mask2 = (param_array<=0)
            param_array[param_log_mask] = numpy.log10(param_array[param_log_mask])
            param_array[param_log_mask2] = numpy.nan
            param_array_nonan[param_log_mask] = numpy.log10(param_array_nonan[param_log_mask]) # NoNaN array
        param_min = numpy.nanmin(param_array)
        param_max = numpy.nanmax(param_array)
        # 
        # set xrange to the user-specified values
        xrange = None
        if 'range' in param_dict:
            if len(param_dict['range'])>=2:
                param_min = param_dict['range'][0]
                param_max = param_dict['range'][1]
                if param_log == True: 
                    param_min = numpy.log10(param_min)
                    param_max = numpy.log10(param_max)
                xrange = [param_min, param_max]
        # 
        param_bin_numb = 50 #<TODO># 
        param_bin_step = (param_max-param_min)/param_bin_numb
        print('Param value min max are %s %s'%(param_min, param_max))
        # 
        # optimize xrange
        xclip_mask = (chisq_array <= chisq_min+Delta_chisq_of_interest)
        if len(numpy.argwhere(xclip_mask)) > 0:
            xclip_min = numpy.nanmin(param_array[xclip_mask])
            xclip_max = numpy.nanmax(param_array[xclip_mask])
            if xclip_max > xclip_min:
                xrange = [xclip_min, xclip_max]
                xrange = [xrange[0]-(xrange[1]-xrange[0])*1.0, xrange[1]+(xrange[1]-xrange[0])*1.0]
                param_bin_step = (xclip_max - xclip_min) / 15.0 #<TODO># bin step
                param_bin_numb = int((param_max-param_min)/param_bin_step)+1
        # 
        # optimize yrange
        yrange = [1/(chisq_min+Delta_chisq_of_interest), 1/(chisq_min)]
        yrange = [yrange[0]-(yrange[1]-yrange[0])*0.25, yrange[1]+(yrange[1]-yrange[0])*0.35]
        # 
        xlog = None
        #if 'Log_plot' in param_dict:
        #    if param_dict['Log_plot'] == True:
        #        xlog = 1 # not working for matplotlib bar plot (i.e., CrabPlot plot_hist)!
        # 
        ylog = None
        # 
        # compute minimum chisq in each parameter bin
        param_bin_x = []
        param_bin_y = []
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
            # 
            param_bin_args = numpy.argwhere(param_bin_mask)
            if verbose>=2:
                print('Binning from param value %s to %s, step %s, count %d'%(param_bin_lower, param_bin_upper, param_bin_step, len(param_bin_args)))
            if len(param_bin_args) > 0:
                chisq_bin_array = chisq_array[param_bin_mask]
                chisq_bin_min = numpy.nanmin(chisq_bin_array)
                chisq_bin_max = numpy.nanmax(chisq_bin_array)
                param_bin_x.append(param_min+i*param_bin_step)
                param_bin_y.append(chisq_bin_min)
        # 
        # Initialize a plot
        if Plot_engine is None:
            Plot_engine = CrabPlot(figure_size=(9.0,5.0))
            Plot_engine.set_margin(panel=0, top=0.96)
        # 
        # Plot xy
        Plot_engine.plot_xy(param_array, 1/numpy.array(chisq_array), overplot = False, 
                                xtitle = param_dict['Par_name'], ytitle = '$1/\chi^2$', useTex = True, 
                                size = 0.12, symbol = 'cross', xrange = xrange, yrange = yrange, xlog = xlog, ylog = ylog)
        # 
        # Plot Cut_chi2
        Plot_engine.plot_line(xrange[0], 1/(chisq_min+Delta_chisq_of_interest), xrange[1], 1/(chisq_min+Delta_chisq_of_interest), overplot = True, linestyle = 'dashed')
        # 
        # Plot histogram
        Plot_engine.plot_hist(param_bin_x, 1/numpy.array(param_bin_y), width = param_bin_step, align = 'edge', overplot = False, 
                                xtitle = param_dict['Par_name'], ytitle = '$1/\chi^2$', useTex = True, 
                                xrange = xrange, yrange = yrange, xlog = xlog, ylog = ylog)
        # 
        # Plot Cut_chi2
        Plot_engine.plot_line(xrange[0], 1/(chisq_min+Delta_chisq_of_interest), xrange[1], 1/(chisq_min+Delta_chisq_of_interest), overplot = True, linestyle = 'dashed')
        # 
    else:
        print('Error! analyze_chisq_distribution() got unaccepted inputs!')
        sys.exit()












##########################################
#               MAIN PROGRAM             #
##########################################

if len(sys.argv) == 1:
    
    print('Usage: michi2_plot_SED_fitting_results.py fit_5.out')
    sys.exit()

else:
    # 
    # Read chi2 table
    DataFile = sys.argv[1]
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
        #    Data_chi2 = DataTable.getColumn('col%d'%(i+1))
        for j in range(int(InfoDict['NLIB'])):
            if DataHeaders[i] == 'i%d'%(j+1):
                DataArray['i%d'%(j+1)] = DataTable.getColumn('col%d'%(i+1))
            elif DataHeaders[i] == 'a%d'%(j+1):
                DataArray['a%d'%(j+1)] = DataTable.getColumn('col%d'%(i+1))
            elif DataHeaders[i] == 'i0':
                DataArray['i0'] = DataTable.getColumn('col%d'%(i+1))
            elif DataHeaders[i] == 'a0':
                DataArray['a0'] = DataTable.getColumn('col%d'%(i+1))
            elif DataHeaders[i] == 'chi2':
                DataArray['chi2'] = DataTable.getColumn('col%d'%(i+1))
    # 
    # Sort chi2 table
    #print(DataTable.TableHeaders)
    #print(DataArray['chi2'])
    SortedIndex = numpy.argsort(DataArray['chi2'])
    #SelectNumber = 10 #<TODO># how many SEDs to show?
    #SelectNumber = len(SortedIndex) if len(SortedIndex)<SelectNumber else SelectNumber
    #SelectIndex = SortedIndex[0:SelectNumber-1]
    ##print(DataTable.TableData[SelectIndex])
    #Arr_chi2 = DataArray['chi2'][SortedIndex] # the sorted chi2 array
    #Min_chi2 = numpy.nanmin(Arr_chi2[0:SelectNumber-1])
    #Max_chi2 = numpy.nanmax(Arr_chi2[0:SelectNumber-1])
    Cut_chi2 = numpy.nanmin(DataArray['chi2']) + Delta_chisq_of_interest # 5 parameter, 68% confidence
    SelectNumber = len(numpy.argwhere(DataArray['chi2']<=Cut_chi2))
    SelectNumber = len(SortedIndex) if len(SortedIndex)<SelectNumber else SelectNumber # do not exceed the total number of chi2
    Arr_chi2 = DataArray['chi2'][SortedIndex[0:SelectNumber]] # the sorted chi2 array, note that when selecting subscript/index with :, the upper index is not included.
    Min_chi2 = numpy.nanmin(Arr_chi2)
    Max_chi2 = numpy.nanmax(Arr_chi2)
    print('Selecting %d chi2 solutions with chi2 <= min(chi2)+%s'%(SelectNumber, Delta_chisq_of_interest))
    # 
    MaxPlotNumber = 50 # Max chi2 solutions to plot, we plot the first MaxPlotNumber/2 and the last MaxPlotNumber/2 solutions, skip solutions in the middle.
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
    if True:
        # 
        # Get SED
        for i in range(SelectNumber):
            # 
            # skip solutions between 11th to last 11th.
            if i > MaxPlotNumber/2 and i<(SelectNumber-1-MaxPlotNumber/2):
                continue
            # 
            # Read SED_LIB
            for j in range(int(InfoDict['NLIB'])):
                if not os.path.isdir('obj_%d'%(i+1)):
                    os.mkdir('obj_%d'%(i+1))
                # 
                if not os.path.isfile('obj_%d/SED_LIB%d'%(i+1,j+1)):
                    #BashCommand = 'cd obj_%d/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_SED ../%s %d %s SED_LIB%d'%\
                    #                    (i+1, \
                    #                        InfoDict['LIB%d'%(j+1)], \
                    #                            DataArray['i%d'%(j+1)][SortedIndex[i]], \
                    #                                DataArray['a%d'%(j+1)][SortedIndex[i]], \
                    #                                    j+1)
                    #print(BashCommand)
                    #os.system(BashCommand)
                    # 
                    # do python way 20180113
                    BashCommand = '%s/michi2_read_lib_SEDs.py %s %d obj_%d | tee obj_%d/log.txt'%\
                                    (os.path.dirname(os.path.realpath(__file__)), \
                                        DataFile, \
                                            SortedIndex[i]+1, \
                                                i+1, \
                                                    i+1)
                    print(BashCommand)
                    os.system(BashCommand)
                    BashCommand = 'echo "%s" > obj_%d/chi2.txt'%\
                                    (DataArray['chi2'][SortedIndex[i]], \
                                                i+1)
                    print(BashCommand)
                    os.system(BashCommand)
                    BashCommand = 'echo "%s" > obj_%d/line_number.txt'%\
                                    (SortedIndex[i]+1, \
                                                i+1)
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
        # Wait for a long time
        # 
        # Then plot SEDs
        Redshift = float(InfoDict['REDSHIFT'])
        Color_list = ['cyan', 'gold', 'red', 'blue', 'purple']
        Plot_engine = CrabPlot(figure_size=(9.0,5.0))
        for i in range(SelectNumber-1,-1,-1):
            # 
            # skip solutions between 11th to last 11th.
            if i > MaxPlotNumber/2 and i<(SelectNumber-1-MaxPlotNumber/2):
                continue
            # 
            # alpha by chi2
            print('Plotting chi2=%s obj_%d'%(Arr_chi2[i], i+1))
            Min_chi2_log = numpy.log10(Min_chi2)
            Max_chi2_log = numpy.log10(Max_chi2)
            Min_chi2_for_plot = numpy.power(10, Min_chi2_log-(Max_chi2_log-Min_chi2_log)*0.8)
            Max_chi2_for_plot = numpy.power(10, Max_chi2_log+(Max_chi2_log-Min_chi2_log)*0.3)
            Alpha_chi2 = Plot_engine.get_color_by_value([Min_chi2_for_plot, Max_chi2_for_plot], 
                                                        input_value=Arr_chi2[i], 
                                                        log=1, 
                                                        cmap=matplotlib.cm.get_cmap('gray_r'))[0]
            print('Alpha_chi2: ', Alpha_chi2)
            # 
            # 
            for j in range(int(InfoDict['NLIB'])):
                linewidth = numpy.sqrt(1.44/float(SelectNumber)) #<TODO># tune line width
                xclip = None
                if j == 0: xclip = [(50,numpy.inf)]
                elif j == 4: xclip = [(-numpy.inf,2e3)]
                Plot_engine.plot_data_file('obj_%d/SED_LIB%d'%(i+1,j+1), xlog=1, ylog=1, xclip=xclip, current=1, \
                                    dataname='obj_%d_SED_LIB%d'%(i+1,j+1), 
                                    redshift = Redshift, 
                                    linestyle='dashed', linewidth=linewidth, color=Color_list[j], alpha=Alpha_chi2)
            # 
            #<20180114><splined else where># obj_SED_1 = Plot_engine.Plot_data['obj_%d_SED_LIB1'%(i+1)]
            #<20180114><splined else where># obj_SED_2 = Plot_engine.Plot_data['obj_%d_SED_LIB2'%(i+1)]
            #<20180114><splined else where># obj_SED_3 = Plot_engine.Plot_data['obj_%d_SED_LIB3'%(i+1)]
            #<20180114><splined else where># obj_SED_4 = Plot_engine.Plot_data['obj_%d_SED_LIB4'%(i+1)]
            #<20180114><splined else where># obj_SED_5 = Plot_engine.Plot_data['obj_%d_SED_LIB5'%(i+1)]
            #<20180114><splined else where># obj_SED_sum_x_lg = numpy.arange(-2,6,0.001) # wavelength_um grid
            #<20180114><splined else where># obj_SED_sum_x = numpy.power(10, obj_SED_sum_x_lg) # make it in linear space
            #<20180114><splined else where># obj_SED_sum_y = []
            #<20180114><splined else where># for j in range(int(InfoDict['NLIB'])):
            #<20180114><splined else where>#     obj_SED_single_x = Plot_engine.Plot_data['obj_%d_SED_LIB%d'%(i+1,j+1)][:,0]
            #<20180114><splined else where>#     obj_SED_single_y = Plot_engine.Plot_data['obj_%d_SED_LIB%d'%(i+1,j+1)][:,1]
            #<20180114><splined else where>#     obj_SED_spline_y = Plot_engine.spline(obj_SED_single_x, obj_SED_single_y, obj_SED_sum_x_lg, xlog=1, ylog=1, outputxlog=0)
            #<20180114><splined else where>#     if j == 0:
            #<20180114><splined else where>#         obj_SED_sum_y = obj_SED_spline_y
            #<20180114><splined else where>#     else:
            #<20180114><splined else where>#         obj_SED_sum_y = obj_SED_sum_y + obj_SED_spline_y
            # 
            #pprint(obj_SED_sum_y)
            #print(numpy.column_stack((obj_SED_sum_x, obj_SED_sum_y)))
            #print(Arr_chi2[i])
            Min_chi2_for_plot = numpy.power(10, Min_chi2_log-(Max_chi2_log-Min_chi2_log)*0.05)
            Max_chi2_for_plot = numpy.power(10, Max_chi2_log+(Max_chi2_log-Min_chi2_log)*0.85)
            Color_chi2 = Plot_engine.get_color_by_value([Min_chi2_for_plot, Max_chi2_for_plot], 
                                                        input_value=Arr_chi2[i], 
                                                        log=1, 
                                                        cmap=matplotlib.cm.get_cmap('gray'))
            #print('Color_chi2: ', Color_chi2)
            #Plot_engine.plot_line(obj_SED_sum_x, obj_SED_sum_y, current=1, color=Color_chi2)
            Plot_engine.plot_data_file('obj_%d/SED_SUM'%(i+1), xlog=1, ylog=1, current=1, \
                                dataname='obj_%d_SED_SUM'%(i+1), 
                                redshift = Redshift, 
                                linestyle='solid', linewidth=1.0, color=Color_chi2, alpha=1.0, zorder=8) # alpha=Alpha_chi2
            # 
            # show chi2 on the figure
            if i == SelectNumber-1:
                Plot_engine.xyouts(0.05, 0.95, '$\chi^2:$', NormalizedCoordinate=True, useTex=True)
            if i >= SelectNumber-6:
                print(0.09, 0.95-0.03*(SelectNumber-1-i), '%.1f'%(Arr_chi2[i]))
                Plot_engine.xyouts(0.09, 0.95-0.03*(SelectNumber-1-i), '%.1f'%(Arr_chi2[i]), NormalizedCoordinate=True, useTex=True, color=Color_chi2)
            if i == 0:
                Plot_engine.xyouts(0.09, 0.95-0.03*(6), '......', NormalizedCoordinate=True, color=Color_chi2)
                Plot_engine.xyouts(0.09, 0.95-0.03*(7), '%.1f'%(Arr_chi2[i]), NormalizedCoordinate=True, useTex=True, color=Color_chi2)
            if i == SelectNumber-1:
                Plot_engine.xyouts(0.15, 0.95, '$z=%s$'%(Redshift), NormalizedCoordinate=True, useTex=True)
            #break
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
            Plot_engine.plot_xy(Wavelength_obs[Detection_mask], Flux_obs[Detection_mask], yerr=FluxErr_obs[Detection_mask], dataname='obs', overplot=True, symbol='open square', symsize=3, thick=1.75, capsize=4, zorder=10)
            Plot_engine.plot_xy(Wavelength_obs[UpperLimits_mask], 3.0*FluxErr_obs[UpperLimits_mask], dataname='upper limits', overplot=True, symbol='upper limits', symsize=3, thick=1.25, alpha=0.5, zorder=9)
        except Exception as err:
            print(err)
        # 
        Plot_engine.set_xrange([0.1,1e6])
        Plot_engine.set_yrange([1e-6,1e4])
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
                    Stellar_mass_dict['Par_name'] = '$M_{*}$ [$\mathrm{M}_{\odot}$]' # Lib_dict[Key_TPAR]
                    Stellar_mass_dict['Col_numb'] = Col_number
                    Stellar_mass_dict['Log_calc'] = True
                    Stellar_mass_dict['range'] = numpy.power(10,[8.0,13.0])
                    Stellar_mass_dict['value'] = DataArray['a%d'%(j+1)] / (3.839e33*1e26/(4*pi*dL**2*9.52140e48)) * DataTable.getColumn('col%d'%(Col_number)) / (1+Redshift)
                    Stellar_mass_dict['chisq'] = DataArray['chi2']
            elif InfoDict[Lib_name].find('DL07.HiExCom') >= 0:
                if 'lgLTIR' == Lib_dict[Key_TPAR]:
                    LTIR_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LTIR_warm_dust_dict['Lib_name'] = Lib_name
                    LTIR_warm_dust_dict['Lib_numb'] = j+1
                    LTIR_warm_dust_dict['Par_name'] = '$L_{\mathrm{IR}}$ (warm) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LTIR_warm_dust_dict['Col_numb'] = Col_number
                    LTIR_warm_dust_dict['Log_calc'] = True
                    LTIR_warm_dust_dict['range'] = numpy.power(10,[9.0,14.0])
                    LTIR_warm_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn('col%d'%(Col_number))) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LTIR_warm_dust_dict['chisq'] = DataArray['chi2']
                elif 'Umin' == Lib_dict[Key_TPAR]:
                    Umin_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Umin_warm_dust_dict['Lib_name'] = Lib_name
                    Umin_warm_dust_dict['Lib_numb'] = j+1
                    Umin_warm_dust_dict['Par_name'] = '$U_{\mathrm{min}}$ (warm)' # Lib_dict[Key_TPAR]
                    Umin_warm_dust_dict['Col_numb'] = Col_number
                    Umin_warm_dust_dict['Log_plot'] = True # 'Log_plot', plot X axis in log scale
                    Umin_warm_dust_dict['range'] = [0.08,30.0]
                    Umin_warm_dust_dict['value'] = DataTable.getColumn('col%d'%(Col_number))
                    Umin_warm_dust_dict['chisq'] = DataArray['chi2']
                    # 
                    Mass_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Mass_warm_dust_dict['Lib_name'] = Lib_name
                    Mass_warm_dust_dict['Lib_numb'] = j+1
                    Mass_warm_dust_dict['Par_name'] = '$M_{\mathrm{dust}}$ (warm) [$\mathrm{M}_{\odot}$]'
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
                    LTIR_cold_dust_dict['Par_name'] = '$L_{\mathrm{IR}}$ (cold) [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    LTIR_cold_dust_dict['Col_numb'] = Col_number
                    LTIR_cold_dust_dict['Log_calc'] = True
                    LTIR_cold_dust_dict['range'] = numpy.power(10,[9.0,14.0])
                    LTIR_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn('col%d'%(Col_number))) * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LTIR_cold_dust_dict['chisq'] = DataArray['chi2']
                elif 'Umin' == Lib_dict[Key_TPAR]:
                    Umin_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Umin_cold_dust_dict['Lib_name'] = Lib_name
                    Umin_cold_dust_dict['Lib_numb'] = j+1
                    Umin_cold_dust_dict['Par_name'] = '$U_{\mathrm{min}}$ (cold)' # Lib_dict[Key_TPAR]
                    Umin_cold_dust_dict['Col_numb'] = Col_number
                    Umin_cold_dust_dict['Log_plot'] = True # 'Log_plot', plot X axis in log scale
                    Umin_cold_dust_dict['range'] = [0.08,30.0]
                    Umin_cold_dust_dict['value'] = DataTable.getColumn('col%d'%(Col_number))
                    Umin_cold_dust_dict['chisq'] = DataArray['chi2']
                    # 
                    Mass_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Mass_cold_dust_dict['Lib_name'] = Lib_name
                    Mass_cold_dust_dict['Lib_numb'] = j+1
                    Mass_cold_dust_dict['Par_name'] = '$M_{\mathrm{dust}}$ (cold) [$\mathrm{M}_{\odot}$]'
                    Mass_cold_dust_dict['Col_numb'] = 2+2*(j+1)
                    Mass_cold_dust_dict['Log_calc'] = True
                    Mass_cold_dust_dict['range'] = numpy.power(10,[7.0,12.0])
                    Mass_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * dL**2 / (1+Redshift) # Mdust # Mdust #NOTE# no need to multiply a '4*pi'!
                    Mass_cold_dust_dict['chisq'] = DataArray['chi2']
            elif InfoDict[Lib_name].find('MullaneyAGN') >= 0:
                # Mullaney AGN
                #   by integrating the AGN template (AGN_TYPE=2) from 1um to 1000um, we get an integration of 5133.913101
                if 'AGN_TYPE' == Lib_dict[Key_TPAR]:
                    Lumin_AGN_dict['Lib_file'] = InfoDict[Lib_name]
                    Lumin_AGN_dict['Lib_name'] = Lib_name
                    Lumin_AGN_dict['Lib_numb'] = j+1
                    Lumin_AGN_dict['Par_name'] = '$L_{\mathrm{AGN}}$ [$\mathrm{L}_{\odot}$]' # Lib_dict[Key_TPAR]
                    Lumin_AGN_dict['Col_numb'] = Col_number
                    Lumin_AGN_dict['Log_calc'] = True
                    Lumin_AGN_dict['range'] = numpy.power(10,[9.0,14.0])
                    Lumin_AGN_dict['value'] = DataArray['a%d'%(j+1)] * 5133.913101 * 4*pi*dL**2 / (1+Redshift) # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    Lumin_AGN_dict['chisq'] = DataArray['chi2']
                    # 
            # 
            # count the column number in the big chi-square data table. All params in each LIB are listed in the big chi-square data table.
            Col_number = Col_number + 1
    # 
    # Total LIR
    LTIR_total_dust_dict = copy(LTIR_warm_dust_dict)
    LTIR_total_dust_dict['value'] = LTIR_warm_dust_dict['value'] + LTIR_cold_dust_dict['value']
    LTIR_total_dust_dict['Lib_file'] = [LTIR_warm_dust_dict['Lib_file'], LTIR_cold_dust_dict['Lib_file']]
    LTIR_total_dust_dict['Lib_name'] = [LTIR_warm_dust_dict['Lib_name'], LTIR_cold_dust_dict['Lib_name']]
    LTIR_total_dust_dict['Col_numb'] = [LTIR_warm_dust_dict['Col_numb'], LTIR_cold_dust_dict['Col_numb']]
    LTIR_total_dust_dict['Par_name'] = '$L_{\mathrm{IR}}$ (total) [$\mathrm{L}_{\odot}$]'
    # 
    # Total Mdust
    Mass_total_dust_dict = copy(Mass_warm_dust_dict)
    Mass_total_dust_dict['value'] = Mass_warm_dust_dict['value'] + Mass_cold_dust_dict['value']
    Mass_total_dust_dict['Lib_file'] = [Mass_warm_dust_dict['Lib_file'], Mass_cold_dust_dict['Lib_file']]
    Mass_total_dust_dict['Lib_name'] = [Mass_warm_dust_dict['Lib_name'], Mass_cold_dust_dict['Lib_name']]
    Mass_total_dust_dict['Col_numb'] = [Mass_warm_dust_dict['Col_numb'], Mass_cold_dust_dict['Col_numb']]
    Mass_total_dust_dict['Par_name'] = '$M_{\mathrm{dust}}$ (total) [$\mathrm{M}_{\odot}$]'
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
    fPDR_total_dust_dict = copy(Umin_warm_dust_dict)
    fPDR_total_dust_dict['value'] = Mass_warm_dust_dict['value'] / (Mass_warm_dust_dict['value'] + Mass_cold_dust_dict['value'])
    fPDR_total_dust_dict['Lib_file'] = [Mass_warm_dust_dict['Lib_file'], Mass_cold_dust_dict['Lib_file']]
    fPDR_total_dust_dict['Lib_name'] = [Mass_warm_dust_dict['Lib_name'], Mass_cold_dust_dict['Lib_name']]
    fPDR_total_dust_dict['Col_numb'] = [Mass_warm_dust_dict['Col_numb'], Mass_cold_dust_dict['Col_numb']]
    fPDR_total_dust_dict['Par_name'] = '$\log \delta_{\mathrm{PDR}}$ (total)'
    fPDR_total_dust_dict['Log_calc'] = True
    fPDR_total_dust_dict['range'] = [1e-4,1.0]
    # 
    Umean_total_dust_dict = copy(Umin_warm_dust_dict)
    Umean_total_dust_dict['value'] = (1-fPDR_total_dust_dict['value']) * Umin_cold_dust_dict['value'] + fPDR_total_dust_dict['value'] * Umin_warm_dust_dict['value'] * (numpy.log(1e6/Umin_warm_dust_dict['value'])/(1-Umin_warm_dust_dict['value']/1e6))
    Umean_total_dust_dict['Lib_file'] = [Mass_warm_dust_dict['Lib_file'], Mass_cold_dust_dict['Lib_file']]
    Umean_total_dust_dict['Lib_name'] = [Mass_warm_dust_dict['Lib_name'], Mass_cold_dust_dict['Lib_name']]
    Umean_total_dust_dict['Col_numb'] = [Mass_warm_dust_dict['Col_numb'], Mass_cold_dust_dict['Col_numb']]
    Umean_total_dust_dict['Par_name'] = '$\\left<U\\right>$ (total)'
    Umean_total_dust_dict['range'] = [0.08,50.0]
    # 
    # analyze 
    print('Num_params', Num_params)
    print('Lib_params', Lib_params)
    Plot_engine = CrabPlot(figure_size=(12.0,10.0))
    Plot_engine.set_margin(panel=0, top=0.96, bottom=0.04)
    analyze_chisq_distribution(Stellar_mass_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(Lumin_AGN_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(Umin_warm_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(Umin_cold_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(LTIR_warm_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(LTIR_cold_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(Mass_warm_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(Mass_cold_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(LTIR_total_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(Mass_total_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(fPDR_total_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(Umean_total_dust_dict, Plot_engine = Plot_engine)
    Plot_engine.savepdf(Output_name+'.chisq.pdf')
    #Plot_engine.show()
    Plot_engine.close()
    print('Output to "%s"!'%(Output_name+'.chisq.pdf'))
    # 
    






