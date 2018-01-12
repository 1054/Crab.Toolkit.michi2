#!/usr/bin/env python3.6
# 

import os
import sys

#sys.path.append(os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabtable')
#sys.path.append(os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabplot')
sys.path.insert(1, '/Users/dzliu/Softwares/Python/lib/crab/crabtable')
sys.path.insert(1, '/Users/dzliu/Softwares/Python/lib/crab/crabplot')

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
#               Functions               #
#########################################

def analyze_chisq_distribution(param_dict, verbose=1, Plot_engine=None):
    if 'Lib_file' in param_dict and \
        'Lib_name' in param_dict and \
        'Lib_numb' in param_dict and \
        'Par_name' in param_dict and \
        'Col_numb' in param_dict and \
        'value' in param_dict and \
        'chisq' in param_dict :
        if verbose>=1:
            print('Analyzing the chi-square distribution for parameter "%s" in library %s from file "%s"'%(param_dict['Par_name'], param_dict['Lib_name'], param_dict['Lib_file']))
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
        param_array = param_dict['value']
        chisq_array = param_dict['chisq']
        # 
        # convert to physical values
        #if param_dict['Par_name']
        # 
        # bin param
        param_log = False
        if 'Log' in param_dict:
            if param_dict['Log'] == True:
                param_log = True
        if param_log == True: 
            param_log_mask = (param_array>0)
            param_log_mask2 = (param_array<=0)
            param_array[param_log_mask] = numpy.log10(param_array[param_log_mask])
            param_array[param_log_mask2] = numpy.nan
        param_min = numpy.nanmin(param_array)
        param_max = numpy.nanmax(param_array)
        param_bin_numb = 15 #<TODO># 
        param_bin_step = (param_max-param_min)/param_bin_numb
        # 
        # compute minimum chisq in each parameter bin
        param_bin_x = []
        param_bin_y = []
        for i in range(param_bin_numb):
            if i == param_bin_numb-1:
                if verbose>=2:
                    print('Binning from param value %s to %s, step %s'%(param_min+i*param_bin_step, param_max, param_bin_step))
                param_bin_mask = (param_array>=param_min+i*param_bin_step) & (param_array<=param_max)
            else:
                if verbose>=2:
                    print('Binning from param value %s to %s, step %s'%(param_min+i*param_bin_step, param_min+(i+1)*param_bin_step, param_bin_step))
                param_bin_mask = (param_array>=param_min+i*param_bin_step) & (param_array<param_min+(i+1)*param_bin_step)
            # 
            param_bin_args = numpy.argwhere(param_bin_mask)
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
        # 
        # Plot histogram
        Plot_engine.plot_hist(param_bin_x, param_bin_y, width = 0.90 * param_bin_step, align = 'edge', overplot = False, 
                                xtitle = param_dict['Par_name'], ytitle = '$\chi^2$', useTex = True)
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
    SelectNumber = 10 #<TODO># how many SEDs to show?
    SelectIndex = SortedIndex[0:(len(SortedIndex) if len(SortedIndex)<SelectNumber else SelectNumber)]
    #print(DataTable.TableData[SelectIndex])
    Arr_chi2 = DataArray['chi2'][SortedIndex]
    Min_chi2 = numpy.nanmin(Arr_chi2[0:SelectNumber-1])
    Max_chi2 = numpy.nanmax(Arr_chi2[0:SelectNumber-1])
    Cut_chi2 = Min_chi2 + 5.89 # 5 parameter, 68% confidence
    SelectNumber = len(numpy.argwhere(DataArray['chi2'][SortedIndex]<=Cut_chi2))
    SelectIndex = SortedIndex[0:(len(SortedIndex) if len(SortedIndex)<SelectNumber else SelectNumber)]
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
        for i in range(len(SelectIndex)):
            for j in range(int(InfoDict['NLIB'])):
                if not os.path.isdir('obj_%d'%(i+1)):
                    os.mkdir('obj_%d'%(i+1))
                # 
                if not os.path.isfile('obj_%d/SED_LIB%d'%(i+1,j+1)):
                    BashCommand = 'cd obj_%d/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_SED ../%s %d %s SED_LIB%d'%\
                                        (i+1, \
                                            InfoDict['LIB%d'%(j+1)], \
                                                DataArray['i%d'%(j+1)][SelectIndex[i]], \
                                                    DataArray['a%d'%(j+1)][SelectIndex[i]], \
                                                        j+1)
                    print(BashCommand)
                    os.system(BashCommand)
        # 
        # Wait for a long time
        # 
        # Then plot SEDs
        Redshift = float(InfoDict['REDSHIFT'])
        Color_list = ['cyan', 'gold', 'red', 'blue', 'purple']
        Plot_engine = CrabPlot(figure_size=(9.0,5.0))
        for i in range(len(SelectIndex)-1,-1,-1):
            # 
            # alpha by chi2
            #print(Arr_chi2[i])
            Alpha_chi2 = Plot_engine.get_color_by_value([Min_chi2-(Max_chi2-Min_chi2)*0.2, Max_chi2+(Max_chi2-Min_chi2)*0.3], 
                                                        input_value=Arr_chi2[i], 
                                                        log=1, 
                                                        cmap=matplotlib.cm.get_cmap('gray_r'))[0]
            #print('Alpha_chi2: ', Alpha_chi2)
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
            obj_SED_1 = Plot_engine.Plot_data['obj_%d_SED_LIB1'%(i+1)]
            obj_SED_2 = Plot_engine.Plot_data['obj_%d_SED_LIB2'%(i+1)]
            obj_SED_3 = Plot_engine.Plot_data['obj_%d_SED_LIB3'%(i+1)]
            obj_SED_4 = Plot_engine.Plot_data['obj_%d_SED_LIB4'%(i+1)]
            obj_SED_5 = Plot_engine.Plot_data['obj_%d_SED_LIB5'%(i+1)]
            obj_SED_sum_x_lg = numpy.arange(-2,6,0.001) # wavelength_um grid
            obj_SED_sum_x = numpy.power(10, obj_SED_sum_x_lg) # make it in linear space
            obj_SED_sum_y = []
            for j in range(int(InfoDict['NLIB'])):
                obj_SED_single_x = Plot_engine.Plot_data['obj_%d_SED_LIB%d'%(i+1,j+1)][:,0]
                obj_SED_single_y = Plot_engine.Plot_data['obj_%d_SED_LIB%d'%(i+1,j+1)][:,1]
                obj_SED_spline_y = Plot_engine.spline(obj_SED_single_x, obj_SED_single_y, obj_SED_sum_x_lg, xlog=1, ylog=1, outputxlog=0)
                if j == 0:
                    obj_SED_sum_y = obj_SED_spline_y
                else:
                    obj_SED_sum_y = obj_SED_sum_y + obj_SED_spline_y
            # 
            #pprint(obj_SED_sum_y)
            #print(numpy.column_stack((obj_SED_sum_x, obj_SED_sum_y)))
            #print(Arr_chi2[i])
            Color_chi2 = Plot_engine.get_color_by_value([Min_chi2, Max_chi2+(Max_chi2-Min_chi2)*0.3], 
                                                        input_value=Arr_chi2[i], 
                                                        log=1, 
                                                        cmap=matplotlib.cm.get_cmap('gray'))
            #print('Color_chi2: ', Color_chi2)
            Plot_engine.plot_line(obj_SED_sum_x, obj_SED_sum_y, current=1, color=Color_chi2)
            # 
            # show chi2 on the figure
            if i == len(SelectIndex)-1:
                Plot_engine.xyouts(0.05, 0.95, '$\chi^2:$', NormalizedCoordinate=True, useTex=True)
                Plot_engine.xyouts(0.15, 0.95, '$z=%s$'%(Redshift), NormalizedCoordinate=True, useTex=True)
            Plot_engine.xyouts(0.09, 0.95-0.03*i, '%.1f'%(Arr_chi2[i]), NormalizedCoordinate=True, useTex=True, color=Color_chi2)
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
        Plot_engine.set_xtitle('Wavelength [um]')
        Plot_engine.set_ytitle('Flux Density [mJy]')
        Plot_engine.show()
        Plot_engine.savepdf(Output_name+'.pdf')
        Plot_engine.close()
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
    LIR_warm_dust_dict = {}
    LIR_cold_dust_dict = {}
    Umin_warm_dust_dict = {}
    Umin_cold_dust_dict = {}
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
                    Stellar_mass_dict['Lib_file'] = InfoDict[Lib_name]
                    Stellar_mass_dict['Lib_name'] = Lib_name
                    Stellar_mass_dict['Lib_numb'] = j+1
                    Stellar_mass_dict['Par_name'] = Lib_dict[Key_TPAR]
                    Stellar_mass_dict['Col_numb'] = Col_number
                    Stellar_mass_dict['Log'] = True
                    Stellar_mass_dict['value'] = DataArray['a%d'%(j+1)] * DataTable.getColumn('col%d'%(Col_number))
                    Stellar_mass_dict['chisq'] = DataArray['chi2']
            elif InfoDict[Lib_name].find('DL07.HiExCom') >= 0:
                if 'lgLTIR' == Lib_dict[Key_TPAR]:
                    LIR_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LIR_warm_dust_dict['Lib_name'] = Lib_name
                    LIR_warm_dust_dict['Lib_numb'] = j+1
                    LIR_warm_dust_dict['Par_name'] = Lib_dict[Key_TPAR]
                    LIR_warm_dust_dict['Col_numb'] = Col_number
                    LIR_warm_dust_dict['Log'] = True
                    LIR_warm_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn('col%d'%(Col_number))) * 4*pi*dL**2 # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LIR_warm_dust_dict['chisq'] = DataArray['chi2']
                elif 'Umin' == Lib_dict[Key_TPAR]:
                    Umin_warm_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Umin_warm_dust_dict['Lib_name'] = Lib_name
                    Umin_warm_dust_dict['Lib_numb'] = j+1
                    Umin_warm_dust_dict['Par_name'] = Lib_dict[Key_TPAR]
                    Umin_warm_dust_dict['Col_numb'] = Col_number
                    Umin_warm_dust_dict['Log'] = True
                    Umin_warm_dust_dict['value'] = DataTable.getColumn('col%d'%(Col_number))
                    Umin_warm_dust_dict['chisq'] = DataArray['chi2']
            elif InfoDict[Lib_name].find('DL07.LoExCom') >= 0:
                if 'lgLTIR' == Lib_dict[Key_TPAR]:
                    LIR_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    LIR_cold_dust_dict['Lib_name'] = Lib_name
                    LIR_cold_dust_dict['Lib_numb'] = j+1
                    LIR_cold_dust_dict['Par_name'] = Lib_dict[Key_TPAR]
                    LIR_cold_dust_dict['Col_numb'] = Col_number
                    LIR_cold_dust_dict['Log'] = True
                    LIR_cold_dust_dict['value'] = DataArray['a%d'%(j+1)] * numpy.power(10,DataTable.getColumn('col%d'%(Col_number))) * 4*pi*dL**2 # Note that we need to carefully convert lgLTIR from log space to LIR in linear space, and apply the normalization.
                    LIR_cold_dust_dict['chisq'] = DataArray['chi2']
                elif 'Umin' == Lib_dict[Key_TPAR]:
                    Umin_cold_dust_dict['Lib_file'] = InfoDict[Lib_name]
                    Umin_cold_dust_dict['Lib_name'] = Lib_name
                    Umin_cold_dust_dict['Lib_numb'] = j+1
                    Umin_cold_dust_dict['Par_name'] = Lib_dict[Key_TPAR]
                    Umin_cold_dust_dict['Col_numb'] = Col_number
                    Umin_cold_dust_dict['Log'] = True
                    Umin_cold_dust_dict['value'] = DataTable.getColumn('col%d'%(Col_number))
                    Umin_cold_dust_dict['chisq'] = DataArray['chi2']
            # 
            # count the column number in the big chi-square data table. All params in each LIB are listed in the big chi-square data table.
            Col_number = Col_number + 1
    # 
    # Total LIR
    LIR_total_dust_dict = copy(LIR_warm_dust_dict)
    LIR_total_dust_dict['value'] = LIR_warm_dust_dict['value'] + LIR_cold_dust_dict['value']
    # 
    # analyze 
    print('Num_params', Num_params)
    print('Lib_params', Lib_params)
    Plot_engine = CrabPlot(figure_size=(9.0,7.0))
    analyze_chisq_distribution(Stellar_mass_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(LIR_warm_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(LIR_cold_dust_dict, Plot_engine = Plot_engine)
    analyze_chisq_distribution(LIR_total_dust_dict, Plot_engine = Plot_engine)
    Plot_engine.savepdf(Output_name+'.chisq.pdf')
    Plot_engine.show()
    # 
    






