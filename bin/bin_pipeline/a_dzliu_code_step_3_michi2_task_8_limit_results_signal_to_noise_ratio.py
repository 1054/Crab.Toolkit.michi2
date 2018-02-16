#!/usr/bin/env python
# 

import os
import sys
import numpy
import astropy
import astropy.io.ascii as asciitable
from copy import copy
import json



if not os.path.isdir('Results'):
    print('The "Results" folder is not found! Maybe you should run the *.sh script instead of this *.py script?')
    sys.exit()




####################################
#               MAIN               #
####################################

for param in ['LAGN', 'LIR_total', 'Mdust_total', 'Mstar', 'Umean_total', 'fPDR_total']:
    data_file = 'Results/best-fit_param_%s.txt'%(param)
    data_table = asciitable.read(data_file)
    # append a column
    if not ('%s_width_PDF'%(param) in data_table.colnames):
        data_table['%s_width_PDF'%(param)] = copy(data_table['%s_sigma'%(param)])
    else:
        print('Found column %s_width_PDF'%(param))
    data_table['SNR_100_3000_um'] = copy(data_table['%s_width_PDF'%(param)]) * 0.0
    print(data_table.colnames)
    #continue
    # 
    # limit param_sigma to be no lower than 1/snr['SNR_100_3000_um']
    if len(data_table['Source']) > 0:
        for i in range(len(data_table['Source'])):
            source_name = data_table['Source'][i]
            with open('%s/SNR.json'%(source_name),'r') as fp:
                snr = json.load(fp)
                fp.close()
            if 'SNR_100_3000_um' in snr:
                data_table['SNR_100_3000_um'][i] = snr['SNR_100_3000_um']
                if data_table['%s_width_PDF'%(param)][i] < 1.0/snr['SNR_100_3000_um']: 
                    data_table['%s_sigma'%(param)][i] = 1.0/snr['SNR_100_3000_um'] 
                    data_table['%s_L68'%(param)][i] = data_table['%s_median'%(param)][i] - data_table['%s_width_PDF'%(param)][i]
                    data_table['%s_H68'%(param)][i] = data_table['%s_median'%(param)][i] + data_table['%s_width_PDF'%(param)][i]
                else: 
                    data_table['%s_sigma'%(param)][i] = data_table['%s_width_PDF'%(param)][i]
                #data_table['%s_sigma'%(param)][i] = 1.0/snr['SNR_100_3000_um'] if data_table['%s_width_PDF'%(param)][i] < 1.0/snr['SNR_100_3000_um'] else data_table['%s_width_PDF'%(param)][i]
            else:
                print('Error! Failed to read SNR from "%s/SNR.json"!'%(source_name))
                sys.exit()
    else:
        print('Error! Failed to read any content from "%s"!'%(data_file))
        sys.exit()
    #print(data_table)
    
    # backup
    os.system('cp "Results/best-fit_param_%s.txt" "Results/best-fit_param_%s.txt.backup"'%(param, param))
    
    # output
    out_file = data_file
    asciitable.write(data_table, out_file, Writer=asciitable.Ipac, delimiter='    ', overwrite=True)
    #asciitable.write(data_table, sys.stdout, Writer=asciitable.Ipac, delimiter='  ')
    with open(out_file, 'r+') as fp:
        out_content = fp.readlines() # read everything in the file
        out_iline = 0
        out_header = [] # Ipac format has multiple comment lines (commented by the char '\\') and 4 header lines.
        fp.seek(0)
        while out_iline < len(out_content):
            if out_content[out_iline][0] == '\\':
                # if his is a commented line, then we change the comment mark to '#'
                out_content[out_iline] = '#' + out_content[out_iline][1:]
                fp.write(out_content[out_iline])
            else:
                if len(out_header) == 0:
                    # if this is the first header line, then replace the first white space by '#', or if there is no white space, preprend '#'.
                    if out_content[out_iline][0] == ' ':
                        out_content[out_iline] = '#' + out_content[out_iline][1:]
                    else:
                        out_content[out_iline] = '#' + out_content[out_iline]
                    # append header to 'out_header' list
                    out_header.append(out_content[out_iline])
                    # write only one header line
                    fp.write(out_content[out_iline])
                    # 
                elif len(out_header) < 4:
                    # append header to 'out_header' list
                    out_header.append(out_content[out_iline])
                    # skip the 2nd to 4th header line
                    pass
                else:
                    # write data line
                    fp.write(out_content[out_iline])
                    # 
            out_iline = out_iline + 1
        fp.truncate()
        fp.close()
    #os.system('sed -i.bak -e "$(grep \"\\\" %s | wc -l)s/^ /#/" "%s"'%(out_file, out_file))
    #os.system('sed -i.bak -e "2d;3d;4d" "%s"'%(out_file))
    #if os.path.isfile(out_file+'.bak'):
    #    os.system('rm "%s"'%(out_file+'.bak'))
    print('Output to "%s"!'%(out_file))




