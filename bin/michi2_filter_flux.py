#!/usr/bin/env python
# 

import os, sys, re
import numpy as np
import astropy.io.ascii as asciitable


# fixing ascii table read write has no Reader Writer in astropy 7.0.0
def asciitable_write(table, output, **kwargs):
    #from distutils.version import LooseVersion
    #if astropy.__version__ > '7.0.0'
    try:
        asciitable.write(table, output=output, **kwargs)
    except TypeError:
        if 'Writer' in kwargs:
            if kwargs['Writer'] is asciitable.FixedWidthTwoLine:
                kwargs['format'] = 'fixed_width_two_line'
            elif kwargs['Writer'] is asciitable.Ipac:
                kwargs['format'] = 'ipac'
            del kwargs['Writer']
        asciitable.write(table, output=output, **kwargs)


# fixing ascii table read write has no Reader Writer in astropy 7.0.0
def asciitable_read(table, **kwargs):
    #from distutils.version import LooseVersion
    #if astropy.__version__ > '7.0.0'
    try:
        table = asciitable.read(table, **kwargs)
    except TypeError:
        if 'Reader' in kwargs:
            if kwargs['Reader'] is asciitable.NoHeader:
                kwargs['format'] = 'no_header'
            del kwargs['Reader']
        table = asciitable.read(table, **kwargs)
    return table




####################################
#               MAIN               #
####################################

data_table_file = ''
out_file = ''
thresh_optical = 0.0
thresh_infrared = 0.0
thresh_millimeter = 0.0
thresh_radio = 0.0
method_optical = 0
method_infrared = 0
method_millimeter = 0
method_radio = 0
low_SNR_options = ['keep', 'discard', 'replace']
iarg = 1
while iarg < len(sys.argv):
    argstr = sys.argv[iarg].lower().replace('--', '-')
    if re.match(r'^[-][a-z]+.*$', argstr):
        if argstr == '-optical' or argstr == '-opt':
            try:   
                thresh_optical = float(sys.argv[iarg+1])
                iarg += 1
            except:
                raise Exception('Error! Option -optical should be followed with a float number!')
            try:   
                method_optical = low_SNR_options.index(str(sys.argv[iarg+1]).lower())
                iarg += 1
            except:
                pass
        elif argstr == '-infrared' or argstr == '-ir':
            try:   
                thresh_infrared = float(sys.argv[iarg+1])
                iarg += 1
            except:
                raise Exception('Error! Option -infrared should be followed with a float number!')
            try:   
                method_infrared = low_SNR_options.index(str(sys.argv[iarg+1]).lower())
                iarg += 1
            except:
                pass
        elif argstr == '-millimeter' or argstr == '-mm':
            try:   
                thresh_millimeter = float(sys.argv[iarg+1])
                iarg += 1
            except:
                raise Exception('Error! Option -millimeter should be followed with a float number!')
            try:   
                method_millimeter = low_SNR_options.index(str(sys.argv[iarg+1]).lower())
                iarg += 1
            except:
                pass
        elif argstr == '-radio':
            try:   
                thresh_radio = float(sys.argv[iarg+1])
                iarg += 1
            except:
                raise Exception('Error! Option -radio should be followed with a float number!')
            try:   
                method_radio = low_SNR_options.index(str(sys.argv[iarg+1]).lower())
                iarg += 1
            except:
                pass
    else:
        if data_table_file == '':
            data_table_file = sys.argv[iarg]
        elif out_file == '':
            out_file = sys.argv[iarg]
    iarg += 1

if data_table_file == '' or out_file == '':
    print('Usage: michi2_filter_flux.py input_flux.txt output_flux.txt -options')
    print("Options:")
    print("  -optical <Nsigma=2> <method=discard>    : filter data points in optical to mid-infrared (<8micron)")
    print("  -infrared <Nsigma=3> <method=replace>   : filter data points in the far-infrared (=8-1000micorn)")
    print("  -millimeter <Nsigma=3> <method=replace> : filter data points in millimeter and radio (>1000micron)")
    print("  -radio <Nsigma=3> <method=replace>      : filter data points in millimeter and radio (>1000micron)")
    print("Notes: ")
    print("  Nsigma = <float>                        : the sigma threshold to filter data points by the method.")
    print("  method = keep|discard|replace           : keep, or discard, or replace lower S/N data points.")
    print("           In the case of replace, we replace S/N~1-2 data with f+1*e, and S/N~0-1 data with f+2*e")
    print("Example:")
    print("  michi2_filter_flux.py input_flux.txt output_flux.txt -optical 2 discard -infrared 3 replace -millimeter 3 replace -radio 2 discard")
    print("")
    sys.exit()

data_table = asciitable_read(data_table_file)

if not len(data_table.colnames) >= 3:
    print('Error! The input flux data table does not have at least three columns: wavelength, flux density and error in flux density.')
    sys.exit()

w = data_table.field(data_table.colnames[0])
f = data_table.field(data_table.colnames[1])
ferr = data_table.field(data_table.colnames[2])

print('Filtering optical with threshold {} sigma and method "{}"'.format(thresh_optical, low_SNR_options[method_optical]))
print('Filtering infrared with threshold {} sigma and method "{}"'.format(thresh_infrared, low_SNR_options[method_infrared]))
print('Filtering millimeter with threshold {} sigma and method "{}"'.format(thresh_millimeter, low_SNR_options[method_millimeter]))
print('Filtering radio with threshold {} sigma and method "{}"'.format(thresh_radio, low_SNR_options[method_radio]))

def process_flux_data(w, f, ferr, mask, method):
    if method == 1:
        rows = np.argwhere(mask).flatten()
        if len(rows) > 0:
            for i in rows:
                print('Discarding upper limit at row %d w %s f %s +- %s'%(i, w[i], f[i], ferr[i]))
                w[i] = -99.
                f[i] = -99.
                ferr[i] = 1e10
    elif method == 2:
        rows = np.argwhere(mask & (f>=1.0*ferr) & (f<2.0*ferr)).flatten()
        if len(rows) > 0:
            for i in rows:
                f[i] = f[i]+1.0*ferr[i]
                print('Fitting upper limit at row %d w %s as f+1*e +- e: %s +- %s'%(i, w[i], f[i], ferr[i]))
        rows = np.argwhere(mask & (f<1.0*ferr)).flatten()
        if len(rows) > 0:
            for i in rows:
                f[i] = f[i]+2.0*ferr[i]
                print('Fitting upper limit at row %d w %s as f+2*e +- e: %s +- %s'%(i, w[i], f[i], ferr[i]))
    # 
    # set zero error data error to 1/3 flux
    mask2 = (ferr==0)
    rows2 = np.argwhere(mask2).flatten()
    if len(rows2) > 0:
        ferr[mask2] = f[mask2] / 3.0
        for i in rows2:
            print('Limited row %d error from zero to 1/3: w = %s, f = %s, ferr = %s'%(i, w[i], f[i], ferr[i]))
    # 
    # 
    # limit S/N to be not larger than 10
    mask3 = (ferr<f/10.0)
    rows3 = np.argwhere(mask3).flatten()
    if len(rows3) > 0:
        ferr[mask3] = f[mask3] / 10.0
        for i in rows3:
            print('Limited row %d S/N no larger than 10: w = %s, f = %s, ferr = %s'%(i, w[i], f[i], ferr[i]))
    # 
    return w, f, ferr

# 
# filter data at optical
print('f', f)
print('ferr', ferr)
mask = (w<8.0) & (f<thresh_optical*ferr)
w, f, ferr = process_flux_data(w, f, ferr, mask, method_optical)

# 
# filter data at infrared
mask = (w>=8.0) & (w<1000.0) & (f<thresh_infrared*ferr)
w, f, ferr = process_flux_data(w, f, ferr, mask, method_infrared)

# 
# filter data at mm
mask = (w>=1000.0) & (w<10000.0) & (f<thresh_millimeter*ferr)
w, f, ferr = process_flux_data(w, f, ferr, mask, method_millimeter)

# 
# filter data at radio
mask = (w>=10000.0) & (f<thresh_radio*ferr)
w, f, ferr = process_flux_data(w, f, ferr, mask, method_radio)


# 
# deal with duplicated w
i = 0
while i < len(data_table):
    # identify duplicated w
    mask2 = (w==w[i]) & (w>0.0)
    isel2 = np.argwhere(mask2).flatten()
    if len(isel2) >= 2:
        # found duplicated w
        print('Found wavelength-duplicated rows: %s'%(isel2))
        print(data_table[isel2])
        f_to_average = f[mask2]
        ferr_to_average = ferr[mask2]
        f_averaged = np.sum(f_to_average/ferr_to_average**2)/np.sum(1/ferr_to_average**2)
        ferr_averaged = np.sqrt(1/np.sum(ferr_to_average**(-2))) # error propagation of weighted mean, see -- http://www.physics.umd.edu/courses/Phys261/F06/ErrorPropagation.pdf
        # limit S/N not larger than 10
        #if ferr_averaged < f_averaged/10.0:
        #    ferr_averaged = f_averaged/10.0
        # store into data_table
        f[i] = f_averaged # change f will directly change data_table!
        ferr[i] = ferr_averaged # change ferr will directly change data_table!
        print('Averaged wavelength-duplicated rows: w = %s, f = %s, ferr = %s'%(w[i], f_averaged, ferr_averaged))
        # remove those duplicated rows, but keep current i row.
        isel3 = isel2[(isel2 != i)]
        for iseli in isel3:
            print('data_table.remove_rows(%d)'%(iseli))
        data_table.remove_rows(isel3)
    i = i+1


# 
# exclude invalid wavelength and invalid flux error data
mask = (w<=0) | (ferr>=1e10) #<20180814> fixed bug, added "| (ferr>=1e10)"
rows = np.argwhere(~mask).flatten()
w = w[rows]
f = f[rows]
ferr = ferr[rows]


# 
# sort
data_table.sort(data_table.colnames[0])


# 
# output
asciitable_write(data_table, out_file, Writer=asciitable.Ipac, delimiter='    ', overwrite=True)
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







