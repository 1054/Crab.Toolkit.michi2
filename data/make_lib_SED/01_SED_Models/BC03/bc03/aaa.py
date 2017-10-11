#!/usr/bin/env python2.7
# 

# Run BC2003 models with a reduced number of options.
# call.: python runBC03reduced.py inputfile
# some hacks from pah@jb.man.ac.uk
# last edited by dzliu



import os
import sys
import string
import shutil


# Check input argument

#if(len(sys.argv)<=1):
#    print("Usage: input argument 1 is the inputfile, see \"inputfile.default\" as an example.")
#    sys.exit()


# Get input file
#inputfile = sys.argv[1]
inputfile = os.path.dirname(sys.argv[0])+os.sep+'inputfile.default'


# Read input file
input0 = []
input1 = []
file = open(inputfile,"r")
while 1:
    line = file.readline()
    if not line:
        break
    else:
        cols = string.split(line)
        input0.append(cols[0])
        input1.append(cols[1])
file.close()

inputs = {}
for i in range(len(input0)):
    inputs[input0[i]]=input1[i]

#quick hack
inputs["ROOTDIR"]=os.path.dirname(sys.argv[0])+ os.sep
inputs["CSPOUTFILE"]="output_csp_file"



if inputs["TRACKS"] == "Padova2000":
    if float(inputs["METALLICITY"]) == 0.0004:
        inputs["mcode"] = "m122"
    elif float(inputs["METALLICITY"]) == 0.001:
        inputs["mcode"] = "m132"
    elif float(inputs["METALLICITY"]) == 0.004:
        inputs["mcode"] = "m142"
    elif float(inputs["METALLICITY"]) == 0.008:
        inputs["mcode"] = "m152"
    elif float(inputs["METALLICITY"]) == 0.019:
        inputs["mcode"] = "m162"
    elif float(inputs["METALLICITY"]) == 0.03:
        inputs["mcode"] = "m172"
elif inputs["TRACKS"] == "Padova1994":
    if float(inputs["METALLICITY"]) == 0.0001:
        inputs["mcode"] = "m22"
    elif float(inputs["METALLICITY"]) == 0.0004:
        inputs["mcode"] = "m32"
    elif float(inputs["METALLICITY"]) == 0.004:
        inputs["mcode"] = "m42"
    elif float(inputs["METALLICITY"]) == 0.008:
        inputs["mcode"] = "m52"
    elif float(inputs["METALLICITY"]) == 0.02:
        inputs["mcode"] = "m62"
    elif float(inputs["METALLICITY"]) == 0.05:
        inputs["mcode"] = "m72"

if inputs["TRACKS"] == 'Padova1994':
    if inputs["IMF"] == 'salpeter':
        shutil.copy(inputs["ROOTDIR"]+"models/Padova1994/salpeter/bc2003_"+inputs["RESOLUTION"]+"_"+inputs["mcode"]+"_salp_ssp.ised_ASCII.gz","TMPinputSSP.ised_ASCII.gz")
    elif inputs["IMF"] == 'chabrier':
        shutil.copy(inputs["ROOTDIR"]+"models/Padova1994/chabrier/bc2003_"+inputs["RESOLUTION"]+"_"+inputs["mcode"]+"_chab_ssp.ised_ASCII.gz","TMPinputSSP.ised_ASCII.gz")
elif inputs["TRACKS"] == 'Padova2000':
    if inputs["IMF"] == 'salpeter':
        shutil.copy(inputs["ROOTDIR"]+"models/Padova2000/salpeter/bc2003_"+inputs["RESOLUTION"]+"_"+inputs["mcode"]+"_salp_ssp.ised_ASCII.gz","TMPinputSSP.ised_ASCII.gz")
    elif inputs["IMF"] == 'chabrier':
        shutil.copy(inputs["ROOTDIR"]+"models/Padova2000/chabrier/bc2003_"+inputs["RESOLUTION"]+"_"+inputs["mcode"]+"_chab_ssp.ised_ASCII.gz","TMPinputSSP.ised_ASCII.gz")

os.system("gunzip TMPinputSSP.ised_ASCII.gz")
os.system(inputs["ROOTDIR"]+"src/bin_ised TMPinputSSP.ised_ASCII")
inputs["SSPmodel"] = "TMPinputSSP.ised"

if inputs["SFH"] == "exponential":
    inputs["SFHcode"] = "1"
elif inputs["SFH"] == "SSP":
    inputs["SFHcode"] = "0"
elif inputs["SFH"] == "singleburst":
    inputs["SFHcode"] = "2"

file = open("TMPinput.cspgal","w")
file.write( inputs["SSPmodel"] +"\n")
file.write(inputs["DUST"]+"\n")
if inputs["DUST"] == "Y":
    file.write(inputs["TAU_V"] +"\n")
    file.write(inputs["MU"] +"\n")

if inputs["SFH"] == "exponential":
    file.write(inputs["SFHcode"] +"\n")
    file.write(inputs["TAU"]  +"\n")
    file.write(inputs["GASRECYCLE"]  +"\n")
    if inputs["GASRECYCLE"] == "Y":
        file.write(inputs["EPSILON"]  +"\n")
    file.write(inputs["TCUTSFR"]  +"\n")
    file.write(inputs["CSPOUTFILE"]  +"\n")
    file.close()
elif inputs["SFH"] == "SSP":
    file.write(inputs["SFHcode"] +"\n")
    file.write(inputs["CSPOUTFILE"]  +"\n")
    file.close()
elif inputs["SFH"] == "singleburst":
    file.write(inputs["SFHcode"] +"\n")
    file.write(inputs["TAU"]  +"\n")
    file.write(inputs["CSPOUTFILE"]  +"\n")
    file.close()
os.system(inputs["ROOTDIR"]+"src/csp_galaxev < TMPinput.cspgal")
os.remove("TMPinput.cspgal")
os.remove("TMPinputSSP.ised_ASCII")
os.remove("TMPinputSSP.ised")

#Extract a spectrum at a given age
if inputs["FNU_LAMBDA"] == "nu":
    wlrangestring = "-"+inputs["SPECRANGEMIN"] +","+inputs["SPECRANGEMAX"]
elif  inputs["FNU_LAMBDA"] == "lambda":
    wlrangestring = inputs["SPECRANGEMIN"] +","+inputs["SPECRANGEMAX"]
if inputs["F0SCALING"] != "N":
    wlrangestring = wlrangestring + "," + inputs["W0SCALING"] +"," + inputs["F0SCALING"]
inputs["extract_spec_name"]  = inputs["CSPOUTFILE"]+"_ages" #+inputs["EXTRACT_AGE"]
file = open("TMPinput.gpl","w")
file.write(inputs["CSPOUTFILE"] +"\n")
file.write(wlrangestring +"\n")
file.write(inputs["EXTRACT_AGE"]   +"\n")
file.write(inputs["extract_spec_name"]   +"\n")
file.close()
os.system(inputs["ROOTDIR"]+"src/galaxevpl < TMPinput.gpl")
os.remove("TMPinput.gpl")

#Redshift evolution
file = open("TMPinput.cmevol","w")
file.write(inputs["H_0"] + ","+inputs["Omega"]+","+inputs["Omega_Lambda"]  +"\n")
file.write(inputs["GALAXY_AGE_TODAY"]   +"\n")
file.write(inputs["CSPOUTFILE"] +"\n")
if inputs["REDSHIFT_EVOL"] == "mag":
    file.write(inputs["FILTER1"]   +"\n")
elif inputs["REDSHIFT_EVOL"] == "col":
    file.write(inputs["FILTER1"] +","+inputs["FILTER2"] +   "\n")
file.close()
os.system(inputs["ROOTDIR"]+"src/cm_evolution < TMPinput.cmevol")
os.remove("TMPinput.cmevol")

#do some file renaming to make it easier to pick up the files
os.rename(inputs["CSPOUTFILE"]+".magnitude_F"+"%03d" % int(inputs["FILTER1"]), inputs["CSPOUTFILE"]+".magnitude_FILTER1")
if inputs["REDSHIFT_EVOL"] == "col":
    os.rename(inputs["CSPOUTFILE"]+".magnitude_F"+"%03d" % int(inputs["FILTER2"]), inputs["CSPOUTFILE"]+".magnitude_FILTER2")







