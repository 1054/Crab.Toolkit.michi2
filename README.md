# Crab.Toolkit.michi2
This is a _Crab_ toolkit for minimumizing chi2 for any 1D model fit, for example, galaxy SED fitting, molecular line ladder fitting. 

## Usage ##

### SED Fitting ###

#### Get the code (with git)
```
cd /some/path/
git clone https://github.com/1054/Crab.Toolkit.michi2
```

#### Source the code (must under BASH shell)
```
bash
source /some/path/Crab.Toolkit.michi2/SETUP.bash
```

#### Prepare photometry data
Assuming we have a photometric catalog, where columns are different photometric bands and rows are different sources. Then we need to prepare one SED fitting input file for each source you would like to fit. Each SED fitting input file should be a text file and has three columns: first column the wavelength in micron-meter unit, second column the flux density in milli-Jansky unit, and the third column the error in flux density in milli-Jansky unit. It should have multiple rows corresponding to each photometric band, and better S/N>3 bands. Columns should be separated by white space. Rows can be commented by # character. 

An example of the SED fitting input file is like: 
```
# wavelength          flux        flux_err
#        um            mJy             mJy
    3.56343     0.12417548     0.012417548
    4.51101     0.09488227     0.009488227
    5.75934    0.075890755    0.0075890755
    7.95949     0.08685837     0.009947749
       24.0      1.3679016      0.13679016
      100.0      41.577702       4.1577702
      160.0        60.0042       6.3393602
      250.0        58.0821         5.80821
      350.0      36.100101       6.0408401
      500.0      18.878901       2.2385001
      850.0       3.866907        1.258796
```

#### Run michi2 ####
If you have already `sourced` the `SETUP.bash`, then just change directory to where you store your SED fitting input file (assuming it's named "extracted_flux.txt"), and run michi2.

```
cd /path/to/your/data/directory/

ls "extracted_flux.txt"

michi2-run-fitting-5-components-applying-evolving-qIR # call it without any argument will print the usage

michi2-run-fitting-5-components-applying-evolving-qIR -redshift 1.5 -flux "extracted_flux.txt" -parallel 2

# Note that for this example we set redshift to 1.5, and fit with 2 CPU cores. 
```

The michi2 SED fitting is currently VERY SLOW. It can take three hours on laptops! It is because it stupidly loops all the models. We will adopt Markov chain Monte Carlo method in the future... For our own, we use it on 100plus-CPU-core machine, so the current code is fine. 

The output of michi2 SED fitting will be: 
```
ls fit_5.out
ls fit_5.out.info
```

Then we make SED plots and chi-square plots. 

#### Plot chi2 distribution ####
Here we make the SED and chi-square plots, assuming that you have already `sourced` the `SETUP.bash`. 
```
michi2-plot-fitting-results # call it without any argument will print the usage

michi2-plot-fitting-results fit_5.out -flux extracted_flux.txt -source YOUR_SOURCE_NAME
```

Then the output files will be:
```
ls fit_5.pdf
ls fit_5.chisq.pdf
ls best-fit*.txt
```



### LVG Fitting ###

#### Prepare CO ladder, run michi2, make plots ####
TODO
```
cat co_ladder.txt
source /some/path/Crab.Toolkit.michi2/SETUP.bash
michi2-deploy-files LVG
michi2_v04 -obs co_ladder.txt -lib lib_z_1.500_co_dvddr_5.0_dVDoW_50.0_Faster.lvg lib_z_1.500_co_dvddr_5.0_dVDoW_50.0_Faster.lvg -out fit_2_lvg.out -constraint LIB1 PAR1 LT LIB2 PAR1
./run_plotting_bestfit_LVG_two_components.sh
```


