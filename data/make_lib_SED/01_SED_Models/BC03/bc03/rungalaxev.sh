#!/bin/csh
#trivial script to run the python script as the env method is not working for me
# expects to find the config file as the first parameter
# note that this needs to be csh for the python part to work for some unknown reason.

set dirname=`dirname $0`
#need to set up the environment.
setenv bc03 $dirname/src/
echo $bc03
source $bc03/.bc_cshrc

python $dirname/runBC03reduced.py $1

if ( -f $dirname/galevot/Galevot.jar ) then
 mv testout_ages testout_ages.txt; java -jar $dirname/galevot/Galevot.jar testout_ages.txt testout_ages;
else
 echo -e "the galevot converter is not present";exit 1;
endif
