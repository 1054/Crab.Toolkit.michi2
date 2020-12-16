#!/bin/bash
# 
# common caller for wcstools:
#=  bin_portable.sh
#=  michi2_v04
#=  michi2_read_lib_SED
#=  michi2_compute_delta_chisq_1sigma
#=  lumdist
# 
# for ff in $(cat bin_portable.sh | grep "^#= " | tr -s ' ' | cut -d ' ' -f 2); do cp bin_portable.sh $ff; done
#   
# 
root_dirname=$(dirname "${BASH_SOURCE[0]}")
root_basename=$(basename "${BASH_SOURCE[0]}")
# Linux
if [[ $(uname -s) == Linux ]]; then
    ldd_versions=($(ldd --version | head -n 1 | perl -p -e 's/.* ([0-9.]+) *$/\1/g' | perl -p -e 's/\./ /g'))
    bin_dirname="bin_linux_glibc_2_22"
    if [[ ${#ldd_versions[@]} -ge 2 ]]; then
        ldd_version_major="${ldd_versions[0]}"
        ldd_version_minor="${ldd_versions[1]}"
        if [[ $(bc <<< "${ldd_version_major}==2") -eq 1 ]]; then
            if [[ $(bc <<< "${ldd_version_minor}<=12") -eq 1 ]]; then
                bin_dirname="bin_linux_glibc_2_12" # For computer with old GLIBC version
            elif [[ $(bc <<< "${ldd_version_minor}<=14") -eq 1 ]]; then
                bin_dirname="bin_linux_glibc_2_14" # For computer with old GLIBC version
            elif [[ $(bc <<< "${ldd_version_minor}<=22") -eq 1 ]]; then
                bin_dirname="bin_linux_glibc_2_22" # For computer with old GLIBC version
            fi
        fi
    fi
    ${root_dirname}/${bin_dirname}/${root_basename}_linux_x86_64 "$@"
    #
    #if [[ $(bc <<< "$ldd_version_number<=12") -eq 1 ]]; then
    #    # the supercomputer planer has an old GLIBC version 2.5
    #    $(dirname $0)/bin_linux_glibc_2_12/$(basename $0)_linux_$(arch) $@
    #elif [[ $(bc <<< "$ldd_version_number<=14") -eq 1 ]]; then
    #    $(dirname $0)/bin_linux_glibc_2_14/$(basename $0)_linux_$(arch) $@
    #elif [[ $(bc <<< "$ldd_version_number<=22") -eq 1 ]]; then
    #    # isaac, aida
    #    $(dirname $0)/bin_linux_glibc_2_22/$(basename $0)_linux_$(arch) $@
    #else
    #    $(dirname $0)/bin_linux_glibc_2_22/$(basename $0)_linux_$(arch) $@
    #fi
    #${root_dirname}/${bin_dirname}/${root_basename}_linux_x86_64 "$@"
fi
# Darwin
if [[ $(uname -s) == Darwin ]]; then
    $(dirname $0)/bin_mac/$(basename $0)_mac $*
fi
# Cygwin
if [[ $(uname -s) == *CYGWIN* ]]; then
    $(dirname $0)/bin_cygwin/$(basename $0)_cygwin_x86.exe $@
fi


