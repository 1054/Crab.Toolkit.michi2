#!/bin/bash
#
# readlink
if [[ $(uname) == *"Darwin"* ]]; then
    function readlink() {
        if [[ $# -gt 1 ]]; then if [[ "$1" == "-f" ]]; then shift; fi; fi
        DIR="$1"; if [[ "$DIR" != *"/"* ]]; then DIR="./$DIR"; fi # 20170228: fixed bug: path without "/"
        DIR=$(echo "${DIR%/*}") # 20160410: fixed bug: source SETUP just under the Softwares dir
        if [[ -d "$DIR" ]]; then cd "$DIR" && echo "$(pwd -P)/$(basename ${1})"; 
        else echo "$(pwd -P)/$(basename ${1})"; fi
    }
fi
CRABTOOLKITMICHI2=$(dirname $(readlink -f "${BASH_SOURCE[0]}"))
export CRABTOOLKITMICHI2
#
# PATH
if [[ x"$PATH" != x*"$CRABTOOLKITMICHI2/bin"* ]]; then
    if [[ ! -z "IDL_PATH" ]]; then
        export PATH="$CRABTOOLKITMICHI2/bin:$PATH"
    else
        export PATH="$CRABTOOLKITMICHI2/bin:/bin:/usr/bin:/usr/local/bin:/opt/local/bin"
    fi
fi
#
# IDL_PATH
if [[ x"$IDL_PATH" != x*"+$CRABTOOLKITMICHI2/lib/idl"* ]]; then
    if [[ ! -z "IDL_PATH" ]]; then
        export IDL_PATH="+$CRABTOOLKITMICHI2/lib/idl:$IDL_PATH"
    else
        export IDL_PATH="+$CRABTOOLKITMICHI2/lib/idl"
    fi
fi
#
# LIST
CRABTOOLKITCMD=("michi2-deploy-files" "michi2-run-SED-fitting-v5")
# 
# CHECK
# -- 20160427 only for interactive shell
# -- http://stackoverflow.com/questions/12440287/scp-doesnt-work-when-echo-in-bashrc
if [[ $- =~ "i" ]]; then 
  for TEMPTOOLKITCMD in ${CRABTOOLKITCMD[@]}; do
    type $TEMPTOOLKITCMD
  done
fi


