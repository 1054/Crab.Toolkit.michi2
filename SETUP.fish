#!/usr/bin/fish
# 
# 
# CRABTOOLKITMICHI2
if contains "Linux" (uname)
    set -x CRABTOOLKITMICHI2 (dirname (readlink -f (status --current-filename)))
end
if contains "Darwin" (uname)
    set -x CRABTOOLKITMICHI2 (dirname (perl -MCwd -e 'print Cwd::abs_path shift' (status --current-filename)))
end
export CRABTOOLKITMICHI2
#<DEBUG># echo "$CRABTOOLKITMICHI2"
# 
# Check
if [ x"$CRABTOOLKITMICHI2" = x"" ]
    exit
end
#
# PATH
if not contains "$CRABTOOLKITMICHI2/bin" $PATH
    set -x PATH "$CRABTOOLKITMICHI2/bin" $PATH
end
if not contains "$CRABTOOLKITMICHI2/lib/idl" $IDL_PATH
    set -x IDL_PATH "$CRABTOOLKITMICHI2/lib/idl" $IDL_PATH
end
#
# LIST
set -x CRABTOOLKITCMD "michi2-deploy-files" "michi2-run-fitting-5-components" "michi2-plot-results-of-fitting-5-components"
# 
# CHECK
# -- 20160427 only for interactive shell
# -- http://stackoverflow.com/questions/12440287/scp-doesnt-work-when-echo-in-bashrc
if status --is-interactive
  for TEMPTOOLKITCMD in {$CRABTOOLKITCMD}
    type $TEMPTOOLKITCMD
  end
end


