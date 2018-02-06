#!/bin/bash
# 

#ls ../../

idl << EOF
.r read_spec
.r makelibSED_FSPS_CSP
resolve_routine, 'read_spec', /COMPILE_FULL_FILE
makelibSED_FSPS_CSP
EOF

