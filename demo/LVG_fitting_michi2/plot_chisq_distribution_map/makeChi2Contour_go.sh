#!/bin/bash
#

idl -quiet << EOF
.compile makeChi2MinArray.pro makeChi2Contour.pro makeChi2ContourTkVnH2.pro
makeChi2ContourTkVnH2
EOF
