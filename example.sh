#!/bin/csh
#
# Example script for the `xplat' program...
#
# As long as we don't invoke the 'g' command, then 'xplat' can be used
# quite well in the background (of course this could be run with plain
# 'plat' as well).
#
# First of all some interesting statistics are selected via the
# commands `I -w <information>'.  These will be written to the file
# `info.pb'.  Then a Voronoi network of 10 cells is created `V 10'.
# The network is then subject to Hencky strain, first squeezing
# horizontally and then vertically.
#
# A number of foam pictures are generated and saved in Postscript format
# in the subdirectory `data'.
#
# Note: If you want to run `xplat' from a script you should generally
# avoid invoking plotting commands as these will cause the program to
# crash.
#

xplat <<!
I -w netenergy
I -w henckyeps
I -w phi
V 10
i
W data/fmraw.pb
e
i
W data/fm0.pb
F 0.98
i
W data/fm_f0.98.pb
{bei} 40
W data/fm_f0.98_b40.pb
q
!
