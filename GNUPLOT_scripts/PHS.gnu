# Generate the pulse height spectrum image using detected_***.dat file
#
# Replace 'filename.dat' with the detected_.. name. 
# Specify the bin size for constant 'B'.
#
#	@file    PHS.gnu
#       @author  Diksha Sharma (Diksha.Sharma@fda.hhs.gov)
#       @date    Apr 13, 2012
#


reset
clear

B=10	# set the bin size here. From hybridMANTIS_input.in calculate, (max photons - min photons detected)/(# of bins).

set style line 1 linetype 1 linewidth 4
set xlabel '# detected photons per primary'
set ylabel 'Frequency'
set si sq
unset key

plot 'filename.dat' u ($0*B):1 w l ls 1

set terminal png enh font "Arial" 26 size 1028,1028 crop
set out 'PHS.png'
rep

set term x11
rep

