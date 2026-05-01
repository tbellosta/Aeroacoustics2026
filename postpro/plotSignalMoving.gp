set terminal pdfcairo size 10cm,5cm enhanced font "Times,12"
set output 'signal.pdf'

set xlabel 'Time [s]'
set ylabel "p' [Pa]"
set grid
set key top right
set samples 5000


plot '../bin/signal_01.dat' using 1:2 with lines lw 3 lc rgb '#2060A0' title "FWH numerical"
