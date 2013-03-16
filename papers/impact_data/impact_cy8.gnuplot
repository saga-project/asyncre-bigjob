#!/usr/local/bin/gnuplot
 
set title "IMPACT cy8 Small System (384 core BJ)"
set xlabel "# of replicas"
set ylabel "simulation speed (ns/day)"
#set ytics "10"
plot "impact_cy8.dat" title "" with lines
pause -1 "Hit any key to continue"