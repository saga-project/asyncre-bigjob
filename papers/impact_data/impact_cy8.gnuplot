#!/usr/local/bin/gnuplot
 
set terminal postscript eps enhanced color lw 3 16
set size 0.65,0.65
set output "impact_cy8.eps"


set title "IMPACT cy8 Small System (384 core BJ)"
set xlabel "# of replicas"
set ylabel "simulation speed (ns/day)"
#set ytics "10"
plot "impact_cy8.dat" ti "" w l lt 1
