

#!/usr/local/bin/gnuplot
 
set terminal postscript eps enhanced color lw 3 16
set size 0.65,0.65


set title "IMPACT TrpCage Medium System (504 cores)"
set xlabel "# of replicas"
set ylabel "simulation speed (ns/day)"
#set ytics "10"
set xrange[50:1100]
plot "impact_trpcage.dat" ti "" pt 6 lc 1, "impact_trpcage.dat" ti "" lt 1 lc 1 w l

set output "impact_trpcage.eps"
replot
