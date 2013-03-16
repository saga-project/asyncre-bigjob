set terminal postscript eps enhanced color lw 3 16
set output "qmmm.eps"
set size 0.70,0.70

set key left
f(x) = a*x + b
g(x) = c*x + d
set format y "%4.2f"

ppn = 12
ns_per_cycle = 0.001 # each cycle 

set xrange[6:150]
set ytics nomirror
set y2tics
set y2label "time per simulation (s)"
set y2range[100:200]
set title "Scaling of AMBER QM/MM MD Simulations"

#
# It is common for chemists to use simulation speed as a metric.
#
set ylabel "simulation speed (ns/day)"
set xlabel "# of concurrent replicas"
fit f(x) 'qmmm.dat' u 1:($2*ns_per_cycle*24) via a,b
fit g(x) 'qmmm.dat' u 1:5:6 via c,d
#plot 'qmmm.dat' u 1:($2*ns_per_cycle*24) ti "", f(x) ti sprintf("y = %.2fx+%.2f",a,b), 'qmmm.dat' u 1:5:6 w yerrorbars axes x1y2 ti ""
plot 'qmmm.dat' u 1:($2*ns_per_cycle*24) pt 6 lc 1 ti "", f(x) lt 2 lc 1 ti "", 'qmmm.dat' u 1:5:6 w yerrorbars axes x1y2 pt 4 lt 1 lc 2 ti "", g(x) lt 2 lc 2 ti "" axes x1y2


#
# Alternatively you can use simulation time if that seems more natural.
#
# set ylabel "simulation time (ns)"
# set xlabel "number of cores"
# fit f(x) 'qmmm.dat' u 1:($2*ns_per_cycle) via a,b
# plot 'qmmm.dat' u 1:($2*ns_per_cycle) ti ""


#
# You might also want to count in terms of compute nodes if that is better.
#
# set ylabel "simulation time (ns)"
# set xlabel "number of nodes"
# fit f(x) 'qmmm.dat' u ($1/ppn):($2*ns_per_cycle) via a,b
# plot 'qmmm.dat' u ($1/ppn):($2*ns_per_cycle) ti ""

replot
