set terminal postscript eps enhanced color lw 3 16
set size 0.65,0.65

set key left
f(x) = a*x + b
set format y %4.2f

ppn = 12
ns_per_cycle = 0.001 # each cycle 

set xrange[6:150]
set ytics nomirror
set y2tics
set y2label "time per simulation (s)"
set y2range[100:200]
set title "Scaling of Parallel Independent AMBER QM/MM MD Simulations"

#
# Iit is common for chemists to use simulation speed as a metric.
#
set ylabel "simulation speed (ns/day)"
set xlabel "number of cores"
fit f(x) 'qmmm.dat' u 1:($2*ns_per_cycle*24) via a,b
plot 'qmmm.dat' u 1:($2*ns_per_cycle*24) ti "", f(x) ti sprintf("y = %.2fx+%.2f",a,b), 'qmmm.dat' u 1:5:6 w yerrorbars axes x1y2 ti ""

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

set output "qmmm.eps"
replot