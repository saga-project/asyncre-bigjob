set terminal postscript eps enhanced color lw 3 16
set size 0.65,0.65

set key left reverse
ideal_2min = 

ns_per_step = 1e-6  

#set ytics nomirror
#set y2tics
#set y2label "time per simulation (s)"
#set y2range[100:200]
# set title "Scaling of Parallel Independent AMBER MM MD Simulations"

# set ylabel "simulation speed (ns/day)"
# set xlabel "# of coordinated simulations"
# set xrange[6:140] 
# plot 'amber_mm_2MinCycles.dat' u 2:($4*$3*ns_per_step*24) ti "2 min.", 'amber_mm_1MinCycles.dat' u 2:($4*$3*ns_per_step*24) ti "1 min.", 'amber_mm_30SCycles.dat' u 2:($4*$3*ns_per_step*24) ti "30 sec."

# 'amber_mm_2MinCycles.dat' u 2:5:6 w yerrorbars axes x1y2 ti "" lc 1, 'amber_mm_1MinCycles.dat' u 2:5:6 w yerrorbars axes x1y2 ti "" lc 2 


set ylabel "simulation time (ns/day)"
set xlabel "cores per simulation"
plot 'amber_mm_2MinCycles.dat' u 1:($4*$3*ns_per_step*24) ti "2 min.", 'amber_mm_1MinCycles.dat' u 1:($4*$3*ns_per_step*24) ti "1 min.", 'amber_mm_30SCycles.dat' u 1:($4*$3*ns_per_step*24) ti "30 sec."


set output "amber_mm.eps"
replot