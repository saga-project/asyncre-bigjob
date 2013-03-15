set terminal postscript eps enhanced color lw 3 16
set size 0.65,0.65

set key left reverse
ideal_2min = 

ns_per_step = 1e-6  

#set ytics nomirror
#set y2tics
#set y2label "time per simulation (s)"
#set y2range[100:200]
set title "AMBER RE-US MD Simulations (576 replicas/720 cores)"

set ylabel "simulation speed (ns/day/replica)"
set xlabel "# of coordinated simulations"
set xrange[6:140] 
set yrange[0:1650]
#plot 'amber_mm_2MinCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "2 min.", 'amber_mm_1MinCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "1 min.", 'amber_mm_30SCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "30 sec.", 'amber_mm_2MinCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "" w l lt 1, 'amber_mm_1MinCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "" w l lt 2, 'amber_mm_30SCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "" w l lt 3

set ylabel "simulation speed (ns/day)"
plot 'amber_mm_2MinCycles.dat' u 2:($4*$3*ns_per_step*24) ti "2 min.", 'amber_mm_1MinCycles.dat' u 2:($4*$3*ns_per_step*24) ti "1 min.", 'amber_mm_30SCycles.dat' u 2:($4*$3*ns_per_step*24) ti "30 sec.", 'amber_mm_2MinCycles.dat' u 2:($4*$3*ns_per_step*24) ti "" w l lt 1, 'amber_mm_1MinCycles.dat' u 2:($4*$3*ns_per_step*24) ti "" w l lt 2, 'amber_mm_30SCycles.dat' u 2:($4*$3*ns_per_step*24) ti "" w l lt 3, 'amber_mm_ideal.dat' u (720/$1):($2*720/$1):(($2+$3)*720/$1) lt 1 lc 7 ti "" w l


#set ylabel "simulation time (ns/day/replica)"
#set xlabel "cores per simulation"
#plot 'amber_mm_2MinCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "2 min.", 'amber_mm_1MinCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "1 min.", 'amber_mm_30SCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "30 sec.", 'amber_mm_2MinCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "" w l lt 1, 'amber_mm_1MinCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "" w l lt 2, 'amber_mm_30SCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "" w l lt 3


set output "amber_mm.eps"
replot