set terminal postscript eps enhanced color lw 3 16
set size 0.65,0.65

set key left #reverse
set key title "{/Symbol w} (s^{-1})"

ns_per_step = 1e-6  

#set ytics nomirror
#set y2tics
#set y2label "time per simulation (s)"
#set y2range[100:200]
set title "AMBER RE-US MD (576 total replicas/720 cores)"

set ylabel "simulation speed (ns/day/replica)"
set xlabel "# of concurrent replicas"
set xrange[6:126] 
set yrange[0:1950]
#plot 'amber_mm_2MinCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "2 min.", 'amber_mm_1MinCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "1 min.", 'amber_mm_30SCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "30 sec.", 'amber_mm_2MinCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "" w l lt 1, 'amber_mm_1MinCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "" w l lt 2, 'amber_mm_30SCycles.dat' u 2:($4*$3*ns_per_step*24/$2) ti "" w l lt 3

set ylabel "simulation speed (ns/day)"

plot 'amber_mm_ideal.dat' u (720/$1):($2*720/$1):(($2+$3)*720/$1) lt 1 lc 7 ti "(ideal) 0" w l, 'amber_mm_2MinCycles.dat' u 2:($4*$3*ns_per_step*24) lc 1 pt 6 ti "", 'amber_mm_1MinCycles.dat' u 2:($4*$3*ns_per_step*24) lc 2 pt 4 ti "", 'amber_mm_30SCycles.dat' u 2:($4*$3*ns_per_step*24) lc 3 pt 8 ti "", 'amber_mm_2MinCycles.dat' u 2:($4*$3*ns_per_step*24) ti "30" w l lt 1 lc 1, 'amber_mm_1MinCycles.dat' u 2:($4*$3*ns_per_step*24) ti "60" w l lt 1 lc 2, 'amber_mm_30SCycles.dat' u 2:($4*$3*ns_per_step*24) ti "120" w l lt 1 lc 3

#set ylabel "simulation time (ns/day/replica)"
#set xlabel "cores per simulation"
#plot 'amber_mm_2MinCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "2 min.", 'amber_mm_1MinCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "1 min.", 'amber_mm_30SCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "30 sec.", 'amber_mm_2MinCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "" w l lt 1, 'amber_mm_1MinCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "" w l lt 2, 'amber_mm_30SCycles.dat' u 1:($4*$3*ns_per_step*24/$2) ti "" w l lt 3


set output "amber_mm.eps"
replot
