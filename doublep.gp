set terminal gif animate delay 1
set output "double_pendulum.gif"

stats 'results'

n = STATS_records
k=500
step=n/k

set xrange [-2.1:2.1]
set yrange [-2.1:2.1]
set size square

do for [i=2:n:step] {
   set multiplot
   plot "results" every ::1::i with points pt 7 ps 0.2  notitle
   unset multiplot
}
