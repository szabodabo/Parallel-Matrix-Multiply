set term postscript eps enhanced color

set xlabel "Number of Cores"
set ylabel "Execution Time (sec)"

set output "graphs/time.eps"
plot "out.csv" using 1:2 title "Execution Time vs. Number of Cores" with points

set ylabel "Maximum Bandwidth (MB/sec)"
set output "graphs/max.eps"
set xrange [1:16]
set yrange [0:500]
plot "out.csv" using 1:5 title "Maximum Send Bandwidth" with points, \
     "out.csv" using 1:6 title "Maximum Receive Bandwidth" with points

set ylabel "Minimum Bandwidth (MB/sec)"
set output "graphs/min.eps"
set xrange [1:16]
set yrange [0:75]
plot "out.csv" using 1:3 title "Minimum Send Bandwidth" with points, \
     "out.csv" using 1:4 title "Minimum Receive Bandwidth" with points

set ylabel "Average Bandwidth (MB/sec)"
set output "graphs/avg.eps"
set xrange [1:16]
set yrange [0:200]
plot "out.csv" using 1:7 title "Average Send Bandwidth" with points, \
     "out.csv" using 1:8 title "Average Receive Bandwidth" with points
