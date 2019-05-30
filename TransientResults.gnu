set terminal pngcairo dashed enhanced size 1920,1080
set output 'TimeResult.png'

set grid ytics mytics  # draw lines for each ytics and mytics
set mytics 5          # set the spacing for the mytics

set grid xtics mxtics  # draw lines for each ytics and mytics
set mxtics 5          # set the spacing for the mytics



set multiplot layout 2,1 rowsfirst title "Equilibrium Iteration Results"
#set yrange [-1.2:1.2]

set title 'Magnetization'
set xlabel 'Time (t)'
set ylabel 'Magnetization (T)'
set grid
set key off
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2
plot 'TimeDependant.dat'using 1:2 with lines ls 6 lc 7 lw 2

set title 'Energy Plot'
set xlabel 'Time (t)'
set ylabel 'Energy (J)'
#set yrange [-2.2:2.2]
set grid
plot 'TimeDependant.dat'using 1:3 with lines ls 8 lc 9 lw 2

unset multiplot
unset output





  
