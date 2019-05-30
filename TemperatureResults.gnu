set terminal pngcairo dashed enhanced size 1920,1080
set output 'TemperatureResult.png'

set title 'Magnetization/Spin'
set grid

set grid ytics mytics  
set mytics 4         

set grid xtics mxtics  
set mxtics 4         


set multiplot layout 2,2 rowsfirst title "MC Transient Temperature Results"
set ylabel 'Magnetization (m)'
set xlabel 'Temperature (T)'
set style line 1 linecolor rgb '#33B2FF' 
set style circle radius 0.025
set style fill solid 1 noborder
plot 'TempDependant.dat'using 1:6 with lines ls 1 lw 2 t "Exact Solution", "" using 1:2 with circles t "MC" lc rgb "blue"

set key off
set title 'Energy/Spin'
set ylabel 'Energy (e)'
set xlabel 'Temperature (T)'
set grid
set style line 1 linecolor rgb '#ba00ff' linetype 1 linewidth 2
plot 'TempDependant.dat'using 1:3 with circles lc 4 lw 2

set title 'Magnetic Susceptibility'
set ylabel 'Susceptibility (χ)'
set xlabel 'Temperature (T)'
set grid
set style line 1 linecolor rgb '#33f500' linetype 1 linewidth 2
plot 'TempDependant.dat'using 1:4 with circles lc 2 lw 2

set title 'Spacific Heat'
set ylabel 'Spacific Heat (Cv)'
set xlabel 'Temperature (T)'
set grid
set style line 1 linecolor rgb '#fea900' linetype 1 linewidth 2
plot 'TempDependant.dat'using 1:5 with circles lc 8 lw 2
set xlabel 'Temperature (K)'

unset multiplot
unset output






  
