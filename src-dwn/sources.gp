set term x11
set size ratio -1
plot 'sources.dat' u 3:2 w p,'fault.dat' u 2:1 w l,'stations.dat' u 2:1  w p
pause 35
