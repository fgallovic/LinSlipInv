DL=1.
DW=1.

set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )
set xlabel 'Along strike (km)'
set ylabel 'Up-dip (km)'
set cblabel 'Mu (Pa)'

plot 'NEZsor.mu' matrix u ($1*DL+DL/2.):($2*DW+DW/2.):3 w image

pause -1
