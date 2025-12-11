source load_intel
ifx -oSlipInv-seisplot SlipInv-seisplot.f90
./SlipInv-seisplot
gnuplot SlipInv-seisplot.gp
