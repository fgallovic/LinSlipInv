source load_intel
ifx -oslipratesonfault slipratesonfault.f90
./slipratesonfault
gnuplot slipratesonfault.plot.gp
