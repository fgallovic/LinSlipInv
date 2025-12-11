source load_intel
ifx -fast -autodouble -ocnv_nez cnv_nez.for
ifx -fast -autodouble -ogr_nez gr_nez.for
ifx -fast -oprepare prepare.f90
ifx -fast -oresort resort.f90

./prepare
rm -fr dat
mkdir dat
