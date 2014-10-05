#!/bin/bash

processors=6

cat sources.dat | xargs -P$processors -n7 ./gr_nez
cat sources.dat | xargs -P$processors -n7 ./cnv_nez
./resort

echo "Done!"
