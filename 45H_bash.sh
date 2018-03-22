cd H40
gfortran -O H40_posinit.f90
mv a.out H40.out
nohup ./H40.out &
cd ..

cd H50
gfortran -O H50_posinit.f90
mv a.out H50.out
nohup ./H50.out &
cd ..