cd H30
cd 12teta
gfortran -O 12teta.f90
mv a.out T12.out
nohup ./T12.out &
cd ..

cd 14teta
gfortran -O 14teta.f90
mv a.out T14.out
nohup ./T14.out &
cd ..

cd 16teta
gfortran -O 16teta.f90
mv a.out T16.out
nohup ./T16.out &
cd ..

cd 20teta
gfortran -O 20teta.f90
mv a.out T20.out
nohup ./T20.out &
cd ..

cd 22teta
gfortran -O 22teta.f90
mv a.out T22.out
nohup ./T22.out &
cd ..

cd 24teta
gfortran -O 24teta.f90
mv a.out T24.out
nohup ./T24.out &
cd ..
cd ..

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