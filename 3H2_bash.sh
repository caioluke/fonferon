cd H30
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