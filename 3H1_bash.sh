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
cd ..