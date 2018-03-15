cd H10
cd teta10
gfortran -O 10teta.f90
mv a.out T10.out
nohup ./T10.out &
cd ..
cd teta15
gfortran -O 15teta.f90
mv a.out T15.out
nohup ./T15.out &
cd ..
cd teta20
gfortran -O 20teta.f90
mv a.out T20.out
nohup ./T20.out &
cd ..
cd ..
cd H10MAG
cd teta10MAG
gfortran -O 10tetaMAG.f90
mv a.out T10MAG.out
nohup ./T10MAG.out &
cd ..
cd teta15MAG
gfortran -O 15tetaMAG.f90
mv a.out T15MAG.out
nohup ./T15MAG.out &
cd ..
cd teta20MAG
gfortran -O 20tetaMAG.f90
mv a.out T20MAG.out
nohup ./T20MAG.out &
cd ..
cd ..