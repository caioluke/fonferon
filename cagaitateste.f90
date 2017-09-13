program testecagaita
implicit none

open (unit=21,file='helloworld.txt',status='unknown')

write(21,*) "Hello world! A vida é bela"

close(unit=21)

end program testecagaita