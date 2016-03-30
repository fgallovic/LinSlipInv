program Resort

implicit none

integer i,k,kk
integer np,nr,rec,dum1,dum2,nfmax,pom,NSeg,gntot
integer,allocatable:: gn1(:),gn2(:)
real dat(7),T,dum
CHARACTER*6 filename

open(3,form='unformatted',file='NEZsor.dat')
open(5,file='input.dat')

read(5,*)
read(5,*) nfmax
read(5,*)
read(5,*) T
read(5,*)
read(5,*) dum, NSeg
allocate(gn1(NSeg),gn2(NSeg))
read(5,*)
read(5,*) nr
read(5,*)
read(5,*) (gn2(kk),gn1(kk),kk=1,NSeg)
do i=1,13
  read(5,*)
enddo
read(5,*) np

gntot=sum(gn1(:)*gn2(:))
do i=1,gntot
  if(i<10)then
    write(filename,'(A5,I1)')'00000',i
  elseif(i<100)then
    write(filename,'(A4,I2)')'0000',i
  elseif(i<1000)then
    write(filename,'(A3,I3)')'000',i
  elseif(i<10000)then
    write(filename,'(A2,I4)')'00',i
  elseif(i<100000)then
    write(filename,'(A1,I5)')'0',i
  else
    write(filename,'(I6)')i
  endif
  open(100+i,FILE='dat/'//filename//'.nez')
enddo

write(*,*) nfmax,nr,NSeg
do kk=1,NSeg
  do rec=1,nr
    write(*,*) rec
    write(3) rec
    do i=sum(gn1(1:kk-1)*gn2(1:kk-1))+1,sum(gn1(1:kk)*gn2(1:kk))
      read(100+i,*)
      write(3) i
      do k=1,nfmax
        read(100+i,*) dat(1:7)
        write(3) dat(2:7)
      enddo
    enddo
  enddo
enddo

end



