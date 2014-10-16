program PREPARE

implicit none

integer, parameter :: ngmax=1000
real, parameter :: pi=3.141592654
real df,aw1,t0
complex rseis,ui,freq
integer ng1,ng2,np
real strike,dip,rake
real leng,widt
real hhypo
real gleng,gwidt
real TM(3,3),ITM(3,3)
real NEZhypo(3),hypo(3),xi(3),sour(3)
real dx1,dx2,x1a(ngmax),x2a(ngmax)
integer i,j,k
interface
  function Transf(NEZ,smer)
    logical smer
    real Transf(3)
    real NEZ(3)
  end function Transf
end interface
common /transform/ TM,NEZhypo,ITM,hypo
integer nc,nfreq,nr,ns,ikmax
real Stat(2)
real tl,aw,xl,uconv,fref
namelist  /input/ nc,nfreq,tl,aw,nr,ns,xl,ikmax,uconv,fref
CHARACTER*6 filename

open(1,file='input.dat')

read(1,*)
read(1,*) nfreq
read(1,*)
read(1,*) tl
read(1,*)
read(1,*) t0
read(1,*)
read(1,*) nr
read(1,*)
read(1,*) ng2,ng1
read(1,*)
read(1,*)
read(1,*)
read(1,*) strike,dip,rake
read(1,*)
read(1,*) hhypo
read(1,*)
read(1,*) leng,widt
read(1,*)
read(1,*) hypo(2),hypo(1)   !WARNING! IN PREVIOUS VERSION THE ORDER WAS hypo(1),hypo(2)!
read(1,*)
read(1,*) np
close(1)

if(ng2>ngmax.or.ng1>ngmax)stop 'Check dimensions!'

hypo(3)=0.

open(2,file='GRDAT.HED')
aw=1.;ns=1;xl=3000000.;ikmax=200000;uconv=1.E-12;fref=1. !Axitra values that does not have to be generally changed
nc=0;  !set up formal value that are actualy not used by Axitra (they are readed from elsewhere)
write(2,input)
close(2)

NEZhypo(1)=0.
NEZhypo(2)=0.
NEZhypo(3)=-hhypo

TM(1,1)=sind(strike)*cosd(dip)
TM(1,2)=-cosd(strike)*cosd(dip)
TM(1,3)=sind(dip)
TM(2,1)=cosd(strike)
TM(2,2)=sind(strike)
TM(2,3)=0.
TM(3,1)=-sind(strike)*sind(dip)
TM(3,2)=cosd(strike)*sind(dip)
TM(3,3)=cosd(dip)
ITM=transpose(TM)

! Change 3.10.2014
dx1=widt/float((ng1))
dx2=leng/float((ng2))

open(3,file='XYGreen.dat')
do i=1,ng1
  x1a(i)=(float(i)-.5)*dx1
  write(3,*) x1a(i)
enddo
do i=1,ng2
  x2a(i)=(float(i)-.5)*dx2
  write(3,*) x2a(i)
enddo

open(1,file='sources.dat')
k=0
do i=1,ng1
  do j=1,ng2
    xi(1)=x1a(i)
    xi(2)=x2a(j)
    xi(3)=0.
    k=k+1
    if(k<10)then
      write(filename,'(A5,I1)')'00000',k
    elseif(k<100)then
      write(filename,'(A4,I2)')'0000',k
    elseif(k<1000)then
      write(filename,'(A3,I3)')'000',k
    elseif(k<10000)then
      write(filename,'(A2,I4)')'00',k
    elseif(k<100000)then
      write(filename,'(A1,I5)')'0',k
    else
      write(filename,'(I6)')k
    endif
    sour=Transf(xi,.FALSE.)
    sour(3)=-sour(3)
    sour=sour/1000.
    write(1,'(A6,6E13.5)') filename,sour,strike,dip,rake
  enddo
enddo
close(1)
open(1,file='fault.dat')
xi=0.
write(1,*) Transf(xi,.FALSE.)
xi(1)=widt
xi(2)=0.
xi(3)=0.
write(1,*) Transf(xi,.FALSE.)
xi(1)=widt
xi(2)=leng
xi(3)=0.
write(1,*) Transf(xi,.FALSE.)
xi(1)=0.
xi(2)=leng
xi(3)=0.
write(1,*) Transf(xi,.FALSE.)
xi=0.
write(1,*) Transf(xi,.FALSE.)

close(1)    

ui=cmplx(0.,1.)
df=1./tl
aw1=-aw/(2.*tl)
open(1,file='dirac.dat')
do i=1,nfreq
  freq=cmplx(df*(i-1),aw1)
  rseis=exp(-ui*2.*pi*t0*freq)
  write(1,*) real(rseis),imag(rseis)
enddo
do i=nfreq+1,np
  rseis=0.
  write(1,*) real(rseis),imag(rseis)
enddo
close(1)        

end

function Transf(NEZ,smer)
implicit none
logical smer
real Transf(3)
real NEZhypo(3)
real TM(3,3),ITM(3,3),hypo(3)
real NEZ(3)
common /transform/ TM,NEZhypo,ITM,hypo

if (smer) then
  Transf=matmul(TM,(NEZ-NEZhypo))+hypo
else
  Transf=matmul(ITM,(NEZ-hypo))+NEZhypo
endif
end function
