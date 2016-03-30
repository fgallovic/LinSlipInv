program PREPARE

implicit none

integer, parameter :: ngmax=10000
real, parameter :: pi=3.141592654
real df,aw1,t0
complex rseis,ui,freq
real,allocatable,dimension(:):: strike,dip,rake,leng,widt,hhypo
real,allocatable:: hypo(:,:),offsets(:,:)
real gleng,gwidt
real TM(3,3),ITM(3,3)
real NEZhypo(3),hypo2(3),xi(3),sour(3)
real dx1,dx2
real,allocatable:: x1a(:),x2a(:)
integer,allocatable:: ng1(:),ng2(:)
integer NSeg,np
integer i,j,k,kk
interface
  function Transf(NEZ,smer)
    logical smer
    real Transf(3)
    real NEZ(3)
  end function Transf
end interface
common /transform/ TM,NEZhypo,ITM,hypo2
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
read(1,*) t0,NSeg
allocate(ng1(NSeg),ng2(NSeg),strike(NSeg),dip(NSeg),rake(NSeg),hhypo(NSeg),leng(NSeg),widt(NSeg),hypo(3,NSeg))
allocate(offsets(2,NSeg))
read(1,*)
read(1,*) nr
read(1,*)
read(1,*) (ng2(kk),ng1(kk),kk=1,NSeg)
read(1,*)
read(1,*)
read(1,*)
read(1,*) (strike(kk),dip(kk),rake(kk),kk=1,NSeg)
read(1,*)
read(1,*) (hhypo(kk),kk=1,NSeg)
read(1,*)
read(1,*) (leng(kk),widt(kk),kk=1,NSeg)
read(1,*)
read(1,*) (hypo(2,kk),hypo(1,kk),kk=1,NSeg)   !WARNING! IN PREVIOUS VERSION THE ORDER WAS hypo(1),hypo(2)!
read(1,*)
read(1,*) np
close(1)

offsets=0.
if(NSeg>1)then
  open(1,file='prepare.offsets')
  read(1,*)
  do kk=2,NSeg
    read(1,*)offsets(:,kk)
  enddo
  close(1)
endif

open(1,file='sources.dat')
open(2,file='fault.dat')
k=0
do kk=1,NSeg
  if(ng2(kk)>ngmax.or.ng1(kk)>ngmax)stop 'Check dimensions!'

  hypo(3,kk)=0.

  allocate(x1a(ng1(kk)),x2a(ng2(kk)))

  NEZhypo(1)=offsets(1,kk)
  NEZhypo(2)=offsets(2,kk)
  NEZhypo(3)=-hhypo(kk)
  hypo2=hypo(:,kk)

  TM(1,1)=sind(strike(kk))*cosd(dip(kk))
  TM(1,2)=-cosd(strike(kk))*cosd(dip(kk))
  TM(1,3)=sind(dip(kk))
  TM(2,1)=cosd(strike(kk))
  TM(2,2)=sind(strike(kk))
  TM(2,3)=0.
  TM(3,1)=-sind(strike(kk))*sind(dip(kk))
  TM(3,2)=cosd(strike(kk))*sind(dip(kk))
  TM(3,3)=cosd(dip(kk))
  ITM=transpose(TM)

  dx1=widt(kk)/float(ng1(kk))
  dx2=leng(kk)/float(ng2(kk))

  do i=1,ng1(kk)
    x1a(i)=(float(i)-.5)*dx1
  enddo
  do i=1,ng2(kk)
    x2a(i)=(float(i)-.5)*dx2
  enddo

  do i=1,ng1(kk)
    do j=1,ng2(kk)
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
      write(1,'(A6,6E13.5)') filename,sour,strike(kk),dip(kk),rake(kk)
    enddo
  enddo
  deallocate(x1a,x2a)

  xi=0.
  write(2,*) Transf(xi,.FALSE.)/1000.
  xi(1)=widt(kk)
  xi(2)=0.
  xi(3)=0.
  write(2,*) Transf(xi,.FALSE.)/1000.
  xi(1)=widt(kk)
  xi(2)=leng(kk)
  xi(3)=0.
  write(2,*) Transf(xi,.FALSE.)/1000.
  xi(1)=0.
  xi(2)=leng(kk)
  xi(3)=0.
  write(2,*) Transf(xi,.FALSE.)/1000.
  xi=0.
  write(2,*) Transf(xi,.FALSE.)/1000.
  write(2,*)
enddo
close(1)    
close(2)

open(2,file='GRDAT.HED')
aw=1.;ns=1;xl=3000000.;ikmax=200000;uconv=1.E-12;fref=1. !Axitra values that does not have to be generally changed
nc=0;  !set up formal value that are actualy not used by Axitra (they are readed from elsewhere)
write(2,input)
close(2)

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
real NEZhypo(3),hypo2(3)
real TM(3,3),ITM(3,3)
real NEZ(3)
common /transform/ TM,NEZhypo,ITM,hypo2

if (smer) then
  Transf=matmul(TM,(NEZ-NEZhypo))+hypo2
else
  Transf=matmul(ITM,(NEZ-hypo2))+NEZhypo
endif
end function
