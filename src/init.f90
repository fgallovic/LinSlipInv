!   Truncated SVD solution of the Linear slip inversion problem
!   Requires as an input the output from SlipInvSVD1
!   AUTHOR:  Frantisek Gallovic

    MODULE SISVDmodule
      REAL*8,PARAMETER:: PI=3.1415926535d0
      INTEGER Msvd,Nsvd    !Msvd - number of model parameters, Nsvd - number of data
      INTEGER Nsmooth,Ssvd,NSTAcomp,Nseis,Ngps,Nslip
      INTEGER NL,NW,nfmax,NRseis,NRgps,np
      REAL*8 T,TS,T1,T2,T0,artifDT,Mfix,strike,dip,hypodepth,leng,widt,epicW,epicL,vr
      INTEGER iT1,iT2,nT,iT0
      REAL*8,ALLOCATABLE,DIMENSION(:):: fc1,fc2
      INTEGER,ALLOCATABLE,DIMENSION(:):: fcsta
      INTEGER nfc
      REAL*8 dt,df,dL,dW,elem
      REAL*8 smoothkoef,smoothkoefGF,relatweightGPS,lambdalim,maxw
      REAL*8 norminput,normdatGPS
      CHARACTER*256 inputtffile
      INTEGER syntdata,syntdatai,syntdataj,compweights,fixM0weight,slipweight
      INTEGER eigsumchoice,lambdafrom,lambdato,lambdanum
      INTEGER minMNsvd,minSVDchoice
      REAL*8,ALLOCATABLE:: staweight(:,:)
      INTEGER,ALLOCATABLE:: stainfo(:,:)
      REAL*8,ALLOCATABLE:: G(:,:),mu(:,:),lambda(:,:)
      REAL*8,ALLOCATABLE:: D(:),normdat(:)
      REAL*8,ALLOCATABLE:: V(:,:),VT(:,:),U(:,:),W(:)
      REAL*8,ALLOCATABLE:: Dout(:,:),M(:,:)
      REAL*4,ALLOCATABLE,DIMENSION(:,:):: sigmaGPS,sourNgps,sourEgps,sourZgps,strikeGPS,dipGPS,rakeGPS
    END MODULE

    SUBROUTINE Init()
    USE SISVDmodule
    IMPLICIT NONE
    INTEGER ndepth
    REAL*8 dum,R,staN,staE
    REAL*8,ALLOCATABLE:: depth(:),vp(:),vs(:),rho(:)
    INTEGER i,j

    write(*,*)'Reading parameters...'

    open(10,file='input.dat',action='read')
    read(10,*)
    read(10,*) nfmax
    read(10,*)
    read(10,*) T,TS,T1,T2
    read(10,*)
    read(10,*) artifDT
    read(10,*)
    read(10,*) NRseis,NRgps
    read(10,*)
    read(10,*) NL,NW
    read(10,*)
    read(10,*) Mfix
    read(10,*)
    read(10,*) strike,dip
    read(10,*)
    read(10,*) hypodepth
    read(10,*)
    read(10,*) leng,widt
    read(10,*)
    write(*,*) '  (Warning! Assumig epicL, epicW in input.dat!)'
    read(10,*) epicL,epicW
    read(10,*)
    read(10,*) np
    read(10,*)
    read(10,*) vr
    read(10,*)
    read(10,*) nfc   !number of frequency bands
    allocate(fc1(nfc),fc2(nfc))
    do i=1,nfc
      read(10,*) fc1(i),fc2(i)
    enddo
    close(10)
    dL=leng/dble(NL)
    dW=widt/dble(NW)
    elem=dL*dW

    if(NRseis>0)then
      dt=T/float(np)
      df=1./T
      iT1=T1/dt+1
      iT2=T2/dt+1
      nT=iT2-iT1+1
      open(10,file='stainfo.dat',action='read')
      allocate(stainfo(3,NRseis),staweight(3,NRseis),fcsta(NRseis))
      do i=1,NRseis
        read(10,*)stainfo(:,i),staweight(:,i),fcsta(i)
      enddo
      close(10)
    endif

    open(10,file='SlipInvSVD.in',action='read')
    read(10,*)
    read(10,*)syntdata
    if(syntdata==1 )read(10,*)syntdatai,syntdataj
    if(syntdata==-1)read(10,*)inputtffile
    read(10,*)
    read(10,*)smoothkoef,smoothkoefGF,relatweightGPS,fixM0weight,slipweight
    read(10,*)
    read(10,*)compweights
    read(10,*)
    read(10,*)eigsumchoice
    read(10,*)
    if(eigsumchoice==1)then
      read(10,*)lambdalim
      lambdanum=1
    else
      read(10,*)lambdafrom,lambdato
      lambdanum=lambdato-lambdafrom+1
    endif
    if(abs(smoothkoef)>0.d0)lambdanum=1
    read(10,*)
    read(10,*)T0
    iT0=T0/dt
	  read(10,*)
	  read(10,*)minSVDchoice
    close(10)

! Evaluating mu
    allocate(mu(NL,NW),lambda(NL,NW))
    open(10,FILE='crustal.dat',ACTION='READ',STATUS='OLD',ERR=181)
    write(*,*)'  (Using mu values from file crustal.dat)'
    read(10,*)
    read(10,*)
    read(10,*)ndepth
    allocate(depth(ndepth),vp(ndepth),vs(ndepth),rho(ndepth))
    read(10,*)
    read(10,*)
    do i=1,ndepth
      read(10,*)depth(i),vp(i),vs(i),rho(i)
    enddo
    close(10)
    do i=1,NW
      dum=(hypodepth+(epicW-dW*(dble(i)-.5d0))*sin(dip/180.d0*PI))/1000.d0
      if(dum>depth(ndepth))then
        mu(:,i)=rho(ndepth)*vs(ndepth)**2*1.d9
        lambda(:,i)=rho(ndepth)*vp(ndepth)**2*1.d9-2*mu(:,i)
      else
        do j=1,ndepth
          if(dum<depth(j))exit
        enddo
        mu(:,i)=rho(j-1)*vs(j-1)**2*1.d9
        lambda(:,i)=rho(j-1)*vp(j-1)**2*1.d9-2*mu(:,i)
      endif
    enddo
    deallocate(depth,vp,vs,rho)
    goto 183

181 open(10,FILE='NEZsor.mu',ACTION='READ',STATUS='OLD',ERR=182)
    write(*,*)'  (Using mu values from file NEZsor.mu)'
    do j=1,NW
      read(10,*)(mu(i,j),i=1,NL)
    enddo
    lambda=0.
    goto 183

182 write(*,*)'ERROR - neither crustal.dat nor NEZsor.mu found!'
    stop

183 norminput=Mfix*dble(NL*NW)/sum(mu(:,:))/elem/dt

! Location of subsources (needed for GPS and SCRMOD output format)
    open(10,FILE='sources.dat')
    allocate(sourNgps(NL,NW),sourEgps(NL,NW),sourZgps(NL,NW),strikeGPS(NL,NW),dipGPS(NL,NW),rakeGPS(NL,NW))
    do j=1,NW
      do i=1,NL
        read(10,*)dum,sourNgps(i,j),sourEgps(i,j),sourZgps(i,j),strikeGPS(i,j),dipGPS(i,j),rakeGPS(i,j)
      enddo
    enddo
    close(10)
    sourNgps=sourNgps*1.e3
    sourEgps=sourEgps*1.e3
    sourZgps=sourZgps*1.e3

!Applying distance dependent weigths
    if(compweights==0)then
       write(*,*)'  (Applying distance-dependent weigths)'
       open(10,FILE='stations.dat')
       open(11,FILE='stainfo.out')
       do i=1,NRseis
         read(10,*)staN,staE
         R=minval(sqrt((sourNgps(:,:)-staN*1.e3)**2+(sourEgps(:,:)-staE*1.e3)**2))
!         dum=max(R,10.e3)/10.e3
         dum=max(R,leng/4.)/(leng/4.)
         staweight(:,i)=staweight(:,i)*dum
         write(11,*)staN,staE,dum
       enddo
       close(10)
       close(11)
    endif

! Main allocations

    if(NRseis>0)then
      Ssvd=int(TS/dt+1.d0)       !number of time samples of the slip velocity model
      if(iT1-iT0-Ssvd<0.or.iT2-iT0>np)then
        write(*,*)'Error! Time window out of range, check input.dat...'
        stop
      endif
    else
      dt=1.d0
      df=1.d0
      np=0
      Ssvd=1                      !Slip value
    endif

    Msvd=NW*NL*Ssvd               !number of model parameters
    if(smoothkoef<0.d0)then       !smoothing by means of first differences
      Nsmooth=(Ssvd-1)*NL*NW+Ssvd*(NL-1)*NW+Ssvd*NL*(NW-1)
    elseif(smoothkoef>0.d0)then   !smoothing by means of a covariance function 
      Nsmooth=Msvd
    else                          !no smoothing
      Nsmooth=0
    endif

    Ngps=NRgps*3
    if(NRseis>0)then
      NSTAcomp=sum(stainfo(:,:))
      Nseis=nT*NSTAcomp
    else
      NSTAcomp=0
      Nseis=0
    endif

    if(slipweight>0.d0)then
      Nslip=Msvd
    else
      Nslip=0
    endif
    
    Nsvd=Nseis+Ngps+1+Nsmooth+Nslip

	if(minSVDchoice==1)then
	  minMNsvd=min(Msvd,Nsvd)
	else
    minMNsvd=Msvd
	endif

    END
