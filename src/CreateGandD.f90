!   AUTHOR:  Frantisek Gallovic
    
    
    SUBROUTINE CreateGandD()
    USE SISVDmodule
    IMPLICIT NONE
    INTEGER i,j,k,l,IRET
    REAL*8 dum
    REAL*8,ALLOCATABLE:: H(:,:),Uvec(:),tfin(:),CM(:,:)
    REAL*4,ALLOCATABLE,DIMENSION(:,:):: stalocGPS,dataGPS
    REAL*4 gpsgfN,gpsgfE,gpsgfZ
    REAL*4 ALPHA,dumGPS,dLgps,dWgps
    REAL*8,ALLOCATABLE:: slipconstraint(:,:)

    if(syntdata==0)then
      write(*,*)'  (reading observed data from files)'
    elseif(syntdata>0)then
      write(*,*)'  (synthetic model data)'
      allocate(tfin(Msvd))
      CALL synthetic_model()
    else
      write(*,*)'  (reading custom target model from file ',trim(inputtffile),')'
      allocate(tfin(Msvd))
      open(10,FILE=trim(inputtffile))
      do i=1,Msvd
        read(10,*)tfin(i)
      enddo
      close(10)
!      CALL synthetic_model()
    endif

    if(NRseis>0)then
      allocate(H(Nseis,Msvd),Uvec(Nseis))
      CALL CreateH()
      CALL CreateUvec()
      write(*,*)'Creating matrix G and vector d...'
! Weights for matrix H and vector u
      normdat(:)=sqrt(sum(Uvec(1:Nseis)**2)) ! unified weights
      k=0
      do i=1,NRseis
        do j=1,3
          if(stainfo(j,i)==0)cycle
          k=k+1
          normdat(k)=normdat(k)/staweight(j,i)
        enddo
      enddo
! Applying weights
      do i=1,NSTAcomp
        G((i-1)*nT+1:i*nT,:)=H((i-1)*nT+1:i*nT,:)/normdat(i)
        D((i-1)*nT+1:i*nT)=Uvec((i-1)*nT+1:i*nT)/normdat(i)
      enddo
      deallocate(H,Uvec)
! If smoothing is on, apply general std.dev for matrix G - same as taking into account diagonal CD
      if(abs(smoothkoef)>0.d0)then
        G(1:Nseis,:)=G(1:Nseis,:)/sum(1.d0/normdat)*sum(staweight(:,:)*dble(stainfo(:,:)))/smoothkoefGF
        D(1:Nseis)=D(1:Nseis)/sum(1.d0/normdat)*sum(staweight(:,:)*dble(stainfo(:,:)))/smoothkoefGF
      endif
!    write(*,*)'  (Mean sigma CD = ',sum(normdat)/dble(NSTAcomp),')'
    endif

! GPS
    if(NRgps>0)then
      open(10,file='stations-GPS.dat',action='read')
      allocate(dataGPS(3,NRgps),stalocGPS(3,NRgps))
      normdatGPS=0.d0
      do i=1,NRgps
        read(10,*)stalocGPS(1:2,i),dataGPS(1:3,i),sigmaGPS(1:3,i)
        stalocGPS(3,i)=0.d0
        sigmaGPS(:,i)=sigmaGPS(:,i)/relatweightGPS
        normdatGPS=normdatGPS+sum((dataGPS(:,i)/sigmaGPS(:,i))**2)
      enddo
      close(10)
      stalocGPS=stalocGPS*1.e3
      normdatGPS=sqrt(normdatGPS)
      dumGPS=100.
      dLgps=real(dL)
      dWgps=real(dW)
!NEFUNGUJE $OMP parallel do private(i,j,k,ALPHA,gpsgfN,gpsgfE,gpsgfZ,IRET) DEFAULT(SHARED) SCHEDULE(DYNAMIC,1)      
      do k=1,NRgps
        do j=1,NW
          do i=1,NL
            ALPHA=(lambda(i,j)+mu(i,j))/(lambda(i,j)+2.*mu(i,j))
            CALL DC3Dmodif(ALPHA,sourNgps(i,j),sourEgps(i,j),sourZgps(i,j),strikeGPS(i,j),dipGPS(i,j),rakeGPS(i,j),dLgps,dWgps,dumGPS,stalocGPS(1,k),stalocGPS(2,k),stalocGPS(3,k),gpsgfN,gpsgfE,gpsgfZ,IRET)
            G(Nseis+(k-1)*3+1,(j-1)*Ssvd*NL+(i-1)*Ssvd+1:(j-1)*Ssvd*NL+i*Ssvd)=gpsgfN*dt/sigmaGPS(1,k)
            G(Nseis+(k-1)*3+2,(j-1)*Ssvd*NL+(i-1)*Ssvd+1:(j-1)*Ssvd*NL+i*Ssvd)=gpsgfE*dt/sigmaGPS(2,k)
            G(Nseis+(k-1)*3+3,(j-1)*Ssvd*NL+(i-1)*Ssvd+1:(j-1)*Ssvd*NL+i*Ssvd)=gpsgfZ*dt/sigmaGPS(3,k)
          enddo
        enddo
        D(Nseis+(k-1)*3+1:Nseis+(k-1)*3+3)=dataGPS(1:3,k)/sigmaGPS(1:3,k)
      enddo
!NEFUNGUJE $OMP end parallel do
      deallocate(dataGPS,stalocGPS)
    endif

! Scalar moment constraint
    if(fixM0weight>0.d0)write(*,*)'  (constraint on M0 applied)'
    G(Nseis+Ngps+1,:)=1.d0/norminput*fixM0weight
    D(Nseis+Ngps+1)=1.d0*fixM0weight

! Smoothing constraints
    if(smoothkoef<0.d0)then
      write(*,*)'  (smoothing by first differences applied ... ABSOLETE!!!)'
      G(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth,:)=0.d0
      l=0
      do i=1,NW
        do j=1,NL
          do k=1,Ssvd-1
            l=l+1
            G(Nseis+Ngps+1+l,(i-1)*NL*Ssvd+(j-1)*Ssvd+(k-1)+1)=+1.d0
            G(Nseis+Ngps+1+l,(i-1)*NL*Ssvd+(j-1)*Ssvd+(k-0)+1)=-1.d0
          enddo
        enddo
      enddo
      do k=1,Ssvd
        do i=1,NW
          do j=1,NL-1
            l=l+1
            G(Nseis+Ngps+1+l,(i-1)*NL*Ssvd+(j-1)*Ssvd+(k-1)+1)=+1.d0
            G(Nseis+Ngps+1+l,(i-1)*NL*Ssvd+(j-0)*Ssvd+(k-1)+1)=-1.d0
          enddo
        enddo
      enddo
      do j=1,NL
        do k=1,Ssvd
          do i=1,NW-1
            l=l+1
            G(Nseis+Ngps+1+l,(i-1)*NL*Ssvd+(j-1)*Ssvd+(k-1)+1)=+1.d0
            G(Nseis+Ngps+1+l,(i-0)*NL*Ssvd+(j-1)*Ssvd+(k-1)+1)=-1.d0
          enddo
        enddo
      enddo
      G(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth,:)=G(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth,:)/norminput*dble(Msvd/2)*abs(smoothkoef)
      D(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth)=0.d0
    elseif(smoothkoef>0.d0)then
      write(*,*)'  (smoothing by model covariance matrix applied)'
      ALLOCATE(CM(Msvd,Msvd))
      CALL fillcovarmatrix()
      G(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth,1:Msvd)=CM(1:Msvd,1:Msvd)
      D(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth)=0.d0
      DEALLOCATE(CM)
    else
      write(*,*)'  (no smoothing)'
    endif

    ! Slip constraint (based on aftershock accurrence in the present version)
    if(slipweight>0.d0)then
      write(*,*)'  (slip constraint based on aftershocks applied)'
      allocate(slipconstraint(NL,NW))
      open(388,FILE='aftershocks-rates.dat')
      do j=1,NW
        read(388,*)(slipconstraint(i,j),i=1,NL)
      enddo
      close(388)
      G(Nseis+Ngps+1+Nsmooth+1:Nseis+Ngps+1+Nsmooth+Nslip,:)=0.d0
      k=0
      do j=1,NW
        do i=1,NL
          do l=1,Ssvd
            k=k+1  !=(j-1)*Ssvd*NL+(i-1)*Ssvd+l
            G(Nseis+Ngps+1+Nsmooth+k,k)=slipconstraint(i,j)/dble(Ssvd)*slipweight
          enddo
        enddo
      enddo
      deallocate(slipconstraint)
      D(Nseis+Ngps+1+Nsmooth+1:Nseis+Ngps+1+Nsmooth+Nslip)=0.d0
    endif
        
    CONTAINS

! NEW - FILTERING AND INTEGRATING IN THE TIME DOMAIN
    SUBROUTINE CreateH()
    USE SISVDmodule
    IMPLICIT NONE
    real*4 dumarr(6),dtr4,f1r4,f2r4
    INTEGER dumi
    COMPLEX*16,DIMENSION(:,:,:),ALLOCATABLE:: cirN,cirE,cirZ
    COMPLEX*16,DIMENSION(:),ALLOCATABLE:: cseis
    REAL*4,DIMENSION(:),ALLOCATABLE:: fltr4
    INTEGER i,j,jj,k,mm,jw,jl

    write(*,*)'Creating matrix H...'

    allocate(cirN(min(nfmax,np),NL*NW,NRseis),cirE(min(nfmax,np),NL*NW,NRseis),cirZ(min(nfmax,np),NL*NW,NRseis))
    open(20,form='unformatted',file='NEZsor.dat')  !if problems with NEZsor.dat file appears, check that resort.f90 writes unformatted file (and not binary)!
    do i=1,NRseis
      j=0
      read(20) dumi
      do jw=1,NW
        do jl=1,NL
          j=j+1
          read(20) dumi
          do k=1,nfmax
            read(20) dumarr
            if(k>np)cycle
            cirN(k,j,i)=cmplx(dumarr(1),dumarr(4))*mu(jl,jw)*elem*dt
            cirE(k,j,i)=cmplx(dumarr(2),dumarr(5))*mu(jl,jw)*elem*dt
            cirZ(k,j,i)=cmplx(dumarr(3),dumarr(6))*mu(jl,jw)*elem*dt
          enddo
        enddo
      enddo
    enddo
    close(20)

    dtr4=dt
    j=0
    write(*,*)'  (Correcting GFs for artifical time delay of',artifDT,'sec.)'
    do jj=1,NRseis
      f1r4=fc2(fcsta(jj));f2r4=fc3(fcsta(jj))
      do k=1,3
        if(stainfo(k,jj)==0)cycle
        j=j+1
!$OMP parallel private(i,mm,cseis,fltr4) DEFAULT(SHARED)
        allocate(cseis(np),fltr4(np))
!$OMP do SCHEDULE(DYNAMIC,1)
        do i=1,NL*NW
          cseis=0.d0
          SELECT CASE (k)
          CASE(1)
            cseis(1:min(nfmax,np))=cirN(1:min(nfmax,np),i,jj)*df
          CASE(2)
            cseis(1:min(nfmax,np))=cirE(1:min(nfmax,np),i,jj)*df
          CASE(3)
            cseis(1:min(nfmax,np))=cirZ(1:min(nfmax,np),i,jj)*df
          END SELECT
          cseis(np/2+2:np)=conjg(cseis(np/2:2:-1))
          call four1(cseis,np,1)
          fltr4=real(cseis)*dt   !dt for time integration
          do mm=1,int(artifDT/dt)  ! REMOVING DWN ARTIFACT FROM THE SEISMOGRAM BEGINNING
            fltr4(mm)=fltr4(mm)*(cos((dt*real(mm-1)-artifDT)/artifDT*PI)/2.+.5);
          enddo
          if(f1r4>0.)then
            CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'BP', f1r4, f2r4, dtr4, 1, np)
          else
            CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'LP', f1r4, f2r4, dtr4, 1, np)
          endif
          do mm=2,np   !time integration
            fltr4(mm)=fltr4(mm)+fltr4(mm-1)
          enddo
          if(i==NL*(NW-1)+1.and.jj==5)then
            do mm=1,np
              write(238,*)dt*(mm-1),fltr4(mm)
            enddo
            write(238,*);write(238,*)
          endif
          do mm=1,nT
            H((j-1)*nT+mm,(i-1)*Ssvd+1:i*Ssvd)=dble(fltr4(iT1-1+mm-iT0:iT1+mm-Ssvd-iT0:-1))
          enddo
        enddo
!$omp end do
        deallocate(cseis,fltr4)
!$omp end parallel
      enddo
    enddo

    deallocate(cirN,cirE,cirZ)
    END SUBROUTINE
    

    SUBROUTINE CreateUvec()
    USE SISVDmodule
    IMPLICIT NONE
    REAL*8 dum
    INTEGER i,j,jj,k,kk,mm

    write(*,*)'Creating vector u...'

    if(syntdata==0)then

      j=0
      do jj=1,NRseis
        do k=1,3
          if(stainfo(k,jj)==0)cycle
          j=j+1
          do i=1,NL*NW
            SELECT CASE (k)
            CASE(1)
              open(290,FILE='rvseisn.dat')
            CASE(2)
              open(290,FILE='rvseise.dat')
            CASE(3)
              open(290,FILE='rvseisz.dat')
            END SELECT
            do kk=1,iT2
              if(kk<iT1)then
                read(290,*)
              else
                read(290,*)(dum,mm=1,jj),Uvec((j-1)*nT+kk-iT1+1)
              endif
            enddo
            close(290)
          enddo
        enddo
      enddo

    else

      Uvec=matmul(H,tfin);   !POZOR, v pripade hlazeni by se melo jeste LP filtrovat, jelikoz Greenovky jsou jen HP filtrovane
      
      open(296,FILE='rvseisnez.dat')
      do i=1,nT
        write(296,'(1000E13.5)')dt*(iT1-1+i-1),(Uvec((j-1)*nT+i),j=1,NSTAcomp)
      enddo
      close(296)

    endif

    END SUBROUTINE


    SUBROUTINE synthetic_model()
    USE SISVDmodule
    IMPLICIT NONE
    REAL*8,ALLOCATABLE:: timefunc(:,:,:),slip1(:,:),slip2(:,:),risetime(:,:),ruptime(:,:)
    REAL*8 posL,posW,time,dum,moment
    INTEGER i,j,k

    allocate(timefunc(Ssvd,NL,NW),slip1(NL,NW),slip2(NL,NW),risetime(NL,NW),ruptime(NL,NW))
    risetime=2.d0
    if(syntdatai==0.and.syntdataj==0)then
      write(*,*)'  (full slip model)'
      slip1=1.d0
      slip2=1.d0
      do j=1,NW
        posW=(dble(j)-.5d0)*widt/dble(NW)
        do i=1,NL
          posL=(dble(i)-.5d0)*leng/dble(NL)
          ruptime(i,j)=sqrt((posL-epicL)**2+(posW-epicW)**2)/vr
! Zpozdena trhlina
          if(i>NL/2.and.j<NW/2)ruptime(i,j)=ruptime(i,j)+3.

!	if(posL>25.d3.and.posL<35.d3.and.posW>=5.d3)then
!	  slip(i,j)=slip(i,j)*3. !*2.
!	elseif(posL>5.d3.and.posL<15.d3.and.posW<=5.d3)then
!	  slip(i,j)=slip(i,j)*3. !;ruptime(i,j)=ruptime(i,j)+2.d0
 !	else
 !	  slip(i,j)=0.d0
 !	endif

!       !hladky skluz (Synteticky model pro clanek GRL):
!        slip1(i,j)=1.d0*exp(-(posL-25.d3)**2/5.d3**2)+1.d0*exp(-(posL-10.d3)**2/5.d3**2)
!        slip2(i,j)=2.d0*exp(-((posL-15.d3)/5.d3)**2)

! Skluz pro synteticke testy do clanku Korelace a Resolmatrix (SVD 1D zlom)
!         slip1(i,j)=1.d0*exp(-(posL-10.d3)**2/5.d3**2)+1.d0*exp(-(posL-30.d3)**2/5.d3**2)

!DVE ASPERITY:
!jednostranne
!         if(    i>=NL/2+1.and.j>=NW/2+1)then
!           slip1(i,j)=exp(-(posL-5.d3)**2/2.5d3**2)+exp(-(posL-15.d3)**2/2.5d3**2)
!           ruptime(i,j)=ruptime(i,j)+2.d0
!         elseif(i<=NL/2+1.and.j<NW/2+1)then
!           slip1(i,j)=exp(-(posL-5.d3)**2/2.5d3**2)+exp(-(posL-15.d3)**2/2.5d3**2)
!         else
!           slip1(i,j)=0.0d0
!         endif
!bilateral:
!         if(j>=NW/4.+1..and.j<=NW/4.*3.+1.)then
           slip1(i,j)=exp(-(posL-5.d3)**2/2.5d3**2)+exp(-(posL-15.d3)**2/2.5d3**2)
!         else
!           slip1(i,j)=0.0d0
!         endif

!Van
!         if(j>=NW/4.+1..and.j<=NW/4.*3.+1.)then
!           slip1(i,j)=exp(-(posL-20.d3)**2/10.d3**2)+exp(-(posL-55.d3)**2/10.d3**2)
!         else
!           slip1(i,j)=0.0d0
!         endif

!POKUS O MODEL MYIAGI:
!         if(    i<=NL/2)then
!           slip1(i,j)=0.d0
!         elseif(i>NL/2.and.j<NW/2)then
!           slip1(i,j)=0.d0
!         elseif(i>NL*7/8.and.j>NW*3/4)then
!           slip1(i,j)=0.d0
!         elseif(i<=NL*7/8.and.j<=NW*3/4)then
!           slip1(i,j)=2.d0*cos(PI*(i-11*NL/16)*8/3*NL)**2
!         else
!           slip1(i,j)=1.d0
!           ruptime(i,j)=sqrt((leng*.75d0-epicL)**2+(widt*.75d0-epicW)**2)/3000.d0+sqrt((posL-leng*.75d0)**2+(posW-widt*.75d0)**2)/2200.d0
!         endif

!         slip1(i,j)=1.d0   !homogenni skluz

        enddo
      enddo
    else
      write(*,*)'  (delta function slip-rate for resolution analysis)'
      slip2=0.d0;slip1=0.d0;
      slip1(syntdatai,syntdataj)=1.d0;ruptime(syntdatai,syntdataj)=5.
    endif

    timefunc=0.d0
    do j=1,NW
      posW=(dble(j)-.5d0)*dW
      do i=1,NL
        posL=(dble(i)-.5d0)*dL
        do k=1,Ssvd
          time=dt*dble(k-1)
          if(time>ruptime(i,j).and.time<ruptime(i,j)+risetime(i,j))then
            timefunc(k,i,j)=(time-ruptime(i,j))*exp(-(time-ruptime(i,j))/1.)*slip1(i,j)
  	      else
            timefunc(k,i,j)=0.d0
          endif
!          if(time>3.d0+ruptime(i,j).and.time<5.d0+ruptime(i,j))then
!            timefunc(k,i,j)=timefunc(k,i,j)+(time-3.d0-ruptime(i,j))*exp(-(time-3.d0-ruptime(i,j))/1.)*slip2(i,j)
!          endif
        enddo

        dum=sum(real(timefunc(:,i,j)))*dt
        if(dum>0.d0)timefunc(:,i,j)=timefunc(:,i,j)/dum*(slip1(i,j)+slip2(i,j))
        if(dum==0.d0)ruptime(i,j)=0.d0
      enddo
    enddo
    moment=0.d0
    do k=1,Ssvd
      moment=moment+sum(timefunc(k,:,:)*mu(:,:))*elem*dt
    enddo
    timefunc=timefunc/moment*Mfix

    write(*,*)'  (saving input model)'
    open(201,FILE='inputtf1d.dat')
    do i=1,Ssvd
      write(201,'(1000E13.5)')(sum(timefunc(i,j,:)),j=1,NL)
    enddo
    close(201)
    open(201,FILE='inputtf.dat')
    write(201,'(1E13.5)')timefunc
    close(201)
    open(296,FILE='inputtfslip2D.dat')
    do j=1,NW
      write(296,'(1000E13.5)')(sum(timefunc(:,i,j))*dt,i=1,NL),sum(timefunc(:,NL,j))*dt
!      write(296,'(1000E13.5)')(sum(timefunc(:,i,j))*dt*mu(i,j)*elem,i=1,NL),sum(timefunc(:,NL,j))*dt*mu(NL,j)*elem
    enddo
    write(296,'(1000E13.5)')(sum(timefunc(:,i,NW))*dt,i=1,NL),sum(timefunc(:,NL,NW))*dt
!    write(296,'(1000E13.5)')(sum(timefunc(:,i,NW))*dt*mu(i,NW)*elem,i=1,NL),sum(timefunc(:,NL,NW))*dt*mu(NL,NW)*elem
    write(296,*);write(296,*)
    do j=1,NW
      write(296,'(1000E13.5)')ruptime(1:NL,j),ruptime(NL,j)
    enddo
    write(296,'(1000E13.5)')ruptime(1:NL,NW),ruptime(NL,NW)
    close(296)

    tfin=RESHAPE(timefunc,(/Msvd/))

    deallocate(timefunc,slip1,slip2,risetime,ruptime)
    END SUBROUTINE


    SUBROUTINE fillcovarmatrix()
    USE SISVDmodule
    IMPLICIT NONE
    INTEGER, PARAMETER:: AA=2 !(AntiAlias for cross-correlation function)
    INTEGER FFT3D(3)
    COMPLEX*16,ALLOCATABLE:: CMtime(:),CMspace(:,:),CMst(:,:,:)
    INTEGER NLFFT,NWFFT,NPFFT
    REAL*8 dkL,dkW,kL,kW
    INTEGER i,j,k,m1,m2,i1,j1,k1,dumi
    REAL*8 sum,freq
    REAL*8 normCM
    
    write(*,*)'  (Preparing the covariance matrix ...)'
    NLFFT=2**(int(log(dble(NL*AA))/log(2.d0)+.99999999d0))
    NWFFT=2**(int(log(dble(NW*AA))/log(2.d0)+.99999999d0))
    if(NRseis>0)then
      NPFFT=np*AA
    else
      NPFFT=1
    endif
    dkL=1.d0/dble(NLFFT)/dL
    dkW=1.d0/dble(NWFFT)/dW
    FFT3D(1)=NPFFT
    FFT3D(2)=NLFFT
    FFT3D(3)=NWFFT
    ALLOCATE(CMtime(NPFFT),CMspace(NLFFT,NWFFT),CMst(NPFFT,NLFFT,NWFFT))

! CMtime
    CMtime(:)=1.d0;write(*,*)'   (k^-2 in space and f^0 in time)'

! CMspace
!    CMspace(:,:)=1.d0   ! spatialy independent
    do i=1,NLFFT
      if(i>NLFFT/2)then
        kL=-dble(NLFFT-i+1)*dkL
      else
        kL=dble(i-1)*dkL
      endif
      do j=1,NWFFT
        if(j>NWFFT/2)then
          kW=-dble(NWFFT-j+1)*dkW
        else
          kW=dble(j-1)*dkW
        endif
        CMspace(i,j)=1.d0/(1.d0+(kL*leng)**2+(kW*widt)**2)
      enddo
    enddo
    
! CMspacetime (CMst)
    do i=1,NLFFT
      do j=1,NWFFT
        CMst(:,i,j)=CMspace(i,j)*CMtime(:)
      enddo
    enddo
    CMst=CMst*conjg(CMst)*df*dkL*dkW
    CALL fourn(CMst(:,:,:),FFT3D(:),3,1)

! Fill-in CM
    m2=0
    do j=1,NW
      do i=1,NL
        do k=1,Ssvd
          m2=m2+1
          m1=0
          do j1=j,1,-1
            do i1=i,1,-1
              do k1=k,1,-1
                m1=m1+1
                CM(m2,m1)=dble(CMst(k1,i1,j1))
              enddo
              do k1=2,Ssvd-k+1
                m1=m1+1
                CM(m2,m1)=dble(CMst(k1,i1,j1))
              enddo
            enddo
            do i1=2,NL-i+1
              do k1=k,1,-1
                m1=m1+1
                CM(m2,m1)=dble(CMst(k1,i1,j1))
              enddo
              do k1=2,Ssvd-k+1
                m1=m1+1
                CM(m2,m1)=dble(CMst(k1,i1,j1))
              enddo
            enddo
          enddo
          do j1=2,NW-j+1
            do i1=i,1,-1
              do k1=k,1,-1
                m1=m1+1
                CM(m2,m1)=dble(CMst(k1,i1,j1))
              enddo
              do k1=2,Ssvd-k+1
                m1=m1+1
                CM(m2,m1)=dble(CMst(k1,i1,j1))
              enddo
            enddo
            do i1=2,NL-i+1
              do k1=k,1,-1
                m1=m1+1
                CM(m2,m1)=dble(CMst(k1,i1,j1))
              enddo
              do k1=2,Ssvd-k+1
                m1=m1+1
                CM(m2,m1)=dble(CMst(k1,i1,j1))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    normCM=0.d0
	do i=1,Msvd
      normCM=normCM+dble(CM(i,i))
    enddo
    normCM=normCM/dble(Msvd)
	CM(1:Msvd,1:Msvd)=CM(1:Msvd,1:Msvd)/normCM*(smoothkoef)**2

!goto 10

    open(111,FILE='CMaprior.dat')
	do i=1,Msvd
	  write(111,*)sqrt(CM(i,i))
	enddo
    close(111)

    open(111,FILE='mtildeslip2D-aprior-stddev.dat')
    write(*,*)sqrt(sum(CM))*dt/dble(Msvd)   !støední odchylka prùmìrného skluzu
    do j=1,NW
      dumi=(j-1)*NL*Ssvd
      write(111,'(1000E13.5)')(sqrt(sum(CM(dumi+(i-1)*Ssvd+1:dumi+i*Ssvd,dumi+(i-1)*Ssvd+1:dumi+i*Ssvd)))*dt,i=1,NL),sqrt(sum(CM(dumi+(NL-1)*Ssvd+1:dumi+NL*Ssvd,dumi+(NL-1)*Ssvd+1:dumi+NL*Ssvd)))*dt  !støední odchylka skluzu
    enddo
    dumi=(NW-1)*NL*Ssvd
    write(111,'(1000E13.5)')(sqrt(sum(CM(dumi+(i-1)*Ssvd+1:dumi+i*Ssvd,dumi+(i-1)*Ssvd+1:dumi+i*Ssvd)))*dt,i=1,NL),sqrt(sum(CM(dumi+(NL-1)*Ssvd+1:dumi+NL*Ssvd,dumi+(NL-1)*Ssvd+1:dumi+NL*Ssvd)))*dt  !støední odchylka skluzu
    close(111)

! Inversion of CMs    
#ifdef MKL
    call dpotrf('U',Msvd,CM,Msvd,i)   ! Upper triangle of CM becomes U
    if (i>0) then
      print *,' the matrix is not positive definite!'
      stop
    endif
    call dtrtri('U','N',Msvd,CM,Msvd,i)  ! Inverse of U
    do i=1,Msvd-1                        ! Transpose and clear upper triangle
      CM(i+1:Msvd,i)=CM(i,i+1:Msvd)
      CM(i,i+1:Msvd)=0.d0
    enddo
#else
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!Not yet coded!!!!!!!!!!!!!!!!!!!!!!!'
!    CMinv=CM
!    CALL cholsl(Msvd,CMinv,CM)
#endif

!    open(333,FILE='CM.dat')
!    do i=1,Msvd
!      write(333,'(10000E13.5)')CM(i,:)
!    enddo
!    close(333)
10  write(*,*)'  (Covariance matrix ready)'

    DEALLOCATE(CMtime,CMspace,CMst)
    END SUBROUTINE

    END



    SUBROUTINE filtconstruct(flt,N,flo,fleft,fright,fro,dt)
    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535d0
    INTEGER N,i
    COMPLEX*16 flt(N)
    REAL*8 flo,fleft,fright,fro,f2,dt,df
    REAL*8 freq
    REAL*4 fltr4(N),f1r4,f2r4,dtr4
    fltr4=0.;fltr4(1)=1.
    dtr4=dt
    f1r4=fleft;f2r4=fright
    CALL XAPIIR(fltr4, N, 'BU', 0.0, 0.0, 4,'BP', f1r4, f2r4, dtr4, 1, N)
    flt=fltr4
    do i=1,N
      write(321,*)dt*(i-1),fltr4(i)
    enddo
    CALL four1(flt,N,-1)
goto 10   !comment for double-cos filter
    df=1.d0/(dt*dble(N))
    flt=0.d0
    do i=1,N/2+1
      freq=df*(i-1)
      if(freq>=fleft.and.freq<=fright)then
        flt(i)=1.d0
      elseif(freq>=flo.and.freq<=fleft.and.flo/=fleft)then
        flt(i)=(.5+.5*cos(PI*(freq-fleft)/(flo-fleft)))
      elseif(freq>=fright.and.freq<=fro.and.fro/=fright)then
        flt(i)=(.5+.5*cos(PI*(fright-freq)/(fright-fro)))
      else
        cycle
      endif
    enddo
10  END


    SUBROUTINE filtconstructHP(flt,N,flo,fleft,fright,fro,dt)
    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535d0
    INTEGER N,i
    COMPLEX*16 flt(N)
    REAL*8 flo,fleft,fright,fro,f2,dt,df
    REAL*8 freq
    REAL*4 fltr4(N),f1r4,f2r4,dtr4
    fltr4=0.;fltr4(1)=1.
    dtr4=dt
    f1r4=fleft;f2r4=fright
    CALL XAPIIR(fltr4, N, 'BU', 0.0, 0.0, 4,'HP', f1r4, f2r4, dtr4, 1, N)
    flt=fltr4
    do i=1,N
      write(321,*)dt*(i-1),fltr4(i)
    enddo
    CALL four1(flt,N,-1)
    END


    SUBROUTINE filtconstructLP(flt,N,flo,fleft,fright,fro,dt)
    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535d0
    INTEGER N,i
    COMPLEX*16 flt(N)
    REAL*8 flo,fleft,fright,fro,f2,dt,df
    REAL*8 freq
    REAL*4 fltr4(N),f1r4,f2r4,dtr4
    fltr4=0.;fltr4(1)=1.
    dtr4=dt
    f1r4=fleft;f2r4=fright
    CALL XAPIIR(fltr4, N, 'BU', 0.0, 0.0, 4,'LP', f1r4, f2r4, dtr4, 1, N)
    flt=fltr4
    do i=1,N
      write(321,*)dt*(i-1),fltr4(i)
    enddo
    CALL four1(flt,N,-1)
goto 10   !comment for double-cos filter
    df=1.d0/(dt*dble(N))
    flt=0.d0
    do i=1,N/2+1
      freq=df*(i-1)
      if(freq>=fleft.and.freq<=fright)then
        flt(i)=1.d0
      elseif(freq>=flo.and.freq<=fleft.and.flo/=fleft)then
        flt(i)=(.5+.5*cos(PI*(freq-fleft)/(flo-fleft)))
      elseif(freq>=fright.and.freq<=fro.and.fro/=fright)then
        flt(i)=(.5+.5*cos(PI*(fright-freq)/(fright-fro)))
      else
        cycle
      endif
    enddo
10  END



Subroutine cholsl(n,A,aa)
  integer n
  real*8 A(0:n-1,0:n-1), aa(0:n-1,0:n-1)
  integer i,j,k

  call choldcsl(n,A,aa)

  do i = 0, n-1
    do j = i + 1, n-1
	  aa(i,j) = 0.d0
    end do
  end do

  do i = 0, n-1
    aa(i,i) = aa(i,i) * aa(i,i)
    do k = i + 1, n-1
      aa(i,i) = aa(i,i) + aa(k,i) * aa(k,i)
    end do

    do j = i + 1, n-1
      do k = j, n-1
        aa(i,j) = aa(i,j) + aa(k,i) * aa(k,j)
      end do
    end do
  end do
  do i = 0,  n-1
    do j = 0, i-1
      aa(i,j) = aa(j,i)
    end do
  end do
  return
End

Subroutine choldcsl(n,A,aa)
  integer n
  real*8 A(0:n-1,0:n-1), aa(0:n-1,0:n-1)
  integer i,j,k, ialloc
  real*8 sum
  real*8, pointer :: p(:)
  allocate(p(0:n-1),stat=ialloc)

  aa = A

  call choldc1(n, aa, p)

  do i = 0, n-1
    aa(i,i) = 1.d0 / p(i)
    do j = i + 1, n-1
      sum = 0.d0
      do k = i, j-1
	    sum = sum - aa(j,k) * aa(k,i)
      end do
      aa(j,i) = sum / p(j)
    end do
  end do
  deallocate(p)
  return
End

Subroutine choldc1(n,a,p)
  integer n
  real*8 a(0:n-1,0:n-1), p(0:n-1)
  integer i,j,k
  real*8 sum
  do i = 0, n-1
    do j = i, n-1
      sum = a(i,j)
      do k = i - 1, 0, -1
        sum = sum - a(i,k) * a(j,k)
      end do
      if (i.eq.j) then
        if (sum <= 0.d0) then
          write(*,*)sum
          print *,' the matrix is not positive definite!'
        endif
        p(i) = dsqrt(sum)
      else
        a(j,i) = sum / p(i)
      end if
	end do
  end do
  return
End

