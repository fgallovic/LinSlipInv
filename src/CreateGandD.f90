!   AUTHOR:  Frantisek Gallovic
    
    SUBROUTINE CreateGandD()
    USE SISVDmodule
    IMPLICIT NONE
    INTEGER i,j,k,kk,l,IRET,SegShift,Msvdkk
    REAL*8 dum
    REAL*8,ALLOCATABLE:: H(:,:),Uvec(:),tfin(:),CM(:,:),CD(:,:),Hnew(:,:)
    REAL*4,ALLOCATABLE,DIMENSION(:,:):: stalocGPS,dataGPS
    REAL*4 gpsgfN,gpsgfE,gpsgfZ
    REAL*4 ALPHA,dumGPS,dLgps,dWgps
    REAL*8,ALLOCATABLE:: slipconstraint(:,:,:)

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
    endif

    if(NRseis>0)then
      allocate(H(Nseis,Msvd),Uvec(Nseis))
      CALL CreateH()
      CALL CreateUvec()
      if(compweights==2)then
        write(*,*)'  (Applying distance-dependent CD - still in testing mode!)'
        allocate(CD(Nseis,Nseis))
        CALL fillCDmatrix(Uvec,CD)   !Following Hallo and Gallovic (2016)
#ifdef MKL
        allocate(Hnew(Nseis,Msvd))
        call dgemm('N','N',Nseis,Msvd,Nseis,1.d0,CD,Nseis,H,Nseis,0.d0,Hnew,Nseis)
        H=Hnew
        deallocate(Hnew)
#else
        allocate(Hnew(Nseis,Msvd))
        Hnew=matmul(CD,H)
        H=Hnew
        deallocate(Hnew)
#endif
        Uvec=matmul(CD,Uvec)    !POZOR, TATO STANDARDIZACE SE UZ NEKORIGUJE, TAKZE POROVNANI SEISMOGRAMU JE STANDARDIZOVANE
        deallocate(CD)
        smoothkoefGF=1.   !POZOR - VYPINAM TIM TEN PREPINAC A ZUSTAVA JEN VLIV TECH CASOVYCH POSUNU
      endif
      
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
        G(1:Nseis,:)=G(1:Nseis,:)/sum(1.d0/normdat(:))*sum(staweight(:,:)*dble(stainfo(:,:)))/smoothkoefGF
        D(1:Nseis)=D(1:Nseis)/sum(1.d0/normdat(:))*sum(staweight(:,:)*dble(stainfo(:,:)))/smoothkoefGF
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
      do kk=1,NSeg
        dLgps=real(dL(kk))
        dWgps=real(dW(kk))
        SegShift=sum(NW(1:kk-1)*NL(1:kk-1))*Ssvd
        do k=1,NRgps
          do j=1,NW(kk)
            do i=1,NL(kk)
              ALPHA=(lambda(i,j,kk)+mu(i,j,kk))/(lambda(i,j,kk)+2.*mu(i,j,kk))
              CALL DC3Dmodif(ALPHA,sourNgps(i,j,kk),sourEgps(i,j,kk),sourZgps(i,j,kk),strikeGPS(i,j,kk),dipGPS(i,j,kk),rakeGPS(i,j,kk),dLgps,dWgps,dumGPS,stalocGPS(1,k),stalocGPS(2,k),stalocGPS(3,k),gpsgfN,gpsgfE,gpsgfZ,IRET)
              G(Nseis+(k-1)*3+1,SegShift+(j-1)*Ssvd*NL(kk)+(i-1)*Ssvd+1:SegShift+(j-1)*Ssvd*NL(kk)+i*Ssvd)=gpsgfN*dt/sigmaGPS(1,k)
              G(Nseis+(k-1)*3+2,SegShift+(j-1)*Ssvd*NL(kk)+(i-1)*Ssvd+1:SegShift+(j-1)*Ssvd*NL(kk)+i*Ssvd)=gpsgfE*dt/sigmaGPS(2,k)
              G(Nseis+(k-1)*3+3,SegShift+(j-1)*Ssvd*NL(kk)+(i-1)*Ssvd+1:SegShift+(j-1)*Ssvd*NL(kk)+i*Ssvd)=gpsgfZ*dt/sigmaGPS(3,k)
            enddo
          enddo
          D(Nseis+(k-1)*3+1:Nseis+(k-1)*3+3)=dataGPS(1:3,k)/sigmaGPS(1:3,k)
        enddo
      enddo
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
      do kk=1,NSeg
        SegShift=sum(NW(1:kk-1)*NL(1:kk-1))*Ssvd
        do i=1,NW(kk)
          do j=1,NL(kk)
            do k=1,Ssvd-1
              l=l+1
              G(Nseis+Ngps+1+l,SegShift+(i-1)*NL(kk)*Ssvd+(j-1)*Ssvd+(k-1)+1)=+1.d0
              G(Nseis+Ngps+1+l,SegShift+(i-1)*NL(kk)*Ssvd+(j-1)*Ssvd+(k-0)+1)=-1.d0
            enddo
          enddo
        enddo
        do k=1,Ssvd
          do i=1,NW(kk)
            do j=1,NL(kk)-1
              l=l+1
              G(Nseis+Ngps+1+l,SegShift+(i-1)*NL(kk)*Ssvd+(j-1)*Ssvd+(k-1)+1)=+1.d0
              G(Nseis+Ngps+1+l,SegShift+(i-1)*NL(kk)*Ssvd+(j-0)*Ssvd+(k-1)+1)=-1.d0
            enddo
          enddo
        enddo
        do j=1,NL(kk)
          do k=1,Ssvd
            do i=1,NW(kk)-1
              l=l+1
              G(Nseis+Ngps+1+l,SegShift+(i-1)*NL(kk)*Ssvd+(j-1)*Ssvd+(k-1)+1)=+1.d0
              G(Nseis+Ngps+1+l,SegShift+(i-0)*NL(kk)*Ssvd+(j-1)*Ssvd+(k-1)+1)=-1.d0
            enddo
          enddo
        enddo
        G(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth,:)=G(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth,:)/norminput*dble(Msvd/2)*abs(smoothkoef)
        D(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth)=0.d0
      enddo
    elseif(smoothkoef>0.d0)then
      write(*,*)'  (smoothing by model covariance matrix applied)'
!      G(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth,1:Msvd)=0.d0
      do kk=1,NSeg
        Msvdkk=NL(kk)*NW(kk)*Ssvd
        ALLOCATE(CM(Msvdkk,Msvdkk))
        CALL fillcovarmatrix(kk)
        SegShift=sum(NW(1:kk-1)*NL(1:kk-1))*Ssvd
        G(SegShift+Nseis+Ngps+2:SegShift+Nseis+Ngps+1+Msvdkk,SegShift+1:SegShift+Msvdkk)=CM(1:Msvdkk,1:Msvdkk)
        DEALLOCATE(CM)
      enddo
      D(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth)=0.d0
    else
      write(*,*)'  (no smoothing)'
    endif

    ! Slip constraint (based on aftershock accurrence in the present version)
    if(slipweight>0.d0)then
      write(*,*)'  (slip constraint based on aftershocks applied)'
      allocate(slipconstraint(maxval(NL),maxval(NW),NSeg))
      open(388,FILE='aftershocks-rates.dat')
      do kk=1,NSeg
        do j=1,NW(kk)
          read(388,*)(slipconstraint(i,j,kk),i=1,NL(kk))
        enddo
      enddo
      close(388)
      G(Nseis+Ngps+1+Nsmooth+1:Nseis+Ngps+1+Nsmooth+Nslip,:)=0.d0
      k=0
      do kk=1,NSeg
        do j=1,NW(kk)
          do i=1,NL(kk)
            do l=1,Ssvd
              k=k+1  !=(j-1)*Ssvd*NL+(i-1)*Ssvd+l
              G(Nseis+Ngps+1+Nsmooth+k,k)=slipconstraint(i,j,kk)/dble(Ssvd)*slipweight
            enddo
          enddo
        enddo
      enddo
      deallocate(slipconstraint)
      D(Nseis+Ngps+1+Nsmooth+1:Nseis+Ngps+1+Nsmooth+Nslip)=0.d0
    endif
        
    CONTAINS

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
    open(20,form='unformatted',file='NEZsor.dat')  !if problems with NEZsor.dat file appears, check that resort.f90 writes unformatted file (and not binary)!
    write(*,*)'  (Correcting GFs for artifical time delay by',artifDT,'sec.)'

    do kk=1,NSeg
      allocate(cirN(min(nfmax,np),NL(kk)*NW(kk),NRseis),cirE(min(nfmax,np),NL(kk)*NW(kk),NRseis),cirZ(min(nfmax,np),NL(kk)*NW(kk),NRseis))
      SegShift=sum(NW(1:kk-1)*NL(1:kk-1))*Ssvd
      do i=1,NRseis
        j=0
        read(20) dumi
        do jw=1,NW(kk)
          do jl=1,NL(kk)
            j=j+1
            read(20) dumi
            do k=1,nfmax
              read(20) dumarr
              if(k>np)cycle
              cirN(k,j,i)=cmplx(dumarr(1),dumarr(4))*mu(jl,jw,kk)*elem(kk)*dt
              cirE(k,j,i)=cmplx(dumarr(2),dumarr(5))*mu(jl,jw,kk)*elem(kk)*dt
              cirZ(k,j,i)=cmplx(dumarr(3),dumarr(6))*mu(jl,jw,kk)*elem(kk)*dt
            enddo
          enddo
        enddo
      enddo

      dtr4=dt
      j=0
      do jj=1,NRseis
        f1r4=fc1(fcsta(jj));f2r4=fc2(fcsta(jj))
        do k=1,3
          if(stainfo(k,jj)==0)cycle
          j=j+1
!$OMP parallel private(i,mm,cseis,fltr4) DEFAULT(SHARED)
          allocate(cseis(np),fltr4(np))
!$OMP do SCHEDULE(DYNAMIC,1)
          do i=1,NL(kk)*NW(kk)
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
            if(kk==1.and.i==NL(kk)*(NW(kk)-1)+1.and.jj==5)then
              do mm=1,np
                write(238,*)dt*(mm-1),fltr4(mm)
              enddo
              write(238,*);write(238,*)
            endif
            do mm=1,nT
              H((j-1)*nT+mm,SegShift+(i-1)*Ssvd+1:SegShift+i*Ssvd)=dble(fltr4(iT1-1+mm-iT0:iT1+mm-Ssvd-iT0:-1))
            enddo
          enddo
!$omp end do
          deallocate(cseis,fltr4)
!$omp end parallel
        enddo
      enddo

      deallocate(cirN,cirE,cirZ)
    enddo

    close(20)
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

    else

      Uvec=matmul(H,tfin);
      
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
    REAL*8,ALLOCATABLE:: timefunc(:,:,:,:)
    REAL*8,ALLOCATABLE,DIMENSION(:,:,:):: slip1,slip2,risetime,ruptime
    REAL*8 posL,posW,time,dum,moment
    INTEGER i,j,k,kk

    allocate(timefunc(Ssvd,maxval(NL),maxval(NW),NSeg),slip1(maxval(NL),maxval(NW),NSeg),slip2(maxval(NL),maxval(NW),NSeg),risetime(maxval(NL),maxval(NW),NSeg),ruptime(maxval(NL),maxval(NW),NSeg))
    risetime=2.d0
    if(syntdatai==0.and.syntdataj==0)then
      write(*,*)'  (full slip model)'
      slip1=1.d0
      slip2=1.d0
      do kk=1,NSeg
        do j=1,NW(kk)
          posW=(dble(j)-.5d0)*widt(kk)/dble(NW(kk))
          do i=1,NL(kk)
            posL=(dble(i)-.5d0)*leng(kk)/dble(NL(kk))
            ruptime(i,j,kk)=sqrt((posL-epicL(kk))**2+(posW-epicW(kk))**2)/vr
! Zpozdena trhlina
            if(kk==1.and.i>NL(1)/2.and.j<NW(1)/2)ruptime(i,j,kk)=ruptime(i,j,kk)+3.

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
             slip1(i,j,kk)=exp(-(posL-5.d3)**2/2.5d3**2)+exp(-(posL-15.d3)**2/2.5d3**2)
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
      enddo
    else
      write(*,*)'  (delta function slip-rate for resolution analysis at segment 1)'
      slip2=0.d0;slip1=0.d0;
      slip1(syntdatai,syntdataj,1)=1.d0;ruptime(syntdatai,syntdataj,1)=5.
    endif

    timefunc=0.d0
    do kk=1,NSeg
      do j=1,NW(kk)
        posW=(dble(j)-.5d0)*dW(kk)
        do i=1,NL(kk)
          posL=(dble(i)-.5d0)*dL(kk)
          do k=1,Ssvd
            time=dt*dble(k-1)
            if(time>ruptime(i,j,kk).and.time<ruptime(i,j,kk)+risetime(i,j,kk))then
              timefunc(k,i,j,kk)=(time-ruptime(i,j,kk))*exp(-(time-ruptime(i,j,kk))/1.)*slip1(i,j,kk)
  	        else
              timefunc(k,i,j,kk)=0.d0
            endif
!            if(time>3.d0+ruptime(i,j).and.time<5.d0+ruptime(i,j))then
!              timefunc(k,i,j)=timefunc(k,i,j)+(time-3.d0-ruptime(i,j))*exp(-(time-3.d0-ruptime(i,j))/1.)*slip2(i,j)
!            endif
          enddo
          if(Ssvd==1)timefunc(1,i,j,kk)=slip1(i,j,kk)    !GPS inversion only
          dum=sum(real(timefunc(:,i,j,kk)))*dt
          if(dum>0.d0)timefunc(:,i,j,kk)=timefunc(:,i,j,kk)/dum*(slip1(i,j,kk)+slip2(i,j,kk))
          if(dum==0.d0)ruptime(i,j,kk)=0.d0
        enddo
      enddo
    enddo
    moment=0.d0
    do kk=1,NSeg
      do k=1,Ssvd
        moment=moment+sum(timefunc(k,1:NL(kk),1:NW(kk),kk)*mu(1:NL(kk),1:NW(kk),kk))*elem(kk)*dt
      enddo
    enddo
    timefunc=timefunc/moment*Mfix

    write(*,*)'  (saving input model)'
    open(201,FILE='inputtf1d.dat')
    do kk=1,NSeg
      do i=1,Ssvd
        write(201,'(1000E13.5)')(sum(timefunc(i,j,1:NW(kk),kk))*dW(kk),j=1,NL(kk))
      enddo
      write(201,*);write(201,*)
    enddo
    close(201)
    open(201,FILE='inputtf.dat')
    do kk=1,NSeg
      write(201,'(1E13.5)')timefunc(:,1:NL(kk),1:NW(kk),kk)
    enddo
    close(201)
    open(296,FILE='inputtfslip2D.dat')
    do kk=1,NSeg
      do j=1,NW(kk)
        write(296,'(1000E13.5)')(sum(timefunc(:,i,j,kk))*dt,i=1,NL(kk)),sum(timefunc(:,NL(kk),j,kk))*dt
!       write(296,'(1000E13.5)')(sum(timefunc(:,i,j,kk))*dt*mu(i,j,kk)*elem(kk),i=1,NL(kk)),sum(timefunc(:,NL(kk),j,kk))*dt*mu(NL(kk),j,kk)*elem(kk)     !Output is scalar moment instead of slip
      enddo
      write(296,'(1000E13.5)')(sum(timefunc(:,i,NW(kk),kk))*dt,i=1,NL(kk)),sum(timefunc(:,NL(kk),NW(kk),kk))*dt
!      write(296,'(1000E13.5)')(sum(timefunc(:,i,NW(kk),kk))*dt*mu(i,NW,kk)*elem(kk),i=1,NL(kk)),sum(timefunc(:,NL(kk),NW(kk),kk))*dt*mu(NL(kk),NW(kk),kk)*elem(kk)   !Output is scalar moment instead of slip
      write(201,*);write(201,*)
    enddo
    write(296,*);write(296,*)
    do kk=1,NSeg
      do j=1,NW(kk)
        write(296,'(1000E13.5)')ruptime(1:NL(kk),j,kk),ruptime(NL(kk),j,kk)
      enddo
      write(296,'(1000E13.5)')ruptime(1:NL(kk),NW(kk),kk),ruptime(NL(kk),NW(kk),kk)
      write(296,*);write(296,*)
    enddo
    close(296)

!    tfin=RESHAPE(timefunc,(/Msvd/))
    do kk=1,NSeg
      SegShift=sum(NW(1:kk-1)*NL(1:kk-1))*Ssvd
      tfin(SegShift+1:SegShift+NW(kk)*NL(kk)*Ssvd)=reshape(timefunc(:,1:NL(kk),1:NW(kk),kk),(/NL(kk)*NW(kk)*Ssvd/))
    enddo
    
    deallocate(timefunc,slip1,slip2,risetime,ruptime)
    END SUBROUTINE

    
    SUBROUTINE fillCDmatrix(Uvec,CD)
    USE SISVDmodule
    IMPLICIT NONE
    REAL*8 Uvec(Nseis),CD(Nseis,Nseis)
    REAL*8 dvec(nT),CDsub(nT,nT),maxdiagCD
    REAL*8 staN,staE
    INTEGER i,j,jj,k,m,L1
    CD=0.
    j=0
!    L1=5.   !Nahodne casove posuny
    do jj=1,NRseis
      L1=int(max(1.,R(jj)/5000.)/dt)
      write(*,*)R(jj),L1
      do k=1,3
        if(stainfo(k,jj)==0)cycle
        j=j+1
        dvec(1:nT)=Uvec((j-1)*nT+1:j*nT)
        CALL xcov2(nT,dvec,dvec,L1,0,CDsub)
        CD((j-1)*nT+1:j*nT,(j-1)*nT+1:j*nT)=CDsub(1:nT,1:nT)
      enddo
    enddo
    maxdiagCD=0.
    do i=1,Nseis
      if(CD(i,i)>maxdiagCD)maxdiagCD=CD(i,i)
    enddo

    do i=1,Nseis
      CD(i,i)=CD(i,i)+maxdiagCD/100.
      write(111,'(100000E13.5)')CD(i,i) !CD(i,:)
    enddo
    do i=1,nT
       write(155,'(1000E13.5)')CD(1:nT,i)
       write(156,*)Uvec(i)
    enddo

    ! Inversion of CD    
#ifdef MKL
    call dpotrf('U',Nseis,CD,Nseis,i)   ! Upper triangle of CD becomes U
    if (i>0) then
      print *,' the matrix is not positive definite!'
      stop
    endif
    call dpotri('U',Nseis,CD,Nseis,i)  ! Inverse of CD
    call dpotrf('U',Nseis,CD,Nseis,i)   ! Upper triangle of CD becomes U
    do i=1,Nseis-1                        ! Transpose and clear upper triangle
!      CD(i+1:Nseis,i)=CD(i,i+1:Nseis)
!      CD(i,i+1:Nseis)=0.d0
      CD(i+1:Nseis,i)=0.d0 ! Clear lower triangle
    enddo
    continue
#else
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!Not yet coded!!!!!!!!!!!!!!!!!!!!!!!'
    stop
!    CDinv=CD
!    CALL cholsl(Nseis,CDinv,CD)
#endif
    END

    
    SUBROUTINE xcov2(N,f1,f2,L1,L12,C)   !just xcovariance for L12=0 !!!
    IMPLICIT NONE
    INTEGER N,L1,L12
    REAL*8 f1(N),f2(N),C(N,N)
    COMPLEX*16,ALLOCATABLE:: smooth1(:),smooth12(:),s1(:),s2(:),s1smooth1(:),s2smooth1(:),s2smooth12(:),XC(:,:)
    INTEGER i,j,FFT
    C=0.
    FFT=2**(int(log(dble(2*N))/log(2.d0)+.99999999d0)+1)
    ALLOCATE(smooth1(FFT),s1(FFT),s2(FFT),s1smooth1(FFT),s2smooth1(FFT),s2smooth12(FFT),XC(FFT,FFT))

    s1=0.;s1(1:N)=f1(1:N)
    s2=0.;s2(1:N)=f2(1:N)
    if(L12>1)then
      allocate(smooth12(FFT))
      smooth12=0.
      smooth12(1:int(dble(L12)/2.+.51))=1./dble(L12*FFT)
      smooth12(FFT:FFT-L12/2+1:-1)=1./dble(L12*FFT)
      CALL four1(smooth12,FFT,1)
      CALL four1(s2,FFT,1)
      s2smooth12=s2*smooth12
      CALL four1(s2smooth12,FFT,-1)
      deallocate(smooth12)
    else
      s2smooth12=s2
    endif
    smooth1=0.
    smooth1(1:int(dble(L1)/2.+.51))=1./dble(L1*FFT)
    smooth1(FFT:FFT-L1/2+1:-1)=1./dble(L1*FFT)
    CALL four1(smooth1,FFT,1)

    do i=1,FFT
      XC(:,i)=s1(:)*cshift(s2smooth12(:),SHIFT=(FFT/2+1)-i)
      CALL four1(XC(:,i),FFT,1)
      XC(:,i)=XC(:,i)*smooth1(:)
      CALL four1(XC(:,i),FFT,-1)
    enddo

    CALL four1(s1,FFT,1)
    s1smooth1=s1*smooth1
    CALL four1(s1smooth1,FFT,-1)
    CALL four1(s2smooth12,FFT,1)
    s2smooth1=s2smooth12*smooth1
    CALL four1(s2smooth1,FFT,-1)
    
    do i=1,FFT
      XC(:,i)=XC(:,i)-s1smooth1(:)*cshift(s2smooth1(:),SHIFT=(FFT/2+1)-i)
    enddo

    do i=1,N
     do j=1,N
       C(i,j)=dble(XC(i,FFT/2+1-(j-i)))
     enddo
    enddo

    DEALLOCATE(smooth1,s1,s2,s2smooth12,s1smooth1,s2smooth1,XC)
    END



    
    SUBROUTINE fillcovarmatrix(kk)
    USE SISVDmodule
    IMPLICIT NONE
    INTEGER, PARAMETER:: AA=2 !(Antialias for cross-correlation function)
    INTEGER FFT3D(3)
    COMPLEX*16,ALLOCATABLE:: CMtime(:),CMspace(:,:),CMst(:,:,:)
    INTEGER NLFFT,NWFFT,NPFFT,Msvdkk
    REAL*8 dkL,dkW,kL,kW
    INTEGER i,j,k,kk,m1,m2,i1,j1,k1,dumi
    REAL*8 sum,freq
    REAL*8 normCM
    
    write(*,*)'  (Preparing the covariance matrix ...)'
    NLFFT=2**(int(log(dble(NL(kk)*AA))/log(2.d0)+.99999999d0))
    NWFFT=2**(int(log(dble(NW(kk)*AA))/log(2.d0)+.99999999d0))
    if(NRseis>0)then
      NPFFT=np*AA
    else
      NPFFT=1
    endif
    dkL=1.d0/dble(NLFFT)/dL(kk)
    dkW=1.d0/dble(NWFFT)/dW(kk)
    FFT3D(1)=NPFFT
    FFT3D(2)=NLFFT
    FFT3D(3)=NWFFT
    ALLOCATE(CMtime(NPFFT),CMspace(NLFFT,NWFFT),CMst(NPFFT,NLFFT,NWFFT))
    Msvdkk=NL(kk)*NW(kk)*Ssvd

! CMtime
    CMtime(:)=1.d0;write(*,*)'    (k^-2 in space and f^0 in time)'

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
        CMspace(i,j)=1.d0/(1.d0+(kL*leng(kk))**2+(kW*widt(kk))**2)
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
    do j=1,NW(kk)
      do i=1,NL(kk)
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
            do i1=2,NL(kk)-i+1
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
          do j1=2,NW(kk)-j+1
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
            do i1=2,NL(kk)-i+1
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
	  do i=1,Msvdkk
      normCM=normCM+dble(CM(i,i))
    enddo
    normCM=normCM/dble(Msvdkk)
  	CM(:,:)=CM(:,:)/normCM*(smoothkoef)**2

!goto 10

    open(111,FILE='CMaprior.dat')   ! Last segment overwrites the others
	  do i=1,Msvdkk
	    write(111,*)sqrt(CM(i,i))
	  enddo
    close(111)

    open(111,FILE='mtildeslip2D-aprior-stddev.dat')   ! Last segment overwrites the others
    do j=1,NW(kk)
      dumi=(j-1)*NL(kk)*Ssvd
      write(111,'(1000E13.5)')(sqrt(sum(CM(dumi+(i-1)*Ssvd+1:dumi+i*Ssvd,dumi+(i-1)*Ssvd+1:dumi+i*Ssvd)))*dt,i=1,NL(kk)),sqrt(sum(CM(dumi+(NL(kk)-1)*Ssvd+1:dumi+NL(kk)*Ssvd,dumi+(NL(kk)-1)*Ssvd+1:dumi+NL(kk)*Ssvd)))*dt  !støední odchylka skluzu
    enddo
    dumi=(NW(kk)-1)*NL(kk)*Ssvd
    write(111,'(1000E13.5)')(sqrt(sum(CM(dumi+(i-1)*Ssvd+1:dumi+i*Ssvd,dumi+(i-1)*Ssvd+1:dumi+i*Ssvd)))*dt,i=1,NL(kk)),sqrt(sum(CM(dumi+(NL(kk)-1)*Ssvd+1:dumi+NL(kk)*Ssvd,dumi+(NL(kk)-1)*Ssvd+1:dumi+NL(kk)*Ssvd)))*dt  !støední odchylka skluzu
    close(111)

! Inversion of CMs    
#ifdef MKL
    call dpotrf('U',Msvdkk,CM,Msvdkk,i)   ! Upper triangle of CM becomes U
    if (i>0) then
      print *,' the matrix is not positive definite!'
      stop
    endif
    call dtrtri('U','N',Msvdkk,CM,Msvdkk,i)  ! Inverse of U
    do i=1,Msvdkk-1                        ! Transpose and clear upper triangle
      CM(i+1:Msvdkk,i)=CM(i,i+1:Msvdkk)
      CM(i,i+1:Msvdkk)=0.d0
    enddo
#else
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!Not yet coded!!!!!!!!!!!!!!!!!!!!!!!'
    stop
!    CMinv=CM
!    CALL cholsl(Msvdkk,CMinv,CM)
#endif

!    open(333,FILE='CM.dat')
!    do i=1,Msvdkk
!      write(333,'(10000E13.5)')CM(i,:)
!    enddo
!    close(333)
10  write(*,*)'  (Covariance matrix ready)'

    DEALLOCATE(CMtime,CMspace,CMst)
    END SUBROUTINE

    END




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

