!   CONVERTS DATA TO SEISMOGRAM FORMAT USED IN LinSlipInv
!   AUTHOR:  Frantisek Gallovic

    IMPLICIT NONE
    REAL,PARAMETER:: PI=3.1415926535
    INTEGER originY,originMONTH,originD,originH,originM,STA,FILTERS,NN,NNout
    REAL originS,DTout
    character*4 staname
    character*20 filename
    INTEGER staY,staMONTH,staD,staH,staM
    real staS,staLAT,staLON,dt,artftimeshift
    integer Istart,FLT,INTEGRATE,detrend
    real,dimension(:),allocatable:: f1,f2
    real,dimension(:),allocatable:: seisn,seise,seisz,tt
    real,dimension(:,:),allocatable:: seisnout,seiseout,seiszout
    real dum,dumN,dumE,dumZ,UNDER,divide
    integer i,j,k,jstart,id,sampl
    REAL aN,bN,aE,bE,aZ,bZ

    open(100,FILE="processseis.in")
    read(100,*)
    read(100,*)originY,originMONTH,originD,originH,originM,originS
    read(100,*)
    read(100,*)STA
    read(100,*)
    read(100,*)DTout,NNout,artftimeshift
    allocate(seisnout(NNout,STA),seiseout(NNout,STA),seiszout(NNout,STA))
    read(100,*)
    read(100,*)FILTERS
    allocate(f1(FILTERS),f2(FILTERS))
    do i=1,FILTERS
      read(100,*)f1(i),f2(i)
    enddo
    read(100,*)
    read(100,*)
    read(100,*)
    
    open(101,FILE='stations.txt')
    do i=1,STA
      read(100,*)staLAT,staLON,staY,staMONTH,staD,staH,staM,staS,dt,divide,detrend,FLT,INTEGRATE,filename
      Istart=int(((staD-originD)*86400.+(staH-originH)*3600.+(staM-originM)*60.+staS-originS+artftimeshift)/dt)
      write(staname,'(A4)')filename
      write(*,*)trim(staname),Istart
      write(101,*)staLAT,staLON,staname
      UNDER=dtout/dt  !undersampling
      NN=int(UNDER*NNout)+1
      allocate(seisn(NN),seise(NN),seisz(NN),tt(NN))
      seisn=0.
      seise=0.
      seisz=0.
      do k=1,NN
        tt(k)=dt*(k-1)
      enddo
      open(201,FILE=trim(filename))
      j=Istart
      if(j>0)then
        jstart=j
      else
        jstart=0
      endif
10    read(201,*,END=11)dum,dumN,dumE,dumZ
      j=j+1
      if(j>0)then
        seisn(j)=dumN/divide
        seise(j)=dumE/divide
        seisz(j)=dumZ/divide
      endif
      if(j<NN)goto 10
11    close(201)
      sampl=j

!detrend
      if(detrend==1)then
        bN=(sum(seisn(jstart+1:sampl)*tt(jstart+1:sampl))-sum(seisn(jstart+1:sampl))*sum(tt(jstart+1:sampl))/dble(sampl-jstart))/(sum(tt(jstart+1:sampl)**2)-sum(tt(jstart+1:sampl))**2/dble(sampl-jstart))
        aN=(sum(seisn(jstart+1:sampl))-bN*sum(tt(jstart+1:sampl)))/dble(sampl-jstart)
        bE=(sum(seise(jstart+1:sampl)*tt(jstart+1:sampl))-sum(seise(jstart+1:sampl))*sum(tt(jstart+1:sampl))/dble(sampl-jstart))/(sum(tt(jstart+1:sampl)**2)-sum(tt(jstart+1:sampl))**2/dble(sampl-jstart))
        aE=(sum(seise(jstart+1:sampl))-bE*sum(tt(jstart+1:sampl)))/dble(sampl-jstart)
        bZ=(sum(seisz(jstart+1:sampl)*tt(jstart+1:sampl))-sum(seisz(jstart+1:sampl))*sum(tt(jstart+1:sampl))/dble(sampl-jstart))/(sum(tt(jstart+1:sampl)**2)-sum(tt(jstart+1:sampl))**2/dble(sampl-jstart))
        aZ=(sum(seisz(jstart+1:sampl))-bZ*sum(tt(jstart+1:sampl)))/dble(sampl-jstart)
        seisn(jstart+1:sampl)=seisn(jstart+1:sampl)-bN*tt(jstart+1:sampl)-aN
        seise(jstart+1:sampl)=seise(jstart+1:sampl)-bE*tt(jstart+1:sampl)-aE
        seisz(jstart+1:sampl)=seisz(jstart+1:sampl)-bZ*tt(jstart+1:sampl)-aZ
      endif

!substract first amplitude
      aN=seisn(jstart+1)
      aE=seise(jstart+1)
      aZ=seisz(jstart+1)
      seisn(jstart+1:sampl)=seisn(jstart+1:sampl)-aN
      seise(jstart+1:sampl)=seise(jstart+1:sampl)-aE
      seisz(jstart+1:sampl)=seisz(jstart+1:sampl)-aZ

      if(f1(FLT)>0.)then
        CALL XAPIIR(seisn(:), NN, 'BU', 0.0, 0.0, 4,'BP', f1(FLT), f2(FLT), dt, 1, NN)
        CALL XAPIIR(seise(:), NN, 'BU', 0.0, 0.0, 4,'BP', f1(FLT), f2(FLT), dt, 1, NN)
        CALL XAPIIR(seisz(:), NN, 'BU', 0.0, 0.0, 4,'BP', f1(FLT), f2(FLT), dt, 1, NN)
      else
        CALL XAPIIR(seisn(:), NN, 'BU', 0.0, 0.0, 4,'LP', f1(FLT), f2(FLT), dt, 1, NN)
        CALL XAPIIR(seise(:), NN, 'BU', 0.0, 0.0, 4,'LP', f1(FLT), f2(FLT), dt, 1, NN)
        CALL XAPIIR(seisz(:), NN, 'BU', 0.0, 0.0, 4,'LP', f1(FLT), f2(FLT), dt, 1, NN)
      endif

      if(INTEGRATE>0)then
        do k=2,NN   !time integration
          seisn(k)=seisn(k)+seisn(k-1)
          seise(k)=seise(k)+seise(k-1)
          seisz(k)=seisz(k)+seisz(k-1)
        enddo
        seisn(:)=seisn(:)*dt
        seise(:)=seise(:)*dt
        seisz(:)=seisz(:)*dt
      endif
      if(INTEGRATE>1)then
        do k=2,NN   !time integration
          seisn(k)=seisn(k)+seisn(k-1)
          seise(k)=seise(k)+seise(k-1)
          seisz(k)=seisz(k)+seisz(k-1)
        enddo
        seisn(:)=seisn(:)*dt
        seise(:)=seise(:)*dt
        seisz(:)=seisz(:)*dt
      endif

      do j=1,NNout !undersampling
        seisnout(j,i)=seisn(int(real(j-1)*UNDER)+1)
        seiseout(j,i)=seise(int(real(j-1)*UNDER)+1)
        seiszout(j,i)=seisz(int(real(j-1)*UNDER)+1)
      enddo
      deallocate(seisn,seise,seisz,tt)
      
    enddo
    close(100)
    close(101)

    open(101,FILE='rvseisn.dat')
    open(102,FILE='rvseise.dat')
    open(103,FILE='rvseisz.dat')
    do i=1,NNout
      write(101,'(100E13.5)')dtout*(i-1),seisnout(i,:)
      write(102,'(100E13.5)')dtout*(i-1),seiseout(i,:)
      write(103,'(100E13.5)')dtout*(i-1),seiszout(i,:)
    enddo
    close(101)
    close(102)
    close(103)
    
    END
