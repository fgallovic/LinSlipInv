!   Plot of seismogram comparison
!   -----------------------------
!   AUTHOR:  Frantisek Gallovic

    IMPLICIT NONE
    REAL*8, PARAMETER:: margin=0.05d0   !relative margin size
    REAL*8,ALLOCATABLE,DIMENSION(:,:):: rseisn,rseise,rseisz,sseisn,sseise,sseisz
    REAL*8,ALLOCATABLE,DIMENSION(:):: stepa,maxampl
    INTEGER,ALLOCATABLE,DIMENSION(:,:):: stainfo
    CHARACTER*4,ALLOCATABLE,DIMENSION(:):: staname
    REAL*8 TW,DT,T1,T2
    INTEGER NS,FFT,NT,NTFROM
    REAL*8 startx,starty,stept,dum
    REAL*8 maxampln,maxample,maxamplz
    INTEGER i,j,k,m
    
    open(100,FILE='input.dat')
    read(100,*)
    read(100,*)
    read(100,*)
    read(100,*)TW,dum,T1,T2
    read(100,*)
    read(100,*)
    read(100,*)
    read(100,*)NS
    do i=1,13
      read(100,*)
    enddo
    read(100,*)FFT
    close(100)
    DT=TW/dble(FFT)
    NT=floor((T2-T1)/DT+1)
    NTfrom=floor(T1/DT+1)
    write(*,*)NT

    allocate(rseisn(FFT,NS),rseise(FFT,NS),rseisz(FFT,NS),sseisn(FFT,NS),sseise(FFT,NS),sseisz(FFT,NS))
    allocate(stepa(NS),maxampl(NS),staname(NS),stainfo(3,NS))

    open(100,FILE='stations.dat')
    do k=1,NS
      read(100,*)dum,dum,dum,staname(k)
    enddo
    close(100)
    open(100,FILE='stainfo.dat')
    do k=1,NS
      read(100,*)stainfo(1:3,k)
    enddo
    close(100)

! In typical applications rvseis[nez].dat should contain the same seismogram as in rvseisnez.dat
!    open(101,FILE='rvseisn.dat')
!    open(102,FILE='rvseise.dat')
!    open(103,FILE='rvseisz.dat')
!    do k=1,FFT
!      read(101,*)dum,rseisn(k,1:NS)
!      read(102,*)dum,rseise(k,1:NS)
!      read(103,*)dum,rseisz(k,1:NS)
!    enddo
!    close(101)
!    close(102)
!    close(103)
    
    open(105,FILE='rvseisnez.dat')
    sseisn=0.d0;sseise=0.d0;sseisz=0.d0
    do k=int(T1/DT)+1,int(T2/DT)+1
      read(105,'(E13.5)',ADVANCE='NO')dum
      do i=1,NS
        if(stainfo(1,i)==1)read(105,'(E13.5)',ADVANCE='NO')sseisn(k,i)
        if(stainfo(2,i)==1)read(105,'(E13.5)',ADVANCE='NO')sseise(k,i)
        if(stainfo(3,i)==1)read(105,'(E13.5)',ADVANCE='NO')sseisz(k,i)
      enddo
      read(105,*)
    enddo
    close(105)

    stept=1./dble(3)*(1.d0-margin)/(NT-1)
    do k=1,NS
      maxampln=maxval(abs(sseisn(NTfrom:NT+NTfrom,k)))
      maxample=maxval(abs(sseise(NTfrom:NT+NTfrom,k)))
      maxamplz=maxval(abs(sseisz(NTfrom:NT+NTfrom,k)))
      maxampl(k)=max(maxampln,max(maxample,maxamplz))
      stepa(k)=.5d0/dble(NS)*(1.d0-margin)/maxampl(k)
    enddo

    open(205,FILE='SlipInv-seisplot.dat')
    !do j=1,NS
    !  starty=((NS-j+1)+margin/2.d0)/dble(NS+1)
    !  startx=((1-1)+margin/2.d0)/3.d0
    !  do k=NTfrom,NT+NTfrom-1
    !    write(205,'(3E13.5)')startx+stept*(k-NTfrom),starty+stepa(j)*rseisn(k,j)
    !  enddo
    !  write(205,*)
    !  startx=((2-1)+margin/2.d0)/3.d0
    !  do k=NTfrom,NT+NTfrom-1
    !    write(205,'(3E13.5)')startx+stept*(k-NTfrom),starty+stepa(j)*rseise(k,j)
    !  enddo
    !  write(205,*)
    !  startx=((3-1)+margin/2.d0)/3.d0
    !  do k=NTfrom,NT+NTfrom-1
    !    write(205,'(3E13.5)')startx+stept*(k-NTfrom),starty+stepa(j)*rseisz(k,j)
    !  enddo
    !  write(205,*)
    !enddo
    !write(205,*)
    !write(205,*)
    
    do j=1,NS
      starty=((NS-j+1)+margin/2.d0)/dble(NS+1)
      if(stainfo(1,j)==1)then
        startx=((1-1)+margin/2.d0)/3.d0
        do k=NTfrom,NT+NTfrom-1
          write(205,'(3E13.5)')startx+stept*(k-NTfrom),starty+stepa(j)*sseisn(k,j)
        enddo
        write(205,*)
      endif
      if(stainfo(2,j)==1)then
        startx=((2-1)+margin/2.d0)/3.d0
        do k=NTfrom,NT+NTfrom-1
          write(205,'(3E13.5)')startx+stept*(k-NTfrom),starty+stepa(j)*sseise(k,j)
        enddo
        write(205,*)
      endif
      if(stainfo(3,j)==1)then
        startx=((3-1)+margin/2.d0)/3.d0
        do k=NTfrom,NT+NTfrom-1
          write(205,'(3E13.5)')startx+stept*(k-NTfrom),starty+stepa(j)*sseisz(k,j)
        enddo
        write(205,*)
      endif
    enddo

    open(105,FILE='svseisnez.dat')
    sseisn=0.d0;sseise=0.d0;sseisz=0.d0
    do k=int(T1/DT)+1,int(T2/DT)+1
      read(105,'(E13.5)',ADVANCE='NO')dum
      do i=1,NS
        if(stainfo(1,i)==1)read(105,'(E13.5)',ADVANCE='NO')sseisn(k,i)
        if(stainfo(2,i)==1)read(105,'(E13.5)',ADVANCE='NO')sseise(k,i)
        if(stainfo(3,i)==1)read(105,'(E13.5)',ADVANCE='NO')sseisz(k,i)
      enddo
      read(105,*)
    enddo
    close(105)
    write(205,*)
    write(205,*)
    do j=1,NS
      starty=((NS-j+1)+margin/2.d0)/dble(NS+1)
      if(stainfo(1,j)==1)then
        startx=((1-1)+margin/2.d0)/3.d0
        do k=NTfrom,NT+NTfrom-1
          write(205,'(3E13.5)')startx+stept*(k-NTfrom),starty+stepa(j)*sseisn(k,j)
        enddo
        write(205,*)
      endif
      if(stainfo(2,j)==1)then
        startx=((2-1)+margin/2.d0)/3.d0
        do k=NTfrom,NT+NTfrom-1
          write(205,'(3E13.5)')startx+stept*(k-NTfrom),starty+stepa(j)*sseise(k,j)
        enddo
        write(205,*)
      endif
      if(stainfo(3,j)==1)then
        startx=((3-1)+margin/2.d0)/3.d0
        do k=NTfrom,NT+NTfrom-1
          write(205,'(3E13.5)')startx+stept*(k-NTfrom),starty+stepa(j)*sseisz(k,j)
        enddo
        write(205,*)
      endif
    enddo
    close(205)

    open(201,FILE='SlipInv-seisplot.gp')
    write(201,*)'set term postscript portrait color solid enh'
    write(201,*)'set output "SlipInv-seisplot.ps"'
    write(201,*)'set xrange [0:1]'
    write(201,*)'set yrange [0:1]'
    write(201,*)'unset xtics'
    write(201,*)'unset ytics'
    write(201,*)'set border 0'

    do j=1,NS
      write(201,'(A11,F5.2,A23,G,A6)')'set label "',real(maxampl(j))*100.,'" at graph 1.09, graph ',((NS-j+1)+margin/2.d0)/dble(NS+1),' right'
      write(201,*)'set label "'//staname(j)//'" at graph -0.01, graph ',((NS-j+1)+margin/2.d0)/dble(NS+1),' right'
    enddo
    write(201,*)'set label "N-S" at .15,1'
    write(201,*)'set label "E-W" at .48,1'
    write(201,*)'set label "Z" at .86,1'

    write(201,*)'plot "SlipInv-seisplot.dat" u 1:2 index 0 notitle w l lt -1 lw 1,\'
    write(201,*)'"SlipInv-seisplot.dat" u 1:2 index 1 notitle w l lt 1 lw 2'

    END
