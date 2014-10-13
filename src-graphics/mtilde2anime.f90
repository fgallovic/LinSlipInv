    IMPLICIT NONE

    INTEGER,PARAMETER:: NM=1  !Column number in the input file (mtilde.dat) to plot
    REAL*8, PARAMETER:: sizex=0.48d0,sizey=0.24d0
    REAL*8 L,W,tw,tw0
    INTEGER NL,NW,NT,NS
    REAL*8 dum,maxik,minik
    REAL*8,ALLOCATABLE:: mtilde(:,:,:)
    REAL*8 DL,DW,DT
    INTEGER i,j,k,m,tfrom,tto,NTP

    open(10,FILE='input.dat')
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)tw,tw0
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)NL,NW
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)L,W
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)NTP
    close(10)
    
    DT=tw/dble(NTP)
    NT=int(tw0/DT+1.d0)
    NS=int(tw0)  !Number of slides equals the slip rate time window
    allocate(mtilde(NT,NL,NW))
    L=L/1.d3
    W=W/1.d3
    DL=L/dble(NL)
    DW=W/dble(NW)

    open(101,FILE='mtilde.dat');write(*,*)'Converting mtilde.dat'
!    open(101,FILE='inputtf.dat');write(*,*)'Converting inputtf.dat'
!    open(101,FILE='GTd.dat');write(*,*)'Converting GTd.dat'
!     open(101,FILE='eigenvectors.dat');write(*,*)'Converting eigenvectors.dat'
    do k=1,NW
      do j=1,NL
        do i=1,NT
          read(101,*)(dum,m=1,NM)
            mtilde(i,j,k)=dum
        enddo
      enddo
    enddo
    close(101)

    maxik=0.d0;minik=0.d0
    open(202,FILE='mtilde2anime.dat')
    do i=1,NS
      do j=1,NW
        tfrom=int(dble((i-1)*NT)/dble(NS)+1)
        tto=int(dble((i)*NT)/dble(NS))
        write(202,'(1000E13.5)')(sum(Mtilde(tfrom:tto,k,j))*DT,k=1,NL),sum(Mtilde(tfrom:tto,NL,j))*DT
      enddo
      write(202,'(1000E13.5)')(sum(Mtilde(tfrom:tto,k,NW))*DT,k=1,NL),sum(Mtilde(tfrom:tto,NL,NW))*DT
      write(202,*);write(202,*)
      dum=maxval(sum(Mtilde(tfrom:tto,:,:),dim=1))*DT
      if(dum>maxik)maxik=dum
      dum=minval(sum(Mtilde(tfrom:tto,:,:),dim=1))*DT
      if(dum<minik)minik=dum
    enddo

    open(201,FILE='mtilde2anime.gp')
    write(201,*)'set term postscript portrait color enh'
    write(201,*)'set output "mtilde2anime.ps"'
!    write(201,*)'set output "inputtf2anime.ps"'
    write(201,*)'set multiplot'
    write(201,*)'set size ',sizex,',',sizey
    write(201,*)'set size ratio -1'
    write(201,*)'set pm3d map corners2color c1'
    write(201,'(80A)')'set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )'
    write(201,*)'DW=',DW
    write(201,*)'DL=',DL
    write(201,*)' W=',W
    write(201,*)'set cbrange [',minik,':',maxik,']'
    write(201,*)'unset colorbox'
    write(201,*)'set xrange [0:',L,']'
    write(201,*)'set yrange [0:',W,']'
    write(201,*)'set xtics out scale 0.5 offset 0,.6'
    write(201,*)'set ytics out scale 0.5 offset 0.6,0'
    write(201,*)'set label 1 "" front at graph .98,graph .88 right'

    do i=1,NS
      write(201,*)'set format x ""'
      write(201,*)'set format y ""'
      write(201,*)'unset xlabel'
      write(201,*)'unset ylabel'
      if(i==NS/2)then
        write(201,*)'set format x "%g"'
        write(201,*)'set format y "%g"'
        write(201,*)'set xlabel "Along strike (km)" offset 0,1'
        write(201,*)'set ylabel "Up-dip (km)"'
      endif
      if(i==NS)then
!        write(201,*)'set colorbox horizontal user origin 0.41,0.05 size 0.3,0.01'
        write(201,*)'set colorbox horizontal user origin 0.41,',real(1.-dble(mod(i-1,NS/2)+1)*sizey/1.7+.04),' size 0.3,0.01'
        write(201,*)'set cblabel "Slip velocity (m/s)" offset 0,.5'
!        write(201,*)'set cbtics 0.1'
      endif
      write(201,'(A13,I2,A1,I2,A2)')'set label 1 "',i-1,'-',i,'s"'
      write(201,*)'set origin ',sizex*((i-1)/(NS/2))/1.5,',',1.-dble(mod(i-1,NS/2)+1)*sizey/1.7

      write(201,'(A58,I5,A20)')'splot "mtilde2anime.dat" matrix u ($1*DL):($2*DW):3 index ',i-1,' notitle w pm3d'  ! bez inv2.dat
    enddo


    END
