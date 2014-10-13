    !Plot of velocity function
    IMPLICIT NONE
    INTEGER,PARAMETER:: NM=1  !Column number in the input file (mtilde.dat) to plot
    REAL*8, PARAMETER:: maxampl=1.5d0   !maximum amplitude of the slip vel funcs
    REAL*8, PARAMETER:: margin=0.2d0   !relative size of the margin of each fault element
    REAL*8 L,W,tw,tw0
    INTEGER NL,NW,NT,NTP
    REAL*8 DT,DL,DW,waterlevel
    REAL*8,ALLOCATABLE:: mtilde(:,:,:),CM(:,:,:)
    REAL*8,ALLOCATABLE:: cumulmtilde(:),rupttime(:,:),risetime(:,:),peaktime(:,:,:)
    REAL*8 startx,starty,stept,stepa,dum
    INTEGER i,j,k,m
    
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
    allocate(mtilde(NT,NL,NW),CM(NT,NL,NW),cumulmtilde(NT),rupttime(NL,NW),risetime(NL,NW),peaktime(1,NL,NW))
    L=L/1.d3
    W=W/1.d3
    DL=L/dble(NL)
    DW=W/dble(NW)
    
    open(101,FILE='mtilde.dat')
!    open(102,FILE='inputtf.dat')
    do k=1,NW
      do j=1,NL
        cumulmtilde=0.d0
        do i=1,NT
          read(101,*)(dum,m=1,NM)
          if(dum>waterlevel)then
            mtilde(i,j,k)=dum
          else
            mtilde(i,j,k)=0.d0
          endif 
          if(i==1)then
            cumulmtilde(1)=mtilde(1,j,k)
          else
            cumulmtilde(i)=cumulmtilde(i-1)+mtilde(i,j,k)
          endif
!          read(102,*)(CM(i,j,k),m=1,NM)
        enddo
        do i=1,NT
          if(cumulmtilde(i)<cumulmtilde(NT)*.1)then
            rupttime(j,k)=dt*(i-1)
            cycle
          elseif(cumulmtilde(i)<cumulmtilde(NT)*.9)then
            risetime(j,k)=dt*(i-1)-rupttime(j,k)
            cycle
          endif
        enddo
        peaktime(:,j,k)=dt*(maxloc(mtilde(:,j,k))-1)
      enddo
    enddo
    close(101)
!    close(102)

    waterlevel=minval(mtilde)
    if(waterlevel>0.)waterlevel=0.
    stept=L/dble(NL+1)*(1.d0-margin)/(NT-1)
    stepa=W/dble(NW+1)*(1.d0-margin)/(maxampl-waterlevel)

    open(205,FILE='slipratesonfault.plot1.dat')
    open(206,FILE='slipratesonfault.plot2.dat')

    do i=1,NL
      startx=((i-1)+margin/2.d0)*L/dble(NL)
      do j=1,NW
        starty=((j-1)+margin/2.d0)*W/dble(NW)
        do k=1,NT
          write(205,'(3E13.5)')startx+stept*(k-1),starty+stepa*(mtilde(k,i,j)-waterlevel)
!          write(206,'(3E13.5)')startx+stept*(k-1),starty+stepa*(CM(k,i,j)-waterlevel)
        enddo
        write(205,*)
!        write(206,*)
      enddo
    enddo
    write(205,'(3E13.5)')0.d0,0.d0
    write(206,'(3E13.5)')0.d0,0.d0
    close(205)
    close(206)
    open(205,FILE='slipratesonfault.grid.dat')
    do i=1,NL-1
      write(205,*)dble(i)*L/dble(NL),0.d0
      write(205,*)dble(i)*L/dble(NL),W
      write(205,*)
      write(205,*)
    enddo
    do i=1,NW-1
      write(205,*)0.d0,dble(i)*W/dble(NW)
      write(205,*)L,dble(i)*W/dble(NW)
      write(205,*)
      write(205,*)
    enddo

    open(201,FILE='slipratesonfault.plot.gp')
    write(201,*)'set term postscript landscape color solid enh 20'
    write(201,*)'set output "slipratesonfault.ps"'
    write(201,*)'set pm3d map corners2color c1'
    write(201,*)'set size ratio -1'
    write(201,'(80A)')'set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )'
    write(201,*)'DL=',DL
    write(201,*)'DW=',DW
    write(201,*)'set xrange [',-L/float(NL)/2.,':',L+L/float(NL)/2.,']'
    write(201,*)'set yrange [',-W/float(NW)/2.,':',W+W/float(NW)/2.,']'
    write(201,*)'set xlabel "Along strike (km)'
    write(201,*)'set ylabel "Up-dip (km)"'
    write(201,*)'set cblabel "Slip (m)"'
    write(201,*)'set xtics 5 scale 0.5 out'
    write(201,*)'set ytics 5 scale 0.5 out'
    write(201,*)'set mxtics 5'
    write(201,*)'set mytics 5'

    write(201,*)'set colorbox'
!    write(201,*)'set cbrange [0:.8]'

    write(201,*)'splot "mtildeslip2D.dat" matrix u ($1*DL):($2*DW):3 index \'
    write(201,*)NM-1,' notitle w pm3d, \'
!    write(201,*)'"slipratesonfault.plot2.dat" u 1:2:(0) notitle w l lt 2 lw 2, \'
    write(201,*)'"slipratesonfault.plot1.dat" u 1:2:(0) notitle w l lt -1 lw 1, \'
    write(201,*)'"slipratesonfault.grid.dat" u 1:2:(0) notitle w l lt 0 lw 1.'

    END
