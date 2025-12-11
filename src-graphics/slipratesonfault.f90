    !Plot of velocity function
    IMPLICIT NONE
    INTEGER,PARAMETER:: NM=1  !Column number in the input file (mtilde.dat) to plot
    REAL*8, PARAMETER:: maxampl=.5d0   !maximum amplitude of the slip vel funcs
    REAL*8, PARAMETER:: margin=0.2d0   !relative size of the margin of each fault element
    REAL*8,ALLOCATABLE,DIMENSION(:):: L,W
    REAL*8 tw,tw0
    INTEGER,ALLOCATABLE,DIMENSION(:):: NL,NW
    INTEGER NT,NTP,NSeg
    REAL*8 DT,DL,DW,waterlevel
    REAL*8,ALLOCATABLE:: mtilde(:,:,:),CM(:,:,:)
    REAL*8,ALLOCATABLE:: cumulmtilde(:),rupttime(:,:),risetime(:,:),peaktime(:,:)
    REAL*8 startx,starty,stept,stepa,dum
    INTEGER i,j,k,kk,m
    
    open(10,FILE='input.dat')
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)tw,tw0
    read(10,*)
    read(10,*)dum,NSeg
    allocate(NL(NSeg),NW(NSeg),L(NSeg),W(NSeg))
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)(NL(kk),NW(kk),kk=1,NSeg)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)(L(kk),W(kk),kk=1,NSeg)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)NTP
    close(10)
    
    DT=tw/dble(NTP)
    NT=int(tw0/DT+1.d0)
    open(101,FILE='mtilde.dat')
!    open(102,FILE='inputtf.dat')
    open(205,FILE='slipratesonfault.plot1.dat')
    open(206,FILE='slipratesonfault.plot2.dat')
    open(215,FILE='slipratesonfault.grid.dat')
    open(211,FILE='slipratesonfault.plot.gp')
    open(292,FILE='kinemati_parameters.dat')
    waterlevel=0.
    
    do kk=1,NSeg
      allocate(mtilde(NT,NL(kk),NW(kk)),CM(NT,NL(kk),NW(kk)),cumulmtilde(NT),rupttime(NL(kk),NW(kk)),risetime(NL(kk),NW(kk)),peaktime(NL(kk),NW(kk)))
      L(kk)=L(kk)/1.d3
      W(kk)=W(kk)/1.d3
      DL=L(kk)/dble(NL(kk))
      DW=W(kk)/dble(NW(kk))
      do k=1,NW(kk)
        do j=1,NL(kk)
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
!            read(102,*)(CM(i,j,k),m=1,NM)
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
          peaktime(j,k)=maxval(dt*(maxloc(mtilde(:,j,k))-1))
        enddo
      enddo

      do k=1,NW(kk)
        do j=1,NL(kk)
          write(292,'(1000E13.5)')DL*(dble(j)-.5),DW*(dble(k)-.5),rupttime(j,k),risetime(j,k),maxval(mtilde(:,j,k)),0.d0,0.d0,sum(mtilde(:,j,k))*dt,peaktime(j,k)
        enddo
      enddo
      write(292,*);write(292,*)

      waterlevel=minval(mtilde)
      if(waterlevel>0.)waterlevel=0.
      stept=L(kk)/dble(NL(kk)+1)*(1.d0-margin)/(NT-1)
      stepa=W(kk)/dble(NW(kk)+1)*(1.d0-margin)/(maxampl-waterlevel)

      do i=1,NL(kk)
        startx=((i-1)+margin/2.d0)*L(kk)/dble(NL(kk))
        do j=1,NW(kk)
          starty=((j-1)+margin/2.d0)*W(kk)/dble(NW(kk))
          do k=1,NT
            write(205,'(3E13.5)')startx+stept*(k-1),starty+stepa*(mtilde(k,i,j)-waterlevel)
!            write(206,'(3E13.5)')startx+stept*(k-1),starty+stepa*(CM(k,i,j)-waterlevel)
          enddo
          write(205,*)
!          write(206,*)
        enddo
      enddo
      write(205,'(3E13.5)')0.d0,0.d0
!      write(206,'(3E13.5)')0.d0,0.d0
      write(205,*)
      write(205,*)

!grid
      do i=1,NL(kk)+1
        write(215,*)dble(i-1)*L(kk)/dble(NL(kk)),0.d0
        write(215,*)dble(i-1)*L(kk)/dble(NL(kk)),W(kk)
        write(215,*)
      enddo
      do i=1,NW(kk)+1
        write(215,*)0.d0,dble(i-1)*W(kk)/dble(NW(kk))
        write(215,*)L(kk),dble(i-1)*W(kk)/dble(NW(kk))
        write(215,*)
      enddo
      write(215,*)
      write(215,*)

      write(211,*)'set term postscript landscape color solid enh 16'
      write(211,'(A28,I1,A4)')'set output "slipratesonfault',kk,'.ps"'
      write(211,*)'set pm3d map corners2color c1'
      write(211,*)'set size ratio -1'
      write(211,'(80A)')'set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )'
      write(211,*)'DL=',DL
      write(211,*)'DW=',DW
      write(211,*)'set xrange [',-L(kk)/float(NL(kk))/2.,':',L(kk)+L(kk)/float(NL(kk))/2.,']'
      write(211,*)'set yrange [',-W(kk)/float(NW(kk))/2.,':',W(kk)+W(kk)/float(NW(kk))/2.,']'
      write(211,*)'set xlabel "Along strike (km)'
      write(211,*)'set ylabel "Up-dip (km)"'
      write(211,*)'set cblabel "Slip (m)"'
      write(211,*)'set xtics 5 scale 0.5 out'
      write(211,*)'set ytics 5 scale 0.5 out'
      write(211,*)'set mxtics 5'
      write(211,*)'set mytics 5'

      write(211,*)'set colorbox'
!      write(211,*)'set cbrange [0:.8]'

      write(211,*)'splot "mtildeslip2D.dat" matrix u ($1*DL):($2*DW):3 index \'
      write(211,*)(NM-1)*NSeg+kk-1,' notitle w pm3d, \'
!      write(211,*)'"slipratesonfault.plot2.dat" u 1:2:(0) index ',kk-1,' notitle w l lt 2 lw 2, \'
      write(211,*)'"slipratesonfault.plot1.dat" u 1:2:(0) index \'
      write(211,*)kk-1,' notitle w l lt -1 lw 1, \'
      write(211,*)'"slipratesonfault.grid.dat" u 1:2:(0) index \'
      write(211,*)kk-1,' notitle w l lt 0 lw 1,\'
      write(211,*)'"epic.dat" u 1:2:(0) index ',kk-1,' notitle w p pt 3 lc 3 ps 2.'

      deallocate(mtilde,CM,cumulmtilde,rupttime,risetime,peaktime)
    enddo
    close(101)
!    close(102)
    close(205)
    close(206)
    close(215)
    close(211)
    close(292)

    END
