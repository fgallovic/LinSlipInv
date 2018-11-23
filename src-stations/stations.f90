!   Converts stations locations from lat,long to X,Y (X points towards north, Y towards east)
!   --------
!   AUTHOR:  Frantisek Gallovic
!   --------

    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535d0
    REAL*8 x,y,lat,long
    REAL*8 reflat,reflong

    OPEN(110,FILE='stations.in')   !reference point (corresponds to the reference point in the input.dat file)
    read(110,*)
    read(110,*)reflat,reflong
    CLOSE(110)

    OPEN(110,FILE='stations.txt')  !INPUT
    OPEN(111,FILE='stations.dat')  !OUTPUT

10  read(110,*,END=12)lat,long
    CALL POISTA(lat,long,reflat,reflong,x,y)
    write(111,*)x,y,0.    !Station co-ordinates  x(N>0,km),y(E>0,km),z(km)
    goto 10
12  continue
    END
    
       
    
    
    SUBROUTINE POISTA(lat,long,lat0,long0,x,y)
    IMPLICIT NONE
    REAL*8 x,y,lat,long,lat0,long0,azi,dist
    REAL*8 glat,glong,geocn,fltdis,fltazi
    azi=FLTAZI(lat0,long0,lat,long)
    dist=FLTDIS(lat0,long0,lat,long)
    x=dist*cos(azi*0.017453292519943296d0)
    y=dist*sin(azi*0.017453292519943296d0)
    END

    FUNCTION FLTAZI(ELAT,ELON,SLAT,SLON)
! Returns azimut between 2 geodetic points
    REAL*8 ELAT,ELON,SLAT,SLON,AZ,DIST
    REAL*8 GELAT,GELON,GSLAT,GSLON,GEOCN,FLTAZI
    GELAT=GEOCN(ELAT)
    GELON=ELON
    GSLAT=GEOCN(SLAT)
    GSLON=SLON
    CALL AZDIST(GELAT,GELON,GSLAT,GSLON,AZ,DIST)
    FLTAZI=AZ
    RETURN
    END

    FUNCTION FLTDIS(RLAT1,RLON1,RLAT2,RLON2)
! Returns the distance in kilometers of 2 geodetic points
    REAL*8 GEOCN,GELAT,ELON,GSLAT,SLON,AZ,DIST,FLTDIS
    REAL*8 RLAT1,RLON1,RLAT2,RLON2
    GSLAT=GEOCN(RLAT1)
    SLON=RLON1
    GELAT=GEOCN(RLAT2)
    ELON=RLON2
    CALL AZDIST(GELAT,ELON,GSLAT,SLON,AZ,DIST)
    FLTDIS=111.1*DIST
    RETURN
    END



!------------------------------------------------------------------------------
      SUBROUTINE AZDIST(ELAT,ELON,SLAT,SLON,AZ,DIST)
! Returns the epicentral distance in degrees and the azimut (on the ellipsoid)
      REAL*8 PI,RTOD,DTOR
      REAL*8 ELAT,ELON,SLAT,SLON,AZ,DIST
      REAL*8 SLA,SLO,ELA,ELO,SLAC,SLAS,SLOC,SLOS,ELAC,ELAS,ELOC,ELOS,AS,BS,CS,DS,ES,GS,HS,SK,AE,BE,CE,DE,EE,GE,HE,EK,CDIST,SDIST,CSDIST,CAZ,SAZ
      DATA PI/3.1415927D0/,RTOD/57.29578D0/,DTOR/0.0174533D0/

      ELAT=ELAT+1.0E-5
      ELON=ELON+1.0E-5
      SLA=SLAT*DTOR
      SLO=SLON*DTOR
      ELA=ELAT*DTOR
      ELO=ELON*DTOR

      SLAC=DCOS(SLA)
      SLAS=DSIN(SLA)
      SLOC=DCOS(SLO)
      SLOS=DSIN(SLO)

      ELAC=DCOS(ELA)
      ELAS=DSIN(ELA)
      ELOC=DCOS(ELO)
      ELOS=DSIN(ELO)

      AS=SLAC*SLOC
      BS=SLAC*SLOS
      CS=SLAS
      DS=SLOS
      ES=-SLOC
      GS=SLAS*SLOC
      HS=SLAS*SLOS
      SK=-SLAC

      AE=ELAC*ELOC
      BE=ELAC*ELOS
      CE=ELAS
      DE=ELOS
      EE=-ELOC
      GE=ELAS*ELOC
      HE=ELAS*ELOS
      EK=-ELAC

      CDIST=AE*AS+BE*BS+CE*CS
      SDIST=DSQRT(1.0-CDIST*CDIST)
      DIST=RTOD*DATAN2(SDIST,CDIST)
      CSDIST=1./SDIST
      SAZ=-(AS*DE+BS*EE)*CSDIST
      CAZ=-(AS*GE+BS*HE+CS*EK)*CSDIST

      AZ=DATAN2(SAZ,CAZ)
      IF(AZ.LT.0.0)AZ=AZ+2.*PI
      AZ=AZ*RTOD

      RETURN
      END
!------------------------------------------------------------------------------
!Converts geodetic latitude to geocentric
      FUNCTION GEOCN(ALAT)
      REAL*8 RTOD,DTOR,ALAT,GCON,GEOCN
      DATA RTOD/57.29578D0/,DTOR/0.0174533D0/
      GCON=0.9932315D0
      GEOCN=RTOD*ATAN(GCON*(SIN(ALAT*DTOR)/COS(ALAT*DTOR)))
      RETURN
      END
!------------------------------------------------------------------------------
!Converts geocentric latitude to geodetic
      FUNCTION CNGEO(ALAT)
      REAL*8 RTOD,DTOR,ALAT,GCON,CNGEO
      DATA RTOD/57.29578D0/,DTOR/0.0174533D0/
      GCON=0.9932315D0
      CNGEO=atan(tan(ALAT/RTOD)/GCON)/DTOR
      RETURN
      END

