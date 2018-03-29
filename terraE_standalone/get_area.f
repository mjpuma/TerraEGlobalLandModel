      module globe_geom
      implicit none
      save

      real*8,parameter :: PI = 3.1415926535897932d0 !@param pi    pi

      contains

!----------------------------------------------------------------------
      SUBROUTINE GEOM_B(im_grid,jm_grid,dxyp)
!@sum  GEOM_B Calculate spherical geometry for B grid
!@auth Original development team (modifications by G. Schmidt)
!@ver  1.0 (B grid version) (further modification by M.J. Puma)
      implicit none
C**** The primary grid is the A grid (including both poles)
C**** The secondary grid is for the B grid velocities, located on the
C**** vertices of the A grid (Note: first velocity point is J=2)
C**** Polar boxes can have different latitudinal size and are treated
C**** as though they were 1/IM_GRID of their actual area
      integer,intent(in) :: im_grid,jm_grid
      real*8,parameter :: pi = 3.145926535897932d0
      real*8,parameter :: TWOPI = 6.283185307179586477D0
      real*8,parameter :: RADIUS = 6375000.d0
      real*8,parameter :: radian = pi/180d0
      real*8,parameter :: sday = 86400.d0
!@param omega earth's rotation rate (7.29 s^-1)
      real*8,parameter :: EDPERD = 1.d0
      real*8,parameter :: EDPERY = 365.d0
      real*8,parameter :: omega = TWOPI*(EDPERD+EDPERY)/
     &                    (EDPERD*EDPERY*sday)
!@param FIM, BYIM real values related to the number of long. grid boxes
      real*8 :: FIM , BYIM
!@param  DLON grid spacing in longitude (deg)
      REAL*8 :: DLON
C**** For the wonderland model set DLON=DLON/3
c      REAL*8, PARAMETER :: DLON=TWOPI/(IM_GRID*3)
!@param  DLAT,DLAT_DG grid spacing in latitude (rad,deg)
      REAL*8  :: DLAT,DLAT_DG
!@param  FJEQ equatorial value of J
      REAL*8 :: FJEQ
!@var  J1U index of southernmost latitude (currently 2, later 1)
      INTEGER, parameter :: J1U = 2
!@var  JRANGE_HEMI lowest,highest lat index for SH,NH for A,B grid
      INTEGER, dimension(2,2,2) :: JRANGE_HEMI
!@var  LAT latitude of mid point of primary grid box (radians)
      REAL*8, DIMENSION(JM_GRID) :: LAT
!@var  LAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, DIMENSION(JM_GRID,2) :: LAT_DG
!@var  LON longitude of mid points of primary grid box (radians)
      REAL*8, DIMENSION(IM_GRID) :: LON
!@var  LON_DG longitude of mid points of prim. and sec. grid boxes (deg)
      REAL*8, DIMENSION(IM_GRID,2) :: LON_DG

!@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
C**** Note that this is not the exact area, but is what is required for
C**** some B-grid conservation quantities
      REAL*8, DIMENSION(JM_GRID) :: DXYP,BYDXYP

!@var AREAG global integral of area (m^2)
      REAL*8 :: AREAG
!@var WTJ area weighting used in JLMAP, JKMAP (for hemispheric means)
      REAL*8, DIMENSION(JM_GRID,2,2) :: WTJ
!@var DXYV,BYDXYV area of grid box around velocity point (recip.)(m^2)
      REAL*8, DIMENSION(JM_GRID) :: DXYV,BYDXYV
!@var  DXP,DYP,BYDXP,BYDYP distance between points on primary grid
!@+     (+inverse)
      REAL*8, DIMENSION(JM_GRID) :: DXP,DYP,BYDXP,BYDYP
!@var  DXV,DYV distance between velocity points (secondary grid)
      REAL*8, DIMENSION(JM_GRID) :: DXV,DYV
!@var  DXYN,DXYS half box areas to the North,South of primary grid point
      REAL*8, DIMENSION(JM_GRID) :: DXYS,DXYN
!@var  SINP sin of latitude at primary grid points
      REAL*8, DIMENSION(JM_GRID) :: SINP
!@var  COSP, COSV cos of latitude at primary, secondary latitudes
      REAL*8, DIMENSION(JM_GRID) :: COSP,COSV
!@var  RAPVS,RAPVN,RAVPS,RAVPN area scalings for primary and sec. grid
      REAL*8, DIMENSION(JM_GRID) :: RAPVS,RAPVN,RAVPS,RAVPN
!@var SINIV,COSIV,SINIP,COSIP longitud. sin,cos for wind,pressure grid
      REAL*8, DIMENSION(IM_GRID) :: SINIV,COSIV,SINIP,COSIP
!@var  RAVJ scaling for A grid U/V to B grid points (func. of lat. j)
!@var  RAPJ scaling for B grid -> A grid conversion (1/4,1/im_grid at poles)
      REAL*8, DIMENSION(IM_GRID,JM_GRID) :: RAPJ,RAVJ
!@var  IDJJ J index of adjacent U/V points for A grid (func. of lat. j)
      INTEGER, DIMENSION(IM_GRID,JM_GRID) :: IDJJ
!@var  IDIJ I index of adjacent U/V points for A grid (func. of lat/lon)
      INTEGER, DIMENSION(IM_GRID,IM_GRID,JM_GRID) :: IDIJ
!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, DIMENSION(JM_GRID) :: KMAXJ
!@var  IMAXJ varying number of used longitudes
      INTEGER, DIMENSION(JM_GRID) :: IMAXJ
!@var  FCOR latitudinally varying coriolis parameter
      REAL*8, DIMENSION(JM_GRID) :: FCOR

      real*8  :: RAVPO,LAT1,COSP1,DXP1
      real*8 :: acor,acor2,polwt
      integer :: I,J,K,IM1  !@var I,J,K,IM1  loop variables
      integer :: JVPO,JMHALF

!     Compute parameters dependent on grid resolution
      FJEQ=.5*(1+JM_GRID)
      FIM = IM_GRID
      BYIM=1./FIM
      DLON = TWOPI*BYIM
      JRANGE_HEMI = reshape(
     &  (/1,JM_GRID/2,  1+JM_GRID/2,JM_GRID, J1U,J1U-1+JM_GRID/2, 
     &  J1U-1+JM_GRID/2,JM_GRID+J1U-2/),(/2,2,2/))

C**** latitudinal spacing depends on whether you have even spacing or
C**** a partial box at the pole
      DLAT_DG=180./JM_GRID                     ! even spacing (default)
      IF (JM_GRID.eq.46) DLAT_DG=180./(JM_GRID-1)   ! 1/2 box at pole for 4x5
cc    IF (JM_GRID.eq.24) DLAT_DG=180./(JM_GRID-1)   ! 1/2 box @pole, orig 8x10
      IF (JM_GRID.eq.24) DLAT_DG=180./(JM_GRID-1.5) ! 1/4 box @pole,'real'8x10
      DLAT=DLAT_DG*radian
      LAT(1)  = -.25*TWOPI
      LAT(JM_GRID) = -LAT(1)
      SINP(1)  = -1.
      SINP(JM_GRID) = 1.
      COSP(1)  = 0.
      COSP(JM_GRID) = 0.
      DXP(1)  = 0.
      DXP(JM_GRID) = 0.
      DO J=2,JM_GRID-1
        LAT(J)  = DLAT*(J-FJEQ)
        SINP(J) = SIN(LAT(J))
        COSP(J) = COS(LAT(J))
        DXP(J)  = RADIUS*DLON*COSP(J)
      END DO
      BYDXP(2:JM_GRID-1) = 1.D0/DXP(2:JM_GRID-1)
      LAT1    = DLAT*(1.-FJEQ)
      COSP1   = COS(LAT1)
      DXP1    = RADIUS*DLON*COSP1
      DO J=2,JM_GRID
        COSV(J) = .5*(COSP(J-1)+COSP(J))
        DXV(J)  = .5*(DXP(J-1)+DXP(J))
        DYV(J)  = RADIUS*(LAT(J)-LAT(J-1))
C**** The following corrections have no effect for half polar boxes
C**** but are important for full and quarter polar box cases.
        IF (J.eq.2) THEN
          polwt = cosv(j)
          COSV(J) = .5*(COSP1+COSP(J))
          DXV(J)  = .5*(DXP1+DXP(J))
        END IF
        IF (J.eq.JM_GRID) THEN
          COSV(J) = .5*(COSP(J-1)+COSP1)
          DXV(J)  = .5*(DXP(J-1)+DXP1)
        END IF
C****
      END DO
      DYP(1)  = RADIUS*(LAT(2)-LAT(1)-0.5*DLAT)
      DYP(JM_GRID) = RADIUS*(LAT(JM_GRID)-LAT(JM_GRID-1)-0.5*DLAT)
      DXYP(1) = .5*DXV(2)*DYP(1)
      BYDXYP(1) = 1./DXYP(1)
      DXYP(JM_GRID)= .5*DXV(JM_GRID)*DYP(JM_GRID)
      BYDXYP(JM_GRID) = 1./DXYP(JM_GRID)
      DXYS(1)  = 0.
      DXYS(JM_GRID) = DXYP(JM_GRID)
      DXYN(1)  = DXYP(1)
      DXYN(JM_GRID) = 0.
      polwt = (cosv(3)-cosv(2))/(cosv(3)-polwt)
      AREAG = DXYP(1)+DXYP(JM_GRID)
      DO J=2,JM_GRID-1
        DYP(J)  =  radius*dlat !.5*(DYV(J)+DYV(J+1))
        DXYP(J) = .5*(DXV(J)+DXV(J+1))*DYP(J)
        BYDXYP(J) = 1./DXYP(J)
        DXYS(J) = .5*DXYP(J)
        DXYN(J) = .5*DXYP(J)
        AREAG = AREAG+DXYP(J)
      END DO
      BYDYP(:) = 1.D0/DYP(:)
      AREAG = AREAG*FIM
      RAVPS(1)  = 0.
      RAPVS(1)  = 0.
      RAVPN(JM_GRID) = 0.
      RAPVN(JM_GRID) = 0.
      DO J=2,JM_GRID
        DXYV(J) = DXYN(J-1)+DXYS(J)
        BYDXYV(J) = 1./DXYV(J)
        RAPVS(J)   = .5*DXYS(J)/DXYV(J)
        RAPVN(J-1) = .5*DXYN(J-1)/DXYV(J)
        RAVPS(J)   = .5*DXYS(J)/DXYP(J)
        RAVPN(J-1) = .5*DXYN(J-1)/DXYP(J-1)
      END DO
      acor = dxyv(2)/(.5*dxp(2)*dyv(2)) ! gridbox area correction factor
      acor2 = dxyv(2)/(dxv(2)*dyv(2))
C**** LONGITUDES (degrees); used in ILMAP
      LON_DG(1,1) = -180.+360./(2.*FLOAT(IM_GRID))
      LON_DG(1,2) = -180.+360./    FLOAT(IM_GRID)
      DO I=2,IM_GRID
        LON_DG(I,1) = LON_DG(I-1,1)+360./FLOAT(IM_GRID)
        LON_DG(I,2) = LON_DG(I-1,2)+360./FLOAT(IM_GRID)
      END DO
C**** LATITUDES (degrees); used extensively in the diagn. print routines
      LAT_DG(1,1:2)=-90.
      LAT_DG(JM_GRID,1)=90.
      DO J=2,JM_GRID-1
        LAT_DG(J,1)=DLAT_DG*(J-FJEQ)    ! primary (tracer) latitudes
      END DO
      DO J=2,JM_GRID
        LAT_DG(J,2)=DLAT_DG*(J-JM_GRID/2-1)  ! secondary (velocity) latitudes
      END DO
C**** WTJ: area weighting for JKMAP, JLMAP hemispheres
      JMHALF= JM_GRID/2
      DO J=1,JM_GRID
        WTJ(J,1,1)=1.
        WTJ(J,2,1)=2.*FIM*DXYP(J)/AREAG
      END DO
      DO J=2,JM_GRID
        WTJ(J,1,2)=1.
        WTJ(J,2,2)=2.*FIM*DXYV(J)/AREAG
      END DO
cgsfc      WTJ(JMHALF+1,1,2)=.5
cgsfc      WTJ(JMHALF+1,2,2)=WTJ(JMHALF+1,2,2)/2.
      WTJ(1,1,2)=0.
      WTJ(1,2,2)=0.
C**** CALCULATE CORIOLIS PARAMETER
c      OMEGA = TWOPI*(EDPERD+EDPERY)/(EDPERD*EDPERY*SDAY)
      FCOR(1)  = -OMEGA*DXV(2)*DXV(2)/DLON
      FCOR(JM_GRID) = OMEGA*DXV(JM_GRID)*DXV(JM_GRID)/DLON
      DO J=2,JM_GRID-1
        FCOR(J) = OMEGA*(DXV(J)*DXV(J)-DXV(J+1)*DXV(J+1))/DLON
      END DO

C**** Set indexes and scalings for the influence of A grid points on
C**** adjacent velocity points

C**** Calculate relative directions of polar box to nearby U,V points
      DO I=1,IM_GRID
        SINIV(I)=SIN((I-1)*DLON)
        COSIV(I)=COS((I-1)*TWOPI*BYIM) ! DLON)
        LON(I)=DLON*(I-.5)
        SINIP(I)=SIN(LON(I))
        COSIP(I)=COS(LON(I))
      END DO

C**** Conditions at the poles
      DO J=1,JM_GRID,JM_GRID-1
        IF(J.EQ.1) THEN
          JVPO=2
          RAVPO=2.*RAPVN(1)
        ELSE
          JVPO=JM_GRID
          RAVPO=2.*RAPVS(JM_GRID)
        END IF
        KMAXJ(J)=IM_GRID
        IMAXJ(J)=1
        RAVJ(1:KMAXJ(J),J)=RAVPO
        RAPJ(1:KMAXJ(J),J)=BYIM
        IDJJ(1:KMAXJ(J),J)=JVPO
        DO K=1,KMAXJ(J)
          IDIJ(K,1:IM_GRID,J)=K
        END DO
      END DO
C**** Conditions at non-polar points
      DO J=2,JM_GRID-1
        KMAXJ(J)=4
        IMAXJ(J)=IM_GRID
        DO K=1,2
          RAVJ(K,J)=RAPVS(J)
          RAPJ(K,J)=RAVPS(J)    ! = .25
          IDJJ(K,J)=J
          RAVJ(K+2,J)=RAPVN(J)
          RAPJ(K+2,J)=RAVPN(J)  ! = .25
          IDJJ(K+2,J)=J+1
        END DO
        IM1=IM_GRID
        DO I=1,IM_GRID
          IDIJ(1,I,J)=IM1
          IDIJ(2,I,J)=I
          IDIJ(3,I,J)=IM1
          IDIJ(4,I,J)=I
          IM1=I
        END DO
      END DO
      
      RETURN
      END SUBROUTINE GEOM_B

!----------------------------------------------------------------------

      end module globe_geom


!======================================================================
! Purpose: Get grid cell areas for different model resolutions        =
!======================================================================

      program get_area
      use globe_geom

      implicit none
      integer :: i,j

      integer, parameter :: im=144, jm=90
      integer, parameter :: im_1deg=7200, jm_1deg=3600
!      integer, parameter :: im_1deg=360, jm_1deg=180
      real*8,dimension(jm) :: lat_area
      real*8,dimension(im,jm) :: area
      real*8,dimension(jm_1deg) :: lat_area1deg
!      real*8,dimension(im_1deg,jm_1deg) :: area1deg
!      call GEOM_B(im,jm,lat_area)
!      do i=1,im
!         do j = 1,jm
!            area(i,j) = lat_area(j)/1000000.d0 !convert m^2 -> km^2
!         enddo
!      enddo
!       write(900) area !note:compile with -convert big_endian

      call GEOM_B(im_1deg,jm_1deg,lat_area1deg)
 !     do i=1,im_1deg
 !        do j = 1,jm_1deg
 !           area1deg(i,j) = lat_area1deg(j)/1000000.d0 !convert m^2 -> km^2
 !        enddo
 !     enddo
!
!       write(900) area1deg !note:compile with -convert big_endian
       write(900) lat_area1deg

      end program get_area
