      module avg_main
      implicit none
      save

      real*8,parameter :: PI = 3.1415926535897932d0 !@param pi    pi

      contains

!----------------------------------------------------------------------

      subroutine check_err(iret)
      implicit none
      include 'netcdf.inc'
      integer                :: iret
      if(iret /= nf_noerr) then
        print*, nf_strerror(iret)
        stop
      endif
      end subroutine check_err

!----------------------------------------------------------------------
C**** HNTRP.FOR   Horizontal Interpolation Program REAL*8   2004/03/17
C****
      SUBROUTINE HNTRP0 (INA,JNA,OFFIA,DIVJA,
     *                   INB,JNB,OFFIB,DIVJB, SKIB)
C****
C**** HNTRP performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value of SKIP.
C**** The area weighted integral of the quantity is conserved.
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      implicit integer (i-k)
      PARAMETER (TWOPI=6.283185307179586477d0)
C     REAL*8 WTA(INA,JNA),A(INA,JNA),B(INB,JNB),
      REAL*8 WTA(*),      A(*),      B(*),
     *       OFFIA,DIVJA, OFFIB,DIVJB, SKIB,SKIP
      COMMON /HNTRCB/ SINA(0:181),SINB(0:181),
     *       FMIN(360),FMAX(360),GMIN(181),GMAX(181),
     *       IMIN(360),IMAX(360),JMIN(182),JMAX(182)
      LOGICAL*4 QMPOLE
      DATA IMA,JMA,IMB,JMB/4*0/, SKIP/0/
C****
      IMA = INA
      JMA = JNA
      IMB = INB
      JMB = JNB
      SKIP = SKIB
      IF(IMA.lt.1 .or. IMA.gt.360 .or. JMA.lt.1 .or. JMA.gt.181 .or.
     *   IMB.lt.1 .or. IMB.gt.360 .or. JMB.lt.1 .or. JMB.gt.181)
     *   GO TO 400
C****
C**** Partitions in the I direction
C**** RIA = longitude in degrees of right edge of grid box IA on grid A
C**** RIB = longitude in degrees of right edge of grid box IB of grid B
C**** IMIN(IB) = box on grid A containing left edge of box IB on B
C**** IMAX(IB) = box on grid A containing right edge of box IB on B
C**** FMIN(IB) = fraction of box IMIN(IB) on A that is left of box IB
C**** FMAX(IB) = fraction of box IMAX(IB) on A that is right of box IB
C****
      DIA = 360d0/IMA
      DIB = 360d0/IMB
      IA  = 1
      RIA = (IA+OFFIA)*DIA - 360
      IB  = IMB
      DO 150 IBP1=1,IMB
      RIB = (IBP1-1+OFFIB)*DIB
  110 IF(RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GO TO 110
C**** Right edge of A box IA and right edge of B box IB coincide
  130 IMAX(IB) = IA
      FMAX(IB) = 0.
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBP1) = IA
      FMIN(IBP1) = 0.
      GO TO 150
C**** A box IA contains right edge of B box IB
  140 IMAX(IB) = IA
      FMAX(IB) = (RIA-RIB)/DIA
      IMIN(IBP1) = IA
      FMIN(IBP1) = 1.-FMAX(IB)
  150 IB = IBP1
      IMAX(IMB) = IMAX(IMB) + IMA
C       WRITE (0,*) 'IMIN=',(IMIN(I),I=1,IMB)
C       WRITE (0,*) 'IMAX=',(IMAX(I),I=1,IMB)
C       WRITE (0,*) 'FMIN=',(FMIN(I),I=1,IMB)
C       WRITE (0,*) 'FMAX=',(FMAX(I),I=1,IMB)
C****
C**** Partitions in the J direction
C****
C**** RJA = latitude in radians at top edge of box JA on grid A
C**** SINA(JA) = sine of latitude of top edge of box JA on grid A
      OFFJA = (DIVJA-JMA)/2.
      DJA   = .5*TWOPI/DIVJA
      DO 210 JA=1,JMA-1
      RJA = (JA+OFFJA)*DJA - .25*TWOPI
  210 SINA(JA) = DSIN(RJA)
      SINA(0)  = -1.
      SINA(JMA)=  1.
C**** RJB = latitude in radians at top edge of box JB on grid B
C**** SINB(JB) = sine of latitude of top edge of box JB on grid B
      OFFJB = (DIVJB-JMB)/2.
      DJB   = .5*TWOPI/DIVJB
      DO 220 JB=1,JMB-1
      RJB = (JB+OFFJB)*DJB - .25*TWOPI
  220 SINB(JB) = DSIN(RJB)
      SINB(0)  = -1.
      SINB(JMB)=  1.
C****
C**** JMIN(JB) = index of box of A that contains bottom edge of box JB
C**** JMAX(JB) = index of box of A that contains top edge of box JB
C**** GMIN(JB) = fraction of box JMIN(JB) on A grid that is below box JB
C**** GMAX(JB) = fraction of box JMAX(JB) on A grid that is above box JB
C****
      JMIN(1) = 1
      GMIN(1) = 0.
      JA = 1
      DO 350 JB=1,JMB-1
  310 IF(SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GO TO 310
C**** Top edge of A box JA and top edge of B box JB coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0.
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0.
      GO TO 350
C**** A box JA contains top edge of B box JB
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA)-SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB)-SINA(JA-1)
  350 CONTINUE
      JMAX(JMB) = JMA
      GMAX(JMB) = 0.
C       WRITE (0,*) 'JMIN=',(JMIN(J),J=1,JMB)
C       WRITE (0,*) 'JMAX=',(JMAX(J),J=1,JMB)
C       WRITE (0,*) 'GMIN=',(GMIN(J),J=1,JMB)
C       WRITE (0,*) 'GMAX=',(GMAX(J),J=1,JMB)
      RETURN
C****
C**** Invalid parameters or dimensions out of range
C****
  400 WRITE (0,940) IMA,JMA,OFFIA,DIVJA, IMB,JMB,OFFIB,DIVJB, SKIP
      STOP 400
  940 FORMAT ('0Arguments received by HNTRP0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid boxes from',
     *                    ' IDL to left edge of grid box I=1'/
     *  E24.8,' = DIVJA   = number of whole grid boxes from SP to NP'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid boxes from',
     *                    ' IDL to left edge of grid box I=1'/
     *  E24.8,' = DIVJB   = number of whole grid boxes from SP to NP'/
     *  E24.8,' = SKIP    = value to be put in B array when B',
     *  ' grid box is subset of A grid boxes with WTA = 0'/
     *  '0These arguments are invalid or out of range.')
C****

      ENTRY HNTRP (WTA,A,B)
C****
C**** HNTRP performs the horizontal interpolation
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
      QMPOLE = .FALSE.
      GO TO 500

      ENTRY HNTRPP (WTA,A,B)
C****
C**** HNTRPP is similar to HNTRP but polar values are replaced by
C**** their longitudinal mean
C****
      QMPOLE = .TRUE.
C****
C**** Interpolate the A grid onto the B grid
C****
  500 DO 520 JB=1,JMB
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      DO 520 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      WEIGHT= 0.
      VALUE = 0.
      IAMIN = IMIN(IB)
      IAMAX = IMAX(IB)
      DO 510 JA=JAMIN,JAMAX
      G = SINA(JA)-SINA(JA-1)
      IF(JA.eq.JAMIN)  G = G - GMIN(JB)
      IF(JA.eq.JAMAX)  G = G - GMAX(JB)
      DO 510 IAREV=IAMIN,IAMAX
      IA  = 1+MOD(IAREV-1,IMA)
      IJA = IA + IMA*(JA-1)
      F   = 1.
      IF(IAREV.eq.IAMIN)  F = F - FMIN(IB)
      IF(IAREV.eq.IAMAX)  F = F - FMAX(IB)
      WEIGHT = WEIGHT + F*G*WTA(IJA)
  510 VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = SKIP
      IF(WEIGHT.ne.0.)  B(IJB) = VALUE/WEIGHT
  520 continue
C****
C**** Replace individual values near the poles by longitudinal mean
C****
      IF(.NOT.QMPOLE)  RETURN
      DO 630 JB=1,JMB,JMB-1
      WEIGHT = 0.
      VALUE  = 0.
      DO 610 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      IF(B(IJB).eq.SKIP)  GO TO 610
      WEIGHT = WEIGHT + 1.
      VALUE  = VALUE  + B(IJB)
  610 continue
      BMEAN = SKIP
      IF(WEIGHT.ne.0.)  BMEAN = VALUE/WEIGHT
      DO 620 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
  620 B(IJB) = BMEAN
  630 continue
      RETURN

      END SUBROUTINE HNTRP0

!======================================================================!


      subroutine avg_gswp_to_giss(iu_var)

      include 'netcdf.inc'
!
      integer, parameter :: long_max=360,lat_max=180 !GSWP2 dimensions
      integer, parameter :: im=144, jm=90

      real*8, parameter :: SKIP = -1.E30
      integer :: iu_var
      character*80 :: infile,varname

      real*8, dimension(long_max,lat_max) :: WTA,gissWTA,arr2d,giss2d
      real*8, dimension(im,jm) :: OUT
      real*8 OFFIA,DIVJA,OFFIB,DIVJB ! grid averaging variables

      ! Local variables
      integer :: fid,varid,status
      integer i,j


      OFFIA=0.0 !# boxes from IDL to left edge of gridbox I=1 (origin grid)
      DIVJA=180.d0 ! # of whole boxes from SP to NP (original grid)
      OFFIB=0.0 ! # boxes from IDL to left edge of gridbox I=1 (new grid)
      DIVJB = 90.d0 ! # of whole grid boxes from SP to NP (new grid)
                    ! = 180/latitutde spacing = 180/2 deg

      arr2d(:,:)= SKIP
      giss2d(:,:) = SKIP
      WTA(:,:) = 0.d0
      gissWTA(:,:) = 0.d0
      OUT(:,:)=0.d0

!     open file
      infile='/Users/mpuma/Documents/SUITABILITY/soil_orgc_180x360.nc'
      status = nf_open(infile,nf_nowrite,fid)
      call check_err(status)
      varname='orgc'
      status = nf_inq_varid(fid,trim(varname),varid)
      call check_err(status)
      status = nf_get_var_double(fid,varid,arr2d)
      call check_err(status)

      do i=1,long_max
         do j = 1,lat_max
            if( arr2d(i,j)> 0.d0 )then
               WTA(i,j) = 1.d0
	    else
               arr2d(i,j)=SKIP
               WTA(i,j) = 0.d0
            endif
         enddo
      enddo

!     flip matrices to align with GISS input
      do j = 1,lat_max
         giss2d(:,j) = arr2d(:,lat_max-j+1)
         gissWTA(:,j) = WTA(:,lat_max-j+1)
      enddo

!     Interpolate observations to GCM spatial resolution
      call HNTRP0(long_max,lat_max,OFFIA,DIVJA,
     &                 im,jm,OFFIB,DIVJB,SKIP)
      call HNTRPP (gissWTA,giss2d,OUT)

!     WRITE FORCINGS 
      write(iu_var) OUT

      end subroutine avg_gswp_to_giss

!-------------------------------------------------------------------

      subroutine open_outfiles(iu_var)

      integer :: iu_var
      character*80 :: infile

      iu_var = 20
      
      infile='/Users/mpuma/Documents/SUITABILITY/orgc30cm_144x90.bi'
      open(iu_var,FILE=infile, STATUS='UNKNOWN',
     &              FORM='UNFORMATTED')
      
      end subroutine open_outfiles

!-------------------------------------------------------------------

      subroutine close_outfiles(iu_var)
      integer :: iu_var
      integer i
      print *, 'in close subroutine'

      close(iu_var)
      
      end subroutine close_outfiles
!-------------------------------------------------------------------

      end module avg_main

!======================================================================
! Purpose: Read & average soil data                                   =
!======================================================================

      program suitability_avg
      use avg_main

      implicit none
      include 'netcdf.inc'

      integer, parameter :: im=144, jm=90
      integer, parameter :: im_1deg=360, jm_1deg=180
      integer :: iu_var

      call open_outfiles(iu_var)

      call avg_gswp_to_giss(iu_var)
            
      call close_outfiles(iu_var)

      
      end program suitability_avg
