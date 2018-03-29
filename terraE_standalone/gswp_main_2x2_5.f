      module drv_main_gswp
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

C**** HNTRPS.FOR   Horizontal Interpolation Program REAL*8   2004/03/17
      SUBROUTINE HNTRP0 (INA,JNA,OFFIA,DIVJA,
     &                   INB,JNB,OFFIB,DIVJB, SKIB)
C**** HNTRP performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value of SKIP.
C**** The area weighted integral of the quantity is conserved.

      IMPLICIT REAL*8 (A-H,O-Z)
      implicit integer (i-k)
      PARAMETER (TWOPI=6.283185307179586477d0)
C     REAL*4 WTA(INA,JNA),A(INA,JNA),B(INB,JNB),
      REAL*4 WTA(*),      A(*),      B(*),
     *       OFFIA,DIVJA, OFFIB,DIVJB, SKIB,SKIP
      COMMON /HNTRCB/ SINA(0:3601),SINB(0:3601),
     *       FMIN(7202),FMAX(7202),GMIN(3601),GMAX(3601),
     *       IMIN(7202),IMAX(7202),JMIN(3602),JMAX(3602)
      LOGICAL*4 QMPOLE
      DATA IMA,JMA,IMB,JMB/4*0/, SKIP/0/
C****
      IMA = INA
      JMA = JNA
      IMB = INB
      JMB = JNB
      SKIP = SKIB
      IF(IMA.lt.1 .or. IMA.gt.7202 .or. JMA.lt.1 .or. JMA.gt.3601 .or.
     *   IMB.lt.1 .or. IMB.gt.7202 .or. JMB.lt.1 .or. JMB.gt.3601)
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


      subroutine avg_gswp_to_giss(n_monthyr,nt,monlen,iu_vector)

      include 'netcdf.inc'
!
      integer,intent(in) :: n_monthyr ! simulation month number
      integer,intent(in) :: nt ! current month timestep (3-hourly) number
      integer,intent(in) :: monlen ! # of timesteps in current month
!
      integer, parameter :: long_max=360,lat_max=150 !GSWP2 dimensions
      integer, parameter :: lat_maxfull=180 ! # of 1 deg latitude divisions
      real*4, parameter :: SKIP = -1.E30
!
      integer, parameter :: im=144, jm=90
      integer, parameter :: num_months = 162
      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer, parameter :: num_outfiles = num_months*varnums
      integer, parameter :: num_landcells = 15238
      integer,dimension(varnums) :: iu_vector

      character*80 :: varname, infile
      character*80 :: basefolder, foldnames(varnums),vnames(varnums)
      character*80 :: v2names(varnums), monthyear(num_months)

      real*4, dimension(long_max,lat_max) :: WTA
      real*4, dimension(long_max,lat_maxfull) :: gswpWTA,gissWTA,giss2d
      real*4, dimension(long_max,lat_max) :: arr2d,lat,long
      real*4, dimension(long_max,lat_maxfull) :: gswp2d
      real*4, dimension(im,jm) :: OUT

      real*4 OFFIA,DIVJA,OFFIB,DIVJB ! grid averaging variables

      ! Local variables
      integer :: fid,varid,status
      integer :: status_alloc1 ! status = 0 if success
      integer i,j,k,kk,land

      integer,dimension(num_landcells) :: xy
      real*4 ,dimension(num_landcells) :: var_nt
      real*4, allocatable, dimension(:,:) :: arrpack

      real*4, dimension(im,jm,varnums) :: gswp2X2_5_all

      DATA foldnames/
     &         'SWdown_srb/','LWdown_srb/','Rainf_gswp/','Snowf_gswp/'
     &        ,'Qair_cru/','Rainf_C_gswp/','Tair_cru/'
     &        ,'PSurf_ecor/','Wind_ncep/'/
      DATA vnames/
     &         'SWdown_srb','LWdown_srb','Rainf_gswp', 'Snowf_gswp'
     &        ,'Qair_cru','Rainf_C_gswp','Tair_cru'
     &        ,'PSurf_ecor','Wind_ncep'/
      DATA v2names/
     &         'SWdown','LWdown','Rainf', 'Snowf'
     &        ,'Qair','Rainf_C','Tair'
     &        ,'PSurf','Wind'/

      DATA monthyear/          '198207.nc','198208.nc',
     & '198209.nc','198210.nc','198211.nc','198212.nc',
     & '198301.nc','198302.nc','198303.nc','198304.nc',
     & '198305.nc','198306.nc','198307.nc','198308.nc',
     & '198309.nc','198310.nc','198311.nc','198312.nc',
     & '198401.nc','198402.nc','198403.nc','198404.nc',
     & '198405.nc','198406.nc','198407.nc','198408.nc',
     & '198409.nc','198410.nc','198411.nc','198412.nc',
     & '198501.nc','198502.nc','198503.nc','198504.nc',
     & '198505.nc','198506.nc','198507.nc','198508.nc',
     & '198509.nc','198510.nc','198511.nc','198512.nc',
     & '198601.nc','198602.nc','198603.nc','198604.nc',
     & '198605.nc','198606.nc','198607.nc','198608.nc',
     & '198609.nc','198610.nc','198611.nc','198612.nc',
     & '198701.nc','198702.nc','198703.nc','198704.nc',
     & '198705.nc','198706.nc','198707.nc','198708.nc',
     & '198709.nc','198710.nc','198711.nc','198712.nc',
     & '198801.nc','198802.nc','198803.nc','198804.nc',
     & '198805.nc','198806.nc','198807.nc','198808.nc',
     & '198809.nc','198810.nc','198811.nc','198812.nc',
     & '198901.nc','198902.nc','198903.nc','198904.nc',
     & '198905.nc','198906.nc','198907.nc','198908.nc',
     & '198909.nc','198910.nc','198911.nc','198912.nc',
     & '199001.nc','199002.nc','199003.nc','199004.nc',
     & '199005.nc','199006.nc','199007.nc','199008.nc',
     & '199009.nc','199010.nc','199011.nc','199012.nc',
     & '199101.nc','199102.nc','199103.nc','199104.nc',
     & '199105.nc','199106.nc','199107.nc','199108.nc',
     & '199109.nc','199110.nc','199111.nc','199112.nc',
     & '199201.nc','199202.nc','199203.nc','199204.nc',
     & '199205.nc','199206.nc','199207.nc','199208.nc',
     & '199209.nc','199210.nc','199211.nc','199212.nc',
     & '199301.nc','199302.nc','199303.nc','199304.nc',
     & '199305.nc','199306.nc','199307.nc','199308.nc',
     & '199309.nc','199310.nc','199311.nc','199312.nc',
     & '199401.nc','199402.nc','199403.nc','199404.nc',
     & '199405.nc','199406.nc','199407.nc','199408.nc',
     & '199409.nc','199410.nc','199411.nc','199412.nc',
     & '199501.nc','199502.nc','199503.nc','199504.nc',
     & '199505.nc','199506.nc','199507.nc','199508.nc',
     & '199509.nc','199510.nc','199511.nc','199512.nc' /

      basefolder = '/discover/nobackup/mjpuma/GSWP/'

      ! intialize variables--------
      gswp2X2_5_all(:,:,:) = 0.d0

      loop_gswp_variables: do kk = 1, varnums

         OFFIA=0.0 !# boxes from IDL to left edge of gridbox I=1 (origin grid)
         DIVJA=180.d0 ! # of whole boxes from SP to NP (original grid)
         OFFIB=0.0 ! # boxes from IDL to left edge of gridbox I=1 (new grid)
         DIVJB = 90.d0 ! # of whole grid boxes from SP to NP (new grid)
                       ! = 180/latitutde spacing = 180/2

         arr2d(:,:)= 0.d0!SKIP
         giss2d(:,:) = SKIP
         gswp2d(:,:) = SKIP
         WTA(:,:) = 0.d0
         gissWTA(:,:) = 0.d0
         gswpWTA(:,:) = 0.d0
         OUT(:,:)=0.d0

         ! open file
         infile=basefolder(1:len_trim(basefolder))//
     &       foldnames(kk)(1:len_trim(foldnames(kk)))//
     &       vnames(kk)(1:len_trim(vnames(kk)))//
     &       monthyear(n_monthyr)(1:len_trim(monthyear(n_monthyr)))
         status = nf_open(infile,nf_nowrite,fid)
         call check_err(status)
!        ! read in x y indices
         varname='land'
         status = nf_inq_varid(fid,trim(varname),varid)
         call check_err(status)
         status = nf_get_var_int(fid,varid,xy)
         call check_err(status)

         allocate(arrpack(1:num_landcells,1:monlen),STAT=status_alloc1)
         allocate_ok: if(status_alloc1 == 0) then

            ! read in data over land points 
            status = nf_inq_varid(fid,trim(v2names(kk)),varid)
            call check_err(status)
            status = nf_get_var_real(fid,varid,arrpack)
            call check_err(status)
            status = nf_close(fid)
            call check_err(status)

            ! separate data for this timestep (nt)
            var_nt = arrpack(:,nt)

            ! spread data to full x y grid
            do k = 1,num_landcells
               j = xy(k)/long_max +1
               i = xy(k) - (j-1)*long_max

               if(var_nt(k) < 1.0E19)then ! missing value = 1.0E20
                   arr2d(i,j) = var_nt(k)
                   WTA(i,j) = 1.d0
               endif      
            enddo

            ! assign to full 360x180 matrices from 360x150
            do j = 1,lat_max
               gswp2d(:,j) = arr2d(:,j)
               gswpWTA(:,j) = WTA(:,j)
            enddo

            ! flip matrices to align with GISS input
            do j = 1,lat_maxfull
              giss2d(:,j) = gswp2d(:,lat_maxfull-j+1)
              gissWTA(:,j) = gswpWTA(:,lat_maxfull-j+1)
            enddo

            !**** Interpolate observations to GCM spatial resolution
            call HNTRP0(long_max,lat_maxfull,OFFIA,DIVJA,
     &                          im,jm,OFFIB,DIVJB,SKIP)
            call HNTRPP (gissWTA,giss2d,OUT)

            do i = 1,im
               do j = 1,jm
                 gswp2X2_5_all(i,j,kk) = OUT(i,j)
               enddo
            enddo

            deallocate(arrpack,STAT=status_alloc1)
         else
            print *, 'allocate failed'

         endif allocate_ok

      enddo loop_gswp_variables !kk

!     WRITE FORCINGS (note: variable order is different for GSWP & GISS)
      ! Incoming shortwave radiation [W/m2]
       write(iu_vector(1)) gswp2X2_5_all(:,:,1)
      ! Incoming longwave radiation [W/m2]
       write(iu_vector(2)) gswp2X2_5_all(:,:,2)
      ! Surface air temperature [K]
       write(iu_vector(3)) gswp2X2_5_all(:,:,7)
      ! Surface air moisture [kg/kg]
       write(iu_vector(4)) gswp2X2_5_all(:,:,5)
      ! Surface pressure [Pa]
       write(iu_vector(5)) gswp2X2_5_all(:,:,8)
      ! Wind speed [m/s]
       write(iu_vector(6)) gswp2X2_5_all(:,:,9)
      ! Rainfall rate (total) [kg/m2/s]
       write(iu_vector(7)) gswp2X2_5_all(:,:,3)
      ! Rainfall rate [kg/m2/s]
       write(iu_vector(8)) gswp2X2_5_all(:,:,6)
      ! Convective rainfall rate [kg/m2/s]
       write(iu_vector(9)) gswp2X2_5_all(:,:,4)

      end subroutine avg_gswp_to_giss

!-------------------------------------------------------------------

      subroutine open_outfiles(n_monthyr,iu_vector)

      integer,intent(in) :: n_monthyr ! simulation month number

      integer, parameter :: num_months = 162
      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer, parameter :: num_outfiles = num_months*varnums
      integer,dimension(varnums) :: iu_vector

      integer :: i
      character*80 :: infile
      character*80 :: monthyearBIN(num_months)
      character*80 :: lsm_vars(varnums),basefold_out
      character*80 :: lsm_fold(varnums)

      DATA lsm_fold/ 'force_srheat/'
     &        ,'force_trheat/', 'force_ts/', 'force_qs/'
     &        ,'force_ps/','force_ws/'
     &        ,'force_rainf/', 'force_rainf_c/', 'force_snowf/'/


      DATA lsm_vars/ 'force_srheat'
     &        ,'force_trheat', 'force_ts', 'force_qs'
     &        ,'force_ps','force_ws'
     &        ,'force_rainf', 'force_rainf_c', 'force_snowf'/


      DATA monthyearBIN/       '198207.bi','198208.bi',
     & '198209.bi','198210.bi','198211.bi','198212.bi',
     & '198301.bi','198302.bi','198303.bi','198304.bi',
     & '198305.bi','198306.bi','198307.bi','198308.bi',
     & '198309.bi','198310.bi','198311.bi','198312.bi',
     & '198401.bi','198402.bi','198403.bi','198404.bi',
     & '198405.bi','198406.bi','198407.bi','198408.bi',
     & '198409.bi','198410.bi','198411.bi','198412.bi',
     & '198501.bi','198502.bi','198503.bi','198504.bi',
     & '198505.bi','198506.bi','198507.bi','198508.bi',
     & '198509.bi','198510.bi','198511.bi','198512.bi',
     & '198601.bi','198602.bi','198603.bi','198604.bi',
     & '198605.bi','198606.bi','198607.bi','198608.bi',
     & '198609.bi','198610.bi','198611.bi','198612.bi',
     & '198701.bi','198702.bi','198703.bi','198704.bi',
     & '198705.bi','198706.bi','198707.bi','198708.bi',
     & '198709.bi','198710.bi','198711.bi','198712.bi',
     & '198801.bi','198802.bi','198803.bi','198804.bi',
     & '198805.bi','198806.bi','198807.bi','198808.bi',
     & '198809.bi','198810.bi','198811.bi','198812.bi',
     & '198901.bi','198902.bi','198903.bi','198904.bi',
     & '198905.bi','198906.bi','198907.bi','198908.bi',
     & '198909.bi','198910.bi','198911.bi','198912.bi',
     & '199001.bi','199002.bi','199003.bi','199004.bi',
     & '199005.bi','199006.bi','199007.bi','199008.bi',
     & '199009.bi','199010.bi','199011.bi','199012.bi',
     & '199101.bi','199102.bi','199103.bi','199104.bi',
     & '199105.bi','199106.bi','199107.bi','199108.bi',
     & '199109.bi','199110.bi','199111.bi','199112.bi',
     & '199201.bi','199202.bi','199203.bi','199204.bi',
     & '199205.bi','199206.bi','199207.bi','199208.bi',
     & '199209.bi','199210.bi','199211.bi','199212.bi',
     & '199301.bi','199302.bi','199303.bi','199304.bi',
     & '199305.bi','199306.bi','199307.bi','199308.bi',
     & '199309.bi','199310.bi','199311.bi','199312.bi',
     & '199401.bi','199402.bi','199403.bi','199404.bi',
     & '199405.bi','199406.bi','199407.bi','199408.bi',
     & '199409.bi','199410.bi','199411.bi','199412.bi',
     & '199501.bi','199502.bi','199503.bi','199504.bi',
     & '199505.bi','199506.bi','199507.bi','199508.bi',
     & '199509.bi','199510.bi','199511.bi','199512.bi' /

      basefold_out = '/discover/nobackup/mjpuma/lsm_gswp_2x2_5/'

      do i = 1,varnums
         iu_vector(i) = 20 + (n_monthyr-1)*varnums + i
      enddo

!     open forcings file
      do i = 1,varnums
          infile=basefold_out(1:len_trim(basefold_out))//
     &      lsm_fold(i)(1:len_trim(lsm_fold(i)))//
     &      lsm_vars(i)(1:len_trim(lsm_vars(i)))//
     &      monthyearBIN(n_monthyr)(1:len_trim(monthyearBIN(n_monthyr)))
          print *, infile
          open(iu_vector(i),FILE=infile, STATUS='UNKNOWN',
     &              FORM='UNFORMATTED')
          print *, 'iu_vector(i)=',iu_vector(i)
      enddo
      
      end subroutine open_outfiles



      subroutine close_outfiles(iu_vector)

      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer,dimension(varnums) :: iu_vector
      integer i
      print *, 'in close subroutine'
      do i = 1,varnums
         close(iu_vector(i))
      enddo
      
      end subroutine close_outfiles

      end module drv_main_gswp


!======================================================================
! Purpose: Read GSWP2 data                                            =
!======================================================================

      program gswp_time
      use drv_main_gswp

      implicit none
      include 'netcdf.inc'
      integer :: i,j,nt

      integer, parameter :: im=144, jm=90
      integer, parameter :: im_1deg=360, jm_1deg=180
      integer, parameter :: num_months = 162
      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer,dimension(varnums) :: iu_vector

      INTEGER, DIMENSION(num_months) :: monlen= ! # of 3-hourly timesteps
     &        (/                         248,248,240,248,240,248,   !1982
     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1983
     &           248,232,248,240,248,240,248,248,240,248,240,248,   !1984
     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1985
     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1986
     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1987
     &           248,232,248,240,248,240,248,248,240,248,240,248,   !1988
     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1989
     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1990
     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1991
     &           248,232,248,240,248,240,248,248,240,248,240,248,   !1992
     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1993
     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1994
     &           248,224,248,240,248,240,248,248,240,248,240,248 /) !1995

      loop_months: do i = 1,num_months

         call open_outfiles(i,iu_vector)

         loop_hours:  do nt = 1,monlen(i)
            print *, 'month =',  i, 'timestep = ' , nt
            call avg_gswp_to_giss(i, nt, monlen(i), iu_vector)
            

         enddo loop_hours

         call close_outfiles(iu_vector)

      enddo loop_months


      end program gswp_time
