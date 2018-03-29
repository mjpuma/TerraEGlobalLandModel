      module drv_main_gswp
      implicit none
      save

      real*8,parameter :: PI = 3.1415926535897932d0 !@param pi    pi

      contains

!----------------------------------------------------------------------

      SUBROUTINE drv_finterp(ip,trp_flag,val0,val1,val2,val3,
     &                       madtt,valtrp)
!======================================================================
! Description: INTERPOLATION OF ATMOSPHERIC FORCING DATA
!
! $Log: drv_finterp.f90,v $
! Revision 1.2  2002/09/20  20:55:32  dirmeyer
! Introduced a new interpolation type "P" which specifies a temporal
! disaggregation of interpolated precipitation to improve partitioning
! of infiltration/runoff when forcing data interval is large compared to
! the typical length of a convective (or total) rain event.
!
! Revision 1.1  2002/08/16  15:55:18  guo
! Initial revision
!
!======================================================================
 
      !** This routine interpolates atmospheric adtt seconds period forcing
      !   data to dtt seconds valtrp time-step for one adtt 
      !   interval. Interpolation is performed based on the value of 
      !   'trp_flag':
      !
      !   "L" or "l" = value represents average over interval ending at current time
      !   "N" or "n" = value represents average over interval beginning at current time
      !   "C" or "c" = value represents average over interval centered on current time
      !   "I" or "i" = instantaneous value at current time (linear interpolation)
      !   "P" or "p" = PDF disaggregation applied in time (for precip)
      !   "0" (zero) = no interpolation, centered on current time
      !   Otherwise  = no interpolation, applied beginning at current time
      !
      !
      !
      !************************************************************************
      !* NOTE!!!!!!!
      !*   For "L", "N", and "C" to conserve the mean over the interval after 
      !*   interpolation, 'madtt' MUST BE A MULTIPLE OF 2!
      !************************************************************************
      !
      !   Algorithm for conserving interpolation scheme - courtesy of 
      !   Dag Lohmann, NOAA/NCEP/EMC
      !
      !   - P. Dirmeyer, November 2001
      !===================================================================
 
      IMPLICIT NONE
  
! Input 
      INTEGER,INTENT(IN)           :: ip       ! Grid point in question
      CHARACTER(len=1),INTENT(IN)  :: trp_flag ! Type of interpolation

      REAL,INTENT(IN)              :: val0   ! last forcing value (for cases N & C)
      REAL,INTENT(IN)              :: val1   ! current forcing value (all cases)
      REAL,INTENT(IN)              :: val2   ! next forcing value (all cases but default)
      REAL,INTENT(IN)              :: val3   ! ueber-next forcing value (for cases C & L)

      INTEGER,INTENT(IN)           :: madtt  ! Number of valtrp timesteps in
                                             ! a forcing timestep.
! Output
      REAL, DIMENSION(madtt)       :: valtrp ! Interpolated forcing data vector

! Local
      REAL               :: fac0 ! Weight for last value 
      REAL               :: fac1 ! Weight for current value
      REAL               :: fac2 ! Weight for next value
      REAL               :: rmadtt ! real madtt
      REAL               :: denom  ! denominator of scaling factor for 
                                   ! conserving interpolation
      REAL               :: numer  ! numerator of scaling factor for 
                                   ! conserving interpolation
      INTEGER            :: i
      INTEGER            :: j
      INTEGER            :: ip1  ! i + 1
      INTEGER            :: im1  ! i - 1

      REAL               :: rtdist(120) ! Weights for precip disag PDF
      REAL               :: factor      ! Scaling factor for brevity/intensity
      REAL               :: expo        ! Exponent to fit slope of log-log relationship
      REAL               :: rtsum       ! Ensure PDF weights add to 1.0
      REAL               :: rtsum0      ! Sum from last time interval
      REAL               :: p0, p1      ! Define edges of each time bin
      REAL               :: rintr1      ! Log-log tail target value

      LOGICAL            :: iniflag
      SAVE iniflag,rtdist
      DATA iniflag/.true./


      !===================================================================
      rmadtt = float(madtt)

!
!>>> Generate disaggregation PDF - Just do this once
!
      IF (iniflag) THEN
         write(6,4030)
 4030    format('Initializing precipitation disaggregation PDF:')
         factor = 66.6666  ! For 33% dry time steps in a rainy forcing interval
!         factor = 29.63  ! For 0% dry time steps in a rainy forcing interval
         expo = -3.0      
         rtsum = 0.0
         rintr1 = factor * float(madtt) ** expo

         DO j = 1, madtt
           rtsum0 = rtsum
           p0 = float(j-1)/float(madtt)
           p1 = float(j)/float(madtt)
           rtdist(j) = rintr1*expo/((expo+1)*(p1-p0)*100)*
     &                  ((100*p1/rintr1)**((expo+1)/expo)-
     &                   (100*p0/rintr1)**((expo+1)/expo))
           rtsum = rtsum + rtdist(j)
           IF (rtsum > 1.0) THEN
             IF (rtsum0 < 1.0) THEN
               rtdist(j) = rtdist(j) + 1.0 - rtsum   ! Stick any remaining weight in last interval
             ELSE
               rtdist(j) = 0.0
             ENDIF
           ENDIF
           write(6,4050) j,rtdist(j)
 4050      format('Time disag: ',i3,': ',f11.6)
         ENDDO
         IF (rtsum < 1.0) THEN
             rtdist(1) = rtdist(1) + 1.0 - rtsum  ! Ensure integral of PDF = 1.0
         ENDIF
      ENDIF
      iniflag = .false.

      IF (trp_flag == 'I' .or. trp_flag == 'i') THEN
      !** Current value valid at midpoint of interval
        DO j = 1, madtt
          fac1     = float(madtt+1-j)/madtt
          fac2     = 1.0-fac1
          valtrp(j)   = val1*fac1+val2*fac2
        ENDDO

      ELSEIF (trp_flag == 'N' .or. trp_flag == 'n') THEN
      !** Current value is average over next interval
        IF (mod(madtt,2) /= 0) STOP 'drv_interp mod(madtt,2) /= 0'
        DO j = 1, madtt
          fac1 = (2.0*rmadtt-abs(float(2*j-madtt-1)))/(rmadtt*2.0)
          fac0 = max(1.0-float(j*2+madtt-1)/(rmadtt*2.0),0.0)
          fac2 = max(1.0-float((madtt+1-j)*2+madtt-1)/(rmadtt*2.0),0.0)
          denom = 0.5*(val0+val2)+3.0*val1
          numer = 4.0*val1
          IF (denom > EPSILON(denom)) THEN
             valtrp(j) = (val0*fac0+val1*fac1+val2*fac2) * numer / denom
          ELSE
             valtrp(j) = 0.0
          ENDIF
        ENDDO

      ELSEIF (trp_flag == 'C' .or. trp_flag == 'c') THEN
      !** Current value is average centered on current time
        IF (mod(madtt,2) /= 0) STOP 'drv_interp mod(madtt,2) /= 0'
        DO j = 1, madtt
          fac1 = (2.0*rmadtt-(float(2*j-1)))/(rmadtt*2.0)
          fac0 = 0.0
          fac2 = 1.0-fac1
          IF (j > madtt/2) THEN
            denom = 0.5*(val1+val3)+3.0*val2
            numer = 4.0*val2
          ELSE
            denom = 0.5*(val0+val2)+3.0*val1
            numer = 4.0*val1
          ENDIF
          IF (denom > EPSILON(denom)) THEN
             valtrp(j) = (val0*fac0+val1*fac1+val2*fac2) * numer / denom
          ELSE
             valtrp(j) = 0.0
          ENDIF
        ENDDO

      ELSEIF (trp_flag == 'L' .or. trp_flag == 'l') THEN
      !** Current value is average for period ending at current time
        IF (mod(madtt,2) /= 0) STOP 'drv_interp mod(madtt,2) /= 0'
        DO j = 1, madtt
          fac1 = (2.0*rmadtt-abs(float(2*j-madtt-1)))/(rmadtt*2.0)
          fac0 = max(1.0-float(j*2+madtt-1)/(rmadtt*2.0),0.0)
          fac2 = max(1.0-float((madtt+1-j)*2+madtt-1)/(rmadtt*2.0),0.0)
          denom = 0.5*(val1+val3)+3.0*val2
          numer = 4.0*val2
          IF (denom > EPSILON(denom)) THEN
             valtrp(j) = (val1*fac0+val2*fac1+val3*fac2) * numer / denom
          ELSE
             valtrp(j) = 0.0
          ENDIF
        ENDDO

      ELSEIF (trp_flag == 'P' .or. trp_flag == 'p') THEN
      !** Current value is from a PDF that disagregates in time over the forcing interval
        DO j = 1, madtt
          valtrp(j) = MAX(val2*rtdist(j)*float(madtt),0.0)
        ENDDO

      ELSEIF (trp_flag == '0') THEN
      !** Current value is applied centered on current time without interpolation
        DO j = 1, madtt
          IF (j > madtt/2) THEN
            valtrp(j) = val2
          ELSE
            valtrp(j) = val1
          ENDIF
        ENDDO

      ELSE
      !** Current value is applied over this interval without interpolation
        DO j = 1, madtt
          valtrp(j) = val1
        ENDDO

      ENDIF

      END SUBROUTINE drv_finterp

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

      function LEAPYR(year)
      integer :: year
      logical :: LEAPYR
      LEAPYR = .false.
      if(mod(year,400) == 0) then
         LEAPYR = .true.
      else
         if(mod(year,4) == 0 .and. mod(year,100) /= 0) then
           LEAPYR = .true.
         endif
      endif
      end function LEAPYR

!----------------------------------------------------------------------

      function solarelevangle(lat,td,th) Result(selev) 
      implicit none
      real*8 :: lat !Latitude (radians)
      real*8 :: td  !Julian day time (includes sub-day fraction of time
      real*8 :: th  !Local hour of day
      real*8 :: selev !Solar elevation angle (radians)
      !----Local-----
      real*8 :: delta
      real*8 :: sinselev !Sine of solar elev angle

      delta=asin(-sin(pi/180.0*23.45)*cos(2*pi*(td+10.0)/365.0))
      sinselev=sin(lat)*sin(delta)+cos(lat)*cos(delta)*
     &     cos(pi/180.0*15.0*(th-12.0)) 
      selev = asin(sinselev)

      end function solarelevangle

!----------------------------------------------------------------------

      function fdiffuse(solarelev,day,Rg) Result(fdif)
      implicit none
      real*8,intent(in) :: solarelev !Solar zenith angle (radians)
      real*8,intent(in) :: day   !Day of year
      real*8,intent(in) :: Rg  !Global incident radiation (W m-2)
      real*8 :: fdif
      !---Local------
      real*8,parameter :: SOLARCONST = 1370.0d0  !Solar constant, (W m-2)
      real*8 :: sbeta !Sine of solar elevation angle
      real*8 :: S0    !Incident rad at top-of-atmosphere at time (W m-2)
      real*8 :: transm !Transmissivity through the atmosphere (fraction)
      real*8 :: R, K   

      sbeta = sin(solarelev)
      
      if(sbeta.gt.0.0d0)then
        S0=SOLARCONST*(1.0d0+0.033d0*cos(2.0d0*pi*day/365.0d0))*sbeta
        !Rg=I0/2.3
        transm=Rg/S0
        R=0.847d0-1.61d0*sbeta+1.04d0*(sbeta**2.d0)
        K=(1.47d0-R)/1.66d0
        if(transm.le.0.22d0)then
          fdif=1.0d0
        elseif(transm.le.0.35)then
          fdif=1.0d0-6.4d0*((transm-0.22d0)**2.d0)
        elseif(transm.le.K)then
          fdif=1.47d0-1.66d0*transm
        else
          fdif=R
        endif
      else
        fdif = 1.d0
      endif
      end function fdiffuse

!----------------------------------------------------------------------

!      SUBROUTINE calc_Tcp1(zs,z0,disp,h,Ta, P_Pa, E, Fh, U,Ustar,RH, 
!     $     rho, Ch, Tcp, TcpC,Qsurf)
!      !* From Andrew Friend's physics.f

!      implicit none
!      !* Input parameters *!
!      real*8,intent(in) :: zs,z0,disp,h,Ta,P_Pa,E,Fh,U,Ustar,RH

      !* Outputs *!
!      real*8, intent(out) :: rho       !air density (kg/m3) 
!      real*8,intent(out) :: Ch         !heat transfer coefficient
!      real*8,intent(out) :: Tcp, TcpC  !Potential temperature of canopy, 
                                       !Kelvin and Celsius
!      real*8,intent(out) :: Qsurf      !canopy air specific humidity

      !* Local variables *!
!      real*8, parameter :: tfrz=273.15d0  !Kelvin
!      real*8, parameter :: rk = 0.41d0 !.41 !von Karman
!      real*8, parameter :: cp=1012.d0 !Heat capacity of dry air (J kg-1 K-1) 
!      real*8, parameter :: R = 8.3144d0 !Ideal gas constant (J mol-1 K-1)
!      real*8, parameter :: grav = 9.8067d0  !gravitational acceleration, m s-2
!      real*8 :: Ts              !Temperature at ?? height, sensor or surface?
!      real*8 :: Tsp, Tgp  !Pot. temperature at the surface layer and ground
!      real*8 :: delt, err, func
!      integer :: iChc
!      integer :: J
!      real*8 :: Ris  !bulk Richardson number for surface layer
!      real*8 :: DM   !variation of drag coeff for non-neutral stability
!      real*8 :: CD   !Drag coefficient
!      real*8 :: CDN  !Drag coefficient for neutral stability, k=0.35

    !----------------------------------------------------------------------!
!      real*8 :: a(10),b(10),c(10),d(10),f(10)
!      data a / 16.60   , 10.40   , 5.24    , 3.13    , 2.32    ,
!     &         1.98    , 1.83    , 1.76    , 1.71    , 1.69      /
!      data b / 3.25    , 8.45E-01, 1.83E-01, 4.62E-02, 1.25E-02,
!     &         3.20E-03, 7.28E-04, 1.47E-04, 2.71E-05, 4.60E-06  /
!      data c / 5.11    , 1.68    , 0.56    , 0.33    , 0.30    ,
!     &         0.33    , 0.39    , 0.45    , 0.51    , 0.56      /
!      data d / 1.24    , 0.81    , 0.66    , 0.80    , 0.99    ,
!     &         1.18    , 1.35    , 1.50    , 1.65    , 1.78      /
!      data f / 0.13    , 0.14    , 0.24    , 0.44    , 0.66    ,
!     &         0.87    , 1.05    , 1.22    , 1.37    , 1.52     /
    !----------------------------------------------------------------------!

!      Ts=Ta !+tfrz  !Ta is already passed in K -PK 11/19/06
!      Tsp=Ts*(101325.d0/P_Pa)**(R/cp)
!      Qsurf=(RH/100.d0)*QSAT(Tsp,(2500800.d0-2360.d0*
!     & (Ts-tfrz)),P_Pa/100.d0)
!!      write(*,*) Qsurf, RH, Tsp, Ts, P_Pa
!      rho=28.964d-3*P_Pa/(R*Ts)
!      iChc=anint(log10((zs-disp)/z0))  !?What's this from?
!      CDN=(rk**2.d0)/((alog((zs-disp)/z0))**2.d0)
!      delt=-2.d0
!      Tgp=Tsp+delt
!      err=0.001d0
!      do J=1,100
!        Ris=((zs-disp)*grav*(Tsp-Tgp))/(Tgp*(U**2))
!        if(Ris.lt.0.d0)then  !unstable case
!          DM=sqrt(((1.d0-a(iChc)*Ris)*(1.d0-b(iChc)*Ris))/ 
!     &         (1.d0-c(iChc)*Ris))
!          CD=CDN*DM
!          Ch=CD*1.35d0*sqrt((1.d0-d(iChc)*Ris)/(1.d0-f(iChc)*Ris))
!        else               !stable case
!          DM=1.d0/(1.d0+(11.2d0+90.d0*Ris)*Ris)
!          CD=CDN*DM
!          Ch=CD*1.35d0/(1.d0+1.93d0*Ris)
!        endif
!        func=abs(Tgp-(Tsp+Fh/(cp*rho*Ch*U)))
!        if(func>err)then
!          delt=-0.5d0*delt
!        endif
!        err=func
!        Tgp=Tgp+delt
!      enddo
!      Tcp=Tgp/((101325.d0/P_Pa)**(R/cp))  !Put back to actual temperature of canopy
!      TcpC=Tcp-tfrz
!      if((TcpC.gt.-99.d0).and.(TcpC.lt.99.d0))then
!        TcpC=TcpC
!      else
!        TcpC=-99.d0
!      endif
!
!      return  !Tcp, Ch, TcpC
!      END SUBROUTINE calc_Tcp1

!======================================================================!

      FUNCTION QSAT (TM,QL,PR) Result(Qsatval)
      implicit none
!@sum  QSAT calculates saturation vapour mixing ratio (kg/kg)
!@auth Gary Russell
!@ver  1.0
!      USE CONSTANT, only : mrat,rvap,tf
!      IMPLICIT NONE
!@var A,B,C   expansion coefficients for QSAT
      real*8 :: Qsatval
      REAL*8, PARAMETER :: A=3.797915e0    !3.797915d0
      REAL*8, PARAMETER :: B=7.93252e-6    !7.93252d-6
      REAL*8, PARAMETER :: C=2.166847e-3         !2.166847d-3
      real*8 :: TM, QL, PR
!**** Note that if QL is considered to be a function of temperature, the
!**** correct argument in QSAT is the average QL from t=0 (C) to TM, ie.
!**** QL = 0.5*(QL(0)+QL(t))
!      REAL*8, INTENT(IN) :: TM  !@var TM   potential temperature (K)
!      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
!     REAL*8, INTENT(IN) :: QL  !@var QL   lat. heat of vap./sub. (J/kg)
!      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
      QSATval = A*EXP(QL*(B-C/TM))/PR
      RETURN
      END function QSAT


!======================================================================!

      subroutine get_gswp_forcings(
     &         force_srheat          ,
     &         force_trheat          ,
     &         force_ts              ,
     &         force_qs              ,
     &         force_ps              ,
     &         force_ws              ,
     &         force_rainf           ,
     &         force_rainf_c         ,
     &         force_snowf   )


      integer, parameter :: im=72, jm=46
      integer, parameter :: num_months = 18  !162
      integer, parameter :: varnums = 9 ! number of GSWP2 variables

      real*4, dimension(im,jm) ::
     &         force_srheat         ,
     &         force_trheat         ,
     &         force_ts             ,
     &         force_qs             ,
     &         force_ps             ,
     &         force_ws             ,
     &         force_rainf          ,
     &         force_rainf_c        ,
     &         force_snowf 


      integer,dimension(varnums) :: iu_vector
      character*80 :: varname, infile
      character*80 :: basefolder, foldnames(varnums),vnames(varnums)
      character*80 :: v2names(varnums), monthyear(num_months)
      integer i,j,k


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

      DATA monthyear/          '198207.bi','198208.bi',
     & '198209.bi','198210.bi','198211.bi','198212.bi',
     & '198301.bi','198302.bi','198303.bi','198304.bi',
     & '198305.bi','198306.bi','198307.bi','198308.bi',
     & '198309.bi','198310.bi','198311.bi','198312.bi'/!,
!     & '198401.bi','198402.bi','198403.bi','198404.bi',
!     & '198405.bi','198406.bi','198407.bi','198408.bi',
!     & '198409.bi','198410.bi','198411.bi','198412.bi',
!     & '198501.bi','198502.bi','198503.bi','198504.bi',
!     & '198505.bi','198506.bi','198507.bi','198508.bi',
!     & '198509.bi','198510.bi','198511.bi','198512.bi',
!     & '198601.bi','198602.bi','198603.bi','198604.bi',
!     & '198605.bi','198606.bi','198607.bi','198608.bi',
!     & '198609.bi','198610.bi','198611.bi','198612.bi',
!     & '198701.bi','198702.bi','198703.bi','198704.bi',
!     & '198705.bi','198706.bi','198707.bi','198708.bi',
!     & '198709.bi','198710.bi','198711.bi','198712.bi',
!     & '198801.bi','198802.bi','198803.bi','198804.bi',
!     & '198805.bi','198806.bi','198807.bi','198808.bi',
!     & '198809.bi','198810.bi','198811.bi','198812.bi',
!     & '198901.bi','198902.bi','198903.bi','198904.bi',
!     & '198905.bi','198906.bi','198907.bi','198908.bi',
!     & '198909.bi','198910.bi','198911.bi','198912.bi',
!     & '199001.bi','199002.bi','199003.bi','199004.bi',
!     & '199005.bi','199006.bi','199007.bi','199008.bi',
!     & '199009.bi','199010.bi','199011.bi','199012.bi',
!     & '199101.bi','199102.bi','199103.bi','199104.bi',
!     & '199105.bi','199106.bi','199107.bi','199108.bi',
!     & '199109.bi','199110.bi','199111.bi','199112.bi',
!     & '199201.bi','199202.bi','199203.bi','199204.bi',
!     & '199205.bi','199206.bi','199207.bi','199208.bi',
!     & '199209.bi','199210.bi','199211.bi','199212.bi',
!     & '199301.bi','199302.bi','199303.bi','199304.bi',
!     & '199305.bi','199306.bi','199307.bi','199308.bi',
!     & '199309.bi','199310.bi','199311.bi','199312.bi',
!     & '199401.bi','199402.bi','199403.bi','199404.bi',
!     & '199405.bi','199406.bi','199407.bi','199408.bi',
!     & '199409.bi','199410.bi','199411.bi','199412.bi',
!     & '199501.bi','199502.bi','199503.bi','199504.bi',
!     & '199505.bi','199506.bi','199507.bi','199508.bi',
!     & '199509.bi','199510.bi','199511.bi','199512.bi' /

      basefolder = '/archive/g03/mjpuma/GSWP/'


      loop_gswp_variables: do k = 1, varnums

         ! open file
         infile=basefolder(1:len_trim(basefolder))//
     &       foldnames(k)(1:len_trim(foldnames(k)))//
     &       vnames(k)(1:len_trim(vnames(k)))//
     &       monthyear(n_monthyr)(1:len_trim(monthyear(n_monthyr)))

         open(iu_vector(i),FILE=infile, STATUS='UNKNOWN',
     &              FORM='UNFORMATTED')


      enddo loop_gswp_variables !k

      print *,"reading forcings"

      read(951)
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_ws (1:im,1:jm)             ,
     &         force_rainf (1:im,1:jm)   ,
     &         force_rainf_c (1:im,1:jm),
     &         force_snowf (1:im,1:jm)

      print *,"forcings ok"


      end subroutine get_gswp_forcings

!======================================================================!


      subroutine assign_gswp_forcings(
     i         force_srheat (1:im,1:jm)         ,
     i         force_trheat (1:im,1:jm)         ,
     i         force_ts (1:im,1:jm)             ,
     i         force_qs (1:im,1:jm)             ,
     i         force_ps (1:im,1:jm)             ,
     i         force_ws (1:im,1:jm)             ,
     i         force_rainf (1:im,1:jm)          ,
     i         force_rainf_c (1:im,1:jm)        ,
     i         force_snowf (1:im,1:jm)          ,
     o         force_Ca (1:im,1:jm)             ,
     o         force_cos_zen_angle (1:im,1:jm)  ,
     o         force_vis_rad (1:im,1:jm)        ,
     o         force_direct_vis_rad (1:im,1:jm) ,
     o         force_prec_ms (1:im,1:jm)        ,
     o         force_eprec_w (1:im,1:jm)        ,
     o         force_sprec_ms (1:im,1:jm)       ,
     o         force_seprec_w (1:im,1:jm)       ,
     o         force_rhosrf (1:im,1:jm)         ,
     o         force_cdh (1:im,1:jm)            ,
     o         force_qm1 (1:im,1:jm)            ,
     o         force_pbl_args_ws0 (1:im,1:jm)   ,
     o         force_pbl_args_tprime (1:im,1:jm),
     o         force_pbl_args_qprime (1:im,1:jm) )


      do i=1,im
         do j = 1,jm
            ! Cosine of solar zenith angle 
            force_cos_zen_angle(i,j)= 1.d0
           ! Visible radiation [W/m2]
            force_vis_rad(i,j)= 1.d0
           ! Direct visible radiation [W/m2]
            force_direct_vis_rad(i,j)= 1.d0
           ! Total precip. (rain+snow) both convective & large-scale [m/s]
            force_prec_ms(i,j)= (puma(i,j,3)+puma(i,j,4)+puma(i,j,6))
     &                         / 1000.d0
           ! Energy of precip.[W/m2]: 0 for rain; lhm [units?]*prec for snow
            force_eprec_w(i,j)= 1.d0
            ! Large-scale precipitation [m/s]: currently set to 0 in modelE runs
            force_sprec_ms(i,j)= 1.d0
            ! Energy of large-scale precipitation [W/m2]
            force_seprec_w(i,j)= 1.d0
            ! Incoming shortwave radiation [W/m2]
            force_srheat(i,j)= puma(i,j,1)
            ! Incoming longwave radiation [W/m2]
            force_trheat(i,j)= puma(i,j,2)
            ! Surface air temperature [???]
            force_ts(i,j)= puma(i,j,7)
            ! Surface air moisture [???]
            force_qs(i,j)= puma(i,j,5)
            ! Surface pressure [???]
            force_ps(i,j)= puma(i,j,8)
            ! Surface air density [???]
            force_rhosrf(i,j)= 1.d0
            ! Turbulent transfer coefficient [???}
            force_cdh(i,j)= 1.d0
            ! Amount of water in the 1st atm layer [???]
            force_qm1(i,j)= 1.d0
            ! Wind speed [m/s]
            force_ws(i,j)= puma(i,j,9)
            ! Rainfall rate [kg/m2/s]
            force_rainf(i,j) = puma(i,j,3) 
            ! Convective rainfall rate [kg/m2/s]            
            force_rainf_c(i,j) = puma(i,j,6)
            ! Snowfall rate [kg/m2/s]
            force_snowf(i,j) = puma(i,j,4)

            !     Assign non-GSWP2 inputs
            force_Ca(i,j) = 0.0127609 !(350 ppm @STP in mol m-3)
            force_pbl_args_ws0 (i,j) = 0.d0
            force_pbl_args_tprime (i,j) = 0.d0
            force_pbl_args_qprime (i,j) = 0.d0


         enddo
      enddo

      end subroutine assign_gswp_forcings

!-------------------------------------------------------------------

      subroutine open_outfiles(n_monthyr,iu_vector)

      integer,intent(in) :: n_monthyr ! simulation month number

      integer, parameter :: num_months = 18  !162
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
     & '198309.bi','198310.bi','198311.bi','198312.bi'/!,
!     & '198401.bi','198402.bi','198403.bi','198404.bi',
!     & '198405.bi','198406.bi','198407.bi','198408.bi',
!     & '198409.bi','198410.bi','198411.bi','198412.bi',
!     & '198501.bi','198502.bi','198503.bi','198504.bi',
!     & '198505.bi','198506.bi','198507.bi','198508.bi',
!     & '198509.bi','198510.bi','198511.bi','198512.bi',
!     & '198601.bi','198602.bi','198603.bi','198604.bi',
!     & '198605.bi','198606.bi','198607.bi','198608.bi',
!     & '198609.bi','198610.bi','198611.bi','198612.bi',
!     & '198701.bi','198702.bi','198703.bi','198704.bi',
!     & '198705.bi','198706.bi','198707.bi','198708.bi',
!     & '198709.bi','198710.bi','198711.bi','198712.bi',
!     & '198801.bi','198802.bi','198803.bi','198804.bi',
!     & '198805.bi','198806.bi','198807.bi','198808.bi',
!     & '198809.bi','198810.bi','198811.bi','198812.bi',
!     & '198901.bi','198902.bi','198903.bi','198904.bi',
!     & '198905.bi','198906.bi','198907.bi','198908.bi',
!     & '198909.bi','198910.bi','198911.bi','198912.bi',
!     & '199001.bi','199002.bi','199003.bi','199004.bi',
!     & '199005.bi','199006.bi','199007.bi','199008.bi',
!     & '199009.bi','199010.bi','199011.bi','199012.bi',
!     & '199101.bi','199102.bi','199103.bi','199104.bi',
!     & '199105.bi','199106.bi','199107.bi','199108.bi',
!     & '199109.bi','199110.bi','199111.bi','199112.bi',
!     & '199201.bi','199202.bi','199203.bi','199204.bi',
!     & '199205.bi','199206.bi','199207.bi','199208.bi',
!     & '199209.bi','199210.bi','199211.bi','199212.bi',
!     & '199301.bi','199302.bi','199303.bi','199304.bi',
!     & '199305.bi','199306.bi','199307.bi','199308.bi',
!     & '199309.bi','199310.bi','199311.bi','199312.bi',
!     & '199401.bi','199402.bi','199403.bi','199404.bi',
!     & '199405.bi','199406.bi','199407.bi','199408.bi',
!     & '199409.bi','199410.bi','199411.bi','199412.bi',
!     & '199501.bi','199502.bi','199503.bi','199504.bi',
!     & '199505.bi','199506.bi','199507.bi','199508.bi',
!     & '199509.bi','199510.bi','199511.bi','199512.bi' /

      basefold_out = '/archive/g03/mjpuma/lsm_gswp/'

      do i = 1,varnums
         iu_vector(i) = 30 + (n_monthyr-1)*varnums + i
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
! Purpose: Read GSWP2 data - put into lsm_standalone.f                =
!======================================================================

      program gswp_time
      use drv_main_gswp

      implicit none
      include 'netcdf.inc'
      integer :: i,nt

      integer, parameter :: im=72, jm=46
      integer, parameter :: num_months = 18  !162
      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer,dimension(varnums) :: iu_vector

      real*4, dimension(im,jm) :: force_Ca
     &        ,force_cos_zen_angle, force_vis_rad
     &        ,force_direct_vis_rad, force_prec_ms 
     &        ,force_eprec_w, force_sprec_ms
     &        ,force_seprec_w, force_srheat
     &        ,force_trheat, force_ts, force_qs
     &        ,force_ps, force_rhosrf, force_cdh
     &        ,force_qm1, force_ws, force_pbl_args_ws0
     &        ,force_pbl_args_tprime, force_pbl_args_qprime


      INTEGER, DIMENSION(num_months) :: monlen= ! # of 3-hourly timesteps
     &        (/                         248,248,240,248,240,248,   !1982
     &           248,224,248,240,248,240,248,248,240,248,240,248 /)!,   !1983
!     &           248,232,248,240,248,240,248,248,240,248,240,248,   !1984
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1985
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1986
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1987
!     &           248,232,248,240,248,240,248,248,240,248,240,248,   !1988
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1989
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1990
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1991
!     &           248,232,248,240,248,240,248,248,240,248,240,248,   !1992
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1993
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1994
!     &           248,224,248,240,248,240,248,248,240,248,240,248 /) !1995


      loop_months: do i = 1,num_months

         call open_outfiles(i,iu_vector)

         loop_hours:  do nt = 1,monlen(i)

            call get_gswp_forcings(
     o         force_srheat          ,
     o         force_trheat          ,
     o         force_ts              ,
     o         force_qs              ,
     o         force_ps              ,
     o         force_ws              ,
     o         force_rainf           ,
     o         force_rainf_c         ,
     o         force_snowf   )

            call assign_gswp_forcings(
     i         force_srheat (1:im,1:jm)         ,
     i         force_trheat (1:im,1:jm)         ,
     i         force_ts (1:im,1:jm)             ,
     i         force_qs (1:im,1:jm)             ,
     i         force_ps (1:im,1:jm)             ,
     i         force_ws (1:im,1:jm)             ,
     i         force_rainf (1:im,1:jm)          ,
     i         force_rainf_c (1:im,1:jm)        ,
     i         force_snowf (1:im,1:jm)          ,
     o         force_Ca (1:im,1:jm)             ,
     o         force_cos_zen_angle (1:im,1:jm)  ,
     o         force_vis_rad (1:im,1:jm)        ,
     o         force_direct_vis_rad (1:im,1:jm) ,
     o         force_prec_ms (1:im,1:jm)        ,
     o         force_eprec_w (1:im,1:jm)        ,
     o         force_sprec_ms (1:im,1:jm)       ,
     o         force_seprec_w (1:im,1:jm)       ,
     o         force_rhosrf (1:im,1:jm)         ,
     o         force_cdh (1:im,1:jm)            ,
     o         force_qm1 (1:im,1:jm)            ,
     o         force_pbl_args_ws0 (1:im,1:jm)   ,
     o         force_pbl_args_tprime (1:im,1:jm),
     o         force_pbl_args_qprime (1:im,1:jm) )


            print *, 'month =',  i, 'timestep = ' , nt

         enddo loop_hours

         call close_outfiles(iu_vector)

      enddo loop_months


      end program gswp_time
