      module drv_main_gswp
      implicit none
      save

      real*8,parameter :: PI = 3.1415926535897932d0 !@param pi    pi

      contains

!======================================================================!

      subroutine check_err(iret)
      implicit none
      include 'netcdf.inc'
      integer                :: iret
      if(iret /= nf_noerr) then
        print*, nf_strerror(iret)
        stop
      endif
      end subroutine check_err


!======================================================================!


      subroutine gswp_timeavg(n_monthyr,nt,monlen,iu_vector,
     &                        gswp1X1_all_avg,num_all)

      include 'netcdf.inc'
!
      integer,intent(in) :: n_monthyr ! simulation month number
      integer,intent(in) :: nt ! current month timestep (3-hourly) number
      integer,intent(in) :: monlen ! # of timesteps in current month
!
      integer, parameter :: long_max=360,lat_max=150 !GSWP2 dimensions
      integer, parameter :: lat_maxfull=180 ! # of 1 deg latitude divisions
      real*4, parameter  :: SKIP = -1.E30
!
      integer, parameter :: num_months = 120
      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer, parameter :: num_outfiles = num_months*varnums
      integer, parameter :: num_landcells = 15238
      integer,dimension(varnums+3) :: iu_vector

      character*80 :: varname, infile
      character*80 :: basefolder, foldnames(varnums),vnames(varnums)
      character*80 :: v2names(varnums), monthyear(num_months)

      real*4, dimension(long_max,lat_maxfull) :: giss2d
      real*4, dimension(long_max,lat_max) :: arr2d,lat,long
      real*4, dimension(long_max,lat_maxfull) :: gswp2d
      real*4, dimension(long_max,lat_maxfull,varnums) :: gswp1X1_all
      real*4, dimension(long_max,lat_maxfull,varnums)::gswp1X1_all_avg
      integer,dimension(long_max,lat_maxfull,varnums)::num_all
      ! Local variables
      integer :: fid,varid,status
      integer :: status_alloc1 ! status = 0 if success
      integer i,j,k,kk,land

      integer,dimension(num_landcells) :: xy
      real*4 ,dimension(num_landcells) :: var_nt
      real*4, allocatable, dimension(:,:) :: arrpack

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

      DATA monthyear/ 
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

      basefolder = '/archive/g03/mjpuma/GSWP/'

      ! intialize variables--------
      gswp1X1_all(:,:,:) = 0.d0

      loop_gswp_variables: do kk = 1, varnums

         arr2d(:,:)= SKIP
         giss2d(:,:) = SKIP
         gswp2d(:,:) = SKIP

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
               else
                   arr2d(i,j) = SKIP
               endif      
            enddo

            ! assign to full 360x180 matrix from 360x150
            do j = 1,lat_max
               gswp2d(:,j) = arr2d(:,j)
            enddo

            ! flip matrices to align with GISS input
            do j = 1,lat_maxfull
              giss2d(:,j) = gswp2d(:,lat_maxfull-j+1)
            enddo

            do i = 1,long_max
               do j = 1,lat_maxfull
                 gswp1X1_all(i,j,kk) = giss2d(i,j)

                 if(giss2d(i,j)>-1.e20)then
                   num_all(i,j,kk) = num_all(i,j,kk) + 1
                 endif 

               enddo
            enddo

            deallocate(arrpack,STAT=status_alloc1)
         else
            print *, 'allocate failed'

         endif allocate_ok

      enddo loop_gswp_variables !kk

!     Write the 3-hrly data for these variables
      ! Rainfall rate (total) [kg/m2/s]
      write(iu_vector(10)) gswp1X1_all_avg(:,:,3)
      ! Rainfall rate [kg/m2/s]
      write(iu_vector(11)) gswp1X1_all_avg(:,:,6)
      ! Snowfall rate [kg/m2/s]
      write(iu_vector(12)) gswp1X1_all_avg(:,:,4)


!     Average forcings
      do i = 1,long_max
         do j = 1,lat_maxfull
            do k = 1,9

               if (gswp1X1_all(i,j,k)>-1.e20)then
                  gswp1X1_all_avg(i,j,k)=
     &               gswp1X1_all_avg(i,j,k)+gswp1X1_all(i,j,k)
               endif

            enddo
!      ! Incoming longwave radiation [W/m2]
!      gswp1X1_all_avg(i,j,2)=gswp1X1_all_avg(i,j,2)+gswp1X1_all(i,j,2)
!      ! Surface air temperature [K]
!      gswp1X1_all_avg(i,j,7)=gswp1X1_all_avg(i,j,7)+gswp1X1_all(i,j,7)
!      ! Surface air moisture [kg/kg]
!      gswp1X1_all_avg(i,j,5)=gswp1X1_all_avg(i,j,5)+gswp1X1_all(i,j,5)
!      ! Surface pressure [Pa]
!      gswp1X1_all_avg(i,j,8)=gswp1X1_all_avg(i,j,8)+gswp1X1_all(i,j,8)
!      ! Wind speed [m/s]
!      gswp1X1_all_avg(i,j,9)=gswp1X1_all_avg(i,j,9)+gswp1X1_all(i,j,9)
!      ! Rainfall rate (total) [kg/m2/s]
!           gswp1X1_all_avg(i,j,3)=gswp1X1_all_avg(i,j,3)+gswp1X1_all(i,j,3)
!      ! Convective rainfall rate [kg/m2/s]
!      gswp1X1_all_avg(i,j,6)=gswp1X1_all_avg(i,j,6)+gswp1X1_all(i,j,6)
!      ! Snowfall rate [kg/m2/s]
!      gswp1X1_all_avg(i,j,4)=gswp1X1_all_avg(i,j,4)+gswp1X1_all(i,j,4)

         enddo
      enddo

      end subroutine gswp_timeavg

!======================================================================!

      subroutine write_outday(iu_vector, gswp1X1_all_avg)

      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer, parameter :: long_max=360,lat_max=180 
      integer,dimension(varnums+3) :: iu_vector
      real*4, dimension(long_max,lat_max,varnums)::gswp1X1_all_avg

!     WRITE FORCINGS (note: variable order is different for GSWP & GISS)
      ! Incoming shortwave radiation [W/m2]
      write(iu_vector(1)) gswp1X1_all_avg(:,:,1)
      ! Incoming longwave radiation [W/m2]
      write(iu_vector(2)) gswp1X1_all_avg(:,:,2)
      ! Surface air temperature [K]
      write(iu_vector(3)) gswp1X1_all_avg(:,:,7)
      ! Surface air moisture [kg/kg]
      write(iu_vector(4)) gswp1X1_all_avg(:,:,5)
      ! Surface pressure [Pa]
      write(iu_vector(5)) gswp1X1_all_avg(:,:,8)
      ! Wind speed [m/s]
      write(iu_vector(6)) gswp1X1_all_avg(:,:,9)
      ! Rainfall rate (total) [kg/m2/s]
      write(iu_vector(7)) gswp1X1_all_avg(:,:,3)
      ! Rainfall rate [kg/m2/s]
      write(iu_vector(8)) gswp1X1_all_avg(:,:,6)
      ! Snowfall rate [kg/m2/s]
      write(iu_vector(9)) gswp1X1_all_avg(:,:,4)

      end subroutine write_outday

!======================================================================!

      subroutine avg_sum(gswp1X1_all_avg,num_all)
      implicit none

      integer, parameter :: long_max    = 360
      integer, parameter :: lat_maxfull = 180 ! # of 1 deg latitude divisions
      real*4, parameter  :: SKIP    = -1.E30
      integer, parameter :: varnums = 9 ! number of GSWP2 variables

      real*4, dimension(long_max,lat_maxfull,varnums)::gswp1X1_all_avg
      integer,dimension(long_max,lat_maxfull,varnums)::num_all
      integer :: i_avg,j_avg,k_avg

      do i_avg = 1,long_max
         do j_avg = 1,lat_maxfull
            do k_avg = 1,varnums 

               if (num_all(i_avg,j_avg,k_avg) > 0) then
                  gswp1X1_all_avg(i_avg,j_avg,k_avg) = 
     &                gswp1X1_all_avg(i_avg,j_avg,k_avg)
     &                / num_all(i_avg,j_avg,k_avg)
               else
                  gswp1X1_all_avg(i_avg,j_avg,k_avg) = SKIP
               endif

            enddo
         enddo
      enddo

!     gswp1X1_all_avg(:,:,:) = gswp1X1_all_avg(:,:,:)/8 

      end subroutine avg_sum

!======================================================================!

      subroutine open_outfiles(n_monthyr,iu_vector)

      integer,intent(in) :: n_monthyr ! simulation month number

      integer, parameter :: num_months = 120
      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer, parameter :: num_outfiles = num_months*varnums
      integer,dimension(varnums+3) :: iu_vector

      integer :: i
      character*80 :: infile
      character*80 :: monthyearBIN(num_months)
      character*80 :: lsm_vars(varnums+3),basefold_out
      character*80 :: lsm_fold(varnums+3)

      DATA lsm_fold/ 'onedeg_srheat/'
     &        ,'onedeg_trheat/', 'onedeg_ts/', 'onedeg_qs/'
     &        ,'onedeg_ps/','onedeg_ws/'
     &        ,'onedeg_rainf/', 'onedeg_rainf_c/', 'onedeg_snowf/'
     &        ,'rain_3hr/', 'rainc_3hr/', 'snow_3hr/'/


      DATA lsm_vars/ 'onedeg_srheat'
     &        ,'onedeg_trheat', 'onedeg_ts', 'onedeg_qs'
     &        ,'onedeg_ps','onedeg_ws'
     &        ,'onedeg_rainf', 'onedeg_rainf_c', 'onedeg_snowf'
     &        ,'rain_3hr', 'rainc_3hr', 'snow_3hr'/


      DATA monthyearBIN/   
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

      basefold_out = '/archive/g03/mjpuma/gswp1x1_avg/'

      do i = 1,(varnums+3)
         iu_vector(i) = 30 + (n_monthyr-1)*(varnums+3) + i
      enddo

!     open output file
      do i = 1,(varnums+3)
          infile=basefold_out(1:len_trim(basefold_out))//
     &      lsm_fold(i)(1:len_trim(lsm_fold(i)))//
     &      lsm_vars(i)(1:len_trim(lsm_vars(i)))//
     &      monthyearBIN(n_monthyr)(1:len_trim(monthyearBIN(n_monthyr)))
          !print *, infile
          open(iu_vector(i),FILE=infile, STATUS='UNKNOWN',
     &              FORM='UNFORMATTED')
      enddo


      end subroutine open_outfiles

!======================================================================!


      subroutine close_outfiles(iu_vector)

      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer,dimension(varnums+3) :: iu_vector
      integer i
      print *, 'in close subroutine'
      do i = 1,(varnums+3)
         close(iu_vector(i))
      enddo
      
      end subroutine close_outfiles

!======================================================================!

      end module drv_main_gswp


!======================================================================
! Purpose: Read GSWP2 data                                            =
!======================================================================

      program gswp_time
      use drv_main_gswp

      implicit none
      include 'netcdf.inc'
      integer :: i,j,nt,n_d,n_h,num_days

      integer, parameter :: num_months = 120
      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer, parameter :: long_max=360,lat_max=180 
      integer,dimension(varnums+3) :: iu_vector
      real*4, dimension(long_max,lat_max,varnums) ::gswp1X1_all_avg
      integer, dimension(long_max,lat_max,varnums)::num_all


      INTEGER, DIMENSION(num_months) :: monlen= ! # of 3-hourly timesteps
     &        (/ 
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

         num_days = monlen(i)/8

         loop_days: do n_d = 1,num_days
            ! set the average matrix to zero each day
            gswp1X1_all_avg(:,:,:) = 0.d0
            ! set number of non-SKIP values to zero each day
            num_all(:,:,:) = 0

            loop_hours: do n_h = 1,8
               nt = (n_d-1)*8 + n_h
               call gswp_timeavg(i,nt,monlen(i),iu_vector,
     &                           gswp1X1_all_avg,num_all)
               print *, 'month =',  i, 'timestep = ' , nt
               !print *, 'n_d   = ',n_d,'n_h      = ', n_h
            enddo loop_hours

            ! Divide by # of non-SKIP values to get the daily averages
            call avg_sum(gswp1X1_all_avg,num_all)
            ! Write the daily averages to file
            call write_outday(iu_vector,gswp1X1_all_avg)

         enddo loop_days

         call close_outfiles(iu_vector)

      enddo loop_months


      end program gswp_time
