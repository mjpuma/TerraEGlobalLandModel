      module gswp_open_routines
      implicit none


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

!======================================================================!

      subroutine get_single_cell(n_monthyr,nt,monlen,iu_num)
      include 'netcdf.inc'
!
      integer, parameter :: long_max=360,lat_max=150 !GSWP2 dimensions
      integer, parameter :: lat_maxfull=180 ! # of 1 deg latitude divisions
      integer, parameter :: num_months = 120!!162
      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer,parameter :: i_site=217,j_site=91 !Mpala

      integer,intent(in) :: n_monthyr ! simulation month number
      integer,intent(in) :: nt ! current month timestep (3-hourly) number
      integer,intent(in) :: monlen ! # of timesteps in current month
      integer,intent(in) :: iu_num

      real*4, parameter :: SKIP = -1.E30
      integer, parameter :: num_landcells = 15238

      character*80 :: varname, infile
      character*80 :: basefolder, foldnames(varnums),vnames(varnums)
      character*80 :: v2names(varnums), monthyear(num_months)

      real*4, dimension(long_max,lat_maxfull) :: giss2d
      real*4, dimension(long_max,lat_max) :: arr2d,lat,long
      real*4, dimension(long_max,lat_maxfull) :: gswp2d
      real*4, dimension(long_max,lat_maxfull,varnums) :: giss2d_all

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
      DATA monthyear/  !        '198207.nc','198208.nc',
!     & '198209.nc','198210.nc','198211.nc','198212.nc',
!     & '198301.nc','198302.nc','198303.nc','198304.nc',
!     & '198305.nc','198306.nc','198307.nc','198308.nc',
!     & '198309.nc','198310.nc','198311.nc','198312.nc',
!     & '198401.nc','198402.nc','198403.nc','198404.nc',
!     & '198405.nc','198406.nc','198407.nc','198408.nc',
!     & '198409.nc','198410.nc','198411.nc','198412.nc',
!     & '198501.nc','198502.nc','198503.nc','198504.nc',
!     & '198505.nc','198506.nc','198507.nc','198508.nc',
!     & '198509.nc','198510.nc','198511.nc','198512.nc',
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

      basefolder = '/discover/nobackup/mpuma/GSWP/'

      ! intialize variables--------
      giss2d_all(:,:,:) = 0.d0

      loop_gswp_variables: do kk = 1, varnums

         arr2d(:,:)= 0.d0!SKIP
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
               endif      
            enddo

            ! assign to full 360x180 matrices from 360x150
            do j = 1,lat_max
               gswp2d(:,j) = arr2d(:,j)
            enddo

            ! flip matrices to align with GISS input
            do j = 1,lat_maxfull
              giss2d(:,j) = gswp2d(:,lat_maxfull-j+1)
            enddo

            do i = 1,long_max
               do j = 1,lat_maxfull
                 giss2d_all(i,j,kk) = giss2d(i,j)
               enddo
            enddo

            deallocate(arrpack,STAT=status_alloc1)
         else
            print *, 'allocate failed'

         endif allocate_ok

      enddo loop_gswp_variables !kk

!     WRITE FORCINGS
      ! Incoming shortwave radiation [W/m2]
      ! Incoming longwave radiation [W/m2]
      ! Surface air temperature [K]
      ! Surface air moisture [kg/kg]
      ! Surface pressure [Pa]
      ! Wind speed [m/s]
      ! Rainfall rate (total) from [kg/m2/s]
      ! Convective rainfall rate [kg/m2/s]
      ! Snowfall rate [kg/m2/s]
       write(iu_num,'(9(1pe16.8))') 
     & giss2d_all(i_site,j_site,1),
     & giss2d_all(i_site,j_site,2),
     & giss2d_all(i_site,j_site,7),
     & giss2d_all(i_site,j_site,5),
     & giss2d_all(i_site,j_site,8),
     & giss2d_all(i_site,j_site,9),
     & giss2d_all(i_site,j_site,3),
     & giss2d_all(i_site,j_site,6),
     & giss2d_all(i_site,j_site,4)

      end subroutine get_single_cell

!-------------------------------------------------------------------

      subroutine open_outfiles(iu_num)
      integer,intent(out) :: iu_num
      character*80 :: infile
      character*80 :: site_nameyr,basefold_out
      DATA site_nameyr/ 'Mpala1986_1995.txt'/
      basefold_out = '/home/mpuma/gswp/'
      iu_num = 35

!     Site meteorological file
      infile=basefold_out(1:len_trim(basefold_out))//
     &      site_nameyr(1:len_trim(site_nameyr))
      print *, infile
      open(iu_num,FILE=infile,STATUS='UNKNOWN',FORM='FORMATTED')

      end subroutine open_outfiles

!-------------------------------------------------------------------

      subroutine close_outfiles(iu_num)
      integer,intent(in) :: iu_num
      print *, 'in close subroutine'
      close(iu_num)
      
      end subroutine close_outfiles

!-------------------------------------------------------------------

      end module gswp_open_routines


!======================================================================
! Purpose: Read GSWP2 data and save data for a single cell            =
!======================================================================

      program gswp_time
      use gswp_open_routines

      implicit none
      include 'netcdf.inc'
      integer :: i,j,nt
      integer :: iu_num
      integer, parameter :: long_max=360,lat_max=150 !GSWP2 dimensions
      integer, parameter :: lat_maxfull=180 ! # of 1 deg latitude divisions
      integer, parameter :: num_months = 120!162
      integer, parameter :: varnums = 9 ! number of GSWP2 variables

      INTEGER, DIMENSION(num_months) :: monlen= ! # of 3-hourly timesteps
     &        (/ !                        248,248,240,248,240,248,   !1982
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1983
!     &           248,232,248,240,248,240,248,248,240,248,240,248,   !1984
!     &           248,224,248,240,248,240,248,248,240,248,240,248,   !1985
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

      call open_outfiles(iu_num)

      loop_months: do i = 1,num_months

         loop_hours:  do nt = 1,monlen(i)

            call get_single_cell(i, nt, monlen(i), iu_num)
            print *, 'month =',  i, 'timestep = ' , nt

         enddo loop_hours

      enddo loop_months

      call close_outfiles(iu_num)

      end program gswp_time
