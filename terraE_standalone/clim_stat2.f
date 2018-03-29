      program clim_stat2

      USE GEOM !,ONLY:DXYP
      implicit none

      character*80 TITLE,title_dum
      integer,parameter :: num_var_pr = 9 ! # of ModelE vars to print
      integer,parameter :: num_veg_typ = 10 ! # of veg types
      real*4,dimension(72,46) :: precip,precip_con,topo,
     &       stress,soil_h2o_1,soil_h2o_2,
     &       ev_pen_pot,ev_bare_soil,ev_dry_canop,
     &       ev_soil,ev_wet_canop,frac_veg_soil,
     &       evap, sa_temp,grow_seas,area_cells,
     &       incid_sol_rad,mat_dum,land_temp

      real*4,dimension(72,46) :: frac_bbare,frac_tund,frac_gras,
     &       frac_shgr,frac_trgr,frac_decfor,frac_evfor,
     &       frac_rainfor, frac_cultiv,frac_dbare,diff_test

      real*4,dimension(1000,num_veg_typ):: precip_dat,cprecip_dat,
     &      SAT_dat,ISOLR_dat,frac_dat    
	 
      real*4,dimension(num_var_pr) :: sum_vals,sum_vals2,min_vals,
     &                                max_vals
      real*4,dimension(num_var_pr,num_veg_typ) :: av_vals,var_vals,
     &                                min_vals_all,max_vals_all 
      real*4,dimension(1000)::vec_hist,frac_hist
      real*4 sum_frac

      integer i,j,k,num_cells,kk,num_temp
      
      open(9,file='ANN1951-1955.ijE8piM20',STATUS = 'unknown',
     &       form='unformatted')
      open(10,file='V72X46.1.cor2_no_crops',STATUS = 'unknown',
     &       form='unformatted')
     
      open(11,file='av_ann.txt',form='formatted')
      open(12,file='var_ann.txt',form='formatted')
      open(41,file='min_ann.txt',form='formatted')
      open(42,file='max_ann.txt',form='formatted')
      
      open(43,file='prec_dat.txt',form='formatted')
      open(44,file='cprec_dat.txt',form='formatted')      
      open(45,file='SAT_dat.txt',form='formatted')
      open(46,file='ISOLR_dat.txt',form='formatted')
      open(47,file='frac_dat.txt',form='formatted')
      open(48,file='grid_areas.txt',form='formatted')
      
      open(13,file='tundra.out',form='unformatted')
      open(14,file='grass.out',form='unformatted')
      open(15,file='shrub_gras.out',form='unformatted')
      open(16,file='tree_gras.out',form='unformatted')
      open(17,file='decfor.out',form='unformatted')
      open(18,file='evgreen.out',form='unformatted')      
      open(19,file='rainfor.out',form='unformatted')
      open(30,file='precip.out',form='unformatted')
      open(31,file='incid_solrad.out',form='unformatted')
      open(32,file='SAT.out',form='unformatted')

!***********************************************************************
      do j=1,num_veg_typ
         do i=1,1000
	    precip_dat(i,j) = -1.d0
            cprecip_dat(i,j) = -1.d0
	    SAT_dat(i,j) = -1.d0
	    ISOLR_dat(i,j) = -1.d0
	    frac_dat(i,j) = -1.d0
	 enddo
      enddo

!***********************************************************************
! READ IN MODELE OUTPUT
      read(9) TITLE,topo ! Topography #1

      do k = 1,5 ! Skip 5
         read(9) title_dum,mat_dum
      enddo

      read(9) TITLE,precip ! Precipitation #7
      write(30) TITLE,precip

      ! Skip 0

      read(9) TITLE,evap ! Evaporation #8
      print *, TITLE

      do k = 1,14 ! Skip 14
         read(9) title_dum,mat_dum 
      enddo	  

      read(9) TITLE,incid_sol_rad ! incident solar radiation # 23
      write(31) TITLE,incid_sol_rad

      do k = 1,59 ! Skip 59
         read(9) title_dum,mat_dum 
      enddo

      read(9) TITLE,sa_temp ! Surface air temperature # 83
      write(32) TITLE,sa_temp
     
      ! Skip 0
	 
      read(9) TITLE,precip_con ! Convective precipitation # 84
      print *, TITLE
	  
      do k = 1,51 ! Skip 51
         read(9) title_dum,mat_dum
      enddo	
	  
      read(9) TITLE,stress ! Plant water stress #136
      print *, TITLE
		  
      do k = 1,10 ! Skip 10
         read(9) title_dum,mat_dum
      enddo	
	  
      read(9) TITLE,ev_soil ! Soil evaporation #147
      print *, TITLE
	  
      do k = 1,37 ! Skip 37
         read(9) title_dum,mat_dum
      enddo	
	  
      read(9) TITLE,soil_h2o_1 ! Vegetated soil water layer 1 #185
      print *, TITLE	  
      
      ! Skip 0

      read(9) TITLE,soil_h2o_2 ! Vegetated soil water layer 2 #186
      print *, TITLE
	  
      do k = 1,4 ! Skip 4
         read(9) title_dum,mat_dum
      enddo	
	  
      read(9) TITLE,ev_pen_pot ! Penman potential evaporation #191
      print *, TITLE
	  
      do k = 1,3 ! Skip 3
         read(9) title_dum,mat_dum
      enddo	
	  
      read(9) TITLE,ev_bare_soil ! Bare soil evaporation #195
      print *, TITLE
	  
      ! Skip 0
	  
      read(9) TITLE,ev_dry_canop ! Dry canopy evaporation #196
      print *, TITLE
	  
      ! Skip 0	
	  
      read(9) TITLE,ev_wet_canop ! Wet canopy evaporation #197
      print *, TITLE
	  
      do k = 1,11 ! Skip 11
         read(9) title_dum,mat_dum
      enddo	
	  
      read(9) TITLE,frac_veg_soil ! Fraction of vegetated soil #209
      print *, TITLE 	  	  

      do k = 1,25 ! Skip 25
         read(9) title_dum,mat_dum
      enddo	
	  
      read(9) TITLE,grow_seas ! Growing season [days] #235
!************************************************************************
! CREATE MATRIX OF GRID CELL AREAS in [km^2]
      call GEOM_B
      do j = 1,46
         do i = 1,72
            area_cells(i,j)=DXYP(j)/1000.d0/1000.d0
         enddo
         print *, DXYP(j)
      enddo 

!************************************************************************
! READ LANDCOVER INPUT ASSIGNED IN MODELE
      read(10) TITLE,frac_bbare
   !   print *, TITLE
      do i = 1,72
         do j = 1,46
            if(frac_bbare(i,j)<10e-6)then
               frac_bbare(i,j) = 0.d0
            endif
          end do
      end do

      read(10) TITLE,frac_tund
      write(13) TITLE,frac_tund
      do i = 1,72
         do j = 1,46
            if(frac_tund(i,j)<10e-6)then
               frac_tund(i,j) = 0.d0
            endif
         enddo
      enddo

      read(10) TITLE,frac_gras
      write(14) TITLE,frac_gras
	  do i = 1,72
         do j = 1,46
            if(frac_gras(i,j)<10e-6)then
               frac_gras(i,j) = 0.d0
            endif
         enddo
      enddo
	  
      read(10) TITLE,frac_shgr
      write(15) TITLE,frac_shgr
      do i = 1,72
         do j = 1,46
            if(frac_shgr(i,j)<10e-6)then
               frac_shgr(i,j) = 0.d0
            endif
         enddo
      enddo	  
      
      read(10) TITLE,frac_trgr
      write(16) TITLE,frac_trgr
      do i = 1,72
         do j = 1,46
            if(frac_trgr(i,j)<10e-6)then
               frac_trgr(i,j) = 0.d0
            endif
         enddo
      enddo	  
	  
      read(10) TITLE,frac_decfor
      write(17) TITLE,frac_decfor
      do i = 1,72
         do j = 1,46
            if(frac_decfor(i,j)<10e-6)then
               frac_decfor(i,j) = 0.d0
            endif
         enddo
      enddo	  
	  
      read(10) TITLE,frac_evfor
      write(18) TITLE,frac_evfor
      do i = 1,72
         do j = 1,46
            if(frac_evfor(i,j)<10e-6)then
               frac_evfor(i,j) = 0.d0
            endif
         enddo
      enddo	 
	  
      read(10) TITLE,frac_rainfor
      write(19) TITLE,frac_rainfor
      do i = 1,72
         do j = 1,46
            if(frac_rainfor(i,j)<10e-6)then
               frac_rainfor(i,j) = 0.d0
            endif
         enddo
      enddo	 	  

      read(10) TITLE,frac_cultiv
!      print *, TITLE
      do i = 1,72
         do j = 1,46
            if(frac_cultiv(i,j)<10e-6)then
               frac_cultiv(i,j) = 0.d0
            endif
         enddo
      enddo	 

      read(10) TITLE,frac_dbare
!      print *, TITLE
      do i = 1,72
         do j = 1,46
            if(frac_dbare(i,j)<10e-6)then
               frac_dbare(i,j) = 0.d0
            endif
         enddo
      enddo	 
	  	  		  	  
      print *, 'eof reached (read completed)'

!***********************************************************************
      do kk = 1,num_veg_typ
         do k = 1,num_var_pr
	    av_vals(k,kk) = 0.d0
	    var_vals(k,kk) = 0.d0
	 enddo
      enddo
      ! Vegetation type loop
      do kk = 1,num_veg_typ

	 if(kk.eq.1)then
            land_temp = frac_bbare
         elseif(kk.eq.2)then
            land_temp = frac_tund
         elseif(kk.eq.3)then
            land_temp = frac_gras    
         elseif(kk.eq.4)then
            land_temp = frac_shgr     
         elseif(kk.eq.5)then
            land_temp = frac_trgr
         elseif(kk.eq.6)then
            land_temp = frac_decfor
         elseif(kk.eq.7)then
            land_temp = frac_evfor
         elseif(kk.eq.8)then
            land_temp = frac_rainfor
         elseif(kk.eq.9)then
            land_temp = frac_cultiv    
         elseif(kk.eq.10)then
            land_temp = frac_dbare
         endif
         
	 do k = 1,num_var_pr
            sum_vals(k) = 0.d0
            sum_vals2(k) = 0.d0
            min_vals(k) = 100000.d0
            max_vals(k) = 0.d0
         enddo	 
         
	 num_cells = 0

         do i = 1,72
            do j = 1,46
               if(land_temp(i,j)>0.d0)then

                  num_cells = num_cells + 1
	          
	          precip_dat(num_cells,kk) = precip(i,j)*365.d0
                  cprecip_dat(num_cells,kk) = precip_con(i,j)*365.d0
	          SAT_dat(num_cells,kk) = sa_temp(i,j)
	          ISOLR_dat(num_cells,kk) = incid_sol_rad(i,j)
		  frac_dat(num_cells,kk) = land_temp(i,j)

                  sum_vals(1) = sum_vals(1) + precip(i,j)
                  sum_vals(2) = sum_vals(2) + precip_con(i,j)
                  sum_vals(3) = sum_vals(3) + sa_temp(i,j)
                  sum_vals(4) = sum_vals(4) + incid_sol_rad(i,j)
                  sum_vals(5) = sum_vals(5) + grow_seas(i,j)
                  sum_vals(6) = sum_vals(6) + stress(i,j)
                  sum_vals(7) = sum_vals(7) + topo(i,j)
                  sum_vals(8) = sum_vals(8) + soil_h2o_1(i,j)
                  sum_vals(9) = sum_vals(9) + soil_h2o_2(i,j)
 
                  sum_vals2(1) = sum_vals2(1) + precip(i,j)**2.d0
                  sum_vals2(2) = sum_vals2(2) + precip_con(i,j)**2.d0
                  sum_vals2(3) = sum_vals2(3) + sa_temp(i,j)**2.d0
                  sum_vals2(4) = sum_vals2(4) + incid_sol_rad(i,j)**2.d0
                  sum_vals2(5) = sum_vals2(5) + grow_seas(i,j)**2.d0
                  sum_vals2(6) = sum_vals2(6) + stress(i,j)**2.d0
                  sum_vals2(7) = sum_vals2(7) + topo(i,j)**2.d0
                  sum_vals2(8) = sum_vals2(8) + soil_h2o_1(i,j)**2.d0
                  sum_vals2(9) = sum_vals2(9) + soil_h2o_2(i,j)**2.d0

                  ! calculate minimums
                  if (precip(i,j)<min_vals(1)) then
                      min_vals(1) = precip(i,j)
                  endif
                  if (precip_con(i,j)<min_vals(2)) then
                      min_vals(2) = precip_con(i,j)
                  endif
                  if (sa_temp(i,j)<min_vals(3)) then
                      min_vals(3) = sa_temp(i,j)
                  endif
                  if (incid_sol_rad(i,j)<min_vals(4)) then
                      min_vals(4) = incid_sol_rad(i,j)
                  endif
                  if (grow_seas(i,j)<min_vals(5)) then
                      min_vals(5) = grow_seas(i,j)
                  endif
                  if (stress(i,j)<min_vals(6)) then
                      min_vals(6) = stress(i,j)
                  endif
                  if (topo(i,j)<min_vals(7)) then
                      min_vals(7) = topo(i,j)
                  endif
                  if (soil_h2o_1(i,j)<min_vals(8)) then
                      min_vals(8) = soil_h2o_1(i,j)
                  endif
                  if (soil_h2o_2(i,j)<min_vals(9)) then
                      min_vals(9) = soil_h2o_2(i,j)
                  endif
                  ! calculate maximums
                  if (precip(i,j)>max_vals(1)) then
                      max_vals(1) = precip(i,j)
                  endif
                  if (precip_con(i,j)>max_vals(2)) then
                      max_vals(2) = precip_con(i,j)
                  endif
                  if (sa_temp(i,j)>max_vals(3)) then
                      max_vals(3) = sa_temp(i,j)
                  endif
                  if (incid_sol_rad(i,j)>max_vals(4)) then
                      max_vals(4) = incid_sol_rad(i,j)
                  endif
                  if (grow_seas(i,j)>max_vals(5)) then
                      max_vals(5) = grow_seas(i,j)
                  endif
                  if (stress(i,j)>max_vals(6)) then
                      max_vals(6) = stress(i,j)
                  endif
                  if (topo(i,j)>max_vals(7)) then
                      max_vals(7) = topo(i,j)
                  endif
                  if (soil_h2o_1(i,j)>max_vals(8)) then
                      max_vals(8) = soil_h2o_1(i,j)
                  endif
                  if (soil_h2o_2(i,j)>max_vals(9)) then
                      max_vals(9) = soil_h2o_2(i,j)
                  endif
		  
               endif
            enddo
         enddo

         do k = 1,num_var_pr
            av_vals(k,kk) = sum_vals(k)/real(num_cells)
            var_vals(k,kk) = (real(num_cells)*sum_vals2(k)-
     &              sum_vals(k)**2.d0)/(real(num_cells)*
     &              real(num_cells-1))
            min_vals_all(k,kk) = min_vals(k)
            max_vals_all(k,kk) = max_vals(k)
         enddo
      end do !veg types

      ! Write statistical output
      do i=1,num_var_pr
         write(11,100) (av_vals(i,j),j=1,num_veg_typ)
         write(12,100) (var_vals(i,j),j=1,num_veg_typ)
         write(41,100) (min_vals_all(i,j),j=1,num_veg_typ)
         write(42,100) (max_vals_all(i,j),j=1,num_veg_typ)
      enddo

      do i=1,1000
         write(43,100) (precip_dat(i,j),j=1,num_veg_typ)
         write(44,100) (cprecip_dat(i,j),j=1,num_veg_typ)
         write(45,100) (SAT_dat(i,j),j=1,num_veg_typ)
         write(46,100) (ISOLR_dat(i,j),j=1,num_veg_typ)
	 write(47,100) (frac_dat(i,j),j=1,num_veg_typ)
      enddo      
      
 100  format(10f20.10)

      do i=1,72
         write(48,120) (area_cells(i,j),j=1,46)
      enddo
 120  format(46f25.5)
 
      ! compute histogram
      do j=1,num_veg_typ
         do i = 1,1000
            vec_hist(i) = 0.d0
         enddo
         k=0
	 sum_frac = 0.d0
         do i = 1,1000
            if (precip_dat(i,j)>-1.d0)then
		k=k+1
	        vec_hist(k) = precip_dat(i,j)
		frac_hist(k) = frac_dat(i,j)
		sum_frac = sum_frac + frac_dat(i,j)
            endif
         enddo
	 num_temp = k
         if (num_temp>0)then
            call histogram(vec_hist,frac_hist,num_temp,sum_frac)
         endif
      enddo

      stop

      end program 
