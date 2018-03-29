      program annual_ij
      implicit none

      character*80 TITLE,title_dum
      integer,parameter :: num_var_pr = 9 ! Number of ModelE variables to print
      real*4 precip(72,46),precip_con(72,46),topo(72,46),
     &       stress(72,46),soil_h2o_1(72,46),soil_h2o_2(72,46),
     &       ev_pen_pot(72,46),ev_bare_soil(72,46),ev_dry_canop(72,46),
     &       ev_soil(72,46),ev_wet_canop(72,46),frac_veg_soil(72,46),
     &       evap(72,46), sa_temp(72,46),grow_seas(72,46),
     &       incid_sol_rad(72,46),mat_dum(72,46),diff_max

      real*4 frac_bbare(72,46),frac_tund(72,46),frac_gras(72,46),
     &       frac_shgr(72,46),frac_trgr(72,46),frac_decfor(72,46),
     &       frac_evfor(72,46),frac_rainfor(72,46), frac_cultiv(72,46),
     &       frac_dbare(72,46),diff_test(72,46),
     &       diff_max
	 
      real*4,dimension(num_var_pr) :: sum_vals,sum_vals2,av_vals,
     &                               var_vals

      integer i,j,k,num_cells
      
      open(9,file='ANN1951-1955.ijE8piM20',STATUS = 'unknown',
     &       form='unformatted')
      open(10,file='V72X46.1.cor.ext',STATUS = 'unknown',
     &       form='unformatted')

 !     open(10,file='gras_4x5.txt',form='formatted')
 !     open(11,file='tund_4x5.txt',form='formatted')
!************************************************************************
! READ IN MODELE OUTPUT
      read(9) TITLE,topo ! Topography #1

      do k = 1,5 ! Skip 5
         read(9) title_dum,mat_dum
      enddo

      read(9) TITLE,precip ! Precipitation #7
      print *, TITLE

      ! Skip 0

      read(9) TITLE,evap ! Evaporation #8
      print *, TITLE

      do k = 1,14 ! Skip 14
         read(9) title_dum,mat_dum 
      enddo	  

      read(9) TITLE,incid_sol_rad ! incident solar radiation # 23
      print *, TITLE

      do k = 1,59 ! Skip 59
         read(9) title_dum,mat_dum 
      enddo	  

      read(9) TITLE,sa_temp ! Surface air temperature # 83
      print *, TITLE
     
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
      print *, TITLE 
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
   !   print *, TITLE	  
      do i = 1,72
         do j = 1,46
            if(frac_tund(i,j)<10e-6)then
               frac_tund(i,j) = 0.d0
            endif
         enddo
      enddo

      read(10) TITLE,frac_gras
!  print *, TITLE
	  do i = 1,72
         do j = 1,46
            if(frac_gras(i,j)<10e-6)then
               frac_gras(i,j) = 0.d0
            endif
         enddo
      enddo
	  
      read(10) TITLE,frac_shgr
!  print *, TITLE
      do i = 1,72
         do j = 1,46
            if(frac_shgr(i,j)<10e-6)then
               frac_shgr(i,j) = 0.d0
            endif
         enddo
      enddo	  
      
      read(10) TITLE,frac_trgr
   !  print *, TITLE
      do i = 1,72
         do j = 1,46
            if(frac_trgr(i,j)<10e-6)then
               frac_trgr(i,j) = 0.d0
            endif
         enddo
      enddo	  
	  
      read(10) TITLE,frac_decfor
!      print *, TITLE
      do i = 1,72
         do j = 1,46
            if(frac_decfor(i,j)<10e-6)then
               frac_decfor(i,j) = 0.d0
            endif
         enddo
      enddo	  
	  
      read(10) TITLE,frac_evfor
!      print *, TITLE
      do i = 1,72
         do j = 1,46
            if(frac_evfor(i,j)<10e-6)then
               frac_evfor(i,j) = 0.d0
            endif
         enddo
      enddo	 
	  
      read(10) TITLE,frac_rainfor
!      print *, TITLE
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
	  	  
!      do j=1,46
!         write(10,100) (matt_veg(i,j), i=1,72)
!      enddo
! 100  format(72f10.5)
!      close(10)
      diff_max = maxval(topo)
!      print *, diff_max
      diff_max = maxval(ev_pen_pot)
!      print *, diff_max
      diff_max = maxval(ev_bare_soil)
!      print *, diff_max
      diff_max = maxval(ev_dry_canop)
!      print *, diff_max
      diff_max = maxval(ev_wet_canop)
!      print *, diff_max	  	  
		  	  
      print *, 'eof reached (read completed)'

!************************************************************************
      do k = 1,num_var_pr
         sum_vals(k) = 0.d0
         av_vals(k) = 0.d0
         var_vals(k) = 0.d0
      enddo
     
      num_cells = 0

      do i = 1,72
         do j = 1,46
            if(frac_bbare(i,j)>0.d0)then
      
               sum_vals(1) = sum_vals(1) + precip(i,j)
               sum_vals(2) = sum_vals(2) + precip_con(i,j)
               sum_vals(3) = sum_vals(3) + sa_temp(i,j)
               sum_vals(4) = sum_vals(4) + incid_sol_rad(i,j)
               sum_vals(5) = sum_vals(5) + grow_seas(i,j)
               sum_vals(6) = sum_vals(6) + stress(i,j)
               sum_vals(7) = sum_vals(7) + topo(i,j)
               sum_vals(8) = sum_vals(8) + soil_h2o_1(i,j)
               sum_vals(9) = sum_vals(9) + soil_h2o_2(i,j)

               sum_vals2(1) = sum_vals2(1) + precip(i,j)
               sum_vals2(2) = sum_vals2(2) + precip_con(i,j)
               sum_vals2(3) = sum_vals2(3) + sa_temp(i,j)
               sum_vals2(4) = sum_vals2(4) + incid_sol_rad(i,j)
               sum_vals2(5) = sum_vals2(5) + grow_seas(i,j)
               sum_vals2(6) = sum_vals2(6) + stress(i,j)
               sum_vals2(7) = sum_vals2(7) + topo(i,j)
               sum_vals2(8) = sum_vals2(8) + soil_h2o_1(i,j)
               sum_vals2(9) = sum_vals2(9) + soil_h2o_2(i,j)
           
               num_cells = num_cells + 1
            endif
         enddo
      enddo

      do k = 1,num_var_pr
         av_vals(k) = sum_vals(k)/real(num_cells)
         var_vals(k) = (real(num_cells)*sum_vals2(k)-sum_vals(k)**2)/
     &                  (real(num_cells)*real(num_cells-1))
      enddo


      stop

      end program 

