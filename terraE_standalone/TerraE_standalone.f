!@sum TerraE Global Land Model, Standalone Version
!@sum Routine to run the standalone version of the TerraE Global Land Model
!@sum Options for meteorological forcings curently include:
!@sum     1) GSWP2 reananalysis data
!%sum     2) GISS GCM output
!@auth I. Aleinov, M.J. Puma

#define USE_GSWP_FORCINGS
#define PRINT_DIAGS

      module TerraE
      use ent_mod
      use ghy_h, only : ngm, imt, nlsn, LS_NFRAC
      use sle001, only : advnc
      implicit none

      integer, parameter :: im=144, jm=90
      integer, parameter :: j_0h=1, j_1h=jm
      type t_terraE_state
        real*8, pointer :: w_ij(:,:,:,:)      ! soil water [m]
        real*8, pointer :: ht_ij(:,:,:,:)     ! soil heat [J]
        integer, pointer :: nsn_ij(:, :, :)   ! number of snow layers [-]     
        real*8, pointer :: dzsn_ij(:, :, :, :)! snow layer thicknesses [m]
        real*8, pointer :: wsn_ij(:, :, :, :) ! snow layer water-equiv.depth[m]
        real*8, pointer :: hsn_ij(:, :, :, :) ! snow layer heat contents [J]
        real*8, pointer :: fr_snow_ij(:, :, :)! fraction of land with snowcover
        real*8, pointer :: Qf_ij(:,:)! Foliage surf. vapor mixing ratio [kg/kg]
        type (entcelltype_public), pointer :: entcells(:,:)
      end type t_terraE_state

      type t_terraE_bc
        real*8, pointer ::
     &       fearth(:,:),
     &       top_index_ij(:,:),                 
     &       top_dev_ij(:,:),                   
     &       dz_ij(:,:,:),                   
     &       q_ij(:,:,:,:),              
     &       qk_ij(:,:,:,:),             
     &       sl_ij(:,:)        
      end type t_terraE_bc

      contains
!-----------------------------------------------------------------------
      subroutine terraE_init_state( s, fearth, jday, jyear )
      type(t_terraE_state) :: s
      real*8 :: fearth(:,:)
      integer jday, jyear
      integer i,j
      integer sumj(jm)

      allocate(
     &     s%w_ij(0:ngm,ls_nfrac,im,j_0h:j_1h),
     &     s%ht_ij(0:ngm,ls_nfrac,im,j_0h:j_1h),
     &     s%nsn_ij(     2,im,j_0h:j_1h),
     &     s%dzsn_ij(nlsn,2,im,j_0h:j_1h),
     &     s%wsn_ij(nlsn,2,im,j_0h:j_1h),
     &     s%hsn_ij(nlsn,2,im,j_0h:j_1h),
     &     s%fr_snow_ij(2,im,j_0h:j_1h),
     &     s%Qf_ij(im,j_0h:j_1h) )


      read(952)
     &       s%w_ij(0:ngm,1:ls_nfrac,:,:),         
     &       s%ht_ij(0:ngm,1:ls_nfrac,:,:),        
     &       s%nsn_ij    (1:2, :,:),              
     &       s%dzsn_ij   (1:nlsn, 1:2, :,:),      
     &       s%wsn_ij    (1:nlsn, 1:2, :,:),      
     &       s%hsn_ij    (1:nlsn, 1:2, :,:),      
     &       s%fr_snow_ij(1:2, :,:)
!     &      ,s%Qf_ij(:,:)

      ! ******FIX:something is wrong with snow 
      ! remove it for now
      s%nsn_ij = 0
      s%dzsn_ij = 0.d0
      s%wsn_ij = 0.d0
      s%hsn_ij = 0.d0
      s%fr_snow_ij = 0.d0
      s%Qf_ij(:,:)=0.1d0
      ! ******

      allocate( s%entcells(im,j_0h:j_1h) )
      call ent_cell_nullify( s%entcells )

      sumj(:) = 0
      do j=1,jm
        do i=1,im
          if ( fearth(i,j) > 0.d0 )
     &         call ent_cell_construct( s%entcells(i,j) )
          if ( fearth(i,j) > 0.d0 ) sumj(j) = sumj(j) + 1
        enddo
      enddo

      print *,"sumj= ", sumj

      call set_vegetation_data( s%entcells,
     &     IM, JM, 1, IM, 1, JM, jday, jyear )

      end subroutine terraE_init_state

!-----------------------------------------------------------------------

      subroutine terraE_init_bc( bc )
      type (t_terraE_bc) :: bc

        allocate(
     &     bc%fearth      (im, j_0h:j_1h),
     &     bc%top_index_ij(im, j_0h:j_1h),
     &     bc%top_dev_ij  (im, j_0h:j_1h),
     &     bc%dz_ij       (im, j_0h:j_1h, ngm),
     &     bc%q_ij        (im, j_0h:j_1h, imt, ngm),
     &     bc%qk_ij       (im, j_0h:j_1h, imt, ngm),
     &     bc%sl_ij       (im, j_0h:j_1h) )

        read(953)
     &       bc%fearth(:,:),
     &       bc%top_index_ij(:,:),                 
     &       bc%top_dev_ij(:,:),                   
     &       bc%dz_ij(:,:,1:ngm),                   
     &       bc%q_ij(:,:,1:imt,1:ngm),              
     &       bc%qk_ij(:,:,1:imt,1:ngm),             
     &       bc%sl_ij(:,:)

      end subroutine terraE_init_bc

!-----------------------------------------------------------------------

      subroutine get_forcings(
     &         force_Ca              ,
     &         force_cos_zen_angle   ,
     &         force_vis_rad         ,
     &         force_direct_vis_rad  ,
     &         force_prec_ms         ,
     &         force_eprec_w         ,
     &         force_sprec_ms        ,
     &         force_seprec_w        ,
     &         force_srheat          ,
     &         force_trheat          ,
     &         force_ts              ,
     &         force_qs              ,
     &         force_ps              ,
     &         force_rhosrf          ,
     &         force_cdh             ,
     &         force_qm1             ,
     &         force_ws              ,
     &         force_pbl_args_ws0    ,
     &         force_pbl_args_tprime ,
     &         force_pbl_args_qprime )

      real*8
     &         force_Ca (1:im,1:jm)             ,
     &         force_cos_zen_angle (1:im,1:jm)  ,
     &         force_vis_rad (1:im,1:jm)        ,
     &         force_direct_vis_rad (1:im,1:jm) ,
     &         force_prec_ms (1:im,1:jm)        ,
     &         force_eprec_w (1:im,1:jm)        ,
     &         force_sprec_ms (1:im,1:jm)       ,
     &         force_seprec_w (1:im,1:jm)       ,
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_rhosrf (1:im,1:jm)         ,
     &         force_cdh (1:im,1:jm)            ,
     &         force_qm1 (1:im,1:jm)            ,
     &         force_ws (1:im,1:jm)             ,
     &         force_pbl_args_ws0 (1:im,1:jm)   ,
     &         force_pbl_args_tprime (1:im,1:jm),
     &         force_pbl_args_qprime (1:im,1:jm)

      print *,"reading forcings"

      read(951)
     &         force_Ca (1:im,1:jm)             ,
     &         force_cos_zen_angle (1:im,1:jm)  ,
     &         force_vis_rad (1:im,1:jm)        ,
     &         force_direct_vis_rad (1:im,1:jm) ,
     &         force_prec_ms (1:im,1:jm)        ,
     &         force_eprec_w (1:im,1:jm)        ,
     &         force_sprec_ms (1:im,1:jm)       ,
     &         force_seprec_w (1:im,1:jm)       ,
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_rhosrf (1:im,1:jm)         ,
     &         force_cdh (1:im,1:jm)            ,
     &         force_qm1 (1:im,1:jm)            ,
     &         force_ws (1:im,1:jm)             ,
     &         force_pbl_args_ws0 (1:im,1:jm)   ,
     &         force_pbl_args_tprime (1:im,1:jm),
     &         force_pbl_args_qprime (1:im,1:jm)

      print *,"forcings ok"

      end subroutine get_forcings

!-----------------------------------------------------------------------

      subroutine terraE_initialize( dt_in , n_month , tbcs_ij, tcanopy_ij
#ifdef PRINT_DIAGS
     &                  , aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m, canopyH2O_m, canopyheat_m
#endif
     &                           )

      use sle001, only : hl0, dt
      real*8,dimension(im,jm), intent(out) :: tbcs_ij, tcanopy_ij

      real*8 :: dt_in
      integer :: n_month

#ifdef PRINT_DIAGS
      real*8, dimension(im,jm,12),intent(out) ::
     &                    aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                  , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                  , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                  , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                  , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                  , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                  , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                  , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m, canopyH2O_m, canopyheat_m
      real*8,dimension(12),intent(out) :: n_count_mnth
#endif

      call ent_initialize(
     &     do_soilresp=.true.
     &     ,do_phenology=.false.
     &     ,do_frost_hardiness=.false.
     &     ,do_patchdynamics=.false.
     &     )

      n_month = 0          ! intialize month number
      dt = dt_in           ! intialize timestep size [sec]

      tbcs_ij(:,:) = -1d30 ! initialize ground surface temperature[C]
      tcanopy_ij(:,:) = -1d30 !initialize canopy temperature[C]

      call hl0

#ifdef PRINT_DIAGS
!     Initialize diagnostic accumulator
      aevap_mnth(:,:,:) = 0.d0
      aevapw_mnth(:,:,:) = 0.d0
      aevapd_mnth(:,:,:) = 0.d0
      aevapb_mnth(:,:,:) = 0.d0
      aintercep_mnth(:,:,:) = 0.d0
      aruns_mnth(:,:,:) = 0.d0
      arunu_mnth(:,:,:) = 0.d0
      agpp_mnth(:,:,:) = 0.d0
      arauto_mnth(:,:,:) = 0.d0
      asoilresp_mnth(:,:,:) = 0.d0
      aevapvg_mnth(:,:,:)=0.d0
      aevapvs_mnth(:,:,:)=0.d0
      aevapbs_mnth(:,:,:)=0.d0
      asrht_mnth(:,:,:)=0.d0
      atrht_mnth(:,:,:)=0.d0
      aalbedo_mnth(:,:,:)=0.d0
      aclab_mnth(:,:,:)=0.d0
      asoilCpoolsum_mnth(:,:,:)=0.d0
      aepp_mnth(:,:,:)=0.d0
      atrg_mnth(:,:,:)=0.d0
      ashg_mnth(:,:,:)=0.d0
      alhg_mnth(:,:,:)=0.d0
      aepc_mnth(:,:,:)=0.d0
      w_b1_mnth(:,:,:)=0.d0
      w_b2_mnth(:,:,:)=0.d0
      w_b3_mnth(:,:,:)=0.d0
      w_b4_mnth(:,:,:)=0.d0
      w_b5_mnth(:,:,:)=0.d0
      w_b6_mnth(:,:,:)=0.d0
      w_v1_mnth(:,:,:)=0.d0
      w_v2_mnth(:,:,:)=0.d0
      w_v3_mnth(:,:,:)=0.d0
      w_v4_mnth(:,:,:)=0.d0
      w_v5_mnth(:,:,:)=0.d0
      w_v6_mnth(:,:,:)=0.d0
      ht_b1_mnth(:,:,:)=0.d0
      ht_b2_mnth(:,:,:)=0.d0
      ht_b3_mnth(:,:,:)=0.d0
      ht_b4_mnth(:,:,:)=0.d0
      ht_b5_mnth(:,:,:)=0.d0
      ht_b6_mnth(:,:,:)=0.d0
      ht_v1_mnth(:,:,:)=0.d0
      ht_v2_mnth(:,:,:)=0.d0
      ht_v3_mnth(:,:,:)=0.d0
      ht_v4_mnth(:,:,:)=0.d0
      ht_v5_mnth(:,:,:)=0.d0
      ht_v6_mnth(:,:,:)=0.d0
      dzsn_b1(:,:,:)=0.d0
      dzsn_b2(:,:,:)=0.d0
      dzsn_b3(:,:,:)=0.d0
      dzsn_v1(:,:,:)=0.d0
      dzsn_v2(:,:,:)=0.d0
      dzsn_v3(:,:,:)=0.d0
      wsn_b1(:,:,:)=0.d0
      wsn_b2(:,:,:)=0.d0
      wsn_b3(:,:,:)=0.d0
      wsn_v1(:,:,:)=0.d0
      wsn_v2(:,:,:)=0.d0
      wsn_v3(:,:,:)=0.d0
      hsn_b1(:,:,:)=0.d0
      hsn_b2(:,:,:)=0.d0
      hsn_b3(:,:,:)=0.d0
      hsn_v1(:,:,:)=0.d0
      hsn_v2(:,:,:)=0.d0
      hsn_v3(:,:,:)=0.d0
      nsn_b(:,:,:)=0.d0
      nsn_v (:,:,:)=0.d0
      frsnow_b(:,:,:)=0.d0
      frsnow_v(:,:,:)=0.d0
      Qf_mnth(:,:,:)=0.d0
      abetad_mnth(:,:,:)=0.d0
      aClivepool_leaf_m(:,:,:)=0.d0
      aClivepool_froot_m(:,:,:)=0.d0
      aClivepool_wood_m(:,:,:)=0.d0
      aCdeadpool_surfmet_m(:,:,:)=0.d0
      aCdeadpool_surfstr_m(:,:,:)=0.d0
      aCdeadpool_soilmet_m(:,:,:)=0.d0
      aCdeadpool_soilstr_m(:,:,:)=0.d0
      aCdeadpool_cwd_m(:,:,:)=0.d0
      aCdeadpool_surfmic_m(:,:,:)=0.d0
      aCdeadpool_soilmic_m(:,:,:)=0.d0
      aCdeadpool_slow_m(:,:,:)=0.d0
      aCdeadpool_passive_m(:,:,:)=0.d0
      alai_m(:,:,:)=0.d0
      canopyH2O_m(:,:,:)=0.d0
      canopyheat_m(:,:,:)=0.d0
      n_count_mnth(:) = 0.d0 !Initialize counter
#endif

      end subroutine terraE_initialize

!-----------------------------------------------------------------------

      subroutine terraE_run(s,bc,jday,jyear,tyr_sec,dtsec,time,j0,j1,k_mnth
#ifdef PRINT_DIAGS
     &                 , aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
#endif
#ifdef USE_GSWP_FORCINGS
     &                 , iu_vector,tbcs_ij,tcanopy_ij
#endif
     &                   )

#ifdef USE_GSWP_FORCINGS
      use drv_gswp_force, only : get_gswp_forcings
#endif
      use sle001, only : tp,tbcs,aevap,aevapw,aevapd,aevapb,aintercep
     &     ,aruns,arunu,agpp,arauto,asoilresp
     &     ,aevapvg,aevapvs,aevapbs,asrht,atrht,aalbedo
     &     ,aclab,asoilCpoolsum,aepp,atrg,ashg,alhg,aepc
     &     ,abetad,alai
     &     ,aClivepool_leaf,aClivepool_froot,aClivepool_wood
     &     ,aCdeadpool_surfmet,aCdeadpool_surfstr,aCdeadpool_soilmet
     &     ,aCdeadpool_soilstr,aCdeadpool_cwd,aCdeadpool_surfmic
     &     ,aCdeadpool_soilmic,aCdeadpool_slow,aCdeadpool_passive
     &     ,alai

      use domain_decomp, only : mype
      type(t_terraE_state) :: s
      type (t_terraE_bc) :: bc

      integer, intent(in) :: jday, jyear, tyr_sec, j0, j1
      integer, intent(in) :: k_mnth !current month number
      real*8, intent(in) :: dtsec
      real*8, intent(in) :: time
      real*8 fb,fv
      integer i,j
      integer, save :: jday_old = -32768

      real*8 ::
     &     force_Ca (1:im,1:jm)             ,
     &     force_cos_zen_angle (1:im,1:jm)  ,
     &     force_vis_rad (1:im,1:jm)        ,
     &     force_direct_vis_rad (1:im,1:jm) ,
     &     force_prec_ms (1:im,1:jm)        ,
     &     force_eprec_w (1:im,1:jm)        ,
     &     force_sprec_ms (1:im,1:jm)       ,
     &     force_seprec_w (1:im,1:jm)       ,
     &     force_srheat (1:im,1:jm)         ,
     &     force_trheat (1:im,1:jm)         ,
     &     force_ts (1:im,1:jm)             ,
     &     force_qs (1:im,1:jm)             ,
     &     force_ps (1:im,1:jm)             ,
     &     force_rhosrf (1:im,1:jm)         ,
     &     force_cdh (1:im,1:jm)            ,
     &     force_qm1 (1:im,1:jm)            ,
     &     force_ws (1:im,1:jm)             ,
     &     force_pbl_args_ws0 (1:im,1:jm)   ,
     &     force_pbl_args_tprime (1:im,1:jm),
     &     force_pbl_args_qprime (1:im,1:jm)

#ifdef PRINT_DIAGS
      real*8, dimension(im,jm,12),intent(out) ::
     &                    aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  ,w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                  ,w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                  ,w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                  ,w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                  ,ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                  ,ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                  ,ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                  ,ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                  ,dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                  ,wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                  ,hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                  ,nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                  ,abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m

      real*8, dimension(12) :: n_count_mnth
#endif

#ifdef USE_GSWP_FORCINGS
      integer, parameter :: varnums = 9      ! # of GSWP2 variables
      integer,dimension(varnums) :: iu_vector
      real*8,dimension(im,jm), intent(inout) :: tbcs_ij,tcanopy_ij
#else
      integer,dimension(varnums) :: iu_vector
      real*8,dimension(im,jm), intent(inout) :: tbcs_ij,tcanopy_ij
#endif

      !call sysusage(mype+4,1)

#ifdef USE_GSWP_FORCINGS
      call get_gswp_forcings(
     &     jday,
     &     jyear,
     &     tyr_sec,
     &     dtsec,
     &     iu_vector,
     &     tbcs_ij,
     &     tcanopy_ij,
#else
          call get_forcings(
#endif
     &         force_Ca (1:im,1:jm)             ,
     &         force_cos_zen_angle (1:im,1:jm)  ,
     &         force_vis_rad (1:im,1:jm)        ,
     &         force_direct_vis_rad (1:im,1:jm) ,
     &         force_prec_ms (1:im,1:jm)        ,
     &         force_eprec_w (1:im,1:jm)        ,
     &         force_sprec_ms (1:im,1:jm)       ,
     &         force_seprec_w (1:im,1:jm)       ,
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_rhosrf (1:im,1:jm)         ,
     &         force_cdh (1:im,1:jm)            ,
     &         force_qm1 (1:im,1:jm)            ,
     &         force_ws (1:im,1:jm)             ,
     &         force_pbl_args_ws0 (1:im,1:jm)   ,
     &         force_pbl_args_tprime (1:im,1:jm),
     &         force_pbl_args_qprime (1:im,1:jm)! xxx
     &         )
      !call sysusage(mype+4,2)
      !call sysusage(mype+12,1)

      ! update vegetation only once per day
      if ( jday .ne. jday_old) then
        call update_vegetation_data( s%entcells(:,j0:j1),
     &       im, jm, 1, im, j0, j1, jday, jyear )
        jday_old = jday
      endif

      !call sysusage(mype+12,2)
      !call sysusage(mype+8,1)

#ifdef PRINT_DIAGS
          n_count_mnth(k_mnth) = n_count_mnth(k_mnth) + 1.d0
#endif

      loop_j:   do j = j0,j1
        loop_i: do i = 1,im
          if ( bc%fearth(i,j) <= 0.d0 ) cycle loop_i
!          write(933,*) "terraE_run: pricessing i,j ", i, j

          call ent_get_exports( s%entcells(i,j),
     &         fraction_of_vegetated_soil=fv
     &         )
          fb = 1.d0 - fv

          call advnc(
!-------------- Ent specific
     &         s%entcells(i,j), force_Ca(i,j),
     &         force_cos_zen_angle(i,j), force_vis_rad(i,j),
     &         force_direct_vis_rad(i,j),
     &         s%Qf_ij(i,j),
!-------------- prognostic vars
     &         s%w_ij(0:ngm,1:LS_NFRAC,i,j),         
     &         s%ht_ij(0:ngm,1:LS_NFRAC,i,j),        
     &         s%nsn_ij    (1:2, i, j),              
     &         s%dzsn_ij   (1:nlsn, 1:2, i, j),      
     &         s%wsn_ij    (1:nlsn, 1:2, i, j),      
     &         s%hsn_ij    (1:nlsn, 1:2, i, j),      
     &         s%fr_snow_ij(1:2, i, j),          
!-------------- BC's    
     &         bc%top_index_ij(i, j),                 
     &         bc%top_dev_ij(i, j),                   
     &         bc%dz_ij(i,j,1:ngm),                   
     &         bc%q_ij(i,j,1:imt,1:ngm),              
     &         bc%qk_ij(i,j,1:imt,1:ngm),             
     &         bc%sl_ij(i,j),                         
     &         fb,                                 
     &         fv,                        
!-------------- forcings         
     &         force_prec_ms (i,j)        ,
     &         force_eprec_w (i,j)        ,
     &         force_sprec_ms (i,j)       ,
     &         force_seprec_w (i,j)       ,
     &         force_srheat (i,j)         ,
     &         force_trheat (i,j)         ,
     &         force_ts (i,j)             ,
     &         force_qs (i,j)             ,
     &         force_ps (i,j)             ,
     &         force_rhosrf (i,j)         ,
     &         force_cdh (i,j)            ,
     &         force_qm1 (i,j)            ,
     &         force_ws (i,j)             ,
     &         force_pbl_args_ws0 (i,j)   ,
     &         force_pbl_args_tprime (i,j),
     &         force_pbl_args_qprime (i,j) )

!         Assign surface temperaute for next timestep's Ch calculation  
          tbcs_ij(i,j)=tbcs
          tcanopy_ij(i,j)=tp(0,2)

!         Assign accumulators to output diagnostocs
#ifdef PRINT_DIAGS
          aevap_mnth(i,j,k_mnth)  = aevap_mnth(i,j,k_mnth) + aevap
          aevapw_mnth(i,j,k_mnth) = aevapw_mnth(i,j,k_mnth) + aevapw
          aevapd_mnth(i,j,k_mnth) = aevapd_mnth(i,j,k_mnth) + aevapd
          aevapb_mnth(i,j,k_mnth) = aevapb_mnth(i,j,k_mnth) + aevapb
          aintercep_mnth(i,j,k_mnth) = aintercep_mnth(i,j,k_mnth) + 
     &                                 aintercep
          aruns_mnth(i,j,k_mnth)  = aruns_mnth(i,j,k_mnth) + aruns
          arunu_mnth(i,j,k_mnth)  = arunu_mnth(i,j,k_mnth) + arunu
          agpp_mnth(i,j,k_mnth)   = agpp_mnth(i,j,k_mnth) + agpp
          arauto_mnth(i,j,k_mnth) = arauto_mnth(i,j,k_mnth) + arauto
          asoilresp_mnth(i,j,k_mnth)= asoilresp_mnth(i,j,k_mnth)+ 
     &                                asoilresp
          aevapvg_mnth(i,j,k_mnth)  = aevapvg_mnth(i,j,k_mnth) + aevapvg
          aevapvs_mnth(i,j,k_mnth)  = aevapvs_mnth(i,j,k_mnth) + aevapvs
          aevapbs_mnth(i,j,k_mnth)  = aevapbs_mnth(i,j,k_mnth) + aevapbs

          asrht_mnth(i,j,k_mnth) = asrht_mnth(i,j,k_mnth) + asrht
          atrht_mnth(i,j,k_mnth) = atrht_mnth(i,j,k_mnth) + atrht
          aalbedo_mnth(i,j,k_mnth) = aalbedo_mnth(i,j,k_mnth) + aalbedo

          aclab_mnth(i,j,k_mnth) = aclab_mnth(i,j,k_mnth) + aclab
          alai_m(i,j,k_mnth) = alai_m(i,j,k_mnth) + alai

          asoilCpoolsum_mnth(i,j,k_mnth) = 
     &          asoilCpoolsum_mnth(i,j,k_mnth) + asoilCpoolsum
          aepp_mnth(i,j,k_mnth) = aepp_mnth(i,j,k_mnth) + aepp
          atrg_mnth(i,j,k_mnth) = atrg_mnth(i,j,k_mnth) + atrg
          ashg_mnth(i,j,k_mnth) = ashg_mnth(i,j,k_mnth) + ashg
          alhg_mnth(i,j,k_mnth) = alhg_mnth(i,j,k_mnth) + alhg
          aepc_mnth(i,j,k_mnth) = aepc_mnth(i,j,k_mnth) + aepc
          abetad_mnth(i,j,k_mnth) = abetad_mnth(i,j,k_mnth) + abetad

          aClivepool_leaf_m(i,j,k_mnth) = 
     &        aClivepool_leaf_m(i,j,k_mnth)  + aClivepool_leaf
          aClivepool_froot_m(i,j,k_mnth) = 
     &        aClivepool_froot_m(i,j,k_mnth) + aClivepool_froot
          aClivepool_wood_m(i,j,k_mnth) =
     &        aClivepool_wood_m(i,j,k_mnth)  + aClivepool_wood 

          aCdeadpool_surfmet_m(i,j,k_mnth) =
     &        aCdeadpool_surfmet_m(i,j,k_mnth) + aCdeadpool_surfmet 
          aCdeadpool_surfstr_m(i,j,k_mnth) =
     &        aCdeadpool_surfstr_m(i,j,k_mnth) + aCdeadpool_surfstr 
          aCdeadpool_soilmet_m(i,j,k_mnth) =
     &        aCdeadpool_soilmet_m(i,j,k_mnth) + aCdeadpool_soilmet 
          aCdeadpool_soilstr_m(i,j,k_mnth) =
     &        aCdeadpool_soilstr_m(i,j,k_mnth) + aCdeadpool_soilstr 
          aCdeadpool_cwd_m(i,j,k_mnth) =
     &        aCdeadpool_cwd_m(i,j,k_mnth) + aCdeadpool_cwd 
          aCdeadpool_surfmic_m(i,j,k_mnth) =
     &        aCdeadpool_surfmic_m(i,j,k_mnth) + aCdeadpool_surfmic 
          aCdeadpool_soilmic_m(i,j,k_mnth) =
     &        aCdeadpool_soilmic_m(i,j,k_mnth) + aCdeadpool_soilmic 
          aCdeadpool_slow_m(i,j,k_mnth) =
     &        aCdeadpool_slow_m(i,j,k_mnth) + aCdeadpool_slow 
          aCdeadpool_passive_m(i,j,k_mnth) =
     &        aCdeadpool_passive_m(i,j,k_mnth) + aCdeadpool_passive 

          ! Accumulate prognostic variables
          w_b1_mnth(i,j,k_mnth)=w_b1_mnth(i,j,k_mnth)+s%w_ij(1,1,i,j)!Bare soil
          w_b2_mnth(i,j,k_mnth)=w_b2_mnth(i,j,k_mnth)+s%w_ij(2,1,i,j)!Bare soil
          w_b3_mnth(i,j,k_mnth)=w_b3_mnth(i,j,k_mnth)+s%w_ij(3,1,i,j)!Bare soil
          w_b4_mnth(i,j,k_mnth)=w_b4_mnth(i,j,k_mnth)+s%w_ij(4,1,i,j)!Bare soil
          w_b5_mnth(i,j,k_mnth)=w_b5_mnth(i,j,k_mnth)+s%w_ij(5,1,i,j)!Bare soil
          w_b6_mnth(i,j,k_mnth)=w_b6_mnth(i,j,k_mnth)+s%w_ij(6,1,i,j)!Bare soil

          w_v1_mnth(i,j,k_mnth)=w_v1_mnth(i,j,k_mnth)+s%w_ij(1,2,i,j)!Veg soil
          w_v2_mnth(i,j,k_mnth)=w_v2_mnth(i,j,k_mnth)+s%w_ij(2,2,i,j)!Veg soil
          w_v3_mnth(i,j,k_mnth)=w_v3_mnth(i,j,k_mnth)+s%w_ij(3,2,i,j)!Veg soil
          w_v4_mnth(i,j,k_mnth)=w_v4_mnth(i,j,k_mnth)+s%w_ij(4,2,i,j)!Veg soil
          w_v5_mnth(i,j,k_mnth)=w_v5_mnth(i,j,k_mnth)+s%w_ij(5,2,i,j)!Veg soil
          w_v6_mnth(i,j,k_mnth)=w_v6_mnth(i,j,k_mnth)+s%w_ij(6,2,i,j)!Veg soil

          ht_b1_mnth(i,j,k_mnth)=ht_b1_mnth(i,j,k_mnth)+s%ht_ij(1,1,i,j)!Bare
          ht_b2_mnth(i,j,k_mnth)=ht_b2_mnth(i,j,k_mnth)+s%ht_ij(2,1,i,j)!Bare
          ht_b3_mnth(i,j,k_mnth)=ht_b3_mnth(i,j,k_mnth)+s%ht_ij(3,1,i,j)!Bare
          ht_b4_mnth(i,j,k_mnth)=ht_b4_mnth(i,j,k_mnth)+s%ht_ij(4,1,i,j)!Bare
          ht_b5_mnth(i,j,k_mnth)=ht_b5_mnth(i,j,k_mnth)+s%ht_ij(5,1,i,j)!Bare
          ht_b6_mnth(i,j,k_mnth)=ht_b6_mnth(i,j,k_mnth)+s%ht_ij(6,1,i,j)!Bare

          ht_v1_mnth(i,j,k_mnth)=ht_v1_mnth(i,j,k_mnth)+s%ht_ij(1,2,i,j)!Veg
          ht_v2_mnth(i,j,k_mnth)=ht_v2_mnth(i,j,k_mnth)+s%ht_ij(2,2,i,j)!Veg
          ht_v3_mnth(i,j,k_mnth)=ht_v3_mnth(i,j,k_mnth)+s%ht_ij(3,2,i,j)!Veg
          ht_v4_mnth(i,j,k_mnth)=ht_v4_mnth(i,j,k_mnth)+s%ht_ij(4,2,i,j)!Veg
          ht_v5_mnth(i,j,k_mnth)=ht_v5_mnth(i,j,k_mnth)+s%ht_ij(5,2,i,j)!Veg
          ht_v6_mnth(i,j,k_mnth)=ht_v6_mnth(i,j,k_mnth)+s%ht_ij(6,2,i,j)!Veg

          dzsn_b1(i,j,k_mnth) = dzsn_b1(i,j,k_mnth)+s%dzsn_ij(1,1,i,j)
          dzsn_b2(i,j,k_mnth) = dzsn_b2(i,j,k_mnth)+s%dzsn_ij(2,1,i,j)
          dzsn_b3(i,j,k_mnth) = dzsn_b3(i,j,k_mnth)+s%dzsn_ij(3,1,i,j)
          dzsn_v1(i,j,k_mnth) = dzsn_v1(i,j,k_mnth)+s%dzsn_ij(1,2,i,j)
          dzsn_v2(i,j,k_mnth) = dzsn_v2(i,j,k_mnth)+s%dzsn_ij(2,2,i,j)
          dzsn_v3(i,j,k_mnth) = dzsn_v3(i,j,k_mnth)+s%dzsn_ij(3,2,i,j)

          wsn_b1(i,j,k_mnth) = wsn_b1(i,j,k_mnth)+s%wsn_ij(1,1,i,j)
          wsn_b2(i,j,k_mnth) = wsn_b2(i,j,k_mnth)+s%wsn_ij(2,1,i,j)
          wsn_b3(i,j,k_mnth) = wsn_b3(i,j,k_mnth)+s%wsn_ij(3,1,i,j)
          wsn_v1(i,j,k_mnth) = wsn_v1(i,j,k_mnth)+s%wsn_ij(1,2,i,j)
          wsn_v2(i,j,k_mnth) = wsn_v2(i,j,k_mnth)+s%wsn_ij(2,2,i,j)
          wsn_v3(i,j,k_mnth) = wsn_v3(i,j,k_mnth)+s%wsn_ij(3,2,i,j)

          hsn_b1(i,j,k_mnth) = hsn_b1(i,j,k_mnth)+s%hsn_ij(1,1,i,j)
          hsn_b2(i,j,k_mnth) = hsn_b2(i,j,k_mnth)+s%hsn_ij(2,1,i,j)
          hsn_b3(i,j,k_mnth) = hsn_b3(i,j,k_mnth)+s%hsn_ij(3,1,i,j)
          hsn_v1(i,j,k_mnth) = hsn_v1(i,j,k_mnth)+s%hsn_ij(1,2,i,j)
          hsn_v2(i,j,k_mnth) = hsn_v2(i,j,k_mnth)+s%hsn_ij(2,2,i,j)
          hsn_v3(i,j,k_mnth) = hsn_v3(i,j,k_mnth)+s%hsn_ij(3,2,i,j)

          nsn_b(i,j,k_mnth) = nsn_b(i,j,k_mnth)+ s%nsn_ij(1,i,j)
          nsn_v(i,j,k_mnth) = nsn_v(i,j,k_mnth)+ s%nsn_ij(2,i,j)

          frsnow_b(i,j,k_mnth)=frsnow_b(i,j,k_mnth)+s%fr_snow_ij(1,i,j)
          frsnow_v(i,j,k_mnth)=frsnow_v(i,j,k_mnth)+s%fr_snow_ij(2,i,j)

          Qf_mnth(i,j,k_mnth)=Qf_mnth(i,j,k_mnth)+ s%Qf_ij(i,j)
          canopyH2O_m(i,j,k_mnth)=canopyH2O_m(i,j,k_mnth) + 
     &                                   s%w_ij(0,2,i,j)
          canopyheat_m(i,j,k_mnth)=canopyheat_m(i,j,k_mnth) +
     &                                   s%ht_ij(0,2,i,j)
#endif

        enddo loop_i
      enddo loop_j

      end subroutine terraE_run

!-----------------------------------------------------------------------

      subroutine set_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, year )
!@sum read standard TerraE vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      use ent_prescribed_drv, only : init_canopy_physical,prescr_vegdata
      use ent_prescr_veg, only : prescr_calc_shc,prescr_calcconst
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      real*8, dimension(N_COVERTYPES) :: dbhdata !Diameter at breast height for woody veg.(cm)
      real*8, dimension(N_COVERTYPES) :: craddata !Crown radius (m)
      real*8, dimension(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1) :: cpooldata !Carbon pools in individuals
      integer, dimension(N_COVERTYPES) :: soildata ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpool_ini
      !-----Local---------
      integer i,j
      real*8 heat_capacity

      call prescr_calcconst()

      call prescr_vegdata(jday, year, 
     &     IM,JM,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
     &     soildata,soil_texture,Tpool_ini,.true.)
      call init_canopy_physical(I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
      call ent_cell_set(entcells, vegdata, popdata, laidata,
     &     hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata, 
     &     soildata, albedodata, soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini)

      ! just in case, do nothing, just set heat capacities
      call ent_prescribe_vegupdate(entcells)

      end subroutine set_vegetation_data

!-----------------------------------------------------------------------

      subroutine update_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, year )
!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      use ent_prescribed_drv, only:
     &     prescr_get_laidata,prescr_veg_albedodata
      !use ent_prescr_veg, only: prescr_get_laidata,prescr_veg_albedodata
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      !-----Local---------
      integer hemi(I0:I1,J0:J1)
      integer i,j
      real*8 , dimension(N_COVERTYPES):: LAI_buf

!!!! HACK : trying to update Ent exactly like in ent_prog
            !* Set hemisphere flags.
      if ( J0<=JM/2 )   hemi(:,J0:min(JM/2,J1))   = -1    ! S.
      if ( J1>=JM/2+1 ) hemi(:,max(JM/2+1,J0):J1) =  1    ! N.
 
          call ent_prescribe_vegupdate(entcells,hemi,jday,year,
     &         do_giss_phenology=.true.,
     &         do_giss_albedo=.true.,
     &         do_giss_lai=.true.,
     &         update_crops=.false. )

      return
!!!!!END HACK

      call prescr_get_laidata(jday,JM,I0,I1,J0,J1,laidata)
      call prescr_veg_albedodata(jday,JM,I0,I1,J0,J1,albedodata)

      call ent_prescribe_vegupdate(entcells
     &     ,laidata=laidata(2:2+N_PFT-1,:,:)
     &     ,albedodata=albedodata(:,2:2+N_PFT-1,:,:)
     &     )

      end subroutine update_vegetation_data

!-----------------------------------------------------------------------

      subroutine heat_to_temperature(tp, fice, ht, w, ht_cap)
      real*8, intent(out) :: tp, fice
      real*8, intent(in) :: ht, w, ht_cap

      real*8,parameter::rhow = 1000.d0
      real*8,parameter::lhm = 334590.d0 !J/kg latent heat of melting at 0 deg C
      real*8,parameter::shw_kg = 4185.d0 !J/kg C heat capacity of water at 20 C
      real*8,parameter::shi_kg = 2060.d0 !J/kg heat capacity of pure ice at 0C
      ! volumetric quantities
      real*8, parameter :: lhmv=lhm*rhow, shwv=shw_kg*rhow,
     &     shiv=shi_kg*rhow

      fice = 0.d0
      if( lhmv*w+ht < 0.d0 ) then ! all frozen
        tp = ( ht + w*lhmv )/( ht_cap + w*shiv )
        fice = 1.d0
      else if( ht > 0.d0 ) then ! all melted
        tp = ht /( ht_cap + w*shwv )
      else  ! part frozen
        tp = 0.d0
        if( w > 1d-12 ) fice = -ht /( lhmv*w )
      endif

      end subroutine heat_to_temperature

      end module terraE


!-----------------------------------------------------------------------

      subroutine terraE_standalone
      use parser
      use param
      use filemanager, only : openunit,closeunit
      use terraE
      use ent_mod
#ifdef USE_GSWP_FORCINGS
      use drv_gswp_force, only : init_forcings, get_month
#endif
      use domain_decomp, only : app_init, array_gather, mype
      implicit none
#ifdef USE_ESMF
#include "mpi_defs.h"
#include "mpif.h"
#endif

      type (t_terraE_bc) :: bc
      type(t_terraE_state) :: terraE_state
      integer jday, jyear, tyr_sec, itime_3hr
      integer :: itime,itime0=1,itime1=17520!spin83_85=52608;86_95=175296
      real*8 ::  dtsec=1800.d0     ! timestep size
      real*8 :: time
      integer :: n_month, i_mnth  ! current month number & month counter
      
#ifdef PRINT_DIAGS  /*  Monthly diagnostic accumulators */
      character*80 :: title       ! title of diagnostic to print to file
      real*8, dimension(im,jm,12) ::
     &                    aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                  , aintercep_mnth,aruns_mnth,arunu_mnth
     &                  , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                  , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                  , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                  , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                  , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                  , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                  , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                  , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                  , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                  , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                  , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                  , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                  , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
      real*8,dimension(12) :: n_count_mnth
#endif

#ifdef USE_GSWP_FORCINGS
      integer, parameter :: varnums = 9        ! # of GSWP2 variables
      integer,parameter :: num_times_83_85 = 2 ! # of times to repeat
      integer :: i_repeat83_85                 ! counter
      real*8,dimension(im,jm), save :: tbcs_ij, tcanopy_ij
      integer,dimension(varnums),save :: iu_vector
#endif

      character*120 :: ifile="ent_input"
      integer iu_IFILE
      integer j0, j1, ierr

      call app_init(jm,j0,j1)


#ifdef USE_GSWP_FORCINGS
      ! set jday, jyear etc. fro GSWP run
      jyear = 1983
      jday = 1
      tyr_sec = (jday-1)*24*3600
      i_repeat83_85 = 1
#else
      ! set jday, jyear as in default GCM setup
      jday = 335
      jyear = 1949
      jday = 1
      tyr_sec = (jday-1)*24*3600
#endif

      ! read input file with run settings
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      call parse_params(iu_IFILE)
      call closeunit(iu_IFILE)

      ! get some run parameters
      call sync_param( "itime1", itime1 )
      call sync_param( "jyear", jyear )
      call sync_param( "jday", jday )
      call sync_param( "tyr_sec", tyr_sec )


!need to initialize GHY and ENT here ...
      call terraE_initialize( dtsec,n_month,tbcs_ij,tcanopy_ij
#ifdef PRINT_DIAGS
     &                 , aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
#endif
     &                     )

      call terraE_init_bc( bc )
      call terraE_init_state( terraE_state, bc%fearth, jday, jyear )

#ifdef USE_GSWP_FORCINGS
      ! initialize GSWP forcings
      call init_forcings(jday,jyear,(jday-1)*24*3600,dtsec,iu_vector
     &                  ,tbcs_ij,tcanopy_ij)
#endif

      call sysusage(mype,0)
      call sysusage(mype,1)

      ! main time loop
      time = 0.d0
      do itime=itime0,itime1
        !print *, 'itime=', itime
        itime_3hr = tyr_sec/10800 + 1 ! 10800 is # of sec in 3 hr
        call get_month(jyear,itime_3hr,n_month) ! month # from start of year

        call terraE_run(terraE_state,bc,jday,jyear,tyr_sec,dtsec,time
     &                 ,j0,j1,n_month
#ifdef PRINT_DIAGS
     &                 , aevap_mnth,aevapw_mnth,aevapd_mnth,aevapb_mnth
     &                 , aintercep_mnth,aruns_mnth,arunu_mnth
     &                 , agpp_mnth,arauto_mnth,asoilresp_mnth
     &                 , aevapvg_mnth,aevapvs_mnth,aevapbs_mnth
     &                 , asrht_mnth,atrht_mnth,aalbedo_mnth
     &                 , aclab_mnth,asoilCpoolsum_mnth,aepp_mnth
     &                 , atrg_mnth,ashg_mnth,alhg_mnth,aepc_mnth
     &                 , n_count_mnth
     &                 , w_b1_mnth,w_b2_mnth,w_b3_mnth
     &                 , w_b4_mnth,w_b5_mnth,w_b6_mnth
     &                 , w_v1_mnth,w_v2_mnth,w_v3_mnth
     &                 , w_v4_mnth,w_v5_mnth,w_v6_mnth
     &                 , ht_b1_mnth,ht_b2_mnth,ht_b3_mnth
     &                 , ht_b4_mnth,ht_b5_mnth,ht_b6_mnth
     &                 , ht_v1_mnth,ht_v2_mnth,ht_v3_mnth
     &                 , ht_v4_mnth,ht_v5_mnth,ht_v6_mnth
     &                 , dzsn_b1,dzsn_b2,dzsn_b3,dzsn_v1,dzsn_v2,dzsn_v3
     &                 , wsn_b1 ,wsn_b2 ,wsn_b3 ,wsn_v1 ,wsn_v2 ,wsn_v3
     &                 , hsn_b1 ,hsn_b2 ,hsn_b3 ,hsn_v1 ,hsn_v2 ,hsn_v3
     &                 , nsn_b,nsn_v,frsnow_b,frsnow_v,Qf_mnth
     &                 , abetad_mnth
     &  , aClivepool_leaf_m,aClivepool_froot_m,aClivepool_wood_m
     &  , aCdeadpool_surfmet_m,aCdeadpool_surfstr_m,aCdeadpool_soilmet_m
     &  , aCdeadpool_soilstr_m,aCdeadpool_cwd_m,aCdeadpool_surfmic_m
     &  , aCdeadpool_soilmic_m,aCdeadpool_slow_m,aCdeadpool_passive_m
     &  , alai_m,canopyH2O_m,canopyheat_m
#endif
#ifdef USE_GSWP_FORCINGS
     &                 , iu_vector, tbcs_ij, tcanopy_ij
#endif
     &                   )

        tyr_sec = tyr_sec + nint(dtsec)
        jday = tyr_sec/(24*3600) + 1
        
        time = time + dtsec
        ! if end of year - reset time variables and print diags
        if (      (jday > 365 .and. mod(jyear,4).ne.0) ! non-leap
     &       .or.  jday > 366 ) then                   ! leap
          jday = 1
          tyr_sec = 0
!         Reset year back to 1983 if repeating spinup years 1983 to 1985
          if (i_repeat83_85<num_times_83_85)then
             jyear = jyear + 1
             if (jyear>1985)then
                i_repeat83_85 = i_repeat83_85 + 1
                jyear = 1983
             endif
          else ! Otherwise just increase year 
             jyear = jyear + 1
          endif

#ifdef PRINT_DIAGS
#ifdef USE_ESMF
          !print *, n_count_mnth
          do i_mnth=1,12
             title = 'TerraE: Total evaporation [kg/m2/month]'
             call array_gather( aevap_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(980) title,
     &                        real(aevap_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Canopy evaporation (evapw)[kg/m2/month]'
             call array_gather( aevapw_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(981) title,
     &                        real(aevapw_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Transpiration (evapd)[kg/m2/month]'
             call array_gather( aevapd_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(982) title,
     &                         real(aevapd_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Bare soil evap. (evapb)[kg/m2/month]'
             call array_gather( aevapb_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(983) title,
     &                         real(aevapb_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Canopy interception'
             call array_gather( aintercep_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(984) title,
     &                        real(aintercep_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Surface runoff [kg/m2/mnth]'
             call array_gather( aruns_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(985) title,
     &                         real(aruns_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Subsurface flow out of soil[kg/m2/mnth]'
             call array_gather( arunu_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(986) title,
     &                         real(arunu_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Gross primary productivity[kgC/m2/mnth]'
             call array_gather( agpp_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(987) title,
     &                         real(agpp_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Autotrophic respiration[kgC/m2/mnth]'
             call array_gather( arauto_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(988) title,
     &                         real(arauto_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil respiration[kgC/m2/mnth]'
             call array_gather( asoilresp_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(989) title,
     &                         real(asoilresp_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Evap. from vegetated soil[kg/m2/month]'
             call array_gather( aevapvg_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(990) title,
     &                         real(aevapvg_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Evap. from snow on veg[kg/m2/month]'
             call array_gather( aevapvs_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(991) title,
     &                         real(aevapvs_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Evap. from snow on bare soil[kg/m2/mn]'
             call array_gather( aevapbs_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(992) title,
     &                         real(aevapbs_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Net shortwave radiation [W/m2]'
             call array_gather( asrht_mnth(:,:,i_mnth) )
             asrht_mnth(:,:,i_mnth) = 
     &                asrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(993) title,
     &                         real(asrht_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Net longwave radiation [W/m2]'
             call array_gather( atrht_mnth(:,:,i_mnth) )
              atrht_mnth(:,:,i_mnth) = 
     &                atrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(994) title,
     &                         real(atrht_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Grid-cell mean albedo'
              aalbedo_mnth(:,:,i_mnth) = 
     &                aalbedo_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             call array_gather( aalbedo_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(995) title,
     &                          real(aalbedo_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Plant labile carbon [kgC/m2]'
             aclab_mnth(:,:,i_mnth) = 
     &                aclab_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             call array_gather( aclab_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(996) title,
     &                         real(aclab_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Sum of soil carbon pools[kgC/m2]'
             asoilCpoolsum_mnth(:,:,i_mnth) = 
     &             asoilCpoolsum_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             call array_gather( asoilCpoolsum_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(997) title,
     &                      real(asoilCpoolsum_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Penman potential evap. [kg/m2/month]'
             call array_gather( aepp_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(998) title,
     &                      real(aepp_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Thermal heat from ground [J/m2/mnth]'
             call array_gather( atrg_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(999) title,
     &                      real(atrg_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Sensible heat from ground [J/m2/mnth]'
             call array_gather( ashg_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(1000) title,
     &                      real(ashg_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Latent heat flux [W/m2]'
             call array_gather( alhg_mnth(:,:,i_mnth) )
             alhg_mnth(:,:,i_mnth) = 
     &                alhg_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1001) title,
     &                      real(alhg_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Pot. evap. from canopy [kg/m2/mnth]'
             call array_gather( aepc_mnth(:,:,i_mnth) )
             if ( mype == 0 ) write(1002) title,
     &                      real(aepc_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Bare layer 1 [m]'
             call array_gather( w_b1_mnth(:,:,i_mnth) )
             w_b1_mnth(:,:,i_mnth) = 
     &                w_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1003) title,
     &                      real(w_b1_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Bare layer 2 [m]'
             call array_gather( w_b2_mnth(:,:,i_mnth) )
             w_b2_mnth(:,:,i_mnth) = 
     &                w_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1004) title,
     &                      real(w_b2_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Bare layer 3 [m]'
             call array_gather( w_b3_mnth(:,:,i_mnth) )
             w_b3_mnth(:,:,i_mnth) = 
     &                w_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1005) title,
     &                      real(w_b3_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Bare layer 4 [m]'
             call array_gather( w_b4_mnth(:,:,i_mnth) )
             w_b4_mnth(:,:,i_mnth) = 
     &                w_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1006) title,
     &                      real(w_b4_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Bare layer 5 [m]'
             call array_gather( w_b5_mnth(:,:,i_mnth) )
             w_b5_mnth(:,:,i_mnth) = 
     &                w_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1007) title,
     &                      real(w_b5_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Bare layer 6 [m]'
             call array_gather( w_b6_mnth(:,:,i_mnth) )
             w_b6_mnth(:,:,i_mnth) = 
     &                w_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1008) title,
     &                      real(w_b6_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Veg layer 1 [m]'
             call array_gather( w_v1_mnth(:,:,i_mnth) )
             w_v1_mnth(:,:,i_mnth) = 
     &                w_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1009) title,
     &                      real(w_v1_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Veg layer 2 [m]'
             call array_gather( w_v2_mnth(:,:,i_mnth) )
             w_v2_mnth(:,:,i_mnth) = 
     &                w_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1010) title,
     &                      real(w_v2_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Veg layer 3 [m]'
             call array_gather( w_v3_mnth(:,:,i_mnth) )
             w_v3_mnth(:,:,i_mnth) = 
     &                w_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1011) title,
     &                      real(w_v3_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Veg layer 4 [m]'
             call array_gather( w_v4_mnth(:,:,i_mnth) )
             w_v4_mnth(:,:,i_mnth) = 
     &                w_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1012) title,
     &                      real(w_v4_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Veg layer 5 [m]'
             call array_gather( w_v5_mnth(:,:,i_mnth) )
             w_v5_mnth(:,:,i_mnth) = 
     &                w_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1013) title,
     &                      real(w_v5_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil water; Veg layer 6 [m]'
             call array_gather( w_v6_mnth(:,:,i_mnth) )
             w_v6_mnth(:,:,i_mnth) = 
     &                w_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1014) title,
     &                      real(w_v6_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Bare layer 1 [J/m2]'
             call array_gather( ht_b1_mnth(:,:,i_mnth) )
             ht_b1_mnth(:,:,i_mnth) = 
     &                ht_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1015) title,
     &                      real(ht_b1_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Bare layer 2 [J/m2]'
             call array_gather( ht_b2_mnth(:,:,i_mnth) )
             ht_b2_mnth(:,:,i_mnth) = 
     &                ht_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1016) title,
     &                      real(ht_b2_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Bare layer 3 [J/m2]'
             call array_gather( ht_b3_mnth(:,:,i_mnth) )
             ht_b3_mnth(:,:,i_mnth) = 
     &                ht_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1017) title,
     &                      real(ht_b3_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Bare layer 4 [J/m2]'
             call array_gather( ht_b4_mnth(:,:,i_mnth) )
             ht_b4_mnth(:,:,i_mnth) = 
     &                ht_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1018) title,
     &                      real(ht_b4_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Bare layer 5 [J/m2]'
             call array_gather( ht_b5_mnth(:,:,i_mnth) )
             ht_b5_mnth(:,:,i_mnth) = 
     &                ht_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1019) title,
     &                      real(ht_b5_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Bare layer 6 [J/m2]'
             call array_gather( ht_b6_mnth(:,:,i_mnth) )
             ht_b6_mnth(:,:,i_mnth) = 
     &                ht_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1020) title,
     &                      real(ht_b6_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Veg layer 1 [J/m2]'
             call array_gather( ht_v1_mnth(:,:,i_mnth) )
             ht_v1_mnth(:,:,i_mnth) = 
     &                ht_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1021) title,
     &                      real(ht_v1_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Veg layer 2 [J/m2]'
             call array_gather( ht_v2_mnth(:,:,i_mnth) )
             ht_v2_mnth(:,:,i_mnth) = 
     &                ht_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1022) title,
     &                      real(ht_v2_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Veg layer 3 [J/m2]'
             call array_gather( ht_v3_mnth(:,:,i_mnth) )
             ht_v3_mnth(:,:,i_mnth) = 
     &                ht_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1023) title,
     &                      real(ht_v3_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Veg layer 4 [J/m2]'
             call array_gather( ht_v4_mnth(:,:,i_mnth) )
             ht_v4_mnth(:,:,i_mnth) = 
     &                ht_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1024) title,
     &                      real(ht_v4_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Veg layer 5 [J/m2]'
             call array_gather( ht_v5_mnth(:,:,i_mnth) )
             ht_v5_mnth(:,:,i_mnth) = 
     &                ht_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1025) title,
     &                      real(ht_v5_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Soil heat; Veg layer 6 [J/m2]'
             call array_gather( ht_v6_mnth(:,:,i_mnth) )
             ht_v6_mnth(:,:,i_mnth) = 
     &                ht_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1026) title,
     &                      real(ht_v6_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow layer thickness; bare layer 1 [m]'
             call array_gather( dzsn_b1(:,:,i_mnth) )
             dzsn_b1(:,:,i_mnth) = 
     &                dzsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1027) title,
     &                      real(dzsn_b1(:,:,i_mnth),kind=4)

             title = 'TerraE:  Snow layer thickness; bare layer 2 [m]'
             call array_gather( dzsn_b2(:,:,i_mnth) )
             dzsn_b2(:,:,i_mnth) = 
     &                dzsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1028) title,
     &                      real(dzsn_b2(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow layer thickness; bare layer 3 [m]'
             call array_gather( dzsn_b3(:,:,i_mnth) )
             dzsn_b3(:,:,i_mnth) = 
     &                dzsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1029) title,
     &                      real(dzsn_b3(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow layer thickness; bare layer 1 [m]'
             call array_gather( dzsn_v1(:,:,i_mnth) )
             dzsn_v1(:,:,i_mnth) = 
     &                dzsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1030) title,
     &                      real(dzsn_v1(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow layer thickness; bare layer 2 [m]'
             call array_gather( dzsn_v2(:,:,i_mnth) )
             dzsn_v2(:,:,i_mnth) = 
     &                dzsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1031) title,
     &                      real(dzsn_v2(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow layer thickness; bare layer 3 [m]'
             call array_gather( dzsn_v3(:,:,i_mnth) )
             dzsn_v3(:,:,i_mnth) = 
     &                dzsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1032) title,
     &                      real(dzsn_v3(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow water equival; bare layer 1 [m]'
             call array_gather( wsn_b1(:,:,i_mnth) )
             wsn_b1(:,:,i_mnth) = 
     &                wsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1033) title,
     &                      real(wsn_b1(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow water equival; bare layer 2 [m]'
             call array_gather( wsn_b2(:,:,i_mnth) )
             wsn_b2(:,:,i_mnth) = 
     &                wsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1034) title,
     &                      real(wsn_b2(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow water equival; bare layer 3 [m]'
             call array_gather( wsn_b3(:,:,i_mnth) )
             wsn_b3(:,:,i_mnth) = 
     &                wsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1035) title,
     &                      real(wsn_b3(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow water equival; veg layer 1 [m]'
             call array_gather( wsn_v1(:,:,i_mnth) )
             wsn_v1(:,:,i_mnth) = 
     &                wsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1036) title,
     &                      real(wsn_v1(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow water equival; veg layer 2 [m]'
             call array_gather( wsn_v2(:,:,i_mnth) )
             wsn_v2(:,:,i_mnth) = 
     &                wsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1037) title,
     &                      real(wsn_v2(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow water equival; veg layer 3 [m]'
             call array_gather( wsn_v3(:,:,i_mnth) )
             wsn_v3(:,:,i_mnth) = 
     &                wsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1038) title,
     &                      real(wsn_v3(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow layer heat; bare layer 1 [J/m2]'
             call array_gather( hsn_b1(:,:,i_mnth) )
             hsn_b1(:,:,i_mnth) = 
     &                hsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1039) title,
     &                      real(hsn_b1(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow layer heat; bare layer 1 [J/m2]'
             call array_gather( hsn_b2(:,:,i_mnth) )
             hsn_b2(:,:,i_mnth) = 
     &                hsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1040) title,
     &                      real(hsn_b2(:,:,i_mnth),kind=4)

             title = 'TerraE:  Snow layer heat; bare layer 1 [J/m2]'
             call array_gather( hsn_b3(:,:,i_mnth) )
             hsn_b3(:,:,i_mnth) = 
     &                hsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1041) title,
     &                      real(hsn_b3(:,:,i_mnth),kind=4)

             title = 'TerraE:  Snow layer heat; veg layer 1 [J/m2]'
             call array_gather( hsn_v1(:,:,i_mnth) )
             hsn_v1(:,:,i_mnth) = 
     &                hsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1042) title,
     &                      real(hsn_v1(:,:,i_mnth),kind=4)

             title = 'TerraE:  Snow layer heat; veg layer 2 [J/m2]'
             call array_gather( hsn_v2(:,:,i_mnth) )
             hsn_v2(:,:,i_mnth) = 
     &                hsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1043) title,
     &                      real(hsn_v2(:,:,i_mnth),kind=4)

             title = 'TerraE:  Snow layer heat; veg layer 3 [J/m2]'
             call array_gather( hsn_v3(:,:,i_mnth) )
             hsn_v3(:,:,i_mnth) = 
     &                hsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1044) title,
     &                      real(hsn_v3(:,:,i_mnth),kind=4)

             title = 'TerraE: # of snow layers on bare land'
             call array_gather( nsn_b(:,:,i_mnth) )
             nsn_b(:,:,i_mnth) = 
     &                nsn_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1045) title,
     &                      real(nsn_b(:,:,i_mnth),kind=4)

             title = 'TerraE: # of snow layers on vegetated land'
             call array_gather( nsn_v(:,:,i_mnth) )
             nsn_v(:,:,i_mnth) = 
     &                nsn_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1046) title,
     &                      real(nsn_v(:,:,i_mnth),kind=4)

             title = 'TerraE: Snow-covered bare fraction [-]'
             call array_gather( frsnow_b(:,:,i_mnth) )
             frsnow_b(:,:,i_mnth) = 
     &                frsnow_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1047) title,
     &                      real(frsnow_b(:,:,i_mnth),kind=4)

             title = 'TerraE:  Snow-covered vegetated fraction [-]'
             call array_gather( frsnow_v(:,:,i_mnth) )
             frsnow_v(:,:,i_mnth) = 
     &                frsnow_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1048) title,
     &                      real(frsnow_v(:,:,i_mnth),kind=4)

             title = 'TerraE:Foliage surf. vapor mixing ratio [kg/kg]'
             call array_gather( Qf_mnth(:,:,i_mnth) )
             Qf_mnth(:,:,i_mnth) = 
     &                Qf_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1049) title,
     &                      real(Qf_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE:root zone betad [-]'
             call array_gather( abetad_mnth(:,:,i_mnth) )
             abetad_mnth(:,:,i_mnth) = 
     &                abetad_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1050) title,
     &                      real(abetad_mnth(:,:,i_mnth),kind=4)

             title = 'TerraE:Live leaf carbon pool [kg/m2]'
             call array_gather( aClivepool_leaf_m(:,:,i_mnth) )
             aClivepool_leaf_m(:,:,i_mnth) = 
     &              aClivepool_leaf_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1051) title,
     &                      real(aClivepool_leaf_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Live froot carbon pool [kg/m2]'
             call array_gather( aClivepool_froot_m(:,:,i_mnth) )
             aClivepool_froot_m(:,:,i_mnth) = 
     &              aClivepool_froot_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1052) title,
     &                    real(aClivepool_froot_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Live wood carbon pool [kg/m2]'
             call array_gather( aClivepool_wood_m(:,:,i_mnth) )
             aClivepool_wood_m(:,:,i_mnth) = 
     &              aClivepool_wood_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1053) title,
     &                    real(aClivepool_wood_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Dead surface metabolic C pool [kg/m2]'
             call array_gather( aCdeadpool_surfmet_m(:,:,i_mnth) )
             aCdeadpool_surfmet_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1054) title,
     &                    real(aCdeadpool_surfmet_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Dead surface structural C pool [kg/m2]'
             call array_gather(aCdeadpool_surfstr_m(:,:,i_mnth) )
             aCdeadpool_surfstr_m(:,:,i_mnth) = 
     &           aCdeadpool_surfstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1055) title,
     &                    real(aCdeadpool_surfstr_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Dead soil metabolic C pool [kg/m2]'
             call array_gather(aCdeadpool_soilmet_m(:,:,i_mnth) )
             aCdeadpool_soilmet_m(:,:,i_mnth) = 
     &           aCdeadpool_soilmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1056) title,
     &                    real(aCdeadpool_soilmet_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Dead soil structural C pool [kg/m2]'
             call array_gather(aCdeadpool_soilstr_m(:,:,i_mnth) )
             aCdeadpool_soilstr_m(:,:,i_mnth) = 
     &           aCdeadpool_soilstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1057) title,
     &                    real(aCdeadpool_soilstr_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Coarse woody debris C pool [kg/m2]'
             call array_gather( aCdeadpool_cwd_m(:,:,i_mnth) )
             aCdeadpool_cwd_m(:,:,i_mnth) = 
     &            aCdeadpool_cwd_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1058) title,
     &                    real(aCdeadpool_cwd_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Dead surface microbial C pool [kg/m2]'
             call array_gather(aCdeadpool_surfmic_m(:,:,i_mnth) )
             aCdeadpool_surfmic_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1059) title,
     &                    real(aCdeadpool_surfmic_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Dead soil microbial C pool [kg/m2]'
             call array_gather( aCdeadpool_soilmic_m(:,:,i_mnth) )
             aCdeadpool_soilmic_m(:,:,i_mnth) = 
     &            aCdeadpool_soilmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1060) title,
     &                    real(aCdeadpool_soilmic_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Dead slow (~10 yr) C pool [kg/m2]'
             call array_gather(aCdeadpool_slow_m(:,:,i_mnth) )
             aCdeadpool_slow_m(:,:,i_mnth) = 
     &            aCdeadpool_slow_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1061) title,
     &                    real(aCdeadpool_slow_m(:,:,i_mnth),kind=4)

             title = 'TerraE:Dead passive (~100 yr) C pool [kg/m2]'
             call array_gather(aCdeadpool_passive_m(:,:,i_mnth) )
             aCdeadpool_passive_m(:,:,i_mnth) = 
     &            aCdeadpool_passive_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1062) title,
     &                    real(aCdeadpool_passive_m(:,:,i_mnth),kind=4)

             title = 'TerraE: Leaf Area Index [-]'
             call array_gather(alai_m(:,:,i_mnth) )
             alai_m(:,:,i_mnth)=alai_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1063) title,
     &                    real(alai_m(:,:,i_mnth),kind=4)

             title = 'TerraE: Water stored on the canopy [m]'
             call array_gather(canopyH2O_m(:,:,i_mnth) )
             canopyH2O_m(:,:,i_mnth)=canopyH2O_m(:,:,i_mnth)
     &                               /n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1064) title,
     &                    real(canopyH2O_m(:,:,i_mnth),kind=4)

             title = 'TerraE: Heat of the canopy [J/m2]'
             call array_gather(canopyheat_m(:,:,i_mnth) )
             canopyheat_m(:,:,i_mnth)=canopyheat_m(:,:,i_mnth)
     &                                /n_count_mnth(i_mnth)
             if ( mype == 0 ) write(1065) title,
     &                    real(canopyheat_m(:,:,i_mnth),kind=4)


          enddo
#else
!     Annual monthly diagnostic accumulators
          do i_mnth=1,12
              title = 'TerraE: Total evaporation'
              write(980) title,real(aevap_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Canopy evaporation (evapw)'
              write(981) title,real(aevapw_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Transpiration (evapd)'
              write(982) title,real(aevapd_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Bare soil evaporation (evapb)'
              write(983) title,real(aevapb_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Canopy interception'
              write(984) title,real(aintercep_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Surface runoff (runs)'
              write(985) title,real(aruns_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Subsurface flow out of soil column'
              write(986) title,real(arunu_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Gross primary productivity'
              write(987) title,real(agpp_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Autotrophic respiration'
              write(988) title,real(arauto_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Soil respiration'
              write(989) title,real(asoilresp_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Evap. from vegetated soil'
              write(990) title,real(aevapvg_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Evap. from snow on vegetation'
              write(991) title,real(aevapvs_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Evap. from snow on bare soil'
              write(992) title,real(aevapbs_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Net shortwave radiation [W/m2]'
              asrht_mnth(:,:,i_mnth) = 
     &                asrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(993) title,real(asrht_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Net longwave radiation [W/m2]'
              atrht_mnth(:,:,i_mnth) = 
     &                atrht_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(994) title,real(atrht_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Grid-cell mean albedo'
              aalbedo_mnth(:,:,i_mnth) = 
     &                aalbedo_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(995) title,real(aalbedo_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Plant labile carbon [kgC/m2]'
              aclab_mnth(:,:,i_mnth) = 
     &                aclab_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(996) title,real(aclab_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Sum of soil carbon pools[kgC/m2]'
              asoilCpoolsum_mnth(:,:,i_mnth) = 
     &          asoilCpoolsum_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(997) title,real(asoilCpoolsum_mnth(:,:,i_mnth)
     &                              ,kind=4)
              title = 'TerraE: Penman potential evap. [kg/m2/month]'
              write(998) title,real(aepp_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Thermal heat from ground [J/m2/mnth]'
              write(999) title,real(atrg_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Sensible heat from ground [J/m2/mnth]'
              write(1000) title,real(ashg_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Latent heat flux [W/m2]'
              alhg_mnth(:,:,i_mnth) = 
     &                alhg_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
              write(1001) title,real(alhg_mnth(:,:,i_mnth),kind=4)
              title = 'TerraE: Pot. evap. from canopy [kg/m2/mnth]'
              write(1002) title,real(aepc_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Bare layer 1 [m]'
             w_b1_mnth(:,:,i_mnth) = 
     &                w_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1003) title,real(w_b1_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Bare layer 2 [m]'
             w_b2_mnth(:,:,i_mnth) = 
     &                w_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1004) title,real(w_b2_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Bare layer 3 [m]'
             w_b3_mnth(:,:,i_mnth) = 
     &                w_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1005) title,real(w_b3_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Bare layer 4 [m]'
             w_b4_mnth(:,:,i_mnth) = 
     &                w_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1006) title,real(w_b4_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Bare layer 5 [m]'
             w_b5_mnth(:,:,i_mnth) = 
     &                w_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1007) title,real(w_b5_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Bare layer 6 [m]'
             w_b6_mnth(:,:,i_mnth) = 
     &                w_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1008) title,real(w_b6_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Veg layer 1 [m]'
             w_v1_mnth(:,:,i_mnth) = 
     &                w_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1009) title,real(w_v1_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Veg layer 2 [m]'
             w_v2_mnth(:,:,i_mnth) = 
     &                w_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1010) title,real(w_v2_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Veg layer 3 [m]'
             w_v3_mnth(:,:,i_mnth) = 
     &                w_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1011) title,real(w_v3_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Veg layer 4 [m]'
             w_v4_mnth(:,:,i_mnth) = 
     &                w_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1012) title,real(w_v4_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Veg layer 5 [m]'
             w_v5_mnth(:,:,i_mnth) = 
     &                w_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1013) title,real(w_v5_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil water; Veg layer 6 [m]'
             w_v6_mnth(:,:,i_mnth) = 
     &                w_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1014) title,real(w_v6_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Bare layer 1 [J/m2]'
             ht_b1_mnth(:,:,i_mnth) = 
     &                ht_b1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1015) title,real(ht_b1_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Bare layer 2 [J/m2]'
             ht_b2_mnth(:,:,i_mnth) = 
     &                ht_b2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1016) title,real(ht_b2_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Bare layer 3 [J/m2]'
             ht_b3_mnth(:,:,i_mnth) = 
     &                ht_b3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1017) title,real(ht_b3_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Bare layer 4 [J/m2]'
             ht_b4_mnth(:,:,i_mnth) = 
     &                ht_b4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1018) title,real(ht_b4_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Bare layer 5 [J/m2]'
             ht_b5_mnth(:,:,i_mnth) = 
     &                ht_b5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1019) title,real(ht_b5_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Bare layer 6 [J/m2]'
             ht_b6_mnth(:,:,i_mnth) = 
     &                ht_b6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1020) title,real(ht_b6_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Veg layer 1 [J/m2]'
             ht_v1_mnth(:,:,i_mnth) = 
     &                ht_v1_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1021) title,real(ht_v1_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Veg layer 2 [J/m2]'
             ht_v2_mnth(:,:,i_mnth) = 
     &                ht_v2_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1022) title,real(ht_v2_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Veg layer 3 [J/m2]'
             ht_v3_mnth(:,:,i_mnth) = 
     &                ht_v3_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1023) title,real(ht_v3_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Veg layer 4 [J/m2]'
             ht_v4_mnth(:,:,i_mnth) = 
     &                ht_v4_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1024) title,real(ht_v4_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Veg layer 5 [J/m2]'
             ht_v5_mnth(:,:,i_mnth) = 
     &                ht_v5_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1025) title,real(ht_v5_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Soil heat; Veg layer 6 [J/m2]'
             ht_v6_mnth(:,:,i_mnth) = 
     &                ht_v6_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1026) title,real(ht_v6_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow layer thickness; bare layer 1 [m]'
             dzsn_b1(:,:,i_mnth) = 
     &                dzsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1027) title,real(dzsn_b1(:,:,i_mnth),kind=4)
             title = 'TerraE:  Snow layer thickness; bare layer 2 [m]'
             dzsn_b2(:,:,i_mnth) = 
     &                dzsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1028) title,real(dzsn_b2(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow layer thickness; bare layer 3 [m]'
             dzsn_b3(:,:,i_mnth) = 
     &                dzsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1029) title,real(dzsn_b3(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow layer thickness; bare layer 1 [m]'
             dzsn_v1(:,:,i_mnth) = 
     &                dzsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1030) title,real(dzsn_v1(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow layer thickness; bare layer 2 [m]'
             dzsn_v2(:,:,i_mnth) = 
     &                dzsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1031) title,real(dzsn_v2(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow layer thickness; bare layer 3 [m]'
             dzsn_v3(:,:,i_mnth) = 
     &                dzsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1032) title,real(dzsn_v3(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow water equival; bare layer 1 [m]'
             wsn_b1(:,:,i_mnth) = 
     &                wsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1033) title,real(wsn_b1(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow water equival; bare layer 2 [m]'
             wsn_b2(:,:,i_mnth) = 
     &                wsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1034) title,real(wsn_b2(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow water equival; bare layer 3 [m]'
             wsn_b3(:,:,i_mnth) = 
     &                wsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1035) title,real(wsn_b3(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow water equival; veg layer 1 [m]'
             wsn_v1(:,:,i_mnth) = 
     &                wsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1036) title,real(wsn_v1(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow water equival; veg layer 2 [m]'
             wsn_v2(:,:,i_mnth) = 
     &                wsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1037) title,real(wsn_v2(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow water equival; veg layer 3 [m]'
             wsn_v3(:,:,i_mnth) = 
     &                wsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1038) title,real(wsn_v3(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow layer heat; bare layer 1 [J/m2]'
             hsn_b1(:,:,i_mnth) = 
     &                hsn_b1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1039) title,real(hsn_b1(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow layer heat; bare layer 1 [J/m2]'
             hsn_b2(:,:,i_mnth) = 
     &                hsn_b2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1040) title,real(hsn_b2(:,:,i_mnth),kind=4)
             title = 'TerraE:  Snow layer heat; bare layer 1 [J/m2]'
             hsn_b3(:,:,i_mnth) = 
     &                hsn_b3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1041) title,real(hsn_b3(:,:,i_mnth),kind=4)
             title = 'TerraE:  Snow layer heat; veg layer 1 [J/m2]'
             hsn_v1(:,:,i_mnth) = 
     &                hsn_v1(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1042) title,real(hsn_v1(:,:,i_mnth),kind=4)
             title = 'TerraE:  Snow layer heat; veg layer 2 [J/m2]'
             hsn_v2(:,:,i_mnth) = 
     &                hsn_v2(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1043) title,real(hsn_v2(:,:,i_mnth),kind=4)
             title = 'TerraE:  Snow layer heat; veg layer 3 [J/m2]'
             hsn_v3(:,:,i_mnth) = 
     &                hsn_v3(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1044) title,real(hsn_v3(:,:,i_mnth),kind=4)
             title = 'TerraE: # of snow layers on bare land'
             nsn_b(:,:,i_mnth) = 
     &                nsn_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1045) title,real(nsn_b(:,:,i_mnth),kind=4)
             title = 'TerraE: # of snow layers on vegetated land'
             nsn_v(:,:,i_mnth) = 
     &                nsn_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1046) title,real(nsn_v(:,:,i_mnth),kind=4)
             title = 'TerraE: Snow-covered bare fraction [-]'
             frsnow_b(:,:,i_mnth) = 
     &                frsnow_b(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1047) title,real(frsnow_b(:,:,i_mnth),kind=4)
             title = 'TerraE:  Snow-covered vegetated fraction [-]'
             frsnow_v(:,:,i_mnth) = 
     &                frsnow_v(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1048) title,real(frsnow_v(:,:,i_mnth),kind=4)
             title = 'TerraE:Foliage surf. vapor mixing ratio [kg/kg]'
             Qf_mnth(:,:,i_mnth) = 
     &                Qf_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1049) title,real(Qf_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE:root zone betad [-]'
             abetad_mnth(:,:,i_mnth) = 
     &                abetad_mnth(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1050) title,real(abetad_mnth(:,:,i_mnth),kind=4)
             title = 'TerraE:Live leaf carbon pool [kg/m2]'
             aClivepool_leaf_m(:,:,i_mnth) = 
     &              aClivepool_leaf_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1051) title,
     &              real(aClivepool_leaf_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Live froot carbon pool [kg/m2]'
             aClivepool_froot_m(:,:,i_mnth) = 
     &              aClivepool_froot_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1052) title,
     &                    real(aClivepool_froot_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Live wood carbon pool [kg/m2]'
             aClivepool_wood_m(:,:,i_mnth) = 
     &              aClivepool_wood_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1053) title,
     &                    real(aClivepool_wood_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Dead surface metabolic C pool [kg/m2]'
             aCdeadpool_surfmet_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1054) title,
     &                    real(aCdeadpool_surfmet_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Dead surface structural C pool [kg/m2]'
             aCdeadpool_surfstr_m(:,:,i_mnth) = 
     &           aCdeadpool_surfstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1055) title,
     &                    real(aCdeadpool_surfstr_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Dead soil metabolic C pool [kg/m2]'
             aCdeadpool_soilmet_m(:,:,i_mnth) = 
     &           aCdeadpool_soilmet_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1056) title,
     &                    real(aCdeadpool_soilmet_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Dead soil structural C pool [kg/m2]'
             aCdeadpool_soilstr_m(:,:,i_mnth) = 
     &           aCdeadpool_soilstr_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1057) title,
     &                    real(aCdeadpool_soilstr_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Coarse woody debris C pool [kg/m2]'
             aCdeadpool_cwd_m(:,:,i_mnth) = 
     &            aCdeadpool_cwd_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1058) title,
     &                    real(aCdeadpool_cwd_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Dead surface microbial C pool [kg/m2]'
             aCdeadpool_surfmic_m(:,:,i_mnth) = 
     &            aCdeadpool_surfmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1059) title,
     &                    real(aCdeadpool_surfmic_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Dead soil microbial C pool [kg/m2]'
             aCdeadpool_soilmic_m(:,:,i_mnth) = 
     &            aCdeadpool_soilmic_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1060) title,
     &                    real(aCdeadpool_soilmic_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Dead slow (~10 yr) C pool [kg/m2]'
             aCdeadpool_slow_m(:,:,i_mnth) = 
     &            aCdeadpool_slow_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1061) title,
     &                    real(aCdeadpool_slow_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Dead passive (~100 yr) C pool [kg/m2]'
             aCdeadpool_passive_m(:,:,i_mnth) = 
     &            aCdeadpool_passive_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1062) title,
     &                    real(aCdeadpool_passive_m(:,:,i_mnth),kind=4)
             title = 'TerraE:Leaf Area Index [-]'
             alai_m(:,:,i_mnth)=alai_m(:,:,i_mnth)/n_count_mnth(i_mnth)
             write(1063) title,real(alai_m(:,:,i_mnth),kind=4)

             title = 'TerraE: Water stored on the canopy [m]'
             canopyH2O_m(:,:,i_mnth)=canopyH2O_m(:,:,i_mnth)
     &                               /n_count_mnth(i_mnth)
             write(1064) title,real(canopyH2O_m(:,:,i_mnth),kind=4)

             title = 'TerraE: Heat of the canopy [J/m2]'
             canopyheat_m(:,:,i_mnth)=canopyheat_m(:,:,i_mnth)
     &                                /n_count_mnth(i_mnth)
             write(1065) title,real(canopyheat_m(:,:,i_mnth),kind=4)

          enddo
#endif
!         Reset annual monthly diagnostic/prognostic accumulators
          aevap_mnth(:,:,:) = 0.d0
          aevapw_mnth(:,:,:) = 0.d0
          aevapd_mnth(:,:,:) = 0.d0
          aevapb_mnth(:,:,:) = 0.d0
          aintercep_mnth(:,:,:) = 0.d0
          aruns_mnth(:,:,:) = 0.d0
          arunu_mnth(:,:,:) = 0.d0
          agpp_mnth(:,:,:) = 0.d0
          arauto_mnth(:,:,:) = 0.d0
          asoilresp_mnth(:,:,:) = 0.d0
          aevapvg_mnth(:,:,:) = 0.d0
          aevapvs_mnth(:,:,:) = 0.d0
          aevapbs_mnth(:,:,:) = 0.d0
          asrht_mnth(:,:,:) = 0.d0
          atrht_mnth(:,:,:) = 0.d0
          aalbedo_mnth(:,:,:) = 0.d0
          aclab_mnth(:,:,:) = 0.d0
          asoilCpoolsum_mnth(:,:,:) = 0.d0
          aepp_mnth(:,:,:) = 0.d0
          atrg_mnth(:,:,:) = 0.d0
          ashg_mnth(:,:,:) = 0.d0
          alhg_mnth(:,:,:) = 0.d0
          aepc_mnth(:,:,:) = 0.d0
          w_b1_mnth(:,:,:)=0.d0
          w_b2_mnth(:,:,:)=0.d0
          w_b3_mnth(:,:,:)=0.d0
          w_b4_mnth(:,:,:)=0.d0
          w_b5_mnth(:,:,:)=0.d0
          w_b6_mnth(:,:,:)=0.d0
          w_v1_mnth(:,:,:)=0.d0
          w_v2_mnth(:,:,:)=0.d0
          w_v3_mnth(:,:,:)=0.d0
          w_v4_mnth(:,:,:)=0.d0
          w_v5_mnth(:,:,:)=0.d0
          w_v6_mnth(:,:,:)=0.d0
          ht_b1_mnth(:,:,:)=0.d0
          ht_b2_mnth(:,:,:)=0.d0
          ht_b3_mnth(:,:,:)=0.d0
          ht_b4_mnth(:,:,:)=0.d0
          ht_b5_mnth(:,:,:)=0.d0
          ht_b6_mnth(:,:,:)=0.d0
          ht_v1_mnth(:,:,:)=0.d0
          ht_v2_mnth(:,:,:)=0.d0
          ht_v3_mnth(:,:,:)=0.d0
          ht_v4_mnth(:,:,:)=0.d0
          ht_v5_mnth(:,:,:)=0.d0
          ht_v6_mnth(:,:,:)=0.d0
          dzsn_b1(:,:,:)=0.d0
          dzsn_b2(:,:,:)=0.d0
          dzsn_b3(:,:,:)=0.d0
          dzsn_v1(:,:,:)=0.d0
          dzsn_v2(:,:,:)=0.d0
          dzsn_v3(:,:,:)=0.d0
          wsn_b1(:,:,:)=0.d0
          wsn_b2(:,:,:)=0.d0
          wsn_b3(:,:,:)=0.d0
          wsn_v1(:,:,:)=0.d0
          wsn_v2(:,:,:)=0.d0
          wsn_v3(:,:,:)=0.d0
          hsn_b1(:,:,:)=0.d0
          hsn_b2(:,:,:)=0.d0
          hsn_b3(:,:,:)=0.d0
          hsn_v1(:,:,:)=0.d0
          hsn_v2(:,:,:)=0.d0
          hsn_v3(:,:,:)=0.d0
          nsn_b(:,:,:)=0.d0
          nsn_v (:,:,:)=0.d0
          frsnow_b(:,:,:)=0.d0
          frsnow_v(:,:,:)=0.d0
          Qf_mnth(:,:,:)=0.d0
          abetad_mnth(:,:,:)=0.d0
          aClivepool_leaf_m(:,:,:)=0.d0
          aClivepool_froot_m(:,:,:)=0.d0
          aClivepool_wood_m(:,:,:)=0.d0
          aCdeadpool_surfmet_m(:,:,:)=0.d0
          aCdeadpool_surfstr_m(:,:,:)=0.d0
          aCdeadpool_soilmet_m(:,:,:)=0.d0
          aCdeadpool_soilstr_m(:,:,:)=0.d0
          aCdeadpool_cwd_m(:,:,:)=0.d0
          aCdeadpool_surfmic_m(:,:,:)=0.d0
          aCdeadpool_soilmic_m(:,:,:)=0.d0
          aCdeadpool_slow_m(:,:,:)=0.d0
          aCdeadpool_passive_m(:,:,:)=0.d0
          alai_m(:,:,:)=0.d0
          canopyH2O_m(:,:,:)=0.d0
          canopyheat_m(:,:,:)=0.d0
          n_count_mnth(:) = 0.d0 !Reset counter
#endif
        endif
      enddo

      call sysusage(mype,2)
      call sysusage(mype,3)


#ifdef USE_ESMF
      call MPI_Finalize(ierr)
#endif

      end subroutine terraE_standalone
