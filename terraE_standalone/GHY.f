c**** SLE001 E001M12 SOMTQ SLB211M9
c**** (same as frank''s soils64+2bare_soils+old runoff)
c**** change to evap calculation to prevent negative runoff
c**** soils45 but with snowmelt subroutine snmlt changed
c**** to melt snow before 1st layer ground ice.
ccc   comments from new soils
c**** 8/11/97 - modified to include snow model
c**** 10/29/96 - aroot and broot back to original values; added
c**** call to cpars to change vegetation parameters for pilps.
c**** 9/7/96 - back to heat capacities/10
c**** 5/14/96 - added soils64 surface runoff changes/corrections
c**** 11/13/95 - changed aroot and broot for rainf for 1.5m root depth
c**** 10/11/95 - back to full heat capacity, to avoid cd oscillations.
c**** changes for pilps: (reversed)
c**** use soils100.com instead of soils45.com
c**** set im=36,jm=24
c**** set sdata,vdata and fdata to real*4
c**** divide canopy heat capacities by 10.d0
c**** change aroot of grass to 1.0d0, to reduce root depth.
c**** end of changes for pilps
c**** changes for pilps: (kept)
c**** modify gdtm surface flux timestep limits
c**** define new diagnostics
c**** zero out diagnostics at start of advnc
c**** end of changes for pilps
c**** modified for 2 types of bare soils
c****
c**** soils62 soils45 soils45          cdfxa 04/27/95
c**** same as soils45 but with snowmelt subroutine snmlt changed
c**** to melt snow before 1st layer ground ice.
ccc   end comments from new soils
c**** also corrects evaps calculation.
c**** also includes masking effects in radation fluxes.
c**** modifies timestep for canopy fluxes.
c**** soils45 10/4/93
c**** uses bedrock as a soil texture.  soil depth of 3.5m
c**** everywhere, where layers can have bedrock.
c**** requires sm693.data instead of sm691.data.
c**** sdata needs to be changed in calling program.
c**** soils44b 8/25/93
c**** uses snow conductivity of .088 w m-1 c-1 instead of .3d0
c**** soils44 8/16/93
c**** adds bedrock for heat calculations, to fill out the
c**** number of layers to ngm.
c**** soils43 6/25/93
c**** comments out call to fhlmt heat flux limits.
c**** uses ghinij to return wfc1, eliminates rewfc.
c**** soils42 6/15/93
c**** adds snow insulation
c**** soils41 5/24/93
c**** uses snow masking depth from vegetation height to determine
c**** fraction of snow that is exposed.
c**** reth must be called prior to retp.
c**** soils40 5/10/93
c**** removes snow from canopy and places it on vegetated soil.
c**** soils39 4/19/93
c**** modifications for real*8 or real*4 runs.  common block
c**** ordering changed for efficient alignment.  sdata,fdata,
c**** and vdata are explicitly real*4.  on ibm rs/6000, should
c**** be compiled with -qdpc=e option for real*8 operation.
c**** to run real*4, change implicit statement in include file.
c**** soils38 2/9/93
c**** adds heat flux correction to handle varying coefficients
c**** of drag.
c**** soils37 1/25/93
c**** changes soil crusting parameter ku/d from .05 per hour to .1d0,
c**** to agree with morin et al.
c**** soils36 11/12/92
c**** calculates heat conductivity of soils using devries method.
c**** changes loam material heat capacity and conductivity
c**** to mineral values.
c**** soils35 10/27/92
c**** includes effect of soil crusting for infiltration by
c**** modifying hydraulic conductivity calculation of layer
c**** 1 in hydra.
c**** soils34 8/28/92
c**** uses effective leaf area index alaie for purposes of
c**** canopy conductance calculation.
c**** soils33 8/9/92
c**** changes canopy water storage capacity to .1mm per lai from 1.d0
c**** soils32 7/15/92
c**** 1) for purposes of infiltration only, reduces soil conductivity
c**** by (1-thetr*fice) instead of (1-fice).
c**** 2) betad is reduced by fraction of ice in each layer.
c**** 3) transpired water is removed by betad fraction in each layer,
c**** instead of by fraction of roots. prevents negative runoff.
c**** 4) speeds up hydra by using do loop instead of if check,
c**** by using interpolation point from bisection instead of logs,
c**** and by avoiding unnecessary calls to hydra.  also elimates call
c**** to hydra in ma89ezm9.f.
c**** soils31 7/1/92
c**** 1) fixes fraction of roots when soil depth is less than root
c**** depth, thus fixing non-conservation of water.
c**** soils30 6/4/92
c**** 1) uses actual final snow depth in flux limit calculations,
c**** instead of upper and lower limits.  fixes spurious drying
c**** of first layer.
c**** Added gpp, GPP terms  4/25/03 nyk
c**** Added decks parameter cond_scheme  5/1/03 nyk
c**** Added decks parameter vegCO2X_off  3/2/04 nyk

#include "rundeck_opts.h"

!#define EVAP_VEG_GROUND
!#define RAD_VEG_GROUND

      module sle001

      use constant, only : stbo,tfrz=>tf,sha,lhe,one,zero,rhow
     &     ,shw_kg=>shw,shi_kg=>shi,lhm
      use ghy_com, only : ngm, imt, nlsn, LS_NFRAC !, wsn_max
#ifdef TRACERS_WATER
      use tracer_com, only : ntm
#endif

      implicit none
      save


      private

c**** public functions:
      public hl0, set_snow, advnc, evap_limits, xklh

ccc   physical constants and global model parameters
ccc   converting constants from 1/kg to 1/m^3
!@var shw heat capacity of water (J/m^3 C)
      real*8, parameter, public :: shw= shw_kg * rhow
!@var shi heat capacity of pure ice (J/m^3 C)
      real*8, parameter, public :: shi= shi_kg * rhow
!@var fsn latent heat of melt (J/m^3)
      real*8, parameter, public :: fsn= lhm * rhow
!@var elh latent heat of evaporation (J/m^3)
      real*8, parameter :: elh= lhe * rhow
!@var prfr fraction (by area) of precipitation
      real*8, parameter :: prfr = 0.1d0

ccc   data needed for debugging

      type, public :: ghy_debug_str
        real*8 water(2), energy(2)
      end type ghy_debug_str

      type ( ghy_debug_str ), public :: ghy_debug
      integer, public :: ijdebug

ccc   public data
ccc   main accumulators
      real*8, public :: atrg,ashg,aevap,alhg,aruns,arunu,aeruns,aerunu
!@var tbcs surface temperature (C)
      real*8, public ::  tbcs

ccc   diagnostics accumulatars
      real*8, public :: aevapw,aevapd,aevapb,aepc,aepb,aepp,af0dt,af1dt
     &      , agpp,aflmlt,aintercep
ccc   some accumulators that are currently not computed:
     &     ,acna,acnc
ccc   beta''s
     &     ,abetad,abetav,abetat,abetap,abetab,abeta
!@var zw water table (not computed ?)
      real*8, public :: zw(2)

ccc   private accumulators:
      real*8 aedifs

ccc   input fluxes
      real*8, public :: pr,htpr,prs,htprs,srht,trht
ccc   input bc''s
      real*8, public :: ts,qs,pres,rho,vsm,ch,qm1
!@var dt earth time step (s)
      real*8, public :: dt

ccc   soil prognostic variables
!!!      integer, parameter, public :: ngm=6, ng=ngm+1, imt=5
      integer, parameter, public :: ng=ngm+1
      real*8, public :: w(0:ngm,LS_NFRAC),ht(0:ngm,LS_NFRAC),dz(ngm)

ccc   soil properties which need to be passed in:
      real*8, public :: q(imt,ngm),qk(imt,ngm),sl,fv,fb

ccc   topographic info which is passed in
      real*8, public :: top_index, top_stdev

ccc   soil internal variables wchich need to be passed from/to ghinij
      real*8, public :: zb(ng),zc(ng),fr(ng),snowm !veg alaie, rs,
     &     ,thets(0:ng,2),thetm(0:ng,2),ws(0:ngm,2),shc(0:ng,2)
     &     ,thm(0:64,imt-1)
!veg     &     ,nm,nf,vh ! added by adf ! ,alai
!@var n the worst choice for global name: number of soil layers
!@var nth some magic number computed in hl0, used in ghinij and hydra
      integer, public :: n,nth

ccc   soil internal variables wchich need to be passed to driver anyway
      real*8, public :: snowd(2),tp(0:ngm,2),fice(0:ngm,2)
      !tp(:,:) is the temperature (C) of soil layers numbered from top
      !to bottom for both bare and vegetated fractions.
      !tp(:,1) - temperature for bare soil fraction
      !tp(:,2) - temperature for vegetated fraction

      !tp(0,2) - canopy
      !tp(0,1) - not used
      !tp(1,:) - upper soil layer
      !tp(2,:) - second soil layer

ccc   snow prognostic variables
!!!      integer, parameter, public :: nlsn=3
      integer, public :: nsn(2)
      real*8, public ::  dzsn(nlsn+1,2),wsn(nlsn,2),hsn(nlsn,2)
     &     ,fr_snow(2)
ccc   tsn1 is private
!@var tsn1 temperature of the upper layer of snow (C)
      real*8 tsn1(2)

      real*8 betat,betad
      real*8 gpp,dts,trans_sw
!xxx      real*8, public :: cnc

ccc fractions of dry,wet,covered by snow canopy
!@var fd effective fraction of dry canopy (=0 if dew)
!@var fv effective fraction of wet canopy (=1 if dew)
!@var fd0 actual fraction of dry canopy
!@var fv0 actual fraction of wet canopy
!@var fm snow masking fraction for canopy (=1 if all snow)
      real*8 fd,fw,fd0,fw0,fm

!@var theta fraction of water in the soil ( 1 = 100% )
!@var f water flux between the layers ( > 0 is up ) (m/s)
!@var fh heat flux between the layers ( > 0 is up ) (W)
!@var rnf surface runoff (m/s)
!@var rnff underground runoff (m/s)
!@var fc water flux at canopy boundaries ( > 0 is up ) (m/s)
!@var fch heat flux at canopy boundaries ( > 0 is up ) (W)
      real*8 theta(0:ng,2),f(1:ng,2),fh(1:ng,2),rnf(2),rnff(ng,2)
     &     ,fc(0:1),fch(0:1)
ccc   the following three declarations are various conductivity
ccc   coefficients (or data needed to compute them)
      real*8 xkh(ng,2),xkhm(ng,2),h(0:ng,2) ,xk(0:ng,2),xinfc(2)
     &     ,d(0:ng,2) ,xku(0:ng,2)
      real*8 hlm(0:64), xklm(0:64,imt-1),dlm(0:64,imt-1)
      real*8 xkus(ng,2), xkusa(2)


ccc   the following are fluxes computed inside the snow model
ccc   they are unit/(m^2 of snow), not multiplied by aby fractions
!@var flmlt flux of water down from snow to the ground (m/m^2)
!@var fhsng flux of heat down from snow to the ground (J/m^2)
!@var thrmsn flux of radiative heat up from snow (J/m^2)
      real*8  flmlt(2),fhsng(2),thrmsn(2)
ccc   similar values defined for large scale, i.e. flux to entire cell
ccc   total flux is val*fr_snow + value_scale
      real*8  flmlt_scale(2),fhsng_scale(2)


ccc   put here some data needed for print-out during the crash
      type :: debug_data_str
        sequence
        real*8 qc,qb
      end type debug_data_str

      type (debug_data_str) debug_data

ccc   evaporation limiting data which one needs to pass the pbl
      real*8, public :: evap_max_sat, evap_max_nsat, fr_sat

ccc   potential evaporation
!@var epb,epbs,epv,epvs potential evaporation (b-bare, v-vege, s-snow)
      real*8 epb, epv, epbs, epvs, epvg

ccc   evaporation fluxes: these are pure fluxes over m^2
ccc   i.e. they are not multiplied by any fr_...
!@var evapb, evapbs, evapvw, evapvd, evapvs, evapvg evaporation flux (m/s)
!@+ (b=bare, v=vege, d=dry, w=wet, s=snow, g=ground)
      real*8 :: evapb, evapbs, evapvw, evapvd, evapvs, evapvg
!@var devapbs_dt, devapvs_dt d(evap)/dT for snow covered soil (bare/veg)
      real*8 :: devapbs_dt, devapvs_dt
c!@var evapor mean (weighted with fr_..) evaporation from bare/vege soil
ccc      real*8 :: evapor(2)

!@var evapdl thanspiration from each soil layer
      real*8 :: evapdl(ngm,2)

ccc sensible heat fluxes: these are pure fluxes over m^2
ccc   i.e. they are not multiplied by any fr_...
!@var snsh sensible heat from soil to air (bare/vege) no snow (J/s)
!@var snshs sensible heat from snow to air (bare/vege)  (J/s)
!@var dsnsh_dt d(sensible heat)/d(soil temp) : same for bare/vege/snow
      real*8 :: snsh(2), snshs(2), dsnsh_dt

!!!! vars below this line were not included into GHYTPC yet !!!

!@var dripw liquid water flux from canopy (similar for bare soil) (m/s)
!@var htdripw heat of drip water (similar for bare soil) (J/s)
!@var drips snow flux from canopy (similar for bare soil) (m/s)
!@var htdrips heat of drip snow (similar for bare soil) (J/s)
      real*8 dripw(2), htdripw(2), drips(2), htdrips(2)

!@var evap_tot total evap. (weighted with fr_snow,fm)(bare/veg)(m/s)
!@var snsh_tot total sens. heat(weighted with fr_snow,fm)(bare/veg)(m/s)
!@var thrm_tot total T^4 heat (weighted with fr_snow,fm)(bare/veg)(m/s)
      real*8 evap_tot(2), snsh_tot(2), thrm_tot(2)
      public evap_tot ! for debug only !

ccc data for tracers
!@var flux_snow water flux between snow layers (>0 down) (m/s)
!@var wsn_for_tr snow water before call to snow_adv
      real*8 flux_snow(0:nlsn,2), wsn_for_tr(nlsn,2)

#ifdef TRACERS_WATER
ccc should be passed from elsewhere
!@var ntgm maximal number of tracers that can by passesd to HGY
      integer, parameter, public :: ntgm = ntm
!@var ntg actual number of tracers passed to ground hydrology
      integer, public :: ntg
!@var trpr flux of tracers in precipitation (?/m^2 s)
!@var trdd flux of tracers as dry deposit (?/m^2 s)
!@var tr_w amount of tracers in the soil (?)
!@var tr_wsn amount of tracers in the snow (?)
      real*8, public :: trpr(ntgm), trdd(ntgm), tr_surf(ntgm),
     &     tr_w(ntgm,0:ngm,2), tr_wsn(ntgm,nlsn,2)
ccc tracers output:
!@var tr_evap flux of tracers to atm due to evaporation (?/m^2 s)
!@var tr_evap flux of tracers due to runoff (?/m^2 s)
      real*8 tr_evap(ntgm,2),tr_rnff(ntgm,2)
!@var tr_evap flux of tracers to atm due to evaporation
!@var tr_evap flux of tracers due to runoff
      real*8, public :: atr_evap(ntgm),atr_rnff(ntgm),atr_g(ntgm)
#ifdef TRACERS_SPECIAL_O18
      character*8, public :: tr_name(ntgm)
#endif
#endif

ccc the data below this line is not in GHYTPC yet !
ccc the following vars control if bare/vegetated fraction has to
ccc be computed (i.e. f[bv] is not zero)
      integer :: i_bare, i_vege
      logical :: process_bare, process_vege

C***
C***   Thread Private Common Block GHYTPC
C***
      COMMON /GHYTPC/
     &     abeta,abetab,abetad,abetap,abetat,abetav,acna,acnc,agpp
     &     ,aedifs,aepb,aepc,aepp,aeruns,aerunu,aevap,aevapb
     &     ,aevapd,aevapw,af0dt,af1dt,alhg,aruns,arunu,aflmlt,aintercep
     &     ,ashg,atrg,betad,betat,ch,gpp,d,devapbs_dt,devapvs_dt
     &     ,drips,dripw,dsnsh_dt,dts,dz,dzsn,epb,epbs,epvs,epvg  ! dt dlm
     &     ,epv,evap_max_nsat,evap_max_sat,evap_tot,evapb
     &     ,evapbs,evapdl,evapvd,evapvs,evapvw,evapvg,f !evapor,
     &     ,fb,fc,fch,fd,fd0,fh,fhsng,fhsng_scale,fice,flmlt,flmlt_scale
     &     ,fm,fr,fr_sat,fr_snow,fv,fw,fw0,h,hsn,ht !hlm
     &     ,htdrips,htdripw,htpr,htprs,pr,pres,prs,q,qk,qm1,qs
     &     ,rho,rnf,rnff,shc,sl,snowd,snowm,snsh,snsh_tot !veg rs,
     &     ,snshs,srht,tbcs,theta,thetm,thets,thrm_tot,thrmsn !thm
     &     ,top_index,top_stdev,tp,trht,ts,tsn1,vsm,w,ws,wsn,xinfc,xk
     &     ,xkh,xkhm,xku,xkus,xkusa,zb,zc,zw ! xklm
     &     ,ijdebug,n,nsn !nth
     &     ,flux_snow,wsn_for_tr,trans_sw

!----------------------------------------------------------------------!
     &     ,i_bare,i_vege,process_bare,process_vege
#ifdef TRACERS_WATER
     &     ,trpr,trdd, tr_surf, tr_w, tr_wsn,tr_evap,tr_rnff ! ntg
     &     ,atr_evap,atr_rnff,atr_g
#endif
c     not sure if it works with derived type. if not - comment the
c     next line out (debug_data used only for debug output)
c    &     ,debug_data           ! needs to go to compile on COMPAQ
!$OMP  THREADPRIVATE (/GHYTPC/)
C***
C***

ccc   external functions
      real*8, external :: qsat,dqsatdt



      contains

      subroutine reth
!@sum computes theta(:,:), snowd(:), fm, fw, fd
c**** revises values of theta based upon w.
c**** input:
c**** w - water depth, m
c**** ws - saturated water depth, m
c**** dz - layer thickness, m
c**** thets - saturated theta
c**** snowd - snow depth, water equivalent m
c**** snowm - snow masking depth, water equivalent m
c**** output:
c**** theta - water saturation
c**** fw - fraction of wet canopy
c**** fd - fraction of dry canopy
c**** fm - fraction of snow that is exposed, or masking.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      integer lsn,ibv,k
      do ibv=i_bare,i_vege
        do  k=1,n
          theta(k,ibv)=w(k,ibv)/dz(k)
        end do
      end do
c**** do canopy layer
c**** here theta is the fraction of canopy covered by water
      if( process_vege .and. ws(0,2).gt.0.d0 )then
        theta(0,2)=(w(0,2)/ws(0,2))**(2.d0/3.d0)
      else
        theta(0,2)=0.d0
      endif
      theta(0,2)=min(theta(0,2),one)
c**** set up snowd variables
      do ibv=i_bare,i_vege
        snowd(ibv)=0.d0
        do lsn=1,nsn(ibv)
ccc    we compute snowd as if all snow was distributed uniformly
ccc    over the cell (i.e. snowd = wsn * sn_frac / 1.d0 )
          snowd(ibv) = snowd(ibv) + wsn(lsn,ibv) * fr_snow(ibv)
        enddo
      enddo
c**** fraction of wet canopy fw
      fw=theta(0,2)
c**** determine fm from snowd depth and masking depth
      fm=1.d0-exp(-snowd(2)/(snowm+1d-12))
!!! testing if setting bounds on fm will make tracers work ...
      if ( fm < 1.d-3 ) fm=0.d0
!      if ( fm > .99d0 ) fm=1.d0
c**** correct fraction of wet canopy by snow fraction
      !fw=fw+fm*(1.d0-fw)  !!! no idea what does this formula mean (IA)
      fd=1.d0-fw
      fw0=fw
      fd0=fd
      return
      end subroutine reth

      subroutine hydra
!@sum computes matr. potential h(k,ibv) and H2O condactivity xk(k,ibv)
c     routine to return the equlibrium value of h in a mixed soil
c     layer.  the h is such that each soil texture has the same
c     value of h, but differing values of theta.
c     hydra also calculates the conductivity xk and diffussivity d.
c**** input:
c**** theta(k,ibv) - volumetric water concentration
c**** thetm(k,ibv) - minimum theta
c**** thets(k,ibv) - maximum theta
c**** nth - number of h0 intervals in table, a power of two.
c**** hlm(j) - table of h values, from 0 (at j=0) to hmin (at j=nth)
c**** thm(j,i) - value of relative theta at hlm(j) in texture i,
c**** ranging between thets(k,ibv)  at j=0 to thetm(k,ibv) at j=nth.
c**** output:
c**** h - potential, m, including both matric and gravitational
c**** d - diffusivity, dl.
c**** xk - conductivity m/s.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c     solve for h using bisection
c     we assume that if j1.lt.j2 then hlm(j1).gt.hlm(j2)
c     and thm(j1,i).gt.thm(j2,i).
c
c     algdel=log(1.d0+alph0)
c
      real*8 d1,d2,dl,hl,temp,thr,thr0,thr1,thr2,xk1,xkl,xklu,xku1,xku2
      real*8 xkud
      integer i,j,ibv,k,ith,j1,j2,jcm,jc
      real*8 dz_total
      xkud=2.78d-5
      jcm=nint(log(float(nth))/log(2.d0))
      do ibv=i_bare,i_vege
        xk(n+1,ibv)=0.0d0
        xku(0,ibv)=0.d0
        do k=1,n
          j1=0
          j2=nth
          thr1=thets(k,ibv)
          thr2=thetm(k,ibv)
          thr0=theta(k,ibv)
          thr0=min(thr1,thr0)
          thr0=max(thr2,thr0)
          do jc=1,jcm
            j=(j1+j2)/2
            thr=0.d0
            do i=1,imt-1
              thr=thr+thm(j,i)*q(i,k)
            end do
            if(thr-thr0 .lt. 0.d0) then
c     here thr is too small, bisect on low j end
              j2=j
              thr2=thr
            else if (thr-thr0 .gt. 0.d0) then
c     here thr is too large, bisect on high j end
              j1=j
              thr1=thr
            else                ! i.e. .eq.
c     here thr is equal to thr0
              hl=hlm(j)
              j1=j
              thr1=thr0
c     the strange value for thr2 below is only for calculating temp
              thr2=-10.d0
              go to 500
            end if
          end do                ! jc
c     here theta is between two adjacent thr''s. interpolate.
          hl=(hlm(j1)*(thr0-thr2)+hlm(j2)*(thr1-thr0))/(thr1-thr2)
 500      continue
c**** only filling hl array with matric potential (gravitational to be
c**** added later)
          h(k,ibv)=hl
c**** calculate diffusivity
          ith=j1
          temp=(thr1-thr0)/(thr1-thr2)
          d1=0.d0
          d2=0.d0
          xku1=0.d0
          xku2=0.d0
          xkus(k,ibv) = 0.d0
          do i=1,imt-1
            d1=d1+q(i,k)*dlm(ith,i)
            d2=d2+q(i,k)*dlm(ith+1,i)
            xku1=xku1+q(i,k)*xklm(ith,i)
            xku2=xku2+q(i,k)*xklm(ith+1,i)
            xkus(k,ibv) = xkus(k,ibv) + q(i,k)*xklm(0,i)
          end do
          dl=(1.d0-temp)*d1+temp*d2
          dl=(1.d0-fice(k,ibv))*dl
          d(k,ibv)=dl
c**** calculate conductivity
          xklu=(1.d0-temp)*xku1+temp*xku2
          xklu=(1.d0-fice(k,ibv))*xklu
          xku(k,ibv)=xklu
          if(k.eq.1) then
            xk1=0.d0
            do i=1,imt-1
              xk1=xk1+qk(i,1)*xklm(0,i)
            end do
            xkl=xk1
            xkl=xkl/(1.d0+xkl/(-zc(1)*xkud))
            xkl=(1.d0-fice(1,ibv)*theta(1,ibv)/thets(1,ibv))*xkl
            xkl=max(zero,xkl)
            xk(1,ibv)=xkl
          else
            xk(k,ibv)=sqrt(xku(k-1,ibv)*xku(k,ibv))
          end if
        end do                  ! k
      end do                    ! ibv
ccc compute conductivity for topmodel (i.e. mean saturated conductivity)
      do ibv=i_bare,i_vege
        xkusa(ibv) = 0.d0
        dz_total = 0.d0
        do k=1,n
          xkusa(ibv) = xkusa(ibv) + xkus(k,ibv)*dz(k)
          dz_total = dz_total + dz(k)
        enddo
        xkusa(ibv) = xkusa(ibv) / dz_total
      enddo
c     add gravitational potential to hl
      do k=1,n
        do ibv=i_bare,i_vege
          h(k,ibv)=h(k,ibv)+zc(k)
        end do
      end do
      return
      end subroutine hydra

      subroutine hl0
c**** hl0 sets up a table of theta values as a function of matric
c**** potential, h.  h is tabulated in a geometric series from
c**** 0 to hmin, with a first step of delh1.  the theta values
c**** depend not only on the matric potential, but also on the
c**** soil texture.  we solve a cubic equation to determine
c**** theta as a function of h.  hl0 also outputs the conductivity
c**** and diffusivity tables.
c**** input:
c**** a - matric potential function parameters
c**** b - hydraulic conductivity function parameters
c**** p - hydraulic diffusivity function parameters
c**** sat - saturated thetas
c**** output:
c**** thm(j,i) - theta at j''th h point for texture i
c**** xklm(j,i) - conductivity at j''th h point for texture i
c**** dlm(j,i) - diffusivity at j''th h point for texture i
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      integer, parameter :: nexp=6
      real*8, parameter :: c=2.3025851d0
      real*8, dimension(4,imt-1), parameter :: a=reshape(
     &     (/                   ! matric potential coefficients for
     &     .2514d0,  0.0136d0, -2.8319d0,  0.5958d0, ! sand
     &     .1481d0,  1.8726d0,  0.1025d0, -3.6416d0, ! loam
     &     .2484d0,  2.4842d0,  0.4583d0, -3.9470d0, ! clay
     &     .8781d0, -5.1816d0, 13.2385d0,-11.9501d0/), ! peat
     &     (/4,imt-1/))
      real*8, dimension(4,imt-1), parameter :: b=reshape(
     &     (/                   ! conductivity coefficients for
     &     -0.4910d0, -9.8945d0,  9.7976d0, -3.2211d0, ! sand
     &     -0.3238d0,-12.9013d0,  3.4247d0,  4.4929d0, ! loam
     &     -0.5187d0,-13.4246d0,  2.8899d0,  5.0642d0, ! clay
     &     -3.0848d0,  9.5497d0,-26.2868d0, 16.6930d0/), ! peat
     &     (/4,imt-1/))
      real*8, dimension(4,imt-1), parameter :: p=reshape(
     &     (/                   ! diffusivity coefficients for
     &     -0.1800d0, -7.9999d0,  5.5685d0, -1.8868d0, ! sand
     &     -0.1000d0,-10.0085d0,  3.6752d0,  1.2304d0, ! loam
     &     -0.1951d0, -9.7055d0,  2.7418d0,  2.0054d0, ! clay
     &     -2.1220d0,  5.9983d0,-16.9824d0,  8.7615d0/), ! peat
     &     (/4,imt-1/))
      real*8, parameter :: sat(imt-1) = (/.394d0,.537d0,.577d0,.885d0/)
      real*8 a1,a2,a3,alph0o,alpls1,arg(imt-1),delh1,delhn,dfunc,alph0
      real*8 diff,func,hmin,hs,s,sxtn,testh,xtol
      integer i,j,k,m,mmax

      sxtn=16.d0
      nth=2**nexp
      hlm(0)=0.0d0
      delh1=-0.00625d0
      hmin=-1000.d0
      delhn=delh1
c     solve for alph0 in s=((1+alph0)**n-1)/alph0
      s=hmin/delh1
      alph0=1.d0/8.d0
 10   alph0o=alph0
      alph0=(s*alph0+1.d0)**(1.d0/nth)-1.d0
      if(abs(alph0o-alph0).ge.1d-8) go to 10
      alpls1=1.0d0+alph0
      ! algdel=log(1.d0+alph0)  ! not used
      do 100 j=1,nth
        hlm(j)=hlm(j-1)+delhn
        delhn=alpls1*delhn
 100  continue
      mmax=100
      xtol=1d-6
      do 200 i=1,imt-1
        thm(0,i)=1.00d0
        do 150 j=1,nth
          hs=-exp(c*(a(1,i)+a(2,i)+a(3,i)+a(4,i)))
          a1=a(3,i)/a(4,i)
          a2=(a(2,i)-(log(-hlm(j)-hs))/c)/a(4,i)
          a3=a(1,i)/a(4,i)
          testh=thm(j-1,i)
          do 130 m=1,mmax
            func=(testh**3)+(a1*(testh**2))+(a2*(testh))+a3
            dfunc=(3*testh**2)+(2*a1*testh)+a2
            diff=func/dfunc
            testh=testh-diff
            if(abs(diff).lt.xtol) go to 140
 130      continue
          print *,'max # iterations:',mmax
 140      thm(j,i)=testh
 150    continue
 200  continue
      do 280 j=0,nth
        do 245 i=1,imt-1
          xklm(j,i)=0.d0
          arg(i)=0.d0
          do 240 k=-1,2
            arg(i)=arg(i)+b(k+2,i)*thm(j,i)**k
 240      continue
          arg(i)=min(arg(i),sxtn)
          arg(i)=max(arg(i),-sxtn)
          xklm(j,i)=exp(c*arg(i))
 245    continue
        !print *, 'xklm ', j, (xklm(j,i), i=1,imt-1)
        do 265 i=1,imt-1
          dlm(j,i)=0.d0
          arg(i)=0.d0
          do 260 k=-1,2
            arg(i)=arg(i)+p(k+2,i)*thm(j,i)**k
 260      continue
          arg(i)=min(arg(i),sxtn)
          arg(i)=max(arg(i),-sxtn)
          dlm(j,i)=exp(c*arg(i))
 265    continue
 280  continue
      do 350 j=0,nth
        do 310 i=1,imt-1
          thm(j,i)=thm(j,i)*sat(i)
 310    continue
 350  continue
      return
      end subroutine hl0


      subroutine evap_limits( compute_evap, evap_max_out, fr_sat_out )
!@sum computes maximal evaporation fluxes for current soil properties
!@calls cond
      use vegetation, only : veg_conductance
!@var compute_evap if .true. compute evap, else just evap_max,fr_sat
!@var evap_max_out max evaporation from unsaturated soil
!@var fr_sat_out fraction of saturated soil
      logical, intent(in) :: compute_evap
      real*8, intent(out) :: evap_max_out, fr_sat_out
ccc   local variables
!@var evap_max evap limits due to total amount of water in the soil
!@var evap_max_snow evap limits due to amount of snow
!@var evap_max_wet evap limit for saturated soil/canopy
!@var evap_max_dry evap limit for unsaturated soil/canopy
      real*8 :: evap_max(2), evap_max_snow(2)
      real*8 :: evap_max_wet(2), evap_max_dry(2)
!@var qb,qbs,qv,qvs specific humidity (b-bare, v-vege, s-snow)
!----------------------------------------------------------------------!
! adf
      real*8 :: qb, qbs, qv, qvs, qvg
!      real*8 :: qb, qbs, qvs
!----------------------------------------------------------------------!
!!!@var epb,epbs,epv,epvs potential evaporation (b-bare, v-vege, s-snow)
!      real*8 :: epbs, epvs, epvg ! , epb, epv
!@var rho3,cna,qm1dt local variable
      real*8 rho3, cna, qm1dt
      integer ibv, k
!@var hw the wilting point (m)
      real*8, parameter :: hw = -100.d0
!@var betadl transpiration efficiency for each soil layer
      real*8 betadl(ngm) ! used in evaps_limits only
      real*8 pot_evap_can
      real*8 cnc         ! local cnc from veg_conductance, nyk
#ifdef EVAP_VEG_GROUND
      real*8 evap_max_vegsoil
#endif
c     cna is the conductance of the atmosphere
      cna=ch*vsm
      rho3=rho/rhow ! i.e divide by rho_water to get flux in m/s

ccc make sure that important vars are initialized (needed for ibv hack)
      evap_max(:) = 0.d0
      evap_max_snow(:) = 0.d0
      evap_max_wet(:) = 0.d0
      evap_max_dry(:) = 0.d0
      betadl(:) = 0.d0
      betad = 0.d0
      abetad = 0.d0
      acna = 0.d0
      acnc = 0.d0

ccc !!! it''s a hack should call it somewhere else !!!
      !call hydra

      ! soil moisture
      do ibv=i_bare,i_vege
        evap_max(ibv) = 0.
        do k=1,n
          evap_max(ibv) = evap_max(ibv) +
     &         (w(k,ibv)-dz(k)*thetm(k,ibv))/dt
        enddo
      enddo

      ! maximal evaporation from the snow fraction
      ! may be too restrictive with resp. to pr, but will leave for now
      do ibv=i_bare,i_vege
        evap_max_snow(ibv) = pr
        do k=1,nsn(ibv)
          evap_max_snow(ibv) = evap_max_snow(ibv) +
     &         wsn(k,ibv)/dt
        enddo
      enddo

      ! evaporation from bare soil
      if ( process_bare ) then
        ibv = 1
        !!! no support for saturated soil yet, setting just in case...
        evap_max_wet(ibv) =  evap_max(ibv) + pr
        ! evap limited by diffusion and precipitation
        evap_max_dry(ibv) = min( evap_max(ibv),
     &       2.467d0*d(1,1)*(theta(1,1)-thetm(1,1))/dz(1) + pr )
      endif

      ! evaporation from the canopy
      if ( process_vege ) then
        ibv = 2
        evap_max_wet(ibv) = w(0,2)/dt !+ pr ! pr doesn''t work for snow
        ! soil under dry canopy
#ifdef EVAP_VEG_GROUND
        evap_max_vegsoil = min( evap_max(ibv),
     &       2.467d0*d(1,2)*(theta(1,2)-thetm(1,2))/dz(1) + pr )
#endif
        ! dry canopy
!!! this needs "qs" from the previous time step
c     betad is the the root beta for transpiration.
c     hw is the wilting point.
c     fr(k) is the fraction of roots in layer k
        betad=0.d0
        do 30 k=1,n
          betadl(k)=(1.d0-fice(k,2))*fr(k)*max((hw-h(k,2))/hw,zero)
          betad=betad+betadl(k)
 30     continue
        if ( betad < 1.d-12 ) betad = 0.d0 ! to avoid 0/0 divisions
        abetad=betad            ! return to old diagnostics
c     Get canopy conductivity cnc and gpp
        qv  = qsat(tp(0,2)+tfrz,lhe,pres) ! for cond_scheme==2
        call veg_conductance(
     &       cnc
     &       ,gpp
     &       ,trans_sw       !nyk
     &       ,betad          ! evaporation efficiency
     &       ,tp(0,2)          ! canopy temperature C
     &       ,qv
     &       ,dts
     &       )

        !print *,"HGY_COND: ",ijdebug, cnc, betadl
        !print *,"GHY_FORCINGS: ", ijdebug, tp(0,2)

!!!! test
 !!!       trans_sw = .1d0
 !!!       print *,'trans_sw = ', trans_sw

 !!!       trans_sw = .1d0
 !       trans_sw = min( trans_sw, .9d0 )

        betat=cnc/(cnc+cna+1d-12)
        abetat=betat            ! return to old diagnostics
        acna=cna                ! return to old diagnostics
        acnc=cnc                ! return to old diagnostics
!        agpp=gpp                !nyk 4/25/03.  Put in subroutine accm.
        pot_evap_can = betat*rho3*cna*(qsat(tp(0,2)+tfrz,lhe,pres) - qs)
        evap_max_dry(ibv) = 0.d0
        if ( betad > 0.d0 .and. pot_evap_can > 0.d0 ) then
          do k=1,n
            evap_max_dry(ibv) = evap_max_dry(ibv)
     &           +  min( pot_evap_can*betadl(k)/betad,
     &         (w(k,ibv)-dz(k)*thetm(k,ibv))/dt )
          enddo
        endif
!        evap_max_dry(ibv) = min( evap_max(ibv),
!     &       betat*rho3*cna*( qsat(tp(0,2)+tfrz,lhe,pres) - qs ) )
                   ! may be pr should be included somehow in e_m_d
      endif  !process_vege

      !! now we have to add the fluxes according to fractions
      ! bare soil
      evap_max_sat = fb*fr_snow(1)*evap_max_snow(1)
      evap_max_nsat = fb*(1.-fr_snow(1))*evap_max_dry(1)
      ! canopy
      evap_max_sat = evap_max_sat +
     &     fv*( fr_snow(2)*fm*evap_max_snow(2) +
     &          (1.-fr_snow(2)*fm)*theta(0,2)*evap_max_wet(2) )
      evap_max_nsat = evap_max_nsat +
     &     fv*( (1. - fr_snow(2)*fm) * ( 1. - theta(0,2) )
     &           * evap_max_dry(2) )

      fr_sat = fb * fr_snow(1) +
     &         fv * ( fr_snow(2)*fm + (1.-fr_snow(2)*fm)*theta(0,2) )

ccc set variables for output
      evap_max_out = evap_max_nsat
      fr_sat_out = fr_sat
ccc not sure if all this should be in one subroutine, but otherwise
ccc one has to pass all these flux limits

      if ( .not. ( evap_max_out > 0. .or.  evap_max_out <= 0. ) )then
        write(99,*)evap_max_out,cnc,evap_max(2)
        call stop_model("GHY::evap_limits: evap_max_out = NaN",255)
      endif

      if ( .not. compute_evap ) return

cccccccccccccccccccccccccccccccccccccccc

ccc the rest of evap_limits does not take into account process_bare,
ccc process_vege , but it should compute ok for dummy values.

c**** qm1 has mass of water vapor in first atmosphere layer, kg m-2
      qm1dt=.001d0*qm1/dt
! need this ?      if(igcm.ge.0 .and. igcm.le.3) xl=eddy/(z1-zs)

c     calculate bare soil, canopy and snow mixing ratios
      qb  = qsat(tp(1,1)+tfrz,lhe,pres)
      qbs = qsat(tsn1(1)+tfrz,lhe,pres)
      qv  = qsat(tp(0,2)+tfrz,lhe,pres)
      qvs = qsat(tsn1(2)+tfrz,lhe,pres)
#ifdef EVAP_VEG_GROUND
      qvg  = qsat(tp(1,2)+tfrz,lhe,pres)
#endif

c     potential evaporation for bare and vegetated soil and snow
      epb  = rho3*cna*(qb-qs)
      epbs = rho3*cna*(qbs-qs)
      epv  = rho3*cna*(qv-qs)
      epvs = rho3*cna*(qvs-qs)
#ifdef EVAP_VEG_GROUND
      epvg  = rho3*cna*(qvg-qs) ! actually not correct !
#endif

c     bare soil evaporation
      if ( process_bare ) then
        evapb = min( epb, evap_max_dry(1) )
        evapb = max( evapb, -qm1dt )
        evapbs = min( epbs, evap_max_snow(1) )
        evapbs = max( evapbs, -qm1dt )
!      evapor(1) = fr_snow(1)*evapbs + (1.-fr_snow(1))*evapb
      else
        evapb = 0.d0; evapbs = 0.d0;
      endif

c     vegetated soil evaporation
      if ( process_vege ) then
c     evapvd is dry evaporation (transpiration) from canopy
c     evapvw is wet evaporation from canopy (from interception)
        evapvg = 0.d0 ! in case EVAP_VEG_GROUND not defined
        evapvw = min( epv, evap_max_wet(2) )
        evapvw = max( evapvw,-qm1dt )
        evapvd = min(epv,evap_max_dry(2)) !evap_max_dry(2) depends on qs
        evapvd = max( evapvd, 0.d0 )
        evapvs = min( epvs, evap_max_snow(2) )
        evapvs = max( evapvs, -qm1dt )
#ifdef EVAP_VEG_GROUND
        evapvg = min(epvg,evap_max_vegsoil)
        ! keep evapv <= epv
        evapvg = min( evapvg, epv - evapvd*fd - evapvw*fw )
        evapvg = max( evapvg, 0.d0 )
#endif
!      evapor(2) = fr_snow(2)*fm*evapvs + (1.-fr_snow(2)*fm)*
!     &     ( theta(0,2)*evapvw + (1.-theta(0,2))*evapvd )
        if ( evapvw < 0.d0 ) then  ! we have dew
          fw = 1.d0 ; fd = 0.d0    ! let dew fall on entire canopy
        endif
      else
        evapvw = 0.d0; evapvd = 0.d0; evapvs = 0.d0
#ifdef EVAP_VEG_GROUND
        evapvg = 0.d0
#endif
      endif

      devapbs_dt = rho3*cna*qsat(tsn1(1)+tfrz,lhe,pres)
     &     *dqsatdt(tsn1(1)+tfrz,lhe)
      devapvs_dt = rho3*cna*qsat(tsn1(2)+tfrz,lhe,pres)
     &     *dqsatdt(tsn1(2)+tfrz,lhe)

c     compute transpiration for separate layers (=0 for bare soil)
      evapdl(1:n,1:2) = 0.d0
      if ( betad > 0.d0 ) then
        do k=1,n
          evapdl(k,2) = evapvd*betadl(k)/betad
        end do
      endif

      return
cccccccccccccccccccccccccccccccccccccccc
      end subroutine evap_limits


      subroutine sensible_heat
!@sum computes sensible heat for bare/vegetated soil and snow
      implicit none
      real*8 cna

      cna=ch*vsm
      snsh(1)=sha*rho*cna*(tp(1,1)-ts+tfrz)     ! bare soil
      snsh(2)=sha*rho*cna*(tp(0,2)-ts+tfrz)     ! canopy
      snshs(1) = sha*rho*cna*(tsn1(1)-ts+tfrz)  ! bare soil snow
      snshs(2) = sha*rho*cna*(tsn1(2)-ts+tfrz)  ! canopy snow
      dsnsh_dt = sha*rho*cna  ! derivative is the same for all above

      end subroutine sensible_heat


c      subroutine qsbal
c**** finds qs that balances fluxes.
c**** obtains qs by successive approximation.
c**** calculates evaporation.
c**** input:
c**** ch - heat conductivity coefficient from ground to surface
c**** vsm - surface layer wind speed, m s-1
c**** rho - air density, kg m-3
c**** eddy - transfer coefficient from surface to first atmosphere
c**** theta - water saturation of layers and canopy
c**** tp - temperatures of layers and canopy, c
c**** pres - atmospheric pressure at ground
c**** z1 - height of first layer, m
c**** zs - height of surface layer, m
c**** dz - layer thicknesses, m
c**** snowd - snow depths, equivalent water m
c**** pr - precipitation, m s-1
c**** q1 - mixing ratio of first layer
c**** fr - fraction of roots in layer
c**** fb - fraction of bare soil
c**** fv - fraction of vegetated soil
c**** hw - wilting point, m
c**** output:
c**** qs - mixing ratio at surface layer
c**** evap - evaporation from bare and vegetated regions, m s-1
c**** evapw - evaporation from wet canopy, m s-1, including from snow
c**** evapd - evaporation from dry canopy, m s-1
c**** evaps - evaporation from snow from canopy, m s-1
c**** betad - dry canopy beta, based on roots
c****


      subroutine drip_from_canopy
!@sum computes the flux of drip water (and its heat cont.) from canopy
!@+   the drip water is split into liquid water dripw() and snow drips()
!@+   for bare soil fraction the canopy is transparent
      use snow_drvm, only : snow_cover_same_as_rad
      real*8 ptmp,ptmps,pfac,pm,pmax
      real*8 snowf,snowfs,dr,drs
      integer ibv
      real*8 melt_w, melt_h
c     calculate snow fall.  snowf is snow fall, m s-1 of water depth.
      snowf=0.d0
      if(htpr.lt.0.d0)snowf=min(-htpr/fsn,pr)
c     snowfs is the large scale snow fall.
      snowfs=0.d0
      if(htprs.lt.0.d0)snowfs=min(-htprs/fsn,prs)
      if ( process_vege ) then
        ptmps=prs-snowfs
        ptmps=ptmps-evapvw*fw
        ptmp=pr-prs-(snowf-snowfs)
c     use effects of subgrid scale precipitation to calculate drip
        pm=1d-6
        pmax=fd0*pm
        drs=max(ptmps-pmax,zero)
        dr=drs
        if(ptmp.gt.0.d0)then
          pfac=(pmax-ptmps)*prfr/ptmp
          if(pfac.ge.0.d0)then
            if(pfac.lt.30.d0) dr=ptmp*exp(-pfac)
          else
            dr=ptmp+ptmps-pmax
          endif
        endif
ccc   make sure "dr" makes sense
        dr = min( dr, pr-snowf-evapvw*fw )
        dr = max( dr, pr-snowf-evapvw*fw - (ws(0,2)-w(0,2))/dts )
        dr = max( dr, 0.d0 )    ! just in case (probably don''t need it)
        dripw(2) = dr
        htdripw(2) = shw*dr*max(tp(0,2),0.d0) !don''t allow it to freeze
        ! snow falls through the canopy
        drips(2) = snowf
        htdrips(2) = min ( htpr, 0.d0 ) ! liquid H20 is 0 C, so no heat
      else  ! no vegetated fraction
        dripw(2) = 0.d0
        htdripw(2) = 0.d0
        drips(2) = 0.d0
        htdrips(2) = 0.d0
      endif
ccc   for bare soil drip is precipitation
      drips(1) = snowf
      htdrips(1) = min( htpr, 0.d0 )
      dripw(1) = pr - drips(1)
      htdripw(1) = htpr - htdrips(1)
      if ( snow_cover_same_as_rad .ne. 0 ) then
#define TRY_TO_MELT_FRESH_SNOW_ON_WARM_GROUND
#ifdef TRY_TO_MELT_FRESH_SNOW_ON_WARM_GROUND
        do ibv=1,2
          if ( tp(1,ibv) > 0.d0 ) then
            melt_w = (1.d0-fr_snow(ibv)) * drips(ibv)
            melt_h = (1.d0-fr_snow(ibv)) * htdrips(ibv)
            drips(ibv) = drips(ibv) - melt_w
            htdrips(ibv) = htdrips(ibv) - melt_h
            dripw(ibv) = dripw(ibv) + melt_w
            htdripw(ibv) = htdripw(ibv) + melt_h
          endif
        enddo
#endif
      endif
      return
      end subroutine drip_from_canopy


      subroutine flg
!@sum computes the water fluxes at the surface

c     bare soil
      if ( process_bare ) then
        f(1,1) = -flmlt(1)*fr_snow(1) - flmlt_scale(1)
     &       - (dripw(1)-evapb)*(1.d0-fr_snow(1))
      endif
      if ( process_vege ) then
c     upward flux from wet canopy
        fc(0) = -pr+evapvw*fw*(1.d0-fm*fr_snow(2))
                                ! snow masking of pr is ignored since
                                ! it is not included into drip
        fc(1) = -dripw(2) - drips(2)
c     vegetated soil
!!! flux down from canopy is not -f(1,2) any more !!
!!! it is =  -dripw(2) - drips(2)
!!! f(1,2) is a flux up from tyhe soil
        f(1,2) = -flmlt(2)*fr_snow(2) - flmlt_scale(2)
     &       - (dripw(2)-evapvg)*(1.d0-fr_snow(2))
      endif

c     compute evap_tot for accumulators
      evap_tot(1) = evapb*(1.d0-fr_snow(1)) + evapbs*fr_snow(1)
      evap_tot(2) = (evapvw*fw + evapvd*fd)*(1.d0-fr_snow(2)*fm)
     &     + evapvs*fr_snow(2)*fm + evapvg*(1.d0-fr_snow(2))

      return
      end subroutine flg


      subroutine flhg
!@sum calculates the ground heat fluxes (to the surface)
c**** input:
c**** output:
c**** fh - heat fluxes from bare soil surface, and from canopy,
c****      and between canopy and vegetated soil.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c**** bare soil fluxes
      real*8 thrm_can, thrm_soil(2)
      thrm_can = stbo * (tp(0,2)+tfrz)**4
      thrm_soil(1:2) = stbo * (tp(1,1:2)+tfrz)**4

!!!! test
 !!!     thrm_can = 0.d0
      ! bare soil
      if ( process_bare ) then
        fh(1,1) = - fhsng(1)*fr_snow(1) - fhsng_scale(1)
     &       + ( -htdripw(1) +
     &       evapb*elh + snsh(1) + thrm_soil(1) - srht - trht)
     &       *(1.d0-fr_snow(1))
      endif
      ! vegetated soil
      if ( process_vege ) then
        fh(1,2) = - fhsng(2)*fr_snow(2) - fhsng_scale(2)
     &       + ( -htdripw(2)
     &       + evapvg*elh        + thrm_soil(2) - thrm_can
#ifdef RAD_VEG_GROUND
     &       * (1.d0-trans_sw)
     &       - trans_sw*(srht + trht)
#endif
     &       )
     &       *(1.d0-fr_snow(2))
        ! canopy
        fch(0) = -htpr +
     &       (evapvw*elh*fw + snsh(2) + thrm_can
#ifdef RAD_VEG_GROUND
     &       * (1.d0-trans_sw)
     &       - (1.d0 - trans_sw)*(srht + trht)
#else
     &       -srht-trht
#endif
     &       + evapvd*elh*fd)
     &       *(1.d0-fm*fr_snow(2))
        fch(1) = -(thrm_can - thrm_soil(2))
#ifdef RAD_VEG_GROUND
     &       *(1.d0 - trans_sw)
#endif
     &       *(1.d0-fr_snow(2)) !rad soil
     &       - (thrm_can - thrmsn(2))*fr_snow(2)*(1.d0-fm)    !rad snow
#ifdef RAD_VEG_GROUND
     &       *(1.d0 - trans_sw)
#endif
     &       - htdripw(2) - htdrips(2)                  !heat of precip
      endif

      ! compute thrm_tot for accumulators
      thrm_tot(1) = thrm_soil(1)*(1.d0-fr_snow(1))
     &     + thrmsn(1)*fr_snow(1)
      thrm_tot(2) = thrm_can*(1.d0-fr_snow(2)*fm)
#ifdef RAD_VEG_GROUND
     &     *(1.d0-trans_sw)
#endif
     &     + thrmsn(2)*fr_snow(2)*fm
#ifdef RAD_VEG_GROUND
     &     + thrmsn(2)*trans_sw*fr_snow(2)*(1.d0-fm)
     &     + trans_sw*thrm_soil(2)*(1.d0-fr_snow(2))
#endif

!XXXXXXXXXXXXXXXX don't forget to add trans_sw stuff to snow !

c     compute total sensible heat
      snsh_tot(1) = snsh(1)*(1.d0-fr_snow(1)) + snshs(1)*fr_snow(1)
      snsh_tot(2) = snsh(2)*(1.d0-fr_snow(2)*fm)
     &     + snshs(2)*fr_snow(2)*fm
      return
      end subroutine flhg


      subroutine runoff
c**** calculates surface and underground runoffs.
c**** input:
c**** xinfc - infiltration capacity, m s-1
c**** prfr - fraction of precipitation
c**** xk - conductivity, m s-1
c**** dz - layer thicknesses, m
c**** sl - slope
c**** sdstnc - interstream distance, m
c**** output:
c**** rnf - surface runoff
c**** rnff - underground runoff, m s-1
c     use effects of subgrid scale rain
c     use precipitation that includes smow melt
ccc surface runoff was rewritten in a more clear way 7/30/02
      real*8 f_k0_exp_k !@var f_k0_exp_k coefficient for the topmodel
!@var water_down flux of water at the soil surface
!@var satfrac fraction of saturated soil
!@var prec_fr soil fraction at which precipitation is falling
      real*8 water_down, satfrac, prec_fr
      integer ibv,k
!@var sdstnc interstream distance (m)
      real*8, parameter :: sdstnc = 100.d0
!@var rosmp used to compute saturated fraction: (w/ws)**rosmp
      real*8, parameter :: rosmp = 8.
      rnff(:,:) = 0.d0
      rnf(:) = pr  ! hack to conserve water (for ibv != 0,1)
                   ! - should be set to 0 after testing
c**** surface runoff
      do ibv=i_bare,i_vege
        water_down = -f(1,ibv)
        water_down = max( water_down, zero ) ! to make sure rnf > 0
        ! everything that falls on saturated fraction goes to runoff
        satfrac = (w(1,ibv)/ws(1,ibv))**rosmp
        rnf(ibv) = satfrac * water_down
        water_down = (1.d0 - satfrac) * water_down
        ! if we introduce large scale precipitation it should be
        ! applied here
        !!! the following line is a hack. in a more precise approach
        !!! one should treat snow-free fraction separately
        prec_fr = max( prfr, fr_snow(ibv) )
        if ( water_down*30.d0 > xinfc(ibv)*prec_fr ) then
          rnf(ibv) = rnf(ibv) +
     &         water_down * exp( -xinfc(ibv)*prec_fr/water_down )
        endif
      enddo

c**** underground runoff
c     sl is the slope, sdstnc is the interstream distance
      do ibv=i_bare,i_vege
ccc this is some rough estimate for the expression f * k0 exp(zbar/f)
ccc in the topmodel expression for the runoff
!        f_k0_exp = 0.d0
!        do k=1,n
!          f_k0_exp = f_k0_exp + xkus(k,ibv)*w(k,ibv)/ws(k,ibv)*dz(k)
!        enddo
        do k=1,n
          rnff(k,ibv)=xku(k,ibv)*sl*dz(k)/sdstnc
!/* #define do_topmodel_runoff */
#ifdef do_topmodel_runoff
          if ( ws(k,ibv) > 1.d-16 ) then
            f_k0_exp_k = (1.d0-fice(k,ibv))
     $           * xkus(k,ibv)*w(k,ibv)/ws(k,ibv)*dz(k)
          else
            f_k0_exp_k = 0.d0
          endif
c         print *,'rnff: ', k, rnff(k,ibv), f_k0_exp_k*exp( -top_index )
c         print *, xkus(k,ibv),w(k,ibv),ws(k,ibv),dz(k),top_index
          rnff(k,ibv)=f_k0_exp_k*exp( -top_index )
#endif
        end do
!        print *,'runoff: ', sl/sdstnc, exp( -top_index )
      end do
      return
      end subroutine runoff


      subroutine fllmt
c**** places limits on the soil water fluxes
c**** input:
c**** w - water in layers, m
c**** ws - saturated water in layers, m
c**** dts - current time step size, s
c**** f - water fluxes, m s-1
c**** snk - water sink from layers, m s-1
c**** rnf - surface runoff, m s-1
c**** snowd - snow depth, equivalent water m
c**** snowf - snow fall, equivalent water m s-1
c**** output:
c**** f - limited water fluxes, m s-1
c**** snk - limited water sinks, m s-1
c**** rnf - limited surface runoff, m s-1
c**** temp variables:
c**** snowdu - the upper bound on the snow depth at end of time step
c**** snowdl - the lower bound on the snow depth at end of time step
c**** trunc - fix for truncation on ibm mainframes
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      real*8 dflux,drnf,trunc
      integer k, ibv
      real*8 wn
      trunc=1d-6
      trunc=1d-12
      trunc=0.d0 ! works better since at some places thetm=thets=0
ccc   prevent over/undersaturation of layers 2-n
      do ibv=i_bare,i_vege
        do k=n,2,-1
          wn = w(k,ibv) + ( f(k+1,ibv) - f(k,ibv)
     &         - rnff(k,ibv) - fd*(1.-fr_snow(2)*fm)*evapdl(k,ibv) )*dts
ccc   compensate oversaturation by increasing flux up
          if (wn - ws(k,ibv) > trunc)
     &         f(k,ibv) = f(k,ibv) + (wn - ws(k,ibv) + trunc)/dts
ccc   compensate undersaturation by decreasing runoff
          if (wn - dz(k)*thetm(k,ibv) < trunc) then
            rnff(k,ibv) = rnff(k,ibv) +
     &         (wn - dz(k)*thetm(k,ibv) - trunc)/dts
            if ( rnff(k,ibv) < 0.d0 ) then ! have to compensate with f
              f(k,ibv) = f(k,ibv) + rnff(k,ibv)
              rnff(k,ibv) = 0.d0
            endif
          endif
        end do
      end do
ccc   prevent over/undersaturation of first layer
      do ibv=i_bare,i_vege
        wn = w(1,ibv) + ( f(2,ibv) - f(1,ibv)
     &       - rnf(ibv) - rnff(1,ibv)
     &       - fd*(1.-fr_snow(2)*fm)*evapdl(1,ibv) )*dts
        if ( wn - ws(1,ibv) > trunc)
     &       rnf(ibv) = rnf(ibv) + (wn - ws(1,ibv) + trunc)/dts
        if ( wn - dz(1)*thetm(1,ibv) < trunc )
     &       rnf(ibv) = rnf(ibv) +
     &       (wn - dz(1)*thetm(1,ibv) - trunc)/dts
      enddo
ccc   now trying to remove negative runoff
      do ibv=i_bare,i_vege
        k = 1
        do while ( rnf(ibv) .lt. 0.d0 .and. k .le. n )
          if ( k > 1 ) then
ccc         this is how much water we can take from layer k
            dflux = f(k+1,ibv) + (w(k,ibv)-dz(k)*thetm(k,ibv))/dts
     &           - f(k,ibv) - rnff(k,ibv)
     &           - fd*(1.-fr_snow(2)*fm)*evapdl(k,ibv)
            f(k,ibv) = f(k,ibv) - rnf(ibv)
            rnf(ibv) = rnf(ibv) + min( -rnf(ibv), dflux )
          endif
ccc    rnff always >= 0, use it also to compensate rnf<0
          if ( rnff(k,ibv) < 0.d0 )
     &         call stop_model('fllmt: negative underground runoff',255)
          drnf = min( -rnf(ibv), rnff(k,ibv) )
          rnf(ibv) = rnf(ibv) + drnf
          rnff(k,ibv) = rnff(k,ibv) - drnf
          k = k + 1
        enddo
ccc    check if rnf==0 up to machine accuracy
        if ( rnf(ibv) .lt. -1d-12 ) then
          print *, 'fllmt: rnf<0, ibv=',ibv,rnf(ibv)
            call stop_model('fllmt: negative runoff',255)
        endif
ccc    if -1d-12 < rnf < 0. put it to 0 to avoid possible problems
ccc    actually for ground hydrology it is not necessary
        rnf(ibv) = max ( rnf(ibv), 0.d0 )
      enddo

!!! hack to prevent taking water from empty layers
!!! this is probably not quite correct so it is disabled by default
cddd      do ibv=i_bare,i_vege
cddd        do k=n,1,-1
cddd          if ( f(k,ibv) > 0.d0 .and. w(k,ibv) < 1.d-16 ) then
cddd            print *,"GHY: corrected flux->0 ",k,ibv,f(k,ibv)
cddd            f(k,ibv) = 0.d0
cddd          endif
cddd        enddo
cddd      enddo
      return
      end subroutine fllmt


      subroutine xklh( xklh0_flag )
c**** evaluates the heat conductivity between layers
c**** uses the method of DeVries.
c**** input:
c**** zb - soil layer boundaries, m
c**** zc - soil layer centers, m
c**** theta - soil water saturation
c**** fice - fraction of ice in layers
c**** alami - ice heat conductivity
c**** alamw - water heat conductivity
c**** tp - temperature of layers, c
c**** shw - specific heat of water
c**** shi - specific heat of ice
c**** shc - heat capacity of soil layers
c**** dz - layer thicknesses
c**** k,ibv - soil layer
c**** dts - the current time step
c**** output:
c**** xkh(k,ibv) - heat conductivities in each of the soil layers
c**** xkhm(k,ibv) - average heat conductivity between layer k and k-1
c****
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c     dimension xsha(ng,2),xsh(ng,2),gabc(3),hcwt(imt-1)
c
c     calculate with changing ga for air. ga is the depolarization
c     factor for air, calculated by linear interpolation from .333d0
c     at saturation to .035 at 0 water, following devries.
      integer, intent(in), optional :: xklh0_flag
      real*8 gaa,gabc(3),gca,xa,xb,xden,xi,xnum,xs,xw
ccc   real*8, save :: ba,hcwt(imt-1),hcwta,hcwtb,hcwti,hcwtw
ccc   real*8, save :: xsha(ng,2),xsh(ng,2)
      real*8  :: ba,hcwt(imt-1),hcwta,hcwtb,hcwti,hcwtw
      real*8  :: xsha(ng,2),xsh(ng,2)
      COMMON /XKLHSAV/ BA, HCWTW, HCWTI, HCWTB, XSHA, XSH
!$OMP  THREADPRIVATE (/XKLHSAV/)
      integer i, j, ibv, k
c the alam''s are the heat conductivities
      real*8, parameter :: alamw = .573345d0
     &     ,alami = 2.1762d0
     &     ,alama = .025d0
     &     ,alambr= 2.9d0
     &     ,alams(imt-1) = (/ 8.8d0, 2.9d0, 2.9d0, .25d0 /)
      if ( present(xklh0_flag) ) goto 777
      do ibv=i_bare,i_vege
        do k=1,n
          gaa=.298d0*theta(k,ibv)/(thets(k,ibv)+1d-6)+.035d0
          gca=1.d0-2.d0*gaa
          hcwta=(2.d0/(1.d0+ba*gaa)+1.d0/(1.d0+ba*gca))/3.d0
c     xw,xi,xa are the volume fractions. don''t count snow in soil lyr 1
          xw=w(k,ibv)*(1.d0-fice(k,ibv))/dz(k)
          xi=w(k,ibv)*fice(k,ibv)/dz(k)
          xa=(thets(k,ibv)-theta(k,ibv))
          xb=q(imt,k)
          xnum=xw*hcwtw*alamw+xi*hcwti*alami+xa*hcwta*alama+xsha(k,ibv)
     &         + xb*hcwtb*alambr
          xden=xw*hcwtw+xi*hcwti+xa*hcwta+xsh(k,ibv)+xb*hcwtb
          xkh(k,ibv)=xnum/xden
          if ( xkh(k,ibv) .lt. 0.d0 )
     &         call stop_model('xklh: heat conductivity<0',255)
        end do
      end do
c     get the average conductivity between layers
      do ibv=i_bare,i_vege
        do k=2,n
          xkhm(k,ibv)=((zb(k)-zc(k-1))*xkh(k,ibv)
     &         + (zc(k)-zb(k))*xkh(k-1,ibv)
     &         )/(zc(k)-zc(k-1))
        end do
      end do
c****
      return
!      entry xklh0
 777  continue
c gabc''s are the depolarization factors, or relative spheroidal axes.
      gabc(1)=.125d0
      gabc(2)=gabc(1)
      gabc(3)=1.d0-gabc(1)-gabc(2)
c hcwt''s are the heat conductivity weighting factors
      hcwtw=1.d0
      hcwti=0.d0
      hcwtb=1.d0
      do i=1,imt-1
      hcwt(i)=0.d0
      end do
      do j=1,3
        hcwti=hcwti+1.d0/(1.d0+(alami/alamw-1.d0)*gabc(j))
        do i=1,imt-1
          hcwt(i)=hcwt(i)+1.d0/(1.d0+(alams(i)/alamw-1.d0)*gabc(j))
        end do
      end do
      hcwti=hcwti/3.d0
      do i=1,imt-1
        hcwt(i)=hcwt(i)/3.d0
      end do
      do ibv=1,2        ! i_bare,i_vege
        do k=1,n
          xsha(k,ibv)=0.d0
          xsh(k,ibv)=0.d0
          do i=1,imt-1
            xs=(1.d0-thm(0,i))*q(i,k)
            xsha(k,ibv)=xsha(k,ibv)+xs*hcwt(i)*alams(i)
            xsh(k,ibv)=xsh(k,ibv)+xs*hcwt(i)
          end do
        end do
      end do
      ba=alama/alamw-1.d0
      return
      end subroutine xklh


      subroutine fl
!@sum computes water flux between the layers
c**** input:
c**** h - soil potential of layers, m
c**** xk - conductivity of layers, m s-1
c**** zc - layer centers, m
c**** output:
c**** f - fluxes between layers, m s-1
c**** xinfc - infiltration capacity, m s-1
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c****
      integer ibv,k
      do ibv=i_bare,i_vege
        f(n+1,ibv)=0.d0
      end do
c****
      do ibv=i_bare,i_vege
        do k=2,n
          f(k,ibv)=-xk(k,ibv)*(h(k-1,ibv)-h(k,ibv))/(zc(k-1)-zc(k))
        end do
      end do
c**** put infiltration maximum into xinfc
      do ibv=i_bare,i_vege
        xinfc(ibv)=xk(1,ibv)*h(1,ibv)/zc(1)
      end do
      return
      end subroutine fl


      subroutine flh
c**** evaluates the heat flux between layers
c**** subroutine fl must be called first
c**** input:
c**** zb - soil layer boundaries, m
c**** zc - soil layer centers, m
c**** theta - soil water saturation
c**** fice - fraction of ice in layers
c**** alami - ice heat conductivity
c**** alamw - water heat conductivity
c**** tp - temperature of layers, c
c**** shw - specific heat of water
c**** output:
c**** fh - heat flux between layers
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c****
      integer ibv, k
      do ibv=i_bare,i_vege
        fh(n+1,ibv)=0.d0
c total heat flux is heat carried by water flow plus heat conduction
        do k=2,n
          fh(k,ibv)=-xkhm(k,ibv)*(tp(k-1,ibv)-tp(k,ibv))/(zc(k-1)-zc(k))
          if(f(k,ibv).gt.0)then
            fh(k,ibv)=fh(k,ibv)+f(k,ibv)*tp(k,ibv)*shw
          else
            fh(k,ibv)=fh(k,ibv)+f(k,ibv)*tp(k-1,ibv)*shw
          endif
        end do
      end do
      return
      end subroutine flh


      subroutine retp
c**** evaluates the temperatures in the soil layers based on the
c**** heat values.  also executes snow melt.
c**** input:
c**** w - water in soil layers, m
c**** ht - heat in soil layers
c**** fsn - heat of fusion of water
c**** shc - specific heat capacity of soil
c**** shi - specific heat capacity of ice
c**** shw - specific heat capcity of water
c**** snowd - snow depth, equivalent water m
c**** fb - fraction of bare soil
c**** fv - fraction of vegetation
c**** fm - snow vegetation masking fraction (requires reth called first)
c**** output:
c**** tp - temperature of layers, c
c**** fice - fraction of ice of layers
ccc   include 'soils45.com'http://www.giss.nasa.gov/internal/
c**** soils28   common block     9/25/90
      integer ibv, k, kk
      tp(:,:) = 0.d0
      fice(:,:) = 0.d0
      do ibv=i_bare,i_vege
        kk=2-ibv
        do k=kk,n
          ! tp(k,ibv)=0.d0
          !if(w(k,ibv).ge.1d-12)then
          !  fice(k,ibv)=-ht(k,ibv)/(fsn*w(k,ibv))
          !else
          !  fice(k,ibv)=0.d0
          !endif
          if( fsn*w(k,ibv)+ht(k,ibv) .lt. 0.d0 ) then ! all frozen
            tp(k,ibv)=(ht(k,ibv)+w(k,ibv)*fsn)/(shc(k,ibv)+w(k,ibv)*shi)
            fice(k,ibv)=1.d0
          else if( ht(k,ibv) .gt. 0.d0 ) then ! all melted
            tp(k,ibv)=ht(k,ibv)/(shc(k,ibv)+w(k,ibv)*shw)
            !shc -- specific heat of canopy
            !shw -- specific heat of water
            ! fice(k,ibv)=0.d0
          else if( w(k,ibv) .ge. 1d-12 )then  ! part frozen
            fice(k,ibv)=-ht(k,ibv)/(fsn*w(k,ibv))
          endif
        end do
      end do
ccc this is a fix for undefined tsn1 at the beginning of soil routines
ccc probably should be moved to some other place
      tsn1(:) = 0.d0
      do ibv=i_bare,i_vege
         ! tsn1(ibv) = 0.d0
         if (  wsn(1,ibv) .gt. 1.d-6 .and.
     &         hsn(1,ibv) + wsn(1,ibv)*fsn .lt. 0.d0  ) then
            tsn1(ibv) = (hsn(1,ibv) + wsn(1,ibv)*fsn)/(wsn(1,ibv)*shi)
         endif
ccc the following is a hack. it is necessary only at the beginning of th
ccc run, when some temperatures are not initialized properly.
ccc should be removed when program is rewritten in a more clean way...
!!! I think it is not needed any more ...
!         if ( wsn(1,ibv) .le. 1.d-6 ) then
!            tsn1(ibv) = tp(2-ibv,ibv)
!         endif
      enddo

      if(tp(1,1).gt.120.d0.or.tp(0,2).gt.120.d0
     &     .or. tp(1,1)<-150.d0 .or. tp(0,2)<-150.d0 )then
        write(6,*)'retp tp bounds error'
        write(6,*)'ijdebug',ijdebug
        call reth
        call hydra
        call outw(1)
        call stop_model(
     &       'retp: tground > 120C - see soil_outw and fort.99',255)
      endif
      return
      end subroutine retp


      subroutine apply_fluxes
!@sum apply computed fluxes to adwance w and h
!@ver 1.0
c**** shw - specific heat of water
c**** tp - temperature of layers, c
      integer ibv, k
ccc   the canopy
      w(0,2) = w(0,2) + ( fc(1) - fc(0) )*dts
      ht(0,2)=ht(0,2) + ( fch(1) - fch(0) )*dts
ccc   the soil
      do ibv=i_bare,i_vege
ccc     surface runoff
        w(1,ibv) = w(1,ibv) - rnf(ibv)*dts
        ht(1,ibv) = ht(1,ibv) - shw*max(tp(1,ibv),0.d0)*rnf(ibv)*dts
ccc     rest of the fluxes
        do k=1,n
          w(k,ibv) = w(k,ibv) +
     &         ( f(k+1,ibv) - f(k,ibv) - rnff(k,ibv)
     &         - fd*(1.-fr_snow(2)*fm)*evapdl(k,ibv)
     &         )*dts
          ht(k,ibv) = ht(k,ibv) +
     &         ( fh(k+1,ibv) - fh(k,ibv) -
     &         shw*max(tp(k,ibv),0.d0)*( rnff(k,ibv)
!     &         + fd*(1.-fr_snow(2)*fm)*evapdl(k,ibv)   !!! hack !!!
     &         ) )*dts
        end do
      end do

ccc   do we need this check ?
      if ( w(0,2) < 0.d0 ) then
        if (w(0,2)<-1.d-12)write(0,*)'GHY:CanopyH2O<0 at',ijdebug,w(0,2)
        w(0,2) = 0.d0
      endif

ccc check for under/over-saturation
      do ibv=i_bare,i_vege
        do k=1,n
          if ( w(k,ibv) < dz(k)*thetm(k,ibv) - 1.d-14 ) then
            print*,"ghy:",k,ibv,w(k,ibv),dz(k),thetm(k,ibv)
            call stop_model("ghy: w < dz*thetm",255)
          end if
          if ( w(k,ibv) > ws(k,ibv) + 1.d-14 )
     &         call stop_model("ghy: w > ws",255)
          w(k,ibv) = max( w(k,ibv), dz(k)*thetm(k,ibv) )
          w(k,ibv) = min( w(k,ibv), ws(k,ibv) )
        enddo
      enddo

      return
      end subroutine apply_fluxes


      subroutine advnc
c**** advances quantities by one time step.
c**** input:
c**** dt - time step, s
c**** dz - layer thickness, m
c**** tp - layer temperatures, c
c**** tfrz - freezing point of water, k
c**** w - soil water in layers, m
c**** snowd - snow depth, m
c**** f - water flux, m s-1
c**** snk - water sinks, m s-1
c**** ht - heat in soil layers
c**** fh - heat flux in soil layers
c**** snkh - heat sink in layers
c**** snowf - snow fall, m s-1 of equivalent water
c**** evap - evaporation, m s-1
c**** output:
c**** w - updater water in soil layers, m s-1
c**** ht - updated heat in soil layers
c**** snowd - updated snow depth, m s-1 of equivalent water
c**** rus - overall surface runoff, m s-1   replaced by aruns
c**** aruns - overall surface runoff, kg m-2
c**** aeruns - overall surface heat runoff, j m-2
c**** aerunu - underground heat runoff, j m-2
c**** uses:
c**** retp,reth,fl,flg,runoff,sink,sinkh,fllmt,flh,flhg.
c**** also uses surf with its required variables.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      use vegetation, only: update_veg_locals

      real*8 dtm,tb0,tc0,dtr,tot_w1
      integer limit,nit
      real*8 dum1, dum2
      limit=300   ! 200 increase to avoid a few more stops
      nit=0
      dtr=dt
      dts=dt ! make sure that dts is always initialized
ccc trying to skip fb==0 and fv==0 fractions of the cell
ccc reset main water/heat fluxes, so they are always initialized
      f(:,:) = 0.d0
      fh(:,:) = 0.d0
      fc(:) = 0.d0
      fch(:) = 0.d0
ccc normal case (both present)
      i_bare = 1; i_vege = 2
      process_bare = .true.; process_vege = .true.
      if ( fb == 0.d0 ) then  ! bare fraction is missing
        i_bare = 2
        process_bare = .false.
      endif
      if ( fv == 0.d0 ) then  ! bare fraction is missing
        i_vege = 1
        process_vege = .false.
      endif

!debug debug!
!      pr = 0.d0
!      htpr = 0.d0
!      trpr = 0.d0
!      srht = 0.d0
!      trht = 0.d0
!!!
      call reth
      call retp
      tb0=tp(1,1)
      tc0=tp(0,2)
cddd      print '(a,10(e12.4))', 'ghy_temp_b ',
cddd     &     tp(1,1),tp(2,1),tp(0,2),tp(1,2),tp(2,2)
ccc accm0 was not called here in older version - check
      call accm(0)
      do while ( dtr > 0.d0 )
        nit=nit+1
        if(nit.gt.limit)go to 900
        call hydra
        !call qsbal
!debug
!        fm = 1.d0
!!!
        call evap_limits( .true., dum1, dum2 )
        call sensible_heat
!debug debug!
!        evapb = 0.d0
!        evapbs = evapbs * .1
!        evapvw = 0.d0
!        evapvd = 0.d0
!        evapdl = 0.d0
 !       evapvs = 0.d0
!        devapbs_dt = 0.d0
!        devapvs_dt =0.d0
!        dsnsh_dt = 0.d0
!         if ( evapvs < 0.d0 ) then
!           evapvs = 0.d0; devapvs_dt =0.d0
!         endif
!!!

        call xklh
        call gdtm(dtm)
        !print *,'dtm ', ijdebug, dtm
        if ( dtm >= dtr ) then
          dts = dtr
          dtr = 0.d0
        else
          dts = min( dtm, dtr*0.5d0 )
          dtr=dtr-dts
        endif
!!!
 !       drips = 0
 !       dripw = 0
        call drip_from_canopy
        call check_water(0)
        call check_energy(0)
        call snow
        call fl
        call flg
         ! call check_f11
        call runoff
!debug
!        rnff(:,:) = 0.d0
 !       rnf(:) = 0.d0
!!!
        !call sink
        call fllmt
!debug
!        rnff(:,:) = 0.d0
!        rnf(:) = 0.d0
!!!
         ! call check_f11
        !call sinkh
        call flh
        call flhg
c     call fhlmt
         ! call check_f11
        !call check_water(0)
#ifdef TRACERS_WATER
        call ghy_tracers
#endif
        call apply_fluxes

cddd      print '(a,10(e12.4))', 'tr_w_1 '
cddd     &     , tr_w(1,:,1) - w(:,1) * 1000.d0
cddd      print '(a,10(e12.4))', 'tr_w_2 '
cddd     &     , tr_w(1,:,2) - w(:,2) * 1000.d0
!!!
        ! if(wsn_max>0) call restrict_snow (wsn_max)
        call check_water(1)
        call check_energy(1)
        call accm
        call reth
        call retp
cddd      print '(a,i6,10(e12.4))', 'ghy_temp ', ijdebug,
cddd     &     tp(1,1),tp(2,1),tp(0,2),tp(1,2),tp(2,2)

        call update_veg_locals(evap_tot(2), rho, rhow, ch, vsm,qs)

#ifdef TRACERS_WATER
C**** finalise surface tracer concentration here
        tot_w1 = fb*( w(1,1)*(1.d0-fr_snow(1))
     &       + wsn(1,1)*fr_snow(1) )
     &       + fv*( w(0,2)*(1.d0-fm*fr_snow(2))
     &       + wsn(1,2)*fm*fr_snow(2) )
     &       + fv*w(1,2)*0.01 ! hack !!! (for 0 H2O in canopy)
        if ( tot_w1 > 1.d-30 ) then
          atr_g(:ntg) =         ! instantaneous
     &         ( fb*( tr_w(:ntg,1,1)*(1.d0-fr_snow(1))
     &         + tr_wsn(:ntg,1,1) ) !*fr_snow(1)
     &         + fv*( tr_w(:ntg,0,2)*(1.d0-fm*fr_snow(2))
     &         + tr_wsn(:ntg,1,2)*fm )
     &         + fv*tr_w(:ntg,1,2)*0.01 ! hack !!! (for 0 H2O in canopy)
     &         ) /              !*fr_snow(2)
     &         (rhow * tot_w1)  ! * dts
        endif
#endif

      enddo

      call accm(1)
      call hydra
      call wtab  ! for gcm diag. only
      return
  900 continue
      write(99,*)'limit exceeded'
      write(99,*)'dtr,dtm,dts',dtr,dtm,dts
      write(99,*)'tb0,tc0',tb0,tc0
      call hydra
      call outw(2)
      call stop_model(
     &     'advnc: time step too short - see soil_outw and fort.99',255)
      end subroutine advnc


      subroutine accm( flag )
c**** accumulates gcm diagnostics
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c**** the following lines were originally called before retp,
c**** reth, and hydra.
      integer, intent(in), optional :: flag
      real*8 qsats
      real*8 cpfac,dedifs,dqdt,el0,epen,h0
      integer k
#ifdef TRACERS_WATER
      real*8 tot_w1
#endif

      if ( present(flag) ) then
        if ( flag == 0 ) goto 778
                         goto 777
      endif

ccc   main fluxes which should conserve water/energy
      atrg = atrg + ( thrm_tot(1)*fb + thrm_tot(2)*fv )*dts
      ashg = ashg + ( snsh_tot(1)*fb + snsh_tot(2)*fv )*dts
      aevap = aevap + ( evap_tot(1)*fb + evap_tot(2)*fv )*dts
      alhg = elh*aevap
      aruns=aruns+(fb*rnf(1)+fv*rnf(2))*dts
      aeruns=aeruns+shw*( fb*max(tp(1,1),0.d0)*rnf(1)
     &                  + fv*max(tp(1,2),0.d0)*rnf(2) )*dts
      do k=1,n
        arunu=arunu+(rnff(k,1)*fb+rnff(k,2)*fv)*dts
        aerunu=aerunu+ shw*( max(tp(k,1),0.d0)*rnff(k,1)*fb
     *                     + max(tp(k,2),0.d0)*rnff(k,2)*fv )*dts
      end do
ccc   end of main fluxes
ccc   the rest of the fluxes (mostly for diagnostics)
      aflmlt = aflmlt + ( fb*(flmlt(1)*fr_snow(1)+flmlt_scale(1))
     $     + fv*(flmlt(2)*fr_snow(2)+flmlt_scale(2)) )*dts
!      if(process_vege) aintercep = aintercep + pr - dripw(2) - drips(2)
ccc   max in the following expression removes extra drip because of dew
      aintercep = aintercep + fv*max(pr - dripw(2) - drips(2),0.d0)*dts
      aevapw = aevapw + evapvw*fw*(1.d0-fr_snow(2)*fm)*fv*dts
      aevapd = aevapd + evapvd*fd*(1.d0-fr_snow(2)*fm)*fv*dts
      aevapb = aevapb + evapb*(1.d0-fr_snow(1))*fb*dts
!!! I guess we need to accumulate evap from snow here
!      aepc=aepc+(epv*fv*dts)
!      aepb=aepb+(epb*fb*dts)
      aepc=aepc+( epv*(1.d0-fr_snow(2)*fm) + epvs*fr_snow(2)*fm )*fv*dts
#ifdef EVAP_VEG_GROUND
      ! next line is probably not correct (need to use something more
      ! sophisticated for epvg)
      aepc=aepc+( epvg*(1.d0-fr_snow(2)) )*fv*dts
#endif
      aepb=aepb+( epb*(1.d0-fr_snow(1)) + epbs*fr_snow(1) )*fb*dts
      !Accumulate GPP, nyk, like evap_tot(2)
      agpp = agpp + gpp*(1.d0-fr_snow(2)*fm)*fv*dts

      dedifs=f(2,1)*tp(2,1)
      if(f(2,1).lt.0.d0) dedifs=f(2,1)*tp(1,1)
      aedifs=aedifs-dts*shw*dedifs*fb
      dedifs=f(2,2)*tp(2,2)
      if(f(2,2).lt.0.d0) dedifs=f(2,2)*tp(1,2)
      aedifs=aedifs-dts*shw*dedifs*fv          ! not used ?
      af0dt=af0dt-dts*(fb*fh(1,1)+fv*fch(0)+htpr)  ! E0 excludes htpr?
      af1dt=af1dt-dts*(fb*fh(2,1)+fv*fh(2,2))
#ifdef TRACERS_WATER
ccc   accumulate tracer fluxes
      atr_evap(:ntg) = atr_evap(:ntg)
     &     + ( tr_evap(:ntg,1)*fb + tr_evap(:ntg,2)*fv )*dts
      atr_rnff(:ntg) = atr_rnff(:ntg)
     &     + ( tr_rnff(:ntg,1)*fb + tr_rnff(:ntg,2)*fv )*dts
C**** no point in setting this here since fm will change before end
c      tot_w1 = fb*( w(1,1)*(1.d0-fr_snow(1))
c     &     + wsn(1,1)*fr_snow(1) )
c     &     + fv*( w(0,2)*(1.d0-fm*fr_snow(2))   ! *fw0  !! remove fw ?
c     &     + wsn(1,2)*fm*fr_snow(2) )
c      if ( tot_w1 > 1.d-30 ) then
c        atr_g(:ntg) =           ! instantaneous       ! atr_g(:ntg) +
c     &       ( fb*( tr_w(:ntg,1,1)*(1.d0-fr_snow(1))
c     &       + tr_wsn(:ntg,1,1) )          !*fr_snow(1)
c     &       + fv*( tr_w(:ntg,0,2)*(1.d0-fm*fr_snow(2)) ! *fw0 !! remove fw?
c     &       + tr_wsn(:ntg,1,2)*fm ) ) /       !*fr_snow(2)
c     &       (rhow * tot_w1)                            ! * dts
c      endif
cddd      print  '(a,100(e12.4))','tr_evap_err',
cddd     &     ( tr_evap(1,1)*fb + tr_evap(1,2)*fv )/1000.d0 -
cddd     &     ( evap_tot(1)*fb + evap_tot(2)*fv ),
cddd     &     ( evap_tot(1)*fb + evap_tot(2)*fv ),
cddd     &     evapvs, fr_snow, fv, fm
!     &     atr_evap(1) - aevap*1000.d0, aevap*1000.d0, evapvs, fr_snow
#endif
      return
!      entry accmf
 777  continue
c provides accumulation units fixups, and calculates
c penman evaporation.  should be called once after
c accumulations are collected.
      aruns=rhow*aruns
      arunu=rhow*arunu
      aflmlt=rhow*aflmlt
      aintercep=rhow*aintercep
      aevapw=rhow*aevapw
      aevapd=rhow*aevapd
      aevapb=rhow*aevapb
      aevap =rhow*aevap
      aepc=rhow*aepc
      aepb=rhow*aepb
      af1dt=af1dt-aedifs
c**** calculation of penman value of potential evaporation, aepp
      h0=fb*(snsh_tot(1)+elh*evap_tot(1))
     &     +fv*(snsh_tot(2)+elh*evap_tot(2))
c     h0=-atrg/dt+srht+trht
ccc   h0=-thrm(2)+srht+trht
      el0=elh*1d-3
      cpfac=sha*rho * ch*vsm
      qsats=qsat(ts,lhe,pres)
      dqdt = dqsatdt(ts,lhe)*qsats
      epen=(dqdt*h0+cpfac*(qsats-qs))/(el0*dqdt+sha)
      aepp=epen*dt
      abetap=1.d0
      if (aepp.gt.0.d0) abetap=(aevapw+aevapd+aevapb)/aepp
      abetap=min(abetap,one)
      abetap=max(abetap,zero)
ccc   computing surface temperature from thermal radiation fluxes
      tbcs = sqrt(sqrt( atrg/(dt*stbo) )) - tfrz
      return
!      entry accm0
 778  continue
c zero out accumulations

      atrg=0.d0                 ! thermal heat from ground
      ashg=0.d0                 ! sensible heat from ground
      aevap=0.d0                ! evaporation from ground
      !alhg=0.d0  ! not accumulated : computed from aevap
      aruns=0.d0
      aeruns=0.d0
      arunu=0.d0
      aerunu=0.d0
      aflmlt=0.d0
      aintercep=0.d0

      abetad=0.d0 ! not accumulated : do we need it?  YES
      abetav=0.d0 ! not accumulated : do we need it?
      abetat=0.d0 ! not accumulated : do we need it?
      abetap=0.d0 ! not sure how it is computed : probably wrong
      abetab=0.d0 ! not accumulated : do we need it?
      abeta=0.d0  ! not accumulated : do we need it?
      acna=0.d0   ! not accumulated : do we need it?
      acnc=0.d0   ! not accumulated : do we need it?
      agpp=0.d0   ! new accumulator, nyk 4/25/03
      aevapw=0.d0              ! evap from wet canopy
      aevapd=0.d0              ! evap from dry canopy
      aevapb=0.d0              ! evap from bare soil (no snow)
      aepc=0.d0                ! potential evap from canopy
      aepb=0.d0                ! potential evap from bare soil (no snow)
      aedifs=0.d0              ! heat transport by water
      af0dt=0.d0               ! heat from ground - htpr
      af1dt=0.d0               ! heat from 2-nd soil layer
      ! aepp=0.d0 ! not accumulated : computed in accmf
#ifdef TRACERS_WATER
ccc   tracers
      atr_evap(:ntg) = 0.d0
      atr_rnff(:ntg) = 0.d0
      atr_g(:ntg) = 0.d0
#endif

      return
      end subroutine accm

      subroutine gdtm(dtm)
c**** calculates the maximum time step allowed by stability
c**** considerations.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      real*8 ak1,ak2(2),betas(2),cna,xk2(2),dldz2,dqdt,rho3,sgmm
      real*8, intent(out) :: dtm
      real*8 dtm1,dtm2,dtm3,dtm4,t450,xk1
      integer ibv,k
      t450=450.d0
      dqdt=dqsatdt(ts,lhe)*qsat(ts,lhe,pres)
c****
c**** first calculate timestep for water movement in soil.
      sgmm=1.0d0
      dldz2=0.d0
      do ibv=i_bare,i_vege
        do k=1,n
          dldz2=max(dldz2,d(k,ibv)/dz(k)**2)
        end do
      end do
      dtm=sgmm/(dldz2+1d-12)
      if(q(4,1).gt.0.d0)dtm=min(dtm,t450)
      dtm1=dtm
      if ( dtm .lt. 0.d0 ) call stop_model('gdtm: dt1_ghy<0',255)
c****
c**** next calculate timestep for heat movement in soil.
      do ibv=i_bare,i_vege
        do k=1,n
          xk1=xkh(k,ibv)
          ak1=(shc(k,ibv)+((1.d0-fice(k,ibv))*shw+fice(k,ibv)*shi)
     &         *w(k,ibv))/dz(k)
          dtm=min(dtm,.5d0*ak1*dz(k)**2/(xk1+1d-12))
        end do
      end do
      dtm2=dtm
      if ( dtm .lt. 0.d0 ) call stop_model('gdtm: dt2_ghy<0',255)
c****
c**** finally, calculate max time step for top layer bare soil
c**** and canopy interaction with surface layer.
c**** use timestep based on coefficient of drag
      cna=ch*vsm
      rho3=.001d0*rho
      betas(1:2) = 1.d0 ! it''s an overkill but it makes the things
                        ! simpler.
      if(epb.le.0.d0)then
       betas(1)=1.0d0
      else
       betas(1)=evapb/epb
      endif
      if(epv.le.0.d0)then
       betas(2)=1.0d0
      else
       betas(2)=(evapvw*fw+evapvd*(1.d0-fw))/epv
      endif
      do ibv=i_bare,i_vege
        k=2-ibv
        xk2(ibv)=sha*rho*cna
     &       + betas(ibv)*rho3*cna*elh*dqdt
     &       + 8.d0*stbo*(tp(k,ibv)+tfrz)**3
        ak2(ibv)=shc(k,ibv)+((1.d0-fice(k,ibv))*shw+fice(k,ibv)*shi)
     &       *w(k,ibv)
        dtm=min(dtm,0.5*ak2(ibv)/(xk2(ibv)+1d-12))
        if(ibv.eq.1)dtm3=dtm
        if(ibv.eq.2)dtm4=dtm
c
c prevent oscillation of top snow layer
c     if(isn(ibv).ne.0.or.snowd(ibv).ne.0.d0)then
c      ak3(ibv)=.05d0*shi*spgsn
c      dtm=min(dtm,ak3(ibv)/(xk2(ibv)+1.d-12))
c      if(ibv.eq.1)dtm5=dtm
c      if(ibv.eq.2)dtm6=dtm
c     endif
      end do
      if(dtm.lt.5.d0)then
       write(99,*) '*********** gdtm: ijdebug,fb,fv',ijdebug,fb,fv
       write(99,*)'dtm',dtm1,dtm2,dtm3,dtm4
       write(99,*)'xk2',xk2
       write(99,*)'ak2',ak2
       write(99,*)'snsh',snsh
       write(99,*)'xlth',elh*evap_tot(1:2)
       write(99,*)'dqdt',dqdt
       write(99,*)'ts,tfrz',ts,tfrz
       write(99,*)'dlt',tp(1,1)-ts+tfrz,tp(0,2)-ts+tfrz
       call stop_model("gdtm: time step < 5 s",255)
      endif
c****
      return
      end subroutine gdtm


      subroutine outw(i)
c**** prints current state of soil for one grid box
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      use filemanager, only: openunit
      integer i,k
      integer, save :: ichn = 0
      call wtab
      if( ichn == 0 ) call openunit("soil_outw", ichn)
      write(ichn,1000)
      write(ichn,*)'general quantities (bare soil or vegetation)'
      write(ichn,*)'dts,tag',dts,i,' after qsbal(0),retp(1),advnc(2)'
cc    write(ichn,1021)
      write(ichn,1045)
      write(ichn,1023)'i= ',ijdebug/1000,'pr= ',pr,'ts= ',ts-tfrz,
     *     'q1= ',0.d0
      write(ichn,1023)'j= ',mod(ijdebug,1000),
     *     'snowf= ',drips(1),'tg= ',0.d0-tfrz,'qs= ',qs
      write(ichn,1043)'t1= ',0.d0-tfrz,'vg= ',0.d0,'ch= ',ch
      write(ichn,1044)'vsm= ',vsm
      write(ichn,1022)
      write(ichn,1021)
      write(ichn,1014)'bare soil   fb = ',fb
 1014 format(1x,a17,f4.2)
cc    write(ichn,1021)
      write(ichn,1025)
      write(ichn,1026)
      write(ichn,1027)
     &     snowd(1),rnf(1),evap_tot(1),xinfc(1),zw(1),debug_data%qb
      write(ichn,1021)
      write(ichn,1030)
      write(ichn,1031)
      do 100 k=1,n
      write(ichn,1040)k,theta(k,1),tp(k,1),fice(k,1),rnff(k,1),f(k,1),
     & h(k,1),xk(k,1),w(k,1),ws(k,1),shc(k,1),fh(k,1),ht(k,1),
     & q(1,k),q(2,k),q(3,k),q(4,k)
  100 continue
cc    write(ichn,1021)
      write(ichn,1022)
      write(ichn,1021)
      write(ichn,1014)'vegetation  fv = ',fv
cc    write(ichn,1021)
      write(ichn,1035)
      write(ichn,1036)
      write(ichn,1037) snowd(2),rnf(2),evap_tot(2),xinfc(2),zw(2),
     &     debug_data%qc,evapvw,
     *     evapvd,dripw(2)+drips(2),fw
      write(ichn,1021)
      write(ichn,1030)
      write(ichn,1031)
      k=0
      write(ichn,1049)k,theta(k,2),tp(k,2),fice(k,2),     0.d0,f(k,2),
     & w(k,2),ws(k,2),shc(k,2),fh(k,2),ht(k,2)
      do 200 k=1,n
      write(ichn,1040)k,theta(k,2),tp(k,2),fice(k,2),rnff(k,2),f(k,2),
     & h(k,2),xk(k,2),w(k,2),ws(k,2),shc(k,2),fh(k,2),ht(k,2),
     & q(1,k),q(2,k),q(3,k),q(4,k)
  200 continue
cc    write(ichn,1021)
      write(ichn,1022)
      write(ichn,1021)
      write(ichn,1055)
 1055 format(1x,3x,9x,'   kgm-2',2x,9x,'    kgm-2',3x,9x,'   kgm-2',
     *     4x,9x,'1e6jm-2',3x,9x,'1e6jm-2')
      write(ichn,1060) aruns,aevapw,aepc,0.d0,af0dt
 1060 format(1x,3x,'aruns = ',f9.4,2x,'aevapw = ',0pf9.4,4x,'aepc = ',
     *     0pf9.4,4x,'afhg = ',-6pf9.4,2x,'af0dt = ',-6pf9.4)
      write(ichn,1065) arunu,aevapd,aepb,atrg,htpr*dts
 1065 format(1x,3x,'arunu = ',f9.4,2x,'aevapd = ',0pf9.4,4x,'aepb = ',
     *     0pf9.4,4x,'atrg = ',-6pf9.4,2x,'aphdt = ',-6pf9.4)
      write(ichn,1070) aevapb,ashg,aeruns
 1070 format(1x,3x,'        ',9x,2x,'aevapb = ',0pf9.4,4x,'       ',
     *     9x,4x,'ashg = ',-6pf9.4,2x,'aerns = ',-6pf9.4)
      write(ichn,1073) alhg,af1dt,aedifs
 1073 format(1x,3x,'        ',9x,2x,'         ',9x,4x,'       ',
     *     9x,4x,'alhg = ',-6pf9.4,2x,'af1dt = ',-6pf9.4/
     *  1x,3x,'        ',9x,2x,'         ',9x,4x,'       ',
     *     9x,4x,'       ',2x,9x,'aedfs = ',-6pf9.4)
c**** more outw outputs
      write(ichn,*)'thrm ',thrm_tot
      write(ichn,*)'xlth ',elh*evap_tot(1:2)
      write(ichn,*)'snsh ',snsh
      write(ichn,*)'htpr,srht,trht ',htpr,srht,trht
      return
1000  format(' ',121('='))
1001  format('1')
1010  format(1x,a20,f10.0)
1020  format(1x,a20,1pe12.4)
1021  format('0')
1022  format(' ',60('. '),'.')
1023  format(1x,a10,i10,a10,6pf10.2,2(a10,0pf8.2),a10,0pf8.4)
1043  format(1x,10x,10x,10x,10x,2(a10,0pf8.2),a10,0pf8.4)
1044  format(1x,10x,10x,38x,1(a10,f8.2))
1045  format(1x,20x,10x,'  1e-6ms-1',10x,4x,'t(c)',10x,4x,'ms-1')
1024  format(1x,6(a10,e10.2))
1015  format(1x,4(a8,f8.2))
1019  format(1x,12x,4(a8,f8.2))
1025  format(' ',5x,'snowd',7x,'rnf',6x,'evap',6x,'xinfc',
     *     8x,'zw',8x,'qb')
1026  format(' ','      mh2o',2x,'1e-6ms-1',2x,'1e-6ms-1',3x,
     *     '1e-6ms-1',3x,'      m',5x,'     ')
1027  format(' ',0pf10.4,6pf10.4,6pf10.4,1x,6pf10.1,0pf10.4,
     *     0pf10.4)
1030  format(' ',5x,'theta',3x,'tp',2x,'fice',4x,'runoff'
     & ,8x,'fl',9x,'h',8x,'xk',6x,'w',5x,'ws',8x,
     & 'shc',8x,'fh',8x,'ht',1x,'sand',1x,'loam',1x,'clay',1x,'peat')
1031  format(' ',5x,5x,2x,'(c)',2x,6x,'1e-6ms-1',2x,'1e-6ms-1',
     &     3x,'      m',2x,'1e-6ms-1',1x,'     m',1x,'     m',
     &     1x,'1e6jm-3c-1',4x,'  wm-2',3x,'1e6jm-2',4x,'%',4x,'%',4x,'%'
     &     ,4x,'%'/1x,125('-'))
1035  format(' ',5x,'snowd',7x,'rnf',6x,'evap',6x,'xinfc',
     *     8x,'zw',8x,'qc',5x,'evapw',5x,'evapd',8x,'dr',8x,'fw')
1036  format(' ','      mh2o',2x,'1e-6ms-1',2x,'1e-6ms-1',3x,
     *     '1e-6ms-1',2x,'       m',5x,'     ',2x,'1e-6ms-1',2x,
     *     '1e-6ms-1',2x,'1e-6ms-1','          ')
1037  format(' ',0pf10.4,6pf10.4,6pf10.4,1x,6pf10.1,0pf10.4,
     *     0pf10.4,6pf10.4,6pf10.4,6pf10.4,0pf10.2)
1040  format(1x,i3,f7.3,f5.1,f6.3,1p,6pf10.4,6pf10.4,0pf10.3,6pf10.4,
     *     0pf7.4,0pf7.4,1x,-6pf10.4,0pf10.4,-6pf10.4,4(2pf5.1))
1049  format(1x,i3,f7.3,f5.1,f6.3,1p,6pf10.4,6pf10.4,10x,10x,
     *     0pf7.4,0pf7.4,1x,-6pf10.4,0pf10.4,-6pf10.4,3(2pf5.1))
      end subroutine outw
c****
      subroutine wtab
c**** returns water table zw for ibv=1 and 2.d0
c**** input:
c**** zb - layer boundaries, m
c**** zc - soil centers, m
c**** dz - layer thicknesses, m
c**** h - soil potential of layers, m
c**** f - fluxes between layers, m s-1
c**** xk - conductivities of layers, m s-1
c**** output:
c**** zw(2) - water table for ibv=1 and 2, m
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      real*8 denom,hmat,tol
      integer ibv,k
      tol=1d-6
      do 100 ibv=1,2
c**** find non-saturated layer
      do 10 k=n,1,-1
      if(w(k,ibv).lt.ws(k,ibv)*(1.d0-tol))go to 20
   10 continue
      k=1
   20 continue
c**** retrieve matric potential
c     write(6,*)'ij,n,k,hmat,ibv,xkl,ibv',ijdebug,n,k,hmat,ibv,xk(k,ibv)
      hmat=h(k,ibv)-zc(k)
c**** calculate denominator, and keep zw above zb(k+1)
      if(xk(k,ibv).le.1d-20) then
           denom=-2.d0*hmat/dz(k)
           go to 90
           end if
      denom=max(f(k,ibv)/xk(k,ibv)+1.d0,-2.d0*hmat/dz(k))
   90 continue
c**** calculate water table
c     write(6,*) 'denom',denom
      zw(ibv)=zb(k)-sqrt(-2.d0*hmat*dz(k)/(denom+1d-20))
  100 continue
      return
      end subroutine wtab


      subroutine snow
      use snow_drvm, only: snow_drv
      implicit none
      real*8 fmask(2), evapsn(2), devapsn_dt(2), canht
      integer ibv

      fmask(1) = 1.d0
      fmask(2) = fm

      evapsn(1) = evapbs
      evapsn(2) = evapvs
      devapsn_dt(1) = devapbs_dt
      devapsn_dt(2) = devapvs_dt

      canht = stbo*(tp(0,2)+tfrz)**4

!!! correction for canopy transmittance
#ifdef RAD_VEG_GROUND
      canht = canht*(1.d0-trans_sw) + (srht+trht)*trans_sw
#endif

ccc reset some fluxes for ibv hack
      flmlt_scale(:) = 0.d0
      fhsng_scale(:) = 0.d0
      flmlt(:) = 0.d0
      fhsng(:) = 0.d0
      thrmsn(:) = 0.d0
      flux_snow(:,:) = 0.d0

ccc remember initial snow water for tracers
      wsn_for_tr(:,1) = wsn(:,1)*fr_snow(1)
      wsn_for_tr(:,2) = wsn(:,2)*fr_snow(2)

      do ibv=i_bare,i_vege
        call snow_drv(
ccc input
     &       fmask(ibv), evapsn(ibv), snshs(ibv), srht, trht, canht,
     &       drips(ibv), dripw(ibv), htdrips(ibv), htdripw(ibv),
     &       devapsn_dt(ibv), dsnsh_dt, dts,
     &       tp(1,ibv), dz(1), nlsn, top_stdev,
ccc updated
     &       dzsn(1,ibv), wsn(1,ibv), hsn(1,ibv), nsn(ibv),
     &       fr_snow(ibv),
ccc output
     &       flmlt(ibv), fhsng(ibv), flmlt_scale(ibv),
     &       fhsng_scale(ibv), thrmsn(ibv), flux_snow(0,ibv)
     &       )
      enddo

      evapbs = evapsn(1)
      evapvs = evapsn(2)

      end subroutine snow

      subroutine set_snow
!@sum set_snow extracts snow from the first soil layer and initializes
!@+   snow model prognostic variables
!@+   should be called when model restarts from the old restart file
!@+   ( which doesn''t contain new snow model (i.e. 3 layer) data )
c
c input:
c snowd(2) - landsurface snow depth
c w(k,2)   - landsurface water in soil layers
c ht(k,2)  - landsurface heat in soil layers
c fsn      - heat of fusion
c shi      - specific heat of ice
c shc(k,2) - heat capacity of soil layers
c
c output:
c dzsn(lsn,2) - snow layer thicknesses
c wsn(lsn,2)  - snow layer water equivalent depths
c hsn(lsn,2)  - snow layer heat contents
c tsn1(2)     - snow top temperature
c isn(2)      - 0 if no snow, 1 if snow
c nsn(2)      - number of snow layers
c snowd(2)
c w(k,2)
c ht(k,2)
c
c calling sequence:
c
c     assignment of w,ht,snowd
c     call ghinij(i,j,wfc1)
c     call set_snow
c note: only to be called when initializing from landsurface
c       prognostic variables without the snow model.
c
      use snow_model, only: snow_redistr, snow_fraction
      integer ibv

c outer loop over ibv
      do ibv=1,2

c initalize all cases to nsn=1
        nsn(ibv)=1

ccc since we don''t know what kind of data we are dealing with,
ccc better check it

        if( snowd(ibv) .gt. w(1,ibv)-dz(1)*thetm(1,ibv)  ) then
          write(99,*) 'snowd corrected: old=', snowd(ibv)
          snowd(ibv) = w(1,ibv)-dz(1)*thetm(1,ibv) - 1.d-10
          write(99,*) '                 new=', snowd(ibv)
          if ( snowd(ibv) .lt. -0.001d0 )
     &         call stop_model('set_snow: neg. snow',255)
          if ( snowd(ibv) .lt. 0.d0 ) snowd(ibv) = 0.d0 ! rounding error
        endif

c if there is no snow, set isn=0.  set snow variables to 0.d0
        if(snowd(ibv).le.0.d0)then
          dzsn(1,ibv)=0.d0
          wsn(1,ibv)=0.d0
          hsn(1,ibv)=0.d0
          tsn1(ibv)=0.d0
          fr_snow(ibv) = 0.d0
        else

c given snow, set isn=1.d0
c!!!        dzsn(1,ibv)=snowd(ibv)/spgsn
c!!!  replacing prev line considering rho_snow = 200
          dzsn(1,ibv)=snowd(ibv) * 5.d0
          wsn(1,ibv)=snowd(ibv)
c!!! actually have to compute fr_snow and modify dzsn ...
          fr_snow(ibv) = 1.d0

c given snow, temperature of first layer can''t be positive.
c the top snow temperature is the temperatre of the first layer.
          if(fsn*w(1,ibv)+ht(1,ibv).lt.0.d0)then
            tsn1(ibv)=(ht(1,ibv)+w(1,ibv)*fsn)/(shc(1,ibv)+w(1,ibv)*shi)
          else
            tsn1(ibv)=0.d0
          endif

c use snow temperature to get the heat of the snow
          hsn(1,ibv)=tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn

c subtract the snow from the landsurface prognositic variables
          w(1,ibv)=w(1,ibv)-wsn(1,ibv)
          ht(1,ibv)=ht(1,ibv)-hsn(1,ibv)

ccc and now limit all the snow to 5cm water equivalent
          if ( snowd(ibv) .gt. 0.05d0 ) then
            snowd(ibv) = 0.05d0
            dzsn(1,ibv)= snowd(ibv) * 5.d0
            wsn(1,ibv)= snowd(ibv)
            hsn(1,ibv)= tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn
          endif
ccc and now limit all the snow to 50cm water equivalent (for debug)
cddd          if ( snowd(ibv) .gt. 0.0005d0 ) then
cddd            snowd(ibv) = 0.5d0
cddd            dzsn(1,ibv)= snowd(ibv) * 5.d0
cddd            wsn(1,ibv)= snowd(ibv)
cddd            hsn(1,ibv)= tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn
cddd          endif

ccc redistribute snow over the layers and recompute fr_snow
ccc (to make the data compatible with snow model)
          if ( .not. ( hsn(1,ibv) > 0. .or.  hsn(1,ibv) <= 0. ) )
     &        call stop_model("ERR in init_snow: NaN", 255)

          call snow_fraction(dzsn(:,ibv), nsn(ibv), 0.d0, 0.d0,
     &         1.d0, fr_snow(ibv) )
          if ( fr_snow(ibv) > 0.d0 ) then
            call snow_redistr(dzsn(:,ibv), wsn(:,ibv), hsn(:,ibv),
     &           nsn(ibv), 1.d0/fr_snow(ibv) )
            if ( .not. ( hsn(1,ibv) > 0. .or.  hsn(1,ibv) <= 0. ) )
     &        call stop_model("ERR in init_snow 2: NaN", 255)
          else
            snowd(ibv) = 0.d0
            dzsn(1,ibv)=0.d0
            wsn(1,ibv)=0.d0
            hsn(1,ibv)=0.d0
            tsn1(ibv)=0.d0
            fr_snow(ibv) = 0.d0
          endif

        endif
!!! debug : check if  snow is redistributed correctly
        if ( dzsn(1,ibv) > 0.d0 .and. dzsn(1,ibv) < .099d0) then
          call stop_model("set_snow: error in dz",255)
        endif

        if ( dzsn(1,ibv) > 0.d0 ) then
          if (  wsn(1,ibv)/dzsn(1,ibv)  < .1d0)
     &         call stop_model("set_snow: error in dz",255)
        endif


      enddo

            if ( .not. ( hsn(1,1) > 0. .or.  hsn(1,1) <= 0. ) )
     &        call stop_model("ERR in init_snow 2: NaN", 255)
            if ( .not. ( hsn(1,2) > 0. .or.  hsn(1,2) <= 0. ) )
     &        call stop_model("ERR in init_snow 2: NaN", 255)



      return
      end subroutine set_snow


#ifdef TRACERS_WATER
      subroutine check_wc(wc)
      real*8 wc
!!! removing test...
      return
!!!! ;-(
      !wc = 1000.d0
      !return

      if ( abs(wc-1000.d0) > 1.d-3 ) then
        print *,"GHYTR: wc= ", wc, (wc-1000.d0)/1000.d0
        !call stop_model("GHY: wc != 1000",255)
      endif

      end subroutine check_wc


      subroutine ghy_tracers_drydep
!@sum applies dry deposit flux to tracers in upper layers of
!@+   canopy, soil and snow
ccc input from outside:
ccc trdd(ntgm) - flux of tracers as dry deposit (kg /m^2 /s)
!@var m number of tracers (=ntg)
      integer m

      m = ntg

      if ( process_vege ) then
        ! +drydep
        tr_w(:m,0,2) = tr_w(:m,0,2) + trdd(:m)*(1-fm*fr_snow(2))*dts
        tr_wsn(:m,1,2) = tr_wsn(:m,1,2) + trdd(:m)*fm*fr_snow(2)*dts
      endif

      if ( process_bare ) then
        ! +drydep
        tr_w(:m,1,1) = tr_w(:m,1,1) + trdd(:m)*(1-fr_snow(1))*dts
        tr_wsn(:m,1,1) = tr_wsn(:m,1,1) + trdd(:m)*fr_snow(1)*dts
      endif

      end subroutine ghy_tracers_drydep


      subroutine ghy_tracers
ccc will need the following data from GHY:
ccc rnff(k,ibv) - underground runoff fro layers 1:n
ccc rnf(ibv) - surface runoff
ccc evapd*betadl(k)/(betad+1d-12) - transpiration from 1:n (ibv=2 only)
ccc f(k,ibv) - flux between layers (+up) (upper edge)

ccc will need from SNOW:
ccc tr_flux(0:nsl) - flux down

ccc will do it as follows:
ccc 1. propagate tracers to canopy layer (no flux up on lower boundary)
ccc 2. propagate tracers through the snow (no flux up from soil)
ccc 3. propagate them through the soil

ccc input from outside:
ccc pr - precipitation (m/s) ~ (10^3 kg /m^2 /s)
ccc trpr(ntgm) - flux of tracers (kg /m^2 /s)
ccc trdd(ntgm) - flux of tracers as dry deposit (kg /m^2 /s)
ccc tr_surf(ntgm) - surface concentration of tracers (kg/m^3 H2O)

ccc **************************************************

ccc !!! test

ccc internal vars:
!@var wi instantaneous amount of water in soil layers (m)
!@var tr_wc ratio tracers/water in soil layers (?/m)
!@var tr_wcc ratio tracers/water in canopy (?/m)
      real*8 wi(0:ngm,1:2),tr_wc(ntgm,0:ngm,1:2),tr_wcc(ntgm,1:2)
!@var wsni instantaneous amount of water in snow layers (m)
!@var tr_wsnc ratio tracers/water in snow layers (?/m)
      real*8 wsni(nlsn,2),tr_wsnc(ntgm,0:nlsn,2)
!@var flux_snow_tot water flux between snow layers * fr_snow (m/s)
      real*8 flux_snow_tot(0:nlsn,2)
!@var flux_dt amount of water moved in current context (m)
      real*8 flux_dt, err, evap_tmp(ntgm)
      integer ibv,i
!@var m = ntg number of ground tracers (to make notations shorter)
      integer m,mc
      real*8, parameter :: EPSF = 1.d-20 ! 1.d-11
#ifdef TRACERS_SPECIAL_O18
      real*8, external :: fracvl
#endif
!@var tol error tolerance for each type of tracer
      real*8 tol(ntg)

ccc for debug
      real*8 tr_w_o(ntgm,0:ngm,2), tr_wsn_o(ntgm,nlsn,2)

      if ( ntg < 1 ) return   ! no water tracers


      !m = 2 !!! testing
      m = ntg

ccc set error tolerance for each time of tracer (can be increased
ccc if model stops due to round-off error)
      do mc=1,m
         tol(mc) = ( sum(tr_w(mc,:,:))/(size(tr_w,2)*size(tr_w,3))
     &       +   sum(tr_wsn(mc,:,:))/(size(tr_wsn,2)*size(tr_wsn,3))
     &       + trpr(mc)*dts + trdd(mc)*dts + tr_surf(mc)*1.d-8*dts )
     &       *1.d-10 + 1.d-32
      enddo

!!! for test only
!      trpr(:) = .5d0 * pr
  !    tr_surf(:) = 1000.d0 !!!???
      !tr_w = 0.d0
      !tr_wsn = 0.d0
#ifdef DEBUG_TR_WITH_WATER
      if ( abs(tr_surf(1)-1000.) > 1.d-5 )
     &     call stop_model("GHY: tr_surf != 1000",255)
#endif

cddd      print *, 'starting ghy_tracers'
cddd      print *, 'trpr, trpr_dt ', trpr(1), trpr(1)*dts
cddd      print *, 'tr_w 1 ', tr_w(1,:,1)
cddd      print *, 'tr_w 2 ', tr_w(1,:,2)
cddd      print *, 'tr_wsn 1 ', tr_wsn(1,:,1)
cddd      print *, 'tr_wsn 2 ', tr_wsn(1,:,2)

!!! test
 !     tr_wsn(1,:,1) = wsn_for_tr(:,1) * fr_snow(1) * 1000.d0
 !     tr_wsn(1,:,2) = wsn_for_tr(:,2) * fr_snow(2) * 1000.d0


      tr_w_o(:m,:,i_bare:i_vege) = tr_w(:m,:,i_bare:i_vege)
      tr_wsn_o(:m,:,i_bare:i_vege) = tr_wsn(:m,:,i_bare:i_vege)

ccc set internal vars
      do ibv=i_bare,i_vege
        wi(0:n,ibv) = w(0:n,ibv)
        wsni(1:nlsn,ibv) = wsn_for_tr(1:nlsn,ibv) ! *fr_snow(ibv)
        flux_snow_tot(0:nlsn,ibv) = flux_snow(0:nlsn,ibv)*fr_snow(ibv)
!!! hack
!!! snow doesnt suck water from ground
       if ( flux_snow_tot(nlsn,ibv) < 0.d0 ) then
         flux_snow_tot(0:nlsn-1,ibv) = flux_snow_tot(0:nlsn-1,ibv)
     &     - flux_snow_tot(nlsn,ibv)
         flux_snow_tot(nlsn,ibv) = 0.d0
       endif
      enddo

ccc reset accumulators to 0
      tr_evap(:m,1:2) = 0.d0
      tr_rnff(:m,1:2) = 0.d0

ccc apply dry deposit tracers here (may need to be removed outside
ccc but will keep it here for now for compatibility with conservation tests
      call ghy_tracers_drydep

ccc canopy
      tr_wcc(:m,:) = 0.d0
 !!!test
 !!!     tr_wcc(:m,:) = 1000.d0

      if ( process_vege ) then
      ! +precip
        tr_w(:m,0,2) = tr_w(:m,0,2) + trpr(:m)*dts
        wi(0,2) = wi(0,2) + pr*dts
        if ( wi(0,2) > 0.d0 ) tr_wcc(:m,2) = tr_w(:m,0,2)/wi(0,2)
        call check_wc(tr_wcc(1,2))

      ! +- evap
      if ( evapvw >= 0.d0 ) then  ! no dew
        evap_tmp(:m) = fc(0)+pr
#ifdef TRACERS_SPECIAL_O18
        if ( evap_tmp(1)*dts < wi(0,2) .and. tp(0,2) > 0.d0 ) then
          do mc=1,m
ccc tr_name - loop over string comparisons deep inside the nested
ccc loops...
            evap_tmp(mc) = evap_tmp(mc) * fracvl( tp(0,2),tr_name(mc) )
          enddo
        endif
#endif
        tr_evap(:m,2) = tr_evap(:m,2) + evap_tmp(:m)*tr_wcc(:m,2)
        tr_w(:m,0,2) = tr_w(:m,0,2) - evap_tmp(:m)*tr_wcc(:m,2)*dts
      else  ! dew adds tr_surf to canopy
        tr_evap(:m,2) = tr_evap(:m,2) + (fc(0)+pr)*tr_surf(:m)
        tr_w(:m,0,2) = tr_w(:m,0,2) - (fc(0)+pr)*tr_surf(:m)*dts
        wi(0,2) = wi(0,2) - (fc(0)+pr)*dts
        tr_wcc(:m,2) = tr_w(:m,0,2)/wi(0,2)
      endif
      call check_wc(tr_wcc(1,2))

      ! -drip
        tr_w(:m,0,2) = tr_w(:m,0,2) + fc(1)*tr_wcc(:m,2)*dts
      endif

      ! trivial value for bare soil
      if ( pr > 0.d0 ) tr_wcc(:m,1) = trpr(:m)/pr
      call check_wc(tr_wcc(1,1))
ccc end canopy


ccc snow
  !!!>>    tr_wsnc(:m,:,:) = 0.d0 !!! was 0
      tr_wsnc(:m,:,:) = 0.d0 !!! was 1000

      ! dew
      if ( process_bare .and. evapbs < 0.d0 ) then
        flux_dt = - evapbs*fr_snow(1)*dts
        tr_evap(:m,1) = tr_evap(:m,1) - flux_dt/dts*tr_surf(:m)
        tr_wsn(:m,1,1) = tr_wsn(:m,1,1) + flux_dt*tr_surf(:m)
        wsni(1,1) = wsni(1,1) + flux_dt
      endif
      if ( process_vege .and. evapvs < 0.d0 ) then
        flux_dt = - evapvs*fm*fr_snow(2)*dts
        tr_evap(:m,2) = tr_evap(:m,2) - flux_dt/dts*tr_surf(:m)
        tr_wsn(:m,1,2) = tr_wsn(:m,1,2) + flux_dt*tr_surf(:m)
        wsni(1,2) = wsni(1,2) + flux_dt
      endif

      do ibv=i_bare,i_vege
        ! init tr_wsnc
        do i=1,nlsn
          if ( wsni(i,ibv) > 0.d0 ) then
            tr_wsnc(:m,i,ibv) = tr_wsn(:m,i,ibv)/wsni(i,ibv)
          else
            tr_wsnc(:m,i,ibv) = 0.d0
          endif
        enddo
        if ( fr_snow(ibv) > 0.d0 ) then  ! process snow
          flux_dt = (drips(ibv) + dripw(ibv)*fr_snow(ibv))*dts
          tr_wsn(:m,1,ibv) = tr_wsn(:m,1,ibv) + tr_wcc(:m,ibv)*flux_dt
          wsni(1,ibv) = wsni(1,ibv) + flux_dt
          if ( wsni(1,ibv) > 0.d0 )
     &         tr_wsnc(:m,1,ibv) = tr_wsn(:m,1,ibv)/wsni(1,ibv)
!!!
       !   tr_wsnc(:m,1,ibv) = 1000.d0
          call check_wc(tr_wsnc(1,1,ibv))
          ! sweep down
          do i=2,nlsn
            if ( flux_snow_tot(i-1,ibv) > EPSF ) then
              flux_dt = flux_snow_tot(i-1,ibv)*dts
              tr_wsn(:m,i,ibv) = tr_wsn(:m,i,ibv)
     &             + tr_wsnc(:m,i-1,ibv)*flux_dt
              wsni(i,ibv) = wsni(i,ibv) + flux_dt
              if ( wsni(i,ibv) > 0.d0 )
     &             tr_wsnc(:m,i,ibv) = tr_wsn(:m,i,ibv)/wsni(i,ibv)
!!!
       !       tr_wsnc(:m,i,ibv) = 1000.d0
              call check_wc(tr_wsnc(1,i,ibv))
              tr_wsn(:m,i-1,ibv) = tr_wsn(:m,i-1,ibv)
     &             - tr_wsnc(:m,i-1,ibv)*flux_dt
              wsni(i-1,ibv) = wsni(i-1,ibv) + flux_dt   ! extra?
            endif
          enddo
          ! sweep up
          do i=nlsn-1,1,-1
            if ( flux_snow_tot(i,ibv) < -EPSF ) then
              flux_dt = - flux_snow_tot(i,ibv)*dts
              tr_wsn(:m,i,ibv) = tr_wsn(:m,i,ibv)
     &             + tr_wsnc(:m,i+1,ibv)*flux_dt
              wsni(i,ibv) = wsni(i,ibv) + flux_dt
              if ( wsni(i,ibv) > 0.d0 )
     &             tr_wsnc(:m,i,ibv) = tr_wsn(:m,i,ibv)/wsni(i,ibv)
!!!
       !       tr_wsnc(:m,i,ibv) = 1000.d0
              call check_wc(tr_wsnc(1,i,ibv))
              tr_wsn(:m,i+1,ibv) = tr_wsn(:m,i+1,ibv)
     &             - tr_wsnc(:m,i+1,ibv)*flux_dt
              wsni(i+1,ibv) = wsni(i+1,ibv) - flux_dt   ! extra?
            endif
          enddo
          !!!flux_dt = flmlt(ibv)*fr_snow(ibv)*dts
          flux_dt = flux_snow_tot(nlsn,ibv)*dts
          if ( abs(flux_dt-flmlt(ibv)*fr_snow(ibv)*dts) > 1.d-12 )
     &         print *,"GHYTR: tracers flux_snow_tot problem",
     &         flux_dt-flmlt(ibv)*fr_snow(ibv)*dts
          tr_w(:m,1,ibv) = tr_w(:m,1,ibv)
     &         + tr_wsnc(:m,nlsn,ibv)*flux_dt
          wi(1,ibv) = wi(1,ibv) + flux_dt
          tr_wsn(:m,nlsn,ibv) = tr_wsn(:m,nlsn,ibv)
     &         - tr_wsnc(:m,nlsn,ibv)*flux_dt
        else  ! no snow
          tr_w(:m,1,ibv) = tr_w(:m,1,ibv)
     &         + sum( tr_wsn(:m,1:nlsn,ibv), 2 )
!     &         + tr_wcc(:m,ibv)*flmlt_scale(ibv)*dts
     &         + tr_wcc(:m,ibv)*drips(ibv)*dts
          tr_wsn(:m,1:nlsn,ibv) = 0.d0
          wi(1,ibv) = wi(1,ibv) + flmlt_scale(ibv)*dts
        endif
        ! evap
        if ( ibv == 1 ) then
          flux_dt = max( evapbs*fr_snow(1)*dts, 0.d0 )
        else
          flux_dt = max( evapvs*fm*fr_snow(2)*dts, 0.d0 )
        endif
!!! test
!        tr_wsnc(:m,1,ibv) = 1000.d0
        tr_wsn(:m,1,ibv) = tr_wsn(:m,1,ibv) - tr_wsnc(:m,1,ibv)*flux_dt
        tr_evap(:m,ibv)= tr_evap(:m,ibv) + tr_wsnc(:m,1,ibv)*flux_dt/dts
      enddo

!!!! for test
  !    tr_wsn(1,:,1) = wsn(:,1) * fr_snow(1) * 1000.d0
  !    tr_wsn(1,:,2) = wsn(:,2) * fr_snow(2) * 1000.d0
cddd      print '(a,10(e12.4))', 'tr_wsn_1 '
cddd     &     , tr_wsn(1,:,1) - wsn(:,1) * fr_snow(1) * 1000.d0
cddd     &     , fr_snow(1), nsn(1)+0.d0 ,flux_snow_tot(0:nlsn,1)*dts
cddd      print '(a,10(e12.4))', 'tr_wsn_2 '
cddd     &     , tr_wsn(1,:,2) - wsn(:,2) * fr_snow(2) * 1000.d0
cddd     &     , fr_snow(2), nsn(2)+0.d0, flux_snow_tot(0:nlsn,2)*dts

ccc soil layers
  !>>    tr_wc(:m,1:n,1:2) = 0.d0
  !>    tr_wc(:m,1:n,1:2) = 1000.d0  !!! was 0

      ! add dew to bare soil
      if ( process_bare .and. evapb < 0.d0 ) then
        flux_dt = - evapb*(1.d0-fr_snow(1))*dts
        tr_evap(:m,1) = tr_evap(:m,1) - flux_dt/dts*tr_surf(:m)
        tr_w(:m,1,1) = tr_w(:m,1,1) + flux_dt*tr_surf(:m)
        wi(1,1) = wi(1,1) + flux_dt
      endif

      do ibv=i_bare,i_vege
        ! initial tr_wc
        do i=1,n
          if ( wi(i,ibv) > 0.d0 ) then
            tr_wc(:m,i,ibv) = tr_w(:m,i,ibv)/wi(i,ibv)
          else ! actually 'else' should never happen
            tr_wc(:m,i,ibv) = 0.d0
          endif
        enddo
        ! precip
        flux_dt = dripw(ibv)*(1.d0-fr_snow(ibv))*dts
        tr_w(:m,1,ibv) = tr_w(:m,1,ibv) +  tr_wcc(:m,ibv)*flux_dt
        wi(1,ibv) = wi(1,ibv) + flux_dt
        if ( wi(1,ibv) > 0.d0 )
     &       tr_wc(:m,1,ibv) = tr_w(:m,1,ibv)/wi(1,ibv)
        call check_wc(tr_wc(1,1,ibv))


        ! sweep down
        do i=2,n
          if ( f(i,ibv) < 0.d0 ) then
            flux_dt = - f(i,ibv)*dts
            tr_w(:m,i,ibv) = tr_w(:m,i,ibv) + tr_wc(:m,i-1,ibv)*flux_dt
            wi(i,ibv) = wi(i,ibv) + flux_dt
            if ( wi(i,ibv) > 0.d0 )
     &           tr_wc(:m,i,ibv) = tr_w(:m,i,ibv)/wi(i,ibv)
            call check_wc(tr_wc(1,i,ibv))
            tr_w(:m,i-1,ibv) = tr_w(:m,i-1,ibv)
     &           - tr_wc(:m,i-1,ibv)*flux_dt
          endif
        enddo


        ! sweep up
        do i=n-1,1,-1
          if ( f(i+1,ibv) > 0.d0 ) then
            flux_dt = f(i+1,ibv)*dts
            tr_w(:m,i,ibv) = tr_w(:m,i,ibv) + tr_wc(:m,i+1,ibv)*flux_dt
            wi(i,ibv) = wi(i,ibv) + flux_dt
            if ( wi(i,ibv) > 0.d0 )
     &           tr_wc(:m,i,ibv) = tr_w(:m,i,ibv)/wi(i,ibv)
            call check_wc(tr_wc(1,i,ibv))
            tr_w(:m,i+1,ibv) = tr_w(:m,i+1,ibv)
     &           - tr_wc(:m,i+1,ibv)*flux_dt
          endif
        enddo
      enddo

ccc evap from bare soil
      if ( process_bare .and. evapb >= 0.d0 ) then
        evap_tmp(:m) = evapb*(1.d0-fr_snow(1))
#ifdef TRACERS_SPECIAL_O18
        if ( evap_tmp(1)*dts < wi(1,1) .and. tp(1,1) > 0.d0 ) then
          do mc=1,m
ccc tr_name - loop over string comparisons deep inside the nested
ccc loops...
            evap_tmp(mc) = evap_tmp(mc) * fracvl( tp(1,1),tr_name(mc) )
          enddo
        endif
#endif
        tr_w(:m,1,1) = tr_w(:m,1,1) - evap_tmp(:m)*tr_wc(:m,1,1)*dts
        tr_evap(:m,1) = tr_evap(:m,1) + evap_tmp(:m)*tr_wc(:m,1,1)
      endif


ccc runoff
      do ibv=i_bare,i_vege
        tr_rnff(:m,ibv) = tr_rnff(:m,ibv) + tr_wc(:m,1,ibv)*rnf(ibv)
        tr_w(:m,1,ibv) = tr_w(:m,1,ibv) - tr_wc(:m,1,ibv)*rnf(ibv)*dts
        do i=1,n
          tr_rnff(:m,ibv) = tr_rnff(:m,ibv)
     &         + tr_wc(:m,i,ibv)*rnff(i,ibv)
          tr_w(:m,i,ibv) = tr_w(:m,i,ibv)
     &         - tr_wc(:m,i,ibv)*rnff(i,ibv)*dts
        enddo
      enddo


ccc !!! include surface runoff !!!

ccc transpiration
      if ( process_vege ) then
        do i=1,n
          flux_dt = fd*(1.-fr_snow(2)*fm)*evapdl(i,2)*dts
          tr_evap(:m,2) = tr_evap(:m,2) + tr_wc(:m,i,2)*flux_dt/dts
          tr_w(:m,i,2) = tr_w(:m,i,2) - tr_wc(:m,i,2)*flux_dt
        enddo
      endif


cddd      print *, 'end ghy_tracers'
cddd      print *, 'trpr, trpr_dt ', trpr(1), trpr(1)*dts
cddd      print *, 'tr_w 1 ', tr_w(1,:,1)
cddd      print *, 'tr_w 2 ', tr_w(1,:,2)
cddd      print *, 'tr_wsn 1 ', tr_wsn(1,:,1)
cddd      print *, 'tr_wsn 2 ', tr_wsn(1,:,2)
cddd      print *, 'evap ', tr_evap(1,:)*dts
cddd      print *, 'runoff ', tr_rnff(1,:)*dts

      do mc=1,m
         do ibv=i_bare,i_vege
            do i=1,n
               if ( tr_w(mc,i,ibv) < 0.d0 ) then
                  if ( tr_w(mc,i,ibv) < -tol(mc) ) then
                     print *,'GHY:tr_w<0 ', tr_w(mc,i,ibv), ibv, i
                     call stop_model("GHY:tr_w<0",255)
                  endif
                  tr_w(mc,i,ibv) = 0.d0
               endif
            enddo
            do i=1,nlsn
               if ( tr_wsn(mc,i,ibv) < 0.d0 ) then
                  if ( tr_wsn(mc,i,ibv) < -tol(mc) ) then
                     print *,'GHY:tr_wsn<0 ', tr_wsn(mc,i,ibv),
     &                    ibv, i, flmlt(ibv), fr_snow(ibv)
                  endif
                  tr_wsn(mc,i,ibv) = 0.d0
               endif
            enddo
            err = trpr(mc)*dts + trdd(mc)*dts
     &           + sum( tr_w_o(mc,:,ibv) ) + sum( tr_wsn_o(mc,:,ibv) )
     &           - sum( tr_w(mc,:,ibv) ) - sum( tr_wsn(mc,:,ibv) )
     &           - tr_evap(mc,ibv)*dts - tr_rnff(mc,ibv)*dts
                                !print *,'err ', err
            if ( abs( err ) > tol(mc) ) then
               write(0,*)
     $              'ghy tracers not conserved at',ijdebug,' err=',err
               write(99,*)
     $              'ghy tracers not conserved at',ijdebug,' err=',err
     $              ,"mc= ",mc, "ibv= ",ibv
     $              ,"value =",trpr(mc)*dts
     &              + sum(tr_w_o(mc,:,ibv) ) + sum( tr_wsn_o(mc,:,ibv))
               write(99,*) "tr_w_o",tr_w_o(mc,:,ibv)
     $              ,"tr_wsn_o",tr_wsn_o(mc,:,ibv)
     $              ,"tr_w",tr_w(mc,:,ibv)
     $              ,"tr_wsn",tr_wsn(mc,:,ibv)
     $              ,"trpr(mc)*dts",trpr(mc)*dts
     $              ,"trdd(mc)*dts",trdd(mc)*dts
     $              ,"tr_evap*tds",tr_evap(mc,ibv)*dts
     $              ,"tr_rnff",tr_rnff(mc,ibv)*dts
!!!            call stop_model("ghy tracers not conserved",255)
            endif
         enddo
      enddo


!!! hack
!!! get rid of possible garbage
      do ibv=i_bare,i_vege
        do i=1,nlsn
          if ( wsni(i,ibv) < 1.d-14 ) then
            wsni(i,ibv) = 0.d0
            tr_wsn(:m,i,ibv) = 0.d0
          endif
        enddo
      enddo


      end subroutine ghy_tracers
#endif


      subroutine check_water( flag )
      integer flag
      real*8 total_water(2), error_water
      real*8 old_total_water(2), old_fr_snow(2) ! save
      COMMON /check_water_tp/ old_total_water, old_fr_snow
!$OMP  THREADPRIVATE (/check_water_tp/)
      integer k, ibv

      total_water(1) = 0.
      total_water(2) = w(0,2)

      do ibv=1,2
        do k=1,nsn(ibv)
          total_water(ibv) = total_water(ibv) + wsn(k,ibv)*fr_snow(ibv)
        enddo
        do k=1,n
          total_water(ibv) = total_water(ibv) + w(k,ibv)
        enddo
      enddo ! ibv

      if ( flag == 0 ) then
        old_total_water(:) = total_water(:)
        old_fr_snow(:) = fr_snow(:)
        return
      endif

      ! bare soil
      if ( process_bare ) then
        error_water = (total_water(1) - old_total_water(1)) / dts
     $       - pr + evap_tot(1) + sum(rnff(1:n,1)) + rnf(1)
        if (abs(error_water)>1.d-15)
     $       write(99,*)'bare',ijdebug,error_water
        if ( abs( error_water ) > 1.d-13 ) then
c    &       call stop_model('GHY: water conservation problem',255)
           write(0,*)'evap_tot(1)',ijdebug,evap_tot(1)
        endif
      endif

      ! vegetated soil
      if ( process_vege ) then
        error_water = (total_water(2) - old_total_water(2)) / dts
     &       - pr + evap_tot(2) + sum(rnff(1:n,2)) + rnf(2)
        if (abs(error_water)>1.d-15)
     $       write(99,*)'vege',ijdebug,error_water
        if ( abs( error_water ) > 1.d-10 ) then
          write(0,*)'evap_tot(2)',ijdebug,evap_tot(2)
c         call stop_model(
c    &       'GHY: water conservation problem in veg. soil',255)
        endif
      endif

      ghy_debug%water(:) = total_water(:)

      end subroutine check_water


      subroutine check_energy( flag )
      integer flag
      real*8 total_energy(2), error_energy
      real*8 old_total_energy(2), old_fr_snow(2) !save
      COMMON /check_energy_tp/ old_total_energy, old_fr_snow
!$OMP  THREADPRIVATE (/check_energy_tp/)
      integer k, ibv

      total_energy(1) = 0.
      total_energy(2) = ht(0,2)

      do ibv=1,2
        do k=1,nsn(ibv)
          total_energy(ibv) = total_energy(ibv)+hsn(k,ibv)*fr_snow(ibv)
        enddo
        do k=1,n
          total_energy(ibv) = total_energy(ibv) + ht(k,ibv)
        enddo
      enddo ! ibv

      if ( flag == 0 ) then
        old_total_energy(:) = total_energy(:)
        old_fr_snow(:) = fr_snow(:)
        return
      endif

      ! bare soil
      if ( process_bare ) then
        error_energy = (total_energy(1) - old_total_energy(1)) / dts
     $       - htpr + elh*evap_tot(1)
     &       + shw*( sum( rnff(1:n,1)*max(tp(1:n,1),0.d0) )
     &       + rnf(1)*max(tp(1,1),0.d0) )
     &       - srht - trht + thrm_tot(1) + snsh_tot(1)

        if ( abs( error_energy ) > 1.d-5 ) then
          write(0,*)'GHY:bare soil error_energy',ijdebug,error_energy
          !call stop_model('GHY: energy conservation problem',255)
        endif
      endif

      ! vegetated soil
      !!! if ( fr_snow(2) > 0.d0 .or. old_fr_snow(2) > 0.d0 ) return
      if ( process_vege ) then
        error_energy = (total_energy(2) - old_total_energy(2)) / dts
     $       - htpr + elh*evap_tot(2)
     &       + shw*( sum( rnff(1:n,2)*max(tp(1:n,2),0.d0) )
     &       + rnf(2)*max(tp(1,2),0.d0) )
     &       - srht - trht + thrm_tot(2) + snsh_tot(2)

        if ( abs( error_energy ) > 1.d-5) then
          write(0,*)'GHY:veg soil error_energy',ijdebug,error_energy
          !call stop_model('GHY: energy cons problem in veg. soil',255)
        endif
      endif

      ghy_debug%energy(:) = total_energy(:)

      end subroutine check_energy


      subroutine restrict_snow (wsn_max)
!@sum remove the snow in excess of WSN_MAX and dump its water into
!@+   runoff, keeping associated heat in the soil
      real*8 :: wsn_max,wsn_tot,d_wsn,eta,dw(2:3),dh(2:3)
      integer :: ibv

      do ibv=i_bare,i_vege

        if ( fr_snow(ibv) <= 0.d0 ) cycle
        wsn_tot = sum( wsn(1:nsn(ibv),ibv) )
        if ( wsn_tot <= WSN_MAX ) cycle

         ! check if snow structure ok for thick snow
        if ( nsn(ibv) < 3)
     &       call stop_model("remove_extra_snow: nsn<3",255)

        d_wsn = wsn_tot - WSN_MAX

        ! do not remove too much at a time (for now 24mm / day)
        d_wsn = min( d_wsn, 1.d-3*dts/3600.d0 )

        !print *,"restricting snow: wsn_tot = ", wsn_tot,ijdebug
        !print *,"before:",hsn(1:3, ibv),fr_snow(ibv) ,ht(1,ibv)
              ! fraction of snow to be removed:
        eta = d_wsn/sum(wsn(2:3,ibv))
        dw(2:3) = eta*wsn(2:3, ibv)
        dh(2:3) = eta*hsn(2:3, ibv)
        !print *,"dw,dh:",dw,dh
        wsn(2:3, ibv)  = wsn(2:3, ibv) - dw(2:3)
        hsn(2:3, ibv)  = hsn(2:3, ibv) - dh(2:3)
        dzsn(2:3, ibv) = (1.d0-eta)*dzsn(2:3, ibv)

        rnf(ibv) = rnf(ibv)   + sum(dw(2:3))*fr_snow(ibv)/dts
        ht(1,ibv) = ht(1,ibv) + sum(dh(2:3))*fr_snow(ibv)
     &       - sum(dw(2:3))*fr_snow(ibv)*max(tp(1,ibv),0.d0)*shw
        !print *,"after:",hsn(1:3, ibv),fr_snow(ibv) ,ht(1,ibv)

#ifdef TRACERS_WATER
        tr_rnff(1:ntg,ibv) = tr_rnff(1:ntg,ibv) +
     &       eta*(tr_wsn(1:ntg,2,ibv)+tr_wsn(1:ntg,3,ibv))
        tr_wsn(1:ntg,2:3,ibv) = (1.d0-eta)*tr_wsn(1:ntg,2:3,ibv)
#endif

      end do

      end subroutine restrict_snow


      end module sle001



ccccccccccc  notes ccccccccccccccccc

ccc just in case, the thickness of layers is:
c     n, dz =  1,  9.99999642372131348E-2
c     n, dz =  2,  0.17254400253295898
c     n, dz =  3,  0.2977144718170166
c     n, dz =  4,  0.51368874311447144
c     n, dz =  5,  0.88633960485458374
c     n, dz =  6,  1.5293264389038086
