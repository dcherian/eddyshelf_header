!-*-fortran-*-
      SUBROUTINE ana_hmixcoef (ng, tile, model)
!
!! svn $Id: ana_hmixcoef.h 645 2013-01-22 23:21:54Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine rescales horizontal mixing coefficients according      !
!  to the grid size.  Also,  if applicable,  increases horizontal      !
!  in sponge areas.                                                    !
!                                                                      !
!  WARNING:   All biharmonic coefficients are assumed to have the      !
!             square root taken and have  m^2 s^-1/2 units.  This      !
!             will allow multiplying the  biharmonic  coefficient      !
!             to harmonic operator.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"

      CALL ana_hmixcoef_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
#ifdef SOLVE3D
# ifdef TS_DIF2
     &                        MIXING(ng) % diff2,                       &
# endif
# ifdef TS_DIF4
     &                        MIXING(ng) % diff4,                       &
# endif
#endif
#ifdef UV_VIS2
     &                        MIXING(ng) % visc2_p,                     &
     &                        MIXING(ng) % visc2_r,                     &
#endif
#ifdef UV_VIS4
     &                        MIXING(ng) % visc4_p,                     &
     &                        MIXING(ng) % visc4_r,                     &
#endif
     &                        GRID(ng) % grdscl,                        &
     &                        GRID(ng) % xr,                            &
     &                        GRID(ng) % yr)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 8)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_hmixcoef
!
!***********************************************************************
      SUBROUTINE ana_hmixcoef_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
#ifdef SOLVE3D
# ifdef TS_DIF2
     &                              diff2,                              &
# endif
# ifdef TS_DIF4
     &                              diff4,                              &
# endif
#endif
#ifdef UV_VIS2
     &                              visc2_p,                            &
     &                              visc2_r,                            &
#endif
#ifdef UV_VIS4
     &                              visc4_p,                            &
     &                              visc4_r,                            &
#endif
     &                              grdscl, xr, yr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_grid
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: grdscl(LBi:,LBj:)
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
!      real(r8), intent(in) :: pm(LBi:,LBj:)
# ifdef SOLVE3D
#  ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:,LBj:,:)
#  endif
#  ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:,LBj:,:)
#  endif
# endif
# ifdef UV_VIS2
      real(r8), intent(inout) :: visc2_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc2_r(LBi:,LBj:)
# endif
# ifdef UV_VIS4
      real(r8), intent(inout) :: visc4_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc4_r(LBi:,LBj:)
# endif
#else
      real(r8), intent(in) :: grdscl(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
#  ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:UBi,LBj:UBj,NT(ng))
#  endif
# endif
# ifdef UV_VIS2
      real(r8), intent(inout) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: visc2_r(LBi:UBi,LBj:UBj)
# endif
# ifdef UV_VIS4
      real(r8), intent(inout) :: visc4_p(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: visc4_r(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: Iwrk, i, j, itrc
      real(r8) :: cff, cff1, cff2, width, eps, XX, YY
      real(r8) :: maxdiff,maxvisc,dx
#ifdef WC13
      real(r8) :: cff_t, cff_s, cff1_t, cff2_t, cff1_s, cff2_s
#endif
#include "set_bounds.h"

#ifdef VISC_GRID
!
!-----------------------------------------------------------------------
!  Scale horizontal viscosity according to the grid size.
!-----------------------------------------------------------------------
!
!! WARNING:  This section is generic for all applications. Please do not
!!           change the code below.
!!
# ifdef UV_VIS2
      cff=visc2(ng)/grdmax(ng)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          visc2_r(i,j)=cff*grdscl(i,j)
        END DO
      END DO
      cff=0.25_r8*cff
      DO j=Jstr,JendR
        DO i=Istr,IendR
          visc2_p(i,j)=cff*(grdscl(i,j  )+grdscl(i-1,j  )+              &
     &                      grdscl(i,j-1)+grdscl(i-1,j-1))
        END DO
      END DO
# endif
# ifdef UV_VIS4
      cff=visc4(ng)/(grdmax(ng)**3)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          visc4_r(i,j)=cff*grdscl(i,j)**3
        END DO
      END DO
      cff=0.25_r8*cff
      DO j=Jstr,JendR
        DO i=Istr,IendR
          visc4_p(i,j)=cff*(grdscl(i,j  )**3+grdscl(i-1,j  )**3+        &
     &                      grdscl(i,j-1)**3+grdscl(i-1,j-1)**3)
        END DO
      END DO
# endif
#endif
#ifdef DIFF_GRID
!
!-----------------------------------------------------------------------
!  Scale horizontal diffusion according to the grid size.
!-----------------------------------------------------------------------
!
!! WARNING:  This section is generic for all applications. Please do not
!!           change the code below.
!!
# ifdef TS_DIF2
      DO itrc=1,NT(ng)
        cff=tnu2(itrc,ng)/grdmax(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            diff2(i,j,itrc)=cff*grdscl(i,j)
          END DO
        END DO
      END DO
# endif
# ifdef TS_DIF4
      DO itrc=1,NT(ng)
        cff=tnu4(itrc,ng)/(grdmax(ng)**3)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            diff4(i,j,itrc)=cff*grdscl(i,j)**3
          END DO
        END DO
      END DO
# endif
#endif
#ifdef SPONGE
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing in the sponge areas.
!-----------------------------------------------------------------------
!
!! User modifiable section.  Please specify the appropiate sponge area
!! by increasing its horizontal mixing coefficients.
!!
#if defined EDDY_SPONGE
!     Convert width (in points)  of sponge to metres first
! assuming square grid
      dx = 1000

      width=USER(1)*dx
      maxvisc=USER(2)
      maxdiff=USER(3)
      XX = xl(ng)
      YY = el(ng)
      eps=1E-9
      
!      print *,'dx = ', dx,' | sponge width = ', width
!      print *,'maxval(xr) = ', XX
!      print *,'maxval(yr) = ', YY

#if defined UV_VIS4
      maxvisc=sqrt(maxvisc)
#endif
#if defined TS_DIF4
      maxdiff=sqrt(maxdiff)
#endif

! cff contains value between 0 and 1 in sponge region
! multiply by maxvisc / maxdiff later and assign to appropriate array
      IF (width .gt. eps) THEN
      DO j=JstrR,JendR
        DO i=IstrR,IendR
           cff = 0.0d0
#if defined EDDY_SPONGE_SOUTH
            IF (yr(i,j) .lt. width) THEN
              cff=(width-yr(i,j))/width
            END IF
#endif ! EDDY_SPONGE_SOUTH

#if defined EDDY_SPONGE_NORTH
            IF (yr(i,j) .gt. YY-width) THEN
              cff=(yr(i,j)+width-YY)/width
            END IF
#endif

#if defined EDDY_SPONGE_EAST
            IF (xr(i,j) .gt. (XX-width)) THEN
              cff=(xr(i,j)+width-XX)/width
            END IF
#endif

#if defined EDDY_SPONGE_WEST
            IF (xr(i,j) .lt. width) THEN
              cff=(width-xr(i,j))/width
            END IF
#endif
!
! deal with corners - needed with parallel
!
#if defined EDDY_SPONGE_NORTH && defined EDDY_SPONGE_EAST
            IF (xr(i,j) .gt. (XX-width)) THEN
               IF (yr(i,j) .gt. (YY-width)) THEN
                  cff = (yr(i,j)+width-YY)/width +              &
     &                  (xr(i,j)+width-XX)/width
               END IF
            END IF
#endif
#if defined EDDY_SPONGE_SOUTH && defined EDDY_SPONGE_EAST
            IF (xr(i,j) .gt. (XX-width)) THEN
               IF (yr(i,j) .lt. (width)) THEN
                  cff = (width-yr(i,j))/width +              &
     &                  (xr(i,j)+width-XX)/width
               END IF
            END IF
#endif
#if defined EDDY_SPONGE_NORTH && defined EDDY_SPONGE_WEST
            IF (xr(i,j) .lt. (width)) THEN
               IF (yr(i,j) .gt. (YY-width)) THEN
                  cff = (width+yr(i,j)-YY)/width +              &
     &                  (width-xr(i,j))/width
               END IF
            END IF
#endif
#if defined EDDY_SPONGE_SOUTH && defined EDDY_SPONGE_WEST
            IF (xr(i,j) .lt. (width)) THEN
               IF (yr(i,j) .lt. (width)) THEN
                  cff = (width-yr(i,j))/width +              &
     &                  (width-xr(i,j))/width
               END IF
            END IF
#endif

! assign to viscosity and diffusivity arrays
#if defined UV_VIS2
            visc2_r(i,j) = cff*maxvisc + visc2(ng)
            visc2_p(i,j) = cff*maxvisc + visc2(ng)
            IF (visc2_r(i,j) .gt. maxvisc) THEN
               visc2_r(i,j) = maxvisc
               visc2_p(i,j) = maxvisc
            END IF
#endif
#if defined UV_VIS4
            visc4_r(i,j) = cff*maxvisc + visc4(ng)
            visc4_p(i,j) = cff*maxvisc + visc4(ng)
            IF (visc4_r(i,j) .gt. maxvisc) THEN
               visc4_r(i,j) = maxvisc
               visc4_p(i,j) = maxvisc
            END IF
#endif
#if defined TS_DIF2
            diff2(i,j,itemp) = cff*maxdiff + tnu2(itemp,ng)
            diff2(i,j,isalt) = cff*maxdiff + tnu2(isalt,ng)
            IF (diff2(i,j,itemp) .gt. maxdiff) THEN
               diff2(i,j,itemp) = maxdiff
               diff2(i,j,isalt) = maxdiff
            END IF
#endif
#if defined TS_DIF4
            diff4(i,j,itemp) = cff*maxdiff + tnu4(itemp,ng)
            diff4(i,j,isalt) = cff*maxdiff + tnu4(isalt,ng)
            IF (diff4(i,j,itemp) .gt. maxdiff) THEN
               diff4(i,j,itemp) = maxdiff
               diff4(i,j,isalt) = maxdiff
            END IF
#endif

         END DO
         END DO
      END IF
! diagnostic output
!#if defined UV_VIS2
!      print *,'tile = ',tile,'maxval(visc2_r) = ', maxval(visc2_r)
!      print *,'tile = ',tile,'maxval(visc2_p) = ', maxval(visc2_p)
!#endif
!#if defined UV_VIS4
!      print *,'tile = ',tile,'maxval(visc4_r) = ', maxval(visc4_r)
!      print *,'tile = ',tile,'maxval(visc4_p) = ', maxval(visc4_p)
!#endif
!#if defined TS_DIF2
!      print *,'tile = ',tile,'maxval(diff2) = ', maxval(diff2)
!#endif
!#if defined TS_DIF4
!      print *,'tile = ',tile,'maxval(diff4) = ', maxval(diff4)
!#endif

#endif ! EDDY_SPONGE

# if defined ADRIA02
!
!  Adriatic Sea southern sponge areas.
!
      fac=4.0_r8
#  if defined UV_VIS2
      DO i=IstrR,IendR
        DO j=JstrR,MIN(6,JendR)
          cff=visc2(ng)+REAL(6-j,r8)*(fac*visc2(ng)-visc2(ng))/6.0_r8
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
        DO j=MAX(JstrR,7),JendR
          visc2_r(i,j)=0.0_r8
          visc2_p(i,j)=0.0_r8
        END DO
      END DO
#  endif
#  if defined TS_DIF2
      DO i=IstrR,IendR
        DO j=JstrR,MIN(6,JendR)
          cff1=tnu2(itemp,ng)+                                          &
     &         REAL(6-j,r8)*(fac*tnu2(itemp,ng)-tnu2(itemp,ng))/6.0_r8
          cff2=tnu2(isalt,ng)+                                          &
     &         REAL(6-j,r8)*(fac*tnu2(isalt,ng)-tnu2(isalt,ng))/6.0_r8
          diff2(i,j,itemp)=cff1
          diff2(i,j,isalt)=cff2
        END DO
        DO j=MAX(JstrR,7),JendR
          diff2(i,j,itemp)=0.0_r8
          diff2(i,j,isalt)=0.0_r8
        END DO
      END DO
#  endif

# elif defined WC13
!
!  US West Coast sponge areas.
!
      Iwrk=INT(user(1))  ! same for sponge and nudging layers

#  if defined UV_VIS2
!
!  Momentum sponge regions:  sponge viscosities as in Marchesiello
!  et al 2003.
!
      cff1=visc2(ng)
      cff2=100.0_r8
!
!  Southern edge.
!
      DO j=JstrR,MIN(Iwrk,JendR)
        cff=cff1+REAL(Iwrk-j,r8)*(cff2-cff1)/REAL(Iwrk,r8)
        DO i=IstrR,IendR
          visc2_r(i,j)=MAX(MIN(cff,cff2),cff1)
          visc2_p(i,j)=MAX(MIN(cff,cff2),cff1)
        END DO
      END DO
!
!  Northern edge.
!
      DO j=MAX(JstrR,Mm(ng)+1-Iwrk),JendR
        cff=cff2-REAL(Mm(ng)+1-j,r8)*(cff2-cff1)/REAL(Iwrk,r8)
        DO i=IstrR,IendR
          visc2_r(i,j)=MAX(MIN(cff,cff2),cff1)
          visc2_p(i,j)=MAX(MIN(cff,cff2),cff1)
        END DO
      END DO
!
!  Western edge.
!
      DO i=IstrR,MIN(Iwrk,IendR)
        DO j=MAX(JstrR,i),MIN(Mm(ng)+1-i,JendR)
          cff=cff1+REAL(Iwrk-i,r8)*(cff2-cff1)/REAL(Iwrk,r8)
          visc2_r(i,j)=MAX(MIN(cff,cff2),cff1)
          visc2_p(i,j)=MAX(MIN(cff,cff2),cff1)
        END DO
      END DO
#  endif
#  if defined TS_DIF2
!
!  Tracer sponge regions: sponge diffusivities as in Marchesiello
!  et al 2003.
!
      cff1_t=tnu2(itemp,ng)
      cff1_s=tnu2(isalt,ng)
      cff2_t=50.0_r8
      cff2_s=50.0_r8
!
!  Southern edge.
!
      DO j=JstrR,MIN(Iwrk,JendR)
        cff_t=cff1_t+REAL(Iwrk-j,r8)*(cff2_t-cff1_t)/REAL(Iwrk,r8)
        cff_s=cff1_s+REAL(Iwrk-j,r8)*(cff2_s-cff1_s)/REAL(Iwrk,r8)
        DO i=IstrR,IendR
          diff2(i,j,itemp)=MAX(MIN(cff_t,cff2_t),cff1_t)
          diff2(i,j,isalt)=MAX(MIN(cff_s,cff2_s),cff1_s)
        END DO
      END DO
!
!  Northern edge.
!
      DO j=MAX(JstrR,Mm(ng)+1-Iwrk),JendR
        cff_t=cff2_t-REAL(Mm(ng)+1-j,r8)*(cff2_t-cff1_t)/REAL(Iwrk,r8)
        cff_s=cff2_s-REAL(Mm(ng)+1-j,r8)*(cff2_s-cff1_s)/REAL(Iwrk,r8)
        DO i=IstrR,IendR
          diff2(i,j,itemp)=MAX(MIN(cff_t,cff2_t),cff1_t)
          diff2(i,j,isalt)=MAX(MIN(cff_s,cff2_s),cff1_s)
        END DO
      END DO
!
!  Western edge.
!
      DO i=IstrR,MIN(Iwrk,IendR)
        DO j=MAX(JstrR,i),MIN(Mm(ng)+1-i,JendR)
          cff_t=cff1_t+REAL(Iwrk-i,r8)*(cff2_t-cff1_t)/REAL(Iwrk,r8)
          cff_s=cff1_s+REAL(Iwrk-i,r8)*(cff2_s-cff1_s)/REAL(Iwrk,r8)
          diff2(i,j,itemp)=MAX(MIN(cff_t,cff2_t),cff1_t)
          diff2(i,j,isalt)=MAX(MIN(cff_s,cff2_s),cff1_s)
        END DO
      END DO
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
!! WARNING:  This section is generic for all applications. Please do not
!!           change the code below.
!!
#ifdef UV_VIS2
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          visc2_r)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          visc2_p)
      END IF
#endif
#ifdef UV_VIS4
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          visc4_r)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          visc4_p)
      END IF
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            diff2(:,:,itrc))
        END DO
      END IF
# endif
# ifdef TS_DIF4
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            diff4(:,:,itrc))
        END DO
      END IF
# endif
#endif

#ifdef DISTRIBUTE
# ifdef UV_VIS2
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    visc2_r, visc2_p)
# endif
# ifdef UV_VIS4
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    visc4_r, visc4_p)
# endif
# ifdef SOLVE3D
#  ifdef TS_DIF2
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    diff2)
#  endif
#  ifdef TS_DIF4
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    diff4)
#  endif
# endif
#endif

      RETURN
      END SUBROUTINE ana_hmixcoef_tile
