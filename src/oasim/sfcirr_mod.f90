! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module sfcirr_mod

! oasim
use oasim_constants_mod
use slingo_mod, only: slingo
use clrtrans_mod, only: clrtrans

implicit none
private
public sfcirr

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine sfcirr(lam, fobar, thray, oza, awv, ao, aco2, asl, bsl, csl, dsl, esl, fsl, ica, &
                  daycor, cosunz, pres, ws, ozone, wvapor, relhum, ta, wa, asym, am, vi, cov, &
                  clwp, re, ed, es)

! Calls clrtrans to get cloud-free transmittance and slingo to
! get cloudy transmittance, then computes total irradiance in
! W/m2/(variable)nm weighted by the cloudiness.

! Tdclr is spectral clear sky direct transmittance
! Tsclr is spectral clear sky diffuse transmittance
! Tdcld is spectral cloudy direct transmittance
! Tscld is spectral cloudy diffuse transmittance

! Arguments
integer,              intent(in)    :: lam(nlt)
real(kind=kind_real), intent(in)    :: Fobar(nlt)
real(kind=kind_real), intent(in)    :: thray(nlt)
real(kind=kind_real), intent(in)    :: oza(nlt)
real(kind=kind_real), intent(in)    :: awv(nlt)
real(kind=kind_real), intent(in)    :: ao(nlt)
real(kind=kind_real), intent(in)    :: aco2(nlt)
real(kind=kind_real), intent(in)    :: asl(ncld)
real(kind=kind_real), intent(in)    :: bsl(ncld)
real(kind=kind_real), intent(in)    :: csl(ncld)
real(kind=kind_real), intent(in)    :: dsl(ncld)
real(kind=kind_real), intent(in)    :: esl(ncld)
real(kind=kind_real), intent(in)    :: fsl(ncld)
integer,              intent(in)    :: ica(nlt)
real(kind=kind_real), intent(in)    :: daycor
real(kind=kind_real), intent(in)    :: cosunz
real(kind=kind_real), intent(in)    :: pres
real(kind=kind_real), intent(in)    :: ws
real(kind=kind_real), intent(in)    :: ozone
real(kind=kind_real), intent(in)    :: wvapor
real(kind=kind_real), intent(inout) :: relhum
real(kind=kind_real), intent(inout) :: ta(nlt)
real(kind=kind_real), intent(inout) :: wa(nlt)
real(kind=kind_real), intent(in)    :: asym(nlt)
real(kind=kind_real), intent(in)    :: am
real(kind=kind_real), intent(in)    :: vi
real(kind=kind_real), intent(in)    :: cov
real(kind=kind_real), intent(in)    :: clwp
real(kind=kind_real), intent(in)    :: re
real(kind=kind_real), intent(out)   :: ed(nlt)
real(kind=kind_real), intent(out)   :: es(nlt)

! Locals
integer :: nl
real(kind=kind_real) :: Tgas(nlt), Tdclr(nlt), Tsclr(nlt), Tdcld(ncld), Tscld(ncld)
real(kind=kind_real) :: Edclr(nlt), Esclr(nlt), Edcld(nlt), Escld(nlt)
real(kind=kind_real) :: sunz, rtmp, rmu0, rm, otmp, rmo, rmp, to, oarg, ag, gtmp, gtmp2, garg, wtmp
real(kind=kind_real) :: wtmp2, warg, foinc, edir, edif, ccov1

! Paramters
real(kind=kind_real), parameter :: ozfac1 = 44.0_kind_real/6370.0_kind_real
real(kind=kind_real), parameter :: ozfac2 = 1.0_kind_real + 22.0_kind_real/6370.0_kind_real
real(kind=kind_real), parameter :: p0 = 1013.25_kind_real

ed = 0.0
es = 0.0

if (pres .lt. 0.0_kind_real .or. ws .lt. 0.0_kind_real .or. relhum .lt. 0.0_kind_real .or. &
    ozone .lt. 0.0_kind_real .or. wvapor .lt. 0.0_kind_real) then
  return
endif

!  Compute atmospheric path lengths (air mass); not pressure-corrected
sunz = acos(cosunz)*57.29578_kind_real
rtmp = (93.885_kind_real-sunz)**(-1.253_kind_real)
rmu0 = cosunz+0.15_kind_real*rtmp
rm = 1.0_kind_real/rmu0
otmp = (cosunz*cosunz+ozfac1)**0.5_kind_real
rmo = ozfac2/otmp

!  Compute pressure-corrected atmospheric path length (air mass)
rmp = pres/p0*rm

! Loop to compute total irradiances at each grid point
! Compute direct and diffuse irradiance for a cloudy and non-cloudy atmosphere
! Account for gaseous absorption

do nl = 1, nlt

  ! Ozone
  to = oza(nl)*ozone*0.001_kind_real   !optical thickness
  oarg = -to*rmo

  ! Oxygen/gases
  ag = ao(nl) + aco2(nl)
  gtmp = (1.0_kind_real + 118.3_kind_real*ag*rmp)**0.45_kind_real
  gtmp2 = -1.41_kind_real*ag*rmp
  garg = gtmp2/gtmp

  ! Water Vapor
  wtmp = (1.0_kind_real+20.07_kind_real*awv(nl)*wvapor*rm)**0.45_kind_real
  wtmp2 = -0.2385_kind_real*awv(nl)*wvapor*rm
  warg = wtmp2/wtmp
  tgas(nl) = exp(oarg+garg+warg)

enddo

! Compute clear sky transmittances
call clrtrans(lam, thray, cosunz, rm, rmp, ws, relhum, am, vi, ta, wa, asym, tdclr, tsclr)

do nl = 1,nlt
  foinc = fobar(nl)*daycor*cosunz
  ! Direct irradiance
  edir = foinc*tgas(nl)*tdclr(nl)
  ! Diffuse irradiance
  edif = foinc*tgas(nl)*tsclr(nl)
  ! Spectral components
  edclr(nl) = edir
  esclr(nl) = edif
enddo    !end clear nl loop

!  Compute cloudy transmittances
call slingo(asl, bsl, csl, dsl, esl, fsl, rmu0, clwp, re, tdcld, tscld)

do nl = 1,nlt
  foinc = fobar(nl)*daycor*cosunz
  ! Direct irradiance
  edir = foinc*tgas(nl)*tdcld(ica(nl))
  ! Diffuse irradiance
  edif = foinc*tgas(nl)*tscld(ica(nl))
  ! Spectral components
  edcld(nl) = edir
  escld(nl) = edif
enddo   !end cloudy nl loop

! Sum clear and cloudy by percent
ccov1 = cov*0.01_kind_real  !convert from percent to fraction
do nl = 1,nlt
  ed(nl) = (1.0_kind_real-ccov1)*edclr(nl) + ccov1*edcld(nl)
  es(nl) = (1.0_kind_real-ccov1)*esclr(nl) + ccov1*escld(nl)
enddo

end subroutine sfcirr

! --------------------------------------------------------------------------------------------------

end module sfcirr_mod
