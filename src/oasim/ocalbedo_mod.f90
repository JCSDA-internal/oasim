! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ocalbedo_mod

! oasim
use oasim_constants_mod

implicit none
private
public ocalbedo

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------


subroutine ocalbedo(rad, wfac, theta, wspd, rod, ros)

! Computes ocean surface albedo from solar zenith angle (sunz) and wind speed (ws, m/s).
! Albedo is provided as direct (albd) and diffuse (albs).

! Arguments
real(kind=kind_real), intent(in)  :: rad
real(kind=kind_real), intent(in)  :: wfac(nlt)
real(kind=kind_real), intent(in)  :: theta
real(kind=kind_real), intent(in)  :: wspd
real(kind=kind_real), intent(out) :: rod(nlt)
real(kind=kind_real), intent(out) :: ros(nlt)

! Derive surface reflectance as a function of sunz and ws for OASIM Bands
call sfcrfl(rad, wfac, theta, wspd, rod, ros)

end subroutine ocalbedo

! --------------------------------------------------------------------------------------------------

subroutine sfcrfl(rad, wfac, theta, wspd, rod, ros)

! Computes surface reflectance for direct (rod) and diffuse (ros) components separately, as a
! function of theta, wind speed or stress. Includes spectral dependence of foam reflectance derived
! from Frouin et al., 1996 (JGR)

! Arguments
real(kind=kind_real), intent(in)  :: rad
real(kind=kind_real), intent(in)  :: wfac(nlt)
real(kind=kind_real), intent(in)  :: theta
real(kind=kind_real), intent(in)  :: wspd
real(kind=kind_real), intent(out) :: rod(nlt)
real(kind=kind_real), intent(out) :: ros(nlt)

! Locals
integer :: nl
real(kind=kind_real) :: a, b, cn, rmin, rof, rospd, rosps, rpls, rtheta, rthetar, sinp, sinrmin
real(kind=kind_real) :: sinrpls, sintr, tanp, tanrmin, tanrpls

! Paramters
real(kind=kind_real), parameter :: rn = 1.341_kind_real
real(kind=kind_real), parameter :: roair = 1.2E3_kind_real
real(kind=kind_real), parameter :: a0 = 0.9976_kind_real
real(kind=kind_real), parameter :: a1 = 0.2194_kind_real
real(kind=kind_real), parameter :: a2 = 5.554E-2_kind_real
real(kind=kind_real), parameter :: a3 = 6.7E-3_kind_real
real(kind=kind_real), parameter :: b0 = 5.026_kind_real
real(kind=kind_real), parameter :: b1 = -0.01138_kind_real
real(kind=kind_real), parameter :: b2 = 9.552E-6_kind_real
real(kind=kind_real), parameter :: b3 = -2.698E-9_kind_real

!  Foam and diffuse reflectance
if (wspd .gt. 4.0_kind_real)then
  if (wspd .le. 7.0_kind_real)then
    cn = 6.2E-4_kind_real + 1.56E-3_kind_real/wspd
    rof = roair*cn*2.2E-5_kind_real*wspd*wspd - 4.0E-4_kind_real
  else
    cn = 0.49E-3 + 0.065E-3*wspd
    rof = (roair*cn*4.5E-5_kind_real - 4.0E-5_kind_real)*wspd*wspd
  endif
  rosps = 0.057_kind_real
else
  rof = 0.0_kind_real
  rosps = 0.066_kind_real
endif

! Direct
! Fresnel reflectance for theta < 40, wspd < 2 m/s
if (theta .lt. 40.0_kind_real .or. wspd .lt. 2.0_kind_real)then
  if (theta .eq. 0.0_kind_real)then
    rospd = 0.0211_kind_real
  else
    rtheta = theta/rad
    sintr = sin(rtheta)/rn
    rthetar = asin(sintr)
    rmin = rtheta - rthetar
    rpls = rtheta + rthetar
    sinrmin = sin(rmin)
    sinrpls = sin(rpls)
    tanrmin = tan(rmin)
    tanrpls = tan(rpls)
    sinp = (sinrmin*sinrmin)/(sinrpls*sinrpls)
    tanp = (tanrmin*tanrmin)/(tanrpls*tanrpls)
    rospd = 0.5_kind_real*(sinp + tanp)
  endif
else
  ! Empirical fit otherwise
  a = 0.0253_kind_real
  b = -7.14E-4_kind_real*wspd + 0.0618_kind_real
  rospd = a*exp(b*(theta-40.0_kind_real))
endif

! Reflectance totals
do nl = 1,nlt
  rod(nl) = rospd + rof*wfac(nl)
  ros(nl) = rosps + rof*wfac(nl)
enddo

end subroutine sfcrfl

! --------------------------------------------------------------------------------------------------

end module ocalbedo_mod
