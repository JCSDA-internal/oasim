! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module setlte_mod

! oasim
use oasim_constants_mod
use lidata_mod

implicit none
private
public setlte

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine setlte(lam, aw, bw, ac, bc, bpic, excdom, exdet, wtoq, wfac, data_directory)

! Arguments
integer,              intent(inout) :: lam(nlt)
real(kind=kind_real), intent(inout) :: aw(nlt)
real(kind=kind_real), intent(inout) :: bw(nlt)
real(kind=kind_real), intent(inout) :: ac(nchl,nlt)
real(kind=kind_real), intent(inout) :: bc(nchl,nlt)
real(kind=kind_real), intent(inout) :: bpic(nlt)
real(kind=kind_real), intent(inout) :: excdom(nlt)
real(kind=kind_real), intent(inout) :: exdet(nlt)
real(kind=kind_real), intent(inout) :: wtoq(nlt)
real(kind=kind_real), intent(inout) :: wfac(nlt)
character(len=*),     intent(in)  :: data_directory

! Locals
integer :: nl
real(kind=kind_real) :: planck, c, hc, oavo, hcoavo, rlamm
real(kind=kind_real) :: sdom, sdet, rlam
real(kind=kind_real) :: t, tlog, fac
real(kind=kind_real), parameter :: a0 = 0.9976_kind_real
real(kind=kind_real), parameter :: a1 = 0.2194_kind_real
real(kind=kind_real), parameter :: a2 = 5.554E-2_kind_real
real(kind=kind_real), parameter :: a3 = 6.7E-3_kind_real
real(kind=kind_real), parameter :: b0 = 5.026_kind_real
real(kind=kind_real), parameter :: b1 = -0.01138_kind_real
real(kind=kind_real), parameter :: b2 = 9.552E-6_kind_real
real(kind=kind_real), parameter :: b3 = -2.698E-9_kind_real


! Sets parameters for ocean irradiance.

! Obtain Light data
!call lidata(lam, aw, bw, ac, bc, bpic, data_directory)

! Quanta conversion
planck = 6.6256E-34_kind_real   !Plancks constant J sec
c = 2.998E8_kind_real      !speed of light m/sec
hc = 1.0_kind_real/(planck*c)
oavo = 1.0_kind_real/6.023E23_kind_real   ! 1/Avogadros number
hcoavo = hc*oavo
do nl = 1, nlt
  rlamm = real(lam(nl), kind=kind_real)*1.0E-9_kind_real  !lambda in m
  wtoq(nl) = rlamm*hcoavo        !Watts to quanta conversion
enddo

! CDOM and detrital absorption exponent
sdom = 0.014_kind_real
sdet = 0.013_kind_real      !from Gallegos et al., 2011, small detritus nm-1
do nl = 1,nlt
  rlam = real(lam(nl), kind=kind_real)
  excdom(nl) = exp(-sdom*(rlam-443.0_kind_real))
  exdet(nl)  = exp(-sdet*(rlam-440.0_kind_real))
enddo

!  Spectral reflectance factor (from Frouin et al., 1996)
do nl = 1,nlt
  rlam = real(lam(nl), kind=kind_real)
  if (lam(nl) .lt. 900)then
    t = exp(-(aw(nl)+0.5_kind_real*bw(nl)))
    tlog = log(1.0e-36_kind_real + t)
    fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
    wfac(nl) = min(fac, 1.0_kind_real)
    wfac(nl) = max(fac, 0.0_kind_real)
  else
    fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
    wfac(nl) = max(fac, 0.0_kind_real)
  endif
enddo

end subroutine setlte

! --------------------------------------------------------------------------------------------------

end module setlte_mod
