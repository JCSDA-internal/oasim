! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module slingo_mod

! oasim
use oasim_constants_mod

implicit none
private
public slingo

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine slingo(asl, bsl, csl, dsl, esl, fsl, rmu0, cldtau, clwp, cre, tcd, tcs)

!  Slingo's (1989) Delta-Eddington approximation for the two-
!  stream equations applied to clouds.
!  Inputs:
!       rmu0    Kasten's approx for cosine of solar zenith angle
!       cldtau  cloud optical thickness (at 0.6 um)
!       clwp    liquid water path in cloud (g/m2)
!       cre     cloud droplet effective radius (um)
!  Outputs
!       ica     index for translating cloud arrays to light arrays
!       Tcd     spectral transmittance for downwelling direct irradiance
!       Tcs     spectral transmittance for downwelling diffuse irradiance

!  This version uses global mean re for the last 2 bands 2.91-3.42 and
!  3.42-4.00 um to account for unknown errors in the GEOS-5/ESMF system

! Arguments
real(kind=kind_real), intent(in)  :: asl(ncld)
real(kind=kind_real), intent(in)  :: bsl(ncld)
real(kind=kind_real), intent(in)  :: csl(ncld)
real(kind=kind_real), intent(in)  :: dsl(ncld)
real(kind=kind_real), intent(in)  :: esl(ncld)
real(kind=kind_real), intent(in)  :: fsl(ncld)
real(kind=kind_real), intent(in)  :: rmu0
real(kind=kind_real), intent(in)  :: cldtau
real(kind=kind_real), intent(in)  :: clwp
real(kind=kind_real), intent(in)  :: cre
real(kind=kind_real), intent(out) :: tcd(ncld)
real(kind=kind_real), intent(out) :: tcs(ncld)

! Locals
integer :: nc
real(kind=kind_real) :: re, remean, tauc, oneomega, omega, g


Tcd = 0.0_kind_real
Tcs = 0.0_kind_real

!  Compute re as funtion of cldtau and LWP according to eq. 1 in
!  Slingo.
!   tau is derived at this wavelength (0.6 um) in the ISCCP data set
!      re = clwp*bsl(9)/(cldtau - clwp*asl(9))
!      re = min(re,15.0)  !block high re -- produces excessive direct
!  Changes to the ISCCP-D2 data set make this relationship untenable
!  (excessive re's are derived).  Instead choose a fixed re of 10 um
!  for ocean (Kiehl et al., 1998 -- J. Clim.)
!       re = 10.0
!  Paper by Han et al., 1994 (J.Clim.) show mean ocean cloud radius
!  = 11.8 um
!       re = 11.8

! Mean of Kiehl and Han
re = (10.0_kind_real+11.8_kind_real)/2.0_kind_real
remean = re

! Compute spectral cloud characteristics
! If MODIS re is available use it; otherwise use parameterized re above
! I comment this part for now as it cause errors at some profiles, will revisit this later
!if (cre .ge. 0.0_kind_real) then   !use modis re
!  re = cre
!endif

do nc = 1,22
  tauc = clwp*(asl(nc)+bsl(nc)/re)
  oneomega = csl(nc) + dsl(nc)*re
  omega = 1.0_kind_real - oneomega
  g = esl(nc) + fsl(nc)*re

  call slingomath(tauc, omega, g, rmu0, tcd(nc), tcs(nc))
enddo

! Slingo bands 23 and 24 fail due to re.  Workaround is to use mean ocean re
do nc = 23,24
  tauc = clwp*(asl(nc)+bsl(nc)/remean)
  oneomega = csl(nc) + dsl(nc)*remean
  omega = 1.0_kind_real - oneomega
  g = esl(nc) + fsl(nc)*remean

  call slingomath(tauc, omega, g, rmu0, tcd(nc), tcs(nc))

enddo

end subroutine slingo

! --------------------------------------------------------------------------------------------------

subroutine slingomath(tauc, omega, g, rmu0, tcds, tcss)

! Mathematical calculations for the Slingo model

! Inputs
real(kind=kind_real), intent(in)  :: tauc
real(kind=kind_real), intent(in)  :: omega
real(kind=kind_real), intent(in)  :: g
real(kind=kind_real), intent(in)  :: rmu0
real(kind=kind_real), intent(out) :: tcds
real(kind=kind_real), intent(out) :: tcss

! Locals
real(kind=kind_real) :: u1, b0, bmu0, f, U2, alpha1, alpha2, alpha3, alpha4, sqarg, eps, rm, e, val1
real(kind=kind_real) :: val2, rnum, rden, gama1, gama2, tdb, esq, rmsq, em, val3, rdif, tdif, tdir

u1 = 7.0_kind_real/4.0_kind_real
b0 = 3.0_kind_real/7.0_kind_real*(1.0_kind_real-g)
bmu0 = 0.5_kind_real - 0.75_kind_real*rmu0*g/(1.0_kind_real+g)
f = g*g

u2 = u1*(1.0_kind_real-((1.0_kind_real-omega)/(7.0_kind_real*omega*b0)))
u2 = max(u2,0.0_kind_real)
alpha1 = u1*(1.0_kind_real-omega*(1.0_kind_real-b0))
alpha2 = u2*omega*b0
alpha3 = (1.0_kind_real-f)*omega*bmu0
alpha4 = (1.0_kind_real-f)*omega*(1.0_kind_real-bmu0)
sqarg = alpha1*alpha1 - alpha2*alpha2
sqarg = max(sqarg,1.0e-17_kind_real)
eps = sqrt(sqarg)

rm = alpha2/(alpha1+eps)

e = exp(-eps*tauc)

val1 = 1.0_kind_real - omega*f
val2 = eps*eps*rmu0*rmu0

rnum = val1*alpha3 - rmu0*(alpha1*alpha3+alpha2*alpha4)

rden = val1*val1 - val2

gama1 = rnum/rden
rnum = -val1*alpha4 - rmu0*(alpha1*alpha4+alpha2*alpha3)

gama2 = rnum/rden
tdb = exp(-val1*tauc/rmu0)
esq = e*e
rmsq = rm*rm
em = esq*rmsq
val3 = 1.0_kind_real - em

rdif = rm*(1.0_kind_real-esq)/val3
tdif = e*(1.0_kind_real-rmsq)/val3
tdir = -gama2*tdif - gama1*tdb*rdif + gama2*tdb
tcds = tdb
tcss = max(tdir,0.0_kind_real)

end subroutine slingomath

! --------------------------------------------------------------------------------------------------

end module slingo_mod
