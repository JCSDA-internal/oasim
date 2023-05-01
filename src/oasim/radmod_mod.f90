! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module radmod_mod

! oasim
use oasim_constants_mod

implicit none
private
public radmod

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

! Reads in radiative transfer data: specifically water data (seawater absorption and total
! scattering coefficients, and chl-specific absorption and total scattering data for several
! phytoplankton groups).  PAR (350-700) begins at index 3, and ends at index 17.

subroutine radmod(zd, edtop, estop, rmud, a, bt, bb, edz, esz, euz, sfceun, spinup_in)

! Model of irradiance in the water column.  Accounts for three
! irradiance streams:

! Edz = direct downwelling irradiance
! Esz = diffuse downwelling irradiance
! Euz = diffuse upwelling irradiance

! Uses Ackelson's (1994, JGR) mod's to the Aas (1987, AO)
! two-stream model.

! Propagation is done in energy units, tests are done in quanta,
! final is quanta for phytoplankton growth.

! Commented out terms produce a max error of
! 0.8% in Esz for a > 0.004 and bb > 0.0001 and
! 3.9% in Euz for a > 0.004 and bb > 0.00063

! Arguments
real(kind=kind_real), intent(in)  :: zd
real(kind=kind_real), intent(in)  :: edtop
real(kind=kind_real), intent(in)  :: estop
real(kind=kind_real), intent(in)  :: rmud
real(kind=kind_real), intent(in)  :: a
real(kind=kind_real), intent(in)  :: bt
real(kind=kind_real), intent(in)  :: bb
real(kind=kind_real), intent(out) :: edz
real(kind=kind_real), intent(out) :: esz
real(kind=kind_real), intent(out) :: euz
real(kind=kind_real), intent(out) :: sfceun
logical, optional,    intent(in)  :: spinup_in

! Locals
logical :: spinup, euz_flag
real(kind=kind_real) :: rmus, rmuu, cd, au, bu, cu, as, bs, cs, bd, fd, bquad, cquad, sqarg, a1, a2
real(kind=kind_real) :: s, sedz, a2ma1, rm, rn, c2, ta2z, ta2, eutmp

! Paramters
real(kind=kind_real), parameter :: rbot = 0.0_kind_real   !bottom reflectance
real(kind=kind_real), parameter :: rd = 1.5_kind_real     !these are taken from Ackleson, et al. 1994 (JGR)
real(kind=kind_real), parameter :: ru = 3.0_kind_real

! Spinup flag
spinup = .false.
if (present(spinup_in)) spinup = spinup_in
if (spinup) then
  euz_flag = .false.
else
  euz_flag = .true.
endif

!  Constants
rmus = 1.0_kind_real/0.83_kind_real            !avg cosine diffuse down
rmuu = 1.0_kind_real/0.4_kind_real             !avg cosine diffuse up

! Downwelling irradiance: Edz, Esz
! Compute irradiance components at depth

cd = (a+bt)*rmud
edz = edtop*exp(-cd*zd)
au = a*rmuu
bu = ru*bb*rmuu
cu = au+bu
as = a*rmus
bs = rd*bb*rmus
cs = as+bs
bd = bb*rmud
fd = (bt-bb)*rmud
bquad = cs - cu
cquad = bs*bu - cs*cu
sqarg = bquad*bquad - 4.0_kind_real*cquad
a1 = 0.5_kind_real*(-bquad + sqrt(sqarg))
a2 = 0.5_kind_real*(-bquad - sqrt(sqarg))
s = -(bu*bd + cu*fd)
sedz = s*edz
a2ma1 = a2 - a1
rm = sedz/(a1*a2ma1)
rn = sedz/(a2*a2ma1)
c2 = estop - rm + rn
ta2z = exp(a2*zd)

esz = c2*ta2z + rm - rn
esz = max(esz,0.0_kind_real)

if (euz_flag) then
  eutmp = ((a2+cs)*c2)*ta2z + cs*rm - cs*rn - fd*edz
  euz = eutmp/bu
  euz = max(euz,0.0_kind_real)
else
  euz = 0.0_kind_real
endif

!  Computes surface normalized upwelling irradiance.      
sfceun = (a2+cd)/bu

end subroutine radmod

! --------------------------------------------------------------------------------------------------

end module radmod_mod
