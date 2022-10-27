! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module glight_mod

! oasim
use oasim_constants_mod
use daysetrad_mod, only: daysetrad
use edeu_mod,      only: edeu

implicit none
private
public glight

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine glight(km, is_midnight, cosz, lam, aw, bw, ac, bc, bpic, excdom, exdet, wtoq, ed, es, &
                  h, phyto, cdet, pic, cdc, tirrq, cdomabsq, avgq, dt)

! Arguments
integer,              intent(in)  :: km
logical,              intent(in)  :: is_midnight
real(kind=kind_real), intent(in)  :: cosz
integer,              intent(in)  :: lam(nlt)
real(kind=kind_real), intent(in)  :: aw(nlt)
real(kind=kind_real), intent(in)  :: bw(nlt)
real(kind=kind_real), intent(in)  :: ac(nchl,nlt)
real(kind=kind_real), intent(in)  :: bc(nchl,nlt)
real(kind=kind_real), intent(in)  :: bpic(nlt)
real(kind=kind_real), intent(in)  :: excdom(nlt)
real(kind=kind_real), intent(in)  :: exdet(nlt)
real(kind=kind_real), intent(in)  :: wtoq(nlt)
real(kind=kind_real), intent(in)  :: ed(nlt)
real(kind=kind_real), intent(in)  :: es(nlt)
real(kind=kind_real), intent(in)  :: h(km)
real(kind=kind_real), intent(in)  :: phyto(km,nchl)
real(kind=kind_real), intent(in)  :: cdet(km)
real(kind=kind_real), intent(in)  :: pic(km)
real(kind=kind_real), intent(in)  :: cdc(km)
real(kind=kind_real), intent(out) :: tirrq(km)
real(kind=kind_real), intent(out) :: cdomabsq(km)
real(kind=kind_real), intent(out) :: avgq(km)
real(kind=kind_real), intent(in)  :: dt

! Locals
integer :: k, n, nl
real(kind=kind_real) :: p(km,ntyp), deltae(km)
real(kind=kind_real) :: acdom(km,nlt)
real(kind=kind_real) :: sza, sinszaw, rmud, rmudl, szaw

! Parameters
real(kind=kind_real), parameter :: rn = 1.341_kind_real


p = 0.0_kind_real
p(:,nnut+1) = phyto(:,1)
p(:,nnut+2) = phyto(:,2)
p(:,nnut+3) = phyto(:,3)
p(:,nnut+4) = phyto(:,4)
p(:,nnut+5) = phyto(:,5)
p(:,nnut+6) = phyto(:,6)
p(:,nds)    = cdet
p(:,ncs+3)  = pic
p(:,ncs+4)  = cdc

!  Set daily parameters for ocean irradiance.
if (is_midnight) then
  call daysetrad(km, avgq, dt)
endif

! Compute average cosine for direct irradiance in the water column given solar zenith angle
! (in radians) at surface.
sza = acos(cosz)
sinszaw = sin(sza)/rn
szaw = asin(sinszaw)
rmudl = 1.0/cos(szaw)   !avg cosine direct (1 over)
rmud = min(rmudl,1.5)
rmud = max(rmud,0.0)

!  Compute underwater irradiance and total quanta
call edeu(km, lam, aw, bw, ac, bc, bpic, WtoQ, Ed, Es, H, P, excdom, exdet, rmud, tirrq, cdomabsq, &
          avgq)

end subroutine glight

! --------------------------------------------------------------------------------------------------

end module glight_mod
