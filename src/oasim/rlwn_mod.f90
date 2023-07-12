! (C) Copyright 2023 UCAR
! (C) Copyright 2023 United States Government as represented by the Administrator of the 
!     National Aeronautics and Space Administration. All Rights Reserved. 
!                                                                                           
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.                      

module rlwn_mod

! oasim                                                                                                                            
use oasim_constants_mod
use edeu_mod, only: edeu

implicit none
private
public rlwn

! --------------------------------------------------------------------------------------------------

contains

! -------------------------------------------------------------------------------------------------- 
subroutine rlwn(km, lam, cosz, l_chan, aw, bw, ac, bc, bpic, WtoQ, Ed, Es, h, phyto, excdom, exdet, cdet, pic, cdc, rlwnref)

! Arguments 
integer,              intent(in)  :: km
real(kind=kind_real), intent(in)  :: cosz
real(kind=kind_real), intent(in)  :: l_chan(:)
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
real(kind=kind_real), intent(out) :: rlwnref(:)

! Locals                                                                                                                           
real(kind=kind_real) :: pi, Q
real(kind=kind_real) :: rn, rn2
real(kind=kind_real) :: rho, rmudl, rmud
real(kind=kind_real) :: tirrq(km),p(km,ntyp)
real(kind=kind_real) :: cdomabsq(km)
real(kind=kind_real) :: avgq(km)
real(kind=kind_real) :: sfceu(nlt)
integer              :: i, ind

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


pi  = dacos(-1.0_kind_real)
rn  = 1.341     !index of refraction
rho = 0.021     !surface reflectance
Q   = pi        !radiance:irradiance distribution function 
rn2 = rn*rn     !index refr squared 

! Compute average cosine for direct irradiance in the water column given 
!solar zenith angle (in radians) at surface.

rmudl = 1.0/cos( asin( sin(acos(cosz))/rn ) )   !avg cosine direct (1 over)                                                                     
rmud = min(rmudl,1.5)
rmud = max(rmud,0.0)

call edeu(km, lam, aw, bw, ac, bc, bpic, WtoQ, Ed, Es, H, p, excdom, exdet, rmud, tirrq, cdomabsq, &
          avgq, sfceu)

do i = 1,size(l_chan)
   ind = MINLOC(abs(lam-l_chan(i)), DIM=1)
   rlwnref(i) = (1.0-rho)*sfceu(ind)/(rn2*Q)
enddo

end subroutine rlwn

! --------------------------------------------------------------------------------------------------                               

end module rlwn_mod
