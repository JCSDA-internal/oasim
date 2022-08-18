! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module oasim_mod

! oasim
use oasim_constants_mod
use glight_mod,    only: glight
use ocalbedo_mod,  only: ocalbedo
use setlte_mod,    only: setlte
use setsfclte_mod, only: setsfclte
use sfcirr_mod,    only: sfcirr

implicit none
private
public oasim

type :: oasim
  integer :: lam(nlt)
  real(kind=kind_real) :: aw(nlt), bw(nlt), excdom(nlt), exdet(nlt), wtoq(nlt), wfac(nlt)
  real(kind=kind_real) :: ac(nchl,nlt), bc(nchl,nlt)
  real(kind=kind_real) :: bpic(nlt)
  real(kind=kind_real) :: rad, pi2
  real(kind=kind_real) :: fobar(nlt), thray(nlt), oza(nlt), awv(nlt), ao(nlt), aco2(nlt)
  real(kind=kind_real) :: am, vi
  real(kind=kind_real) :: asl(ncld), bsl(ncld), csl(ncld), dsl(ncld), esl(ncld), fsl(ncld)
  integer :: ica(nlt)
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: run
    final :: dummy_final
endtype oasim

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, data_directory)

! Arguments
class(oasim),     intent(inout) :: self
character(len=*), intent(in)    :: data_directory  ! Path to where datafiles are stored

! Construct the oasim structures
! ------------------------------
call setlte(self%lam, self%aw, self%bw, self%ac, self%bc, self%bpic, self%excdom, self%exdet, &
            self%WtoQ, self%wfac, data_directory)

call setsfclte(self%rad, self%pi2, self%lam, self%fobar, self%thray, self%oza, self%awv, self%ao, &
               self%aco2, self%am, self%vi, self%asl, self%bsl, self%csl, self%dsl, self%esl, &
               self%fsl, self%ica, data_directory)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(oasimobj)

! Arguments
class(oasim), intent(inout) :: oasimobj

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine run(self)

! Arguments
class(oasim), intent(in) :: self

! Locals
integer :: lm
logical :: is_midnight
real(kind=kind_real) :: daycor, cosz, slp, wspd, ozone, wvapor, relhum
real(kind=kind_real) :: ta(nlt), wa(nlt), asym(nlt), ed(nlt), es(nlt)
real(kind=kind_real) :: cov, cldtau, clwp, cldre, sunz, rod(nlt), ros(nlt), dt
real(kind=kind_real), allocatable :: dh(:), cdet(:), pic(:), cdc(:), tirrq(:), cdomabsq(:), avgq(:)
real(kind=kind_real), allocatable :: phyto(:,:)

!
!! Use Beer's Law to compute flux divergence
!!------------------------------------------
!uvr = (pruvf+pruvr)*fr
!par = (prpaf+prpar)*fr
!z   = 0.0
!
!if ( associated(qsw) ) then
!  qsw(:,:,1) = uvr + par
!
!  do l=2,lm
!    z = z + dh(:,:,l-1)
!    qsw(:,:,l  ) = uvr*exp(-kuvr*z) + par*exp(-kpar*z)
!    qsw(:,:,l-1) = qsw(:,:,l-1) - qsw(:,:,l)
!  enddo
!  z = z + dh(:,:,lm)
!  qsw(:,:,lm) = qsw(:,:,lm) - (uvr*exp(-kuvr*z) + par*exp(-kpar*z))
!end if
!
!if( .not. associated(kparx) ) deallocate(kpar)
!deallocate(z   )
!deallocate(par )
!deallocate(uvr )
!
!
!hr=1.0
!! Obtain Earth-Sun distance
!    rday = float(DOY) + hr/24.0
!    daycor = (1.0+1.67E-2*cos(self%pi2*(rday-3.0)/365.0))**2
!!   if (Is_Leap)daycor = (1.0+1.67E-2*cos(self%pi2*(rday-3.0)/366.0))**2
!!
!do j = 1, JM
!  do i = 1, IM
!      if (DH(i,j,1) < 1.0E10 .and. COSZ(i,j) > 0.0)then
!       slp = PS(i,j)*0.01  ! convert from Pa to mbar
!       wspd = WSM(i,j)
!       ozone = OZ(i,j)
!       wvapor = WV(i,j)
!       relhum = RH(i,j)
!       do nl = 1,nlt
!        ta(nl) = TAUA(nl)%b(i,j)
!        wa(nl) = SSALB(nl)%b(i,j)
!        asym(nl) = ASYMP(nl)%b(i,j)
!       enddo
!       cov = CCOV(i,j)
!       cldtau = CLDTC(i,j)
!       clwp = RLWP(i,j)
!       cldre = CDRE(i,j)
!
!      ! There are mismatches between the ocean, land and atmosphere in GEOS-5
!      ! and live ocean points for MOM that do not have a corresponding
!      ! atmosphere.  These are called "grottoes" because they are assumed
!      ! to have ocean underneath with land overhead.  Set irradiance to 0
!      ! to represent this condition.
!      if (slp < 0.0 .or. slp >1.0E10)then
!        Ed = 0.0
!        Es = 0.0
!        TIRRQ(i,j,:) = 0.0
!        CDOMABSQ(i,j,:) = 0.0
!        AVGQ(i,j,:) = 0.0
!      else
!        ! Spectral irradiance just above surface
        call sfcirr(self%lam, self%fobar, self%thray, self%oza, self%awv, self%ao, self%aco2, &
                    self%asl, self%bsl, self%csl, self%dsl, self%esl, self%fsl, self%ica, daycor, &
                    cosz, slp, wspd, ozone, wvapor, relhum, ta, wa, asym, self%am, self%vi, &
                    cov, cldtau, clwp, cldre, ed, es)
!
!!   Spectral irradiance just below surface
!        sunz = acos(COSZ(i,j))*self%rad
        call ocalbedo(self%rad, self%wfac, sunz, wspd, rod, ros)
!        do nl = 1,nlt
!         Ed(nl) = Ed(nl)*(1.0-rod(nl))
!         Es(nl) = Es(nl)*(1.0-ros(nl))
!        enddo
!
!        PHYTO(:,1) = DIATOM(i,j,:)
!        PHYTO(:,2) = CHLORO(i,j,:)
!        PHYTO(:,3) = CYANO(i,j,:)
!        PHYTO(:,4) = COCCO(i,j,:)
!        PHYTO(:,5) = DINO(i,j,:)
!        PHYTO(:,6) = PHAEO(i,j,:)
!!   Spectral irradiance in the water column
        call glight(lm, is_midnight, cosz, self%lam, self%aw, self%bw, self%ac, self%bc, &
                    self%bpic, self%excdom, self%exdet, self%wtoq, ed, es, dh, phyto, &
                    cdet, pic, cdc, tirrq, cdomabsq, avgq, dt)
!!TIRRQ(i,j,:) = 1.0
!!CDOMABSQ(i,j,:) = 0.0
!       endif
!      endif
!     enddo
!    enddo

end subroutine run

! --------------------------------------------------------------------------------------------------

! Not really needed but prevents gnu compiler bug
subroutine dummy_final(self)
  type(oasim), intent(inout) :: self
  end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

end module oasim_mod
