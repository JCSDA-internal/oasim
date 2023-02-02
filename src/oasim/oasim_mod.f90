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
  integer, allocatable :: lam(:),  ica(:)
  real(kind=kind_real) :: rad, pi2, am, vi
  real(kind=kind_real), allocatable :: aw(:), bw(:), excdom(:), exdet(:), wtoq(:), wfac(:)
  real(kind=kind_real), allocatable :: fobar(:), thray(:), oza(:), awv(:), ao(:), aco2(:), bpic(:)
  real(kind=kind_real), allocatable :: ac(:,:), bc(:,:)
  real(kind=kind_real), allocatable :: asl(:), bsl(:), csl(:), dsl(:), esl(:), fsl(:)
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

! Allocate the oasim structures
! -----------------------------
allocate(self%lam(nlt))
allocate(self%ica(nlt))
allocate(self%aw(nlt))
allocate(self%bw(nlt))
allocate(self%excdom(nlt))
allocate(self%exdet(nlt))
allocate(self%wtoq(nlt))
allocate(self%wfac(nlt))
allocate(self%fobar(nlt))
allocate(self%thray(nlt))
allocate(self%oza(nlt))
allocate(self%awv(nlt))
allocate(self%ao(nlt))
allocate(self%aco2(nlt))
allocate(self%bpic(nlt))
allocate(self%ac(nchl,nlt))
allocate(self%bc(nchl,nlt))
allocate(self%asl(ncld))
allocate(self%bsl(ncld))
allocate(self%csl(ncld))
allocate(self%dsl(ncld))
allocate(self%esl(ncld))
allocate(self%fsl(ncld))

! Construct the oasim structures
! ------------------------------
call setlte(self%lam, self%aw, self%bw, self%ac, self%bc, self%bpic, self%excdom, self%exdet, &
            self%WtoQ, self%wfac, data_directory)

call setsfclte(self%rad, self%pi2, self%lam, self%fobar, self%thray, self%oza, self%awv, self%ao, &
               self%aco2, self%am, self%vi, self%asl, self%bsl, self%csl, self%dsl, self%esl, &
               self%fsl, self%ica, data_directory)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(oasim), intent(inout) :: self

! Deallocate the oasim structures
! -------------------------------
deallocate(self%lam)
deallocate(self%aw)
deallocate(self%bw)
deallocate(self%excdom)
deallocate(self%exdet)
deallocate(self%wtoq)
deallocate(self%wfac)
deallocate(self%ac)
deallocate(self%bc)
deallocate(self%bpic)
deallocate(self%fobar)
deallocate(self%thray)
deallocate(self%oza)
deallocate(self%awv)
deallocate(self%ao)
deallocate(self%aco2)
deallocate(self%asl)
deallocate(self%bsl)
deallocate(self%csl)
deallocate(self%dsl)
deallocate(self%esl)
deallocate(self%fsl)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine run(self, km, dt, is_midnight, day_of_year, cosz, slp, wspd, ozone, wvapor, rh, cov, &
               cldtau, clwp, cldre, ta_in, wa_in, asym, dh, cdet, pic, cdc, diatom, chloro, cyano, &
               cocco, dino, phaeo, tirrq, cdomabsq, avgq)

! Arguments
class(oasim),         intent(in)  :: self         ! oasim object
integer,              intent(in)  :: km           ! Number of model levels
real(kind=kind_real), intent(in)  :: dt           ! Time interval
logical,              intent(in)  :: is_midnight  ! Number of model levels
integer,              intent(in)  :: day_of_year  ! Day of the year
real(kind=kind_real), intent(in)  :: cosz         ! Cosine of the solar zenith angle (1)
real(kind=kind_real), intent(in)  :: slp          ! Sea level pressure (hPa)
real(kind=kind_real), intent(in)  :: wspd         ! Surface_wind_speed (m s-1)
real(kind=kind_real), intent(in)  :: ozone        ! Ozone thickness (dobson units)
real(kind=kind_real), intent(in)  :: wvapor       ! Water vapor (cm)
real(kind=kind_real), intent(in)  :: rh           ! Relative humidity (percent)
real(kind=kind_real), intent(in)  :: cov          ! Cloud cover (percent)
real(kind=kind_real), intent(in)  :: cldtau       ! Cloud optical thickness (dimensionless)
real(kind=kind_real), intent(in)  :: clwp         ! Cloud liquid water path (dimensionless)
real(kind=kind_real), intent(in)  :: cldre        ! Cloud droplet effective radius (dimensionless)
real(kind=kind_real), intent(in)  :: ta_in(nlt)   ! Aerosol optical thickness (dimensionless)
real(kind=kind_real), intent(in)  :: wa_in(nlt)   ! Single scattering albedo (dimensionless)
real(kind=kind_real), intent(in)  :: asym(nlt)    ! Asymmetry parameter (dimensionless)
real(kind=kind_real), intent(in)  :: dh(km)       ! Layer thickness (m)
real(kind=kind_real), intent(in)  :: cdet(km)     ! Carbon/nitrogen_detritus_concentration ()
real(kind=kind_real), intent(in)  :: pic(km)      ! Particulate inorganic carbon (ugC l-1)
real(kind=kind_real), intent(in)  :: cdc(km)      ! Colored dissolved organic carbon (uM)
real(kind=kind_real), intent(in)  :: diatom(km)   ! Diatom concentration (mg m-3)
real(kind=kind_real), intent(in)  :: chloro(km)   ! chlorophyte_concentration (mg m-3)
real(kind=kind_real), intent(in)  :: cyano(km)    ! Cyano-bacteria concentration (mg m-3)
real(kind=kind_real), intent(in)  :: cocco(km)    ! Coccolithophore concentration (mg m-3)
real(kind=kind_real), intent(in)  :: dino(km)     ! Dinoflagellate concentration (mg m-3)
real(kind=kind_real), intent(in)  :: phaeo(km)    ! Phaeocystis concentration (mg m-3)
real(kind=kind_real), intent(out) :: tirrq(km)    ! Total irradiance (umol quanta m-2 s-1)
real(kind=kind_real), intent(out) :: cdomabsq(km) ! Absorption of quanta by CDOM (umol quanta m-2 s-1)
real(kind=kind_real), intent(out) :: avgq(km)     ! Average quantum irradiance (Average quantum irradiance)

! Locals
integer :: nl
real(kind=kind_real) :: relhum, daycor, rday, sunz
real(kind=kind_real) :: ed(nlt), es(nlt), rod(nlt), ros(nlt), ta(nlt), wa(nlt)
real(kind=kind_real), allocatable :: phyto(:,:)


! Set the phyto variable
! ----------------------
allocate(phyto(km,6))
phyto(:,1) = diatom
phyto(:,2) = chloro
phyto(:,3) = cyano
phyto(:,4) = cocco
phyto(:,5) = dino
phyto(:,6) = phaeo

! Obtain Earth-Sun distance
! -------------------------
rday = real(day_of_year, kind=kind_real) + 1.0_kind_real/24.0_kind_real
daycor = (1.0_kind_real + 1.67e-2_kind_real * cos(self%pi2*(rday-3.0_kind_real)/365.0_kind_real))**2

if (dh(1) < 1.0e10_kind_real .and. cosz > 0.0_kind_real) then

  ! There are mismatches between the ocean, land and atmosphere in GEOS-5 and live ocean points for
  ! MOM that do not have a corresponding atmosphere.  These are called "grottoes" because they are
  ! assumed to have ocean underneath with land overhead.  Set irradiance to 0 to represent this
  ! condition.

  if (slp < 0.0_kind_real .or. slp > 1.0e10_kind_real) then

    ed = 0.0_kind_real
    es = 0.0_kind_real
    tirrq(:) = 0.0_kind_real
    cdomabsq(:) = 0.0_kind_real
    avgq(:) = 0.0_kind_real

  else

    ! Copy the variables that are modified for internal purposes
    relhum = rh
    ta = ta_in
    wa = wa_in

    ! Spectral irradiance just above surface
    call sfcirr(self%lam, self%fobar, self%thray, self%oza, self%awv, self%ao, self%aco2, &
                self%asl, self%bsl, self%csl, self%dsl, self%esl, self%fsl, self%ica, daycor, &
                cosz, slp, wspd, ozone, wvapor, relhum, ta, wa, asym, self%am, self%vi, cov, &
                cldtau, clwp, cldre, ed, es)

    ! Spectral irradiance just below surface
    sunz = acos(cosz)*self%rad
    call ocalbedo(self%rad, self%wfac, sunz, wspd, rod, ros)

    do nl = 1,nlt
      ed(nl) = ed(nl) * (1.0_kind_real - rod(nl))
      es(nl) = es(nl) * (1.0_kind_real - ros(nl))
    enddo

    ! Spectral irradiance in the water column
    call glight(km, is_midnight, cosz, self%lam, self%aw, self%bw, self%ac, self%bc, self%bpic, &
                self%excdom, self%exdet, self%wtoq, ed, es, dh, phyto, cdet, pic, cdc, tirrq, &
                cdomabsq, avgq, dt)
  endif

endif

end subroutine run

! --------------------------------------------------------------------------------------------------

! Not really needed but prevents gnu compiler bug
subroutine dummy_final(self)
  type(oasim), intent(inout) :: self
  end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

end module oasim_mod
