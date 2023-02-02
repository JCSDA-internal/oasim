! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module daysetrad_mod

! oasim
use oasim_constants_mod

implicit none
private
public daysetrad

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine daysetrad(km, avgq, dt)

! Sets daily parameters for ocean irradiance.

! Arguments
integer,              intent(in)  :: km
real(kind=kind_real), intent(out) :: avgq(km)
real(kind=kind_real), intent(in)  :: dt

! Locals
integer :: k

! Parameters
real(kind=kind_real), parameter :: sday=86400.0_kind_real  !seconds per day

!Initialize avgq 
avgq(:) = 1.0_kind_real

! Compute average quanta; Initialize light history arrays
do k = 1,km
  avgq(k) = avgq(k)*sday/dt
enddo

end subroutine daysetrad

! --------------------------------------------------------------------------------------------------

end module daysetrad_mod
