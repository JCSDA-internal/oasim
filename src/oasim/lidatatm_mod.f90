! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module lidatatm_mod

! oasim
use oasim_constants_mod

implicit none
private
public lidatatm

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine lidatatm(fobar, thray, oza, awv, ao, aco2, data_directory)

! Reads in radiative transfer data: specifically water data (seawater absorption and total
! scattering coefficients, and chl-specific absorption and total scattering data for several
! phytoplankton groups).  PAR (350-700) begins at index 3, and ends at index 17.

! Arguments
real(kind=kind_real), intent(out) :: fobar(nlt)
real(kind=kind_real), intent(out) :: thray(nlt)
real(kind=kind_real), intent(out) :: oza(nlt)
real(kind=kind_real), intent(out) :: awv(nlt)
real(kind=kind_real), intent(out) :: ao(nlt)
real(kind=kind_real), intent(out) :: aco2(nlt)
character(len=*),     intent(in)  :: data_directory

! Locals
integer :: nl, ilam
character(len=2048) :: atmo25b
character(len=80) :: title
real(4) :: sfobar, sthray, soza, sawv, sao, saco2

! File
atmo25b = trim(data_directory)//'/atmo25b.dat'

! Atmospheric data file
open(4, file=trim(atmo25b), status='old', form='formatted')
read(4,'(a80)') title
read(4,'(a80)') title
read(4,'(a80)') title
read(4,'(a80)') title
do nl = 1, nlt
  read(4,'(i5,6f11.4)') ilam, sfobar, sthray, soza, sawv, sao, saco2
  fobar(nl) = real(sfobar, kind=kind_real)
  thray(nl) = real(sthray, kind=kind_real)
  oza(nl)   = real(soza,   kind=kind_real)
  awv(nl)   = real(sawv,   kind=kind_real)
  ao(nl)    = real(sao,    kind=kind_real)
  aco2(nl)  = real(saco2,  kind=kind_real)
enddo
close(4)

end subroutine lidatatm

! --------------------------------------------------------------------------------------------------

end module lidatatm_mod
