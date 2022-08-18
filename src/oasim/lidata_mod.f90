! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module lidata_mod

! oasim
use oasim_constants_mod

implicit none
private
public lidata

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

! Reads in radiative transfer data: specifically water data (seawater absorption and total
! scattering coefficients, and chl-specific absorption and total scattering data for several
! phytoplankton groups).  PAR (350-700) begins at index 3, and ends at index 17.

subroutine lidata(lam, aw, bw, ac, bc, bpic, data_directory)

! Arguments
integer,              intent(out) :: lam(nlt)
real(kind=kind_real), intent(out) :: aw(nlt)
real(kind=kind_real), intent(out) :: bw(nlt)
real(kind=kind_real), intent(out) :: ac(nchl,nlt)
real(kind=kind_real), intent(out) :: bc(nchl,nlt)
real(kind=kind_real), intent(out) :: bpic(nlt)
character(len=*),     intent(in)  :: data_directory

! Locals
integer :: i, nl, n
character(len=80) :: title
integer :: lambda
real(4) :: saw, sbw, sac, sbc, sbpic
character(len=2048) :: abw25_morel, acbc25_6, pic_sigma

! Filenames
abw25_morel = trim(data_directory)//'/abw25_morel.dat'
acbc25_6 = trim(data_directory)//'/acbc25_6.dat'
pic_sigma = trim(data_directory)//'/pic_sigma.dat'


! Water data files
open(4, file=trim(abw25_morel), status='old', form='formatted')
do i = 1,5
  read(4,'(a80)')title
enddo
do nl = 1,nlt
  read(4,'(i5, f15.4, f10.4)') lambda, saw, sbw
  lam(nl) = lambda
  aw(nl) = real(saw, kind_real)
  bw(nl) = real(sbw, kind_real)
enddo
close(4)

!  Modify for Lee et al. 2015
aw(3)  = 0.0071
aw(4)  = 0.0041
aw(5)  = 0.0034
aw(6)  = 0.0033
aw(7)  = 0.0068
aw(8)  = 0.0099
aw(9)  = 0.0187
aw(10) = 0.0400
aw(11) = 0.0652


! Phytoplankton group chl-specific absorption and total scattering data. Chl-specific absorption
! data is normalized to 440 nm; convert here to actual ac*(440)
open(4, file=trim(acbc25_6), status='old', form='formatted')
do i = 1,6
  read(4,'(a80)') title
enddo
do n = 1,nchl
  read(4, '(a80)') title
  do nl = 1,19
    read(4, '(i4, 2f10.4)') lambda, sac, sbc
    ac(n,nl) = real(sac, kind_real)
    bc(n,nl) = real(sbc, kind_real)
  enddo
  do nl = 20,nlt
    ac(n,nl) = 0.0_kind_real
    bc(n,nl) = 0.0_kind_real
  enddo
enddo
close(4)


! PIC scattering cross section m2/ugCaCO3 (PIC)  (=PIC-specific scattering coefficient
! (from Gordon et al., 2009)
open(4, file=trim(pic_sigma), status='old', form='formatted')
do i = 1, 5
  read(4, '(a80)') title
enddo
do nl = 1, nlt
  read(4, '(i4, f12.6)') lambda, sbpic
  bpic(nl) = real(sbpic, kind_real)
enddo
close(4)

end subroutine lidata

! --------------------------------------------------------------------------------------------------

end module lidata_mod
