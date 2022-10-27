! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module setsfclte_mod

! oasim
use oasim_constants_mod
use lidatatm_mod, only: lidatatm

implicit none
private
public setsfclte

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine setsfclte(rad, pi2, lam, fobar, thray, oza, awv, ao, aco2, am, vi, asl, bsl, csl, dsl, &
                     esl, fsl, ica, data_directory)

! Computes constants for global irradiance calculations, reads in required data files, and otherwise
! obtains one-time-only information necessary for the run.

! Locals
real(kind=kind_real), intent(out) :: rad
real(kind=kind_real), intent(out) :: pi2
integer,              intent(in)  :: lam(nlt)
real(kind=kind_real), intent(out) :: fobar(nlt)
real(kind=kind_real), intent(out) :: thray(nlt)
real(kind=kind_real), intent(out) :: oza(nlt)
real(kind=kind_real), intent(out) :: awv(nlt)
real(kind=kind_real), intent(out) :: ao(nlt)
real(kind=kind_real), intent(out) :: aco2(nlt)
real(kind=kind_real), intent(out) :: am
real(kind=kind_real), intent(out) :: vi
real(kind=kind_real), intent(out) :: asl(ncld), bsl(ncld), csl(ncld), dsl(ncld), esl(ncld), fsl(ncld)
integer,              intent(out) :: ica(nlt)
character(len=*),     intent(in)  :: data_directory

! Locals
real(kind=kind_real) :: pi

!  Degrees to radians conversion
pi = acos(-1.0_kind_real)
pi2 = pi*2.0_kind_real
rad = 180.0_kind_real/pi

!  Obtain Light data
call lidatatm(fobar, thray, oza, awv, ao, aco2, data_directory)
am = 1.0_kind_real
vi = 25.0_kind_real
call rdslingo(lam, asl, bsl, csl, dsl, esl, fsl, ica, data_directory)

end subroutine setsfclte

! --------------------------------------------------------------------------------------------------

subroutine rdslingo(lam, asl, bsl, csl, dsl, esl, fsl, ica, data_directory)

! Reads cloud parameters by Slingo (1989).

! Arguments
integer,              intent(in)  :: lam(nlt)
real(kind=kind_real), intent(out) :: asl(ncld)
real(kind=kind_real), intent(out) :: bsl(ncld)
real(kind=kind_real), intent(out) :: csl(ncld)
real(kind=kind_real), intent(out) :: dsl(ncld)
real(kind=kind_real), intent(out) :: esl(ncld)
real(kind=kind_real), intent(out) :: fsl(ncld)
integer,              intent(out) :: ica(nlt)
character(len=*),     intent(in)  :: data_directory

! Locals
integer :: n, nl, nc, lamcld
character(len=50) :: title
real(kind=kind_real) :: rnl1(ncld), rnl2(ncld)
real(4) :: rn1, rn2, a4, b4, e4, f4, c4, d4
character(len=2048) :: slingo

slingo = trim(data_directory)//'/slingo.dat'

open(15, file=trim(slingo), status='old', form='formatted')
do n = 1,3
  read(15,'(a50)')title
enddo
do n = 1,ncld
  read(15,'(2f5.2,3x,2f6.3,2f6.3,1pe9.2,1pe8.2)') rn1, rn2, a4, b4, e4, f4, c4, d4
  rnl1(n) = real(rn1, kind=kind_real)
  rnl2(n) = real(rn2, kind=kind_real)
  asl(n)  = real(a4,  kind=kind_real) * 0.01_kind_real
  bsl(n)  = real(b4,  kind=kind_real)
  csl(n)  = real(c4,  kind=kind_real)
  dsl(n)  = real(d4,  kind=kind_real)
  esl(n)  = real(e4,  kind=kind_real)
  fsl(n)  = real(f4,  kind=kind_real) * 0.001_kind_real
enddo
close(15)

! Indices to relate cloud parameters to other light parameters
do nl = 1,nlt
  do nc = 1,ncld
    lamcld = nint(rnl2(nc)*1000.0_kind_real)
    if (lam(nl) .lt. lamcld)then
      ica(nl) = nc
      go to 5
    endif
  enddo
5 continue
enddo

end subroutine rdslingo

! --------------------------------------------------------------------------------------------------

end module setsfclte_mod
