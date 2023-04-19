! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module edeu_mod

! oasim
use oasim_constants_mod
use radmod_mod, only: radmod

implicit none
private
public edeu

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine edeu(km, lam, aw, bw, ac, bc, bpic, WtoQ, Ed, Es, H, P, excdom, exdet, rmud, tirrq, &
                cdomabsq, avgq, sfceu)

! Model of irradiance in the water column.  Accounts for three irradiance streams:

! Edz = direct downwelling irradiance
! Esz = diffuse downwelling irradiance
! Euz = diffuse upwelling irradiance

! Propagation is done in energy units, tests are done in quanta, final is quanta for phytoplankton
! growth.

! Arguments
integer,              intent(in)    :: km
integer,              intent(in)    :: lam(nlt)
real(kind=kind_real), intent(in)    :: aw(nlt)
real(kind=kind_real), intent(in)    :: bw(nlt)
real(kind=kind_real), intent(in)    :: ac(nchl,nlt)
real(kind=kind_real), intent(in)    :: bc(nchl,nlt)
real(kind=kind_real), intent(in)    :: bpic(nlt)
real(kind=kind_real), intent(in)    :: WtoQ(nlt)
real(kind=kind_real), intent(in)    :: Ed(nlt)
real(kind=kind_real), intent(in)    :: Es(nlt)
real(kind=kind_real), intent(in)    :: H(km)
real(kind=kind_real), intent(in)    :: P(km,ntyp)
real(kind=kind_real), intent(in)    :: excdom(nlt)
real(kind=kind_real), intent(in)    :: exdet(nlt)
real(kind=kind_real), intent(in)    :: rmud
real(kind=kind_real), intent(out)   :: tirrq(km)
real(kind=kind_real), intent(out)   :: cdomabsq(km)
real(kind=kind_real), intent(inout) :: avgq(km)
real(kind=kind_real), intent(out)   :: sfceu(nlt)

! Locals
integer :: k, n, nl
real(kind=kind_real) :: acdom(nlt), Edtop(nlt),Estop(nlt)
real(kind=kind_real) :: deltaE(km)
real(kind=kind_real) :: Edz(nlt,km), Esz(nlt,km), Euz(nlt,km), sfceun(nlt,km)
real(kind=kind_real) :: fchl(nchl)
real(kind=kind_real) :: a, actot, adet, bb, bbctot, bbrw, bctot, bdet, bt, ebot, ebotq, etop, etopq
real(kind=kind_real) :: Plte, plte3, rmus, sumq, zd, zirr, zirrq

! Parameters
real(kind=kind_real), parameter, dimension(6) :: bbrc = [0.002_kind_real,  0.00071_kind_real, &
                                                         0.0032_kind_real, 0.00071_kind_real, &
                                                         0.0029_kind_real, 0.002_kind_real]
real(kind=kind_real), parameter :: bbrpic    = 0.01_kind_real      !Balch et al., 1996, low end
real(kind=kind_real), parameter :: bbrd      = 0.005_kind_real     !Gallegos et al., 2011 small detritus
real(kind=kind_real), parameter :: adstar    = 8.0E-5_kind_real    !Gallegos et al., 2011 m2/mg
real(kind=kind_real), parameter :: bdstar    = 0.00115_kind_real   !Gallegos et al., 2011 m2/mg
real(kind=kind_real), parameter :: acdomstar = 2.98E-4_kind_real   !Yacobi et al., 2003 m2/mg
real(kind=kind_real), parameter :: dmax      = 500.0_kind_real     !depth at which Ed = 0


!  Constants and initialize
rmus = 1.0_kind_real/0.83_kind_real            !avg cosine diffuse down
tirrq(:) = 0.0_kind_real
cdomabsq(:) = 0.0_kind_real
deltaE(:) = 0.0_kind_real
acdom(:) = 0.0_kind_real
avgq(:) = 0.0_kind_real
sfceu(:)= 0.0_kind_real
Ebot = 0.0_kind_real
bbrw    = 0.5_kind_real ! we need to confirm with Cecile

do nl = 1,nlt
!do nl = npst,npnd
  Edtop(nl) = Ed(nl)
  Estop(nl) = Es(nl)
  Ebot = Ebot + (Ed(nl)+Es(nl))
enddo

!  Convert to quanta: divide by Avos # to get moles quanta; then mult by
!  1E6 to get uM or uEin
Ebotq = 0.0_kind_real
do nl = 1,nlt
!do nl = npst,npnd   !PAR range only 350-700nm
  Ebotq = Ebotq + (Edtop(nl)+Estop(nl))*WtoQ(nl)*1.0E6
enddo
do k = 1,km
  if (H(k) < 1.0E10_kind_real)then
    Etop = Ebot
    Etopq = Ebotq
    zd = min(Dmax,H(k))
    zirr = 0.0_kind_real
    zirrq = 0.0_kind_real
    do nl = 1,nlt
    !do nl = npst,npnd
      Edz(nl,k) = 0.0_kind_real
      Esz(nl,k) = 0.0_kind_real
      Euz(nl,k) = 0.0_kind_real
      sfceun(nl,k) = 0.0_kind_real
      actot = 0.0_kind_real
      bctot = 0.0_kind_real
      bbctot = 0.0_kind_real
      bbctot = 0.0_kind_real
      Plte = max(P(k,ncs+4),0.0_kind_real)
      acdom(nl) = Plte*12.0_kind_real*acdomstar*excdom(nl) !cdoc in units of uM
      ! No detritus or PIC optics
      do n = 1,nchl
        Plte = max(P(k,nnut+n),0.0_kind_real)
        actot  = actot  + Plte*ac(n,nl)
        bctot  = bctot  + Plte*bc(n,nl)
        bbctot = bbctot + Plte*bbrc(n)*bc(n,nl)
      enddo
      Plte = max(P(k,nds),0.0_kind_real)
      adet = Plte*adstar*exdet(nl)
      bdet = Plte*bdstar*(555.0_kind_real/real(lam(nl), kind=kind_real)**0.5_kind_real)
      a  = aw(nl) + acdom(nl) + actot + adet
      Plte3 = max(P(k,ncs+3),0.0_kind_real)
      bt = bw(nl) + bctot + bdet + bpic(nl)*Plte3
      bb = bbrw*bw(nl) + bbctot + bbrd*bdet + bbrpic*bpic(nl)*Plte3
      bb = max(bb,0.0002_kind_real)
      if (Edtop(nl) .ge. 1.0E-4_kind_real .or. Estop(nl) .ge. 1.0E-4_kind_real) then
        call radmod(zd, Edtop(nl), Estop(nl), rmud, a, bt, bb, Edz(nl,k), Esz(nl,k), Euz(nl,k), sfceun(nl,k))
      endif
      Edtop(nl) = Edz(nl,k)
      Estop(nl) = Esz(nl,k)
      zirr = zirr + (Edz(nl,k)+Esz(nl,k)+Euz(nl,k))
      ! surface normalized upwelling irradiance
      sfceu(nl) = Edtop(nl)*sfceun(nl,k)
    enddo
    Ebot = zirr
    deltaE(k) = Etop - Ebot
    do nl = 1,nlt
    !do nl = npst,npnd
      sumQ = (Edz(nl,k)+Esz(nl,k)+Euz(nl,k))*WtoQ(nl)*1.0E6_kind_real
      cdomabsq(k) = cdomabsq(k) + acdom(nl)*sumQ
      zirrq = zirrq + sumQ
    enddo
    Ebotq = zirrq

    ! tirrq is the average of the natural log between Etopq and Ebotq converted to sqrt using the
    ! identity relationships for speedup

    tirrq(k) = sqrt(Etopq*Ebotq)*rmus
  endif
enddo

! Irradiance summary loops
do k = 1,km
  avgq(k) = avgq(k) + tirrq(k)
enddo

end subroutine edeu

! --------------------------------------------------------------------------------------------------

end module edeu_mod
