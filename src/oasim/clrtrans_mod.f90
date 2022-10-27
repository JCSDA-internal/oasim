! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module clrtrans_mod

! oasim
use oasim_constants_mod

implicit none
private
public clrtrans

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine clrtrans(lam, thray, cosunz, rm, rmp, wspd, relhum, am, vi, ta, wa, asym, td, ts)

! Model for atmospheric transmittance of solar irradiance through a cloudless maritime atmosphere.
! Computes direct and diffuse separately. From Gregg and Carder (1990) Limnology and Oceanography
! 35(8): 1657-1675.

! Td is spectral clear sky direct transmittance
! Ts is spectral clear sky diffuse transmittance

! Arguments
integer,              intent(in)    :: lam(nlt)
real(kind=kind_real), intent(in)    :: thray(nlt)
real(kind=kind_real), intent(in)    :: cosunz
real(kind=kind_real), intent(in)    :: rm
real(kind=kind_real), intent(in)    :: rmp
real(kind=kind_real), intent(in)    :: wspd
real(kind=kind_real), intent(inout) :: relhum
real(kind=kind_real), intent(in)    :: am
real(kind=kind_real), intent(in)    :: vi
real(kind=kind_real), intent(inout) :: ta(nlt)
real(kind=kind_real), intent(inout) :: wa(nlt)
real(kind=kind_real), intent(in)    :: asym(nlt)
real(kind=kind_real), intent(out)   :: Td(nlt)
real(kind=kind_real), intent(out)   :: Ts(nlt)

! Locals
integer :: nl
integer, save :: ifst = 0
real(kind=kind_real), save :: rlamu(nlt)
real(kind=kind_real) :: beta, eta, wa1, afs, bfs, rtra, omegaa, alg, fa, tarm, atra, taa, tas
real(kind=kind_real) :: dray, daer

! First time only
if (ifst .eq. 0)then
  do nl = 1,nlt
    rlamu(nl) = real(lam(nl), kind=kind_real)*1.0e-3_kind_real    !lambda in um
    Td(nl) = 0.0_kind_real
    Ts(nl) = 0.0_kind_real
  enddo
  ifst = 1
endif

!  Obtain aerosol parameters; simplified Navy aerosol model
call navaer(relhum, am, Vi, wspd, beta, eta, wa1, afs, bfs)

! Compute spectral transmittance
do nl = 1,nlt
  ! Rayleigh
  rtra = exp(-thray(nl)*rmp)       !transmittance
  ! Aerosols
  if (ta(nl) .lt. 0.0_kind_real)then
    ta(nl) = beta*rlamu(nl)**eta
  endif
  if (wa(nl) .lt. 0.0_kind_real)then
    omegaa = wa1
  else
    omegaa = wa(nl)
  endif
  if (asym(nl) .ge. 0.0_kind_real)then
    alg = log(1.0_kind_real-asym(nl))
    afs = alg*(1.459_kind_real+alg*(.1595_kind_real+alg*.4129_kind_real))
    bfs = alg*(.0783_kind_real+alg*(-.3824_kind_real-alg*.5874_kind_real))
  endif
  !if (ta(nl) .lt. 0.0_kind_real .or. omegaa .lt. 0.0_kind_real)then
  !  write(6,*)'ERROR in ta or omegaa'
  !  write(6,*)'nl,ta,wa,asym = ', nl,ta(nl),wa(nl),asym(nl)
  !endif
  Fa = 1.0_kind_real - 0.5_kind_real*exp((afs+bfs*cosunz)*cosunz)
  !if (Fa .lt. 0.0_kind_real)then
  !  write(6,*)'ERROR in Fa'
  !  write(6,*)'nl,ta,wa,asym = ', nl,ta(nl),wa(nl),asym(nl)
  !endif
  tarm = ta(nl)*rm
  atra = exp(-tarm)
  taa = exp(-(1.0_kind_real-omegaa)*tarm)
  tas = exp(-omegaa*tarm)
  ! Direct transmittance
  Td(nl) = rtra*atra

  ! Diffuse transmittance
  dray = taa*0.5_kind_real*(1.0_kind_real-rtra**.95_kind_real)
  daer = rtra**1.5_kind_real*taa*Fa*(1.0_kind_real-tas)

  ! Total diffuse
  Ts(nl) = dray + daer

enddo

end subroutine clrtrans

! --------------------------------------------------------------------------------------------------

subroutine navaer(relhum, am, vi, wspd, beta, eta, wa, afs, bfs)

!  Computes aerosol parameters according to a simplified version of the Navy marine aerosol model.

! Arguments
real(kind=kind_real), intent(inout) :: relhum
real(kind=kind_real), intent(in)    :: am
real(kind=kind_real), intent(in)    :: vi
real(kind=kind_real), intent(in)    :: wspd
real(kind=kind_real), intent(out)   :: beta
real(kind=kind_real), intent(out)   :: eta
real(kind=kind_real), intent(inout) :: wa
real(kind=kind_real), intent(out)   :: afs
real(kind=kind_real), intent(out)   :: bfs

! Locals
integer :: i, n
real(kind=kind_real) :: a(3), dndr(3)
real(kind=kind_real) :: wsm, rnum, rden, frh, arg, rval, sumx, sumy, sumxy, sumx2, rlrn, rldndr
real(kind=kind_real) :: gama, rlogc, alpha, cext, asymp, alg

! Paramters
real(kind=kind_real), parameter, dimension(3) :: ro = [0.03, 0.24, 2.0]
real(kind=kind_real), parameter, dimension(3) :: r = [0.1, 1.0, 10.0]
real(kind=kind_real), parameter :: rlam = 0.55

wsm = wspd

! Relative humidity factor
! if (relhum .ge. 100.0_kind_real)relhum = 99.9_kind_real
relhum = min(99.9_kind_real,relhum)
rnum = 2.0_kind_real - relhum/100.0_kind_real
rden = 6.0_kind_real*(1.0_kind_real-relhum/100.0_kind_real)
frh = (rnum/rden)**0.333_kind_real

! Size distribution amplitude components
a(1) = 2000.0_kind_real*am*am
a(2) = 5.866_kind_real*(wsm-2.2_kind_real)
! if (a(2) .lt. 0.5_kind_real)a(2) = 0.5_kind_real
a(2) = max(0.5_kind_real,a(2))
a(3) = 0.01527_kind_real*(wspd-2.2)*0.05_kind_real        !from Hughes 1987
! if (a(3) .lt. 1.4E-5_kind_real)a(3) = 1.4E-5_kind_real
a(3) = max(1.4E-5_kind_real,a(3))

! Compute size distribution at three selected radii according to Navy method
do n = 1,3
  dndr(n) = 0.0_kind_real
  do i = 1,3
    rden = frh*ro(i)
    arg = log(r(n)/rden)*log(r(n)/rden)
    rval = a(i)*exp(-arg)/frh
    dndr(n) = dndr(n) + rval
  enddo
enddo

! Least squares approximation
sumx = 0.0_kind_real
sumy = 0.0_kind_real
sumxy = 0.0_kind_real
sumx2 = 0.0_kind_real
do n = 1,3
  rlrn = log10(r(n))
  rldndr = log10(dndr(n))
  sumx = sumx + rlrn
  sumy = sumy + rldndr
  sumxy = sumxy + rlrn*rldndr
  sumx2 = sumx2 + rlrn*rlrn
enddo
gama = sumxy/sumx2
rlogc = sumy/3.0_kind_real - gama*sumx/3.0_kind_real
alpha = -(gama+3.0_kind_real)
eta = -alpha

! Compute beta
cext = 3.91_kind_real/Vi
beta = cext*rlam**alpha

!  Compute asymmetry parameter -- a function of alpha
if (alpha .gt. 1.2_kind_real)then
  asymp = 0.65_kind_real
else if (alpha .lt. 0.0_kind_real)then
  asymp = 0.82_kind_real
else
  asymp = -0.14167_kind_real*alpha + 0.82_kind_real
endif

! Forward scattering coefficients
alg = log(1.0_kind_real-asymp)
afs = alg*(1.459_kind_real+alg*(.1595_kind_real+alg*.4129_kind_real))
bfs = alg*(.0783_kind_real+alg*(-.3824_kind_real-alg*.5874_kind_real))

! Single scattering albedo at 550; function of RH
wa = (-0.0032_kind_real*am + 0.972_kind_real)*exp(3.06E-4_kind_real*relhum)

end subroutine navaer

! --------------------------------------------------------------------------------------------------

end module clrtrans_mod
