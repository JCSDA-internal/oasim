! (C) Copyright 2022 UCAR
! (C) Copyright 2022 United States Government as represented by the Administrator of the
!     National Aeronautics and Space Administration. All Rights Reserved.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module oasim_constants_mod

use, intrinsic :: iso_c_binding

implicit none
public

! Type
integer, parameter :: kind_int = c_int
integer, parameter :: kind_real = c_double

! comlte
integer, parameter :: nlt = 33
integer, parameter :: npr = 15
integer, parameter :: nhn = 12
integer, parameter :: npar = npr
integer, parameter :: ncld = 24

integer, parameter :: nnut_defined = 4
integer, parameter :: nchl_defined = 6
integer, parameter :: nzoo_defined = 1
integer, parameter :: ndet_defined = 3
integer, parameter :: ncar_defined = 5
integer, parameter :: nnut = nnut_defined
integer, parameter :: nchl = nchl_defined
integer, parameter :: nzoo = nzoo_defined
integer, parameter :: ndet = ndet_defined
integer, parameter :: ncar = ncar_defined
integer, parameter :: ntyp = nnut+nchl+nzoo+ndet+ncar
integer, parameter :: npe = nnut+nchl
integer, parameter :: nds = nnut+nchl+nzoo+1
integer, parameter :: nde = nds+(ndet-1)
integer, parameter :: ncs = ntyp-(ncar-1)

end module oasim_constants_mod
