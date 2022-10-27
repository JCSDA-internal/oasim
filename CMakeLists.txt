# (C) Copyright 2022 UCAR.
# (C) Copyright 2022 United States Government as represented by the Administrator of the
#     National Aeronautics and Space Administration. All Rights Reserved.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Preamble
# --------
cmake_minimum_required( VERSION 3.12 )
project( oasim VERSION 1.0.0 LANGUAGES Fortran )

## Ecbuild integration
set( ECBUILD_DEFAULT_BUILD_TYPE Release )
find_package( ecbuild 3.3.2 REQUIRED )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Dependencies
find_package( jedicmake QUIET )  # Prefer find modules from jedi-cmake

## Sources
add_subdirectory( src )
add_subdirectory( test )

## Export package config
ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
