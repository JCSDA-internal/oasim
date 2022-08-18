# (C) Copyright 2022 UCAR.
# (C) Copyright 2022 United States Government as represented by the Administrator of the
#     National Aeronautics and Space Administration. All Rights Reserved.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Include the current source directory in the path
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# List of source files to compile
list( APPEND oasim_src_files
  clrtrans_mod.f90
  daysetrad_mod.f90
  edeu_mod.f90
  glight_mod.f90
  lidata_mod.f90
  lidatatm_mod.f90
  oasim_constants_mod.f90
  oasim_mod.f90
  ocalbedo_mod.f90
  radmod_mod.f90
  setlte_mod.f90
  setsfclte_mod.f90
  sfcirr_mod.f90
  slingo_mod.f90
)

ecbuild_add_library( TARGET   ${PROJECT_NAME}
                     SOURCES  ${oasim_src_files}
                    )

# Fortran module output directory for build and install interfaces
set(MODULE_DIR module/${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${MODULE_DIR}>)

