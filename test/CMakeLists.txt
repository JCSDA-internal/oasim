# (C) Copyright 2022 UCAR.
# (C) Copyright 2022 United States Government as represented by the Administrator of the
#     National Aeronautics and Space Administration. All Rights Reserved.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


# Macro to link list of files from source to destination
macro( LINK_FILES filelist src_dir dst_dir )
  foreach(FILENAME ${filelist})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${src_dir}/${FILENAME}
      ${dst_dir}/${FILENAME}
    )
  endforeach(FILENAME)
endmacro()

# List of data files
list( APPEND data_files
  data/abw25_morel.dat
  data/acbc25_6.dat
  data/atmo25b.dat
  data/pic_sigma.dat
  data/slingo.dat
)

# Create a directory to hold data files in the testing context
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data/)

# Link data to test directory
LINK_FILES( "${data_files}" ${CMAKE_CURRENT_SOURCE_DIR}/../ ${CMAKE_CURRENT_BINARY_DIR})

# Compile executable
ecbuild_add_executable( TARGET  test_oasim.x
                        SOURCES test_oasim.f90
                        LIBS    oasim )

# Add the test
ecbuild_add_test( TARGET   oasim_test
                  COMMAND  test_oasim.x )
