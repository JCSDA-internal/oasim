# (C) Copyright 2022 UCAR.
# (C) Copyright 2022 United States Government as represented by the Administrator of the
#     National Aeronautics and Space Administration. All Rights Reserved.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Compile executable
ecbuild_add_executable( TARGET  test_oasim.x
                        SOURCES test_oasim.f90
                        LIBS    oasim )

# Add a test
ecbuild_add_test( TARGET   oasim_test
                  COMMAND  test_oasim.x )
