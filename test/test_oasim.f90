program oasim_test

use oasim_mod, only: oasim

! Oasim oject
type(oasim) :: oasim_obj

! Data directory
character(len=2048) :: data_directory

! Write message
print*, 'Running oasim test'

! Run oasim constructor
call oasim_obj%create(data_directory)

end program oasim_test
