program oasim_test

use oasim_mod, only: oasim

! Oasim oject
type(oasim) :: oasim_obj


! Write message
print*, 'Running oasim test'

! Run oasim constructor
call oasim_obj%create('data')

end program oasim_test
