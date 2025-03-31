module kind_defs
    implicit none

    integer, parameter :: sp = selected_real_kind(6,37)   ! single precision
    integer, parameter :: dp = selected_real_kind(15,307) ! double precision

end module kind_defs

