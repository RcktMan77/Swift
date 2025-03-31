module write_output
    use iso_c_binding,        only : c_ptr, c_null_ptr, c_char, c_int32_t,     &
                                     c_int64_t, c_double
    use kind_defs,            only : dp
    use grid_types,           only : grid
    use flow_types,           only : primitive_state
    use namelist_definitions, only : turbulence_model
    use initialization,       only : nturb

    implicit none

    private
    public :: write_szplt

    ! Include the TecIO interface
    include 'tecio.f90'

contains

subroutine write_szplt( mesh, w, filename )
    implicit none

    type(grid), intent(in) :: mesh

    type(primitive_state), intent(in) :: w(:)

    character(len=*), intent(in) :: filename

    type(c_ptr) :: file_handle = c_null_ptr

    integer(c_int64_t) :: num_nodes, total_nodes, num_elems
    integer(c_int64_t) :: num_face_connections = 0_c_int64_t

    integer(c_int32_t) :: file_format         = 1  ! 1 = SZPLT
    integer(c_int32_t) :: file_type           = 0  ! 0 = Full file 
    integer(c_int32_t) :: default_var_type    = 2  ! Double precison data
    integer(c_int32_t) :: share_connectivity  = 0
    integer(c_int32_t) :: face_neighbor_mode  = 0
    integer(c_int32_t) :: nodes_are_one_based = 0  ! Inherited SU2 0-based
!   integer(c_int32_t) :: zone_type        = 5  ! FEMixed 
    integer(c_int32_t) :: zone_type        = 3  ! FEQuadrilateral

    integer(c_int32_t) :: num_tris, num_quads, num_vars, ierr, zone

    integer(c_int64_t), allocatable :: connectivity(:)

    integer(c_int32_t), allocatable :: var_types(:), share_var(:),             &
                                       value_loc(:), passive_var(:)

    real(dp), allocatable :: var_data(:)

    character(len=32) :: title
    
    character(len=32), allocatable :: var_names(:)

    integer :: i, offset, elem_count


    ! Grid & state info
    num_nodes = int(mesh%num_pts, c_int64_t)    ! Convert to 64-bit
    num_elems = int(mesh%num_tris + mesh%num_quads, c_int64_t)
    num_tris  = mesh%num_tris
    num_quads = mesh%num_quads
    num_vars  = 6 + nturb

    ! Calculate total connectivity nodes (3 per triangle, 4 per quadrilateral)
    total_nodes = 3_c_int64_t * num_tris + 4_c_int64_t * num_quads
   
    ! Variable names
    title = 'Flow Solution'
    allocate( var_names(num_vars) )
    var_names(1:6) = ['x  ', 'y  ', 'rho', 'u  ', 'v  ', 'p  ']
    
    if ( nturb == 1 ) then 
        var_names(7) = 'nu '
    else if ( nturb == 2 ) then
        var_names(7) = 'k    '
        var_names(8) = 'omega'
    end if

    ! Open the .szplt file
    ierr = tecFileWriterOpen(trim(filename)//char(0),                          &
                             trim(title)//char(0),                             &
                             trim(adjustl(join_var_names(var_names)))//char(0),&
                             file_format, file_type, default_var_type,         &
                             c_null_ptr, file_handle)

    if ( ierr /= 0 ) stop 'Error opening TecIO file.'

    ! Enable diagnostics
    ierr = tecFileSetDiagnosticsLevel( file_handle, 1 )
    if ( ierr /= 0 ) then
        write(*,*) 'Warning: Could not set diagnostics level, ierr = ', ierr
    end if

    ! Allocate arrays
    allocate( connectivity(total_nodes) )
    allocate( var_types(num_vars), share_var(num_vars),                        &
              value_loc(num_vars), passive_var(num_vars) )
    allocate( var_data(num_nodes) )

    ! Define FEMixed zone
    var_types   = 2     ! Double precision for all variables
    share_var   = 0     ! Not sharing
    value_loc   = 1     ! Node-centered
    passive_var = 0     ! All active

    ! Populate connectivity (0-based indexing for Tecplot)
    offset = 0
    elem_count = 0
    do i = 1, mesh%num_elems
        if ( mesh%elems(i)%num_nodes == 2 ) then
            cycle   ! Skip boundary elements
        else if ( mesh%elems(i)%num_nodes == 3 ) then
             connectivity(offset + 1:offset + 3) =                             &
                mesh%elems(i)%node_ids(1:3) - 1

             offset = offset + 3
             elem_count = elem_count + 1
        else if ( mesh%elems(i)%num_nodes == 4 ) then
             connectivity(offset + 1:offset +4) =                              &
                 mesh%elems(i)%node_ids(1:4) - 1

             offset = offset + 4
             elem_count = elem_count + 1
         else
             write(*,*) 'Unsupported element type with ',                      &
                 mesh%elems(i)%num_nodes, ' nodes.'
             stop 'Error in mesh data.'
        end if
    end do

    ! Debug output
    if ( offset /= total_nodes ) then
        write(*,*) 'Connectivity mismatch: expected ', total_nodes, ', got ',  &
            offset
        stop 'Error in connectivity.'
    end if

    if ( elem_count /= num_elems ) then
        write(*,*) 'Element count mismatch: expected ', num_elems, ' got ',    &
            elem_count
        stop 'Error in element count.'
    end if

    ierr = tecZoneCreateFE( file_handle, 'Zone1'//char(0), zone_type,          &
                            num_nodes, num_elems, var_types, share_var,        &
                            value_loc, passive_var, share_connectivity,        &
                            num_face_connections, face_neighbor_mode, zone )

    if ( ierr /= 0 ) then
        write(*,*) 'tecZoneCreateFE failed with ierr = ', ierr
        stop 'Error writing FE zone.'
    end if

    ! Write connectivities
    ierr = tecZoneNodeMapWrite64( file_handle, zone, 0, nodes_are_one_based,   &
                                  total_nodes, connectivity )
    if ( ierr /= 0 ) then
        write(*,*) 'tecZoneNodeMapWrite64 failed with ierr = ', ierr
        stop 'Error writing node map.'
    end if


    ! Write node-centered data

    ! X-coordinate (x)
    do i = 1, num_nodes
        var_data(i) = mesh%points(i)%x
    end do
    ierr = tecZoneVarWriteDoubleValues( file_handle, zone, 1, 1,               &
                                        num_nodes, var_data )

    ! Y-coordinate (y)
    do i = 1, num_nodes
        var_data(i) = mesh%points(i)%y
    end do
    ierr = tecZoneVarWriteDoubleValues( file_handle, zone, 2, 1,               &
                                       num_nodes, var_data )

    ! Density (rho)
    do i = 1, num_nodes
        var_data(i) = w(i)%rho
    end do
    ierr = tecZoneVarWriteDoubleValues( file_handle, zone, 3, 1,               &
                                        num_nodes, var_data )

    ! X-velocity (u)
    do i = 1, num_nodes
        var_data(i) = w(i)%vel(1)
    end do
    ierr = tecZoneVarWriteDoubleValues( file_handle, zone, 4, 1,               &
                                        num_nodes, var_data )

    ! Y-velocity (v)
    do i = 1, num_nodes
        var_data(i) = w(i)%vel(2)
    end do
    ierr = tecZoneVarWriteDoubleValues( file_handle, zone, 5, 1,               &
                                        num_nodes, var_data )

    ! Pressure (p)
    do i = 1, num_nodes
        var_data(i) = w(i)%p
    end do
    ierr = tecZoneVarWriteDoubleValues( file_handle, zone, 6, 1,               &
                                        num_nodes, var_data )

    ! Turbulence variables
    if ( nturb >= 1 ) then
        do i = 1, num_nodes
            var_data(i) = w(i)%turb(1)      ! ν for SA, κ for SST
        end do
        ierr = tecZoneVarWriteDoubleValues( file_handle, zone, 7, 1,           &
                                            num_nodes, var_data )
    end if

    if ( nturb >= 2 ) then
        do i = 1, num_nodes
            var_data(i) = w(i)%turb(2)      ! ω for SST
        end do
        ierr = tecZoneVarWriteDoubleValues( file_handle, zone, 8, 1,           &
                                            num_nodes, var_data )
    end if

    if ( ierr /= 0 ) stop 'Error writing data.'

    ! Close the file
    ierr = tecFileWriterClose( file_handle )
    if ( ierr /= 0 ) stop 'Error closing TecIO file.'

    ! Clean up
    deallocate( connectivity, var_types, share_var,                            &
                value_loc, passive_var, var_data )


    write(*,*) 'Wrote grid & state to ', trim(filename)
end subroutine write_szplt


! The output variable string list passed to the tecFileWriterOpen procedure
! needs to be dynamically updated to accommodate a dynamic selection of 
! output variables due to turbulence model selection, dimension of simulation,
! etc. A dynamically allocated array containing the variable names is
! passed to the join_var_names function as an argument, and a concatenated
! string of output variable names suitable for the tecFileWriterOpen procedure 
! is then returned.

function join_var_names( var_names ) result( str )
    implicit none

    character(len=32), intent(in) :: var_names(:)

    character(len=256) :: str

    integer :: i

    str = trim( var_names(1) )

    do i = 2, size( var_names )
        str = trim( str ) // ',' // trim( var_names(i) )
    end do
end function join_var_names

end module write_output

