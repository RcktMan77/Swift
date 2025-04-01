program main
    use grid_types, only : grid
    use flow_types, only : primitive_state
    use namelist_definitions
    use initialization
    use grid_properties
    use dual_mesh
    use write_output

    implicit none

    integer :: i, j

    type(grid) :: mesh

    type(primitive_state), allocatable :: w(:)

    character(len=20) :: grid_filename, input_filename, output_filename

    real :: start_time, end_time

    ! Specify the grid file to read
    grid_filename = 'naca0012.su2'

    ! Specify the input file to read
    input_filename = 'swift.nml'

    ! Specify the solution output file
    output_filename = 'swift_tec_volume'

    ! Initialize number of edges within the grid
    mesh%num_edges = 0

    ! Read the grid file
    call cpu_time( start_time )
    call read_grid( grid_filename, mesh )
    call cpu_time( end_time )

    print *, 'Time Elapsed: '
    print *, 'Read Grid: ', end_time - start_time

    ! TODO: Ideally the program should check that the namelist input file
    ! exists and verify that all parameters being passed into the program
    ! have valid input values. This should be done prior to reading the grid
    ! file in order to terminate the program more quickly if needed. 

    ! Read the input file
    call cpu_time( start_time )
    call check_and_read_namelist( mesh, input_filename )
    call cpu_time( end_time )

    print *, 'Read Input File: ', end_time - start_time

    ! Initialize the flow field
    call cpu_time( start_time )
    call initialize_flow( mesh, w )
    call cpu_time( end_time )

    print *, 'Initialize Flow Field: ', end_time - start_time

    ! Process elements to extract edges
    call cpu_time( start_time )
    call extract_edges( mesh )
    call cpu_time( end_time )

    print *, 'Extract Edges: ', end_time - start_time

    ! Calculate edge normal vectors & orient the edge normals for each
    ! edge of an element correctly (facing outward) based on its computed
    ! direction relative to the element's centroid.
    call cpu_time( start_time )
    call calculate_edge_normals( mesh )
    call cpu_time( end_time )

    print *, 'Calculate Edge Normal Vectors: ', end_time - start_time

    ! Process elements to calculate element areas
    call cpu_time( start_time )
    call calculate_element_areas( mesh )
    call cpu_time( end_time )

    print *, 'Calculate Element Areas: ', end_time - start_time

    ! Process edges to calculate edge lengths
    call cpu_time( start_time )
    call calculate_edge_lengths( mesh )
    call cpu_time( end_time )

    print *, 'Calculate Edge Lengths: ', end_time - start_time

    ! Identify element neighbors for each element
    call cpu_time( start_time )
    call identify_neighbors( mesh )
    call cpu_time( end_time )

    print *, 'Identify Neighbors: ', end_time - start_time

    ! Generate median dual mesh
    call cpu_time( start_time )
    call create_dual_mesh( mesh )
    call cpu_time( end_time )

    print *, 'Generate Median Dual Mesh: ', end_time - start_time

    ! Apply boundary conditions
!   call cpu_time( start_time )
!   call apply_boundary_conditions( mesh, q )
!   call cpu_time( end_time )

!   print *, 'Apply Boundary Conditions: ', end_time - start_time

    ! Calculate element area-weighted characteristic length at each node
    call cpu_time( start_time )
    call update_dx_array( mesh )
    call cpu_time( end_time )

    print *, 'Calculate Characteristic Lengths: ', end_time - start_time

    ! Write output
    call cpu_time( start_time )
    call write_szplt( mesh, w, output_filename )
    call cpu_time( end_time )

    print *, 'Write output: ', end_time - start_time

    ! Output for verification
    print *, ' '
    print *, 'Grid Dimension:            ', mesh%grid_dim
    print *, 'Number of points:          ', mesh%num_pts
    print *, 'Number of triangles:       ', mesh%num_tris
    print *, 'Number of quadrilaterals:  ', mesh%num_quads
    print *, 'Number of elements:        ', mesh%num_elems
    print *, 'Number of boundaries:      ', mesh%num_bndrys
    print *, 'Number of edges:           ', mesh%num_edges
    
    print *, ' '
    
    print *, 'Boundary 1:                ', mesh%bndrys(1)%name
    print *, 'Number of elements:        ', mesh%bndrys(1)%num_elems
 
    print *, ' '
 
    print *, 'Boundary 2:                ', mesh%bndrys(2)%name
    print *, 'Number of elements:        ', mesh%bndrys(2)%num_elems
 
    print *, ' '

    print *, 'Element 1:'
    
    do j = 1, mesh%elems(1)%num_nodes
        write(*, '(A, I0, A, I5, A, F12.6, A, F12.6)')                         &
            ' Node ', j, ': node_id = ', mesh%elems(1)%node_ids(j),            &
            ' x = ', mesh%points(mesh%elems(1)%node_ids(j))%x,                 &
            ' y = ', mesh%points(mesh%elems(1)%node_ids(j))%y                 
    end do

    print *, ' '

    do j = 1, mesh%elems(1)%num_edges
        write(*, '(A, I0, A, F12.6, A, F12.6, A, F12.6)')                      &
            ' Edge ', j, ': nx = ', mesh%elems(1)%edge_normals(1,j),           &
            ' ny = ', mesh%elems(1)%edge_normals(2,j),                         &
            ' length = ', mesh%edges(mesh%elems(1)%node_ids(j))%length
    end do

    print *, ' '

    print *, 'Area: ', mesh%elems(1)%area
    print *, 'Element 1 has neighbors: ', mesh%elems(1)%neighbors
 
    print *, ' '

    write(*,'("Reynolds Number:    ", ES12.6E2)') reynolds_number
    write(*,'("Mach Number:        ", ES12.6E2)') mach_number
    write(*,'("Temperature:        ", ES12.6E2)') temperature
    write(*,'("Inviscid Flux Method:  ", A)') flux_construction

    print *, ' '

    do i = 1, mesh%num_bndrys
        print *, 'Boundary ', i, ':'
        print *, '   Name: ', trim( mesh%bndrys(i)%name )
        print *, '   Type: ', mesh%bndrys(i)%bc_type
    end do

    print *, ' '

    print*, 'Dual Face 1'
    write(*, '(A, F8.6, A, F8.6)')                                             &
        ' Start: x = ', mesh%dual_faces(1)%start%x,                            &
        ' y = ', mesh%dual_faces(1)%start%y
    write(*, '(A, F8.6, A, F8.6)')                                             &
        ' End:   x = ', mesh%dual_faces(1)%end%x,                              &
        ' y = ', mesh%dual_faces(1)%end%y

contains


end program main
