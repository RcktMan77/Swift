module grid_properties
    use kind_defs,  only : dp
    use grid_types, only : grid, point, edge

    implicit none

    private
    public :: read_grid, extract_edges, calculate_edge_normals,                &
              calculate_element_areas, identify_neighbors,                     &
              calculate_edge_lengths, update_dx_array

contains


! The read_grid subroutine parses a given grid file for relevant grid
! information including interior element connectivities, node coordinates,
! and boundaries. This grid information is stored in a grid object with
! derived data types for points, elements, & boundaries. Currently only 2-D
! SU2-formatted grids are supported.

subroutine read_grid( filename, mesh )
    implicit none

    character(len=*), intent(in) :: filename

    type(grid), intent(out) :: mesh

    integer, allocatable :: int_elem_types(:), int_elem_data(:,:),             &
                            int_elem_ids(:), bndry_elem_types(:),              &
                            bndry_elem_data(:,:), temp_types(:), temp_data(:,:)

    integer :: unit_num, ios, i, j, num_dim, num_pts, num_int_elems,           &
               num_bndrys, num_bndry_elems, num_elems, bndry_elem_count

    logical :: found_ndime, found_npoin, found_nelem, found_nmark

    character(len=80) :: line

    character(len=20) :: bc_name, keyword

    ! Initialize flags & counters
    found_ndime = .false. ; found_npoin = .false.
    found_nelem = .false. ; found_nmark = .false.

    num_dim = 0 ; num_pts = 0 ;  num_bndrys = 0
    num_elems = 0 ; bndry_elem_count = 0

    ! Open SU2 formatted grid file
    open( newunit=unit_num, file=filename, status='old', action='read',        &
          iostat=ios )

    if ( ios /= 0 ) then
        print *, 'Error opening file: ', filename
        stop
    end if

    ! Parse the grid file for the grid dimension, number of points,
    ! number of elements,, number of boundaries, & number of elements for
    ! each boundary.

    do
        read(unit_num, '(A)', iostat=ios) line

        if ( ios /=0 ) exit

        if ( index(line, 'NDIME=') > 0 ) then
            read(line, *) keyword, num_dim

            found_ndime = .true.
            
            if ( num_dim /= 2 ) then
                print *, 'Error: Invalid grid dimension. Only 2-D grids are ', &
                         'supported.'
                stop
            end if

            mesh%grid_dim = num_dim
        else if ( index(line, 'NPOIN=') > 0 ) then
            read(line, *) keyword, num_pts

            found_npoin = .true.

            mesh%num_pts = num_pts

            if ( .not. allocated(mesh%points) ) then
                allocate( mesh%points(1:num_pts) )
            end if

            i = 1
            do while ( i <= num_pts )
                read(unit_num, '(A)', iostat=ios) line

                if ( ios /= 0 ) stop 'Error reading node coordinates.'

                if ( len_trim(line) > 0 .and. index(line, '=') == 0 ) then
                    read(line, *, iostat=ios) mesh%points(i)%x,                &
                                              mesh%points(i)%y,                &
                                              mesh%points(i)%node_id
                    if ( ios == 0 ) then
                        mesh%points(i)%node_id = i  ! Shift to 1-based index
                        i = i + 1
                    else
                        print *, 'Error reading node coordinates.'
                        stop
                    end if
                end if
            end do
        else if ( index(line, 'NELEM=') > 0 ) then
            read(line, *) keyword, num_int_elems

            found_nelem = .true.

            num_elems = num_int_elems

            allocate( int_elem_types(num_int_elems) )
            allocate( int_elem_data(4,num_int_elems) )
            allocate( int_elem_ids(num_int_elems) )

            int_elem_data = 0

            i = 1
            do while ( i <= num_int_elems )
                read(unit_num, '(A)', iostat=ios) line

                if ( ios /= 0 ) stop 'Error reading element connectivities.'

                if ( len_trim(line) > 0 .and. index(line, '=') == 0 ) then
                    read(line, *) int_elem_types(i)

                    select case ( int_elem_types(i) )
                    case (5) ! Triangle
                        read(line, *) int_elem_types(i),                       &
                                      int_elem_data(1:3,i),                    &
                                      int_elem_ids(i)
                    case (9) ! Quadrilateral
                        read(line, *) int_elem_types(i),                       &
                                      int_elem_data(1:4,i),                    &
                                      int_elem_ids(i)
                    case default
                        stop 'Error: Unsupported interior element type.'
                    end select

                    i = i + 1
                end if
            end do
        else if ( index(line, 'NMARK=') > 0 ) then
            read(line, *) keyword, num_bndrys

            found_nmark = .true.

            mesh%num_bndrys = num_bndrys

            if ( .not. allocated(mesh%bndrys) ) then
                allocate( mesh%bndrys(num_bndrys) )
            end if

            i = 1
            do while ( i <= num_bndrys )
                read(unit_num, '(A)', iostat=ios) line

                if ( ios /= 0 ) stop 'Error: MARKER_TAG= not found.'

                if ( index(line, 'MARKER_TAG=') > 0 ) then
                    read(line, *) keyword, bc_name

                    mesh%bndrys(i)%name = bc_name

                    read(unit_num, '(A)', iostat=ios) line

                    if ( ios /= 0 .or. index(line, 'MARKER_ELEMS=') == 0 ) then
                        print *, 'Error reading number of boundary elements ', &
                                 'for boundary: ', bc_name
                        stop
                    end if

                    read(line, *) keyword, num_bndry_elems

                    mesh%bndrys(i)%num_elems = num_bndry_elems

                    num_elems = num_elems + num_bndry_elems

                    if ( i == 1 ) then
                        allocate( bndry_elem_types(num_bndry_elems) )
                        allocate( bndry_elem_data(2,num_bndry_elems) )
                    else
                        call move_alloc( bndry_elem_types, temp_types )
                        call move_alloc( bndry_elem_data, temp_data )

                        allocate( bndry_elem_types(bndry_elem_count +          &
                                                   num_bndry_elems) )
                        allocate( bndry_elem_data(2,bndry_elem_count +         &
                                                  num_bndry_elems) )

                        bndry_elem_types(1:bndry_elem_count) = temp_types(     &
                            1:bndry_elem_count)
                        bndry_elem_data(:,1:bndry_elem_count) = temp_data(:,   &
                            1:bndry_elem_count)

                        deallocate( temp_types, temp_data )
                    end if

                    do j = 1, num_bndry_elems
                        read(unit_num, '(A)', iostat=ios) line

                        if ( ios /= 0 ) then
                            print *, 'Error reading boundary element ',        &
                                     'connectivities for boundary: ', bc_name
                            stop
                        end if

                        if ( len_trim(line) > 0 .and. index(line, '=') == 0 )  &
                            then
                            read(line, *) bndry_elem_types(bndry_elem_count    &
                                + j), bndry_elem_data(1:2,bndry_elem_count + j)
                        end if
                    end do

                    bndry_elem_count = bndry_elem_count + num_bndry_elems

                    i = i + 1
                end if
            end do
        end if
    end do

    if ( .not. found_ndime ) stop 'NDIME= not found.'
    if ( .not. found_npoin ) stop 'NPOIN= not found.'
    if ( .not. found_nelem ) stop 'NELEM= not found.'
    if ( .not. found_nmark ) stop 'NMARK= not found.'

    ! Populate all elements (interior & boundary)
    mesh%num_elems = num_elems

    if ( .not. allocated( mesh%elems ) ) then
        allocate( mesh%elems(1:num_elems) )
    end if

    mesh%num_tris = 0
    mesh%num_quads = 0

    ! Interior elements
    do i = 1, num_int_elems
        select case ( int_elem_types(i) )
            case ( 5 )
                mesh%elems(i)%num_nodes = 3
                mesh%elems(i)%num_edges = 3
                mesh%num_tris = mesh%num_tris + 1
            case ( 9 )
                mesh%elems(i)%num_nodes = 4
                mesh%elems(i)%num_edges = 4
                mesh%num_quads = mesh%num_quads + 1
            case default
                print *, 'Unsupported interior element type: ',                &
                    int_elem_types(i)
                stop
        end select

        allocate( mesh%elems(i)%node_ids(1:mesh%elems(i)%num_nodes) )
        allocate( mesh%elems(i)%edge_ids(1:mesh%elems(i)%num_edges) )
        allocate( mesh%elems(i)%edge_normals(2, mesh%elems(i)%num_edges) )
        allocate( mesh%elems(i)%neighbors(1:mesh%elems(i)%num_edges) )

        ! Shift node_ids to 1-based
        mesh%elems(i)%node_ids = int_elem_data(1:mesh%elems(i)%num_nodes,i) + 1
        mesh%elems(i)%edge_ids = 0
        mesh%elems(i)%is_bndry = .false.
        mesh%elems(i)%bndry_id = 0
    end do

    ! Boundary elements
    do i = 1, bndry_elem_count
        select case ( bndry_elem_types(i) )
            case ( 3 )
                mesh%elems(num_int_elems + i )%num_nodes = 2
                mesh%elems(num_int_elems + i )%num_edges = 1
            case default
                print *, 'Unsupported boundary element type: ',                &
                    bndry_elem_types(i)
                stop
        end select

        allocate( mesh%elems(num_int_elems + i)%node_ids(1:2) )
        allocate( mesh%elems(num_int_elems + i)%edge_ids(1:1) )
        allocate( mesh%elems(num_int_elems + i)%edge_normals(2,1) )
        allocate( mesh%elems(num_int_elems + i)%neighbors(1:1) )

        ! Shift node_ids to 1-based
        mesh%elems(num_int_elems + i)%node_ids = bndry_elem_data(1:2,i) + 1
        mesh%elems(num_int_elems + i)%edge_ids = 0
        mesh%elems(num_int_elems + i)%is_bndry = .true.
        mesh%elems(num_int_elems + i)%bndry_id = 0
    end do

    deallocate( int_elem_types, int_elem_data, int_elem_ids )

    if ( allocated( bndry_elem_types ) ) then
        deallocate( bndry_elem_types, bndry_elem_data )
    end if

    close( unit_num )
end subroutine read_grid


! The extract_edges subroutine processes interior & boundary elements to
! extract and index unique edges for each element within a given grid. Edges
! are indexed using a sparse hash table. An edge_id (a unique, order-
! independent key for edges) is generated for each pair of nodes, and added
! to the edges array attribute for each element.

! The element IDs shared by each unique edge are then added to the elem_ids
! attribute for each edge.

subroutine extract_edges( mesh )
    implicit none

    type(grid), intent(inout) :: mesh

    integer :: num_edges, n1, n2, edge_id, i, j, k

    integer, allocatable :: edge_map(:,:)

    ! Map SU2's 0-based node_ids to Fortran's native 1-based indexing
    ! for arrays.
    allocate( edge_map(1:mesh%num_pts, 1:mesh%num_pts) )
    edge_map = 0
    num_edges = 0

    ! Note: The mod( j, 3 ) + 1 statement in the following do loop is used
    ! to cycle through the nodes of a triangle in a circular manner.
    ! mod( j, 3 ) calculates the remainder when j is divided by 3. This
    ! operation always results in 0, 1, or 2. Adding 1 to the result gives
    ! 1, 2, or 3. The formula is particularly useful when dealing with
    ! triangular elements because:

    ! when j = 1, mod( 1, 3 ) + 1 = 2,
    ! when j = 2, mod( 2, 3 ) + 1 = 3, &
    ! when j = 3, mod( 3, 3 ) + 1 = 1

    ! The cycling behavior is perfect for connecting the nodes of a triangle
    ! in order, including the wrap-around from the last node back to the
    ! first. It allows the program to:

    ! connect node 1 to node 2,
    ! connect node 2 to node 3, &
    ! connect node 3 back to node 1.

    ! This technique is commonly used in mesh-related algorithms to efficiently
    ! handle the circular nature of element connectivity without using
    ! conditional statements. This technique is also adapted for quadrilateral
    ! elements.


    ! Note: An edge is defined by two node IDs (e.g., n1 = 0, n2 = 1
    ! or n1 = 1, n2 = 0 for the same edge 0-1)

    ! In a mesh, edges are unidirected (i.e., 0-1 is the same as 1-0). Without
    ! consistent ordering, the same edge might be counted twice or assigned 
    ! different IDs depending on the element's ordering (e.g., clockwise vs.
    ! counterclockwise).

    ! Using min(n1,n2) and max(n1,n2) ensures a unique, order-independent key:

    ! For n1 = 0, n2 = 1: min(0,1) = 0, max(0,1) = 1 -> edge_map(0,1).
    ! For n1 = 1, n2 = 0: min(1,0) = 0, max(1,0) = 1 -> edge_map(0,1).

    ! The canonical form (smaller_node,larger_node) guarantees that each edge
    ! maps to a single location in edge_map.

    ! Count unique edges from all elements
    do i = 1, mesh%num_elems
        do j = 1, mesh%elems(i)%num_edges
            n1 = mesh%elems(i)%node_ids(j)
            n2 = mesh%elems(i)%node_ids(mod(j,mesh%elems(i)%num_nodes) + 1)

            if ( edge_map(min(n1,n2),max(n1,n2)) == 0 ) then
                num_edges = num_edges + 1
                edge_map(min(n1,n2), max(n1,n2)) = num_edges
            end if
        end do
    end do

    ! Allocate edges array
    mesh%num_edges = num_edges

    if ( .not. allocated( mesh%edges ) ) then
        allocate( mesh%edges(1:num_edges) )
    end if

    do i = 1, num_edges
        allocate( mesh%edges(i)%elem_ids(1:2) )
        mesh%edges(i)%elem_ids = -1
    end do
    
    ! Populate edges & edge_ids from interior elements only
    !$omp parallel do private(i, j, n1, n2, edge_id, k)
    do i = 1, mesh%num_elems
        do j = 1, mesh%elems(i)%num_edges
            n1 = mesh%elems(i)%node_ids(j)
            n2 = mesh%elems(i)%node_ids(mod(j,mesh%elems(i)%num_nodes) + 1)

            edge_id = edge_map(min(n1,n2), max(n1,n2))

            mesh%elems(i)%edge_ids(j) = edge_id
            !$omp critical
            mesh%edges(edge_id)%node_ids = [min(n1,n2), max(n1,n2)]

            if ( .not. mesh%elems(i)%is_bndry ) then
                k = count( mesh%edges(edge_id)%elem_ids >= 0 )
                
                ! 1-based interior elem_id
                if ( k < 2 ) then
                    mesh%edges(edge_id)%elem_ids(k + 1) = i
                end if
            end if
            !$omp end critical
        end do
    end do
    !$omp end parallel do

    ! Count interior & boundary edges
    mesh%num_interior_edges = 0
    mesh%num_bndry_edges = 0

    do i = 1, mesh%num_edges
        if ( count(mesh%edges(i)%elem_ids >= 0) == 2 ) then
            mesh%num_interior_edges = mesh%num_interior_edges + 1
        else if ( count(mesh%edges(i)%elem_ids >=0) ==1 ) then
            mesh%num_bndry_edges = mesh%num_bndry_edges + 1
        end if
    end do

    deallocate( edge_map )
end subroutine extract_edges


! The calculate_edge_normals subroutine loops ovr all of the edges for each
! element and computes the edge normal direction. The orientation for each
! edge normal is compared to the element's centroid to ensure the edge normal
! direction is outward facing for each element, and finally the normal
! vector is stored in the edge_normals array attribute for each element.

subroutine calculate_edge_normals( mesh )
    implicit none

    type(grid), intent(inout) :: mesh

    integer :: n1, n2, edge_id, i, j

    real(dp) :: dx, dy, nx, ny, cx, cy, edge_length, dot_prod

    real(dp), parameter :: TOL = 1e-10_dp   ! Tolerance for small edge lengths

    do i = 1, mesh%num_elems
        ! Calculate the centroid for the given element
        cx = sum( mesh%points(mesh%elems(i)%node_ids)%x ) /                    &
            mesh%elems(i)%num_nodes

        cy = sum( mesh%points(mesh%elems(i)%node_ids)%y ) /                    &
            mesh%elems(i)%num_nodes

        do j = 1, mesh%elems(i)%num_edges
            edge_id = mesh%elems(i)%edge_ids(j)

            if ( edge_id == 0 ) then
                ! Handle uninitialized edge_id
                nx = 0.0_dp
                ny = 0.0_dp

                print *, 'Warning: Unassinged edge_id at element: ',           &
                    i, ' edge: ', j
            else
                n1 = mesh%edges(edge_id)%node_ids(1)
                n2 = mesh%edges(edge_id)%node_ids(2)

                dx = mesh%points(n2)%x - mesh%points(n1)%x
                dy = mesh%points(n2)%y - mesh%points(n1)%y

                edge_length = dx**2 + dy**2 

                ! Check for small dx & dy before calculating edge length
                if ( edge_length < TOL ) then
                    nx = 0.0_dp
                    ny = 0.0_dp

                    print *, 'Warning: Degenerate edge detected at edge_id: ', &
                        edge_id, ' nodes: ', n1, n2
                else
                    nx = dy
                    ny = -dx
                    nx = nx / sqrt( nx**2 + ny**2 )
                    ny = ny / sqrt( nx**2 + ny**2 )

                    ! Check whether the normal points towards the centroid & 
                    ! if so, then reverse the normal direction
                    dot_prod = ( nx * ( mesh%points(n1 + 1)%x - cx ) ) +       &
                               ( ny * ( mesh%points(n1 + 1)%y - cy ) )

                    if ( dot_prod < 0 ) then
                        nx = -nx
                        ny = -ny
                    end if
                end if
            end if

            mesh%elems(i)%edge_normals(:,j) = [nx, ny]
        end do
    end do
end subroutine calculate_edge_normals


! The calculate_element_areas subroutine calculates the area for each element
! within a given grid derived type object.

subroutine calculate_element_areas( mesh )
    implicit none

    type(grid), intent(inout) :: mesh

    integer :: i, j

    real(dp) :: area

    ! For each element, the area is calculated using the Shoelace (Gauss's
    ! area) formula via cross terms:

    ! x_j * y{j+1} - (y_j * x_{j+1})

    ! mod(j, num_nodes) + 1 wraps around to the first node (e.g., j = 3 -> 1
    ! for a triangle), closing the polygon.

    ! abs(area) * 0.5_dp converts the doubled signed area to the actual
    ! positive area

    do i = 1, mesh%num_elems
        area = 0.0_dp
        do j = 1, mesh%elems(i)%num_nodes
            area = area + ( mesh%points(mesh%elems(i)%node_ids(j))%x *         &
                mesh%points(mesh%elems(i)%node_ids(mod(j,                      &
                mesh%elems(i)%num_nodes) + 1))%y ) -                           &
                ( mesh%points(mesh%elems(i)%node_ids(j))%y *                   &
                mesh%points(mesh%elems(i)%node_ids(mod(j,                      &
                mesh%elems(i)%num_nodes) + 1))%x )
        end do

        mesh%elems(i)%area = abs( area ) * 0.5_dp
    end do
end subroutine calculate_element_areas


! The identify_neighbors subroutine loops over each edge in a given grid 
! type derived object. In 2-D, each interior edge is shared by two unique
! elements. For each element, the subroutine then loops over each of the 
! element's edges. If one of the element's edges matches the current
! edge_id, then the other element shared with the edge is added to the
! current element's neighbors array. Then this is repeated for the second
! element that shares the edge_id.

subroutine identify_neighbors( mesh )
    implicit none

    type(grid), intent(inout) :: mesh

    integer :: edge_id, elem1, elem2, i, j

    do i = 1, mesh%num_edges
        if ( count(mesh%edges(i)%elem_ids >= 0) == 2 ) then
            elem1 = mesh%edges(i)%elem_ids(1)
            elem2 = mesh%edges(i)%elem_ids(2)

            do j = 1, mesh%elems(elem1)%num_edges
                if ( mesh%elems(elem1)%edge_ids(j) == i ) then
                    mesh%elems(elem1)%neighbors(j) = elem2
                end if
            end do

            do j = 1, mesh%elems(elem2)%num_edges
                if ( mesh%elems(elem2)%edge_ids(j) == i ) then
                    mesh%elems(elem2)%neighbors(j) = elem1
                end if
            end do
        end if
    end do
end subroutine identify_neighbors


! The calculate_edge_lengths subroutine loops over each edge within a given
! grid derived type object and calculates the Euclidean distance between
! node pairs corresponding to each edge. The edge length is then stored as
! an attribute for each edge.

subroutine calculate_edge_lengths( mesh )
    implicit none

    type(grid), intent(inout) :: mesh

    integer :: n1, n2, i

    real(dp) :: dx, dy

    do i = 1, mesh%num_edges
        n1 = mesh%edges(i)%node_ids(1)
        n2 = mesh%edges(i)%node_ids(2)

        dx = mesh%points(n2)%x - mesh%points(n1)%x
        dy = mesh%points(n2)%y - mesh%points(n1)%y

        mesh%edges(i)%length = sqrt( dx**2 + dy**2 )
    end do
end subroutine calculate_edge_lengths


! The update_dx_array subroutine updates the dx array attribute for a given
! grid derived type object with an area-weighted characteristic length for
! each node in the mesh. At each node, the sum of the edge lengths for edges
! that share a given node are calculated and weighted by the areas of elements
! that share the edge. These area-weighted characteristic length values are
! used downstream to calculate a time step using the velocity at each node
! and CFL parameter for each solver iteration.

subroutine update_dx_array( mesh )
    implicit none

    type(grid), intent(inout) :: mesh

    integer :: edge_id, i, j, k, n

    real(dp) :: total_length, total_area

    logical, allocatable :: edge_has_node(:)

    if ( .not. allocated( mesh%dx ) ) then
        allocate( mesh%dx(mesh%num_pts) )
        mesh%dx = 0.0_dp
    end if

    !$omp parallel                                                             &
    !$omp private(n, i, j, k, edge_id)                                         &
    !$omp private(total_length, total_area, edge_has_node)
    allocate( edge_has_node(mesh%num_edges) )

    !$omp do
    do n = 1, mesh%num_pts
        total_length = 0.0_dp
        total_area = 0.0_dp

        do i = 1, mesh%num_edges
            edge_id = i
            edge_has_node(i) = any( mesh%edges(edge_id)%node_ids == n )

            if ( edge_has_node(i) ) then
                total_length = total_length + mesh%edges(edge_id)%length

                do j = 1, count( mesh%edges(edge_id)%elem_ids >= 0 )
                    k = mesh%edges(edge_id)%elem_ids(j)

                    if ( any(mesh%elems(k)%node_ids == n) ) then
                        total_area = total_area + mesh%elems(k)%area
                    end if
                end do
            end if
        end do

        if ( total_area > 0.0_dp ) then
            mesh%dx(n) = total_length * total_area /                           &
                real( count(edge_has_node), dp )
        else
            mesh%dx(n) = 0.0_dp
        end if
    end do
    !$omp end do
    deallocate( edge_has_node )
    !$omp end parallel
end subroutine update_dx_array

end module grid_properties
