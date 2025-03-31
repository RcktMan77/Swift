module dual_mesh
    use kind_defs,  only : dp
    use grid_types, only : point, edge, dual_face, element, grid

    implicit none

    private :: get_centroid, get_midpoint, get_normal_vector,                  &
               distance_between_points
    public :: create_dual_mesh

contains

! The get_centroid function returns the centroid for a given element 
! by computing the average x- and y-coordinates over the element's set
! of nodes.

function get_centroid( elem, mesh ) result( centroid )
    implicit none

    type(element), intent(in) :: elem
    
    type(grid), intent(in) :: mesh

    type(point) :: centroid

    integer :: i, ni

    centroid%x = 0.0_dp
    centroid%y = 0.0_dp

    do i = 1, elem%num_nodes
        ni = elem%node_ids(i)
        
        centroid%x = centroid%x + mesh%points(ni)%x
        centroid%y = centroid%y + mesh%points(ni)%y
    end do

    centroid%x = centroid%x / real( elem%num_nodes, dp )
    centroid%y = centroid%y / real( elem%num_nodes, dp )
end function get_centroid


! The get_midpoint function returns the midpoint for a given edge argument
! by calculating the half of the distance between the corresponding nodes
! for the edge.

function get_midpoint( primal_edge, mesh ) result( midpoint )
    implicit none

    type(edge), intent(in) :: primal_edge

    type(grid), intent(in) :: mesh

    type(point) :: midpoint

    integer :: n1, n2

    n1 = primal_edge%node_ids(1)
    n2 = primal_edge%node_ids(2)

    midpoint%x = ( mesh%points(n1)%x + mesh%points(n2)%x ) / 2.0_dp
    midpoint%y = ( mesh%points(n1)%y + mesh%points(n2)%y ) / 2.0_dp
end function get_midpoint


! The get_normal_vector function loops over each edge of an element for a
! given element ID, and if the current edge matches the edge_id argument,
! the normal vector for the edge is returned.

function get_normal_vector( elem_id, edge_id, mesh ) result( normal )
    implicit none

    integer, intent(in) :: elem_id, edge_id

    type(grid), intent(in) :: mesh

    real(dp) :: normal(2)

    integer :: j

    do j = 1, mesh%elems(elem_id)%num_edges
        if ( mesh%elems(elem_id)%edge_ids(j) == edge_id ) then
            normal = mesh%elems(elem_id)%edge_normals(:,j)
            return
        end if
    end do

    print *, 'Error: Edge ', edge_id, 'not found in element ', elem_id
    stop
end function get_normal_vector


! The distance_between_points function simply calulates & returns the
! distance between two given point arguments.

function distance_between_points( p1, p2 ) result( dist )
    implicit none

    type(point), intent(in) :: p1, p2

    real(dp) :: dist

    dist = sqrt( ( p2%x - p1%x )**2 + ( p2%y - p1%y )**2 )
end function distance_between_points


! The create_dual_mesh subroutine generates a median dual mesh where each
! primal edge is bisected at its midpoint, and the midpoints are connected
! to the centroids of adjacent elements. The input mesh argument has its
! dual_face attribute updated to store the dual mesh attributes.

subroutine create_dual_mesh( mesh )
    implicit none

    type(grid), intent(inout) :: mesh

    type(point) :: c1, c2, midpoint

    real(dp) :: normal(2), dist

    integer :: elem_id1, elem_id2, face_count, i

    face_count = 2 * mesh%num_interior_edges + mesh%num_bndry_edges

    if ( face_count <= 0 ) then
        print *, 'Error: No valid dual faces to create.'
        stop
    end if

    ! Allocate dual_faces
    if ( .not. allocated( mesh%dual_faces ) ) then
        allocate( mesh%dual_faces(face_count) )
    end if

    ! Populate dual_faces
    face_count = 0

    do i = 1, mesh%num_edges
        midpoint = get_midpoint( mesh%edges(i), mesh )

        if ( count(mesh%edges(i)%elem_ids >= 0) == 2 ) then
            ! Internal edge: two dual faces
            elem_id1 = mesh%edges(i)%elem_ids(1)
            elem_id2 = mesh%edges(i)%elem_ids(2)

            ! Face 1: midpoint to centroid 1
            face_count = face_count + 1
            c1 = get_centroid( mesh%elems(elem_id1), mesh )
            normal = get_normal_vector( elem_id1, i, mesh )
            dist = distance_between_points( midpoint, c1 )

            mesh%dual_faces(face_count)%edge_id = i
            mesh%dual_faces(face_count)%node_ids = mesh%edges(i)%node_ids
            mesh%dual_faces(face_count)%start = midpoint
            mesh%dual_faces(face_count)%end = c1
            mesh%dual_faces(face_count)%normal = normal
            mesh%dual_faces(face_count)%area = dist

            ! Face 2: midpoint to centroid 2
            face_count = face_count + 1
            c2 = get_centroid( mesh%elems(elem_id2), mesh )
            normal = get_normal_vector( elem_id2, i, mesh )
            dist = distance_between_points( midpoint, c2 )

            mesh%dual_faces(face_count)%edge_id = i
            mesh%dual_faces(face_count)%node_ids = mesh%edges(i)%node_ids
            mesh%dual_faces(face_count)%start = midpoint
            mesh%dual_faces(face_count)%end = c2
            mesh%dual_faces(face_count)%normal = normal
            mesh%dual_faces(face_count)%area = dist
        else if ( count(mesh%edges(i)%elem_ids >=0) == 1 ) then
            ! Boundary edge: one dual face
            face_count = face_count + 1

            elem_id1 = mesh%edges(i)%elem_ids(1)
            c1 = get_centroid( mesh%elems(elem_id1), mesh )
            normal = get_normal_vector( elem_id1, i, mesh )
            dist = distance_between_points( midpoint, c1 )

            mesh%dual_faces(face_count)%edge_id = i
            mesh%dual_faces(face_count)%node_ids = mesh%edges(i)%node_ids
            mesh%dual_faces(face_count)%start = midpoint
            mesh%dual_faces(face_count)%end = c1
            mesh%dual_faces(face_count)%normal = normal
            mesh%dual_faces(face_count)%area = dist
        else
            print *, 'Error: Edge ', i, 'has no associated elements.'
            stop
        end if
    end do
end subroutine create_dual_mesh

end module dual_mesh
