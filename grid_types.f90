module grid_types
    use kind_defs, only : dp

    implicit none

    type :: point
        real(dp) :: x, y
        integer :: node_id
    end type point

    type :: edge
        real(dp) :: length
        integer :: bndry_id = 0
        logical :: shares_bndry_node  ! Used to identify interior edges that
                                      ! share a node with boundary elements
        integer :: node_ids(2)
        integer, allocatable :: elem_ids(:)
        real(dp), allocatable :: flux(:,:)
    end type edge

    type :: dual_face
        integer :: edge_id            ! Primal edge index
        integer :: node_ids(2)        ! Node indices from primal edge
        type(point) :: start          ! Starting point (e.g., edge midpoint)
        type(point) :: end            ! Ending point (e.g., element centroid)
        real(dp) :: normal(2)         ! Normal vector of the dual face
        real(dp) :: area              ! Length of the dual face (area in 2-D)
    end type dual_face

    type :: element
        integer :: num_nodes
        integer :: num_edges
        integer, allocatable :: node_ids (:)
        integer, allocatable :: edge_ids(:)
        integer, allocatable :: neighbors(:) ! Store neighboring element IDs
        real(dp), allocatable :: edge_normals(:,:) ! Store the oriented normal
        real(dp) :: area                           ! for each edge
        logical :: is_bndry
        integer :: bndry_id ! 0 for interior elements
                            ! booundary index for boundary elements
    end type element

    type :: boundary
        character(len=20) :: name
        integer :: num_elems, bc_type
        real(dp) :: ptsp
        real(dp) :: ttsp
        real(dp) :: q(5)
        real(dp) :: psp
        real(dp) :: mach
        real(dp) :: wall_temperature
        logical :: wall_temp_flag
    end type boundary

    type :: grid
        integer :: grid_dim, num_pts, num_elems, num_tris, num_quads,          &
                   num_interior_edges, num_bndry_edges, num_edges, num_bndrys
        type(point), allocatable :: points(:)
        type(edge), allocatable :: edges(:)
        type(dual_face), allocatable :: dual_faces(:)
        type(element), allocatable :: elems(:)
        type(boundary), allocatable :: bndrys(:)
        real(dp), allocatable :: dx(:)
        integer, allocatable :: node_degree(:)  ! Number of neighbors per node
        integer, allocatable :: node_neighbors(:,:) ! Neighbor IDs
    end type grid

end module grid_types
