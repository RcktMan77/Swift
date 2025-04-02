module gradients
      use kind_defs,      only : dp
      use grid_types,     only : grid
      use flow_types,     only : primitive_state
      use initialization, only : nturb

      implicit none

      private
      public :: compute_least_squares_gradients

contains


! The compute_least_squares_gradients procedure updates the grad array
! attribute in primitive_state for each node to store gradients used for
! unstructured Monotonic Upstream-centered Scheme for Conservation Laws
! (UMUSCL) reconstruction in inviscid fluxes. The Unweighted Least-Squares
! gradient scheme gathers neighboring nodes for each node, i, and forms
! the following linear system:

! ∇φ = [∂φ/∂x,∂φ/∂y]

! This system is Ax = b, where A is an n x 2 matrix with rows
! [xj - xi, yj - yi] for n neighbors and b is a vector φj - φi. Solve
! AT * Ax = AT * b (2 x 2) system using explicit inversion for efficiency,
! as it is small and unweighted (no distance-based weights).

subroutine compute_least_squares_gradients( mesh, w )
    implicit none

    type(grid), intent(in) :: mesh

    type(primitive_state), intent(inout) :: w(:)

    real(dp), allocatable :: neighbors(:)

    real(dp) :: dx, dy, ATA(2,2), ATb(2), invATA(2,2), det

    integer :: n1, nvars, i, j, k

    nvars = 4 + nturb

    !$omp parallel do                                                          &
    !$omp private(i, j, k, n1, dx, dy, ATA, ATb, invATA, det)
    do i = 1, mesh%num_pts
        do j = 1, nvars
            ATA = 0.0_dp
            ATb = 0.0_dp
            do k = 1, mesh%node_degree(i)
                n1 = mesh%node_neighbors(k, i)

                dx = mesh%points(n1)%x - mesh%points(i)%x
                dy = mesh%points(n1)%y - mesh%points(i)%y

                ATA(1,1) = ATA(1,1) + dx * dx
                ATA(1,2) = ATA(1,2) + dx * dy
                ATA(2,1) = ATA(1,2)
                ATA(2,2) = ATA(2,2) + dy * dy

                if ( j == 1 ) then
                    ATb(1) = ATb(1) + dx * ( w(n1)%rho - w(i)%rho )
                    ATb(2) = ATb(2) + dy * ( w(n1)%rho - w(i)%rho )
                else if ( j <= 3 ) then
                    ATb(1) = ATb(1) + dx * ( w(n1)%vel(j-1) - w(i)%vel(j-1) )
                    ATb(2) = ATb(2) + dy * ( w(n1)%vel(j-1) - w(i)%vel(j-1) )
                else if ( j == 4 ) then
                    ATb(1) = ATb(1) + dx * ( w(n1)%p - w(i)%p )
                    ATb(2) = ATb(2) + dy * ( w(n1)%p - w(i)%p )
                else
                    ATb(1) = ATb(1) + dx * ( w(n1)%turb(j-4) - w(i)%turb(j-4) )
                    ATb(2) = ATb(2) + dy * ( w(n1)%turb(j-4) - w(i)%turb(j-4) )
                end if
            end do

            det = ATA(1,1) * ATA(2,2) - ATA(1,2) * ATA(2,1)

            if ( det  > 1.0e-10 ) then
                invATA(1,1) = ATA(2,2) / det
                invATA(1,2) = -ATA(1,2) / det
                invATA(2,1) = -ATA(2,1) / det
                invATA(2,2) = ATA(1,1) / det

                w(i)%grad(j,:) = matmul( invATA, ATb )
            else
                w(i)%grad(j,:) = 0.0_dp
            end if
        end do
    end do
    !$omp end parallel do
end subroutine compute_least_squares_gradients

end module gradients
