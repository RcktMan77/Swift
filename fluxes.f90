module fluxes
      use kind_defs,            only : dp
      use flow_types,           only : primitive_state, conservative_state,    &
                                       flux_state
      use grid_types,           only : grid, dual_face
      use namelist_definitions, only : gamma, flux_construction, flux_limiter, &
                                       kappa_umuscl, prandtl_number
      use initialization,       only : nturb

      implicit none

      private
      public :: primitive_to_conservative, conservative_to_primitive,          &
                compute_flux, lax_friedrichs_flux, roe_flux, hlle_flux,        &
                compute_limiter, compute_inviscid_fluxes

contains


function primitive_to_conservative( w ) result( q )
    implicit none

    type(primitive_state), intent(in) :: w

    type(conservative_state) :: q

    real(dp) :: rho, u, v, p, E

    integer :: i

    rho    = w%rho
    u      = w%vel(1)
    v      = w%vel(2)
    p      = w%p

    E = p / ( rho * ( gamma - 1.0_dp ) ) + 0.5 * ( u**2 + v**2 )

    q%rho       = rho
    q%rho_u     = rho * u
    q%rho_v     = rho * v
    q%rho_E     = rho * E

    if ( nturb > 0 ) then
        allocate( q%rho_turb(nturb) )
        do i = 1, nturb
            q%rho_turb(i) = rho * w%turb(i)
        end do
    end if
end function primitive_to_conservative


function conservative_to_primitive( q ) result( w )
    implicit none

    type(conservative_state), intent(in) :: q

    type(primitive_state) :: w

    real(dp) :: rho, u, v, E, p

    integer :: i

    rho   = q%rho
    u     = q%rho_u / rho
    v     = q%rho_v / rho
    E     = q%rho_E / rho

    p = (gamma - 1.0_dp) * rho * (E - 0.5_dp * (u**2 + v**2))

    w%rho    = rho
    w%vel(1) = u
    w%vel(2) = v
    w%p      = p

    if ( nturb > 0 ) then
        allocate( w%turb(nturb) )
        do i = 1, nturb
            w%turb(i) = q%rho_turb(i) / rho
        end do
    end if
end function conservative_to_primitive


! The Finite Volume Method (FVM) discretizes the domain into cells & solves
! the Euler equations by integrating over each cell's volume Ω. F·n
! represents the flux normal to the cell boundary, where n = (nx,ny) is the
! outward unit normal. The boundary integral is approximated by summing
! fluxes across each face.

! At each face, we need a numerical flux to approximate F·n based on the states
! qL (left) and qR (right) of the adjacent cells. The compute_flux function
! calculates the exact physical flux F·n at a face fiven the state q and the
! face normal n. It is used by numerical flux functions to compute left &
! right fluxes, which are then combined to approximate the interface flux.

function compute_flux( q, n ) result( F )
    implicit none

    type(conservative_state), intent(in) :: q

    real(dp), intent(in) :: n(2)

    type(flux_state) :: F

    real(dp) :: rho, u, v, p, E, turb1, turb2, un

    type(primitive_state) :: w

    integer :: i

    w = conservative_to_primitive( q )

    rho = w%rho
    u   = w%vel(1)
    v   = w%vel(2)
    p   = w%p
    E   = q%rho_E / rho

    un = u * n(1) + v * n(2)

    F%rho_u = rho * un
    F%m_flux(1) = rho * u * un + p * n(1)
    F%m_flux(2) = rho * v * un + p * n(2)
    F%e_flux = ( rho * E + p ) * un

    if ( nturb > 0 ) then
        allocate( F%turb_flux(nturb) )
        do i = 1, nturb
            F%turb_flux(i) = rho * w%turb(i) * un
        end do
    end if
end function compute_flux
 

! In FVM, the flux at a cell interface must account for wave propagation &
! ensure stability. Numerical flux functions, F_num(qₗ,qᵣ,n) achieve this
! by blending physical fluxes with dissipation. They must be:

! Consistent: Fnum(q,q,n) = F(q)·n
! Stable: Adding upwinding or dissipation to handle discontinuities like
!         shocks.

! In the Lax Friedrichs method:

! F_LF = 0.5 * ( F(qₗ)·n + F(qᵣ)·n ) - λ/2 * ( qᵣ - qₗ )

! where λ is the maximum wave speed, typically:

! λ = max( |uₙ|ₗ + cₗ, |uₙ|ᵣ + cᵣ ) where c = sqrt( γp/ρ ) is the speed
! of sound. The central term, 0.5 *( Fₗ + Fᵣ ) averages the left & right
! physical fluxes. The dissipation term -λ/2(qᵣ - qₗ) adds numerical viscosity
! proportional to the state jump, scaled by the wave speed.

! The central term alone is unstable for hyperbolic systems like the
! Euler equations. The dissipation term mimics upwinding, stabilizing the
! scheme by smearing discontinuities.

! Pros & Cons
!-------------
! Pros: Simple to implement, robust for any hyperbolic system.
! Cons: Highly dissipative, leading to smeared shocks &
!       contact discontinuities.

function lax_friedrichs_flux( q_left, q_right, n ) result( F_num )
    implicit none

    type(conservative_state), intent(in) :: q_left, q_right

    real(dp), intent(in) :: n(2)

    type(flux_state) :: F_left, F_right, F_num

    real(dp) :: lambda, c_left, c_right, un_left, un_right

    type(primitive_state) :: w_left, w_right

    integer :: i

    F_left = compute_flux( q_left, n )
    F_right = compute_flux( q_right, n )

    w_left = conservative_to_primitive( q_left )
    w_right = conservative_to_primitive( q_right )

    c_left = sqrt( gamma * w_left%p / w_left%rho )
    c_right = sqrt( gamma * w_right%p / w_right%rho )

    un_left = w_left%vel(1) * n(1) + w_left%vel(2) * n(2)
    un_right = w_right%vel(1) * n(1) + w_right%vel(2) * n(2)

    lambda = max( abs(un_left) + c_left, abs(un_right) + c_right )

    F_num%rho_u = 0.5_dp * ( F_left%rho_u + F_right%rho_u ) -                &
                  0.5_dp * lambda * ( q_right%rho - q_left%rho )
    F_num%m_flux(1) = 0.5_dp * ( F_left%m_flux(1) + F_right%m_flux(1) ) -    &
                      0.5_dp * lambda * ( q_right%rho_u - q_left%rho_u )
    F_num%m_flux(2) = 0.5_dp * ( F_left%m_flux(2) + F_right%m_flux(2) ) -    &
                      0.5_dp * lambda * ( q_right%rho_v - q_left%rho_v )
    F_num%e_flux = 0.5_dp * ( F_left%e_flux - F_right%e_flux ) -             &
                   0.5_dp * lambda * ( q_right%rho_E - q_left%rho_E )

    if ( nturb > 0 ) then
        allocate( F_num%turb_flux(nturb) )
        do i = 1, nturb
            F_num%turb_flux(i) = 0.5_dp * ( F_left%turb_flux(i) -            &
                F_right%turb_flux(i) ) - 0.5_dp * lambda *                   &
                ( q_right%rho_turb(i) - q_left%rho_turb(i) )
        end do
    end if
end function lax_friedrichs_flux


! The Roe flux is an approximate Riemann solver that linearizes the Euler
! equations around a Roe-averaged state. The numerical flux is given by:

! F_roe = 0.5 * ( F(qₗ)·n + F(qᵣ)·n ) - 0.5 * ∑|λ|αr for k = (1,2,...,m)

! where:
! λ corresponds to the eigenvalues of the Roe-averaged Jacobian,
! r represents the corresponding eigenvectors,
! α are the wave strengths computed from qᵣ - qₗ, &
! k = 4 for the 2-D Euler equations.

! Combines the central flux with wave-specific dissipation. Upwinding is
! tailored to the physics (e.g., acoustic waves, contact waves), reducing
! dissipation compared to Lax-Friedrichs.

! Numerics:
! Roe-Averaging: Compute an average state q_bar (e.g., using Roe averages
! for ρ, u, v, p) such that:

! A(q_bar) * ( qᵣ - qₗ ) = F(qᵣ) - F(qₗ)

! where A = ∂(F·n)/∂q, is the Roe-averaged flux Jacobian, which is linearized
! around the Roe state.

! Eigen-decomposition: Diagonalize A to get λ (wave speeds, e.g., uₙ - c, 
! uₙ + c ) and r.

! Dissipation: Apply dissipation to each wave, proportional to its speed.

! Pros & Cons
!-------------
! Pros: Captures shocks & discontinuities sharply (less dissipative)
! Cons: Complex to implement, can fail for rare cases (e.g., entropy
!       violations), requiring fixes like Harten's entropy fix.

function roe_flux( q_left, q_right, n ) result( F_num )
    implicit none

    type(conservative_state), intent(in) :: q_left, q_right

    real(dp), intent(in) :: n(2)

    type(flux_state) :: F_num, F_left, F_right

    real(dp) :: rho_l, rho_r, u_l, u_r, v_l, v_r, p_l, p_r,                    &
                rho_roe, u_roe, v_roe, H_roe, c_roe, un_roe,                   &
                delta_rho, delta_p, delta_un, delta_ut

    real(dp), allocatable :: t_l(:), t_r(:), t_roe(:)

    real(dp), allocatable :: delta_U(:), lambda(:), abs_lambda(:), alpha(:),   &
                             tmp(:)

    real(dp), allocatable :: R(:,:), K(:,:)

    type(primitive_state) :: w_left, w_right

    integer :: i, j, nvars

    ! Compute physical fluxes for left & right states
    F_left = compute_flux( q_left, n )
    F_right = compute_flux( q_right, n )

    ! Convert conservative states to primitive states
    w_left = conservative_to_primitive( q_left )
    w_right = conservative_to_primitive( q_right )

    ! Extract primitive variables
    rho_l = w_left%rho
    u_l   = w_left%vel(1)
    v_l   = w_left%vel(2)
    p_l   = w_left%p

    rho_r = w_right%rho
    u_r   = w_right%vel(1)
    v_r   = w_right%vel(2)
    p_r   = w_right%p

    if ( nturb > 0 ) then
        allocate( t_l(nturb), t_r(nturb), t_roe(nturb) )
        t_l = w_left%turb
        t_r = w_right%turb
    end if

    ! Compute Rho-averaged state
    rho_roe = sqrt( rho_l * rho_r )
    u_roe = ( u_l * sqrt(rho_l) + u_r * sqrt(rho_r) ) /                        &
            ( sqrt(rho_l) + sqrt(rho_r) )
    v_roe = ( v_l * sqrt(rho_l) + v_r * sqrt(rho_r) ) /                        &
            ( sqrt(rho_l) + sqrt(rho_r) )
    H_roe = ( (p_l / (gamma - 1.0_dp) + 0.5_dp * rho_l * (u_l**2 + v_l**2)) /  &
              sqrt(rho_l) +                                                    &
              (p_r / (gamma - 1.0_dp) + 0.5_dp * rho_r * (u_r**2 + v_r**2)) /  &
              sqrt(rho_r) ) / ( sqrt(rho_l) + sqrt(rho_r) )

    if ( nturb > 0 ) then
        do i = 1, nturb
            t_roe(i) = ( t_l(i) * sqrt(rho_l) + t_r(i) * sqrt(rho_r) ) /       &
                       ( sqrt(rho_l) + sqrt(rho_r) )
        end do
    end if

    c_roe = sqrt( (gamma - 1.0_dp) * (H_roe - 0.5_dp * (u_roe**2 + v_roe**2)) )
    un_roe = u_roe * n(1) + v_roe * n(2)

    nvars = 4 + nturb
    allocate( delta_U(nvars), lambda(nvars), R(nvars,nvars),                   &
              abs_lambda(nvars), alpha(nvars), tmp(nvars), K(nvars,nvars) )

    ! delta_U is the difference between left & right conservative states. This
    ! difference is the input to the dissipation term, telling us how much
    ! the state changes across the interface.
    delta_U(1) = q_right%rho - q_left%rho
    delta_U(2) = q_right%rho_u - q_left%rho_u
    delta_U(3) = q_right%rho_v - q_left%rho_v
    delta_U(4) = q_right%rho_E - q_left%rho_E

    if ( nturb > 0 ) then
        do i = 1, nturb
            delta_U(4 + i) = q_right%rho_turb(i) - q_left%rho_turb(i)
        end do
    end if

    ! Define eigenvalues (wave speeds). These determine how fast information
    ! travels across the interface. The absolute values |λ| control the
    ! amount of dissipation applied to each wave.
    lambda(1) = un_roe - c_roe      ! Left-moving acoustic wave
    lambda(2) = un_roe              ! Entropy wave
    lambda(3) = un_roe              ! Shear wave
    lambda(4) = un_roe + c_roe      ! Right-moving acoustic wave

    if ( nturb > 0 ) then
        do i = 1, nturb
            lambda(4 + i) = un_roe  ! Turbulence variables (advection speed)
        end do
    end if

    abs_lambda = abs( lambda )

    ! Construct the right eigenvector matrix (R). This eigenvector matrix
    ! allows us to decompose the jump Δq into contributions from each
    ! wave type (acoustic, entropy, shear).
    R = 0.0_dp
    R(:,1) = [1.0_dp, u_roe - c_roe * n(1), v_roe - c_roe * n(2),              &
              H_roe - un_roe * c_roe]
    R(:,2) = [1.0_dp, u_roe, v_roe, 0.5_dp * (u_roe**2 + v_roe**2)]
    R(:,3) = [0.0_dp, -n(2), n(1), -(u_roe * n(2) - v_roe * n(1))]
    R(:,4) = [1.0_dp, u_roe + c_roe * n(1), v_roe + c_roe * n(2),              &
              H_roe + un_roe * c_roe]

    if ( nturb > 0 ) then
        do i = 1, nturb
            R(1:4, 4 + i) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
            R(4 + i, 4 + i) = 1.0_dp
            R(4 + i, 1:4) = [t_roe(i), t_roe(i), 0.0_dp, t_roe(i)]
        end do
    end if

    ! Compute wave strengths (α) for Euler variables
    delta_rho = rho_r - rho_l
    delta_p = p_r - p_l
    delta_un = ( u_r - u_l ) * n(1) + ( v_r - v_l ) * n(2)
    delta_ut = -( u_r - u_l ) * n(2) + ( v_r - v_l ) * n(1)

    alpha(1) = ( delta_p - rho_roe * c_roe * delta_un ) / ( 2.0_dp * c_roe**2 )
    alpha(2) = delta_rho  - delta_p / ( c_roe**2 )
    alpha(3) = rho_roe * delta_ut
    alpha(4) = ( delta_p + rho_roe * c_roe * delta_un ) / (2.0_dp * c_roe**2 )

    if ( nturb > 0 ) then
        do i = 1, nturb
            alpha(4 + i) = delta_U(4 + i) / rho_roe
        end do
    end if

    ! Compute the dissipation term (|A| * delta_U) wave-by-wave
    tmp = 0.0_dp
    do i = 1, nvars
        tmp = tmp + abs_lambda(i) * alpha(i) * R(:,i)
    end do

    ! Compute the final numerical flux
    F_num%rho_u = 0.5_dp * ( F_left%rho_u + F_right%rho_u ) - 0.5_dp * tmp(1)
    F_num%m_flux(1) = 0.5_dp * ( F_left%m_flux(1) + F_right%m_flux(1) ) -      &
                      0.5_dp * tmp(2)
    F_num%m_flux(2) = 0.5_dp * ( F_left%m_flux(2) + F_right%m_flux(2) ) -      &
                      0.5_dp * tmp(3)
    F_num%e_flux = 0.5_dp * ( F_left%e_flux + F_right%e_flux ) - 0.5_dp * tmp(4)

    if ( nturb > 0 ) then
        do i = 1, nturb
            F_num%turb_flux(i) = 0.5_dp * ( F_left%turb_flux(i) +              &
                F_right%turb_flux(i) ) - 0.5_dp * tmp(4 + i)
        end do
    end if
end function roe_flux


! The hlle_flux function is a relatively simple, robust Riemann solver that
! approximates the flux at a cell interface by considering only the fastest
! left- & right-moving waves, ignoring intermediate waves (e.g., contact
! discontinuities). It's less dissipative then Lax Friedrichs but coarser
! than Roe, making it a good middle-ground choise. The formulation ensures
! positivity (e.g., ρ > 0) and stability, but can overly smear contact
! discontinuities.

! This function implements the HLLE++ scheme that includes enhanced features
! that build on HLLE by using more sophisticated bounds on wave speed estimates
! to reduce dissipation, and adds corrections to prevent unphysical solutions
! (e.g. expansion shocks) similar to Harten's entropy fix in Roe's solvers.
function hlle_flux( q_left, q_right, n ) result( F_num )
    implicit none
    
    type(conservative_state), intent(in) :: q_left, q_right

    real(dp), intent(in) :: n(2)

    type(flux_state) :: F_num, F_left, F_right

    real(dp) :: rho_l, rho_r, u_l, u_r, v_l, v_r, p_l, p_r,                    &
                rho_roe, u_roe, v_roe, H_roe, c_roe, un_roe,                   &
                lambda_min, lambda_max, lambda_min_l, lambda_max_r,            &
                c_l, c_r, un_l, un_r

    type(primitive_state) :: w_left, w_right

    integer :: i

    ! Compute physical fluxes
    F_left = compute_flux( q_left, n )
    F_right = compute_flux( q_right, n )

    ! Convert to primitive states
    w_left = conservative_to_primitive( q_left )
    w_right = conservative_to_primitive( q_right )

    ! Extract primitive variables
    rho_l = w_left%rho
    u_l   = w_left%vel(1)
    v_l   = w_left%vel(2)
    p_l   = w_left%p

    rho_r = w_right%rho
    u_r   = w_right%vel(1)
    v_r   = w_right%vel(2)
    p_r   = w_right%p

    ! Compute local wave speeds
    un_l = u_l * n(1) + v_l * n(2)
    un_r = u_r * n(1) + v_r * n(2)
    c_l = sqrt( gamma * p_l / rho_l )
    c_r = sqrt( gamma * p_r / rho_r )

    ! Compute Roe-averaged state for wave speed estimates
    rho_roe = sqrt( rho_l * rho_r )
    u_roe = ( u_l * sqrt(rho_l) + u_r * sqrt(rho_r) ) /                        &
            ( sqrt(rho_l) + sqrt(rho_r) )
    v_roe = ( v_l * sqrt(rho_l) + v_r * sqrt(rho_r) ) /                        &
            ( sqrt(rho_l) + sqrt(rho_r) )
    H_roe = ( (p_l / (gamma - 1.0_dp) + 0.5_dp * rho_l * (u_l**2 + v_l**2)) /  &
              sqrt(rho_l) +                                                    &
              (p_r / (gamma - 1.0_dp) + 0.5_dp * rho_r * (u_r**2 + v_r**2)) /  &
              sqrt(rho_r) ) / ( sqrt(rho_l) + sqrt(rho_r) )
    c_roe = sqrt( (gamma - 1.0_dp) * (H_roe - 0.5_dp * (u_roe**2 + v_roe**2)) )
    un_roe = u_roe * n(1) + v_roe * n(2)

    ! Wave speed estimates (improved over standard HLLE)
    lambda_min_l = un_l - c_l
    lambda_max_r = un_r + c_r
    lambda_min = min( lambda_min_l, un_roe - c_roe )
    lambda_max = max( lambda_max_r, un_roe + c_roe )

    ! Entropy fix: Widen wave speeds near sonic points
    if ( lambda_min_l < 0.0_dp .and. un_roe - c_roe > 0.0_dp ) then
        lambda_min = lambda_min_l - 0.2_dp * ( un_roe - c_roe - lambda_min_l )
    end if

    if ( lambda_max_r > 0.0_dp .and. un_roe + c_roe < 0.0_dp ) then
        lambda_max = lambda_max_r + 0.2_dp * ( un_roe + c_roe - lambda_max_r )
    end if

    ! Compute HLLE++ Flux
    if ( lambda_min >= 0.0_dp ) then
        F_num = F_left
    else if ( lambda_max <= 0.0_dp ) then
        F_num = F_right
    else
        F_num%rho_u = ( lambda_max * F_left%rho_u - lambda_min *               &
                        F_right%rho_u + lambda_max * lambda_min *              &
                        (q_right%rho - q_left%rho) ) /                         &
                        ( lambda_max - lambda_min )

        F_num%m_flux(1) = ( lambda_max * F_left%m_flux(1) - lambda_min *       &
                            F_right%m_flux(1) + lambda_max * lambda_min *      &
                            (q_right%rho_u - q_left%rho_u) ) /                 &
                            ( lambda_max - lambda_min )

        F_num%m_flux(2) = ( lambda_max * F_left%m_flux(2) - lambda_min *       &
                            F_right%m_flux(2) + lambda_max * lambda_min *      &
                            (q_right%rho_v - q_left%rho_v) ) /                 &
                            ( lambda_max - lambda_min )

        F_num%e_flux = ( lambda_max * F_left%e_flux - lambda_min *             &
                         F_right%e_flux + lambda_max * lambda_min *            &
                         (q_right%rho_E - q_left%rho_E) ) /                    &
                         ( lambda_max - lambda_min )

        if ( nturb > 0 ) then
            allocate( F_num%turb_flux(nturb) )

            do i = 1, nturb
                F_num%turb_flux(i) = ( lambda_max * F_left%turb_flux(i) -      &
                    lambda_min * F_right%turb_flux(i) + lambda_max *           &
                    lambda_min * (q_right%rho_turb(i) - q_left%rho_turb(i)) ) /&
                    ( lambda_max - lambda_min )
            end do
        end if
    end if
end function hlle_flux


! This function computes the limiter value ψ, based on the ratio, r and
! the limiter type.
function compute_limiter( r, lim_type ) result( psi )
    implicit none

    real(dp), intent(in) :: r

    character(len=*), intent(in) :: lim_type

    real(dp) :: psi

    select case( trim(lim_type) )
    case ( 'vanalbada' )
        psi = ( r**2 + r ) / ( r**2 + 1.0_dp )
    case ( 'vanleer' )
        psi = ( r + abs(r) ) / ( 1.0_dp + abs(r) )
    case ( 'barth' )
        ! Barth requires pre-computation; return 1.0 here, adjust in subroutine
        psi = 1.0_dp
    case ( 'none' )
        psi = 1.0_dp
    case default
        psi = 1.0_dp ! Fallback to no limiting
    end select
end function compute_limiter


! The compute_inviscid_fluxes subroutine reconstructs left & right states
! at each dual face, converts to conservative variables, and computes fluxes
! using an existing flux scheme, configurable via flux_construction

! Numerics: For each dual face (stored in mesh%dual_faces):

! Use node gradients (w(i)%grad) to exrapolate primitive variables (ρ, u, v,
! p, turb) to the face midpoint from both nodes.

! Extrapolation: φₗ = φn1 + ∇φn1 · rn1->f, φᵣ = φn2 + ∇φn2 · rn2->f,
! where r uses dual_face%area and direction from node to midpoint.

! Applies user-selected limiter, ψ, based on r where:

! r = (φn2 - φn1) / (∇φn1 · r),

! to limit the gradient contribution, ensuring monotonicity without excessive
! dissipation.

subroutine compute_inviscid_fluxes( mesh, w )
    implicit none

    type(grid), intent(inout) :: mesh

    type(primitive_state), intent(in) :: w(:)

    type(conservative_state) :: q_left, q_right

    type(primitive_state) :: w_left, w_right

    type(flux_state) :: F
    
    real(dp) :: r_vec(2), r, psi, dx, dy, grad_contrib

    integer :: n1, n2, nvars, i, j

    nvars = 4 + nturb

    do i = 1, mesh%num_edges
        if ( .not. allocated( mesh%edges(i)%flux ) ) then
            allocate( mesh%edges(i)%flux(nvars) )
        end if

        n1 = mesh%edges(i)%node_ids(1)
        n2 = mesh%edges(i)%node_ids(2)

        dx = ( mesh%points(n2)%x - mesh%points(n1)%x ) * 0.5_dp
        dy = ( mesh%points(n2)%y - mesh%points(n1)%y ) * 0.5_dp

        r_vec = [dx, dy]

        ! Initialize left & right states
        w_left = w(n1)
        w_right = w(n2)

        if ( nturb > 0 ) then
            allocate( w_left%turb(nturb), w_right%turb(nturb) )
            w_left%turb = w(n1)%turb
            w_right%turb = w(n2)%turb
        end if

        ! Reconstruct with selected limiter
        do j = 1, nvars
            if ( j == 1 ) then
                ! Left state
                grad_contrib = dot_product( w(n1)%grad(j,:), r_vec )
                r = ( w(n2)%rho - w(n1)%rho ) / ( grad_contrib + 1.0e-10_dp )
                psi = compute_limiter( r, flux_limiter )
                w_left%rho = w(n1)%rho + kappa_umuscl * psi * grad_contrib

                ! Right state
                grad_contrib = dot_product( w(n2)%grad(j,:), -r_vec )
                r = ( w(n1)%rho - w(n2)%rho ) / ( grad_contrib + 1.0e-10_dp )
                psi = compute_limiter( r, flux_limiter )
                w_right%rho = w(n2)%rho + kappa_umuscl * psi * grad_contrib
            else if ( j <= 3 ) then
                grad_contrib = dot_product( w(n1)%grad(j,:), r_vec )
                r = ( w(n2)%vel(j-1) - w(n1)%vel(j-1) ) /                      &
                    ( grad_contrib + 1.0e-10_dp )
                psi = compute_limiter( r, flux_limiter )
                w_left%vel(j-1) = w(n1)%vel(j-1) + kappa_umuscl * psi *        &
                    grad_contrib

                grad_contrib = dot_product( w(n2)%grad(j,:), -r_vec )
                r = ( w(n1)%vel(j-1) - w(n2)%vel(j-1) ) /                      &
                    ( grad_contrib + 1.0e-10_dp )
                psi = compute_limiter( r, flux_limiter )
                w_right%vel(j-1) = w(n2)%vel(j-1) + kappa_umuscl * psi *       &
                    grad_contrib
            else if ( j == 4 ) then
                grad_contrib = dot_product( w(n1)%grad(j,:), r_vec )
                r = ( w(n2)%p - w(n1)%p ) / ( grad_contrib + 1.0e-10_dp )
                psi = compute_limiter( r, flux_limiter )
                w_left%p = w(n1)%p + kappa_umuscl * psi * grad_contrib

                grad_contrib = dot_product( w(n2)%grad(j,:), -r_vec )
                r = ( w(n1)%p - w(n2)%p ) / ( grad_contrib + 1.0e-10_dp )
                psi = compute_limiter( r, flux_limiter )
                w_right%p = w(n2)%p + kappa_umuscl * psi * grad_contrib
            else
                grad_contrib = dot_product( w(n1)%grad(j,:), r_vec )
                r = ( w(n2)%turb(j-4) - w(n1)%turb(j-4) ) /                    &
                    ( grad_contrib + 1.0e-10 )
                psi = compute_limiter( r, flux_limiter )
                w_left%turb(j-4) = w(n1)%turb(j-4) + kappa_umuscl * psi *      &
                    grad_contrib

                grad_contrib = dot_product( w(n2)%grad(j,:), -r_vec )
                r = ( w(n1)%turb(j-4) - w(n2)%turb(j-4) ) /                    &
                    ( grad_contrib + 1.0e-10_dp )
                psi = compute_limiter( r, flux_limiter )
                w_right%turb(j-4) = w(n2)%turb(j-4) + kappa_umuscl * psi *     &
                    grad_contrib
            end if
        end do

        ! Compute flux
        q_left = primitive_to_conservative( w_left )
        q_right = primitive_to_conservative( w_right )

        select case ( trim(flux_construction) )
        case ( 'lax-friedrichs' )
            F = lax_friedrichs_flux( q_left, q_right,                          &
                                     mesh%dual_faces(i)%normal )
        case ( 'roe' )
            F = roe_flux( q_left, q_right, mesh%dual_faces(i)%normal )
        case ( 'hlle' )
            F = hlle_flux( q_left, q_right, mesh%dual_faces(i)%normal )
        case default
            F = roe_flux( q_left, q_right, mesh%dual_faces(i)%normal )
        end select

        mesh%edges(i)%flux(1) = F%rho_u
        mesh%edges(i)%flux(2:3) = F%m_flux
        mesh%edges(i)%flux(4) = F%e_flux

        if ( nturb > 0 ) then
            mesh%edges(i)%flux(5:4+nturb) = F%turb_flux
        end if

        if ( allocated( w_left%turb ) ) then
            deallocate( w_left%turb, w_right%turb )
        end if
    end do
end subroutine compute_inviscid_fluxes

end module fluxes
