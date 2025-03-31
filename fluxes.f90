module fluxes
      use kind_defs,            only : dp
      use flow_types,           only : primitive_state, conservative_state,    &
                                       flux_state
      use grid_types,           only : dual_face
      use namelist_definitions, only : gamma
      use initialization,       only : nturb

      implicit none

      private
      public :: primitive_to_conservative, conservative_to_primitive,          &
                compute_flux, lax_friedrichs_flux, roe_flux, hlle_flux

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

end module fluxes
