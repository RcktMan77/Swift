module flow_types
    use kind_defs, only : dp

    implicit none

    type :: primitive_state
        real(dp) :: rho
        real(dp) :: vel(2)  ! [u, v]
        real(dp) :: p
        real(dp), allocatable :: turb(:)
        real(dp), allocatable :: residual(:)
        real(dp), allocatable :: grad(:,:)
    end type primitive_state

    type :: conservative_state
        real(dp) :: rho
        real(dp) :: rho_u
        real(dp) :: rho_v
        real(dp) :: rho_E
        real(dp) :: rho_nu
        real(dp), allocatable :: rho_turb(:)
    end type conservative_state

    type :: flux_state
        real(dp) :: rho_u
        real(dp) :: m_flux(2)   ! [ρ × u² + p, ρ × u × v]
        real(dp) :: e_flux
        real(dp), allocatable :: turb_flux(:)  ! ρ × νₜ
    end type flux_state
end module flow_types
