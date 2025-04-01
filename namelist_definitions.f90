module namelist_definitions
    use kind_defs, only : dp

    implicit none
    
    ! Governing Equations
    character(len=15) :: eqn_type = 'compressible', &
                         viscous_terms = 'turbulent'

    real(dp) :: prandtl_number = 0.72

    ! Reference Physical Properties
    real(dp) :: reynolds_number = 0.0_dp, &
                temperature = 273.0_dp, &
                angle_of_attack = 0.0_dp, &
                mach_number = 0.0_dp, &
                gamma = 1.4_dp, &
                sutherland_constant = -1.0_dp

    character(len=10) :: temperature_units = 'Kelvin'

    ! Code Run Control
    integer :: steps = 500

    real(dp) :: stopping_tolerance = 1.0e-15_dp

    ! Nonlinear Solver Parameters
    character(len=15) :: time_accuracy = 'steady'

    integer :: schedule_iteration(2) = [1,50]

    real(dp) :: schedule_cfl(2) = [200.0_dp, 200.0_dp], &
                schedule_cflturb(2) = [50.0_dp, 50.0_dp]

    ! Inviscid Flux Method
    character(len=15) :: flux_construction = 'roe', &
                         flux_construction_lhs = 'consistent', &
                         flux_limiter = 'none'

    real(dp) :: kappa_umuscl = 0.5

    ! Turbulent Diffusion Models
    character(len=15) :: turbulence_model = 'sa'

    real(dp) :: turb_intensity = -0.001, &
                turb_viscosity_ratio = -0.001, &
                prandtlnumber_turbulent = 0.9

    ! Boundary Conditions
    real(dp), allocatable :: total_pressure_ratio(:), &
                             total_temperature_ratio(:), &
                             q_set(:,:), &
                             static_pressure_ratio(:), &
                             mach_bc(:), &
                             wall_temperature(:)

    integer, allocatable :: boundary_type(:)

    character(len=20), allocatable :: boundary_name(:)

    logical, allocatable :: wall_temp_flag(:)

    
    ! Declare namelists
    namelist /boundary_conditions/ &
        boundary_name, &
        boundary_type, &
        total_pressure_ratio, &
        total_temperature_ratio, &
        q_set, &
        static_pressure_ratio, &
        mach_bc, &
        wall_temp_flag, &
        wall_temperature

    namelist /code_run_control/ &
        steps, &
        stopping_tolerance

    namelist /governing_equations/ &
        eqn_type, &
        viscous_terms, &
        prandtl_number

    namelist /inviscid_flux_method/ &
        flux_construction, &
        flux_construction_lhs, &
        flux_limiter, &
        kappa_umuscl

    namelist /nonlinear_solver_parameters/ &
        time_accuracy, &
        schedule_iteration, &
        schedule_cfl, &
        schedule_cflturb

    namelist /reference_physical_properties/ &
        reynolds_number, &
        temperature, &
        temperature_units, &
        angle_of_attack, &
        mach_number, &
        gamma, &
        sutherland_constant

    namelist /turbulent_diffusion_models/ &
        turbulence_model, &
        turb_intensity, &
        turb_viscosity_ratio, &
        prandtlnumber_turbulent


contains

end module namelist_definitions
