module initialization
    use kind_defs,  only: dp
    use grid_types, only: grid
    use flow_types, only: primitive_state
    use namelist_definitions

    implicit none

    private
    public :: check_and_read_namelist, initialize_flow, nturb

    integer :: nturb        ! Number of turbulence transport equations

contains


! The check_and_read_namelist subroutine checks that an input namelist file
! exists, and if found reads user-specified parameters needed for a
! simulation.

subroutine check_and_read_namelist( mesh, filename )
    implicit none

    type(grid), intent(inout) :: mesh

    character(len=*), intent(in) :: filename

    integer :: unit_num, line_num, ios, i, j

    logical :: file_exists, found_ref_phys, found_code_run, found_nonlin,      &
               found_invis, found_bndry, found_gov, found_turb

    character(len=80) :: line


    ! Check if the file exists
    inquire( file=filename, exist=file_exists )

    if ( .not. file_exists ) then
        stop 'Error: Namelist input file does not exist: ' // trim( filename )
    end if

    ! Open namelist file
    open( newunit=unit_num, file=filename, status='old', action='read',        &
          iostat=ios )

    if ( ios /= 0 ) then
        stop 'Error: Failed to open namelist file: ' // trim( filename )
    end if

    ! Allocate boundary condition arrays before reading
    call allocate_boundary_arrays( mesh )

    ! Initialize flags for required namelists
    found_ref_phys = .false.
    found_code_run = .false.
    found_nonlin   = .false.
    found_bndry    = .false.
    found_gov      = .false.     ! Optional namelist
    found_invis    = .false.     ! Optional namelist
    found_turb     = .false.     ! Optional namelist

    line_num = 0

    ! Read the namelist file line-by-line
    do
        read(unit_num, '(A)', iostat=ios) line

        line_num = line_num + 1

        if ( ios /= 0 ) then
            if ( ios > 0 ) write(*,*) 'EOF or error at line ', line_num
            exit
        end if

        ! Skip blank lines or comments
        line = adjustl( trim(line) )
        if ( len_trim(line) == 0 .or. line(1:1) == '!' ) cycle

        ! Detect namelist start
        if ( line(1:1) == '&' ) then
            ! Backspace to reposition the file pointer to the namelist start
            backspace( unit_num)
            
            line_num = line_num - 1
            select case ( trim(line(2:)) )

            case ( 'reference_physical_properties' )
                read(unit_num, nml=reference_physical_properties, iostat=ios)

                if ( ios /= 0 ) then
                    write(*,*) 'Error parsing reference_physical_properties ', &
                        'namelist. ios = ', ios
                    stop
                end if

                found_ref_phys = .true.

            case ( 'code_run_control' )
                read(unit_num, nml=code_run_control, iostat=ios)

                if ( ios /= 0 ) then
                    write(*,*) 'Error parsing code_run_control namelist. ',    &
                        'ios = ', ios
                    stop
                end if

                found_code_run = .true.

            case ( 'nonlinear_solver_parameters' )
                read(unit_num, nml=nonlinear_solver_parameters, iostat=ios)

                if ( ios /= 0 ) then
                    write(*,*) 'Error parsing nonlinear_solver_parameters ',   &
                        'namelist. ios = ', ios
                    stop
                end if

                found_nonlin = .true.

            case ( 'inviscid_flux_method' )
                read(unit_num, nml=inviscid_flux_method, iostat=ios)

                if ( ios /= 0 ) then
                    write(*,*) 'Error parsing inviscid_flux_method ',   &
                        'namelist. ios = ', ios
                    stop
                end if

                found_invis = .true.

            case ( 'governing_equations' )
                read(unit_num, nml=governing_equations, iostat=ios)

                if ( ios /= 0 ) then
                    write(*,*) 'Error parsing governing_equations namelist. ', &
                        'ios = ', ios
                    stop
                end if

                found_gov = .true.

            case ( 'turbulent_diffusion_models' )
                read(unit_num, nml=turbulent_diffusion_models, iostat=ios)

                if ( ios /= 0 ) then
                    write(*,*) 'Error parsing turbulent_diffusion_models ',    &
                        'namelist. ios = ', ios
                    stop
                end if

                found_turb = .true.

            case ( 'boundary_conditions' )
                read(unit_num, nml=boundary_conditions, iostat=ios)

                if ( ios /=0 ) then
                    write(*,*) 'Error parsing boundary_conditions namelist. ', &
                        'ios = ', ios
                    stop
                end if

                found_bndry = .true.

            case default
                write(*,*) 'Warning: Unknown namelist found: ', trim(line)
                
                ! Skip to next namelist be reading until '/'
                do
                    read(unit_num, '(A)', iostat=ios) line
                    line_num = line_num + 1
                    if ( ios /= 0 .or. trim(line) == '/' ) exit
                end do
            end select
        end if
    end do

    ! Close file
    close( unit_num )

    ! Check for required namelists
    if ( .not. found_ref_phys ) then
        stop 'Error: Required namelist reference_physical_properties not found.'
    end if
    if ( .not. found_code_run ) then
        stop 'Error: Required namelist code_run_control not found.'
    end if
    if ( .not. found_nonlin ) then
        stop 'Error: Required namelist nonlinear_solver_parameters not found.'
    end if
    if ( .not. found_bndry ) then
        stop 'Error: Required namelist boundary_conditions not found.'
    end if
    if ( .not. found_gov ) then
        ! Optional namelist
    end if
    if ( .not. found_invis ) then
        ! Optional namelist
    end if
    if ( trim(viscous_terms) == 'turbulent' .and. .not. found_turb ) then
        stop 'Error: turbulent_diffusion_models namelist required when ' //    &
             'viscous_terms = "turbulent".'
    end if

    ! Process boundary conditions
    if ( .not. allocated( boundary_name ) .or. size( boundary_name ) /=        &
         mesh%num_bndrys ) then
        write(*,*) 'Error: Boundary conditions in namelist do not match ',     &
             'grid boundaries.'
        stop
    end if

    do i = 1, mesh%num_bndrys
        do j = 1, size( boundary_name )
            if ( trim( boundary_name(j) ) == trim( mesh%bndrys(i)%name ) ) then
                mesh%bndrys(i)%bc_type = boundary_type(j)

                select case ( boundary_type(j) )
                    case ( 3000 ) ! Tangency (zero normal velocity via fluxes)
                        ! No additional parameters
                    case ( 4000, 4110 ) ! Viscous (explicit/implicit no-slip)
                        mesh%bndrys(i)%wall_temp_flag = wall_temp_flag(j)
                        if ( wall_temp_flag(j) ) then
                            mesh%bndrys(i)%wall_temperature =                  &
                                wall_temperature(j)
                        else
                            ! Adiabatic default
                            mesh%bndrys(i)%wall_temperature = -1.0_dp
                        end if
                    case ( 5000 ) ! Farfield (Riemann invariants)
                        ! No additional parameters
                    case ( 5050 ) ! Freestream (external freestream via fluxes)
                        ! No additonal parameters
                    case ( 5051 ) ! Subsonic outflow (back pressure)
                        mesh%bndrys(i)%psp = static_pressure_ratio(j)
                    case ( 5052 ) ! Subsonic outflow (Mach number)
                        mesh%bndrys(i)%mach = mach_bc(j)
                    case ( 7011 ) ! Subsonic inflow
                        mesh%bndrys(i)%ptsp = total_pressure_ratio(j)
                        mesh%bndrys(i)%ttsp = total_temperature_ratio(j)
                    case ( 7012 ) ! Subsonic outflow
                        mesh%bndrys(i)%psp = static_pressure_ratio(j)
                    case default
                        write(*,*) 'Error: Unknown boundary_type ',            &
                            boundary_type(j), ' for ', trim( boundary_name(j) )
                        stop
                end select
                    
                exit
            end if

            if ( j == size( boundary_name ) ) then
                write(*,*) 'Error: Boundary name not found in grid: ',         &
                    trim( mesh%bndrys(i)%name )
                stop
            end if
        end do
    end do

    ! Deallocate boundary condition arrays after processing
    call deallocate_boundary_arrays()
end subroutine check_and_read_namelist


! The allocate_boundary_arrays subroutine allocates temporary arrays for
! each paramter of the boundary derived type, which corresponds to required
! input parameters associated with various boundary condition types. The
! given grid derived type object is used to determine the number of
! boundaries that exist in the grid, and this value is then used for
! allocating these temporary arrays. The program associates boundary condition
! information related to a unique boundary with a common array index.

subroutine allocate_boundary_arrays( mesh )
    implicit none

    type(grid), intent(in) :: mesh

    integer :: n

    n = mesh%num_bndrys

    allocate( boundary_name(n) )
    allocate( boundary_type(n) )
    allocate( total_pressure_ratio(n) )
    allocate( total_temperature_ratio(n) )
    allocate( q_set(5,n) ) ! Assuming 5 components for state vector
    allocate( static_pressure_ratio(n) )
    allocate( mach_bc(n) )
    allocate( wall_temperature(n) )
    allocate( wall_temp_flag(n) )

    ! Initialize arrays
    boundary_name = ''
    boundary_type = 0
    total_pressure_ratio = 1.0_dp
    total_temperature_ratio = 1.0_dp
    q_set = 0.0_dp
    static_pressure_ratio = 1.0_dp
    mach_bc = 0.0_dp
    wall_temperature = -1.0_dp ! Default for adiabatic wall
    wall_temp_flag = .false.
end subroutine allocate_boundary_arrays


! The deallocate_boundary_arrays subroutine deallocates the temporary arrays
! that were allocated in the allocate_boundary_arrays routine one the user-
! specified boundary information is mapped to the boundary derived type for
! each boundary.

subroutine deallocate_boundary_arrays()
    implicit none

    if ( allocated( boundary_name ) ) deallocate( boundary_name )
    if ( allocated( boundary_type ) ) deallocate( boundary_type )
    if ( allocated( total_pressure_ratio ) ) deallocate( total_pressure_ratio )
    if ( allocated( q_set ) ) deallocate( q_set )
    if ( allocated( mach_bc ) ) deallocate( mach_bc )
    if ( allocated( wall_temperature ) ) deallocate( wall_temperature )
    if ( allocated( wall_temp_flag ) ) deallocate( wall_temp_flag )

    if ( allocated( static_pressure_ratio ) ) then
        deallocate( static_pressure_ratio )
    end if

    if ( allocated( total_temperature_ratio ) ) then
        deallocate( total_temperature_ratio )
    end if
end subroutine deallocate_boundary_arrays


! The sutherland_viscosity function calculates the laminar viscosity for a
! given temperature in degrees Rankine using Sutherland's Law.

function sutherland_viscosity( tsp ) result( mu )
    implicit none

    real(dp), intent(in) :: tsp     ! Temperature in Rankine
    real(dp) :: mu

    real(dp), parameter :: mu0 = 3.5839e-7_dp ! Reference viscosity (lbf-s/ft²)
    real(dp), parameter :: T0 = 518.67_dp     ! Reference temperature (°R)
    real(dp), parameter :: S = 198.72_dp      ! Sutherland constant (°R)

    if ( sutherland_constant > 0.0_dp ) then
        mu = mu0 * ( tsp / T0 )**1.5_dp * ( T0 + sutherland_constant ) /       &
            ( tsp + sutherland_constant )
    else
        mu = mu0 * ( tsp / T0 )**1.5_dp * ( T0 + S ) / ( tsp + S )
    end if
end function sutherland_viscosity


! The initialize_flow subroutine updates the given primitive state derived
! type object at every node in the given mesh with non-dimensional freestream
! values based on those specified in the reference_physical_properties
! namelist.

subroutine initialize_flow( mesh, w )
    implicit none

    type(grid), intent(inout) :: mesh

    type(primitive_state), allocatable, intent(inout) :: w(:)

    real(dp) :: tsp, u_bar, v_bar, rho_bar, p_bar, turb_bar
    real(dp) :: rho_inf, a_inf, mu_inf, pi, aoa_rad

    integer :: i

    real(dp), parameter :: rgas = 53.3523_dp    ! [ft-lbf/(lbm-°R)]
    real(dp), parameter :: g = 32.174_dp        ! [lbm-ft/(lbf-s²)]

    ! Determine nturb based on viscous_terms and turbulence model
    select case ( trim(viscous_terms) )
    case ( 'laminar' )
        nturb = 0
    case ( 'turbulent' )
        select case ( trim(turbulence_model) )
        case ( 'sa', 'sa-neg' )
            nturb = 1
        case ( 'sst' )
            nturb = 2
        case default
            print *, 'Error: Unknown turbulence_model: ', trim(turbulence_model)
            stop
        end select
    case default
        print *, 'Error: Unknown viscous_terms value: ', trim(viscous_terms)
        stop
    end select

    ! Allocate q array
    allocate( w(mesh%num_pts) )

    pi = acos( -1.0_dp )

    ! Convert temperature to Rankine (consistent with namelist)
    if ( trim( temperature_units )  == 'Kelvin' ) then
        tsp = temperature * 9.0_dp / 5.0_dp     ! K to °R: T(°R) = T(K) * 9/5
    else if ( trim ( temperature_units )  == 'Rankine' ) then
        tsp = temperature
    else
        write(*,*) 'Error: Unknown temperature_units specified.'
        stop
    end if

    ! Convert angle-of-attack from degrees to radians
    aoa_rad = angle_of_attack * pi / 180.0_dp

    ! Dimensional free stream reference values (EEU)
    rho_inf = 2116.216_dp / ( rgas * tsp )  ! ρ∞ = p∞ / (R * T∞) [slug/ft³]
    a_inf = sqrt( gamma * rgas * g *  tsp ) ! Speed of sound [ft/s]
    mu_inf = sutherland_viscosity( tsp )    ! Dimensional viscosity [lbf-s/ft²]

    ! Non-dimensional free stream values
    rho_bar = 1.0                           ! ρ / ρ∞
    u_bar = mach_number * cos( aoa_rad )    ! u / a∞
    v_bar = mach_number * sin( aoa_rad )    ! v / a∞
    p_bar = 1.0_dp / gamma                  ! p∞ / (ρ∞ * a∞²) = 1/γ
    turb_bar = mu_inf / ( rho_inf * a_inf ) ! μ∞ / (ρ∞ * a∞ * L), L = 1 ft

    ! Initialize flow variables at each node
    do i = 1, mesh%num_pts
        w(i)%rho = rho_bar
        w(i)%vel(1) = u_bar                 ! Non-dimensional u
        w(i)%vel(2) = v_bar                 ! Non-dimensional v
        w(i)%p = p_bar                      ! Non-dimensional ν

        if ( nturb > 0 ) then
            allocate( w(i)%turb(nturb) )
            w(i)%turb = turb_bar
        end if

        allocate( w(i)%residual(4 + nturb) )  ! [ρ, u, v, p, turb(1:nturb)]
        allocate( w(i)%grad(4 + nturb, 2 ) )  ! [ρ, u, v, p, turb(1:nturb)]
        w(i)%residual = 0.0_dp
        w(i)%grad = 0.0_dp
    end do
end subroutine initialize_flow

end module initialization
