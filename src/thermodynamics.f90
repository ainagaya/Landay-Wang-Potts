
! Compute several thermodynamic variables coming from a Landau-Wang simulation with the Potts model
program thermodynamics
    Implicit none
    ! Initial parameters
    ! q : number of states a spin can take
    ! L : lattice size
    ! n_iter : number of Monte Carlo steps
    ! z : number of nearest neighbors
    ! num_beta : number of beta values
    integer :: q, L, seed
    real(8) :: flatness
    integer, parameter :: z = 4
    integer ::  num_T
    ! beta_min : minimum value of the inverse temperature
    ! beta_max : maximum value of the inverse temperature
    real(8):: T_min, T_max
    ! Variables
    ! N : number of spins
    ! i, k : loop variables
    ! Emax, Emin : maximum and minimum energy
    ! E, : store energy values
    integer :: N, i, k, Emax, Emin, num_E, E
    ! ln_n_density : logarithm of the density of states (g in the paper)
    ! ln_n_norm : normalized density of states
    ! energy_density : energy density values
    ! energy_array : array to store energy values
    real(8), allocatable :: ln_n_density(:), ln_n_norm(:), energy_density(:), energy_array(:)
    ! A : constant to normalize the density of states
    ! beta : inverse temperature
    ! internal_en : internal energy value
    ! free_en : free energy value
    ! S : entropy value
    ! T : temperature value
    ! T_interval : interval between temperatures
    ! specific_heat : specific heat value
    ! beta_min : minimum value of beta
    ! beta_interval : interval between beta values
    real(8) :: A, beta, internal_en, free_en, S, T, T_interval, specific_heat, beta_min, beta_interval
    ! current_run : string to store the current run parameters
    character(len=8) :: current_run
    ! strq, strL : strings to store the values of q and L
    character(len=2) :: strq, strL

    ! Read parameters from namelist file
    namelist /LWparams/ q, L, seed, flatness
    namelist /Thermo/ T_min, T_max, num_T

    open(unit=10, file='LWparams.nml', status='old')
    read(10, nml=LWparams)
    close(10)

    open(unit=10, file='LWparams.nml', status='old')
    read(10, nml=Thermo)
    close(10)

    N = L*L

    Emin = -z*N/2
    Emax = int(Emin*0.2)

    num_E = abs(Emax - Emin) + 1

    allocate(ln_n_density(abs(Emax):abs(Emin)))
    allocate(ln_n_norm(abs(Emax):abs(Emin)))
    allocate(energy_density(abs(Emax):abs(Emin)))
    allocate(energy_array(abs(Emax):abs(Emin)))

    print*, "Parameters of the run: "
    print*, "q: ", q
    print*, "L: ", L
    print*, "T_min: ", T_min
    print*, "T_max: ", T_max

    write (strq, "(I2)")  q
    write (strL, "(I2)") L

    current_run = "_q" // trim(strq) // "_L" // trim(strL)
    open(20, file='results/ln_n_density' // trim(current_run) // '.dat', status='old')

    do i = 1, num_E
        E = Emin + (i-1)
        read(20,*) energy_array(abs(E)), ln_n_density(abs(E))
    end do

    print*, "energy array:", energy_array
  
    print*, "Normalizing density of states..."

    A = log(real(q)) - ln_n_density(abs(Emin))

    print*, "A: ", A

    do i = 1, num_E
        E = Emin + (i-1)

        if (ln_n_density(abs(E)) > 0) then
            ln_n_norm(abs(E)) = ln_n_density(abs(E)) + A
        else
            ln_n_norm(abs(E)) = 0
        end if
    
    end do

    open(20, file='ln_n_density' // trim(current_run) // 'norm.dat')

    print*, "Writing normalized density of states to file..."

    do i = 1, num_E
        E = Emin + (i-1)
        write(20,*) real(E)/N, ln_n_norm(abs(E))
    end do
    close(20)

    print*, "Calculating internal energy, free energy and entropy... (temperatures equiespaiades)"

    T_interval = abs(T_max - T_min) / num_T

    open(10, file='res' // current_run//'.dat')
    open(11, file='energy_density' // current_run//'.dat')
!    write(10, '(A)') "#beta    internal_energy/N   free_energy/N   entropy/N"
    T = T_min
    do i = 0, num_T
        T = T_min + T_interval*i
        beta = 1/T
        print*, "T", T
        call internal_energy_specific_heat(ln_n_norm, beta, internal_en, Emin, Emax, specific_heat)
        print*, "Internal energy at beta = ", beta, " is: ", internal_en, "per site", internal_en/N
        call free_energy(ln_n_norm, beta, free_en, Emin, Emax)
        print*, "Free energy at beta = ", beta, " is: ", free_en, "per site", free_en/N
        call entropy(internal_en, free_en, beta, S)
        print*, "Entropy at beta = ", beta, " is: ", S, "per site", S/N
        write(10,*) beta, internal_en/N, free_en/N, S/N, specific_heat/N
        if ((T<1).and.(T>0.5)) then
            call compute_energy_density(ln_n_norm, beta, E, Emin, Emax, energy_density)
            do k = 1, num_E
                E = Emin + (k-1)
                write(11,*) 1.d0/beta, energy_density(abs(E)), real(E)/N
            end do
        end if

    end do
    close(10)
    close(11)

    stop
    print*, "Calculating internal energy, free energy and entropy... (betes equiespaiades)"

    beta_interval = (1/T_min - 1/T_max) / num_T

    open(12, file='res_ising' // current_run//'beta.dat')
    open(13, file='energy_density' // current_run//'beta.dat')
!    write(10, '(A)') "#beta    internal_energy/N   free_energy/N   entropy/N"
    beta_min = 1/T_max
    beta = beta_min
    do i = 0, num_T
        beta = beta_min + beta_interval*i
        call internal_energy_specific_heat(ln_n_norm, beta, internal_en, Emin, Emax, specific_heat)
        print*, "Internal energy at beta = ", beta, " is: ", internal_en, "per site", internal_en/N
        call free_energy(ln_n_norm, beta, free_en, Emin, Emax)
        print*, "Free energy at beta = ", beta, " is: ", free_en, "per site", free_en/N
        call entropy(internal_en, free_en, beta, S)
        print*, "Entropy at beta = ", beta, " is: ", S, "per site", S/N
        write(12,*) beta, internal_en/N, free_en/N, S/N, specific_heat/N

    end do
    close(12)

    print*, "Done!"

    deallocate(ln_n_norm)
    deallocate(energy_density)
    deallocate(energy_array)



contains   



subroutine internal_energy_specific_heat(ln_n_norm, beta, internal_en, Emin, Emax, specific_heat)
    implicit none
    integer :: Emin, Emax, i
    real(8), intent(in) :: ln_n_norm(abs(Emax):abs(Emin))
    real(8), dimension(abs(Emax):abs(Emin)) :: ln_n_norm_minus
    real(8), intent(in) :: beta
    real(8), intent(out) :: internal_en, specific_heat
    real(8) :: numerator, denominator, lambda, numerator_2
    integer :: E

    numerator = 0
    numerator_2 = 0
    denominator = 0

    do i = 1, size(ln_n_norm_minus)
        E = Emin + (i-1) 
        ln_n_norm_minus(abs(E)) = ln_n_norm(abs(E)) - beta*E
    end do

    lambda = maxval(ln_n_norm_minus)

    do i = 1, size(ln_n_norm_minus)
        E = Emin + (i-1) 

        numerator = numerator + exp(ln_n_norm_minus(abs(E)) - lambda) * E 
        numerator_2 = numerator_2 + exp(ln_n_norm_minus(abs(E)) - lambda) * E**2
        denominator = denominator + exp(ln_n_norm_minus(abs(E)) - lambda)
    end do   

    internal_en = numerator / denominator

    specific_heat = (numerator_2 / denominator - (numerator / denominator)**2) * beta**2

    return

end subroutine

subroutine free_energy(ln_n_norm, beta, free_en, Emin, Emax)
    implicit none
    integer :: Emin, Emax, i
    real(8), intent(in) :: ln_n_norm(abs(Emax):abs(Emin))
    real(8), dimension(abs(Emax):abs(Emin)) :: ln_n_norm_minus
    real(8), intent(in) :: beta
    real(8), intent(out) :: free_en
    real(8) :: numerator, lambda
    integer :: E

    numerator = 0

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1)
        ln_n_norm_minus(abs(E)) = ln_n_norm(abs(E)) - beta*E
    end do

    lambda = maxval(ln_n_norm_minus)

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) 
        numerator = numerator + exp ( ln_n_norm_minus(abs(E)) - lambda) 
    end do    

    free_en = -1/beta * (lambda + log (numerator))

    return

end subroutine

subroutine entropy(U, F, beta, S)
    implicit none
    real(8), intent(in) :: U, F, beta
    real(8), intent(out) :: S

    S = (U - F) * beta

end subroutine

subroutine compute_energy_density(ln_n_norm, beta, E, Emin, Emax, energy_density)
    implicit none
    real(8), intent(in) :: beta
    integer, intent(out) :: E
    integer :: Emin, Emax, i
    real(8), intent(out) :: energy_density(abs(Emax):abs(Emin))
    real(8), intent(in) :: ln_n_norm(abs(Emax):abs(Emin))
    real(8) :: lambda, ln_n_norm_minus(abs(Emax):abs(Emin))


    do i = 1, size(ln_n_norm)
        E = Emin + (i-1)
        ln_n_norm_minus(abs(E)) = ln_n_norm(abs(E)) - beta*E
    end do

    lambda = maxval(ln_n_norm_minus)

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1)
        if (ln_n_norm(abs(E)).ne.0) then
            energy_density(abs(E)) = exp( ln_n_norm(abs(E)) - beta*E - lambda)
!           print*, "energy_density(abs(E))", energy_density(abs(E)), "ln_n_norm(abs(E)- betaE)", ln_n_norm(abs(E)) - beta*E
!            print*, "ln_n_norm(abs(E))", ln_n_norm(abs(E)), "beta*E", -beta*E
        end if
    end do

end subroutine


end program thermodynamics

