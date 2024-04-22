
! 1. Implement the Landau-Wang algorithm for the 2D q-state Potts model as described in the references
!    above. Algorithm:
!    1. set f = finitial and n(E) = 1 for E ∈ [Emin, Emax]
!    2. propose random spin flips, accepting them with probability min(1, n(E)/n(E0)) where
!       E, E0 are the energies before and after the spin flip (reject if E0 falls outside the range
!       [Emin, Emax]).
!    3. each time a move is accepted, modify n(E0) by multiplying it by f (in practice use the
!       logarithm), and update the energy histogram.
!    4. every niter MCS (1 MCS = N attempted spin flips), if the energy histogram is “sufficiently flat”,
!       reduce f as f → √f and reset the histogram.
!    5. stop when f has reached a predetermined value such as 1.0000001, or when the computational time exceeds
!       some reasonable limit.

program thermodynamics
    Implicit none
    ! Initial parameters
    ! q : number of states a spin can take
    ! L : lattice size
    ! n_iter : number of Monte Carlo steps
    ! z : number of nearest neighbors
    ! num_beta : number of beta values
    integer, parameter ::  q = 2, L = 10,  n_iter = 100000, z = 4, seed = 11012000,  num_T = 500
    ! f_initial : initial value of the modification factor f
    ! f_final : final value of the modification factor f
    ! E_min : minimum energy
    ! E_max : maximum energy
    real(8), parameter :: f_initial = exp(1.0), f_final = exp(10e-7)
    ! beta_min : minimum value of the inverse temperature
    ! beta_max : maximum value of the inverse temperature
    real(8), parameter :: T_min = 1/2.d0, T_max = 1.d0/0.002d0, flatness=0.95

    ! Variables
    ! N : number of spins
    ! i, j : loop variables
    ! Emax, Emin : maximum and minimum energy
    integer :: N, i, j, Emax, Emin, k, num_E, E1, E2, E
    ! ln_n_density : logarithm of the density of states (g in the paper)
    ! ln_n_norm : normalized density of states
    real(8), allocatable :: ln_n_density(:), ln_n_norm(:), energy_density(:), energy_array(:)
    ! f : modification factor
    ! p : probability of accepting a move
    real(8) :: f, p
    ! E, E1, E2 : store energy values
    ! ln_A : constant to normalize the density of states
    ! beta : inverse temperature
    ! internal_en : internal energy value
    ! h : step size
    ! free_en : free energy value
    ! S : entropy value
    real(8) :: ln_A, beta, internal_en, h, free_en, S, T, T_interval, specific_heat, beta_min, beta_interval

    ! flat : flag to check if the histogram is flat
    logical :: flat

    ! verbose output
    logical, parameter :: debug = .false.
    ! current_run : string to store the current run parameters
    character(len=20) :: current_run

    character(len=2) :: strq, strL
    character(len=8) :: strMCS

    N = L*L

    ! define the energy range
    ! TO-DO: drop out energies out of range, diminish Emax.
    Emin = -z*N/2
    Emax = 0

    num_E = abs(Emax - Emin) + 1
    ! h : step size
    h = 1

    allocate(ln_n_density(num_E))
    allocate(ln_n_norm(num_E))
    allocate(energy_density(num_E))
    allocate(energy_array(num_E))

    ! Build energy array
    do i = 1, num_E 
        energy_array(i) = Emin + (i-1) * h
    end do

    print*, "energy_array: ", energy_array

    print*, "Initializing system..."
    print*, "Parameters of the run: "
    print*, "q: ", q
    print*, "L: ", L
    print*, "n_iter: ", n_iter
    print*, "T_min: ", T_min
    print*, "T_max: ", T_max

    write (strq, "(I2)")  q
    write (strL, "(I2)") L
    write (strMCS, "(I8)") n_iter

    current_run = "_q" // trim(strq) // "_L" // trim(strL) // "_niter" // trim(strMCS)

    open(20, file='ln_n_density' // trim(current_run) // '.dat', status='old')

    do i = 1, num_E
        read(20,*) energy_array(i), ln_n_density(i)
    end do
  
    print*, "Normalizing density of states..."

    ln_A = log(real(q)) - ln_n_density(num_E)

    print*, log(real(q))
    print*, (ln_n_density(num_E))
    print*, "ln_A: ", ln_A

    do i = 1, num_E

        if (ln_n_density(i) > 0) then
            ln_n_norm(i) = ln_n_density(i) + ln_A
        else
            ln_n_norm(i) = 0
        end if
    
    end do

    open(20, file='ln_n_density' // trim(current_run) // 'norm.dat')

    print*, "Writing normalized density of states to file..."

    do i = 1, num_E
        E = Emin + (i-1) * h
        call e_index(E, energy_array, j)
        write(20,*) E/N, ln_n_norm(j)
    end do
    close(20)

    print*, "Calculating internal energy, free energy and entropy... (temperatures equiespaiades)"

    T_interval = (T_max - T_min) / num_T

    open(10, file='res_ising' // current_run//'.dat')
!    write(10, '(A)') "#beta    internal_energy/N   free_energy/N   entropy/N"
    T = T_min
    do i = 0, num_T
        T = T_min + T_interval*i
        beta = 1/T
        call internal_energy_specific_heat(ln_n_norm, beta, internal_en, Emin, h, specific_heat)
        print*, "Internal energy at beta = ", beta, " is: ", internal_en, "per site", internal_en/N
        call free_energy(ln_n_norm, beta, free_en, Emin, h)
        print*, "Free energy at beta = ", beta, " is: ", free_en, "per site", free_en/N
        call entropy(internal_en, free_en, beta, S)
        print*, "Entropy at beta = ", beta, " is: ", S, "per site", S/N
        call compute_energy_density(ln_n_norm, beta, E, Emin, h, energy_density)
        write(10,*) beta, internal_en/N, free_en/N, S/N, specific_heat/N
        do k = 1, z*N/2
            E = Emin + (k-1) * h
            write(11,*) 1.d0/beta, energy_density(k), E
        end do
        write(11,*) ""
        write(11,*) ""
    end do
    close(10)

    print*, "Calculating internal energy, free energy and entropy... (betes equiespaiades)"

    beta_interval = (1/T_min - 1/T_max) / num_T

    open(12, file='res_ising' // current_run//'beta.dat')
!    write(10, '(A)') "#beta    internal_energy/N   free_energy/N   entropy/N"
    beta_min = 1/T_max
    beta = beta_min
    do i = 0, num_T
        beta = beta_min + beta_interval*i
        call internal_energy_specific_heat(ln_n_norm, beta, internal_en, Emin, h, specific_heat)
        print*, "Internal energy at beta = ", beta, " is: ", internal_en, "per site", internal_en/N
        call free_energy(ln_n_norm, beta, free_en, Emin, h)
        print*, "Free energy at beta = ", beta, " is: ", free_en, "per site", free_en/N
        call entropy(internal_en, free_en, beta, S)
        print*, "Entropy at beta = ", beta, " is: ", S, "per site", S/N
        call compute_energy_density(ln_n_norm, beta, E, Emin, h, energy_density)
        write(12,*) beta, internal_en/N, free_en/N, S/N, specific_heat/N
        do k = 1, z*N/2
            E = Emin + (k-1) * h
            write(11,*) 1.d0/beta, energy_density(k), E
        end do
        write(11,*) ""
        write(11,*) ""
    end do
    close(12)

    close(11)





contains   


subroutine e_index(E, energy_array, i)
    implicit none
    integer, intent(in) :: E
    integer, intent(out) :: i
    real(8), intent(in) :: energy_array(:)

    ! REVISIT THIS
    do i = 1, num_E
        if (energy_array(i).eq.E) then
            return
        end if
    end do

    print*, "did not found index for E: ", E
    print*, "f", f
    stop

end subroutine


subroutine internal_energy_specific_heat(ln_n_norm, beta, internal_en, Emin, h, specific_heat)
    implicit none
    real(8), intent(in) :: ln_n_norm(:)
    real(8), dimension(size(ln_n_norm)) :: ln_n_norm_minus
    real(8), intent(in) :: beta
    integer :: Emin, i ,k
    real(8), intent(out) :: internal_en, specific_heat
    real(8) :: numerator, denominator, lambda, h, numerator_2
    integer :: E

    numerator = 0
    numerator_2 = 0
    denominator = 0

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        ln_n_norm_minus(i) = ln_n_norm(i) - beta*E
    end do

    lambda = maxval(ln_n_norm_minus)

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        call e_index(E, energy_array, k)

        if (debug) then
            print*, E
            print*, "ln_n_norm(E)", ln_n_norm(i), E, beta, lambda
        end if

        numerator = numerator + exp(ln_n_norm_minus(k) - lambda) * E 
        numerator_2 = numerator_2 + exp(ln_n_norm_minus(k) - lambda) * E**2
        denominator = denominator + exp(ln_n_norm_minus(k) - lambda)
    end do

    print*, "numerator: ", numerator
    print*, "denominator: ", denominator    

    internal_en = numerator / denominator

    specific_heat = (numerator_2 / denominator - (numerator / denominator)**2) * beta**2

    return

end subroutine

subroutine free_energy(ln_n_norm, beta, free_en, Emin, h)
    implicit none
    real(8), intent(in) :: ln_n_norm(:)
    real(8), dimension(size(ln_n_norm)) :: ln_n_norm_minus
    real(8), intent(in) :: beta
    integer :: Emin, i, k
    real(8), intent(out) :: free_en
    real(8) :: numerator, h, lambda
    integer :: E

    numerator = 0

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        ln_n_norm_minus(i) = ln_n_norm(i) - beta*E
    end do

    lambda = maxval(ln_n_norm_minus)

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        call e_index(E, energy_array, k)
        numerator = numerator + exp ( ln_n_norm_minus(k) - lambda) 
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

subroutine compute_energy_density(ln_n_norm, beta, E, Emin, h, energy_density)
    implicit none
    real(8), intent(in) :: ln_n_norm(:)
    real(8), intent(in) :: beta
    integer, intent(out) :: E
    integer :: Emin, i, k
    real(8) :: h
    real(8), intent(out) :: energy_density(:)

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        call e_index(E, energy_array, k)
        if (ln_n_norm(k).ne.0) then
            energy_density(i) = exp(ln_n_norm(k)) * exp(-beta*E)
        end if
    end do

end subroutine


end program thermodynamics

