
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

program LandauWangPotts
    Implicit none
    ! Initial parameters
    ! q : number of states a spin can take
    ! L : lattice size
    ! n_iter : number of Monte Carlo steps
    ! z : number of nearest neighbors
    ! num_beta : number of beta values
    integer, parameter ::  q = 2, L = 20,  n_iter = 100000, z = 4, seed = 11012000,  num_T = 500
    ! f_initial : initial value of the modification factor f
    ! f_final : final value of the modification factor f
    ! E_min : minimum energy
    ! E_max : maximum energy
    real(8), parameter :: f_initial = exp(1.0), f_final = 1.01, E_min = 0, E_max = 100
    ! beta_min : minimum value of the inverse temperature
    ! beta_max : maximum value of the inverse temperature
    real(8), parameter :: T_min = 1.d0, T_max = 1.d0/0.004d0, flatness=0.8

    ! Variables
    ! N : number of spins
    ! i, j : loop variables
    ! Emax, Emin : maximum and minimum energy
    integer :: N, i, j, Emax, Emin, k
    ! spins : array of spins
    ! hist : energy histogram. Number of times each energy has been visited
    ! nbr : array of nearest neighbors
    ! in : array to apply periodic boundary conditions
    integer, allocatable :: spins(:), hist(:), nbr(:,:), in(:,:)
    ! ln_n_density : logarithm of the density of states (g in the paper)
    ! ln_n_norm : normalized density of states
    real(8), allocatable :: ln_n_density(:), ln_n_norm(:), energy_density(:)
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
    real(8) :: E, E1, E2, ln_A, beta, internal_en, h, free_en, S, T

    ! flat : flag to check if the histogram is flat
    logical :: flat

    ! verbose output
    logical, parameter :: debug = .false.
    ! current_run : string to store the current run parameters
    character(len=20) :: current_run

    character(len=2) :: strq, strL
    character(len=8) :: strMCS

    ! random numbers generator
    real :: r1279

    N = L*L
    allocate(spins(N))
    allocate(hist(z*N/2))
    allocate(ln_n_density(z*N/2))
    allocate(ln_n_norm(z*N/2))
    allocate(energy_density(z*N/2))
    allocate(nbr(z, N))
    allocate(in(0:1, N))

    ! call random_init(.true., .true.)
    call setr1279(seed)

    print*, "Initializing system..."
    print*, "Parameters of the run: "
    print*, "q: ", q
    print*, "L: ", L
    print*, "n_iter: ", n_iter
    print*, "f_initial: ", f_initial
    print*, "f_final: ", f_final
    print*, "E_min: ", E_min
    print*, "E_max: ", E_max
    print*, "T_min: ", T_min
    print*, "T_max: ", T_max

    write (strq, "(I2)")  q
    write (strL, "(I2)") L
    write (strMCS, "(I8)") n_iter

    print*, "strq: ", strq
    current_run = "_q" // trim(strq) // "_L" // trim(strL) // "_niter" // trim(strMCS)

    call initialize(spins, hist, f, f_initial, ln_n_density, in, L)

    call generate_nbr(L, in, nbr)

    ! define the energy range
    ! TO-DO: drop out energies out of range, diminish Emax.
    Emin = -z*N/2
    Emax = 0

    if (debug) then
        print*, "in1:", in(1,:)
        print*, "in0:", in(0,:)
        print*, "nbr1", nbr(1,:)
        print*, "nbr2", nbr(2,:)
        print*, "nbr3", nbr(3,:)
        print*, "nbr4", nbr(4,:)
    end if

    call energy(spins, nbr, L, z, E)

    ! Main loop
    do i = 1, n_iter
        ! Each step attemps N spin flips
        do j = 1, N
            ! propose random spin flip
            call propose_flip(spins, N, L, hist, ln_n_density, f, Emax, Emin)
        end do

        if (mod(i, 1000).eq.0) then
            print*, "------------------------------------"
            print*, "Iteration: ", i
            print*, "spins: ", spins   
            print*, "hist:", hist
            print*, "ln_n_density:", int(ln_n_density) 
            print*, "f: ", f
        
            call check_histogram(hist, flatness, flat)
            if ( flat ) then
                print*, "Flat histogram, reducing f and resetting histogram..."
                print*, "N MCS", i
                print*, "hist", hist
                f = sqrt(f)
                print*, "f: ", f
                ! Reset histogram
                hist = 0

                if (f < f_final) then
                    print*, "f < f_final"
                    print*, "Achieved the desired f in", i, "MCS"
                    exit
                end if
            else
                print*, "Not flat histogram"
                print*, "Current f: ", f
            end if
        end if

    end do

    print*, "Final density of states: ", ln_n_density

!   From the final density of states, estimate the internal energy for different temperatures. By rescaling n(E) so that n(Egroundstate) = q, 
!   we can also obtain the free energy F = and thus the entropy. Check that internal energy, entropy, and free
!   energy densities agree (for different temperatures) with the exact results of Ferdinand and
!   Fisher, that are computed by the program ferdinand.f written by Bernd Berg.   
    
    print*, "Normalizing density of states..."

    ln_A = log(real(q)) - ln_n_density(z*N/2)

    print*, log(real(q))
    print*, (ln_n_density(z*N/2))
    print*, "ln_A: ", ln_A

    do i = 1, z*N/2

        if (ln_n_density(i) > 0) then
            ln_n_norm(i) = ln_n_density(i) + ln_A
        else
            ln_n_norm(i) = 0
        end if
    
    end do

    print*, "ln_n_norm: ", ln_n_norm

    print*, "Writing output to file..."

    open(20, file='ln_n_density' // trim(current_run) // '.dat')
    write(20, '(A)') "# E/N ln_n_density    hist"
    h = abs(Emax - Emin) / size(ln_n_norm)

    do i = 0, z*N/2
        E = Emin + i * h
        call e_index(E, j)
        write(20,*) E/N, ln_n_norm(j), hist(j)
    end do
    close(20)

    print*, "Calculating internal energy, free energy and entropy..."

    h = (T_max - T_min) / num_T

    open(10, file='res_ising' // current_run//'.dat')
    write(10, '(A)') "#beta    internal_energy/N   free_energy/N   entropy/N"
    T = T_min
    do i = 0, num_T
        T = T_min + h*i
        beta = 1/T
        call internal_energy(ln_n_norm, beta, internal_en, Emin, Emax)
        print*, "Internal energy at beta = ", beta, " is: ", internal_en, "per site", internal_en/N
        call free_energy(ln_n_norm, beta, free_en, Emin, Emax)
        print*, "Free energy at beta = ", beta, " is: ", free_en, "per site", free_en/N
        call entropy(internal_en, free_en, beta, S)
        print*, "Entropy at beta = ", beta, " is: ", S, "per site", S/N
        call compute_energy_density(ln_n_norm, beta, E, energy_density)
        write(10,*) beta, internal_en/N, free_en/N, S/N
        do k = 1, z*N/2
            write(10,*) k, energy_density(k)
        end do
    end do
    close(10)


contains   

subroutine initialize(spins, hist, f, f_initial, ln_n_density, in, L)
    implicit none
    integer, intent(inout) :: spins(:), hist(:), in(:,:)
    real(8), intent(inout) :: ln_n_density(:)
    real(8), intent(out) :: f
    real(8), intent(in) :: f_initial
    integer, intent(in) :: L

    ! initialize the lattice
    spins = 1

    ! initialize the energy histogram
    hist = 0

    f = f_initial

    ! initialize the density of states
    ln_n_density = 0

    ! initialize the array to apply periodic boundary conditions
    do i = 1, L
		in(0,i) = i-1
		in(1,i) = i+1
	end do
	
	in(0,1) = L
	in(1,L) = 1


end subroutine

subroutine propose_flip(spins, N, L, hist, ln_n_density, f, Emax, Emin)
    implicit none
    integer, intent(inout) :: spins(:), hist(:)
    integer, intent(in) :: Emax, Emin 
    real(8) :: i
    real :: r1279
    real(8), intent(inout) :: f, ln_n_density(:)
    integer :: x, y, pos, energy_index
    integer :: N, L, q_new, k
    real(8) :: delta_E, E1, E2, p

    ! propose a random spin flip
    ! call random_number(i)
    i = r1279()

    ! Mapping: fent servir un sol random number decidim x i y NOOOOOOO
    pos = int(i*N) + 1
    if ( debug ) then
        print*, "i", i, "pos: ", pos
    end if

    ! no cal perque ja ho tenim en pla, vector in
    ! x = (pos - 1) / L + 1
    ! y = mod(pos - 1, L) + 1

    ! biased?
    q_new = spins(pos)

    do while (q_new.eq.spins(pos))
        ! call random_number(i)
        i = r1279()
        q_new = int(i*q) + 1
    end do


    call energy(spins, nbr, L, z, E1)

    if ( debug ) then
        print*, "current energy: ", E1
    end if

    delta_E = 0
    do k = 1, z
        if (spins(nbr(k,pos)).eq.q_new) then
            delta_E = delta_E - 1
        else if (spins(nbr(k,pos)).eq.spins(pos)) then
            delta_E = delta_E + 1
        end if
    end do

    if ( debug ) then
        print*, "spins", spins
    end if

    if ( debug ) then
        print*, "delta_E: ", delta_E
    end if

    E2 = E1 + delta_E

    call transition_probability(E1, E2, ln_n_density, p)

    if ( debug ) then
        print*, "p: ", p
    end if    

    if (p.ge.1) then
        spins(pos) = q_new
        if ( debug ) then 
            print*, "move accepted"
        end if 
        call e_index(E2, energy_index)
    else if (p.lt.1) then
        ! call random_number(i)
        i = r1279()
        if (i.lt.p) then
            spins(pos) = q_new
            if ( debug ) then
                print*, "move accepted by luck"
            end if
            call e_index(E2, energy_index)
        else
            if ( debug ) then
                print*, "move NOT accepted"
            end if
            call e_index(E1, energy_index)
        end if
    end if

    ! actualiza ln_n_density 
    ln_n_density(energy_index) = ln_n_density(energy_index) + log(f)
    ! actualiza histograma
    hist(energy_index) = hist(energy_index) + 1



end subroutine

subroutine transition_probability(E1, E2, ln_n_density, p)
    implicit none
    real(8), intent(in) :: E1, E2
    real(8), intent(in) :: ln_n_density(:)
    real(8), intent(out) :: p
    integer :: i1, i2

    if ( debug ) then 
        print*, "Looking for index of E1: ", E1
    end if    
    call e_index(E1, i1)

    if ( debug ) then 
        print*, "Looking for index of E2: ", E2
    end if
    call e_index(E2, i2)

    if ( debug ) then
        print*, "E1: ", E1, i1, ln_n_density(i1)
        print*, "E2: ", E2, i2, ln_n_density(i2)
    end if


    p = min(1.d0, exp(ln_n_density(i1) - ln_n_density(i2)))

    if ( debug ) then
        print*, "p: ", p
    end if


end subroutine

subroutine generate_nbr(L, in, nbr)
    implicit none
    integer, intent(in) :: L, in(:,:)
    integer, intent(out) :: nbr(:, :)
    integer :: x, y, i

    i = 0
    do y = 1, L
		do x = 1, L
			i=i+1
			nbr(1,i) = in(1,x) + L*(y-1) ! ‘‘right’’
			nbr(2,i) = in(0,x) + L*(y-1) ! ‘‘left’’
			nbr(3,i) = x + L*(in(1,y)-1) ! ‘‘up’’
			nbr(4,i) = x + L*(in(0,y)-1) ! ‘‘down’’
		end do
	end do


end subroutine



Subroutine energy(s, nbr, L, z, E)
	Implicit none
	integer, intent(in) :: L, z
	integer, intent(in) :: s(:)
	integer, intent(in) :: nbr(:,:)
	integer :: i, k
	real(8), intent(out) :: E

	E = 0
	do i = 1, L*L
		do k = 1, z
            if (s(i).eq.s(nbr(k,i))) then
                ! Puting 1/2 because we are counting twice the energy of each pair
                E = E - 1.d0/2.d0
            end if
		end do
	end do


End Subroutine


Subroutine magnetization(s, L, M)
	Implicit none
	integer :: L, i
	integer, dimension(L*L), intent(in) :: s
	real(8), intent(out) :: M

	M = 0
	do i = 1, L*L
		M = M + s(i)
	end do

End Subroutine

subroutine e_index(E, i)
    implicit none
    real(8), intent(in) :: E
    integer, intent(out) :: i

    ! REVISIT THIS
    do i = 0, int(z*N/2.d0)
!        print*, "E", E
        if (int(E).eq.(-i)) then
            return
        end if
    end do

    print*, "did not found index for E: ", E
    print*, "spins", spins(:)
    print*, "f", f
    stop

end subroutine

subroutine check_histogram(hist, x, flat)
    implicit none
    integer, intent(in) :: hist(:)
    integer :: k, err
    real(8), intent(in) :: x
    logical, intent(out) :: flat
    real(8) :: sum_hist
    real(8) :: avg_hist

    flat = .false.

    sum_hist = sum(hist)

    avg_hist = sum_hist / size(hist > 0)

    err = 0
    do k = 1, size(hist)
        ! si algun dels valors no compleix la condicio
        if ((hist(k) < (x * avg_hist)).and.(hist(k) > 0)) then
            err = err + 1
        end if
    end do

    if (err.eq.0) then
        print*, "Flat is true"
        print*, "hist: ", hist
        print*, "x * avg_hist: ", x * avg_hist
        print*, "array: ", hist >= x * avg_hist
        flat = .true.
    else
        flat = .false.
        print*, "Percentage of flatness: ", err / real(size(hist > 0 )) * 100.d0
    end if

!    if (debug) then
        print*, "sum_hist: ", sum_hist
        print*, "avg_hist: ", avg_hist
        print*, "x * avg_hist: ", x * avg_hist
        print*, "flat: ", flat
        print*, "hist: ", hist
        print*, "array: ", hist >= x * avg_hist
!    end if

    return
end subroutine


subroutine internal_energy(ln_n_norm, beta, internal_en, Emin, Emax)
    implicit none
    real(8), intent(in) :: ln_n_norm(:)
    real(8), dimension(size(ln_n_norm)) :: ln_n_norm_minus
    real(8), intent(in) :: beta
    integer :: Emin, Emax, i ,k
    real(8), intent(out) :: internal_en
    real(8) :: numerator, denominator, lambda, h, E

    numerator = 0
    denominator = 0

    h = abs(Emax - Emin) / size(ln_n_norm)
    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        ln_n_norm_minus(i) = ln_n_norm(i) - beta*E
    end do

    lambda = maxval(ln_n_norm_minus)

    h = abs(Emax - Emin) / size(ln_n_norm)

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        call e_index(E, k)

        if (debug) then
            print*, E
            print*, "ln_n_norm(E)", ln_n_norm(i), E, beta, lambda
        end if

        numerator = numerator + exp(ln_n_norm(k) - beta*E - lambda) * E 
        denominator = denominator + exp(ln_n_norm(k) - beta*E - lambda)
    end do

    print*, "numerator: ", numerator
    print*, "denominator: ", denominator    

    internal_en = numerator / denominator

    return

end subroutine

subroutine free_energy(ln_n_norm, beta, free_en, Emin, Emax)
    implicit none
    real(8), intent(in) :: ln_n_norm(:)
    real(8), dimension(size(ln_n_norm)) :: ln_n_norm_minus
    real(8), intent(in) :: beta
    integer :: Emin, Emax, i, k
    real(8), intent(out) :: free_en
    real(8) :: numerator, h, E, lambda

    numerator = 0

    h = abs(Emax - Emin) / size(ln_n_norm)

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        ln_n_norm_minus(i) = ln_n_norm(i) - beta*E
    end do

    lambda = maxval(ln_n_norm_minus)

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        call e_index(E, k)
        numerator = numerator + exp ( ln_n_norm(k) - beta*E - lambda) 
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

subroutine compute_energy_density(ln_n_norm, beta, E, energy_density)
    implicit none
    real(8), intent(in) :: ln_n_norm(:)
    real(8), intent(in) :: beta
    real(8), intent(out) :: E
    integer :: Emin, Emax, i, k
    real(8) :: h
    real(8), intent(out) :: energy_density(:)

    h = abs(Emax - Emin) / size(ln_n_norm)

    do i = 1, size(ln_n_norm)
        E = Emin + (i-1) * h
        call e_index(E, k)
        energy_density(i) = ln_n_norm(k) * exp(-beta*E)
    end do

end subroutine

end program LandauWangPotts

