
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
    integer, parameter ::  q = 2, L = 20,  n_iter = 1000000, z = 4, seed = 11012000,  num_T = 500
    ! f_initial : initial value of the modification factor f
    ! f_final : final value of the modification factor f
    ! E_min : minimum energy
    ! E_max : maximum energy
    real(8), parameter :: f_initial = exp(1.0), f_final = exp(10e-7)
    ! beta_min : minimum value of the inverse temperature
    ! beta_max : maximum value of the inverse temperature
    real(8), parameter :: T_min = 1.d0, T_max = 1.d0/0.004d0, flatness=0.8

    ! Variables
    ! N : number of spins
    ! i, j : loop variables
    ! Emax, Emin : maximum and minimum energy
    integer :: N, i, j, Emax, Emin, k, num_E, E
    ! spins : array of spins
    ! hist : energy histogram. Number of times each energy has been visited
    ! nbr : array of nearest neighbors
    ! in : array to apply periodic boundary conditions
    integer, allocatable :: spins(:), hist(:), nbr(:,:), in(:,:)
    ! ln_n_density : logarithm of the density of states (g in the paper)
    ! ln_n_norm : normalized density of states
    real(8), allocatable :: ln_n_density(:), ln_n_norm(:), energy_density(:), energy_array(:)
    ! f : modification factor
    ! p : probability of accepting a move
    real(8) :: f
    ! E, E1, E2 : store energy values
    ! ln_A : constant to normalize the density of states
    ! beta : inverse temperature
    ! internal_en : internal energy value
    ! h : step size
    ! free_en : free energy value
    ! S : entropy value
    real(8) :: h

    ! flat : flag to check if the histogram is flat
    logical :: flat

    ! verbose output
    logical, parameter :: debug = .false.
    ! current_run : string to store the current run parameters
    character(len=23) :: current_run

    character(len=2) :: strq, strL
    character(len=8) :: strMCS

    ! random numbers generator
    real :: r1279

    N = L*L

    ! define the energy range
    ! TO-DO: drop out energies out of range, diminish Emax.
    Emin = -z*N/2
    Emax = 0

    num_E = abs(Emax - Emin) + 1
    ! h : step size
    h = 1

    allocate(spins(N))
    allocate(hist(num_E))
    allocate(ln_n_density(num_E))
    allocate(ln_n_norm(num_E))
    allocate(energy_density(num_E))
    allocate(nbr(z, N))
    allocate(in(0:1, N))
    allocate(energy_array(num_E))

    ! Build energy array
    do i = 1, num_E 
        energy_array(i) = Emin + (i-1) * h
    end do

    print*, "energy_array: ", energy_array
    

    ! call random_init(.true., .true.)
    call setr1279(seed)

    print*, "Initializing system..."
    print*, "Parameters of the run: "
    print*, "q: ", q
    print*, "L: ", L
    print*, "n_iter: ", n_iter
    print*, "f_initial: ", f_initial
    print*, "f_final: ", f_final
    print*, "E_min: ", Emin
    print*, "E_max: ", Emax
    print*, "T_min: ", T_min
    print*, "T_max: ", T_max

    write (strq, "(I2)")  q
    write (strL, "(I2)") L
    write (strMCS, "(I8)") n_iter

    print*, "q: ", strq
    current_run = "_q" // trim(strq) // "_L" // trim(strL) // "_niter" // trim(strMCS)

    call initialize(spins, hist, f, f_initial, ln_n_density, in, L)
    call generate_nbr(L, in, nbr)
    call energy(spins, nbr, L, z, E)

    ! Main loop
    do i = 1, n_iter
        ! Each step attemps N spin flips
        do j = 1, N
            ! propose random spin flip
            call propose_flip(spins, N, L, hist, ln_n_density, f, Emax, Emin, E)
            if (mod(E,2).ne.0) then
                print*, "ERROR: E is not even"
                print*, "spins", spins
                print*, "E", E
                stop
            end if
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

    if (i.ge.n_iter) then
        print*, "Did not achieve the desired f in", n_iter, "MCS"
        stop
    end if

    print*, "Final density of states: ", ln_n_density

!   From the final density of states, estimate the internal energy for different temperatures. By rescaling n(E) so that n(Egroundstate) = q, 
!   we can also obtain the free energy F = and thus the entropy. Check that internal energy, entropy, and free
!   energy densities agree (for different temperatures) with the exact results of Ferdinand and
!   Fisher, that are computed by the program ferdinand.f written by Bernd Berg.   


    open(20, file='ln_n_density' // trim(current_run) // '.dat')

    print*, "Writing density of states to file..."

    do i = 1, num_E
        E = Emin + (i-1) * h
        call e_index(E, energy_array, j)
        write(20,*) E, ln_n_density(j)
    end do
    close(20)

    print*, "Freeing memory: spins"
    deallocate(spins)
    print*, "Freeing memory: hist"
    deallocate(hist)
    print*, "Freeing memory: ln_n_density"
    deallocate(ln_n_density)
    print*, "Freeing memory: energy_density"
    deallocate(energy_density)
    print*, "Freeing memory: nbr"
!    deallocate(nbr)
    print*, "Freeing memory: in"
!    deallocate(in)
    print*, "Freeing memory: energy_array"
    deallocate(energy_array)

    print*, "Program finished successfully"
    print*, "Result is in file: ", 'ln_n_density' // trim(current_run) // '.dat'

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

subroutine propose_flip(spins, N, L, hist, ln_n_density, f, Emax, Emin, E)
    implicit none
    integer, intent(inout) :: spins(:), hist(:)
    integer, intent(in) :: Emax, Emin 
    real(8) :: i
    real :: r1279
    real(8), intent(inout) :: f, ln_n_density(:)
    integer :: x, y, pos, energy_index
    integer :: N, L, q_new, k
    real(8) :: p
    integer, intent(inout) :: E
    integer :: E1, E2, delta_E

    ! propose a random spin flip between 0 and 1
    i = r1279()

    ! convert the random number to a spin
    pos = int(i*N) + 1

    q_new = spins(pos)

    do while (q_new.eq.spins(pos))
        i = r1279()
        q_new = int(i*q) + 1
    end do


    !!!!! PROBLEMS HERE FOR SMALL F !!!!!
    call energy(spins, nbr, L, z, E1)

    if (E1.ne.E) then
        print*, "spins", spins
        print*, "q_new", q_new
        print*, "pos", pos
        print*, "ln_n_density", ln_n_density
        call e_index(E1, energy_array, energy_index)
        print*, "energy_index", energy_index
        print*, "H(i)", hist(energy_index)
        print*, "ln_n_density(i)", ln_n_density(energy_index)

        print*, "ERROR: E does not match E1" 
        stop
    end if

    E1 = E

    delta_E = 0
    do k = 1, z
        if (spins(nbr(k,pos)).eq.q_new) then
            delta_E = delta_E - 1
        else if (spins(nbr(k,pos)).eq.spins(pos)) then
            delta_E = delta_E + 1
        end if
    end do

    E2 = E1 + delta_E

    call transition_probability(E1, E2, ln_n_density, p) 

    if (p.ge.1) then
        ! Move accepted
        spins(pos) = q_new
        call e_index(E2, energy_array, energy_index)
        E = E2
    else if (p.lt.1) then
        i = r1279()
        if (i.lt.p) then
            ! Move accepted by luck
            spins(pos) = q_new
            call e_index(E2, energy_array, energy_index)
            E = E2
        else
            ! Move rejected
            call e_index(E1, energy_array, energy_index)
            E = E1
        end if
    end if

    ! actualiza ln_n_density 
    ln_n_density(energy_index) = ln_n_density(energy_index) + log(f)
    ! actualiza histograma
    hist(energy_index) = hist(energy_index) + 1

end subroutine

subroutine transition_probability(E1, E2, ln_n_density, p)
    implicit none
    integer, intent(in) :: E1, E2
    real(8), intent(in) :: ln_n_density(:)
    real(8), intent(out) :: p
    integer :: i1, i2
  
    call e_index(E1, energy_array, i1)
    call e_index(E2, energy_array, i2)

    p = min(1.d0, exp(ln_n_density(i1) - ln_n_density(i2)))


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
	integer, intent(out) :: E

	E = 0
	do i = 1, L*L
		do k = 1, z
            if (s(i).eq.s(nbr(k,i))) then
                E = E - 1
            end if
		end do
	end do

 ! Puting 1/2 because we are counting twice the energy of each pair
    E = E / 2

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
        flat = .true.
    else
        print*, "Flat is false"
        flat = .false.
        print*, "Percentage of flatness: ", ( 1- err / real(size(hist > 0 ))) * 100.d0
    end if

    return
end subroutine

end program LandauWangPotts

