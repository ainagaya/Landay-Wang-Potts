
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
    integer, parameter ::  q = 10, L = 40,  n_iter = 10000000, z = 4, seed = 11012000
    ! f_initial : initial value of the modification factor f
    ! f_final : final value of the modification factor f
    ! E_min : minimum energy
    ! E_max : maximum energy
    real(8), parameter :: f_initial = exp(1.0), f_final = exp(10e-8)
    ! beta_min : minimum value of the inverse temperature
    ! beta_max : maximum value of the inverse temperature
    real(8), parameter :: T_min = 1.d0, T_max = 1.d0/0.004d0, flatness=0.8

    ! Variables
    ! N : number of spins
    ! i, j : loop variables
    ! Emax, Emin : maximum and minimum energy
    integer :: N, i, j, Emax, Emin, num_E, E
    ! spins : array of spins
    ! hist : energy histogram. Number of times each energy has been visited
    ! nbr : array of nearest neighbors
    ! in : array to apply periodic boundary conditions
    integer, allocatable :: spins(:), hist(:), nbr(:,:)
    integer :: in(0:1, L)
    ! ln_n_density : logarithm of the density of states (g in the paper)
    ! ln_n_norm : normalized density of states
    real(8), allocatable :: ln_n_density(:), ln_n_norm(:), energy_density(:)
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

    ! flat : flag to check if the histogram is flat
    logical :: flat

    ! verbose output
    logical, parameter :: debug = .false.
    ! current_run : string to store the current run parameters
    character(len=8) :: current_run
    character(len=2) :: strq, strL
    character(len=8) :: strMCS

    ! random numbers generator
    real :: r1279

    N = L*L

    ! define the energy range
    ! TO-DO: drop out energies out of range, diminish Emax.
    Emin = -z*N/2
    Emax = int(Emin*0.2)
    ! Emax = 0


    num_E = abs(Emax - Emin) + 1

    allocate(spins(N))
    allocate(hist(abs(Emax):abs(Emin)))
    allocate(ln_n_density(abs(Emax):abs(Emin)))
    allocate(ln_n_norm(abs(Emax):abs(Emin)))
    allocate(energy_density(abs(Emax):abs(Emin)))
    allocate(nbr(z, N))


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

    write (strq, "(I1)")  q
    write (strL, "(I2)") L
    write (strMCS, "(I8)") n_iter

    print*, "q: ", strq
    current_run="_q"//trim(strq)//"_L"//trim(strL)

    call initialize(spins, hist, f, f_initial, ln_n_density, in, L, Emin, Emax)

    print*, "in", in(0,:)

    call generate_nbr(L, in, nbr)
    call energy(spins, nbr, L, z, E)

    print*, "Energy: ", E

    ! Main loop
    do i = 1, n_iter
     !   print*, i
        ! Each step attemps N spin flips
        do j = 1, N
            ! propose random spin flip
            call propose_flip(spins, N, L, hist, ln_n_density, f, Emax, Emin, E)
            if ((mod(E,2).ne.0).and.(q.eq.2)) then
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
        
            call check_histogram(hist, flatness, flat, Emin, Emax, num_E)
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
        E = Emin + (i-1)
        write(20,*) E, ln_n_density(abs(E))
    end do
    close(20)

    print*, "Freeing memory: spins"
    deallocate(spins)
    print*, "Freeing memory: hist"
!    deallocate(hist)
    print*, "Freeing memory: ln_n_density"
!    deallocate(ln_n_density)
    print*, "Freeing memory: energy_density"
!    deallocate(energy_density)
    print*, "Freeing memory: nbr"
!    deallocate(nbr)
    print*, "Freeing memory: in"
!    deallocate(in)


    print*, "Program finished successfully"
    print*, "Result is in file: ", 'ln_n_density' // trim(current_run) // '.dat'

contains   

subroutine initialize(spins, hist, f, f_initial, ln_n_density, in, L, Emin, Emax)
    implicit none
    integer, intent(in) :: L, Emin, Emax
    integer, intent(inout) :: spins(:), hist(abs(Emax):abs(Emin)), in(0:1,L)
    real(8), intent(inout) :: ln_n_density(abs(Emax):abs(Emin))
    real(8), intent(out) :: f
    real(8), intent(in) :: f_initial


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
    integer, intent(in) :: Emax, Emin 
    integer, intent(inout) :: spins(:), hist(abs(Emax):abs(Emin))
    real(8) :: i
    real :: r1279
    real(8), intent(inout) :: f, ln_n_density(abs(Emax):abs(Emin))
    integer :: pos
    integer :: N, L, q_new, k
    real(8) :: p
    integer, intent(inout) :: E
    integer :: E1, E2, delta_E

    ! propose a random spin flip between 0 and 1. convert the random number to a spin
    pos = mod(int(N*r1279()),N)+1

    q_new = spins(pos)

    do while (q_new.eq.spins(pos))
        q_new = mod(int(q*r1279()),q) + 1
    end do

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

    if (E2.gt.Emax) then
        ! Move rejected
    !    print*, "Im exiting because E2 > Emax"
    !    print*, E2, Emax
    !    print*, spins
  !      stop
        return
    end if

  !  print*, delta_E
    call transition_probability(E1, E2, ln_n_density, p, Emin, Emax) 

    if (p.ge.1) then
        ! Move accepted
        spins(pos) = q_new
        E = E2
    else if (p.lt.1) then
        i = r1279()
        if (i.lt.p) then
            ! Move accepted by luck
!            print*, "Move accepted by luck", i, p
            spins(pos) = q_new
            E = E2
        else
            ! Move rejected
!            print*, "Move rejected", i, p  
            E = E1
        end if

  !      if (p.eq.0) then
 !           stop
  !      end if
    end if

    ! actualiza ln_n_density 
    ln_n_density(abs(E)) = ln_n_density(abs(E)) + log(f)
    ! actualiza histograma
    hist(abs(E)) = hist(abs(E)) + 1

end subroutine

subroutine transition_probability(E1, E2, ln_n_density, p, Emin, Emax)
    implicit none
    integer, intent(in) :: E1, E2, Emin, Emax
    real(8), intent(in) :: ln_n_density(abs(Emax):abs(Emin))
    real(8), intent(out) :: p

!    print*, "E1", E1
!    print*, "E2", E2
!    print*, "ln_n_density(E1)", ln_n_density(abs(E1))
!    print*, "ln_n_density(E2)", ln_n_density(abs(E2))


    if (ln_n_density(abs(E1)) - ln_n_density(abs(E2)).gt.600) then
        p = 1.d0
    else if (ln_n_density(abs(E1)) - ln_n_density(abs(E2)).lt.-600) then
        p = 0.d0
    else
  !      print*, ln_n_density(abs(E1)) - ln_n_density(abs(E2))
        p = min(1.d0, exp(ln_n_density(abs(E1)) - ln_n_density(abs(E2))))
        
!    else
!        print*, "ln_n_density(E1)", ln_n_density(abs(E1))
!        print*, "ln_n_density(E2)", ln_n_density(abs(E2))
!        print*, "E1", E1
!        print*, "E2", E2
!        p = 0.d0
     !   print*, "p = 0"
    end if


end subroutine

subroutine generate_nbr(L, in, nbr)
    implicit none
    integer, intent(in) :: L, in(0:1,L)
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


subroutine check_histogram(hist, x, flat, Emin, Emax, num_E)
    implicit none
    integer, intent(in) :: Emin, Emax, num_E
    integer, intent(in) :: hist(abs(Emax):abs(Emin))
    integer :: k, err
    integer :: E_indx
    real(8), intent(in) :: x
    logical, intent(out) :: flat
    real(8) :: sum_hist
    real(8) :: avg_hist

    flat = .false.

    sum_hist = sum(hist)

    avg_hist = sum_hist / num_E

    err = 0
    do k = 1, num_E
        E_indx = Emin + (k-1)
        ! si algun dels valors no compleix la condicio
        if ((hist(abs(E_indx)) < (x * avg_hist)).and.(hist(abs(E_indx)) > 0)) then
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

