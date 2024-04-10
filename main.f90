
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
    integer, parameter ::  q = 4, L = 4,  n_iter = 10000000, z = 4
    ! f_initial : initial value of the modification factor f
    ! f_final : final value of the modification factor f
    logical, parameter :: debug = .false.
    ! E_min : minimum energy
    ! E_max : maximum energy
    real, parameter :: f_initial = exp(1.0), f_final = 1.0000001, E_min = 0, E_max = 100

    ! Variables
    ! N : number of spins
    integer :: N, i, j
    ! spins : array of spins
    ! hist : energy histogram. Number of times each energy has been visited
    ! f : modification factor
    ! n_density : density of states (g in the paper)
    integer, allocatable :: spins(:), hist(:), nbr(:,:), in(:,:)
    real, allocatable :: ln_n_density(:)
    real :: f, p
    real :: E, E1, E2
    logical :: flat

    N = L*L
    allocate(spins(N))
    allocate(hist(z*N/2))
    allocate(ln_n_density(z*N/2))
    allocate(nbr(z, N))
    allocate(in(0:1, N))

    call random_init(.true., .true.)

    call initialize(spins, hist, f, f_initial, ln_n_density, in, L)

    call generate_nbr(L, in, nbr)

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
            call propose_flip(spins, N, L, hist, ln_n_density, f)
        end do

        if (mod(i, 100).eq.0) then
            print*, "spins: ", spins   
            print*, "hist:", hist
            print*, "ln_n_density:", int(ln_n_density) 
            print*, "f: ", f
        
            call check_histogram(hist, 0.8, flat)
            if ( flat ) then
                print*, "Flat histogram, reducing f"
                print*, "N MCS", i
                print*, "hist", hist
                f = sqrt(f)
                hist = 0
            end if
        end if

    end do

contains   

subroutine initialize(spins, hist, f, f_initial, ln_n_density, in, L)
    implicit none
    integer, intent(inout) :: spins(:), hist(:), in(:,:)
    real, intent(inout) :: ln_n_density(:)
    real, intent(out) :: f
    real, intent(in) :: f_initial
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

subroutine propose_flip(spins, N, L, hist, ln_n_density, f)
    implicit none
    integer, intent(inout) :: spins(:), hist(:)
    real :: i
    real, intent(inout) :: f, ln_n_density(:)
    integer :: x, y, pos, energy_index
    integer :: N, L, q_new, k
    real :: delta_E, E1, E2, p

    ! propose a random spin flip
    call random_number(i)

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
        call random_number(i)
        q_new = int(i*q) + 1
    end do


    call energy(spins, nbr, L, z, E1)

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
        ! Fer que depengui de la temperatura?
        call random_number(i)
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
    real, intent(in) :: E1, E2
    real, intent(in) :: ln_n_density(:)
    real, intent(out) :: p
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
	real, intent(out) :: E

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
    real, intent(in) :: E
    integer, intent(out) :: i

    ! REVISIT THIS
    do i = 0, int(z*N/2.d0)
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
    real, intent(in) :: x
    logical, intent(out) :: flat
    real :: sum_hist
    real :: avg_hist

    flat = .false.

    sum_hist = sum(hist)

    avg_hist = sum_hist / size(hist)
    flat = all(hist >= x * avg_hist)


!    if ( debug ) then
        print*, "sum_hist: ", sum_hist
        print*, "avg_hist: ", avg_hist
        print*, "x * avg_hist: ", x * avg_hist
        print*, "flat: ", flat
        print*, "hist: ", hist
        print*, "array: ", hist >= x * avg_hist

!    end if

    if (.not. flat) then
        print*, "Percentage of flatness: ", count(hist >= x * avg_hist) / real(size(hist)) * 100.d0
    end if


end subroutine

end program LandauWangPotts
