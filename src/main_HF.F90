program HartreeFock

   ! Program modified from main_demo as starting point
   ! Yante Gerbers, March 2026

   use molecular_structure
   use ao_basis
   use compute_integrals
   use diagonalization
   implicit none

     ! Variable containing the molecular structure
     type(molecular_structure_t) :: molecule
     ! Variable containing the atomic orbital basis
     type(basis_set_info_t) :: ao_basis

     ! Variable naming as in the description of the exercise
     integer  :: n_AO, n_occ, n_e
     integer  :: kappa, lambda
     integer  :: mu, nu
     real(8)  :: E_old, E_new, V_nn, threshold
     integer  :: iter, max_iter 
     real(8)  :: E_HF
     real(8), allocatable :: F(:,:),V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:)
     real(8), allocatable :: H_core(:,:)
     !!!!! The following large array can be eliminated when Fock matrix contruction is implemented
     real(8), allocatable :: ao_integrals (:,:,:,:) 

     ! Define molecule
     call define_molecule(molecule)
     ! Define GTOs/basis functions
     call define_basis(ao_basis, molecule)
     ! Define number of atomic orbitals
     n_AO = ao_basis%nao
     ! Define number of electrons
     n_e = int(sum(molecule%charge))
     ! Defie number of occupied orbitals
     n_occ = n_e / 2
     
     ! Compute the overlap matrix       (eq. 9)
     allocate (S(n_AO,n_AO))          
     call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)
     ! Compute the kinetic matrix       (eq. 7)
     allocate (T(n_AO,n_AO))
     call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)
     ! Compute the potential matrix     (eq. 8)
     allocate (V(n_AO,n_AO))
     call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)
     ! Compute all 2-electron integrals (eq. 5)
     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
     call generate_2int (ao_basis,ao_integrals)
     ! Compute the core Hamiltonian matrix (the potential is positive, we scale with -e = -1 to get to the potential energy matrix)
     allocate (F(n_AO,n_AO))
     allocate(H_core(n_AO,n_AO))
     H_core = T - V                   ! (eq. 6)
     F = H_core                       ! Initial guess for F in first SCF iteration (equivalent to D=0)
     
     !!!!! Initial diagonalization (with D=0 and F=h_core)
     ! Form the density matrix D and set to 0
     allocate (D(n_AO,n_AO))
     D = 0.0D0
     allocate (C(n_AO,n_AO))
     allocate (eps(n_AO))
     call solve_genev (F,S,C,eps)
     print*, "Orbital energies for the core Hamiltonian:",eps
     
     !!!!! SCF LOOP
     max_iter = 50
     threshold = 1.0D-6
     E_old = 0.0D0 
     do iter = 1, max_iter
        ! Density matrix D (with C, eq. 2)
        do lambda = 1, n_AO
           do kappa = 1, n_AO
              D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
           end do
        end do
        ! Fock matrix F (with h_core, D, g. eq. 3)
        do lambda = 1, n_AO
           do kappa = 1, n_AO
              F(kappa,lambda) = H_core(kappa,lambda)
              do mu = 1, n_AO
                 do nu = 1, n_AO
                    F(kappa,lambda) = F(kappa,lambda) + D(mu,nu) * (2.0D0 * ao_integrals(kappa,lambda,mu,nu) - ao_integrals(kappa,nu,mu,lambda))
                 end do
              end do
           end do
        end do
        ! Compute the Hartree-Fock energy
        E_new = sum((H_core + F) * D)      ! (eq. 4)
        ! Convergence check
        if (abs(E_new - E_old) < threshold) exit
        E_old = E_new
        ! Diagonalize new F
        call solve_genev(F,S,C,eps)
    end do


    ! Compute final HF energy incl. nuclear repulsion Vnn (eq. 4)
    !!!!! Calculate V_nn
    V_nn = 0.0D0
    !
    E_HF = E_new + V_nn
    
   
    ! Print results
    print*, "Final orbital energies:"                   ,eps
    print*, "Final Hartree-Fock energy:                ", E_HF
    print*, "SCF iterations:"                           , iter
    if (iter == max_iter) print*, "problem: no convergence"

end program HartreeFock