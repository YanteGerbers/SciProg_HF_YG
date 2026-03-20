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
     
     ! Compute the overlap matrix
     allocate (S(n_AO,n_AO))
     call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)

     ! Compute the kinetic matrix
     allocate (T(n_AO,n_AO))
     call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)

     ! Compute the potential matrix
     allocate (V(n_AO,n_AO))
     call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)

     ! Compute the core Hamiltonian matrix (the potential is positive, we scale with -e = -1 to get to the potential energy matrix)
     allocate (F(n_AO,n_AO))
     allocate(H_core(n_AO,n_AO))
     H_core = T - V
     F = H_core ! for first SCF iteration 

     ! Diagonalize the Fock matrix
     allocate (C(n_AO,n_AO))
     allocate (eps(n_AO))
     call solve_genev (F,S,C,eps)
     print*, "Orbital energies for the core Hamiltonian:",eps

     ! Form the density matrix
     allocate (D(n_AO,n_AO))
     do lambda = 1, n_ao
        do kappa = 1, n_ao
           D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
       end do
     end do

     !!!!! Compute the Hartree-Fock energy (this should be modified, see the notes)
     E_HF = 2.D0 * sum(F*D)
     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
     ! Compute all 2-electron integrals
     call generate_2int (ao_basis,ao_integrals)
     do lambda = 1, n_ao
        do kappa = 1, n_ao
           E_HF = E_HF + 2.D0 *  D(kappa,lambda) * sum(D*ao_integrals(:,:,kappa,lambda))
           E_HF = E_HF - 1.D0 *  D(kappa,lambda) * sum(D*ao_integrals(:,lambda,kappa,:))
       end do
     end do
   
     print*, "The Hartree-Fock energy:    ", E_HF

   end



   
end program HartreeFock