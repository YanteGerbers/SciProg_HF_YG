program HartreeFock

   ! Program modified from main_demo as starting point
   ! Yante Gerbers, March 2026
   
   use molecular_structure
   use ao_basis
   use system_setup
   use compute_integrals
   use diagonalization 
   use scf_convergence

   implicit none

     ! Variable containing the molecular structure
     type(molecular_structure_t) :: molecule
     ! Variable containing the atomic orbital basis
     type(basis_set_info_t) :: ao_basis

     ! Variable naming as in the description of the exercise
     integer  :: n_AO, n_occ, n_e
     integer  :: kappa, lambda, mu, nu
     integer  :: iter, max_iter
     real(8), allocatable :: H_core(:,:)
     real(8), allocatable :: F(:,:), V(:,:), T(:,:), S(:,:), C(:,:), eps(:), D(:,:)
     real(8), allocatable :: eps_old(:)
     real(8)  :: threshold_eps
     integer  :: a, b
     real(8)  :: d_ab, V_nn, E_HF_elec, E_HF
     integer  :: i

     !!!!! The following large array can be eliminated when Fock matrix contruction is implemented
     real(8), allocatable :: ao_integrals (:,:,:,:)

     ! Specify system size
     call define_molecule(molecule)
     call define_basis(ao_basis, molecule)
     n_AO = ao_basis%nao
     n_e = int(sum(molecule%charge))
     n_occ = n_e / 2
     
     ! Compute integrals
     allocate (S(n_AO,n_AO))                       ! Overlap matrix S          
     call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)
     allocate (T(n_AO,n_AO))                       ! Kinetic matrix T
     call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)
     allocate (V(n_AO,n_AO))                       ! Potential matrix V 
     call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)
     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))  ! 2-electron integrals
     call generate_2int (ao_basis,ao_integrals)
       
     ! Compute core Hamiltonian matrix
     allocate(H_core(n_AO,n_AO))
     H_core = T - V               ! V is electrostatic potential; scaled by electron charge (-1) gives potential energy -V  
     
     ! Allocate matrices
     allocate (F(n_AO,n_AO))      ! Fock matrix 
     allocate (C(n_AO,n_AO))      ! MO coefficients (eigenvectors)
     allocate (eps(n_AO))         ! MO energies (eigenvalues)
     allocate (D(n_AO,n_AO))      ! Density matrix      
     
     ! SCF LOOP
     max_iter = 50
     threshold_eps = 1.0D-6
     ! Initial guess
     F = H_core                   ! equivalent to D=0  
     call solve_genev (F,S,C,eps)
     eps_old = eps 
     
     do iter = 1, max_iter

        ! Density matrix D
        do kappa = 1, n_AO
           do lambda = 1, kappa                                          ! Symmetry 
              D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))  ! 1:n_occ --> n_occ eigenvectors with lowest energy
              D(lambda,kappa) = D(kappa,lambda)                          ! Symmetry
           end do
        end do

        ! Fock matrix F
        do kappa = 1, n_AO
          do lambda = 1, kappa                                           ! Symmetry
            F(kappa,lambda) = H_core(kappa,lambda)
              do mu = 1, n_AO
                do nu = 1, n_AO                         
                  F(kappa,lambda) = F(kappa,lambda) + D(mu,nu) * (2.0D0 * ao_integrals(kappa,lambda,mu,nu) - ao_integrals(kappa,nu,mu,lambda))
                end do
              end do
          F(lambda,kappa) = F(kappa,lambda)                              ! Symmetry
          end do
        end do

        ! Diagonalize F --> new orbitals (C, eigenvectors) and orbital energies (eps, eigenvalues)
        call solve_genev(F,S,C,eps)

        ! Convergence check: until MO energies (ε) are stable between iterations  
        if (conv_check_eps(eps, eps_old, threshold_eps)) exit
        eps_old = eps

    end do

    ! Compute final HF energy incl. nuclear repulsion (Vnn)
    E_HF_elec = sum((H_core + F) * D)
    V_nn = 0.0D0
    do a = 1, molecule%num_atoms
       do b = a+1, molecule%num_atoms
          d_ab = sqrt(sum((molecule%coord(:,a) - molecule%coord(:,b))**2))  ! Euclidian norm
          V_nn = V_nn + molecule%charge(a) * molecule%charge(b) / d_ab
       end do
    end do
    E_HF = E_HF_elec + V_nn
   

    ! Print results
    write(*,*) '--------------------------------------'
    write(*,*) 'Hartree-Fock Results'
    write(*,*) '--------------------------------------'

    write(*,'(a30,i5)') 'Number of electrons: ', n_e
    write(*,'(a30,i5,/)') 'Number of basis functions: ', n_AO
    write(*,*) 'Orbital energies:'
      do i = 1, n_AO
        write(*,'(i3,2x,f12.6)') i, eps(i)
      end do

    write(*,'(/,a,f12.6,a)') 'Final Hartree-Fock energy: ', E_HF, ' Ha'
    write(*,'(a,i5)') 'SCF iterations: ', iter
    if (iter == max_iter) write(*,*) 'problem: no convergence'

end program HartreeFock