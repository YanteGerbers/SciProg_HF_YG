! Hartree-Fock program
! Computes integrals over Gaussian basis functions to solve matrix eigenvalue problems and compute Hartree Fock energy
! Yante Gerbers, March 2026

program HartreeFock
   
   use molecular_structure
   use ao_basis
   use system_setup
   use compute_integrals 
   use scf_loop

   implicit none

     ! System 
     type(molecular_structure_t) :: molecule
     type(basis_set_info_t) :: ao_basis
     integer :: n_AO, n_occ, n_e
    
     ! Matrices 
     real(8), allocatable :: H_core(:,:)  
     real(8), allocatable :: S(:,:), T(:,:), V(:,:)
     real(8), allocatable :: ao_integrals (:,:,:,:)

     ! Results
     real(8), allocatable :: eps(:)
     real(8)  :: E_HF

     integer  :: i

     ! Define system
     call define_molecule(molecule)
     call define_basis(ao_basis, molecule)
     call system_size(ao_basis, molecule, n_AO, n_e, n_occ)

     ! One-electron integrals (overlap-, kinetic- and potential matrix)
     allocate (S(n_AO,n_AO))                               
     call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)
     
     allocate (T(n_AO,n_AO))                   
     call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)

     allocate (V(n_AO,n_AO))                        
     call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)

     ! Two-electron integrals
     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
     call generate_2int (ao_basis,ao_integrals)
       
     ! Core Hamiltonian
     allocate(H_core(n_AO,n_AO))
     H_core = T - V               ! V is electrostatic potential; scaled by electron charge (-1) gives potential energy -V  
       
     ! Run SCF 
     call run_scf(molecule, H_core, S, ao_integrals, n_occ, eps, E_HF)

     ! Output
     write(*,*) '--------------------------------------'
     write(*,*) 'Hartree-Fock Results'
     write(*,*) '--------------------------------------'

     write(*,'(a30,i5)') 'Number of electrons: ', n_e
     write(*,'(a30,i5,/)') 'Number of basis functions: ', n_AO

     write(*,*) 'Orbital energies (Ha):'
       do i = 1, n_AO
         write(*,'(i3,2x,f12.6)') i, eps(i)
       end do

     write(*,'(/,a,f12.6,a)') 'Hartree-Fock energy: ', E_HF, ' Ha'


end program HartreeFock