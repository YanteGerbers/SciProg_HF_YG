! SCF loop module
! Contains routines for SCF iterations in the HartreeFock program
! Yante Gerbers, March 2026

module scf_loop
    
use compute_integrals
use diagonalization

implicit none
private
public :: run_scf, build_density, build_fock, check_conv_eps, nuclear_repulsion



contains

  subroutine run_scf(molecule, H_core, S, ao_integrals, n_occ, eps, E_HF)
    use molecular_structure

  ! Self consistent field (SCF) procedure:
  ! iteratively update Density matrix and Fock matrix until convergence
  !
  ! Input:
  ! molecule      - molecular structure
  ! H_core        - core Hamiltonian
  ! S             - overlap matrix
  ! ao_integrals  - two-electron integrals
  ! n_occ         - number of occupied orbitals
  !
  ! Output:
  ! eps           - orbital energies (eigenvalues)
  ! E_HF          - total Hartree-Fock energy

    type(molecular_structure_t), intent(in) :: molecule
    real(8), intent(in) :: H_core(:,:), S(:,:), ao_integrals(:,:,:,:)
    integer, intent(in) :: n_occ
    real(8), allocatable, intent(out) :: eps(:)
    real(8), intent(out) :: E_HF
    
    integer :: n_AO, iter, max_iter
    real(8) :: threshold_eps, E_HF_elec, V_nn, alfa
    real(8), allocatable :: Fock(:,:), Density(:,:), Density_old(:,:), C(:,:), eps_old(:)
        
    n_AO = size(H_core,1)
    alfa = 0.25D0                   ! Damping parameter

    allocate(Fock(n_AO,n_AO))     
    allocate(Density(n_AO,n_AO))  
    allocate(Density_old(n_AO,n_AO))
    allocate(C(n_AO,n_AO))          
    allocate(eps(n_AO))              
    allocate(eps_old(n_AO))
    
    max_iter = 50
    threshold_eps = 1.0D-6
    Density_old = Density

    ! Initial guess: start with core Hamiltionion (equivalent to Density = 0)
    Fock = H_core              
    call solve_genev (Fock,S,C,eps)
    eps_old = eps 
    
    ! SCF iteration loop
    do iter = 1, max_iter  

      Density_old = Density
      call build_density(C, n_occ, Density)
      ! Density mixing (damping) to improve SCF convergence and stability
      if (iter > 1) then
        Density = alfa*Density + (1.0D0 - alfa)*Density_old
      endif
  
      call build_fock(H_core, Density, ao_integrals, Fock)

      ! Orbitals (C) and energies (eps) are updated
      call solve_genev(Fock,S,C,eps)  

      ! Convergence based on change in MO energies
      write(*,'(a,i3,2x,a,f12.6)') 'iter: ', iter, 'delta eps', sqrt(sum((eps-eps_old)**2))
      if (check_conv_eps(eps, eps_old, threshold_eps)) exit
      eps_old = eps

    end do

    ! Report convergence status
    if (iter == max_iter) then
      write(*,*) 'WARNING: SCF did not converge'
    else
      write(*,'(a,i5)') 'SCF iterations: ', iter
    end if

    ! Final Hartee Fock energy: electronic + nuclear repulsion
    E_HF_elec = sum((H_core + Fock) * Density)
    call nuclear_repulsion(molecule, V_nn)
    E_HF = E_HF_elec + V_nn
    
    deallocate(Fock, Density, C, eps_old)

  end subroutine


  subroutine build_density(C, n_occ, Density)
  ! Builds density matrix from occupied orbitals

    real(8), intent(in) :: C(:,:)
    integer, intent(in) :: n_occ
    real(8), intent(out) :: Density(:,:)

    integer :: kappa, lambda
    integer :: n_AO

    n_AO = size(Density,1)

    do kappa = 1, n_AO
      do lambda = 1, kappa     
        Density(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))  ! only occupied orbitals with lowest energy contribute
        Density(lambda,kappa) = Density(kappa,lambda)                    ! Symmetry
      end do
    end do
  
  end subroutine


  subroutine build_fock(H_core, Density, ao_integrals, Fock)
  ! Builds Fock matrix: core hamiltonial + electron-electron interactions

    real(8), intent(in) :: H_core(:,:), Density(:,:), ao_integrals(:,:,:,:)
    real(8), intent(out) :: Fock(:,:)

    integer :: kappa, lambda, mu, nu
    integer :: n_AO

    n_AO = size(Fock,1)
    
    do kappa = 1, n_AO
      do lambda = 1, kappa                                      
        Fock(kappa,lambda) = H_core(kappa,lambda)
        do mu = 1, n_AO
          do nu = 1, n_AO                         
            Fock(kappa,lambda) = Fock(kappa,lambda) + Density(mu,nu) * (2.0D0 * ao_integrals(kappa,lambda,mu,nu) - ao_integrals(kappa,nu,mu,lambda))
          end do
        end do
        Fock(lambda,kappa) = Fock(kappa,lambda)        ! Symmetry
        end do
    end do

  end subroutine

  
  Function check_conv_eps(eps, eps_old, threshold_eps) result(converged)
  ! Checks SCF convergence based on change in orbital energies (Euclidean norm)
  
    real(8), intent(in) :: eps(:), eps_old(:)
    real(8), intent(in) :: threshold_eps

    logical :: converged
    real(8) :: delta_eps

    delta_eps = sqrt(sum((eps - eps_old)**2))
    
    if (delta_eps < threshold_eps) then
      converged = .true.
    else
      converged = .false.
    end if

  end Function  


  subroutine nuclear_repulsion(molecule, V_nn)
  ! Computes nuclear-nuclear repulsion energy from atomic positions and charges 

    use molecular_structure

    type(molecular_structure_t), intent(in) :: molecule
    real(8), intent(out) :: V_nn
    
    integer :: a,b
    real(8) :: distance_ab
   
    V_nn = 0.0D0
    
    do a = 1, molecule%num_atoms
      do b = a+1, molecule%num_atoms
        distance_ab = sqrt(sum((molecule%coord(:,a) - molecule%coord(:,b))**2)) 
        V_nn = V_nn + molecule%charge(a) * molecule%charge(b) / distance_ab
      end do
    end do

  end subroutine

end module