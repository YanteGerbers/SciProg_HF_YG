module system_setup

use molecular_structure
use ao_basis
implicit none

contains

    subroutine define_molecule(molecule)
     ! This routine should be improved such that an arbitrary molecule can be given as input
     ! the coordinates below are for a be-he dimer oriented along the x-axis with a bond length of 2 au
     type(molecular_structure_t), intent(inout) :: molecule
     real(8) :: charge(2),coord(3,2)
     charge(1)   = 4.D0
     charge(2)   = 2.D0
     coord       = 0.D0
     coord(1,2)  = 2.D0
     call add_atoms_to_molecule(molecule,charge,coord)
   end subroutine

   subroutine define_basis(ao_basis)
    ! This routine can be extended to use better basis sets 
    ! The coordinates of the shell centers are the nuclear coordinates
    ! Think of a refactoring of define_molecule and define_basis to ensure consistency 
     type(basis_set_info_t), intent(inout) :: ao_basis
     type(basis_func_info_t) :: gto
     ! Be:  2 uncontracted s-funs:    l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),4.D0)
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),1.D0)
     ! He:  1 uncontracted s-fun:     l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/2.D0,0.D0,0.D0/),1.D0)
   end subroutine

   subroutine initialize_system(molecule, ao_basis, Ne, Nocc, NAO)
    ! This routine .................. 
    type(molecular_structure_t), intent(out) :: molecule
    type(basis_set_info_t), intent(out) :: ao_basis
    integer, intent(out) :: Ne, Nocc, NAO
    ! Define molecule (hardcoded)
    call define_molecule(molecule)
    ! Define basis
    call define_basis(ao_basis)
    ! Define number of basis functions
    NAO = ao_basis%nao
    ! Define number of electrons
    Ne = int(sum(molecule%charge))
    ! Define number of occupied orbitals
    Nocc = Ne / 2
  end subroutine


end module system_setup