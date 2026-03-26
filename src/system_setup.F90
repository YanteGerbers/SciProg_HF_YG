module system_setup

use molecular_structure
use ao_basis
implicit none

contains

    subroutine define_molecule(molecule)
    ! This subroutine read data from the input.xyz file
      type(molecular_structure_t), intent(inout) :: molecule
      integer :: i, natoms, Z 
      character(len=2) :: element
      real(8), parameter :: ang_to_bohr = 1.8897D0
      real(8), allocatable :: charge(:), coord(:,:)

      open(unit=10, file='input.xyz', status='old')
      read(10,*) ! skip comment lines
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) natoms
      read(10,*) 
      allocate(charge(natoms))
      allocate(coord(3,natoms))
      do i=1,natoms
        read(10,*) element, Z, coord(1,i), coord(2,i), coord(3,i)
        coord(:,i) = coord(:,i) * ang_to_bohr
        charge(i) = Z 
      end do
      close(10)
      call add_atoms_to_molecule(molecule, charge, coord)
   end subroutine


   subroutine define_basis(ao_basis, molecule)
    ! This routine was extended to use more extensive basis sets, according to assigment,
    ! and modified to be consistent with define_molecule.
    ! The coordinates of the shell centers are the nuclear coordinates (of the atoms)
     type(basis_set_info_t), intent(inout) :: ao_basis
     type(molecular_structure_t), intent(in) :: molecule
     !!!!! type(basis_func_info_t) :: gto
     integer :: i
     do i = 1, molecule%num_atoms
       if (molecule%charge(i) == 1) then
         ! H: 3 uncontracted s-shells:    l          coord      exp     3x1=3funcs  (1func/shell)
         call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),0.1D0)
         call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),1.0D0)
         call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),3.0D0)
       else
         ! non-H:
         !    5 uncontracted s-shells:    l          coord      exp      5x1=5funcs (1func/shell)
         call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),0.10D0)
         call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),0.35D0)
         call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),1.0D0)
         call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),3.0D0)
         call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),10.0D0)
         !    3 uncontracted p-shells:    l          coord      exp      3x3=9funcs (3func/shell, px,py,pz)
         call add_shell_to_basis(ao_basis,1,molecule%coord(:,i),0.2D0)
         call add_shell_to_basis(ao_basis,1,molecule%coord(:,i),1.0D0)
         call add_shell_to_basis(ao_basis,1,molecule%coord(:,i),5.0D0)
         !    1 uncontracted d-shell:     l          coord      exp      1x6=6funcs (6funcs/shell, xx, yy, zz, xy, xz, yz) 
         call add_shell_to_basis(ao_basis,2,molecule%coord(:,i),1.0D0)
       end if
     end do
   end subroutine

end module system_setup