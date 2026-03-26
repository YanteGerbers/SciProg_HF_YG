module scf_convergence

implicit none
private
public :: conv_check_eps

contains

  function conv_check_eps (eps, eps_old, threshold) result(converged)

    real(8), intent(in) :: eps(:), eps_old(:)
    real(8), intent(in) :: threshold
    logical :: converged
    real(8) :: delta_eps
    delta_eps = sqrt(sum((eps - eps_old)**2))       ! Euclidian norm
    if (delta_eps < threshold) then
      converged = .true.
    else
      converged = .false.
    end if

  end function  
   

end module 

