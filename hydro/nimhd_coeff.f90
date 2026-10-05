module nimhd_coeff_module
  use amr_parameters, only: dp
  use hydro_parameters, only: nvar
  implicit none
contains
  !-----------------------------------------------------------------
  ! Local resistivities of the split non-ideal MHD solver, given the
  ! conservative cell state uu and the cell-centred magnetic field bb.
  ! The Ohmic EMF is -eo*J and the ambipolar EMF is
  ! ea*[(J.B)B-B^2*J]. By default they are the namelist values etamag
  ! (or eta_o) and eta_a. Edit this routine to make the coefficients
  ! depend on density, temperature, ionization or field strength.
  !-----------------------------------------------------------------
#ifdef _CUDA
  attributes(host,device) &
#endif
  subroutine nimhd_coeff(uu,bb,gamma,etamag,eta_a,eo,ea)
    real(dp),dimension(1:nvar),intent(in)::uu
    real(dp),dimension(1:3),intent(in)::bb
    real(dp),value::gamma,etamag,eta_a
    real(dp),intent(out)::eo,ea

    eo=etamag
    ea=eta_a

  end subroutine nimhd_coeff
end module nimhd_coeff_module
