module nimhd_coeff_module
  use amr_parameters, only: dp
  use hydro_parameters, only: nvar
  implicit none
contains
  !-----------------------------------------------------------------
  ! Local resistivities of the split non-ideal MHD solver, given the
  ! conservative cell state uu and the cell-centred magnetic field bb.
  ! The Ohmic EMF is -eta_o*J and the ambipolar EMF is
  ! eta_a*[(J.B)B-B^2*J]. By default they are the namelist values
  ! etamag (or eta_ohm) and eta_ad. Edit this routine to make the
  ! coefficients depend on density, temperature, ionization or field
  ! strength. They are evaluated once per split step and held over its
  ! substeps, which is exact for coefficients that are constant or depend
  ! on density; the B^2 of the ambipolar EMF follows the current field.
  !-----------------------------------------------------------------
#ifdef _CUDA
  attributes(host,device) &
#endif
  subroutine nimhd_coeff(uu,bb,gamma,etamag,eta_ad,eta_o,eta_a)
    real(dp),dimension(1:nvar),intent(in)::uu
    real(dp),dimension(1:3),intent(in)::bb
    real(dp),value::gamma,etamag,eta_ad
    real(dp),intent(out)::eta_o,eta_a

    eta_o=etamag
    eta_a=eta_ad

  end subroutine nimhd_coeff
end module nimhd_coeff_module
