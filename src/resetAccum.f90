subroutine resetAccum( aver, flx, switch)

  !Purpose: Reset the accumulated averages

  use global
  implicit none

  type(averages), intent(out) :: aver
  type(fluxes), intent(out) :: flx
  integer, intent(in) :: switch

  ! Reset observables
  aver%totEnergy = 0.0d0 
  aver%potEnergy = 0.0d0
  aver%kinEnergyFluid = 0.0d0
  aver%kinEnergyWall = 0.0d0
  aver%tempFluid = 0.0d0
  aver%tempWall  = 0.0d0
  aver%counter = 0.0d0

  ! Reset profiles and gradients
  aver%velProfile  = 0.0d0
  aver%tempProfile = 0.0d0
  aver%rhoProfile  = 0.0d0
  aver%counterProfile = 0
  aver%totalCounterProfile = 0

  aver%tempGradient = 0.0d0
  aver%rhoGradient  = 0.0d0
  aver%counterGradient = 0
  aver%totalCounterGradient = 0

  flx%counter = 0.0d0
  flx%leftFlux = 0.0d0
  flx%rightFlux = 0.0d0

  if (switch == 1) then
    open(unit=7,file='output.dat',status='old',action='write',position='append')
    write(unit=7,fmt="(//,5x, 'Equilibration period has finished')")
    close(unit=7)
  end if

end subroutine resetAccum
