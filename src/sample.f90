 subroutine sample( step, time, aver, flx, timeStep, switch )

  ! Purpose: Write accumulated averages to output.dat

  use global
  implicit none

  type(averages), intent(in) :: aver
  type(fluxes), intent(in) :: flx


  integer, intent(in) :: step, switch

  real(8) :: totEnergy, potEnergy, kinEnergyWall, kinEnergyFluid
  real(8) :: tempFluid, tempWall, flux1!, flux2
  real(8), intent(in) :: time, timeStep


  open(unit=7,file='output.dat',status='old',position='append')

  if( switch == 0 ) then

    write(unit=7,fmt="(//'       step',2x,'       time',2x,'  totEnergy',2x,'  potEnergy',2x,&
                     'kinEnergyFl',2x,'kinEnergyWa',2x,'  tempFluid',2x,'   tempWall',       &
                     2x,'     Flux_1',2x,'     Flux_2')") 
    write(unit=7,fmt="('===================================================================&
    &========================================================',/)") 

    close(unit=7)
    return
  
  end if

  totEnergy = aver%totEnergy(1)/aver%counter
  potEnergy = aver%potEnergy(1)/aver%counter
  kinEnergyFluid = aver%kinEnergyFluid(1)/aver%counter
  kinEnergyWall  = aver%kinEnergyWall(1)/aver%counter
  tempFluid = aver%tempFluid(1)/aver%counter
  tempWall  = aver%tempWall(1)/aver%counter
  flux1 =  (flx%rightFlux(1) -flx%leftFlux(1))/timeStep/flx%counter
!  flux2 =  (flx%rightFlux(2) -flx%leftFlux(2))/timeStep/flx%counter



  write(unit=7,fmt="(i11,2x,f11.5,7(2x,e11.4))") step, time, totEnergy, potEnergy,&
                    kinEnergyFluid, kinEnergyWall, tempFluid, tempWall, flux1!, flux2

  close(unit=7)

 end subroutine sample
