subroutine evolutionProperties( time, prop, switch )

  ! Purpose: Write instantaneous system properties to stats.dat

  use global
  implicit none

  type(properties), intent(in) :: prop

  integer, intent(in) :: switch

  real(8), intent(in) :: time


  select case( switch ) 

  case( 0 )


    open(unit=10,file='stats.dat',status='unknown',action='write')

    write(unit=10,fmt="('#      time',2x,'  totEnergy',2x,'  potEnergy',2x,&
                       ' kinEneryFl',2x,'kinEnergyWa',2x,'  tempFluid',2x,'   tempWall')") 

    close(unit=10)

  case( 1 )

    open(unit=10,file='stats.dat',status='old',action='write',position='append')
  
    write(unit=10,fmt="(f11.5,6(2x,e11.4))") time, prop%totEnergy, prop%potEnergy,&
                     prop%kinEnergyFluid, prop%kinEnergyWall, prop%tempFluid, prop%tempWall

    close(unit=10)

  end select



end subroutine evolutionProperties
