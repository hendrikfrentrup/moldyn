subroutine evolutionProperties( time, prop, switch )

  use global
  implicit none

  type(properties), intent(in) :: prop

  integer, intent(in) :: switch

  real(8), intent(in) :: time


  select case( switch ) 

  case( 0 )


    open(unit=10,file='stats.dat',status='unknown',action='write')

    write(unit=10,fmt="('#      time',2x,'  totEnergy',2x,'  potEnergy',2x,&
                       '   kinEnery',2x,'temperature')") 

    close(unit=10)

  case( 1 )

    open(unit=10,file='stats.dat',status='old',action='write',position='append')
  
    write(unit=10,fmt="(f11.5,4(2x,e11.4))") time, prop%totEnergy, prop%potEnergy,&
                     prop%kinEnergy, prop%temp

    close(unit=10)

  end select



end subroutine evolutionProperties
