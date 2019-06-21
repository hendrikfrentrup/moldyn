subroutine resetAccum( aver )

  use global
  implicit none

  type(averages), intent(out) :: aver

  aver%totEnergy(:) = 0.0d0 
  aver%potEnergy(:) = 0.0d0
  aver%kinEnergy(:) = 0.0d0
  aver%temp  = 0.0d0
  aver%counter = 0.0d0

end subroutine resetAccum
