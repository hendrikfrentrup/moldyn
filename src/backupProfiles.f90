subroutine profiles( str, aver, nBins, switch )

  use global
  implicit none

  type(averages), intent(inout) :: aver
  type(structure), intent(in)   :: str

  integer :: i, bin
  integer, save :: lowIndex
  integer, intent(in) :: nBins, switch

  real(8) :: lBox(3)
  real(8), save :: deltaBin
  real(8), dimension(:), allocatable :: vel, temp, nPart


  select case (switch )

  case(0)

    ! Initialization of the profiles

    allocate( aver%velProfile(nBins) ) 
    allocate( aver%tempProfile(nBins) ) 
    allocate( aver%rhoProfile(nBins) ) 

    aver%velProfile  = 0.0d0
    aver%tempProfile = 0.0d0
    aver%rhoProfile  = 0.0d0
    aver%counterProfile = 0.0d0

    deltaBin = str%lBox(1)/dble(nBins)

    lowIndex = 0
    lowIndex = sum(str%nMol(1:str%nSpecies-1) )

  case(1)

    ! Getting profiles

    allocate(vel(nBins))
    allocate(temp(nBins))
    allocate(nPart(nBins))

    vel=0.0d0
    temp=0.0d0
    nPart=0.0d0

    ! Calculate lower index in order to avoid the wall particles
    lowIndex = sum( str%nMol(1:str%nSpecies-1) )

    aver%counterProfile = aver%counterProfile + 1.0

    do i=1, lowIndex
     
      bin = int( (str%rPart(1,i) + str%lBox(1)*0.5) / deltaBin ) + 1

      ! stream velocity of the bin
      vel(bin) = vel(bin) + str%vPart(3,i)

      ! number of the particles in the bin
      nPart(bin) = nPart(bin) + 1.0

    end do

    where(nPart/=0) vel(:)=vel(:)/nPart(:)


    do i=1, lowIndex
     
      bin = int( (str%rPart(1,i) + str%lBox(1)*0.5) / deltaBin) + 1

      ! temperature
      temp(bin) = temp(bin) + &
                  str%mPart(i)*( str%vPart(1,i)**2.0 + str%vPart(2,i)**2.0 + ( str%vPart(3,i) - vel(bin) )**2.d0)

    end do

    where(nPart/=0) temp(:) = temp(:)/nPart(:)/3.0d0


    aver%rhoProfile(:)  = aver%rhoProfile(:)  + nPart(:)
    aver%velProfile(:)  = aver%velProfile(:)  + vel(:)
    aver%tempProfile(:) = aver%tempProfile(:) + temp(:)


    deallocate(vel)
    deallocate(temp)
    deallocate(nPart)


  case(2)

    ! Finalization of the profiles
    aver%velProfile(:)  = aver%velProfile(:)/aver%counterProfile
    aver%rhoProfile(:)  = aver%rhoProfile(:)/aver%counterProfile/str%lBox(2)/str%lBox(3)/deltaBin
    aver%tempProfile(:) = aver%tempProfile(:)/aver%counterProfile

    open(unit=8,file='profiles.dat',status='unknown',action='write')
    write(unit=8,fmt="('#   x-distance',11x,'rho',11x,'vel',10x,'temp')")
    do bin=1, nBins
    write(unit=8,fmt="(4(2x,f12.6))") (dble(bin)*deltaBin-str%lBox(1)*0.5d0), aver%rhoProfile(bin), &
                                                                              aver%velProfile(bin), &
                                                                              aver%tempProfile(bin)
    end do


  end select

end subroutine profiles
