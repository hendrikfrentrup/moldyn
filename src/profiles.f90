subroutine profiles( str, aver, nBins, switch )

  !Purpose: Calculate profiles (velocity, density and temperature) within the pore

  use global
  implicit none

  type(averages), intent(inout) :: aver
  type(structure), intent(inout)   :: str

  integer :: ind
  integer :: i, bin
  integer, save :: lowIndex
  integer, intent(in) :: nBins, switch
  integer, dimension(:,:), allocatable :: nPart

  real(8), save :: deltaBin
  real(8), dimension(:,:), allocatable :: vel, avVel, temp
  real(8), dimension(:,:), allocatable :: temp1, temp2, temp3

  integer, dimension(:), allocatable :: adsEvent

  character(len=7)  :: name1
  character(len=4)  :: name2
  character(len=13) :: name3

  select case (switch)

  case(0)

    ! Initialization of the profiles

    allocate( aver%velProfile(nBins,0:str%nSpecies-1) ) 
    allocate( aver%tempProfile(nBins,0:str%nSpecies-1) ) 
    allocate( aver%rhoProfile(nBins,0:str%nSpecies-1) )
!   NOTE: counterProfile currently only goes to up to 2 species right now!!!
    allocate( aver%counterProfile(nBins,0:str%nSpecies-1) ) 
!    allocate( aver%counterProfile(nBins,0:2) ) 

    aver%velProfile  = 0.0d0
    aver%tempProfile = 0.0d0
    aver%rhoProfile  = 0.0d0
    aver%counterProfile = 0
    aver%totalCounterProfile = 0

    aver%poreTotal   = 0
    aver%poreAdsor   = 0
    aver%Adsorptions = 0

    deltaBin = str%lBox(3)/dble(nBins)
    lowIndex = 0
    lowIndex = sum(str%nMol(1:str%nSpecies-1) )

    ! Initialization of adsorption status
    allocate( str%poreStat(lowIndex,2) )
    allocate( str%adsStat(lowIndex,2) ) 
    str%poreStat=0
    str%adsStat=0

  case(1)

    ! Getting profiles

    allocate( vel(nBins,0:str%nSpecies-1) )
    allocate( avVel(nBins,0:str%nSpecies-1) )
    allocate( temp(nBins,0:str%nSpecies-1) )
    allocate( nPart(nBins,0:str%nSpecies-1) )

    allocate( adsEvent(lowIndex) )
    str%poreStat(:,2) = str%poreStat(:,1)
    str%poreStat(:,1) = 0
    str%adsStat(:,2)  = str%adsStat(:,1)
    str%adsStat(:,1) = 0
    adsEvent=0

    vel=0.0d0
    avVel=0.0d0
    temp=0.0d0
    nPart=0

    ! Calculate lower index in order to avoid the wall particles
    lowIndex = sum( str%nMol(1:str%nSpecies-1) )


    do i=1, lowIndex

      if ( abs(str%rPart(1,i)) < str%lPore(1) ) then

        str%poreStat(i,1)=1.0d0

        if ( abs(str%rPart(3,i)) > 0.5d0*(str%lPore(3)-1) ) then
          str%adsStat(i,1)=1.0d0
        end if

      	bin = int( (str%rPart(3,i) + str%lBox(3)*0.5) / deltaBin ) + 1

      	! number of the particles in the bin per specie
      	nPart(bin,str%iPart(i)) = nPart(bin,str%iPart(i)) + 1

      	! number of the particles in the bin for total fluid
      	nPart(bin,0) = nPart(bin,0) + 1
      
      	! stream velocity of the bin per specie
      	vel(bin,str%iPart(i)) = vel(bin,str%iPart(i)) + str%vPart(1,i)
     
      	! stream velocity of the bin for total fluid
      	vel(bin,0) = vel(bin,0) + str%vPart(1,i)

      end if

    end do

    do i=0, str%nSpecies-1

      do bin=1, nBins

        if ( nPart(bin,i) > 0 ) then
        
          ! profile bin counter
          aver%counterProfile(bin,i) = aver%counterProfile(bin,i) + 1

          ! average streaming velocity
          avVel(bin,i) = aver%velProfile(bin,i) + vel(bin,i) / dble(nPart(bin,i))

        end if

      end do

    end do

    ! total counter
    aver%totalCounterProfile = aver%totalCounterProfile + 1

    do i=1, lowIndex

      	bin = int( (str%rPart(3,i) + str%lBox(3)*0.5) / deltaBin) + 1

      	if ( nPart(bin,0) > 0 .and. abs(str%rPart(1,i)) < str%lPore(1) ) then
          ! temperature of total fluid
          temp(bin,0) = temp(bin,0) + str%mPart(i)*( str%vPart(1,i)**2 + str%vPart(2,i)**2 + str%vPart(3,i)**2 )

          if ( nPart(bin,str%iPart(i)) > 0 ) then
            ! temperature of each species
            temp(bin,str%iPart(i)) = temp(bin,str%iPart(i)) + &
                                     str%mPart(i)*( str%vPart(1,i)**2.0 + str%vPart(2,i)**2.0 + str%vPart(3,i)**2.0 )
          end if

        end if

    end do

    ! Divide velocity and temperature through number of particles in each bin
    ! AND update averages

    do i=0, str%nSpecies-1

      do bin=1, nBins

        if( nPart(bin,i) > 0 ) then

          vel(bin,i)  = vel(bin,i) / dble(nPart(bin,i))
          temp(bin,i) = temp(bin,i) / dble(nPart(bin,i)) / 3.0d0

          aver%velProfile(bin,i)  = aver%velProfile(bin,i)  + vel(bin,i)        
          aver%rhoProfile(bin,i)  = aver%rhoProfile(bin,i)  + dble( nPart(bin,i) )
          aver%tempProfile(bin,i) = aver%tempProfile(bin,i) + temp(bin,i)
        
        end if

      end do

    end do

    ! Calculate the number of adsorption/desorption events

    ! Calculate total number of particles in pore and number of adsorbed particles
    where( str%poreStat(:,1) .eq. str%poreStat(:,2) .and. str%adsStat(:,1) .ne. str%adsStat(:,2) ) adsEvent(:) = 1
    aver%poreTotal = aver%poreTotal + sum(str%poreStat(:,1))
    aver%poreAdsor = aver%poreAdsor + sum(str%adsStat(:,1))
    aver%Adsorptions = aver%Adsorptions + sum(adsEvent)

!    write(*,*) aver%poreAdsor, ' / ', aver%poreTotal, ' = ', dble(aver%poreAdsor)/dble(aver%poreTotal)

    deallocate(vel)
    deallocate(avVel)
    deallocate(temp)
    deallocate(nPart)

    deallocate(adsEvent)

  case(2)

    ! Finalization of the profiles

    allocate( temp1(nBins,0:str%nSpecies-1) ) ! used for density profile
    allocate( temp2(nBins,0:str%nSpecies-1) ) ! used for velocity profile
    allocate( temp3(nBins,0:str%nSpecies-1) ) ! used for temperature profile

    temp1=0.d0
    temp2=0.d0
    temp3=0.d0

    do i=0, str%nSpecies-1

      where(aver%counterProfile(:,i) .ne. 0) temp2(:,i) = aver%velProfile(:,i)/dble(aver%counterProfile(:,i))
      where(aver%counterProfile(:,i) .ne. 0) temp3(:,i) = aver%tempProfile(:,i)/dble(aver%counterProfile(:,i))

    end do

    temp1(:,:)  = aver%rhoProfile(:,:)/dble(aver%totalCounterProfile)

    name1='profile'
    name2='.dat'
    ind=48

    do i=0, str%nSpecies-1

      name3 = name1//char(ind)//name2

      open(unit=8,file=name3,status='unknown',action='write')

      write(unit=8,fmt="('#   z-distance',11x,'  N',10x,'vel',10x,'temp')")

      do bin=1, nBins

        write(unit=8,fmt="(4(2x,f12.6))") (dble(bin)*deltaBin-str%lBox(3)*0.5d0), temp1(bin,i), temp2(bin,i), temp3(bin,i)

      end do

      close(8)

      ind = ind + 1

    end do


    deallocate( temp1 )
    deallocate( temp2 )
    deallocate( temp3 )


  end select

end subroutine profiles
