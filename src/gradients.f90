subroutine gradients( str, aver, nBins, switch )

  !Purpose: Calculate the gradients (density and temperature) along the x axis

  use global
  implicit none

  type(averages), intent(inout) :: aver
  type(structure), intent(in)   :: str

  integer :: ind
  integer :: i, bin
  integer, save :: lowIndex
  integer, intent(in) :: nBins, switch
  integer, dimension(:,:), allocatable :: nPart

  real(8) :: lBox(3)
  real(8), save :: deltaBin
  real(8), dimension(:,:), allocatable :: vel, temp
  real(8), dimension(:,:), allocatable :: aux1, aux2, aux3

  character(len=8)  :: name1
  character(len=4)  :: name2
  character(len=13) :: name3

  select case (switch )

  case(0)

    ! Initialization of the gradients

    allocate( aver%tempGradient(nBins,0:str%nSpecies-1) ) 
    allocate( aver%rhoGradient(nBins,0:str%nSpecies-1) ) 
!   allocate( aver%counterGradient(nBins,0:str%nSpecies-1) ) 
    allocate( aver%counterGradient(nBins,0:2) ) 

    aver%tempGradient = 0.0d0
    aver%rhoGradient  = 0.0d0
    aver%counterGradient = 0
    aver%totalCounterGradient = 0

    deltaBin = str%lBox(1)/dble(nBins)

    lowIndex = 0
    lowIndex = sum(str%nMol(1:str%nSpecies-1) )

  case(1)

    ! Getting gradients

    allocate(temp(nBins,0:str%nSpecies-1))
    allocate(nPart(nBins,0:str%nSpecies-1))

    temp=0.0d0
    nPart=0

    ! Calculate lower index in order to avoid the wall particles
    lowIndex = sum( str%nMol(1:str%nSpecies-1) )


    do i=1, lowIndex
     
      bin = int( (str%rPart(1,i) + str%lBox(1)*0.5) / deltaBin ) + 1

      ! number of the particles in the bin per specie
      nPart(bin,str%iPart(i)) = nPart(bin,str%iPart(i)) + 1

      ! number of the particles in the bin for total particles in the fluid
      nPart(bin,0) = nPart(bin,0) + 1

    end do


    do bin=1, nBins

      do i=0, str%nSpecies-1

        if( nPart(bin,i) > 0 ) then
        
          ! gradient bin counter
          aver%counterGradient(bin,i) = aver%counterGradient(bin,i) + 1
        
        end if

      end do

    end do

    ! total counter
    aver%totalCounterGradient = aver%totalCounterGradient + 1

    do i=1, lowIndex
     
      bin = int( (str%rPart(1,i) + str%lBox(1)*0.5) / deltaBin) + 1

      if( nPart(bin,str%iPart(i)) > 0 ) then
        ! temperature
        temp(bin,str%iPart(i)) = temp(bin,str%iPart(i)) + &
                    str%mPart(i)*( str%vPart(1,i)**2.0 + str%vPart(2,i)**2.0 + str%vPart(3,i)**2.0 )
      end if

      if( nPart(bin,0) > 0 ) then
        ! temperature
        temp(bin,0) = temp(bin,0) + str%mPart(i)*( str%vPart(1,i)**2 + str%vPart(2,i)**2 + str%vPart(3,i)**2 )
      end if

    end do


    do bin=1, nBins

      do i=0, str%nSpecies-1

        if( nPart(bin,i) > 0 ) then
        
          temp(bin,i) = temp(bin,i) / dble(nPart(bin,i)) / 3.0d0
        
        end if

      end do

    end do


    do bin=1, nBins

      do i=0, str%nSpecies-1

        if( nPart(bin,i) > 0 ) then
        
          aver%rhoGradient(bin,i)  = aver%rhoGradient(bin,i)  + dble( nPart(bin,i) )
          aver%tempGradient(bin,i) = aver%tempGradient(bin,i) + temp(bin,i)
        
        end if

      end do

    end do


    deallocate(temp)
    deallocate(nPart)


  case(2)

    ! Finalization of the gradients

    allocate( aux1(nBins,0:str%nSpecies-1) )
    allocate( aux3(nBins,0:str%nSpecies-1) ) 


    do bin=1, nBins

      do i=0, str%nSpecies-1

        if( aver%counterGradient(bin,i) > 0 ) then
          aux3(bin,i) = aver%tempGradient(bin,i)/dble(aver%counterGradient(bin,i))
        end if


      end do

    end do

    aux1(:,:)  = aver%rhoGradient(:,:)/dble(aver%totalCounterGradient)

    name1='gradient'
    name2='.dat'
    ind=48

    do i=0, str%nSpecies-1

      name3 = name1//char(ind)//name2

      open(unit=8,file=name3,status='unknown',action='write')

      write(unit=8,fmt="('#   x-distance',11x,'  N',10x,'temp')")

      do bin=1, nBins
        write(unit=8,fmt="(3(2x,f12.6))") (dble(bin)*deltaBin-str%lBox(1)*0.5d0), aux1(bin,i), aux3(bin,i)
      end do

      close(8)

      ind = ind + 1

    end do


    deallocate( aux1 )
    deallocate( aux3 )


  end select

end subroutine gradients
