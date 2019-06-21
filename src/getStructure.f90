subroutine getStructure( cnt, str, prop )

  ! Purpose: Read system structure

  use global
  use random
  use gConstants
  implicit none

  type(control),    intent(inout) :: cnt
  type(structure),  intent(inout) :: str
  type(properties), intent(inout) :: prop

  integer :: i, j, k, l, m, fileOption

  real(8) :: comVel(3), temp, comMass

  character(len=2) :: nameSite


  ! Allocate memory for positions, velocities and accelerations
  allocate( str%rPart(3,str%nPart) )
  allocate( str%vPart(3,str%nPart) )
  allocate( str%fPart(3,str%nPart) )
  allocate( str%iPart(str%nPart) )
  allocate( str%mPart(str%nPart) )

  open(unit=300,file='config.inp',status='old',action='read')

  read(unit=300,fmt=*) 
  read(unit=300,fmt=*) str%lBox(:) 
  read(unit=300,fmt=*) 
  read(unit=300,fmt=*) fileOption 
  read(unit=300,fmt=*) 


  l=0

  do i=1, str%nSpecies

    do j=1, str%nMol(i)

      do k=1, str%nSites(i)

        l = l+1
      
        if( fileOption == 0 ) then
           read(unit=300,fmt=*) nameSite, str%rPart(:,l)
        else
           read(unit=300,fmt=*) nameSite, str%rPart(:,l), str%vPart(:,l)
        end if

        ! mass of each site
        str%mPart(l) = str%massSites(i)
        str%iPart(l) = i

        if( nameSite /= str%nameSites(i) ) then
          stop 'Fatal: name of sites in the CONFIG file are not consistent'
        end if


      end do
    
    end do

  end do



  close(unit=300)


   ! Setup number of degrees of fredom
  if( cnt%ensemble == 'nve_res' .or. cnt%ensemble == 'nvt_ber' .or. cnt%ensemble == 'nvt_nos') then
    str%dFreedom = 3.0d0*(dble(str%nPart) - 1.0d0)
  else if( cnt%ensemble == 'nvt_gau' ) then
    str%dFreedom = 3.0d0*(dble(str%nPart) - 1.0d0) - 1.0d0
  end if



  ! Velocities of each molecule at random
  comVel=0.0d0
  temp=0.0d0
  comMass=0.0d0

  do i=1, str%nPart

    if( fileOption == 0 ) then

      do k=1, 3
        str%vPart(k,i) = ran(seed) - 0.5d0
      end do

    end if

    comVel = comVel + str%mPart(i)*str%vPart(:,i)
    temp = temp + str%mPart(i)*dot_product(str%vPart(:,i),str%vPart(:,i))
    comMass = comMass + str%mPart(i)

  end do


  ! calculate velocity of system-center-of-mass and system temperature
  comVel = comVel/comMass
  temp = temp/str%dFreedom
  prop%temp = temp


  ! Rescale velocities to get desired temperature and zero velocity of the
  ! system center of mass.
  do i=1, str%nPart
    str%vPart(:,i) = (str%vPart(:,i) - comVel(:))*dsqrt(cnt%reqTemp/temp)
  end do




end subroutine getStructure
