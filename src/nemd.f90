program nemd

! Purpose: Non-equilibrium Molecular dynamics program for Poiseuille flow
! History: October 2010 by H. Frentrup, C. Avendano and E. A. Muller

  use global
  use gConstants
  implicit none

  type(control)    :: cnt
  type(structure)  :: str
  type(field)      :: fld
  type(properties) :: prop
  type(averages)   :: aver
  type(fluxes)     :: flx

  integer :: i, step, lowIndex

  real(8), dimension(:), pointer :: signOld
  real(8), dimension(:), pointer :: signNew
  real(8) :: time, timeNow, dispHis

  logical :: newList


  ! Get a new seed from system clock
  call system_clock( seed )

  ! Read control variables
  call getControl( cnt )


  ! Read the force field
  call setup( str, fld )


  ! Read the system confguration
  call getStructure( cnt, str, prop )

  
  if( cnt%lTraj .and. cnt%nTraj(1) == 0 ) then
    ! Write initial configuration in the history file
    timeNow = 0.
    call writeConfig( 0, timeNow, str, cnt%nTraj(1:3), 0 )
  end if


  ! Printing the header of the results in the output file
  call sample( step, timeNow, aver, flx, cnt%timeStep, 0 )


  ! Printing the header of the file for history of properties
  call evolutionProperties( 0.0, prop, 0 )

  ! Initializing density, velocity and temperature profiles/gradients
  if( cnt%lProfile ) then
    call gradients( str, aver, cnt%nBinsProfile, 0 )
    call profiles( str, aver, cnt%nBinsProfile, 0 )
  end if


  ! Calculate lower index in order to avoid the wall in the calculation of
  ! fluxes
  lowIndex = sum( str%nMol(1:str%nSpecies-1) )

  ! allocate and set to zero accumulators for fluxes
  allocate( flx%leftFlux(str%nSpecies-1) )
  allocate( flx%rightFlux(str%nSpecies-1) )
  flx%counter = 0.0d0
  flx%leftFlux(:)  = 0.0d0
  flx%rightFlux(:) = 0.0d0

  ! auxiliar vector to decide if a particle has crossed the plane x=0
  allocate( signOld(lowIndex) )
  allocate( signNew(lowIndex) )


  ! Set accumulators to zero
  call resetAccum( aver, flx, 0 )


  !--------------------------------------------!
  ! Starting main MD loop                      !
  !--------------------------------------------!

  newList =.true.
  dispHis=0.0d0

  do step=1, cnt%nSteps


    ! Actual time step
    timeNow= dble(step)*cnt%timeStep


    ! calculate if particles are at right or left hand side of the box before
    ! integrating the equations of motion. This will be use for the flux calculation
    signOld(:) = 1.0
    signNew(:) = 1.0
    signOld(:) = sign( signOld(:), str%rPart(1,1:lowIndex) )



    ! Computation of intra and intermolecular forces
    call computeForces( str, fld, prop, newList )


    ! Integration of equationes of motion
    call integrate_nvtGauss( cnt, str, prop, dispHis )


    ! calculate if particles are at right or left hand side of the box after
    ! integrating the equations of motion. This will be use for the flux calculation
      signNew(:) = sign( signNew(:), str%rPart(1,1:lowIndex) )
      flx%counter = flx%counter + 1.0d0

      do i=1, lowIndex
        ! check net move from left to right or right to left
        if( abs(str%rPart(1,i) ) < 2.0d0 ) then
          if( signOld(i) < 0.0d0 .and. signNew(i) > 0.0d0 ) then
            flx%rightFlux( str%iPart(i) ) = flx%rightFlux( str%iPart(i) ) + 1.0d0
          else if( signOld(i) > 0.0d0 .and. signNew(i) < 0.0d0 ) then
            flx%leftFlux( str%iPart(i) ) = flx%leftFlux( str%iPart(i) ) + 1.0d0
          end if
        end if
      end do


    
    ! Checking the Verl list
    if( dispHis > 0.5d0*fld%delr) then
      dispHis=0.0
      newList=.true.
    end if

    ! Update statistics
    aver%totEnergy(1) = aver%totEnergy(1) + prop%totEnergy 
    aver%totEnergy(2) = aver%totEnergy(2) + prop%totEnergy**2. 
    aver%potEnergy(1) = aver%potEnergy(1) + prop%potEnergy
    aver%potEnergy(2) = aver%potEnergy(2) + prop%potEnergy**2.
    aver%kinEnergyFluid(1) = aver%kinEnergyFluid(1) + prop%kinEnergyFluid
    aver%kinEnergyFluid(2) = aver%kinEnergyFluid(2) + prop%kinEnergyFluid**2.
    aver%kinEnergyWall(1) = aver%kinEnergyWall(1) + prop%kinEnergyWall
    aver%kinEnergyWall(2) = aver%kinEnergyWall(2) + prop%kinEnergyWall**2.
    aver%tempFluid(1) = aver%tempFluid(1) + prop%tempFluid
    aver%tempFluid(2) = aver%tempFluid(2) + prop%tempFluid**2.
    aver%tempWall(1)  = aver%tempWall(1)  + prop%tempWall
    aver%tempWall(2)  = aver%tempWall(2)  + prop%tempWall**2.
    aver%counter = aver%counter + 1.

    
    ! Caclulate velocity profile
    if( mod(step,cnt%nProfile)==0 .and. cnt%lProfile .and. step > cnt%nEquil ) then
      call gradients( str, aver, cnt%nBinsProfile, 1 ) 
      call profiles( str, aver, cnt%nBinsProfile, 1 ) 
    end if


    ! Report partial averages and save actual configuration
    if( mod(step,cnt%nSample) == 0 ) then
      call sample( step, timeNow, aver, flx, cnt%timeStep, 1 )
      call writeConfig( step, timeNow, str, cnt%nTraj(:), 1 )
      ! Finish and print profiles and gradients
      if(  cnt%lProfile .and. step > cnt%nEquil ) then
        call gradients( str, aver, cnt%nBinsProfile, 2 )
        call profiles( str, aver, cnt%nBinsProfile, 2 )  
      end if
    end if


    ! Write history file
    if( cnt%lTraj .and. step >= cnt%nTraj(1) .and. mod(step,cnt%nTraj(2)) == 0 ) then
      call writeConfig( step, timeNow, str, cnt%nTraj(:), 0 )
    end if


    ! Write time evolution of the properties
    if( mod(step,cnt%nProperties) == 0 ) then
      call evolutionProperties( timeNow, prop, 1 )
    end if

    ! End of equilibration - Set accumulators to zero
    if( step == cnt%nEquil ) then
      call resetAccum( aver, flx, 1 )
      ! Printing the header of the results in the output file
      call sample( step, timeNow, aver, flx, cnt%timeStep, 0 )
    end if

  end do

  ! Finish and print profiles gradients
  if(  cnt%lProfile ) then
    call gradients( str, aver, cnt%nBinsProfile, 2 ) 
    call profiles( str, aver, cnt%nBinsProfile, 2 ) 
  end if

  open(unit=7,file='output.dat',status='old',action='write',position='append')

  ! report final value for flux

  write(unit=7,fmt="(//a,f15.8)") 'Total Flux*Area of species 1', &
  -( flx%leftFlux(1) - flx%rightFlux(1) )/cnt%timeStep/flx%counter
!  write(unit=7,fmt="(  a,f15.8)") 'Total Flux*Area of species 2', &
!  -( flx%leftFlux(2) - flx%rightFlux(2) )/cnt%timeStep/flx%counter

  write(unit=7,fmt="(//a,f15.8)") 'Rel. residence time', dble(aver%poreAdsor)/dble(aver%poreTotal)
  write(unit=7,fmt="(//a,f15.8)") 'Ad-/Desorption events per total T_pore', dble(aver%Adsorptions)/dble(aver%poreTotal)
  write(unit=7,fmt="(//a,f15.8)") 'Ad-/Desorption events per T_ads', dble(aver%Adsorptions)/dble(aver%poreAdsor)

  ! Request total cpu_time
  call cpu_time(time)


  write(unit=7,fmt="(//'The total CPU time was: ',e11.4,' min',5x,'(',e11.4,' hrs)')") time/60.00, time/3600.0

  close(unit=7)

end program nemd
