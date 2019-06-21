program md

! Purpose: Molecular dynamics program for LJ fluid
! History: October 2010 by H. Frentrup, C. Avendano and E.A. Muller

  use global
  implicit none

  type(control)    :: cnt
  type(structure)  :: str
  type(field)      :: fld
  type(properties) :: prop
  type(averages)   :: aver

  integer :: step

  real(8) :: time, timeNow, dispHis

  logical :: newList


  ! Read control variables
  call getControl( cnt )


  ! Read the force field
  call setup( str, fld )


  ! Read the system confguration
  call getStructure( cnt, str, prop )

  
  if( cnt%lTraj ) then
    ! Write initial configuration in the history file
    timeNow = 0.
    call writeConfig( 0, timeNow, str, cnt%nTraj(3), 0 )
  end if


  ! Set accumulators to zero
  call resetAccum( aver )


  ! Printing the header of the results in the output file
  call sample( step, timeNow, aver, 0 )


  ! Printing the header of the file for history of properties
  call evolutionProperties( 0.0, prop, 0 )


  !-----------------------!
  ! Starting main MD loop !
  !-----------------------!

  newList =.true.
  dispHis=0.0d0


  do step=1, cnt%nSteps


    ! Actual time step
    timeNow= dble(step)*cnt%timeStep


    ! Computation of intra and intermolecular forces
    call computeForces( str, fld, prop, newList )


    ! Integration of equationes of motion
    if(cnt%ensemble == 'nve_res') then
      if( step .lt. cnt%nEquil ) then
        call integrate_nvtGauss( cnt, str, prop, dispHis, step)
      else
       call integrate_nveLeapFrog( cnt, str, prop, dispHis, step)
      end if
    else if( cnt%ensemble == 'nvt_gau' ) then
      call integrate_nvtGauss( cnt, str, prop, dispHis, step)
    end if
    
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
    aver%kinEnergy(1) = aver%kinEnergy(1) + prop%kinEnergy
    aver%kinEnergy(2) = aver%kinEnergy(2) + prop%kinEnergy**2.
    aver%temp(1) = aver%temp(1) + prop%temp
    aver%temp(2) = aver%temp(2) + prop%temp**2.
    aver%counter = aver%counter + 1.

    

    ! Report partial averages and save actual configuration
    if( mod(step,cnt%nSample) == 0 ) then
      call sample( step, timeNow, aver, 1 )
      call writeConfig( step, timeNow, str, cnt%nTraj(3), 1 )
      if( step < cnt%nEquil ) call resetAccum( aver )
    end if


    ! Write history file
    if( cnt%lTraj .and. step >= cnt%nTraj(1) .and. mod(step,cnt%nTraj(2)) == 0 ) then
      call writeConfig( step, timeNow, str, cnt%nTraj(3), 0 )
    end if


    ! Write time evolution of the properties
    if( mod(step,cnt%nProperties) == 0 ) then
      call evolutionProperties( timeNow, prop, 1 )
    end if

    ! Set accumulators to zero
    if( step == cnt%nEquil ) call resetAccum( aver )

  end do


  ! Request total cpu_time
  call cpu_time(time)

  open(unit=7,file='output.dat',status='old',action='write',position='append')
  write(unit=7,fmt="(//'The total CPU time was: ',e11.4,' min',5x,'(',e11.4,' hrs)')") time/60.00, time/3600.0
  close(unit=7)

end program md
