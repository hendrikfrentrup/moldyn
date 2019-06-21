subroutine integrate_nvtGauss( cnt, str, prop, dispHis )

  ! Purpose: Integration of equations of motions using the Leap-Frog integrator.
  !          The temperature (kinetic energy) is fixed using the Evans
  !          thermostat. Note that the thermostat only applies to the molecules
  !          in the wall

  use global
  implicit none

  type(control),    intent(in) :: cnt
  type(structure),  intent(inout) :: str
  type(properties), intent(inout) :: prop

  integer :: i, j, k, l, lowIndex

  real(8) :: velSq, speedSq, speedMax, eta, unTemp, disp(3), lBoxI(3)
  real(8), dimension(:,:), allocatable :: vvv
  real(8), intent(inout) :: dispHis


  allocate( vvv(3,str%nPart) )
  

  ! saving actual velocities
  vvv(:,:) = str%vPart(:,:)


  ! Integration of the equation of motions of the fluid particles

  ! Calculate lower index in order to avoid the wall particles, which are
  ! integrate in a different way
  lowIndex = sum(str%nMol(1:str%nSpecies-1) )


  outer: do i=1, lowIndex

    ! New velocities at time 't+delta/2'
    
    str%vPart(:,i) = vvv(:,i) + cnt%timeStep*str%fPart(:,i)/str%mPart(i)
    
    ! Particles displacement
    
    disp(:) = cnt%timeStep*str%vPart(:,i)
    
    ! New position of particles at time 't+delta'
    
    str%rPart(:,i) = str%rPart(:,i) + disp(:)

  end do outer


  ! Calculate temperature an kinetic energy of the fluid at full time step and
  ! accumulating maximum speed to check the updating of the Verlet list

  velSq  = 0.0d0
  speedMax=0.0d0

  do i=1, lowIndex

    ! velocity at full time step
    vvv(:,i) = 0.5d0*( vvv(:,i) + str%vPart(:,i) )
    
    ! Square of velocities
    speedSq = dot_product(vvv(:,i),vvv(:,i))

    velSq = velSq + str%mPart(i)*speedSq

    speedMax = max(speedMax,speedSq)
     
  end do

  dispHis=dispHis+dsqrt(speedMax)*cnt%timeStep


  ! Kinetic Energy
  prop%kinEnergyFluid = 0.5d0*velSq


  ! System Temperature
  prop%tempFluid = velSq/(3.0d0*dble(lowIndex)-3.0d0)






  ! Integration of the equation of motions of the wall particles that are
  ! coupled to the Evans thermostat

  ! Calculate unconstrained velocities at new time t 
  velSq = 0.0d0

  do i=lowIndex+1, str%nPart

    str%vPart(:,i) = vvv(:,i) +  0.5d0*cnt%timeStep*str%fPart(:,i)/str%mPart(i)

    velSq  = velSq +  dot_product(str%vPart(:,i),str%vPart(:,i))*str%mPart(i)

  end do

  
  ! Unconstrained system Temperature 
  unTemp = velSq/(3.0d0*dble(str%nMol(str%nSpecies))-1)


  ! Calculate the rescaling velocity factor, eta
  eta = dsqrt( cnt%reqTemp/unTemp )


  ! Constrained wall Temperature 
  prop%tempWall = unTemp*eta*eta


  ! Constrained Kinetic Energy
  prop%kinEnergyWall = 0.5d0*velSq*eta*eta



  ! Calculate new velocities (constrained) and positions using leap-frog algorithm

  do i=lowIndex+1, str%nPart

    ! New velocity
    str%vPart(:,i) = vvv(:,i)*(2.0d0*eta-1.0d0) + eta*cnt%timeStep*str%fPart(:,i)/str%mPart(i)
    
    ! Particles displacement
    disp(:) = cnt%timeStep*str%vPart(:,i)
    
    ! New position
    str%rPart(:,i) = str%rPart(:,i) + disp(:)

  end do
      
      
               



  ! Periodic boundary conditions

  lBoxI = 1.0d0/str%lBox

  do i=1, str%nPart

    str%rPart(:,i) = str%rPart(:,i) - str%lBox(:)*nint(str%rPart(:,i)*lBoxI(:))

  end do


  ! Total system energy
  prop%totEnergy = prop%kinEnergyWall + prop%kinEnergyFluid + prop%potEnergy



  deallocate( vvv )




end subroutine integrate_nvtGauss
