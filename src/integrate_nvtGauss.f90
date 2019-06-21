subroutine integrate_nvtGauss( cnt, str, prop, dispHis, step )

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
  integer, intent(in) :: step


  allocate( vvv(3,str%nPart) )
  

  ! saving actual velocities
  vvv(:,:) = str%vPart(:,:)


  ! Integration of the equation of motions of the wall particles that are
  ! coupled to the Evans thermostat

  ! Calculate unconstrained velocities at new time ( t + dt/2 )
  velSq = 0.0d0

  do i=1, str%nPart

    str%vPart(:,i) = vvv(:,i) +  0.5d0*cnt%timeStep*str%fPart(:,i)/str%mPart(i)

    velSq  = velSq +  dot_product(str%vPart(:,i),str%vPart(:,i))*str%mPart(i)

  end do


  ! Unconstrained system Temperature 
  unTemp = velSq/str%dFreedom           ! check str%dFreedom


  ! Calculate the rescaling velocity factor, eta
  eta = dsqrt( cnt%reqTemp/unTemp )


  ! Constrained wall Temperature 
  prop%temp = unTemp*eta*eta


  ! Constrained Kinetic Energy
  prop%kinEnergy = 0.5d0*velSq*eta*eta


  ! Calculate new velocities (constrained) and positions using leap-frog algorithm
  velSq = 0.0d0
  speedMax=0.0d0

  do i=1, str%nPart

    ! New velocity
    str%vPart(:,i) = vvv(:,i)*(2.0d0*eta-1.0d0) + eta*cnt%timeStep*str%fPart(:,i)/str%mPart(i)
    
    ! Particles displacement
    disp(:) = cnt%timeStep*str%vPart(:,i)
    
    ! New position
    str%rPart(:,i) = str%rPart(:,i) + disp(:)
    
    velSq  = velSq +  dot_product(str%vPart(:,i),str%vPart(:,i))*str%mPart(i)

    ! Square of velocities
    speedSq = dot_product(str%vPart(:,i),str%vPart(:,i))
    speedMax = max(speedMax,speedSq)

  end do

  dispHis=dispHis+dsqrt(speedMax)*cnt%timeStep

  ! Periodic boundary conditions

  lBoxI = 1.0d0/str%lBox

  do i=1, str%nPart

    str%rPart(:,i) = str%rPart(:,i) - str%lBox(:)*nint(str%rPart(:,i)*lBoxI(:))

  end do


  ! Total system energy
  prop%totEnergy = prop%kinEnergy + prop%potEnergy



  deallocate( vvv )




end subroutine integrate_nvtGauss
