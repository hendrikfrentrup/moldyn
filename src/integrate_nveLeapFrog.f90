subroutine integrate_nveLeapFrog( cnt, str, prop, dispHis, step )

  ! Purpose: Integration of equations of motions using the Leap-Frog integrator.

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


  ! Integration of the equation of motions using the leap-frog algorithm

  do i=1, str%nPart
    
    ! New velocities at time 't+delta/2'
    
    str%vPart(:,i) = vvv(:,i) + cnt%timeStep*str%fPart(:,i)/str%mPart(i)
    
    ! Particles displacement
    
    disp(:) = cnt%timeStep*str%vPart(:,i)
    
    ! New position of particles at time 't+delta'
    
    str%rPart(:,i) = str%rPart(:,i) + disp(:)

  end do

  ! Calculate temperature and kinetic energy of the fluid at full time step and
  ! accumulating maximum speed to check the updating of the Verlet list

  velSq  = 0.0d0
  speedMax=0.0d0

  do i=1, str%nPart

    ! velocity at full time step
    vvv(:,i) = 0.5d0*( vvv(:,i) + str%vPart(:,i) )
    
    ! Square of velocities
    speedSq = dot_product(vvv(:,i),vvv(:,i))

    velSq = velSq + str%mPart(i)*speedSq

    speedMax = max(speedMax,speedSq)
     
  end do

  dispHis=dispHis+dsqrt(speedMax)*cnt%timeStep


  ! Kinetic Energy
  prop%kinEnergy = 0.5d0*velSq


  ! System Temperature
  prop%temp = velSq/(3.0d0*dble(str%nPart)-3.0d0)



  ! Periodic boundary conditions

  lBoxI = 1.0d0/str%lBox

  do i=1, str%nPart

    str%rPart(:,i) = str%rPart(:,i) - str%lBox(:)*nint(str%rPart(:,i)*lBoxI(:))

  end do


  ! Total system energy
  prop%totEnergy = prop%kinEnergy + prop%potEnergy



  deallocate( vvv )




end subroutine integrate_nveLeapFrog
