Subroutine writeConfig( actualStep, timeNow, str, trajOption, switch )

  ! Purpose: Write configurations or history or config file

  use global   
  implicit none

  type(structure), intent(in) :: str

  integer :: i, j, k, l
  integer, intent(in) :: actualStep, switch, trajOption(3)

  real(8), intent(in) :: timeNow



  select case(switch)

  case(0)
    ! ----------------------------------------------------------------------
    ! write history file

    if( actualStep == trajOption(1) ) then

      open(unit=1,file='traj.xyz',status='unknown',action='write')

    else

      open(unit=1,file='traj.xyz',status='old',position='append',action='write')

    end if

    write(unit=1,fmt="(i12)") str%nPart
    write(unit=1,fmt="(i12,4(2x,f12.6))") actualStep, timeNow, str%lBox  

    i = 0
    do j=1,str%nSpecies
      do k=1, str%nMol(j)
        do l=1, str%nSites(j)

          i = i + 1
          if( trajOption(3) == 0 ) then
            write(unit=1,fmt="(a3,3(1x,f15.8))")  str%nameSites(j), str%rPart(:,i)
          else if( trajOption(3) == 1 ) then
            write(unit=1,fmt="(a3,6(1x,f15.8))")  str%nameSites(j), str%rPart(:,i), str%vPart(:,i)
          else
            write(unit=1,fmt="(a3,9(1x,f15.8))")  str%nameSites(j), str%rPart(:,i), str%vPart(:,i), str%fPart(:,i)
          end if

        end do
      end do
    end do


    close(unit=1)


  case(1)

    ! ----------------------------------------------------------------------
    ! write config file

    open(unit=1,file='config.new',status='unknown',action='write')

    write(unit=1,fmt="(' Simulation box length')")
    write(unit=1,fmt="(3(2x,f12.6))") str%lBox
    write(unit=1,fmt="(' File option: 0-->coordinates, 1--> coordinates and velocities')")
    write(unit=1,fmt="(2x,i5)") 1
    write(unit=1,fmt="(' Coordinates and velocties of each site')")

    i = 0
    do j=1,str%nSpecies
      do k=1, str%nMol(j)
        do l=1, str%nSites(j)

          i = i + 1
          write(unit=1,fmt="(2x,a2,6(3x,f12.6))")  str%nameSites(j), str%rPart(:,i), str%vPart(:,i)

        end do
      end do
    end do

    write(unit=1,fmt="(' Lattice sites of particles the wall')")
    do i=1, str%nMol(str%nSpecies)
          write(unit=1,fmt="(4x,3(3x,f12.6))")  str%wPart(:,i)
    end do

    close(unit=1)

  end select

end subroutine writeConfig
