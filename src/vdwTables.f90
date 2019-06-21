subroutine vdwTables( vdwType, vdwEpsilon, vdwSigma, vdwRc, vvv, ggg, maxGrid, deltaGrid )

  implicit none

  integer ::  i
  integer, intent(in) ::  maxGrid

  real(8) :: rij, s
  real(8), intent(in) :: vdwEpsilon, vdwSigma, vdwRc, deltaGrid
  real(8), dimension(maxGrid) :: vvv, ggg

  character(len=3) :: vdwType


  vvv = 0.0d0
  ggg = 0.0d0

  do i=1, maxGrid

    s=dble(i)*deltaGrid
    rij=dsqrt(s)


    if( vdwType == 'lj ' ) then

      if(rij > 0.10*vdwSigma) then

        if( rij <= vdwRc ) then
          vvv(i) = 4.0d0*vdwEpsilon*((vdwSigma/rij)**12.0d0-(vdwSigma/rij)**6.0d0)-&
                   4.0d0*vdwEpsilon*((vdwSigma/vdwRc)**12.0d0-(vdwSigma/vdwRc)**6.0d0)
          ggg(i) = 48.0d0*vdwEpsilon/vdwSigma*((vdwSigma/rij)**12.0d0-0.5d0*(vdwSigma/rij)**6.0d0)/rij**2.0d0
                   
        end if

      else

        if( rij <= vdwRc ) then
          vvv(i) = 4.0d0*vdwEpsilon*((1.0d0/0.10d0)**12.0d0-(1.0d0/0.10d0)**6.0d0)-&
                   4.0d0*vdwEpsilon*((vdwSigma/vdwRc)**12.0d0-(vdwSigma/vdwRc)**6.0d0)
          ggg(i) = 48.0d0*vdwEpsilon/vdwSigma*((1.0d0/0.10d0)**12.0d0-0.5d0*(1.0d0/0.10d0)**6.0d0)&
                   /(0.10d0*vdwSigma)**2.0d0
                   
        end if

      end if

    else if( vdwType == 'wca' ) then
      
      if(rij > 0.10*vdwSigma) then

        if( rij <= vdwSigma*2.0d0**(1.0d0/6.0d0) ) then
          vvv(i) = 4.0d0*vdwEpsilon*((vdwSigma/rij)**12.0d0-(vdwSigma/rij)**6.0d0)+vdwEpsilon
          ggg(i) = 48.0d0*vdwEpsilon/vdwSigma*((vdwSigma/rij)**12.0d0-0.5d0*(vdwSigma/rij)**6.0d0)/rij**2.0d0
        end if

      else

        if( rij <= vdwSigma*2.0d0**(1.0d0/6.0d0) ) then
          vvv(i) = 4.0d0*vdwEpsilon*((1.0d0/0.10d0)**12.0d0-(1.0d0/0.10d0)**6.0d0)+vdwEpsilon
          ggg(i) = 48.0d0*vdwEpsilon/vdwSigma*((1.0d0/0.10d0)**12.0d0-0.5d0*(1.0d0/0.10d0)**6.0d0)&
                   /(vdwSigma*0.10d0)**2.0d0
        end if

      end if

    end if

  end do



end subroutine vdwTables
