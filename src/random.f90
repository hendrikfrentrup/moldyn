module random


  contains

  function ran( idum )


  ! Purpose: Generator of random numbers “Minimal” random number generator of Park
  ! and Miller 
  !          combined with a Marsaglia shift sequence. Returns a uniform random
  !          deviate between 
  !          0.0 and 1.0 (exclusive of the endpoint values). This fully portable,
  !          scalar generator 
  !          has the “traditional” (not Fortran 90) calling sequence with a random
  !          deviate as the 
  !          returned function value: call with idum a negative integer to
  !          initialize; thereafter, 
  !          do not alter idum except to reinitialize. The period of this
  !          generator is about 3.1 × 1018 .

  implicit none

  integer, parameter :: k4b = selected_int_kind(9)

  integer (k4b), save :: ix = -1
  integer (k4b), save :: iy = -1
  integer (k4b), save :: k

  integer (k4b), parameter :: ia = 16807
  integer (k4b), parameter :: im = 2147483647
  integer (k4b), parameter :: iq = 127773
  integer (k4b), parameter :: ir = 2836

  integer (k4b), intent(inout) :: idum

  real(8) :: ran

  real(8), save :: am


  ! ... initialize.
  if ( idum <= 0 .or. iy < 0 ) then
      am = nearest( 1.0, -1.0 )/im                 
      iy = ior(ieor ( 888889999, abs( idum )),1)
      ix = ieor(777755555,abs( idum ))
      ! Set idum positive.
      idum = abs( idum ) + 1    
  end if

  ! ... marsaglia shift sequence with period 232 − 1.
  ix = ieor(ix, ishft(ix,13) )
  ix = ieor(ix, ishft(ix,-17) )
  ix = ieor(ix, ishft(ix,5) )

  ! ...Park-Miller sequence by Schrage’s method, period 231 − 2.
  k = iy/iq 
  iy = ia*(iy - k*iq) - ir*k

  if (iy < 0) iy = iy + im

  ! Combine the two generator swith masking to ensure nonzero value.
  ran = am*ior( iand (im, ieor( ix, iy )),1 )      

  return
  end function ran

end module random
