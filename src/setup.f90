subroutine setup( str, fld )

  ! Purpose: Get the force field and tabulation of the intermolecular potentials


  use global
  use gConstants
  implicit none

  type(structure), intent(out) :: str
  type(field),     intent(out) :: fld

  integer :: i, j, k, l
  integer :: nvdw

  real(8) :: lBox(3)
  real(8), dimension(:), allocatable :: vdwSigma, vdwEpsilon, vdwRc

  character(len=2), dimension(:), allocatable :: vdwIndex1, vdwIndex2
  character(len=3), dimension(:), allocatable :: vdwType




  open(unit=200,file='field.inp',status='old',action='read')

  ! Number of different species
  read(unit=200,fmt=*) str%nSpecies

  allocate( str%nMol(str%nSpecies) )
  allocate( str%nSites(str%nSpecies) )
  allocate( str%nameSites(str%nSpecies) )
  allocate( str%massSites(str%nSpecies) )
  allocate( str%hrmR(str%nSpecies) )
  allocate( str%hrmK(str%nSpecies) )

  ! Reads molecular details of the different species
  do k=1, str%nSpecies
    read(unit=200,fmt=*) 
    read(unit=200,fmt=*) str%nMol(k)
    read(unit=200,fmt=*) str%nSites(k)
    read(unit=200,fmt=*) str%nameSites(k)
    read(unit=200,fmt=*) str%massSites(k)
    read(unit=200,fmt=*) str%hrmR(k),str%hrmK(k)
  end do


  ! Calculates total number of sites
  str%nPart= sum(str%nMol*str%nSites)


  ! Read interaction parameters
  nvdw=(str%nSpecies)*(str%nSpecies-1)/2+str%nSpecies  ! number of different type of interactions
                                                     ! that should appear in the 'field' file

  allocate(vdwIndex1(nvdw))
  allocate(vdwIndex2(nvdw))
  allocate(vdwType(nvdw))
  allocate(vdwEpsilon(nvdw))
  allocate(vdwSigma(nvdw))
  allocate(vdwRc(nvdw))


  read(unit=200,fmt=*)
  do k=1,nvdw
    read(unit=200,fmt=*,end=100,err=100) vdwIndex1(k), vdwIndex2(k), vdwType(k), vdwEpsilon(k), vdwSigma(k), vdwRc(k)
  end do

  ! Read global cut-off and verlet list width and verlet option
  read(unit=200,fmt=*)
  read(unit=200,fmt=*) fld%rc
  read(unit=200,fmt=*) fld%delr
  read(unit=200,fmt=*) fld%lVerlet


  ! Checking that the global cut-off is big enough for the system
    ! Getting the simulation box dimensions
  open(unit=300,file='config.inp',status='old',action='read')
  read(unit=300,fmt=*) 
  read(unit=300,fmt=*) str%lBox(:)
  close(unit=300)


  if( fld%rc > minval(str%lBox)*0.5d0 ) then
    open(unit=7,file='output.dat',status='old',action='write',position='append')
    write(unit=7,fmt="(//3x,'Fatal: The cut-off is too large compared with the simulation box') ")
    write(unit=7,fmt="(3x,'The actual global cut-off is: ', f12.6 ) ") fld%rc
    write(unit=7,fmt="(3x,'The largest cut-off that can be used is: ', f12.6 ) ") minVal(str%lBox)*0.5d0
    write(unit=7,fmt="(3x,'The shortest box-lenth is: ', f12.6 ) ") minVal(str%lBox)
    stop
  end if


  if( fld%lVerlet ) then

    ! Scane for the simulation box lenth in the config.inp file
    open(unit=1,file='config.inp',status='old',action='read')
    read(unit=1,fmt=*)
    read(unit=1,fmt=*) lBox(:)
    close(unit=1)

    ! Maximum size of the Verlet list
    fld%nMaxList = int( 1.5*4.0*pi*(fld%delr+fld%rc)**3.0d0*dble(str%nPart**2)/6.0/product(lBox)) + 10 


    allocate(fld%vPoint(str%nPart))
    allocate(fld%vList(fld%nMaxList))
    fld%vPoint=0
    fld%vList=0

  end if



  ! Construct matrix interaction
  allocate( fld%vdwMatrix(str%nSpecies,str%nSpecies) )

  fld%vdwMatrix(:,:)=0
  do i=1, str%nSpecies

    do j=i,str%nSpecies

      do k=1, nvdw

         if( vdwIndex1(k) == str%nameSites(i) .and. vdwIndex2(k) == str%nameSites(j) &
           .or. vdwIndex1(k) == str%nameSites(j) .and. vdwIndex2(k) == str%nameSites(i) ) then
           fld%vdwMatrix(i,j)=k
           fld%vdwMatrix(j,i)=k
         end if

      end do

    end do

  end do

  ! Checking that all interactions have been specified

  do i=1, str%nSpecies
    do j=1,str%nSpecies
      if( fld%vdwMatrix(i,j) == 0 ) then
        write(unit=6,fmt=*) 'Error: Some interactions have not been specified'
        stop
      end if
    end do
  end do

  ! Look up tables for intermolecular interactions
  
  ! The resolution will be such that, using the square of the smallest value of 'sigma,' the
  ! size of the grid will be given by approximate deltaGrid=sigma^2*0.01

  fld%deltaGrid=0.005d0*minVal(vdwSigma)**2.0d0
  fld%maxGrid=int(fld%rc*fld%rc/fld%deltaGrid)
  fld%deltaGrid=fld%rc*fld%rc/dble(fld%maxGrid)




  allocate( fld%vvv(fld%maxGrid+2,nvdw) )
  allocate( fld%ggg(fld%maxGrid+2,nvdw) )
  
  do k=1, nvdw

    call vdwTables( vdwType(k), vdwEpsilon(k), vdwSigma(k), vdwRc(k), fld%vvv(:,k), &
                    fld%ggg(:,k), fld%maxGrid, fld%deltaGrid )

  end do



  close(unit=200)

  open(unit=7,file='output.dat',status='old',action='write',position='append')


  write(unit=7,fmt="(//3x,'--Particles details')                            ")
  do k=1, str%nSpecies
  write(unit=7,fmt="(/3x,'Details of species                         ',i12)")  k
  write(unit=7,fmt="(3x,'Number of molecules                        ',i12)")   str%nMol(k) 
  write(unit=7,fmt="(3x,'Number of sites                            ',i12)")   str%nSites(k) 
  write(unit=7,fmt="(3x,'Name of sites                              ',a12)")   str%nameSites(k) 
  write(unit=7,fmt="(3x,'Mass of the sites                          ',f12.4)") str%massSites(k)
  write(unit=7,fmt="(3x,'Equilibrium distance for harmonic bond     ',f12.4)") str%hrmR(k) 
  write(unit=7,fmt="(3x,'Spring constant distance for harmonic bond ',f12.4)") str%hrmK(k) 
  end do

  write(unit=7,fmt="(//3x,'-- vdw interactions ')                           ")
  write(unit=7,fmt="(3x,'  Site1',3x,'  Site2',3x,'vdwType',3x,'epsilon',3x,'  sigma')")
  do k=1,nvdw
  write(unit=7,fmt="(3x,a7,3x,a7,3x,a7,3x,f7.5,3x,f7.5)") &
       vdwIndex1(k), vdwIndex2(k), vdwType(k), vdwEpsilon(k), vdwSigma(k)
  end do

  write(unit=7,fmt="(/3x,'Global cut-off radius                      ',f12.4)") fld%rc
  write(unit=7,fmt="(3x,'Verlet list width                          ',f12.4)") fld%delr
  write(unit=7,fmt="(3x,'Verlet list ?                              ',l12  )") fld%lVerlet

  close(unit=7)

  return

  100 continue
  write(unit=6,fmt=*) 'Error: wrong number of interaction'
  stop


end subroutine setup
