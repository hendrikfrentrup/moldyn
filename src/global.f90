module global

! Module defining the principal variables involve during the simulation




type control
  ! Control variables
  integer :: nSteps              ! Total number of time steps
  integer :: nEquil              ! Equilibration period
  integer :: nSample             ! Partial averages sampling period
  integer :: nProperties         ! Frecuency of writing history of properties

  ! Control variables of history file
  integer :: nTraj(3)            ! History file options
  logical :: lTraj               ! History file?

  real(8) :: timeStep            ! Simulation time step
  real(8) :: reqTemp             ! System temperature
  character(len=7) :: ensemble   ! Statistical ensemble: nvt_res, nvt_gau, nvt_ber, nvt_nos

end type control

type structure
  ! Structural variables
  integer :: nPart                                 ! Total number of sites in the system
  integer :: nSpecies                              ! Total number of species
  integer, dimension(:),   pointer :: nMol         ! Number of molecules of each species
  integer, dimension(:),   pointer :: nSites       ! Number of sites in each molecule
  integer, dimension(:),   pointer :: iPart        ! Type of site of the sites
  real(8), dimension(:),   pointer :: massSites    ! Equilibrium distances of bond potentials
  real(8), dimension(:),   pointer :: hrmR         ! Equilibrium distances of bond potentials
  real(8), dimension(:),   pointer :: hrmK         ! Sprink constants of bond potentials
  real(8), dimension(:,:), pointer :: rPart        ! Position of the sites
  real(8), dimension(:,:), pointer :: vPart        ! Velocities of the sites
  real(8), dimension(:,:), pointer :: fPart        ! Forces of the sites
  real(8), dimension(:),   pointer :: mPart        ! masses of the sites
  real(8), dimension(3)            :: lBox         ! Simulation box length
  real(8)                          :: dFreedom     ! Degrees of freedom   
  character(len=2), dimension(:), pointer :: nameSites(:)  ! Name the sites
end type structure

type field
  ! Force field variables
  integer, dimension(:,:), pointer :: vdwMatrix    ! Interaction matrix between different species
  integer :: maxGrid                               ! Number of bins in the lookup tables
  real(8) :: deltaGrid                             ! Size of the grid for the lookup tables 
  real(8) :: rc                                    ! Value of the global (the larger) cut-off radius
  real(8), dimension(:,:), pointer :: vvv          ! Lookup table for the intermolecular potentials
  real(8), dimension(:,:), pointer :: ggg          ! Lookup table for the virial
  ! Verlet neighbour list
  integer :: nMaxList                              ! Maximum number of bins in the Verlet list
  integer, dimension(:), pointer   :: vList        ! Verlet list 
  integer, dimension(:), pointer   :: vPoint       ! Verlet pointer
  real(8) :: delr                                  ! Verlet list width
  logical :: lVerlet                               ! Verlet list option
end type field

type properties
  ! Properties variables
  real(8) :: totEnergy
  real(8) :: potEnergy
  real(8) :: kinEnergy
  real(8) :: virial
  real(8) :: temp
end type properties

type averages
  ! Average variables
  real(8) :: totEnergy(2)
  real(8) :: potEnergy(2)
  real(8) :: kinEnergy(2)
  real(8) :: temp(2)
  real(8) :: counter
end type averages


end module global

module gConstants
  integer, save :: seed = 73498403
  real(8), parameter :: pi = 3.14159265358979323846
end module gConstants
