subroutine getControl( cnt )

  ! Purpose: Get control variables from file 

  use global
  implicit none

  type(control),   intent(out) :: cnt

  integer   :: date(8)
  character :: rubbish

  open(unit=100,file='control.inp',status='old',action='read')

  read(unit=100,fmt=*) rubbish
  read(unit=100,fmt=*) rubbish, cnt%ensemble
  read(unit=100,fmt=*) rubbish, cnt%nSteps
  read(unit=100,fmt=*) rubbish, cnt%nEquil
  read(unit=100,fmt=*) rubbish, cnt%nSample
  read(unit=100,fmt=*) rubbish, cnt%nProperties
  read(unit=100,fmt=*) rubbish, cnt%timeStep
  read(unit=100,fmt=*) rubbish, cnt%reqTemp
  read(unit=100,fmt=*) rubbish, cnt%lTraj
  read(unit=100,fmt=*) rubbish, cnt%nTraj(:)

  close(unit=100)

  ! Retrieve actual time and day
  call date_and_time(values=date)

  ! Initial Print out
  open(unit=7,file='output.dat',status='unknown',action='write')

  write(unit=7,fmt="(2x,'----------------------------------------------------------',/,    &
                     2x,'    Equilibrium Molecular Dynamics for Poiseuille flow    ',/,    &
                     2x,'    Last update: December 2010                            ',/,    &
                     2x,'                                                          ',/,    &
                     2x,'    Date : ',i2,'/',i2,'/',i4,'   (d/m/y)                 ',/,    &
                     2x,'    Time : ',i2,'/',i2,'   (h/m)                          ',/,    &
                     2x,'----------------------------------------------------------',//)") &
                     date(3), date(2), date(1), date(5:6)

  write(unit=7,fmt="(3x,'--Control variables')                            ")
  write(unit=7,fmt="(/3x,'Statistical ensemble                       ',a12)")     cnt%ensemble
  write(unit=7,fmt="(3x,'Number of time steps                       ',i12)")     cnt%nSteps
  write(unit=7,fmt="(3x,'Equilibration period                       ',i12)")     cnt%nEquil
  write(unit=7,fmt="(3x,'Frecuency to report partial averages       ',i12)")     cnt%nSample
  write(unit=7,fmt="(3x,'Frecuency to write history of properties   ',i12)")     cnt%nProperties
  write(unit=7,fmt="(3x,'History file?                              ',l12)")     cnt%lTraj
  write(unit=7,fmt="(3x,'  Starts at                                ',i12)")     cnt%nTraj(1)
  write(unit=7,fmt="(3x,'  Frecuency of writting                    ',i12)")     cnt%nTraj(2)
  write(unit=7,fmt="(3x,'  File option content                      ',i12)")     cnt%nTraj(3)
  write(unit=7,fmt="(3x,'Time step                                  ',f12.4)")   cnt%timeStep
  write(unit=7,fmt="(3x,'Target temperature                         ',f12.4)")   cnt%reqTemp


  close(unit=7)




end subroutine getControl
