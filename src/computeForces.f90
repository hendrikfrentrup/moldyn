subroutine computeForces( str, fld, prop, newList )

  ! Purpose: Calculation of inter and intra molecular forces

  use global
  implicit none

  type(structure),  intent(inout) :: str
  type(field),      intent(inout)    :: fld 
  type(properties), intent(out)   :: prop

  integer :: i, j, jj, k, a, b, ab, nWall, nList
  integer :: index1, index2


  real(8) :: sds, xi, vk0, vk1, vk2, t1, t2, uij, wij, lx
  real(8) :: lBoxI(3), rij(3), rr, fVal, rcSq, rListSq

  logical, intent(inout) :: newList


  str%fPart(:,:) = 0.0d0
  prop%potEnergy = 0.0d0
  prop%virial = 0.0d0
  lBoxI(:) = 1.0d0/str%lBox(:)
  rcSq = fld%rc*fld%rc

  
  If( fld%lVerlet .and. newList ) then

    ! Non-bonded interactions
    rListSq = (fld%rc+fld%delr)**2.0d0
    nList=0
  
    do i=1, str%nPart-1

      fld%vPoint(i)=nList+1
  
      a = str%iPart(i)
  
      do j=i+1, str%nPart
  
  
        ! interparticle vector
        rij(:) = str%rPart(:,i) - str%rPart(:,j)
  
        ! periodic boundary conditions
        rij(:) = rij(:) - str%lBox(:)*nint(rij(:)*lBoxI(:))
  
        ! square of the separation distance
        rr = dot_product(rij,rij)
  
        outer: if( rr < rListSq ) then

          nList = nList+1

          if(nList==fld%nMaxList) then
            open(unit=7,file='output.dat',status='unknown',action='write',position='append')
            write(unit=7,fmt="('Fatal: Verlet list too small')")
            close(unit=7)
            stop
          end if

          fld%vList(nList)=j

          inner: if( rr < rcSq ) then
          
            b = str%iPart(j)
          
            ! Type of interaction
            ab=fld%vdwMatrix(b,a)
          
            ! Calculate forces using the Newton-gregory interpolation (see Allen & Tildesley )
            sds = rr/fld%deltaGrid
            k   = int(sds)
            xi  = sds-dble(k)
            vk0 = fld%ggg(k,ab)
            vk1 = fld%ggg(k+1,ab)
            vk2 = fld%ggg(k+2,ab)
            t1  = vk0+(vk1-vk0)*xi
            t2  = vk1+(vk2-vk1)*(xi-1.0d0)
            wij = t1+(t2-t1)*xi*0.5d0
            str%fPart(:,i) = str%fPart(:,i) + wij*rij(:)
            str%fPart(:,j) = str%fPart(:,j) - wij*rij(:)
          
            ! virial
            prop%virial = prop%virial + wij*rr
          
          
            !
            vk0 = fld%vvv(k,ab)
            vk1 = fld%vvv(k+1,ab)
            vk2 = fld%vvv(k+2,ab)
            t1  = vk0+(vk1-vk0)*xi
            t2  = vk1+(vk2-vk1)*(xi-1.0d0)
            uij = t1+(t2-t1)*xi*0.5d0
          
            prop%potEnergy = prop%potEnergy + uij
  
          end if inner

        end if outer
  
      end do

    end do

    fld%vPoint(str%nPart)=nList+1

    newList=.false.

  else if( fld%lVerlet ) then

    ! Non-bonded interactions using Verlet list
  
    do i=1, str%nPart-1
  
      a = str%iPart(i)

      index1=fld%vPoint(i)
      index2=fld%vPoint(i+1)-1

      if( index1 <= index2 ) then
  
        do jj=index1, index2

          j=fld%vList(jj)
        
          ! interparticle vector
          rij(:) = str%rPart(:,i) - str%rPart(:,j)
        
          ! periodic boundary conditions
          rij(:) = rij(:) - str%lBox(:)*nint(rij(:)*lBoxI(:))
        
          ! square of the separation distance
          rr = dot_product(rij,rij)
        
          if( rr < rcSq ) then
        
            b = str%iPart(j)
        
            ! Type of interaction
            ab=fld%vdwMatrix(b,a)
        
            ! Calculate forces using the Newton-gregory interpolation (see Allen & Tildesley )
            sds = rr/fld%deltaGrid
            k   = int(sds)
            xi  = sds-dble(k)
            vk0 = fld%ggg(k,ab)
            vk1 = fld%ggg(k+1,ab)
            vk2 = fld%ggg(k+2,ab)
            t1  = vk0+(vk1-vk0)*xi
            t2  = vk1+(vk2-vk1)*(xi-1.0d0)
            wij = t1+(t2-t1)*xi*0.5d0
            str%fPart(:,i) = str%fPart(:,i) + wij*rij(:)
            str%fPart(:,j) = str%fPart(:,j) - wij*rij(:)
        
            ! virial
            prop%virial = prop%virial + wij*rr
        
        
            !
            vk0 = fld%vvv(k,ab)
            vk1 = fld%vvv(k+1,ab)
            vk2 = fld%vvv(k+2,ab)
            t1  = vk0+(vk1-vk0)*xi
            t2  = vk1+(vk2-vk1)*(xi-1.0d0)
            uij = t1+(t2-t1)*xi*0.5d0
        
            prop%potEnergy = prop%potEnergy + uij
        
          end if
        
        end do

      end if
  
    end do

  else

    ! Non-bonded interactions without Verlet list
  
    do i=1, str%nPart-1
  
      a = str%iPart(i)
  
      do j=i+1, str%nPart
  
  
        ! interparticle vector
        rij(:) = str%rPart(:,i) - str%rPart(:,j)
  
        ! periodic boundary conditions
        rij(:) = rij(:) - str%lBox(:)*nint(rij(:)*lBoxI(:))
  
        ! square of the separation distance
        rr = dot_product(rij,rij)
  
        if( rr < rcSq ) then
  
          b = str%iPart(j)
  
          ! Type of interaction
          ab=fld%vdwMatrix(b,a)
  
          ! Calculate forces using the Newton-gregory interpolation (see Allen & Tildesley )
          sds = rr/fld%deltaGrid
          k   = int(sds)
          xi  = sds-dble(k)
          vk0 = fld%ggg(k,ab)
          vk1 = fld%ggg(k+1,ab)
          vk2 = fld%ggg(k+2,ab)
          t1  = vk0+(vk1-vk0)*xi
          t2  = vk1+(vk2-vk1)*(xi-1.0d0)
          wij = t1+(t2-t1)*xi*0.5d0
          str%fPart(:,i) = str%fPart(:,i) + wij*rij(:)
          str%fPart(:,j) = str%fPart(:,j) - wij*rij(:)
  
          ! virial
          prop%virial = prop%virial + wij*rr
  
  
          !
          vk0 = fld%vvv(k,ab)
          vk1 = fld%vvv(k+1,ab)
          vk2 = fld%vvv(k+2,ab)
          t1  = vk0+(vk1-vk0)*xi
          t2  = vk1+(vk2-vk1)*(xi-1.0d0)
          uij = t1+(t2-t1)*xi*0.5d0
  
          prop%potEnergy = prop%potEnergy + uij
  
        end if
  
      end do
  
    end do
  end if
  

  ! Bonded interactions

  do i=1, str%nSpecies

    do j=1, str%nMol(i)

      do k=1, str%nSites(i)-1

        ! interparticle vector
        rij(:) = str%rPart(:,k) - str%rPart(:,k+1)
        
        ! periodic boundary conditions
        rij(:) = rij(:) - str%lBox(:)*nint(rij(:)*lBoxI(:))
        
        ! square of the separation distance
        rr = dot_product(rij,rij)

        fVal = -str%hrmK(i)/(1.0d0-rr/str%hrmR(i)/str%hrmR(i))
        str%fPart(:,i) = str%fPart(:,i) + fval*rij(:)
        str%fPart(:,j) = str%fPart(:,j) - fval*rij(:)
        prop%potEnergy = prop%potEnergy - & 
                         0.50d0*str%hrmK(i)*str%hrmR(i)*str%hrmR(i)*dlog(1.0d0-rr/str%hrmR(i)/str%hrmR(i))
        prop%virial    = prop%virial + fval*rr

      end do

    end do

  end do


  ! Restoring spring for particles in the wall

  nWall = sum(str%nMol(1:str%nSpecies-1) )+1  ! First index of the wall

  k=0
  do i=nWall, str%nPart

    k=k+1

    rij= str%rPart(:,i)-str%wPart(:,k)

    rij(:) = rij(:) - str%lBox(:)*nint(rij(:)*lBoxI(:))

    rr = dot_product(rij,rij)

    prop%potEnergy = prop%potEnergy + 0.5d0*str%hrmk(str%nSpecies)*rr

    str%fPart(:,i) = str%fPart(:,i) - str%hrmk(str%nSpecies)*rij(:)
    
    prop%virial = prop%virial + wij*rr

  end do


  ! Applied external field to the fluid particles

  lx = str%lBox(1)*0.5d0-fld%rangeExt

  do i=1, nWall-1

     if( str%rPart(1,i) <= - lx .and.  str%iPart(i) == 1 ) then
       str%fPart(1,i) = str%fPart(1,i) - str%mPart(i)*fld%strengthExt
     else if ( str%rPart(1,i) >= lx .and.  str%iPart(i) == 2 ) then
       str%fPart(1,i) = str%fPart(1,i) + str%mPart(i)*fld%strengthExt
     end if

  end do



end subroutine computeForces
