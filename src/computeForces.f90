subroutine computeForces( str, fld, prop, newList )

  ! Purpose: Calculation of inter and intra molecular forces

  use global
  implicit none

  type(structure),  intent(inout) :: str
  type(field),      intent(inout)    :: fld 
  type(properties), intent(out)   :: prop

  integer :: i, j, jj, k, a, b, ab, nList
  integer :: index1, index2


  real(8) :: sds, xi, vk0, vk1, vk2, t1, t2, uij, wij
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
            open(unit=7,file='output.dat',status='old',position='append')
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

  jj=0

  do i=1, str%nSpecies

    do j=1, str%nMol(i)

      do k=1, str%nSites(i)-1

        jj = jj + 1

        ! interparticle vector
        rij(:) = str%rPart(:,jj) - str%rPart(:,jj+1)
        
        ! periodic boundary conditions
        rij(:) = rij(:) - str%lBox(:)*nint(rij(:)*lBoxI(:))
        
        ! square of the separation distance
        rr = dot_product(rij,rij)

        fVal = -str%hrmK(i)/(1.0d0-rr/str%hrmR(i)/str%hrmR(i))
        str%fPart(:,jj) = str%fPart(:,jj) + fval*rij(:)
        str%fPart(:,jj+1) = str%fPart(:,jj+1) - fval*rij(:)
        prop%potEnergy = prop%potEnergy - & 
                         0.50d0*str%hrmK(i)*str%hrmR(i)*str%hrmR(i)*dlog(1.0d0-rr/str%hrmR(i)/str%hrmR(i))
        prop%virial    = prop%virial + fval*rr

      end do

      jj = jj + 1

    end do

  end do



end subroutine computeForces
