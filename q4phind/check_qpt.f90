!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! q4phind. Copyright (C) Antonio Cammarata                            
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Program to calculate the atomic character
! Order a q-point set as GM, H, S and remove possible complex-conjugated couples
! to generate qmatrix.nd and freq.nd compatible with phind v >= 4.3
!
! If used for production, you should cite
! A. Cammarata, M. Dasic, P. Nicolini, J. Chem. Phys. 161, 084111 (2024)
! https://doi.org/10.1063/5.0224108
!
!    This file is part of q4phind.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine check_qpt
use functions, only: i2a
use var, only: npG, npH, npS, nptot, vec_inp, vec_red, npinp

  integer :: i, k
  real(8), parameter :: toldelta = 1.d-4
  real(8) :: qHcheck(3)
  logical :: wrongq_fl
  logical :: qrem_fl(npinp)
!  logical, allocatable :: qrem_fl(:)
  
  write(*,'(a)') '  Checking vector components...'

  wrongq_fl = .false.
  do k = 1, npinp
    if ( (vec_inp(k,1)+0.5d0<=tiny(1.d0)) .or. (vec_inp(k,1)-0.5d0>tiny(1.d0)) .or. &
         (vec_inp(k,2)+0.5d0<=tiny(1.d0)) .or. (vec_inp(k,2)-0.5d0>tiny(1.d0)) .or. &
         (vec_inp(k,3)+0.5d0<=tiny(1.d0)) .or. (vec_inp(k,3)-0.5d0>tiny(1.d0)) ) then
      write(*,'(a,a,3(1x,f10.5))') '   components out of (-1/2,1/2] range: q',i2a(k),vec_inp(k,:)
      wrongq_fl = .true.
    end if
  end do

  if ( wrongq_fl ) then
    write(*,*)
    write(*,'(a)') '   ERROR: remove the vectors with components out of range and rerun.'
    write(*,*)
    stop
  end if

  qrem_fl(:) = .false.

  do i = 1, npinp-1
    do k = i+1, npinp

      ! check for conjugated couples
      if ( ( abs(vec_inp(k,1)+vec_inp(i,1)) < toldelta ) .and. ( abs(vec_inp(k,2)+vec_inp(i,2)) < toldelta ) .and. &
           ( abs(vec_inp(k,3)+vec_inp(i,3)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qrem_fl(k) = .true.
      ! check for conjugated couples at zone border:
      ! if one of the component is 1/2 and the other corresponding ones are opposite, then the two vectors are a conjugated couple;
      ! for example, q1(1/2,1/4,0) is conjugated couple with q2=-q1=(-1/2,-1/4,0). However, q3=(1/2,-1/4,0) is equivalent to q2 because of
      ! the periodicity of the Brillouin zone, as 1/2=-1/2, so q1 and q3 are conjugated couple too.
      else if ( ( abs(abs(vec_inp(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,1)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(vec_inp(k,2)+vec_inp(i,2)) < toldelta ) .and. &
           ( abs(vec_inp(k,3)+vec_inp(i,3)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qrem_fl(k) = .true.
      else if ( ( abs(abs(vec_inp(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(vec_inp(k,1)+vec_inp(i,1)) < toldelta ) .and. &
           ( abs(vec_inp(k,3)+vec_inp(i,3)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qrem_fl(k) = .true.
      else if ( ( abs(abs(vec_inp(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(vec_inp(k,2)+vec_inp(i,2)) < toldelta ) .and. &
           ( abs(vec_inp(k,1)+vec_inp(i,1)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qrem_fl(k) = .true.

      else if ( ( abs(abs(vec_inp(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,1)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(abs(vec_inp(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_inp(k,3)+vec_inp(i,3)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qrem_fl(k) = .true.
      else if ( ( abs(abs(vec_inp(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,1)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(abs(vec_inp(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_inp(k,2)+vec_inp(i,2)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qrem_fl(k) = .true.
      else if ( ( abs(abs(vec_inp(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(abs(vec_inp(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_inp(k,1)+vec_inp(i,1)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qrem_fl(k) = .true.


      ! check for identical vectors
      else if ( ( abs(vec_inp(k,1)-vec_inp(i,1)) < toldelta ) .and. ( abs(vec_inp(k,2)-vec_inp(i,2)) < toldelta ) .and. &
           ( abs(vec_inp(k,3)-vec_inp(i,3)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qrem_fl(k) = .true.
      ! check for identical vectors at zone border
      else if ( ( abs(abs(vec_inp(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,1)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(vec_inp(k,2)-vec_inp(i,2)) < toldelta ) .and. &
           ( abs(vec_inp(k,3)-vec_inp(i,3)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qrem_fl(k) = .true.
      else if ( ( abs(abs(vec_inp(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(vec_inp(k,1)-vec_inp(i,1)) < toldelta ) .and. &
           ( abs(vec_inp(k,3)-vec_inp(i,3)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qrem_fl(k) = .true.
      else if ( ( abs(abs(vec_inp(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(vec_inp(k,2)-vec_inp(i,2)) < toldelta ) .and. &
           ( abs(vec_inp(k,1)-vec_inp(i,1)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qrem_fl(k) = .true.

      else if ( ( abs(abs(vec_inp(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,1)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(abs(vec_inp(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_inp(k,3)-vec_inp(i,3)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qrem_fl(k) = .true.
      else if ( ( abs(abs(vec_inp(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,1)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(abs(vec_inp(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_inp(k,2)-vec_inp(i,2)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qrem_fl(k) = .true.
      else if ( ( abs(abs(vec_inp(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(abs(vec_inp(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_inp(k,1)-vec_inp(i,1)) < toldelta ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qrem_fl(k) = .true.

      else if ( ( abs(abs(vec_inp(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,1)) - 0.5d0) < tiny(1.d0) ) .and. &
      ( abs(abs(vec_inp(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(abs(vec_inp(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_inp(i,3)) - 0.5d0) < tiny(1.d0) ) ) then
        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qrem_fl(k) = .true.
      end if
    end do
  end do
  write(*,'(a)') '  ...done.'
  
  nptot = 0
  do i = 1, npinp
    if ( .not. qrem_fl(i) ) then
      nptot = nptot + 1
    end if
  end do

  write(*,'(*(a))') '  Number of unique q-points: ', i2a(nptot)

  allocate (vec_red(nptot,3), stat=i)
  if (i/=0) stop 'Allocation failed for vec_red'

  nptot = 0
  do i = 1, npinp
    if ( .not. qrem_fl(i) ) then
      nptot = nptot + 1
      vec_red(nptot,:) = vec_inp(i,:)
!      write(*,'(i2,1x,3(f15.10,1x))') i, vec_red(nptot,:)
    end if
  end do

  npG = 0
  npH = 0
  npS = 0
  do i = 1, nptot

    ! q at the border of the Brillouin zone belongs to H only if its components are 0 or 1/2
    ! if q belongs to H then qHcheck components are integer numbers
    qHcheck(:) = 2.d0 * vec_red(i,:)

    if ( (abs(vec_red(i,1))<tiny(1.d0)) .and. (abs(vec_red(i,2))<tiny(1.d0)) .and. (abs(vec_red(i,3))<tiny(1.d0)) ) then
      npG = 1
    else if ( (abs(qHcheck(1)-floor(qHcheck(1)))<tiny(1.d0)) .and. (abs(qHcheck(2)-floor(qHcheck(2)))<tiny(1.d0)) &
              .and. (abs(qHcheck(3)-floor(qHcheck(3)))<tiny(1.d0)) ) then
      npH = npH + 1
    else
      npS = npS + 1
    end if
  end do

  write(*,'(*(a))') '  Q-set composition: G ',i2a(npG),', H ',i2a(npH),', S ',i2a(npS)

end subroutine check_qpt
