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

subroutine check_qopt
use functions, only: i2a
use var, only: npG, npH, npS, npopt, vec_opt, vec_red_opt, nptot

  integer :: i, k
  real(8), parameter :: toldelta = 1.d-4
  logical :: wrongq_fl
  logical :: qopt_rem_fl(npopt)
!  logical, allocatable :: qopt_rem_fl(:)
  
!  write(*,'(a)') '  Checking vector components...'

  wrongq_fl = .false.
  do k = 1, npopt
    if ( (vec_opt(k,1)+0.5d0<=tiny(1.d0)) .or. (vec_opt(k,1)-0.5d0>tiny(1.d0)) .or. &
         (vec_opt(k,2)+0.5d0<=tiny(1.d0)) .or. (vec_opt(k,2)-0.5d0>tiny(1.d0)) .or. &
         (vec_opt(k,3)+0.5d0<=tiny(1.d0)) .or. (vec_opt(k,3)-0.5d0>tiny(1.d0)) ) then
      write(*,'(a,a,3(1x,f10.5))') '   components out of (-1/2,1/2] range: q',i2a(k),vec_opt(k,:)
      wrongq_fl = .true.
    end if
  end do

  if ( wrongq_fl ) then
    write(*,*)
    write(*,'(a)') '   ERROR: remove the vectors with components out of range and rerun.'
    write(*,*)
    stop
  end if

!  allocate (qopt_rem_fl(npopt), stat=i)
!  if (i/=0) stop 'Allocation failed for qopt_rem_fl'
  qopt_rem_fl(:) = .false.

  do i = 1, npopt-1
    do k = i+1, npopt

      ! check for conjugated couples
      if ( ( abs(vec_opt(k,1)+vec_opt(i,1)) < toldelta ) .and. ( abs(vec_opt(k,2)+vec_opt(i,2)) < toldelta ) .and. &
           ( abs(vec_opt(k,3)+vec_opt(i,3)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qopt_rem_fl(k) = .true.
      ! specific check for zone border
      else if ( ( abs(abs(vec_opt(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(vec_opt(k,2)+vec_opt(i,2)) < toldelta ) .and. &
           ( abs(vec_opt(k,3)+vec_opt(i,3)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qopt_rem_fl(k) = .true.
      else if ( ( abs(abs(vec_opt(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(vec_opt(k,1)+vec_opt(i,1)) < toldelta ) .and. &
           ( abs(vec_opt(k,3)+vec_opt(i,3)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qopt_rem_fl(k) = .true.
      else if ( ( abs(abs(vec_opt(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(vec_opt(k,2)+vec_opt(i,2)) < toldelta ) .and. &
           ( abs(vec_opt(k,1)+vec_opt(i,1)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qopt_rem_fl(k) = .true.

      else if ( ( abs(abs(vec_opt(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_opt(k,3)+vec_opt(i,3)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qopt_rem_fl(k) = .true.
      else if ( ( abs(abs(vec_opt(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_opt(k,2)+vec_opt(i,2)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qopt_rem_fl(k) = .true.
      else if ( ( abs(abs(vec_opt(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_opt(k,1)+vec_opt(i,1)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qopt_rem_fl(k) = .true.

      ! check for identical vectors
      else if ( ( abs(vec_opt(k,1)-vec_opt(i,1)) < toldelta ) .and. ( abs(vec_opt(k,2)-vec_opt(i,2)) < toldelta ) .and. &
           ( abs(vec_opt(k,3)-vec_opt(i,3)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qopt_rem_fl(k) = .true.
      ! specific check for zone border
      else if ( ( abs(abs(vec_opt(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(vec_opt(k,2)-vec_opt(i,2)) < toldelta ) .and. &
           ( abs(vec_opt(k,3)-vec_opt(i,3)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qopt_rem_fl(k) = .true.
      else if ( ( abs(abs(vec_opt(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(vec_opt(k,1)-vec_opt(i,1)) < toldelta ) .and. &
           ( abs(vec_opt(k,3)-vec_opt(i,3)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qopt_rem_fl(k) = .true.
      else if ( ( abs(abs(vec_opt(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(vec_opt(k,2)-vec_opt(i,2)) < toldelta ) .and. &
           ( abs(vec_opt(k,1)-vec_opt(i,1)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qopt_rem_fl(k) = .true.

      else if ( ( abs(abs(vec_opt(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_opt(k,3)-vec_opt(i,3)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qopt_rem_fl(k) = .true.
      else if ( ( abs(abs(vec_opt(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_opt(k,2)-vec_opt(i,2)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qopt_rem_fl(k) = .true.
      else if ( ( abs(abs(vec_opt(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,3)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(vec_opt(k,1)-vec_opt(i,1)) < toldelta ) ) then
!        write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qopt_rem_fl(k) = .true.

      else if ( ( abs(abs(vec_opt(k,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,1)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(k,2)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,2)) - 0.5d0) < tiny(1.d0) ) .and. &
           ( abs(abs(vec_opt(k,3)) - 0.5d0) < tiny(1.d0) ) .and. ( abs(abs(vec_opt(i,3)) - 0.5d0) < tiny(1.d0) ) ) then
 !       write(*,'(*(a))') '   vectors ',i2a(i),' and ',i2a(k), ' are identical: the latter will be removed'
        qopt_rem_fl(k) = .true.
      end if
    end do
  end do
!  write(*,'(a)') '  ...done.'
  
  nptot = 0
  do i = 1, npopt
    if ( .not. qopt_rem_fl(i) ) then
      nptot = nptot + 1
    end if
  end do

  write(*,'(*(a))') '  Number of unique optimal q-points: ', i2a(nptot)

  allocate (vec_red_opt(nptot,3), stat=i)
  if (i/=0) stop 'Allocation failed for vec_red_opt'

  nptot = 0
  do i = 1, npopt
    if ( .not. qopt_rem_fl(i) ) then
      nptot = nptot + 1
      vec_red_opt(nptot,:) = vec_opt(i,:)
    end if
  end do

  npG = 0
  npH = 0
  npS = 0
  do i = 1, nptot
    if ( (abs(vec_red_opt(i,1))<tiny(1.d0)) .and. (abs(vec_red_opt(i,2))<tiny(1.d0)) .and. (abs(vec_red_opt(i,3))<tiny(1.d0)) ) then
      npG = 1
    else if ( (abs(abs(vec_red_opt(i,1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red_opt(i,2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red_opt(i,3))-0.5d0)<tiny(1.d0)) ) then
      npH = npH + 1
    else
      npS = npS + 1
    end if
  end do

  write(*,'(*(a))') '  Q-set composition: G ',i2a(npG),', H ',i2a(npH),', S ',i2a(npS)

!  deallocate (qopt_rem_fl, stat=i)
!  if (i/=0) stop 'Deallocation failed for qopt_rem_fl'

end subroutine check_qopt
