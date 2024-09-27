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

subroutine write_qopt
use functions, only: i2a
use var, only: vec_opt, vec_red, npG, npH, npS, nptot, &
        npopt, vec_red_opt

  integer :: i, j, k, l, iord
  logical, allocatable :: vec_ord_fl(:)

  npopt = nptot
  do i = 1, nptot
    do k = i, nptot
      npopt = npopt + 1
    end do
  end do

  allocate ( vec_opt(npopt,3), stat = i )
  if ( i /= 0 ) stop ' ERROR: allocation failed for vec_opt'

  do i = 1, nptot
    vec_opt(i,:) = vec_red(i,:)
  end do

  write(*,'(a)') ' Calculating optimal q-set:'
  k = nptot
  do i = 1, nptot
    do j = i, nptot
      k = k + 1
      vec_opt(k,:) = -vec_red(i,:)-vec_red(j,:)
!  write(*,'(i3,9(f10.6,1x))') k, vec_opt(k,:), vec_red(i,:), vec_red(j,:)
!      vec_opt(K,:) = 1.d0 - sign(vec_opt(k,:),vec_opt(k,:))
      do l = 1, 3
        if ( vec_opt(k,l) > 0.5d0 ) vec_opt(k,l) = 1.d0 - vec_opt(k,l)
        if ( vec_opt(k,l) <= -0.5d0 ) vec_opt(k,l) = vec_opt(k,l) + 1.d0
      end do
    end do
  end do

!open(200,file='qopt_test.nd',action='write')
!
!do i = 1, npopt
!  write(200,'(i3,3(f10.6,1x))') i, vec_opt(i,:)
!end do
!close(200)
!stop

  call check_qopt

  open(unit=20,file='qoptimal.nd',action='write')

  write(20,'(a)') i2a(nptot)
  allocate ( vec_ord_fl(nptot), stat = i )
  if ( i /= 0 ) stop ' ERROR: allocation failed for vec_opt'
  vec_ord_fl(:) = .false.
  i = 0
  iord = 0
  do 
    i = i + 1
    if ( i > nptot ) i = 1
    ! look for GM
    if ( npG == 1 .and. (vec_ord_fl(1) .eqv. .false.) ) then
      if ( (abs(vec_red_opt(i,1))<tiny(1.d0)) .and. (abs(vec_red_opt(i,2))<tiny(1.d0)) .and. (abs(vec_red_opt(i,3))<tiny(1.d0)) ) then
        do j = 1, 3
          if ( abs(vec_red_opt(i,j)) < tiny(1.d0) ) vec_red_opt(i,j) = 0.d0
        end do
        write(20,'(3(f20.15,1x))') vec_red_opt(i,:) 
        iord = iord + 1
        vec_ord_fl(1) = .true.
      else
        cycle
      end if
    ! look for H
    else if ( npH >= 1 .and. (.not. all(vec_ord_fl(npG+1:npG+npH) .eqv. .true.)) ) then
      if ( (abs(abs(vec_red_opt(i,1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red_opt(i,2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red_opt(i,3))-0.5d0)<tiny(1.d0)) ) then
        do j = 1, 3
          if ( abs(vec_red_opt(i,j)) < tiny(1.d0) ) vec_red_opt(i,j) = 0.d0
        end do
        write(20,'(3(f20.15,1x))') vec_red_opt(i,:) 
        iord = iord + 1
        vec_ord_fl(iord) = .true.
      else
        cycle
      end if
    ! look for S
    else if ( npS >= 1 .and. (.not. all(vec_ord_fl(npG+npH+1:npG+npH+npS) .eqv. .true.)) ) then
      if ( .not. ((abs(vec_red_opt(i,1))<tiny(1.d0)) .and. (abs(vec_red_opt(i,2))<tiny(1.d0)) .and. (abs(vec_red_opt(i,3))<tiny(1.d0))) .and. &
           .not. ((abs(abs(vec_red_opt(i,1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red_opt(i,2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red_opt(i,3))-0.5d0)<tiny(1.d0))) ) then
        do j = 1, 3
          if ( abs(vec_red_opt(i,j)) < tiny(1.d0) ) vec_red_opt(i,j) = 0.d0
        end do
        write(20,'(3(f20.15,1x))') vec_red_opt(i,:) 
        iord = iord + 1
        vec_ord_fl(iord) = .true.
      end if
    else
      exit
    end if
  end do

  close(20)

  write(*,'(a)') '  Optimal q-set written in qoptimal.nd'
  write(*,*)

  deallocate ( vec_ord_fl, stat = i )
  if ( i /= 0 ) stop ' ERROR: Deallocation failed for vec_ord_fl'

end subroutine write_qopt
