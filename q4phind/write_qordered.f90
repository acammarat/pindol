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

subroutine write_qordered
use functions, only: i2a
use var, only: npG, npH, npS, nptot, vec_red

  integer :: i, iord
  real(8) :: qHcheck(3)
  logical, allocatable :: vec_ord_fl(:)
  
  allocate (vec_ord_fl(nptot), stat=i)
  if (i/=0) stop 'Allocation failed for vec_ord'
  open(unit=20,file='qordered.nd',action='write')

  write(20,'(a)') i2a(nptot)
  vec_ord_fl(:) = .false.
  i = 0
  iord = 0
  do 
    i = i + 1
    if ( i > nptot ) i = 1

    ! q at the border of the Brillouin zone belongs to H only if its components are 0 or 1/2
    ! if q belongs to H then qHcheck components are integer numbers
    qHcheck(:) = 2.d0 * vec_red(i,:)

    ! look for GM
    if ( npG == 1 .and. (vec_ord_fl(1) .eqv. .false.) ) then
      if ( (abs(vec_red(i,1))<tiny(1.d0)) .and. (abs(vec_red(i,2))<tiny(1.d0)) .and. (abs(vec_red(i,3))<tiny(1.d0)) ) then
        write(20,'(3(f20.15,1x))') vec_red(i,:) 
        iord = iord + 1
        vec_ord_fl(1) = .true.
      else
        cycle
      end if
    ! look for H
    else if ( npH >= 1 .and. (.not. all(vec_ord_fl(npG+1:npG+npH) .eqv. .true.)) ) then
      if ( (abs(qHcheck(1)-floor(qHcheck(1)))<tiny(1.d0)) .and. (abs(qHcheck(2)-floor(qHcheck(2)))<tiny(1.d0)) &
           .and. (abs(qHcheck(3)-floor(qHcheck(3)))<tiny(1.d0)) ) then
        write(20,'(3(f20.15,1x))') vec_red(i,:) 
        iord = iord + 1
        vec_ord_fl(iord) = .true.
      else
        cycle
      end if
    ! look for S
    else if ( npS >= 1 .and. (.not. all(vec_ord_fl(npG+npH+1:npG+npH+npS) .eqv. .true.)) ) then
      if ( .not. ((abs(vec_red(i,1))<tiny(1.d0)) .and. (abs(vec_red(i,2))<tiny(1.d0)) .and. (abs(vec_red(i,3))<tiny(1.d0))) .and. &
           .not. ((abs(qHcheck(1)-floor(qHcheck(1)))<tiny(1.d0)) .and. (abs(qHcheck(2)-floor(qHcheck(2)))<tiny(1.d0)) .and. &
           (abs(qHcheck(3)-floor(qHcheck(3)))<tiny(1.d0))) &
           ) then
        write(20,'(3(f20.15,1x))') vec_red(i,:) 
        iord = iord + 1
        vec_ord_fl(iord) = .true.
      end if
    else
      exit
    end if
  end do

  close(20)

  write(*,'(a)') ' Ordered q-set written in qordered.nd'
  write(*,*)


end subroutine write_qordered
