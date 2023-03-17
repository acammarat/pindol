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
!
! v 0.1
! 
! Antonio Cammarata
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Format of the input file
!
! double double double  reduced components of the first q-point
! ...
! double double double  reduced components of the last q-point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module var

  ! global parameters
  character(3), parameter :: version='0.1'
  character(7), parameter :: progname = 'Q4PHIND'
end module var

module functions
contains
  
  function i2a(i) result(out)
    character(:), allocatable :: out
    integer, intent(in) :: i
    character(range(i)+2) :: x

    write(x,'(i0)') i

    out = trim(x)
    
  end function i2a

end module functions

program q4phind
use functions, only: i2a

  integer :: i, iord
  integer :: npG, npH, npS, nptot, nptmp
  real(8), parameter :: toldelta = 1.d-4
  real(8), allocatable :: vec_red(:,:), vec_tmp(:,:)
  real(8) :: dum
  character(256) :: infile
  logical :: file_exists, wrongq_fl
  logical, allocatable :: vec_ord(:), qconjfl(:)
  
  call show_logo

  call get_command_argument(1,infile)
  if ( (infile == '-h') .or. (infile == '') ) then
    write(*,'(a)') ' Syntax: q4phind <q-point set>'
    write(*,*)
    stop
  end if

  inquire(file=infile,exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if

  write(*,'(*(a))') ' Reading file: ', trim(infile)

  open(unit=10,file=infile,action='read')

  nptmp = 0
  do
    read(10,*,iostat=i) dum
    if ( i == 0 ) then
      nptmp = nptmp + 1
    else if ( i < 0 ) then
      exit
    end if
  end do
  rewind(10)
  write(*,'(*(a))') ' Number of q-points found in input file: ', i2a(nptmp)

  allocate (vec_tmp(nptmp,3), stat=i)
  if (i/=0) stop 'Allocation failed for vec_tmp'
  do i = 1, nptmp
    read(10,*) vec_tmp(i,:)
  end do
  close(10)

  write(*,'(a)') ' Checking vector components...'

  wrongq_fl = .false.
  do k = 1, nptmp
    if ( (vec_tmp(k,1)+0.5d0<=tiny(1.d0)) .or. (vec_tmp(k,1)-0.5d0>tiny(1.d0)) .or. &
         (vec_tmp(k,2)+0.5d0<=tiny(1.d0)) .or. (vec_tmp(k,2)-0.5d0>tiny(1.d0)) .or. &
         (vec_tmp(k,3)+0.5d0<=tiny(1.d0)) .or. (vec_tmp(k,3)-0.5d0>tiny(1.d0)) ) then
      write(*,'(a,a,3(1x,f10.5))') '  components out of (-1/2,1/2] range: q',i2a(k),vec_tmp(k,:)
      wrongq_fl = .true.
    end if
  end do

  if ( wrongq_fl ) then
    write(*,*)
    write(*,'(a)') '  ERROR: remove the vectors with components out of range and rerun.'
    write(*,*)
    stop
  end if

  allocate (qconjfl(nptmp), stat=i)
  if (i/=0) stop 'Allocation failed for qconjfl'
  qconjfl(:) = .false.
  do i = 1, nptmp-1
    do k = i+1, nptmp
      if ( ( abs(vec_tmp(k,1)+vec_tmp(i,1)) < toldelta ) .and. ( abs(vec_tmp(k,2)+vec_tmp(i,2)) < toldelta ) .and. &
           ( abs(vec_tmp(k,3)+vec_tmp(i,3)) < toldelta ) ) then
        write(*,'(*(a))') '  vectors ',i2a(i),' and ',i2a(k), ' are a conjugated couple: the latter will be removed'
        qconjfl(k) = .true.
      end if
    end do
  end do
  write(*,'(a)') ' ...done.'
  
  nptot = 0
  do i = 1, nptmp
    if ( .not. qconjfl(i) ) then
      nptot = nptot + 1
    end if
  end do

  write(*,'(*(a))') ' Number of unique q-points: ', i2a(nptot)

  allocate (vec_red(nptot,3), stat=i)
  if (i/=0) stop 'Allocation failed for vec_red'
  allocate (vec_ord(nptot), stat=i)
  if (i/=0) stop 'Allocation failed for vec_red'

  nptot = 0
  do i = 1, nptmp
    if ( .not. qconjfl(i) ) then
      nptot = nptot + 1
      vec_red(nptot,:) = vec_tmp(i,:)
    end if
  end do

  npG = 0
  npH = 0
  npS = 0
  do i = 1, nptot
    if ( (abs(vec_red(i,1))<tiny(1.d0)) .and. (abs(vec_red(i,2))<tiny(1.d0)) .and. (abs(vec_red(i,3))<tiny(1.d0)) ) then
      npG = 1
    else if ( (abs(abs(vec_red(i,1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(i,2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(i,3))-0.5d0)<tiny(1.d0)) ) then
      npH = npH + 1
    else
      npS = npS + 1
    end if
  end do

  write(*,'(*(a))') ' Q-set composition: G ',i2a(npG),', H ',i2a(npH),', S ',i2a(npS)

  open(unit=20,file='qordered.nd',action='write')

  write(20,'(a)') i2a(nptot)
  vec_ord(:) = .false.
  i = 0
  iord = 0
  do 
    i = i + 1
    if ( i > nptot ) i = 1
    ! look for GM
    if ( npG == 1 .and. (vec_ord(1) .eqv. .false.) ) then
      if ( (abs(vec_red(i,1))<tiny(1.d0)) .and. (abs(vec_red(i,2))<tiny(1.d0)) .and. (abs(vec_red(i,3))<tiny(1.d0)) ) then
        write(20,'(3(f20.15,1x))') vec_red(i,:) 
        iord = iord + 1
        vec_ord(1) = .true.
      else
        cycle
      end if
    ! look for H
    else if ( npH >= 1 .and. (.not. all(vec_ord(npG+1:npG+npH) .eqv. .true.)) ) then
      if ( (abs(abs(vec_red(i,1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(i,2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(i,3))-0.5d0)<tiny(1.d0)) ) then
        write(20,'(3(f20.15,1x))') vec_red(i,:) 
        iord = iord + 1
        vec_ord(iord) = .true.
      else
        cycle
      end if
    ! look for S
    else if ( npS >= 1 .and. (.not. all(vec_ord(npG+npH+1:npG+npH+npS) .eqv. .true.)) ) then
      if ( .not. ((abs(vec_red(i,1))<tiny(1.d0)) .and. (abs(vec_red(i,2))<tiny(1.d0)) .and. (abs(vec_red(i,3))<tiny(1.d0))) .and. &
           .not. ((abs(abs(vec_red(i,1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(i,2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(i,3))-0.5d0)<tiny(1.d0))) ) then
        write(20,'(3(f20.15,1x))') vec_red(i,:) 
        iord = iord + 1
        vec_ord(iord) = .true.
      end if
    else
      exit
    end if
  end do

  close(20)

  write(*,'(a)') ' Ordered q-set written in qordered.nd'
  write(*,*)

  call credits

  stop
end program q4phind
