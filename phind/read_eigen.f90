!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! phind. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Phi for Normal Dynamics
! 
! Calculates the first order anharmonic interaction strength Phi(l,l',l'')
! to be used for the PINDOL code.
! The implemented selection rule is Delta(q+q'+q''); 
! the selection rules described in RSC Advances, 9, 37491 (2019)
! https://doi.org/10.1039/C9RA08294H
! are not implemented.
! Phi(l,l',l'') is not written for any triplet involving at least 
! one of the skip modes as specified in the header of freq.nd
!
! If used for production, you should cite
! A. Cammarata, M. Dasic, P. Nicolini, J. Chem. Phys. 161, 084111 (2024)
! https://doi.org/10.1063/5.0224108
!
!    This file is part of phind.
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

subroutine read_eigen
  use var, only: npunique, nq, eigen, vec, vec_red, &
        toldelta, nskip, skipmod, npG, npH, npS, nptot, &
        npGq, npHq, npSq
  use integer_to_string, only: i2a
  implicit none
  integer :: i, j, k, indpq, indpq1
  real(8) :: tmp_e(2*nq)
  character(1) :: comment
  logical :: file_exists, wrongq_fl
 
  write(*,'(a)') ' Reading polarization vectors...'

  inquire(file='qmatrix.nd',exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(a)') ' ERROR: eigenvector file qmatrix.nd not found.'
    write(*,*)
    stop
  end if

  inquire(file='freq.nd',exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(a)') ' ERROR: frequency file freq.nd not found.'
    write(*,*)
    stop
  end if

  open(unit=20,file='freq.nd',action='READ')
  read(20,*)
  read(20,*) comment, comment, npG, comment, npH, comment, npS, comment, i

  if ( i /= nq ) then
    write(*,'(a)') ' ERROR: the number of modes read in freq.nd is not consistent with the one calculated from the geometry file.'
    write(*,*)
    stop
  end if

  ! number of "unique" q-points (i.e. Gamma + set H + set S) 
  npunique = npG + npH + npS
  write(*,'(*(a))') '  q-set composition: G ',i2a(npG),', H ',i2a(npH),', S ',i2a(npS)
  write(*,'(*(a))') '  Number of q-points: ', i2a(npunique)

  ! set total number of q-points (i.e. Gamma + set H + set S + set S*)
  nptot = npG + npH + 2*npS

  ! set auxiliary numbers of degrees of freedom for the normal dynamics
  npGq = npG * nq
  npHq = npH * nq
  npSq = npS * nq

  allocate ( vec(nptot,3), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for VEC'
  allocate ( vec_red(nptot), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for VEC_RED'
  allocate ( eigen(nptot*nq,nq), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for EIGEN'

  read(20,*) comment, comment, comment, nskip
  allocate ( skipmod(nskip,2), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for skipmod'
  backspace(20)
  read(20,*) comment, comment, comment, nskip, (skipmod(i,1),skipmod(i,2),i=1,nskip)

  write(*,'(*(a))',advance='no') '  Skipping ',i2a(nskip),' modes:'
  do i = 1, nskip-1
    write(*,'(*(a))', advance='no')  '  ',i2a(skipmod(i,1)),' ',i2a(skipmod(i,2))
  end do
  write(*,'(*(a))')  '  ',i2a(skipmod(nskip,1)),' ',i2a(skipmod(nskip,2))

  read(20,*)

  ! Check that the q-set is ordered as GM, H, S, according to npG, npH and npS.
  ! The working assumption in PINDOL is that all the q-components belong to the interval (-1/2,1/2]

  wrongq_fl = .false.
  ! check eigen(GM)
  do k = 1, npG ! The do loop is here used to preserve the general structure of the eigenreading
    read(20,*) comment, (vec(k,i),i=1,3) ! Ang^-1
    read(20,*) comment, indpq1, (vec_red(k)%ei(i),i=1,3) ! red. coord. - indpq1 here is a dummy var
    vec_red(k)%set = 'G'
    if ( .not. (abs(vec_red(k)%ei(1))<tiny(1.d0) .and. abs(vec_red(k)%ei(2))<tiny(1.d0) .and. abs(vec_red(k)%ei(3))<tiny(1.d0)) ) then
      write(*,'(a,a,a,3(f10.5,1x))') ' ERROR: Gamma point expected but found q',i2a(k),' ',vec_red(k)%ei(:)
      wrongq_fl = .true.
    end if
    if ( (vec_red(k)%ei(1)+0.5d0<=tiny(1.d0)) .or. (vec_red(k)%ei(1)-0.5d0>tiny(1.d0)) .or. & 
         (vec_red(k)%ei(2)+0.5d0<=tiny(1.d0)) .or. (vec_red(k)%ei(2)-0.5d0>tiny(1.d0)) .or. &
         (vec_red(k)%ei(3)+0.5d0<=tiny(1.d0)) .or. (vec_red(k)%ei(3)-0.5d0>tiny(1.d0)) ) then
      write(*,'(a,a,a,3(f10.5,1x))') ' ERROR: q-components out of (-1/2,1/2] range: q',i2a(k),' ',vec_red(k)%ei(:)
      wrongq_fl = .true.
    end if
    do i = 1, nq
      read(20,*)
    end do
  end do

  ! check eigen(H)
  do k = npG+1, npG+npH
    read(20,*) comment, (vec(k,i),i=1,3) ! Ang^-1
    read(20,*) comment, indpq1, (vec_red(k)%ei(i),i=1,3) ! red. coord. - indpq1 here is a dummy var
    vec_red(k)%set = 'H'
    if ( .not. ((abs(abs(vec_red(k)%ei(1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(k)%ei(2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(k)%ei(3))-0.5d0)<tiny(1.d0)) ) ) then
      write(*,'(a,a,a,3(f10.5,1x))') ' ERROR: q-point in set H expected but found q',i2a(k),' ',vec_red(k)%ei(:)
      wrongq_fl = .true.
    end if
    if ( (vec_red(k)%ei(1)+0.5d0<=tiny(1.d0)) .or. (vec_red(k)%ei(1)-0.5d0>tiny(1.d0)) .or. & 
         (vec_red(k)%ei(2)+0.5d0<=tiny(1.d0)) .or. (vec_red(k)%ei(2)-0.5d0>tiny(1.d0)) .or. &
         (vec_red(k)%ei(3)+0.5d0<=tiny(1.d0)) .or. (vec_red(k)%ei(3)-0.5d0>tiny(1.d0)) ) then
      write(*,'(a,a,a,3(f10.5,1x))') ' ERROR: q-components out of (-1/2,1/2] range: q',i2a(k),' ',vec_red(k)%ei(:)
      wrongq_fl = .true.
    end if
    do i = 1, nq
      read(20,*)
    end do
  end do

  ! check eigen(S)
  do k = npG+npH+1, npG+npH+npS
    read(20,*) comment, (vec(k,i),i=1,3) ! Ang^-1
    read(20,*) comment, indpq1, (vec_red(k)%ei(i),i=1,3) ! red. coord. - indpq1 here is a dummy var
    vec_red(k)%set = 'S'
    if ( (abs(vec_red(k)%ei(1))<tiny(1.d0) .and. abs(vec_red(k)%ei(2))<tiny(1.d0) .and. abs(vec_red(k)%ei(3))<tiny(1.d0)) .or. &
         ((abs(abs(vec_red(k)%ei(1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(k)%ei(2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(vec_red(k)%ei(3))-0.5d0)<tiny(1.d0))) ) then
      write(*,'(a,a,a,3(f10.5,1x))') ' ERROR: q-point in set S expected but found q',i2a(k),' ',vec_red(k)%ei(:)
      wrongq_fl = .true.
    end if
    if ( (vec_red(k)%ei(1)+0.5d0<=tiny(1.d0)) .or. (vec_red(k)%ei(1)-0.5d0>tiny(1.d0)) .or. & 
         (vec_red(k)%ei(2)+0.5d0<=tiny(1.d0)) .or. (vec_red(k)%ei(2)-0.5d0>tiny(1.d0)) .or. &
         (vec_red(k)%ei(3)+0.5d0<=tiny(1.d0)) .or. (vec_red(k)%ei(3)-0.5d0>tiny(1.d0)) ) then
      write(*,'(a,a,3(1x,f10.5))') ' ERROR: q-components out of (-1/2,1/2] range: q',i2a(k),vec_red(k)%ei(:)
      wrongq_fl = .true.
    end if
    do i = 1, nq
      read(20,*)
    end do
  end do

  if (wrongq_fl) then
    write(*,*)
    stop
  end if
  
  rewind(20)
  do i = 1, 4
    read(20,*)
  end do

  ! read eigenvectors and frequencies
  open(unit=10,file='qmatrix.nd',action='READ')
  read(10,*)

  write(*,'(a)') '  Reading eigenvectors and frequencies at'
  do k = 1, npunique 
    read(20,*) comment, (vec(k,i),i=1,3) ! Ang^-1
    read(20,*) comment, indpq1, (vec_red(k)%ei(i),i=1,3) ! red. coord. - indpq1 here is a dummy var
    write(*,'(a,a,a,3(f10.5,1x),a,a)') '   q',i2a(k), ' ',vec_red(k)%ei(:), ' ',vec_red(k)%set
    do i = 1, nq
      read(10,*) tmp_e(:)
      indpq = (k-1)*nq + i
      indpq1 = 0
      do j = 1, 2*nq-1, 2
        indpq1 = indpq1 + 1
        eigen(indpq,indpq1) = dcmplx( tmp_e(j), tmp_e(j+1) ) ! adim
      end do
      read(20,*)
    end do
  end do

  close(20)
  close(10)
  write(*,'(a)') ' done.'

  ! check if the provided q-point set
  ! contains any complex conjugated couple

  do k = 1 , npunique-1
    do i = k+1, npunique
      if ( ( abs(vec_red(k)%ei(1)+vec_red(i)%ei(1)) < toldelta ) .and. ( abs(vec_red(k)%ei(2)+vec_red(i)%ei(2)) < toldelta ) .and. &
           ( abs(vec_red(k)%ei(3)+vec_red(i)%ei(3)) < toldelta ) ) then
        write(*,'(*(a))') ' ERROR: vectors ',i2a(k),' and ',i2a(i), ' are a conjugated couple.'
        wrongq_fl = .true.
      end if
    end do
  end do

  if (wrongq_fl) then
    write(*,*)
    stop
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! creation of complex conjugated eigenvectors

  write(*,'(a)',advance='no') ' Creating conjugated eigenvectors...'
  do k = npunique+1, npunique+npS
    vec(k,:) = -vec(k-npS,:) ! Ang^-1
    vec_red(k)%ei(:) = -vec_red(k-npS)%ei(:) ! Ang^-1
!    write(*,*) vec_red(k-npS)%ei(:), k-npS
    vec_red(k)%set = 'C'
    do i = 1, nq
       indpq = (k-1)*nq + i
       indpq1 = (k-npS-1)*nq + i
       eigen(indpq,:) = conjg(eigen(indpq1,:)) ! adim
!       write(*,*) eigen(indpq,:)
    end do
  end do
  write(*,'(a)') ' done.'

!temporary check

!  open(100,file='test_eigenstar.tmp',action='write')
!  do k = (npunique)*nq+1, (npunique+npS)*nq
!    write(100,'(*(f17.14,1x))') eigen(k,:)
!  end do
!  close(100)


  return
end subroutine read_eigen
