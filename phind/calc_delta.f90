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
! Phys. Rev. B XX, XXXXX (XXXX)
! https://doi.org/10.1103/xxx
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

subroutine calc_delta
  use integer_to_string
  use var, only: nq, vec_red, toldelta, &
       delta, ll, nskip, skipmod, &
       npGq, npHq, npSq, nptot, npS

  implicit none
  integer :: i, j1, j2, j3, jj3, k3c, l3c
  integer :: k1, k2, k3, l, l1, l2
  integer(8) :: il
  real(8) :: vecsum(3) 
  logical :: fl_qallow_written(nptot,nptot,nptot)
  character(256) :: command

  fl_qallow_written(:,:,:) = .false.

  write(*,'(a)', advance='no') ' Writing Delta-allowed q-triplets in phi.qallowed.nd ...'
  open(unit=10,file='phi.qallowed.nd',action='write')

  write(10,'(a)') '# Allowed q-triplets:'
  write(10,'(a)') '# q1   q2   q3            vecsum'

  ! Counting the unique Phi_(l,l',l'') elements to allocate delta.
  ! Only the possible allowed (q,q',q'') triplets are evaluated.
  ! The working assumption is that all the q-components belong to 
  ! the interval (-1/2,1/2]; this is checked in read_eigen.

  open(unit=11,file='tmp_delta.nd',action='write')
  ll = 0

  ! GM, GM, GM
  ql1: do l = 1, npGq
    k1 = ceiling(dble(l)/dble(nq))
    j1 = l - (k1-1)*nq
    do i = 1, nskip
      if ( k1 == skipmod(i,1) .and. j1 == skipmod(i,2) ) cycle ql1
    end do
    ql11: do l1 = l, npGq
      k2 = ceiling(dble(l1)/dble(nq))
      j2 = l1 - (k2-1)*nq
      do i = 1, nskip
        if ( k2 == skipmod(i,1) .and. j2 == skipmod(i,2) ) cycle ql11
      end do
      ql21: do l2 = l1, npGq
        k3 = ceiling(dble(l2)/dble(nq))
        j3 = l2 - (k3-1)*nq
        do i = 1, nskip
          if ( k3 == skipmod(i,1) .and. j3 == skipmod(i,2) ) cycle ql21
        end do

        ! in this case the check is not necessary as Delta(GM,GM,GM) = 1
        vecsum(1:3) = vec_red(k1)%ei(1:3) + vec_red(k2)%ei(1:3) + vec_red(k3)%ei(1:3)
        !if ( (abs(vecsum(1) - nint(vecsum(1))) < toldelta) .and. (abs(vecsum(2) - nint(vecsum(2))) < toldelta) .and. (abs(vecsum(3) - nint(vecsum(3))) < toldelta) ) then
          if ( .not. fl_qallow_written(k1,k2,k3) ) then
            write(10,'(3(i4,1x),3(f8.5,1x))') k1, k2, k3, vecsum(:)
            fl_qallow_written(k1,k2,k3) = .true.
          end if
          ll = ll + 1
          write(11,*) l, l1, l2
        !end if
      end do ql21
    end do ql11
  end do ql1

  ! GM, H, H
  ql2: do l = 1, npGq
    k1 = ceiling(dble(l)/dble(nq))
    j1 = l - (k1-1)*nq
    do i = 1, nskip
      if ( k1 == skipmod(i,1) .and. j1 == skipmod(i,2) ) cycle ql2
    end do
    ql12: do l1 = npGq + 1, npGq + npHq
      k2 = ceiling(dble(l1)/dble(nq))
      j2 = l1 - (k2-1)*nq
      do i = 1, nskip
        if ( k2 == skipmod(i,1) .and. j2 == skipmod(i,2) ) cycle ql12
      end do
      ql22: do l2 = l1, npGq + npHq
        k3 = ceiling(dble(l2)/dble(nq))
        j3 = l2 - (k3-1)*nq
        do i = 1, nskip
          if ( k3 == skipmod(i,1) .and. j3 == skipmod(i,2) ) cycle ql22
        end do
        vecsum(1:3) = vec_red(k1)%ei(1:3) + vec_red(k2)%ei(1:3) + vec_red(k3)%ei(1:3)
        if ( (abs(vecsum(1) - nint(vecsum(1))) < toldelta) .and. (abs(vecsum(2) - nint(vecsum(2))) < toldelta) .and. (abs(vecsum(3) - nint(vecsum(3))) < toldelta) ) then
          if ( .not. fl_qallow_written(k1,k2,k3) ) then
            write(10,'(3(i4,1x),3(f8.5,1x))') k1, k2, k3, vecsum(:)
            fl_qallow_written(k1,k2,k3) = .true.
          end if
          ll = ll + 1
          write(11,*) l, l1, l2
        end if
      end do ql22
    end do ql12
  end do ql2

  ! GM, S, S*
  ql3: do l = 1, npGq
    k1 = ceiling(dble(l)/dble(nq))
    j1 = l - (k1-1)*nq
    do i = 1, nskip
      if ( k1 == skipmod(i,1) .and. j1 == skipmod(i,2) ) cycle ql3
    end do
    ql13: do l1 = npGq + npHq + 1, npGq + npHq + npSq
      k2 = ceiling(dble(l1)/dble(nq))
      j2 = l1 - (k2-1)*nq
      do i = 1, nskip
        if ( k2 == skipmod(i,1) .and. j2 == skipmod(i,2) ) cycle ql13
      end do
      ql23: do l2 = npGq + npHq + npSq + 1, npGq + npHq + 2*npSq
        k3 = ceiling(dble(l2)/dble(nq))
        j3 = l2 - (k3-1)*nq
        do i = 1, nskip
          if ( k3 == skipmod(i,1) .and. j3 == skipmod(i,2) ) cycle ql23
        end do
        vecsum(1:3) = vec_red(k1)%ei(1:3) + vec_red(k2)%ei(1:3) + vec_red(k3)%ei(1:3)
        if ( (abs(vecsum(1) - nint(vecsum(1))) < toldelta) .and. (abs(vecsum(2) - nint(vecsum(2))) < toldelta) .and. (abs(vecsum(3) - nint(vecsum(3))) < toldelta) ) then
          if ( .not. fl_qallow_written(k1,k2,k3) ) then
            write(10,'(3(i4,1x),3(f8.5,1x))') k1, k2, k3, vecsum(:)
            fl_qallow_written(k1,k2,k3) = .true.
          end if
          jj3 = l2 - (k3-1)*nq
          k3c = k3 - npS
          l3c = jj3 + (k3c-1)*nq
          if ( l3c < l1 ) cycle ql23
          ll = ll + 1
          write(11,*) l, l1, l2
        end if
      end do ql23
    end do ql13
  end do ql3

  ! (H,H,H) and (H,H,S)
  ql4: do l = npGq + 1, npGq + npHq
    k1 = ceiling(dble(l)/dble(nq))
    j1 = l - (k1-1)*nq
    do i = 1, nskip
      if ( k1 == skipmod(i,1) .and. j1 == skipmod(i,2) ) cycle ql4
    end do
    ql14: do l1 = l, npGq + npHq
      k2 = ceiling(dble(l1)/dble(nq))
      j2 = l1 - (k2-1)*nq
      do i = 1, nskip
        if ( k2 == skipmod(i,1) .and. j2 == skipmod(i,2) ) cycle ql14
      end do
      ql24: do l2 = l1, npGq + npHq + npSq
        k3 = ceiling(dble(l2)/dble(nq))
        j3 = l2 - (k3-1)*nq
        do i = 1, nskip
          if ( k3 == skipmod(i,1) .and. j3 == skipmod(i,2) ) cycle ql24
        end do
        vecsum(1:3) = vec_red(k1)%ei(1:3) + vec_red(k2)%ei(1:3) + vec_red(k3)%ei(1:3)
        if ( (abs(vecsum(1) - nint(vecsum(1))) < toldelta) .and. (abs(vecsum(2) - nint(vecsum(2))) < toldelta) .and. (abs(vecsum(3) - nint(vecsum(3))) < toldelta) ) then
          if ( .not. fl_qallow_written(k1,k2,k3) ) then
            write(10,'(3(i4,1x),3(f8.5,1x))') k1, k2, k3, vecsum(:)
            fl_qallow_written(k1,k2,k3) = .true.
          end if
          ll = ll + 1
          write(11,*) l, l1, l2
        end if
      end do ql24
    end do ql14
  end do ql4
        
  ! (H,S,S) and (H,S,S*)
  ql5: do l = npGq + 1, npGq + npHq
    k1 = ceiling(dble(l)/dble(nq))
    j1 = l - (k1-1)*nq
    do i = 1, nskip
      if ( k1 == skipmod(i,1) .and. j1 == skipmod(i,2) ) cycle ql5
    end do
    ql15: do l1 = npGq + npHq + 1, npGq + npHq + npSq
      k2 = ceiling(dble(l1)/dble(nq))
      j2 = l1 - (k2-1)*nq
      do i = 1, nskip
        if ( k2 == skipmod(i,1) .and. j2 == skipmod(i,2) ) cycle ql15
      end do
      ql25: do l2 = l1, npGq + npHq + 2*npSq
        k3 = ceiling(dble(l2)/dble(nq))
        j3 = l2 - (k3-1)*nq
        do i = 1, nskip
          if ( k3 == skipmod(i,1) .and. j3 == skipmod(i,2) ) cycle ql25
        end do
        vecsum(1:3) = vec_red(k1)%ei(1:3) + vec_red(k2)%ei(1:3) + vec_red(k3)%ei(1:3)
        if ( (abs(vecsum(1) - nint(vecsum(1))) < toldelta) .and. (abs(vecsum(2) - nint(vecsum(2))) < toldelta) .and. (abs(vecsum(3) - nint(vecsum(3))) < toldelta) ) then
          if ( .not. fl_qallow_written(k1,k2,k3) ) then
            write(10,'(3(i4,1x),3(f8.5,1x))') k1, k2, k3, vecsum(:)
            fl_qallow_written(k1,k2,k3) = .true.
          end if
          ll = ll + 1
          write(11,*) l, l1, l2
        end if
      end do ql25
    end do ql15
  end do ql5
        
  ! (S,S,S) and (S,S,S*)
  ql6: do l = npGq + npHq + 1, npGq + npHq + npSq
    k1 = ceiling(dble(l)/dble(nq))
    j1 = l - (k1-1)*nq
    do i = 1, nskip
      if ( k1 == skipmod(i,1) .and. j1 == skipmod(i,2) ) cycle ql6
    end do
    ql16: do l1 = l, npGq + npHq + npSq
      k2 = ceiling(dble(l1)/dble(nq))
      j2 = l1 - (k2-1)*nq
      do i = 1, nskip
        if ( k2 == skipmod(i,1) .and. j2 == skipmod(i,2) ) cycle ql16
      end do
      ql26: do l2 = l1, npGq + npHq + 2*npSq
        k3 = ceiling(dble(l2)/dble(nq))
        j3 = l2 - (k3-1)*nq
        do i = 1, nskip
          if ( k3 == skipmod(i,1) .and. j3 == skipmod(i,2) ) cycle ql26
        end do
        vecsum(1:3) = vec_red(k1)%ei(1:3) + vec_red(k2)%ei(1:3) + vec_red(k3)%ei(1:3)
        if ( (abs(vecsum(1) - nint(vecsum(1))) < toldelta) .and. (abs(vecsum(2) - nint(vecsum(2))) < toldelta) .and. (abs(vecsum(3) - nint(vecsum(3))) < toldelta) ) then
          if ( .not. fl_qallow_written(k1,k2,k3) ) then
            write(10,'(3(i4,1x),3(f8.5,1x))') k1, k2, k3, vecsum(:)
            fl_qallow_written(k1,k2,k3) = .true.
          end if
          ll = ll + 1
          write(11,*) l, l1, l2
        end if
      end do ql26
    end do ql16
  end do ql6
        
  close(10)
  close(11)
  write(*,'(a)') ' done.'

  write(*,'(*(a))') ' Number of unique elements: ',i82a(ll)

  allocate ( delta(ll,3), stat = i )
  if ( i /= 0 ) stop ' ERROR: allocation failed for DELTA'

  open(unit=11,file='tmp_delta.nd',action='read')
  do il = 1, ll
    read(11,*) delta(il,:)
  end do
  close(11)

  write(command,'(a)') 'rm -f tmp_delta.nd'
  call system(command, i)
  if ( i /= 0 ) then
    write(*,'(a)') ' ERROR: unable to remove tmp_delta.nd'
    write(*,*)
    stop
  end if

  return
end subroutine calc_delta


