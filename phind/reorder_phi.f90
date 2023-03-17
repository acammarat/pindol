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

subroutine reorder_phi
  use integer_to_string
  use var, only: nq, ll, phi, version, vec_red, npSq

  implicit none
  integer :: i, j, k1, k2, k3
  integer(8) :: il, jl(24), l, l1, l2
  character(256), parameter :: format_phi='(3(a,1x),G22.14,1x,G22.14)'
  character(256) :: command

  write(*,'(a)',advance='no') ' Reordering phi...'

  do i = 1, 24
    open(unit=20+i,file='tmp_'//i2a(i)//'.nd',action='write')
    jl(i) = 0
  end do

  do il = 1, ll

    l  = phi(il)%kk(1)
    l1 = phi(il)%kk(2)
    l2 = phi(il)%kk(3)
    k1 = ceiling(dble(l)/dble(nq))
    k2 = ceiling(dble(l1)/dble(nq))
    k3 = ceiling(dble(l2)/dble(nq))

    ! case I) -> l in GAMMA, l1 in GAMMA, l2 in GAMMA
    if ( vec_red(k1)%set == 'G' .and. vec_red(k2)%set == 'G' .and. vec_red(k3)%set == 'G' ) then 
      ! --- sub-case A) -> l == l1 == l2 (this sub-case is unique, i.e. the permutation property does not apply)
      if ( l == l1 .and. l1 == l2 ) then 
        jl(1) = jl(1) + 1
        write(21,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case B) -> l == l1 /= l2 (in this sub-case we need to update sum_phi for l and l2)
      else if ( l == l1 .and. l1 /= l2 ) then
        jl(2) = jl(2) + 1
        write(22,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case C) -> l /= l1 == l2 (need to update sum_phi for l and l1; in the latter, one permutation is further possible)
      else if ( l /= l1 .and. l1 == l2 ) then
        jl(3) = jl(3) + 1
        write(23,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case D) -> l /= l1 /= l2 (need to update sum_phi for all l, l1 and l2)
      else if ( l /= l1 .and. l1 /= l2 ) then
        jl(4) = jl(4) + 1
        write(24,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      end if
  
    ! case VI) -> l in Gamma, l1 in H(q), l2 in H(q)
    ! this covers also:
    ! case XVIII) -> l in H(q), l1 in Gamma, l2 in H(q) via permutation
    ! case XXI) -> l in H(q), l1 in H(q), l2 in Gamma via permutation
    else if ( vec_red(k1)%set == 'G' .and. vec_red(k2)%set == 'H' .and. vec_red(k3)%set == 'H' ) then
      ! --- sub-case A) -> l1 == l2
      if ( l1 == l2 ) then
        jl(5) = jl(5) + 1
        write(25,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case B) -> l1 /= l2
      else if ( l1 /= l2 ) then
        jl(6) = jl(6) + 1
        write(26,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      end if

    ! case XII) -> l in Gamma, l1 in S(q), l2 in S*(q)
    ! this covers also:
    ! case XV) -> l in Gamma, l1 in S*(q), l2 in S(q) via permutation AND (possibly) complex-conjugation
    ! case XXXVI) -> l in S(q), l1 in Gamma, l2 in S*(q) via permutation AND (possibly) complex-conjugation
    ! case XLV) -> l in S(q), l1 in S*(q), l2 in Gamma via permutation AND (possibly) complex-conjugation
    ! --- sub-case A) -> l1 == l2-nrq (i.e. l1 and l2 are complex-conjugates)
    else if ( vec_red(k1)%set == 'G' .and. vec_red(k2)%set == 'S' .and. vec_red(k3)%set == 'C' ) then
      ! --- sub-case A) -> l1 == l2-nrq (i.e. l1 and l2 are complex-conjugates)
      if ( l1 == l2-npSq ) then
        jl(7) = jl(7) + 1
        write(27,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case B) -> l1 /= l2-nrq
      else if ( l1 /= l2-npSq ) then
        jl(8) = jl(8) + 1
        write(28,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      end if

    ! case XXII) -> l in H(q), l1 in H(q), l2 in H(q)
    else if ( vec_red(k1)%set == 'H' .and. vec_red(k2)%set == 'H' .and. vec_red(k3)%set == 'H' ) then
      ! --- sub-case A) -> l == l1 == l2 (this sub-case is unique, i.e. the permutation property does not apply)
      if ( l == l1 .and. l1 == l2 ) then
        jl(9) = jl(9) + 1
        write(29,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case B) -> l == l1 /= l2 (in this sub-case we need to update sum_phi for l and l2)
      else if ( l == l1 .and. l1 /= l2 ) then
        jl(10) = jl(10) + 1
        write(30,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case C) -> l /= l1 == l2 (need to update sum_phi for l and l1; in the latter, one permutation is further possible)
      else if ( l /= l1 .and. l1 == l2 ) then
        jl(11) = jl(11) + 1
        write(31,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case D) -> l /= l1 /= l2 (need to update sum_phi for all l, l1 and l2)
      else if ( l /= l1 .and. l1 /= l2 ) then
        jl(12) = jl(12) + 1
        write(32,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      end if

    ! case XXIII) -> l in H(q), l1 in H(q), l2 in S(q)
    ! this covers also:
    ! case XXIV) -> l in H(q), l1 in H(q), l2 in S*(q) via complex-conjugation 
    ! case XXVI) -> l in H(q), l1 in S(q), l2 in H(q) via permutation 
    ! case XXX) -> l in H(q), l1 in S*(q), l2 in H(q) via permutation + complex-conjugation 
    ! case XXXVIII) -> l in S(q), l1 in H(q), l2 in H(q) via permutation 
    else if ( vec_red(k1)%set == 'H' .and. vec_red(k2)%set == 'H' .and. vec_red(k3)%set == 'S' ) then
      ! --- sub-case A) -> l == l1
      if ( l == l1 ) then
        jl(13) = jl(13) + 1
        write(33,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case B) -> l /= l1
      else if ( l /= l1 ) then
        jl(14) = jl(14) + 1
        write(34,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      end if

    ! case XXVII) -> l in H(q), l1 in S(q), l2 in S(q)
    ! this covers also:
    ! case XXXII) -> l in H(q), l1 in S*(q), l2 in S*(q) via complex-conjugation
    ! case XXXIX) -> l in S(q), l1 in H(q), l2 in S(q) via permutation
    ! case XLII) -> l in S(q), l1 in S(q), l2 in H(q) via permutation
    else if ( vec_red(k1)%set == 'H' .and. vec_red(k2)%set == 'S' .and. vec_red(k3)%set == 'S' ) then
      ! --- sub-case A) -> l1 == l2
      if ( l1 == l2 ) then
        jl(15) = jl(15) + 1
        write(35,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case B) -> l1 /= l2
      else if ( l1 /= l2 ) then
        jl(16) = jl(16) + 1
        write(36,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      end if

    ! case XXVIII) -> l in H(q), l1 in S(q), l2 in S*(q)
    ! this covers also:
    ! case XXXI) -> l in H(q), l1 in S*(q), l2 in S(q) via permutation AND (possibly) complex-conjugation
    ! case XL) -> l in S(q), l1 in H(q), l2 in S*(q) via permutation AND (possibly) complex-conjugation
    ! case XLVI) -> l in S(q), l1 in S*(q), l2 in H(q) via permutation AND (possibly) complex-conjugation
    else if ( vec_red(k1)%set == 'H' .and. vec_red(k2)%set == 'S' .and. vec_red(k3)%set == 'C' ) then
      ! --- sub-case A) -> l1 == l2-nrq (i.e. l1 and l2 are complex-conjugates)
      if ( l1 == l2-npSq ) then
        jl(17) = jl(17) + 1
        write(37,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case B) -> l1 /= l2-nrq
      else if ( l1 /= l2-npSq ) then
        jl(18) = jl(18) + 1
        write(38,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      end if

    ! case XLIII) -> l in S(q), l1 in S(q), l2 in S(q)
    else if ( vec_red(k1)%set == 'S' .and. vec_red(k2)%set == 'S' .and. vec_red(k3)%set == 'S' ) then
      ! --- sub-case A) -> l == l1 == l2 (this sub-case is unique, i.e. the permutation property does not apply)
      if ( l == l1 .and. l1 == l2 ) then
        jl(19) = jl(19) + 1
        write(39,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case B) -> l == l1 /= l2 (in this sub-case we need to update sum_phi for l and l2)
      else if ( l == l1 .and. l1 /= l2 ) then
        jl(20) = jl(20) + 1
        write(40,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case C) -> l /= l1 == l2 (need to update sum_phi for l and l1; in the latter, one permutation is further possible)
      else if ( l /= l1 .and. l1 == l2 ) then
        jl(21) = jl(21) + 1
        write(41,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case D) -> l /= l1 /= l2 (need to update sum_phi for all l, l1 and l2)
      else if ( l /= l1 .and. l1 /= l2 ) then
        jl(22) = jl(22) + 1
        write(42,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      end if

    ! case XLIV) -> l in S(q), l1 in S(q), l2 in S*(q)
    ! this covers also:
    ! case XLVII) -> l in S(q), l1 in S*(q), l2 in S(q) via permutation
    ! case XLVIII) -> l in S(q), l1 in S*(q), l2 in S*(q) via permutation + complex-conjugation 
    else if ( vec_red(k1)%set == 'S' .and. vec_red(k2)%set == 'S' .and. vec_red(k3)%set == 'C' ) then
      ! --- sub-case A) -> l == l1
      if ( l == l1 ) then
        jl(23) = jl(23) + 1
        write(43,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      ! --- sub-case B) -> l /= l1
      else if ( l /= l1 ) then
        jl(24) = jl(24) + 1
        write(44,format_phi) i82a(l), i82a(l1), i82a(l2), phi(il)%el(:) ! [phi%el]=[uma^(-3/2).eV.Ang^-3]
      end if

    end if

  end do

  write(*,'(a)') ' done.'

  do i = 1, 24
    close(unit=20+i)
  end do

  write(*,'(a)', advance='no') ' Calculating md5sum for qmatrix.nd and freq.nd...'

  write(command,'(a)') "md5sum qmatrix.nd|awk '{print $1}' > tmp.nd"
  call system(command, i)
  if ( i /= 0 ) then
    write(*,'(a)') ' ERROR: md5sum exited with non-zero status.'
    write(*,*)
    stop
  end if

  write(command,'(a)') "md5sum freq.nd|awk '{print $1}' >> tmp.nd"
  call system(command, i)
  if ( i /= 0 ) then
    write(*,'(a)') ' ERROR: md5sum exited with non-zero status.'
    write(*,*)
    stop
  end if

  write(*,'(a)') ' done.'

  write(*,'(a)',advance='no') ' Writing reordered phi in file phi.nd ...'

  open(unit=10,file='phi.nd')

  open(unit=20,file='tmp.nd')
  write(10,'(*(a))') '# phind v. ', version
  read(20,*) command
  write(10,'(*(a))') '# md5sum qmatrix.nd ', trim(command)
  read(20,*) command
  write(10,'(*(a))') '# md5sum freq.nd ', trim(command)
  close(20)
  write(command,'(a)') 'rm -f tmp.nd'
  call system(command, i)
  if ( i /= 0 ) then
    write(*,'(a)') ' ERROR: unable to remove tmp.nd'
    write(*,*)
    stop
  end if

  write(10,'(a)') i82a(ll)
  do i = 1, 23
    write(10,'(a,1x)',advance='no') i82a(jl(i))
  end do
  write(10,'(a)',advance='no') i82a(jl(24))

  close(10)

  write(command,'(a)') 'cat tmp_1.nd tmp_2.nd tmp_3.nd tmp_4.nd tmp_5.nd tmp_6.nd tmp_7.nd tmp_8.nd tmp_9.nd tmp_10.nd tmp_11.nd tmp_12.nd tmp_13.nd tmp_14.nd tmp_15.nd tmp_16.nd tmp_17.nd tmp_18.nd tmp_19.nd tmp_20.nd tmp_21.nd tmp_22.nd tmp_23.nd tmp_24.nd >> phi.nd'
  call system(command, j)
  if ( j /= 0 ) then
    write(*,'(a)') ' ERROR: unable to cat tmp_*.nd files into phi.nd'
    write(*,*)
    stop
  end if

  do i = 1, 24
    write(command,'(a)') 'rm -f tmp_'//i2a(i)//'.nd'
    call system(command, j)
    if ( j /= 0 ) then
      write(*,'(a)') ' ERROR: unable to remove tmp_'//i2a(i)//'.nd'
      write(*,*)
      stop
    end if
  end do

  write(*,'(a)') ' done.'

  deallocate ( phi, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for PHI'

  return
end subroutine reorder_phi


