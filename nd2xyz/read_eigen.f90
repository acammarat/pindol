!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nd2xyz. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Converts an ND trajectory into xyz format
! 
! If used for production, you should cite
! A. Cammarata, M. Dasic, P. Nicolini, J. Chem. Phys. 161, 084111 (2024)
! https://doi.org/10.1063/5.0224108
!
!    This file is part of nd2xyz.
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

! this subroutine reads eigenvectors and q-points from input
! and calculates the corresponding supercell
! note on units:
!   [eigen] = .
!   [vec] = Ang^-1
!   [vec_red] = .
!   [freq] = ps^-1
!   [freq2] = fs^-2
subroutine read_eigen
  use var, only: error_string, eigen, vec, vec_red,  &
       ncells, nq_chk, atoms_UC, ncells_tot, pos_eq_EC, at_EC, &
       side_eq_EC, side_eq_UC, ncells, at_UC, pos_eq_UC, atoms_EC, toldelta, &
       npG, npH, npS, nq, npGq, npHq, npSq, npdiffq, npunique, nptot, npuniqueq, nptotq
  use functions, only: i2a
  implicit none
  integer :: i, j, k, kk, l, p(3), mcm, lcm, sign
  integer :: i1, i2, i3, n
  integer, allocatable :: q(:,:)
  real(8) :: xt, yt, zt
  real(8), allocatable :: tmp_r(:), tmp_i(:)
  character(1) :: comment
  character(1000) :: frac(3)
  
  write(*,'(a)') ' Reading polarization vectors:'

  ! open frequency file
  open( unit=20, file='freq.nd', action='read', status='old' )
  read(20,*)
  read(20,*) comment, comment, npG, comment, npH, comment, npS, comment, nq
  read(20,*)
  read(20,*)

  if ( nq_chk /= nq ) then
    write(*,'(a)') error_string, ' The number of modes read in qmatrix.nd is not consistent with the one calculated from the geometry file.'
    write(*,*)
    stop
  end if

  ! set auxiliary numbers of degrees of freedom for the normal dynamics
  npGq = npG * nq
  npHq = npH * nq
  npSq = npS * nq
  npdiffq = npHq + npSq

  ! set number of unique q-points (i.e. Gamma + set H + set S)
  npunique = npG + npH + npS

  ! set total number of q-points (i.e. Gamma + set H + set S + set S*)
  nptot = npG + npH + 2*npS

  ! set number of degrees of freedom for the normal dynamics (i.e. the number of dynamical variables, note that some of them are complex!)
  npuniqueq = npGq + npHq + npSq

  ! set extended number of degrees for the normal dynamics (i.e. the actual number of real values stored in the configurational part)
  nptotq = npGq + 2*npHq + 2*npSq

  allocate ( tmp_r(nq), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'ALLOC tmp_r'
     stop
  end if

  allocate ( tmp_i(nq), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'ALLOC tmp_i'
     stop
  end if

  ! allocate eigenvectors, frequencies and q-points
  allocate ( eigen(npuniqueq,nq), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'ALLOC EIGEN'
     stop
  end if
  allocate ( vec(npunique,3), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'ALLOC VEC'
     stop
  end if
  allocate ( vec_red(npunique,3), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'ALLOC VEC_RED'
     stop
  end if

  ! allocate temporary array
  allocate ( q(npunique,3), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'ALLOC Q'
     stop
  end if

  ! open eigenvector file
  open( unit=10, file='qmatrix.nd', action='read', status='old' )
  read(10,*)

  ! read
  do k = 1, npunique

     ! read eigenvectors
     ! eigen(1:nq,:) -> eigendisplacements for all atoms in UC in all cartesian directions at Gamma
     ! ...
     ! eigen((nunique-1)*nq+1:nq,:) -> eigendisplacements for all atoms in UC in all cartesian directions at the last q-point
     do i = 1, nq
        l = (k-1) * nq + i
        read(10,*) (tmp_r(j),tmp_i(j),j=1,nq)
        eigen(l,:) = dcmplx( tmp_r(:), tmp_i(:) )
     end do

     ! read only q-vector components
     read(20,*) comment, (vec(k,i),i=1,3) ! Ang^-1 ! pi2 checked!
     read(20,*) comment, comment, (vec_red(k,i),i=1,3) ! red. units
     do i = 1, nq
        read(20,*) ! frequencies are not used
     end do

     ! print the q-point as a fraction
     do i = 1, 3
        if ( abs(vec_red(k,i)) > tiny(1.d0) ) then
           sign = 1
           if ( vec_red(k,i) < 0.d0 ) sign = -1
           call real_to_rational ( dble(sign)*vec_red(k,i), p(i), q(k,i) )
           p(i) = sign*p(i)
        else
           p(i) = 0
           q(k,i) = 1
        end if
        call print_fraction ( p(i), q(k,i), frac(i) )
      end do
      write(*,'(2x,8a)') i2a(k), ' q-point found: (', trim(frac(1)), ',', trim(frac(2)), ',', trim(frac(3)), ')'
  end do

  close(10)
  close(20)

  ! check if the provided q-point set is S, that is,
  ! if there is no complex conjugated couple

  do k = 1 , npunique-1
    do i = k+1, npunique
      if ( ( abs(vec_red(k,1)+vec_red(i,1)) < toldelta ) .and. ( abs(vec_red(k,2)+vec_red(i,2)) < toldelta ) .and. &
           ( abs(vec_red(k,3)+vec_red(i,3)) < toldelta ) ) then
        write(*,'(a16,i4,a5,i4,a24)') ' ERROR: vectors ', k,' and ', i, ' are a conjugated couple'
        stop
      end if
    end do
  end do
  
  ! calculate sizes of the supercell compatible with the q-points
  do i = 1, 3    
     if (npunique == 1) then
        mcm = q(1,i)
     else        
        mcm = 1
        do k = 1, npunique-1
           do kk = k+1, npunique
              mcm = lcm(mcm, lcm(q(k,i), q(kk,i)))
           end do
        end do
     end if         
     ncells(i) = abs(mcm)
  end do
  write(*,'(6a)') '  Direct supercell size: ', i2a(ncells(1)),' x ', i2a(ncells(2)),' x ', i2a(ncells(3))

  ! free temporary array
  deallocate ( q, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'DEALLOC Q'
     stop
  end if
  
  ncells_tot = ncells(1)*ncells(2)*ncells(3)
  atoms_EC = atoms_UC * ncells_tot

  allocate ( pos_eq_EC(atoms_UC,3,ncells_tot), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for POS_EQ_EC'
  allocate ( at_EC(atoms_UC,ncells_tot), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for AT_EC'

  do i = 1, 3
    do k = 1, 3
      side_eq_EC(i,k) = side_eq_UC(i,k)*dble(ncells(i))
    end do
  end do

  n = 0
  do i3 = 1, ncells(3)
     do i2 = 1, ncells(2)
        do i1 = 1, ncells(1)
          n = n + 1
          xt = dble(i1-1)*side_eq_UC(1,1) + dble(i2-1)*side_eq_UC(2,1) + dble(i3-1)*side_eq_UC(3,1)
          yt = dble(i1-1)*side_eq_UC(1,2) + dble(i2-1)*side_eq_UC(2,2) + dble(i3-1)*side_eq_UC(3,2)
          zt = dble(i1-1)*side_eq_UC(1,3) + dble(i2-1)*side_eq_UC(2,3) + dble(i3-1)*side_eq_UC(3,3)
          pos_eq_EC(:,1,n) = pos_eq_UC(:,1) + xt
          pos_eq_EC(:,2,n) = pos_eq_UC(:,2) + yt
          pos_eq_EC(:,3,n) = pos_eq_UC(:,3) + zt
          at_EC(:,n) = at_UC(:)
        end do
     end do
  end do

  
  write(*,'(a)') ' done.'

  return
end subroutine read_eigen


! print the q-point in string format
subroutine print_fraction ( n, d, string )
  use functions, only: i2a
  implicit none
  integer, intent(in) :: n, d
  character(1000), intent(out) :: string

  if ( n == 0 ) then
     string = '0'
  else
     string = i2a(n)//'/'//i2a(d)
  end if

  return
end subroutine print_fraction


! determines p/q for a given real number x
subroutine real_to_rational ( x, p, q )
  implicit none
  real(8), intent(in) :: x
  integer, intent(out) :: p, q
  integer :: f, gcd
  real(8) :: r, e, best 
  
  p = 1 
  q = 1         
  best = x * 6.d0    

  do 
     r = dble(p) / dble(q)                
     e = x - r                  
     if ( abs(e) <= best ) then 
        best = abs(e) * 0.125d0             
        f = gcd(p,q)                    
        if ( abs(e) < 0.000001d0 ) exit                
     end if
     if ( e > 0.d0 ) then 
        p = p + ceiling( e * q )    
     else if ( e < 0.d0 ) then    
        q = q + 1                       
     end if
  end do

  return        
end subroutine real_to_rational


! returns the least common multiplier (lcm) of a pair of integers (a,b)
integer function lcm(a,b)
  implicit none
  integer, intent(in) :: a, b
  integer :: gcd
  
  lcm = a * b / gcd(a,b)
  
  return
end function lcm


! returns the greatest common divisor (gcd) of a pair of integers (a,b)
integer function gcd(a,b)
  implicit none
  integer, intent(in) :: a, b
  integer :: aa, bb, t
  
  aa = a
  bb = b
  do while ( bb /= 0 )
     t = bb
     bb = mod(aa,bb)
     aa = t
  end do
  gcd = abs(aa)
  
  return
end function gcd
