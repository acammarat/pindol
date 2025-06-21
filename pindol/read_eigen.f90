! pindol version 1.0.1, Copyright (C) 2023 P. Nicolini, A. Cammarata, M. Dašić
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
! This subroutine reads eigenvectors and eigenfrequencies from input
! and calculates the corresponding supercell
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
!   - Antonio Cammarata (Czech Technical University in Prague), cammaant@fel.cvut.cz
subroutine read_eigen
  ! Note on units:
  !   [eigen_real] = .
  !   [eigen_imag] = .
  !   [vec] = Ang^-1
  !   [vec_red] = .
  !   [freq] = ps^-1
  !   [freq2*] = fs^-2
  use pars, only: error_string, warning_string, tiny
  use eigen
  use refconf, only: atoms_UC
  use int2str
  use io_units, only: qmatrix_unit, freq_unit
  implicit none
  real(8), parameter :: pi2 = 2.d0 * acos(-1.d0)
  integer :: i, j, k, kk, l, p(3), mcm, lcm, sign_vec
  integer :: nq_tmp, k_tmp
  real(8) :: freq
  character(1) :: comment
  character(6) :: tag1, tag2, tag3, tag4
  character(1000) :: frac(3)
  logical :: skip
  integer, allocatable :: q(:,:), skipmodes(:,:)
  
  write(*,'(a)') 'Reading polarization vectors...'

  ! set number of degrees of freedom for the unit cell
  nq = 3 * atoms_UC 

  ! open frequency file
  open( unit=freq_unit, file='freq.nd', action='read', status='old' )
  read(freq_unit,*) 

  ! read number of "unique" q-points (i.e. Gamma + set H + set S)
  read(freq_unit,*) comment, tag1, npG, tag2, npH, tag3, npS, tag4, nq_tmp
  if ( tag1 /= 'Gamma:' .or. tag2 /= 'H:' .or. tag3 /= 'S:' .or. tag4 /= 'nq:' ) then
     write(0,*) error_string, 'Wrong tags in the line about number of q-points.'
     stop
  end if
  if ( npG /= 0 .and. npG /= 1 ) then
     write(0,*) error_string, 'Wrong number of Gamma-points.'
     stop
  end if
  if ( npH < 0 ) then
     write(0,*) error_string, 'Wrong number of q-points in the set H.'
     stop
  end if
  if ( npS < 0 ) then
     write(0,*) error_string, 'Wrong number of q-points in the set S.'
     stop
  end if
  if ( nq_tmp /= nq ) then
     write(0,*) error_string, 'NQ reported in "qmatrix.nd" differs from what provided in REFCONF.'
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

  ! allocate eigenvectors, frequencies and q-points
  allocate ( eigen_real(npuniqueq,nq), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array EIGEN_REAL.'
     stop
  end if
  allocate ( eigen_imag(npuniqueq,nq), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array EIGEN_IMAG.'
     stop
  end if
  allocate ( vec(npunique,3), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array VEC.'
     stop
  end if
  allocate ( vec_red(npunique,3), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array VEC_RED.'
     stop
  end if
  allocate ( freq2(npuniqueq), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array FREQ2.'
     stop
  end if
  allocate ( freq2_nd(npuniqueq), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array FREQ2_ND.'
     stop
  end if

  ! allocate temporary array
  allocate ( q(npunique,3), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for temporary array Q.'
     stop
  end if
  
  ! reading the number of degrees of freedom that needs to be removed from the computation of temperature and thermostat-related properties
  ! which is 3 for periodic (bulk) system, 5 or 6 for isolated molecules
  read(freq_unit,*) comment, tag1, tag2, ndof_sub

  ! allocate temporary array
  allocate ( skipmodes(ndof_sub,2), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for temporary array SKIPMODES.'
     stop
  end if

  ! read modes to be skipped
  backspace(freq_unit)
  read(freq_unit,*) comment, tag1, tag2, i, (skipmodes(j,1),skipmodes(j,2),j=1,ndof_sub)
  if ( tag1 /= 'Skip' .or. tag2 /= 'modes:' ) then
     write(0,*) error_string, 'Wrong tags in the line about number of acoustic and rotational modes.'
     stop
  end if
  if ( ndof_sub /= 0 .and. ndof_sub /= 3 .and. ndof_sub /= 5 .and. ndof_sub /= 6 ) then
     write(0,*) error_string, 'Wrong number of skip modes.'
     stop
  end if
  read(freq_unit,*) 
  
  ! open eigenvector file
  open( unit=qmatrix_unit, file='qmatrix.nd', action='read', status='old' )
  read(qmatrix_unit,*) !# qpoints v. 1.4

  ! read
  do k = 1, npunique

     ! read eigenvectors
     ! eigen(1:nq,:) -> eigendisplacements for all atoms in UC in all cartesian directions at the first q-point (Gamma or H or S)
     ! ...
     ! eigen(npuniqueq-nq+1:nq,:) -> eigendisplacements for all atoms in UC in all cartesian directions at the last q-point
     do i = 1, nq
        l = (k-1) * nq + i
        read(qmatrix_unit,*) (eigen_real(l,j),eigen_imag(l,j),j=1,nq)
     end do

     ! read frequencies and q-point
     read(freq_unit,*) comment, (vec(k,i),i=1,3) ! Ang^-1, pi2 checked
     read(freq_unit,*) comment, k_tmp, (vec_red(k,i),i=1,3) ! red. units
     if ( k_tmp /= k ) then
        write(0,*) error_string, 'Invalid frequency file.'
        stop
     end if
     do i = 1, nq
        l = (k-1) * nq + i
        read(freq_unit,*) freq ! hard-coded THz
        freq = freq * pi2 * 1.d-3 ! internal units fs^-1, pi2 checked
        freq2_nd(l) = freq * freq
        freq2(l) = sign(1.d0,freq) * freq2_nd(l)
     end do

     ! various checks that the q-points are consistent with the code:
     !
     ! - in general, the code accepts only points that have components in the interval (-1/2,1/2]
     if ( vec_red(k,1) < -0.5d0 .or. vec_red(k,2) < -0.5d0 .or. vec_red(k,3) < -0.5d0 .or. &
          abs(vec_red(k,1)+0.5d0) < tiny .or. abs(vec_red(k,2)+0.5d0) < tiny .or. abs(vec_red(k,3)+0.5d0) < tiny ) then
        write(0,*) error_string, 'A q-point has at least a component equal or smaller than -1/2.'
        stop
     end if
     if ( vec_red(k,1) > 0.5d0 .or. vec_red(k,2) > 0.5d0 .or. vec_red(k,3) > 0.5d0 ) then
        write(0,*) error_string, 'A q-point has at least a component bigger than 1/2.'
        stop
     end if
     !
     ! - if present, the Gamma-point must be the first...
     if ( k == npG .and. ( abs(vec_red(k,1)) > tiny .or. abs(vec_red(k,2)) > tiny .or. abs(vec_red(k,3)) > tiny ) ) then
        write(0,*) error_string, 'Gamma-point should be provided, but the first point is not Gamma.'
        stop
     end if
     !
     ! - ...and not in another position
     if ( k /= npG .and. abs(vec_red(k,1)) < tiny .and. abs(vec_red(k,2)) < tiny .and. abs(vec_red(k,3)) < tiny ) then
        write(0,*) error_string, 'Gamma-point (if present) must be provided as the first point.'
        stop
     end if
     !
     ! - if present, q-points in the set H must have at least one component equal to 1/2
     if ( k /= npG .and. k <= npG+npH .and. abs(vec_red(k,1)-0.5d0) > tiny .and. &
          abs(vec_red(k,2)-0.5d0) > tiny .and. abs(vec_red(k,3)-0.5d0) > tiny ) then
        write(0,*) error_string, 'A q-point in the set H has no component equal to 1/2.'
        stop
     end if
     !
     ! for the set H:
     if ( k /= npG .and. k <= npG+npH ) then
        ! - check that q-points are not duplicated
        do kk = npG+1, k-1
           if ( ( abs(abs(vec_red(k,1)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,1)) - 0.5d0) < tiny ) .and. &
              ( abs(vec_red(k,2)-vec_red(kk,2)) < tiny ) .and. ( abs(vec_red(k,3)-vec_red(kk,3)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H is duplicated.'
              stop
           else if ( ( abs(abs(vec_red(k,2)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,2)) - 0.5d0) < tiny ) .and. &
              ( abs(vec_red(k,1)-vec_red(kk,1)) < tiny ) .and. ( abs(vec_red(k,3)-vec_red(kk,3)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H is duplicated.'
              stop
           else if ( ( abs(abs(vec_red(k,3)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,3)) - 0.5d0) < tiny ) .and. &
              ( abs(vec_red(k,2)-vec_red(kk,2)) < tiny ) .and. ( abs(vec_red(k,1)-vec_red(kk,1)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H is duplicated.'
              stop
           else if ( ( abs(abs(vec_red(k,1)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,1)) - 0.5d0) < tiny ) .and. &
              ( abs(abs(vec_red(k,2)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,2)) - 0.5d0) < tiny ) .and. &
              ( abs(vec_red(k,3)-vec_red(kk,3)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H is duplicated.'
              stop
           else if ( ( abs(abs(vec_red(k,1)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,1)) - 0.5d0) < tiny ) .and. &
               ( abs(abs(vec_red(k,3)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,3)) - 0.5d0) < tiny ) .and. &
               ( abs(vec_red(k,2)-vec_red(kk,2)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H is duplicated.'
              stop
           else if ( ( abs(abs(vec_red(k,2)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,2)) - 0.5d0) < tiny ) .and. &
              ( abs(abs(vec_red(k,3)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,3)) - 0.5d0) < tiny ) .and. &
              ( abs(vec_red(k,1)-vec_red(kk,1)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H is duplicated.'
           else if ( ( abs(abs(vec_red(k,1)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,1)) - 0.5d0) < tiny ) .and. &
              ( abs(abs(vec_red(k,2)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,2)) - 0.5d0) < tiny ) .and. &
              ( abs(abs(vec_red(k,3)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,3)) - 0.5d0) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H is duplicated.'
              stop
           end if
        end do
        ! - check that q-points are unique, i.e. no complex-conjugates
        do kk = npG+1, k-1
           ! if one of the component is 1/2 and the other corresponding ones are opposite, then the two vectors are a conjugated couple;
           ! for example, q1(1/2,1/4,0) is conjugated couple with q2=-q1=(-1/2,-1/4,0). However, q3=(1/2,-1/4,0) is equivalent to q2 because of
           ! the periodicity of the Brillouin zone, as 1/2=-1/2, so q1 and q3 are conjugated couple too.
           if ( ( abs(abs(vec_red(k,1)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,1)) - 0.5d0) < tiny ) .and. &
                ( abs(vec_red(k,2)+vec_red(kk,2)) < tiny ) .and. ( abs(vec_red(k,3)+vec_red(kk,3)) < tiny ) ) then
                write(0,*) error_string, 'A q-point in the set H and its complex-conjugated are provided.'
                stop
           else if ( ( abs(abs(vec_red(k,2)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,2)) - 0.5d0) < tiny ) .and. &
              ( abs(vec_red(k,1)+vec_red(kk,1)) < tiny ) .and. ( abs(vec_red(k,3)+vec_red(kk,3)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H and its complex-conjugated are provided.'
              stop
           else if ( ( abs(abs(vec_red(k,3)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,3)) - 0.5d0) < tiny ) .and. &
              ( abs(vec_red(k,2)+vec_red(kk,2)) < tiny ) .and. ( abs(vec_red(k,1)+vec_red(kk,1)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H and its complex-conjugated are provided.'
              stop

           else if ( ( abs(abs(vec_red(k,1)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,1)) - 0.5d0) < tiny ) .and. &
              ( abs(abs(vec_red(k,2)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,2)) - 0.5d0) < tiny ) .and. &
              ( abs(vec_red(k,3)+vec_red(kk,3)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H and its complex-conjugated are provided.'
              stop
           else if ( ( abs(abs(vec_red(k,1)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,1)) - 0.5d0) < tiny ) .and. &
              ( abs(abs(vec_red(k,3)) - 0.5d0) < tiny ) .and. ( abs(abs(vec_red(kk,3)) - 0.5d0) < tiny ) .and. &
              ( abs(vec_red(k,2)+vec_red(kk,2)) < tiny ) ) then
              write(0,*) error_string, 'A q-point in the set H and its complex-conjugated are provided.'
              stop
           else if ( abs(vec_red(k,1)+vec_red(kk,1)) < tiny .and. abs(vec_red(k,2)+vec_red(kk,2)) < tiny .and. &
              abs(vec_red(k,3)+vec_red(kk,3)) < tiny ) then
              write(0,*) error_string, 'A q-point in the set H and its complex-conjugated are provided.'
              stop
           end if
        end do
     end if
     !
     ! for the set S:
     if ( k > npG+npH ) then
        ! - check that q-points are not duplicated
        do kk = npG+1, k-1
           if ( abs(vec_red(k,1)-vec_red(kk,1)) < tiny .and. abs(vec_red(k,2)-vec_red(kk,2)) < tiny .and. &
                abs(vec_red(k,3)-vec_red(kk,3)) < tiny ) then
              write(0,*) error_string, 'A q-point in the set S is duplicated.'
              stop
           end if
        end do
        ! - check that q-points are unique, i.e. no complex-conjugates
        do kk = npG+1, k-1
           if ( abs(vec_red(k,1)+vec_red(kk,1)) < tiny .and. abs(vec_red(k,2)+vec_red(kk,2)) < tiny .and. &
                abs(vec_red(k,3)+vec_red(kk,3)) < tiny ) then
              write(0,*) error_string, 'A q-point in the set S and its complex-conjugated are provided.'
              stop
           end if
        end do
     end if
     
     ! print the q-point as a fraction
     do i = 1, 3
        if ( abs(vec_red(k,i)) > tiny ) then
           sign_vec = 1
           if ( vec_red(k,i) < 0.d0 ) sign_vec = -1
           call real_to_rational ( dble(sign_vec)*vec_red(k,i), p(i), q(k,i) )
           p(i) = sign_vec*p(i)
        else
           p(i) = 0
           q(k,i) = 1
        end if
        call print_fraction ( p(i), q(k,i), frac(i) )
      end do
      write(*,'(2x,8a)') i2a(k), ' q-point found: (', trim(frac(1)), ',', trim(frac(2)), ',', trim(frac(3)), ')'
      
  end do

  close( qmatrix_unit )
  close( freq_unit )
  
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
  write(*,'(6a)') '  Supercell size: ', i2a(ncells(1)),' x ', i2a(ncells(2)),' x ', i2a(ncells(3))
  
  ! printout unstable modes (if any)
  do k = 1, npunique
     do i = 1, nq
        skip = .false.
        do j = 1, ndof_sub
           if ( skipmodes(j,1) == k .and. skipmodes(j,2) == i ) skip = .true.
        end do
        l = (k-1) * nq + i
        if ( .not. skip .and. freq2_nd(l) < 0.d0 ) write(*,'(6a,g12.5,a)') warning_string, 'Imaginary frequency found for point ', trim(i2a(k)), ' mode ', trim(i2a(i)), ' = ', sqrt(freq2(l))/pi2*1.d3, ' THz'
     end do
  end do
  
  ! free temporary array
  deallocate ( q, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for temporary array Q.'
     stop
  end if
  deallocate ( skipmodes, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for temporary array SKIPMODES.'
     stop
  end if

  write(*,'(a)') 'Reading polarization vectors...DONE.'

  return
end subroutine read_eigen


! This subroutine returns a string in order to print the q-point as a fraction
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
!   - Miljan Dašić (Czech Technical University in Prague, Czech Republic & Institute of Physics Belgrade, Serbia), miljan.dasic@scl.rs
subroutine print_fraction ( n, d, string )
  use int2str
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


! This subroutine finds the best integer p and q such that p/q is close to a given real number x
! Contributors:
!   - Miljan Dašić (Czech Technical University in Prague, Czech Republic & Institute of Physics Belgrade, Serbia), miljan.dasic@scl.rs
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


! This subroutine returns the least common multiplier (lcm) of a pair of integers (a,b)
! Contributors:
!   - Miljan Dašić (Czech Technical University in Prague, Czech Republic & Institute of Physics Belgrade, Serbia), miljan.dasic@scl.rs
integer function lcm ( a, b )
  implicit none
  integer, intent(in) :: a, b
  integer :: gcd
  
  lcm = a * b / gcd(a,b)
  
  return
end function lcm


! This subroutine returns the greatest common divisor (gcd) of a pair of integers (a,b)
! Contributors:
!   - Miljan Dašić (Czech Technical University in Prague, Czech Republic & Institute of Physics Belgrade, Serbia), miljan.dasic@scl.rs
integer function gcd ( a, b )
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
