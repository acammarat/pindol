! pindol version 1.0, Copyright (C) 2023 P. Nicolini, A. Cammarata, M. Dašić
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
! This subroutine applies random distortions on the atomic positions with a specified maximum amplitude
! the distortions can be applied to all atoms, or only following selected phonon mode(s)  
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine generate_distortions
  ! Note on units:
  !   [pos_EC] = Ang
  !   [distort_amp] = Ang
  !   [disp] = Ang
  !   [Q] = amu^1/2 Ang
  use pars, only: error_string, screen_format
  use eigen, only: nq, npunique, nptotq, npGq, npdiffq, npuniqueq
  use input, only: seed_distort, ndistort, distort, distort_amp, distort_mode
  use supercell, only: pos_EC, ncells_tot
  use refconf, only: atoms_UC
  implicit none
  integer :: i, j, k, l, n
  integer, allocatable :: seed(:)
  real(8) :: fact, disp(atoms_UC,3,ncells_tot), Q(nptotq)
  logical :: check(npuniqueq)
  
  ! init RNG
  call random_seed ( size = n )
  allocate ( seed(n), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for temporary array SEED.'
     stop
  end if
  seed(1:n) = seed_distort
  call random_seed ( put = seed(1:n) )
  deallocate ( seed, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for temporary array SEED.'
     stop
  end if

  ! generate random displacements with a uniform distribution
  do i = 1, atoms_UC
     do n = 1, ncells_tot
        do j = 1, 3
           call random_number ( disp(i,j,n) )
           disp(i,j,n) = disp(i,j,n) - 0.5d0
        end do
        ! scale the displacement according to the user-provided amplitude
        fact = sqrt( dot_product(disp(i,:,n),disp(i,:,n)) )
        disp(i,:,n) = disp(i,:,n) * distort_amp / fact
     end do
  end do

  ! in case the user provided the mode list
  ! note that in this case the final displacements can be smaller than the value set by the user
  if ( distort_mode == 'mode' ) then

     ! project the displacements to normal coordinates
     call cartesian_to_normal ( disp, Q, 'dis' )

     ! find the normal coordinates which are in the mode list...
     check(:) = .false.
     do i = 1, ndistort
        do k = 1, npunique
           do j = 1, nq
              if ( distort(i,1) == k .and. distort(i,2) == j ) check((k-1)*nq+j) = .true.
           end do
        end do
     end do

     ! ... and set them to zero (real and imaginary parts for points in the sets S and H)
     do l = 1, npuniqueq
        if ( .not. check(l) ) then
           Q(l) = 0.d0
           if ( l > npGq ) Q(l+npdiffq) = 0.d0
        end if
     end do
     
     ! project back to Cartesian
     call normal_to_cartesian ( Q, disp, 'dis' )

  end if

  ! apply the displacements
  pos_EC(:,:,:) = pos_EC(:,:,:) + disp(:,:,:)
  
  write(*,'(a,'//screen_format//',a)') 'Atomic positions are distorted up to a maximum of ', distort_amp, ' Ang.'
  
  return
end subroutine generate_distortions
