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
! This subroutine initializes the atomic velocities to values sampled randomly from a Gaussian distribution
! and then rescale them to match the user-given temperature
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine initialize_velocities
  ! Note on units:
  !   [vel_EC] = Ang fs^-1
  !   [mass_EC] = amu
  !   [mom_tot] = amu Ang fs^-1
  !   [mass_tot] = amu
  !   [ekin2] = amu Ang^2 fs^-2
  !   [kB] = J K^-1
  !   [init_temperature] = K 
  use pars, only: error_string, screen_format, kB, amu_to_kg
  use supercell, only: vel_EC, mass_EC, ncells_tot
  use eigen, only: ndof_sub
  use refconf, only: atoms_UC
  use input, only: seed_initvel, init_temperature
  implicit none
  integer :: i, n, j, k
  integer, allocatable :: seed(:)
  real(8) :: mom_tot(3), r, mass_tot, fact, ekin2
  
  ! init RNG
  call random_seed ( size = n )
  allocate ( seed(n), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for temporary array SEED.'
     stop
  end if
  seed(1:n) = seed_initvel
  call random_seed ( put = seed(1:n) )
  deallocate ( seed, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for temporary array SEED.'
     stop
  end if
  
  ! generate Gaussian distribution and calculate total linear momentum
  mom_tot(:) = 0.d0
  mass_tot = 0.d0
  do i = 1, atoms_UC
     do n = 1, ncells_tot
        mass_tot = mass_tot + mass_EC(i,n)
        do j = 1, 3
           ! see Sec. 4.2 of Oleg M. Braun, Computer Modeling in Physics
           ! downloaded at http://www.iop.kiev.ua/~obraun/book_md/book_md_0.pdf
           vel_EC(i,j,n) = 0.d0
           do k = 1, 12
              call random_number ( r )
              vel_EC(i,j,n) = vel_EC(i,j,n) + r
           end do
           vel_EC(i,j,n) = vel_EC(i,j,n) - 6.d0
           mom_tot(j) = mom_tot(j) + mass_EC(i,n) * vel_EC(i,j,n)
        end do
     end do
  end do

  ! shift the distribution
  do j = 1, 3
     vel_EC(:,j,:) = vel_EC(:,j,:) - mom_tot(j) / mass_tot
  end do

  ! compute the kinetic energy
  ekin2 = 0.d0
  do i = 1, atoms_UC
     do n = 1, ncells_tot
        do j = 1, 3
           ekin2 = ekin2 + mass_EC(i,n) * vel_EC(i,j,n)*vel_EC(i,j,n)
        end do
     end do
  end do
  
  ! scale velocities to match the given temperature
  ! note that Boltzmann constant must be used in amu Ang^2 fs^2 K^-1 in order to obtain velocities in Ang/fs
  fact = sqrt( dble(3*atoms_UC*ncells_tot-ndof_sub) * kB / amu_to_kg * 1.d-10 * init_temperature / ekin2 )
  vel_EC(:,:,:) = vel_EC(:,:,:) * fact

  write(*,'(a,'//screen_format//',a)') 'Atomic velocities initialized at ', init_temperature, ' K.'
  
  return
end subroutine initialize_velocities
