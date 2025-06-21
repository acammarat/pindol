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
! This subroutine define the supercell according to the number of replicas found in the read_eigen subroutine
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine setup_supercell
  ! Note on units:
  !   [side*,pos*] = Ang
  !   [mass*] = amu
  !   [vel*] = Ang fs^-1
  use pars, only: error_string, warning_string
  use supercell
  use eigen, only: ncells, nptot
  use refconf, only: side_UC, atoms_UC, pos_eq_UC, mass_UC, at_UC
  use int2str
  implicit none
  integer :: i, i1, i2, i3, k, n
  real(8) :: xt, yt, zt
  
  write(*,'(a)') 'Setting up the supercell...'

  ! calculate supercell size
  ncells_tot = ncells(1) * ncells(2) * ncells(3)
  atoms_EC = atoms_UC * ncells_tot

  ! check that the number of replicas in the supercell is equal to the total number of q-points
  if ( ncells_tot /= nptot ) then
     write(0,*) warning_string, 'The number of replicas in the supercell does not match the total number of q-points.'
     write(0,*) warning_string, 'The set of provided q-points is likely to be not complete.'
  end if
    
  ! calculate supercell box
  do k = 1, 3
     side_EC(:,k) = side_UC(:,k) * dble(ncells(:)) 
  end do

  ! allocate supercell arrays
  allocate ( pos_eq_EC(atoms_UC,3,ncells_tot), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array POS_EQ_EC.'
     stop
  end if
  allocate ( mass_EC(atoms_UC,ncells_tot), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array MASS_EC.'
     stop
  end if
  allocate ( at_EC(atoms_UC,ncells_tot), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array AT_EC.'
     stop
  end if
  allocate ( pos_EC(atoms_UC,3,ncells_tot), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array POS_EC.'
     stop
  end if
  allocate ( vel_EC(atoms_UC,3,ncells_tot), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array VEL_EC.'
     stop
  end if

  ! generate supercell atomic positions according to the reference configuration
  do i = 1, atoms_UC
     n = 0
     do i3 = 1, ncells(3)
        do i2 = 1, ncells(2)
           do i1 = 1, ncells(1)
              n = n + 1
              xt = dble(i1-1)*side_UC(1,1) + dble(i2-1)*side_UC(2,1) + dble(i3-1)*side_UC(3,1)
              yt = dble(i1-1)*side_UC(1,2) + dble(i2-1)*side_UC(2,2) + dble(i3-1)*side_UC(3,2)
              zt = dble(i1-1)*side_UC(1,3) + dble(i2-1)*side_UC(2,3) + dble(i3-1)*side_UC(3,3)
              pos_eq_EC(i,1,n) = pos_eq_UC(i,1) + xt
              pos_eq_EC(i,2,n) = pos_eq_UC(i,2) + yt
              pos_eq_EC(i,3,n) = pos_eq_UC(i,3) + zt
              mass_EC(i,n) = mass_UC(i)
              at_EC(i,n) = at_UC(i)
           end do
        end do
     end do
  end do

  ! initialize pos_EC to pos_eq_EC and vel_EC to zero (in case no initconf nor restart are provided)
  pos_EC(:,:,:) = pos_eq_EC(:,:,:)
  vel_EC(:,:,:) = 0.d0
    
  ! printout
  write(*,'(2x,2a)') i2a(atoms_EC), ' atoms in the supercell'
  
  write(*,'(a)') 'Setting up the supercell...DONE.'

  return
end subroutine setup_supercell


