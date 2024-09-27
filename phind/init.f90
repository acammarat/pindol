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

subroutine init
  use var, only: ncells, natom_types, natoms_UC, at_pertype, &
       atoms_UC, at_UC, mass_UC, mass_pertype, pos_eq_UC, &
       side_eq_UC, pos_eq_EC, ncells_tot, at_EC, &
       side_eq_EC, atoms_EC, nq, &
       phimw_fl, phical_fl
  implicit none
  integer :: i, j, k
  integer :: i1, i2, i3, n
  real(8) :: xt, yt, zt, x_tmp, y_tmp, z_tmp
  character(20) :: word
  character(2), allocatable :: at_tmp(:)
  character(256) :: infile
  logical :: file_exists

  call show_logo

  call get_command_argument(1,infile)
  if ( (infile == '-h') .or. (infile == '') ) then
    write(*,'(a)') ' Syntax: phind <setting file>'
    write(*,'(a)') ' Export the number of threads before executing phind:'
    write(*,'(a)') ' export OMP_NUM_THREADS='
    write(*,*)
    stop
  end if

  inquire(file=infile,exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if

  write(*,*) 'Initializing...'

  open(unit=10,file=infile,action='READ')
  write(*,'(*(a))') '  Reading settings from file: ', trim(infile)
  read(10,*) infile
  inquire(file=infile,exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(*(a))') ' ERROR: POSCAR file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  Reading reference structure from file: ', trim(infile)

  read(10,*) (ncells(i),i=1,3)
  ncells_tot = ncells(1)*ncells(2)*ncells(3)
  
  read(10,*) natom_types
  allocate ( mass_pertype(natom_types), stat = i )
  if ( i /= 0 ) stop ' ERROR: allocation failed for MASS_PERTYPE'
  allocate ( at_pertype(natom_types), stat = i )
  if ( i /= 0 ) stop ' ERROR: allocation failed for AT_PERTYPE'

  do i = 1, natom_types
     read(10,*) at_pertype(i), mass_pertype(i) ! uma
  end do

  read(10,*) phimw_fl     ! phi map write
  if ( (phimw_fl<0) .or. (phimw_fl>1) ) then
    write(*,'(a)') ' ERROR: phi map option not valid.'
    write(*,*)
    stop
  end if

  read(10,*) phical_fl    ! phi_calc_option
  if ( (phical_fl<0) .or. (phical_fl>1) ) then
    write(*,'(a)') ' ERROR: phi calculation option not valid.'
    write(*,*)
    stop
  end if

  close(10)

  if ( (phimw_fl == 0) .and. (phical_fl == 0) ) then
    write(*,'(a)') ' ERROR: No option has been activated in the input file.'
    write(*,*)
    stop
  end if

  open(unit=15,file=infile,action='READ') 
  read(15,*)
  read(15,*)
  do i = 1, 3
    read(15,*) side_eq_UC(i,:)
  end do

  allocate ( at_tmp(natom_types), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for AT_TMP'
  read(15,*) (at_tmp(i),i=1,natom_types)
  do i = 1, natom_types
    if ( at_pertype(i) /= at_tmp(i) ) then
      write(*,'(a)') ' ERROR: atomic symbols in input and geometry file do not match.'
      write(*,*)
      stop
    end if
  end do

  allocate ( natoms_UC(natom_types), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for NATOMS_UC'
  read(15,*) (natoms_UC(i),i=1,natom_types)

  atoms_UC = 0
  do i = 1, natom_types
     atoms_UC = atoms_UC + natoms_UC(i)
  end do
  nq = 3 * atoms_UC

  allocate ( pos_eq_UC(atoms_UC,3), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for POS_EQ_UC'

  read(15,*) word
  do i = 1, atoms_UC
     read(15,*) (pos_eq_UC(i,j),j=1,3) ! Ang or adim
     if ( word == 'Direct' ) then
        x_tmp = side_eq_UC(1,1)*pos_eq_UC(i,1) + side_eq_UC(2,1)*pos_eq_UC(i,2) + side_eq_UC(3,1)*pos_eq_UC(i,3)
        y_tmp = side_eq_UC(1,2)*pos_eq_UC(i,1) + side_eq_UC(2,2)*pos_eq_UC(i,2) + side_eq_UC(3,2)*pos_eq_UC(i,3)
        z_tmp = side_eq_UC(1,3)*pos_eq_UC(i,1) + side_eq_UC(2,3)*pos_eq_UC(i,2) + side_eq_UC(3,3)*pos_eq_UC(i,3)
        pos_eq_UC(i,1) = x_tmp
        pos_eq_UC(i,2) = y_tmp
        pos_eq_UC(i,3) = z_tmp
     end if
  end do
  close(15)

  allocate ( mass_UC(atoms_UC), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for MASS_UC'
  allocate ( at_UC(atoms_UC), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for AT_UC'

  k = 0
  do i = 1, natom_types
     do j = 1, natoms_UC(i)
        k = k + 1
        mass_UC(k) = mass_pertype(i) ! Kg 10^-27
        mass_UC(k) = sqrt(1.d0/mass_UC(k))
        at_UC(k) = at_pertype(i)
     end do
  end do

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

 
  deallocate ( at_tmp, stat = i )
  if ( i /= 0 ) stop 'Deallocation failed for AT_TMP'

  write(*,*) 'Initialization done.'
  
  return
end subroutine init
