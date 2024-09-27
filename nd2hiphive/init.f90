!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nd2hiphive. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Processes an ND trajectory and Creates input files for hiPhive 
! (https://hiphive.materialsmodeling.org/)
! to calculate the effective interatomic force constants
! 
! If used for production, you should cite
! A. Cammarata, M. Dasic, P. Nicolini, J. Chem. Phys. 161, 084111 (2024)
! https://doi.org/10.1063/5.0224108
!
!    This file is part of nd2hiphive.
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
  use functions, only: f2a
  use var, only: t1, t2, &
       natom_types, natoms_UC, at_pertype, &
       atoms_UC, at_UC, mass_UC, mass_pertype, pos_eq_UC, &
       side_eq_UC, infile_qdd, skip, &
       nq_chk, infile_q
  implicit none
  integer :: i, j, k
  real(8) :: x_tmp, y_tmp, z_tmp, volfact
  character(20) :: word
  character(2), allocatable :: at_tmp(:)
  character(256) :: infile
  logical :: file_exists

  call show_logo

  call get_command_argument(1,infile)
  if ( (infile == '-h') .or. (infile == '') ) then
    write(*,'(*(a))') ' Syntax: nd2hiPhive <setting file>'
    write(*,*)
    stop
  end if

  inquire(file=infile,exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(a19,a15,a11)') ' ERROR: input file ',infile,' not found.'
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

  read(10,*) infile_q
  inquire(file=infile_q,exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(*(a))') ' ERROR: ND trajectory file ',trim(infile_q),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  Reading ND trajectory from file: ', trim(infile_q)

  read(10,*) infile_qdd
  inquire(file=infile_qdd,exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(*(a))') ' ERROR: ND acceleration file ',trim(infile_qdd),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  Reading ND accelerations from file: ', trim(infile_qdd)

  read(10,*) natom_types
  allocate ( mass_pertype(natom_types), stat = i )
  if ( i /= 0 ) stop ' ERROR: allocation failed for MASS_PERTYPE'
  allocate ( at_pertype(natom_types), stat = i )
  if ( i /= 0 ) stop ' ERROR: allocation failed for AT_PERTYPE'

  do i = 1, natom_types
     read(10,*) at_pertype(i), mass_pertype(i) ! uma
  end do

  read(10,*) t1, t2, skip
  write(*,'(*(a))') '  Time window: [',trim(f2a(t1)),', ',trim(f2a(t2)),'], ', trim(f2a(t2-t1))
  write(*,'(*(a))') '  Max skip time: ', trim(f2a(skip))

  close(10)

  open(unit=15,file=infile,action='READ') 
  read(15,*)
  read(15,*) volfact
  do i = 1, 3
    read(15,*) side_eq_UC(i,:)
  end do
  side_eq_UC(:,:) = side_eq_UC(:,:)*volfact

  allocate ( at_tmp(natom_types), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for AT_TMP'
  read(15,*) (at_tmp(i),i=1,natom_types)
  do i = 1, natom_types
     if ( at_pertype(i) /= at_tmp(i) ) stop ' ERROR: atomic symbols in input and geometry file do not match.'
  end do

  allocate ( natoms_UC(natom_types), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for NATOMS_UC'
  read(15,*) (natoms_UC(i),i=1,natom_types)

  atoms_UC = 0
  do i = 1, natom_types
     atoms_UC = atoms_UC + natoms_UC(i)
  end do
  nq_chk = 3 * atoms_UC

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
!        mass_UC(k) = sqrt(1.d0/mass_UC(k))
        at_UC(k) = at_pertype(i)
     end do
  end do

  deallocate ( at_tmp, stat = i )
  if ( i /= 0 ) stop 'Deallocation failed for AT_TMP'

  write(*,*) 'Initialization done.'
  
  return
end subroutine init
