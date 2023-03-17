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
! This subroutine reads the initial configuration from a file
! currently LAMMPS data (atomic style only) and POSCAR formats are supported
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_initconf
  use input, only: initconf_filetype, initconf_filename
  use io_units, only: initconf_unit
  implicit none
  
  write(*,'(a)',advance='no') 'Reading the initial configuration...'

  ! open the initconf file
  open( unit=initconf_unit, file=initconf_filename, action='read', status='old' )

  if ( initconf_filetype == 'poscar' ) then

     call read_initconf_poscar
     
  else if ( initconf_filetype == 'lammpsdata' ) then

     call read_initconf_lammpsdata

  end if

  close( initconf_unit )

  write(*,'(a)') 'DONE.'
  
  return
end subroutine read_initconf

! This subroutine reads the initial configuration in LAMMPS data format (atomic style only) 
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_initconf_lammpsdata
  ! Note on units:
  !   [xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz] = Ang for real, metal
  !   [side_EC] = Ang
  !   [mass_tmp] = amu for real, metal
  !   [mass_pertype] = amu
  !   [pos_tmp] = Ang for real, metal
  !   [pos_EC] = Ang
  !   [vel_tmp] = Ang fs^-1 for real, Ang ps^-1 for metal
  !   [vel_EC] = Ang fs^-1
  use pars, only: error_string, tiny
  use input, only: atom_types, mass_pertype, at_pertype, units
  use supercell, only: pos_EC, vel_EC, atoms_EC, ncells_tot, at_EC, side_EC
  use io_units, only: initconf_unit
  implicit none
  integer :: i, j, k, n, kkk, atom_types_tmp, istat, type_tmp, atoms_tmp, ix(3)
  real(8) :: mass_tmp, xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, xxx, pos_tmp(3), vel_tmp(3)
  character(20) :: word, word2, word3

  ! skip two lines
  read(initconf_unit,*)
  read(initconf_unit,*)

  ! read total number of atoms (it must be the 3rd line)
  read(initconf_unit,*) atoms_tmp, word
  if ( word /= 'atoms' ) then
     write(0,*) error_string, 'Wrong format of the data file (the third line must contain the number of atoms).'
     stop
  end if

  ! check total number of atoms
  if ( atoms_tmp /= atoms_EC ) then
     write(0,*) error_string, 'Number of atoms in the supercell not consistent with the reference cell provided.'
     stop
  end if

  ! check number of atom types
  do 
     read(initconf_unit,*) i, word
     if ( word == 'atom' ) then
        backspace(initconf_unit)
        read(initconf_unit,*) i, word, word2
        if ( word2 == 'types' ) then
           backspace(initconf_unit)
           read(initconf_unit,*) atom_types_tmp
           if ( atom_types_tmp /= atom_types ) then
              write(0,*) error_string, 'Number of atomic types not consistent with what is specified in the ATOMS section.'
              stop
           end if
           exit
        end if
     end if
  end do

  ! check box
  do 
     read(initconf_unit,*,iostat=istat) xxx, xxx, word, word2
     if ( istat > 0 ) then ! header is over
        backspace(initconf_unit)
        exit
     end if
     if ( word == 'xlo' ) then
        backspace(initconf_unit)
        read(initconf_unit,*) xlo, xhi
     else if ( word == 'ylo' ) then
        backspace(initconf_unit)
        read(initconf_unit,*) ylo, yhi
     else if ( word == 'zlo' ) then
        backspace(initconf_unit)
        read(initconf_unit,*) zlo, zhi
     end if
     if ( word2 == 'xy' ) then
        backspace(initconf_unit)
        read(initconf_unit,*) xy, xz, yz
     end if
  end do
  if ( abs(xhi-xlo-side_EC(1,1)) > tiny .or. abs(yhi-ylo-side_EC(2,2)) > tiny .or. abs(zhi-zlo-side_EC(3,3)) > tiny .or. &
       abs(xy-side_EC(2,1)) > tiny .or. abs(xz-side_EC(3,1)) > tiny .or. abs(yz-side_EC(3,2)) > tiny ) then
     write(0,*) error_string, 'Sizes of the supercell not consistent with the reference cell provided.'
     stop
  end if
  if ( abs(side_EC(1,2)) > tiny .or. abs(side_EC(1,3)) > tiny .or. abs(side_EC(2,3)) > tiny ) then
     write(0,*) error_string, 'Settings of the supercell not consistent with the LAMMPS data format.'
     stop
  end if

  ! find masses, atoms and velocities sections
  do 
     read(initconf_unit,*,iostat=istat) word
     if ( istat < 0 ) exit ! end of file

     if ( word == 'Masses' ) then
        ! read & check masses
        read(initconf_unit,*) 
        do i = 1, atom_types
           read(initconf_unit,*) j, mass_tmp
           if ( abs(mass_tmp-mass_pertype(j)) > tiny ) then
              write(0,*) error_string, 'Value of the atomic mass not consistent with what is specified in the ATOMS section.'
              stop
           end if
        end do

     else if ( word == 'Atoms' ) then
        backspace(initconf_unit)
        read(initconf_unit,*) word, word2, word3
        if ( word2 == '#' .and. word3 == 'atomic' ) then ! for now it works with atomic style only
           ! read atomic coordinates
           read(initconf_unit,*)
           do j = 1, atoms_EC
              read(initconf_unit,*) kkk, type_tmp, (pos_tmp(k),k=1,3), (ix(k),k=1,3)
              i = ceiling( dble(kkk) / dble(ncells_tot) ) ! number identifying the atom in the unit cell 
              n = kkk - (i-1) * ncells_tot ! number identifying the unit cell in the supercell
              if ( at_pertype(type_tmp) /= at_EC(i,n) ) then
                 write(0,*) error_string, 'Wrong atomic type.'
                 write(0,*) error_string, 'Found atom ', j, ' (', at_pertype(type_tmp), ') in the supercell.'
                 write(0,*) error_string, 'Which should correspond to atom ', i, ' (', at_EC(i,n), ') of the primitive cell.'
                 stop
              end if
              ! apply translation and generate the "unwrapped" positions
              pos_EC(i,:,n) = pos_tmp(:) + side_EC(1,:)*dble(ix(1)) + side_EC(2,:)*dble(ix(2)) + side_EC(3,:)*dble(ix(3))
           end do
        end if

     else if ( word == 'Velocities' ) then
        ! read atomic velocities
        read(initconf_unit,*)
        do j = 1, atoms_EC
           read(initconf_unit,*) kkk, (vel_tmp(k),k=1,3)
           i = ceiling( dble(kkk) / dble(ncells_tot) ) ! number identifying the atom in the unit cell 
           n = kkk - (i-1) * ncells_tot ! number identifying the unit cell in the supercell
           vel_EC(i,:,n) = vel_tmp(:)
        end do

        ! convert into internal units (if needed)
        if ( units == 'metal' ) then
           vel_EC(:,:,:) = vel_EC(:,:,:) * 1.d-3
        end if
 
     end if
     
  end do

  return
end subroutine read_initconf_lammpsdata

! This subroutine reads the initial configuration in POSCAR format
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_initconf_poscar
  ! Note on units:
  !   [side_EC] = Ang
  !   [mass_pertype] = amu
  !   [pos_EC] = Ang
  !   [vel_EC] = Ang fs^-1
  use pars, only: error_string, tiny
  use input, only: atom_types, at_pertype
  use supercell, only: pos_EC, vel_EC, atoms_EC, ncells_tot, side_EC
  use io_units, only: initconf_unit
  implicit none
  integer :: i, j, n, kkk, atoms_tmp
  real(8) :: scale, side_tmp(3,3), x_tmp, y_tmp, z_tmp
  integer, allocatable :: natoms_tmp(:)
  character(2), allocatable :: at_tmp(:)
  character(1) :: word
  
  ! comment
  read(initconf_unit,*) 

  ! scale factor
  read(initconf_unit,*) scale
  if ( abs(scale-1.d0) > tiny ) then
     write(0,*) error_string, 'POSCAR scale factor not implemented yet (it cannot be anything other than 1).'
     stop
  end if

  ! check box
  do i = 1, 3
     read(initconf_unit,*) (side_tmp(i,j),j=1,3) ! Ang
     do j = 1, 3
        if ( abs(side_tmp(i,j)-side_EC(i,j)) > tiny ) then
           write(0,*) error_string, 'Sizes of the supercell not consistent with the reference cell provided.'
           stop
        end if
     end do
  end do

  ! allocate temporary arrays
  allocate ( at_tmp(atom_types), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for temporary array AT_TMP.'
     stop
  end if
  allocate ( natoms_tmp(atom_types), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for temporary array NATOMS_TMP.'
     stop
  end if
  
  ! read atomic types
  read(initconf_unit,*) (at_tmp(i),i=1,atom_types) ! atom symbols
  
  ! check if they match the input
  do i = 1, atom_types
     if ( at_pertype(i) /= at_tmp(i) ) then
        write(0,*) error_string, 'Wrong atomic type.'
        write(0,*) error_string, 'Found atomic type ', i, ' (', at_tmp(i), ') in the supercell;'
        write(0,*) error_string, 'Which should correspond to atom ', i, ' (', at_pertype(i), ') of the primitive cell.'
        stop
     end if
  end do
  
  ! read number of atoms for each type
  read(initconf_unit,*) (natoms_tmp(i),i=1,atom_types) 
  
  ! calculate total number of atoms
  atoms_tmp = 0
  do i = 1, atom_types
     atoms_tmp = atoms_tmp + natoms_tmp(i) 
  end do

  ! check total number of atoms
  if ( atoms_tmp /= atoms_EC ) then
     write(0,*) error_string, 'Number of atoms in the supercell not consistent with the reference cell provided.'
     stop
  end if

  ! free temporary arrays
  deallocate ( at_tmp, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for temporary array AT_TMP.'
     stop
  end if
  deallocate ( natoms_tmp, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for temporary array NATOMS_TMP.'
     stop
  end if
  
  ! read keyword
  read(initconf_unit,*) word
  if ( word == 'S' .or. word == 's' ) read(initconf_unit,*) word
  
  ! read atomic positions and store them in cartesian
  do kkk = 1, atoms_EC
     i = ceiling( dble(kkk) / dble(ncells_tot) ) ! number identifying the atom in the unit cell 
     n = kkk - (i-1) * ncells_tot ! number identifying the unit cell in the supercell
     read(initconf_unit,*) (pos_EC(i,j,n),j=1,3) ! Ang or adim
     if ( word /= 'C' .and. word /= 'c' .and. word /= 'K' .and. word /= 'k' ) then
        x_tmp = side_EC(1,1)*pos_EC(i,1,n) + side_EC(2,1)*pos_EC(i,2,n) + side_EC(3,1)*pos_EC(i,3,n) 
        y_tmp = side_EC(1,2)*pos_EC(i,1,n) + side_EC(2,2)*pos_EC(i,2,n) + side_EC(3,2)*pos_EC(i,3,n) 
        z_tmp = side_EC(1,3)*pos_EC(i,1,n) + side_EC(2,3)*pos_EC(i,2,n) + side_EC(3,3)*pos_EC(i,3,n) 
        pos_EC(i,1,n) = x_tmp
        pos_EC(i,2,n) = y_tmp
        pos_EC(i,3,n) = z_tmp
     end if
  end do

  ! read keyword
  read(initconf_unit,*,iostat=i) word

  ! exit if velocities are not provided in the file
  if ( i /= 0 ) return
  
  ! read atomic velocities and store them in cartesian
  do kkk = 1, atoms_EC
     i = ceiling( dble(kkk) / dble(ncells_tot) ) ! number identifying the atom in the unit cell 
     n = kkk - (i-1) * ncells_tot ! number identifying the unit cell in the supercell
     read(initconf_unit,*) (vel_EC(i,j,n),j=1,3) ! Ang/fs or 1/fs
     if ( word /= 'C' .and. word /= 'c' .and. word /= 'K' .and. word /= 'k' ) then
        x_tmp = side_EC(1,1)*vel_EC(i,1,n) + side_EC(2,1)*vel_EC(i,2,n) + side_EC(3,1)*vel_EC(i,3,n) 
        y_tmp = side_EC(1,2)*vel_EC(i,1,n) + side_EC(2,2)*vel_EC(i,2,n) + side_EC(3,2)*vel_EC(i,3,n) 
        z_tmp = side_EC(1,3)*vel_EC(i,1,n) + side_EC(2,3)*vel_EC(i,2,n) + side_EC(3,3)*vel_EC(i,3,n) 
        vel_EC(i,1,n) = x_tmp
        vel_EC(i,2,n) = y_tmp
        vel_EC(i,3,n) = z_tmp
     end if
  end do
  
  return
end subroutine read_initconf_poscar
