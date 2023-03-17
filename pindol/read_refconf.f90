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
! This subroutine reads the reference (equilibrium and unit cell) configuration from a file
! currently LAMMPS data (atomic style only) and POSCAR formats are supported
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_refconf
  use input, only: refconf_filetype, refconf_filename
  use refconf, only: atoms_UC
  use int2str
  use io_units, only: refconf_unit
  implicit none
  
  write(*,'(a)') 'Reading reference configuration...'

  open( unit=refconf_unit, file=refconf_filename, action='read', status='old' ) 

  if ( refconf_filetype == 'poscar' ) then

     call read_refconf_poscar
     
  else if ( refconf_filetype == 'lammpsdata' ) then

     call read_refconf_lammpsdata

  end if

  close( refconf_unit )

  ! printout
  write(*,'(2x,2a)') i2a(atoms_UC), ' atoms found in the reference cell'
  
  write(*,'(a)') 'Reading reference configuration...DONE.'
  
  return
end subroutine read_refconf

! This subroutine reads the reference (equilibrium and unit cell) configuration in LAMMPS data format (atomic style only)
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_refconf_lammpsdata
! Note on units:
!   [xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz] = Ang for real, metal
!   [side_UC] = Ang
!   [mass_tmp] = amu for real, metal
!   [mass_pertype,mass_UC] = amu
!   [pos_tmp] = Ang for real, metal
!   [pos_eq_UC] = Ang
  use pars, only: error_string, tiny
  use input, only: atom_types, mass_pertype, at_pertype
  use refconf, only: side_UC, natoms_UC, pos_eq_UC, atoms_UC, mass_UC, at_UC
  use io_units, only: refconf_unit
  implicit none
  integer :: i, j, k, atom_types_tmp, istat, type_tmp, ix(3)
  real(8) :: mass_tmp, xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, xxx, pos_tmp(3)
  character(20) :: word, word2, word3


  ! skip two lines
  read(refconf_unit,*)
  read(refconf_unit,*)

  ! read total number of atoms (it must be the 3rd line)
  read(refconf_unit,*) atoms_UC, word
  if ( word /= 'atoms' ) then
     write(0,*) error_string, 'Wrong format of the data file (the third line must contain the number of atoms).'
     stop
  end if

  ! find number of atom types
  do 
     read(refconf_unit,*) i, word
     if ( word == 'atom' ) then
        backspace(refconf_unit)
        read(refconf_unit,*) i, word, word2
        if ( word2 == 'types' ) then
           backspace(refconf_unit)
           read(refconf_unit,*) atom_types_tmp
           if ( atom_types_tmp /= atom_types ) then
              write(0,*) error_string, 'The number of atomic types in the data file does not match with what is reported in the ATOMS section.'
              stop
           end if
           exit
        end if
     end if
  end do

  ! find box and set up box
  do 
     read(refconf_unit,*,iostat=istat) xxx, xxx, word, word2
     if ( istat > 0 ) then ! header is over
        backspace(refconf_unit)
        exit
     end if
     if ( word == 'xlo' ) then
        backspace(refconf_unit)
        read(refconf_unit,*) xlo, xhi
     else if ( word == 'ylo' ) then
        backspace(refconf_unit)
        read(refconf_unit,*) ylo, yhi
     else if ( word == 'zlo' ) then
        backspace(refconf_unit)
        read(refconf_unit,*) zlo, zhi
     end if
     if ( word2 == 'xy' ) then
        backspace(refconf_unit)
        read(refconf_unit,*) xy, xz, yz
     end if
  end do

  ! calculate UC box variables
  side_UC(1,1) = xhi - xlo
  side_UC(1,2:3) = 0.d0
  side_UC(2,1) = xy
  side_UC(2,2) = yhi - ylo
  side_UC(2,3) = 0.d0
  side_UC(3,1) = xz
  side_UC(3,2) = yz
  side_UC(3,3) = zhi - zlo
  
  ! find masses and atoms sections
  do 
     read(refconf_unit,*,iostat=istat) word
     if ( istat < 0 ) exit ! end of file

     if ( word == 'Masses' ) then
        ! read & check masses
        read(refconf_unit,*) 
        do i = 1, atom_types
           read(refconf_unit,*) j, mass_tmp
           if ( abs(mass_tmp-mass_pertype(i)) > tiny ) then
              write(0,*) error_string, 'The atomic masses in the data file do not match with what is reported in the ATOMS section.'
              stop
           end if
        end do

     else if ( word == 'Atoms' ) then
        backspace(refconf_unit)
        read(refconf_unit,*) word, word2, word3
        if ( word2 == '#' .and. word3 == 'atomic' ) then ! for now it works with atomic style only
           ! allocate pertype arrays
           allocate ( mass_UC(atoms_UC), stat = i )
           if ( i /= 0 ) then
              write(0,*) error_string, 'Allocation failed for array MASS_UC.'
              stop
           end if
           allocate ( at_UC(atoms_UC), stat = i )
           if ( i /= 0 ) then
              write(0,*) error_string, 'Allocation failed for array AT_UC.'
              stop
           end if
           allocate ( natoms_UC(atom_types), stat = i )
           if ( i /= 0 ) then
              write(0,*) error_string, 'Allocation failed for array NATOMS_UC.'
              stop
           end if
           natoms_UC(:) = 0
           allocate ( pos_eq_UC(atoms_UC,3), stat = i )
           if ( i /= 0 ) then
              write(0,*) error_string, 'Allocation failed for array POS_EQ_UC.'
              stop
           end if
           ! read atomic coordinates and set masses and atom types
           read(refconf_unit,*) 
           do i = 1, atoms_UC
              read(refconf_unit,*) k, type_tmp, (pos_tmp(j),j=1,3), (ix(j),j=1,3)
              mass_UC(k) = mass_pertype(type_tmp)
              at_UC(k) = at_pertype(type_tmp)
              natoms_UC(type_tmp) = natoms_UC(type_tmp) + 1
              ! apply translation and generate the "unwrapped" positions
              pos_eq_UC(k,:) = pos_tmp(:) + side_UC(1,:)*dble(ix(1)) + side_UC(2,:)*dble(ix(2)) + side_UC(3,:)*dble(ix(3))
           end do
        end if
     end if
     
  end do

  return
end subroutine read_refconf_lammpsdata

! This subroutine reads the reference (equilibrium and unit cell) configuration in POSCAR format
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_refconf_poscar
  ! Note on units:
  !   [side_UC] = Ang
  !   [mass_pertype,mass_UC] = amu
  !   [pos_eq_UC] = Ang
  use pars, only: error_string, tiny
  use input, only: atom_types, at_pertype, mass_pertype
  use refconf, only: side_UC, natoms_UC, pos_eq_UC, atoms_UC, mass_UC, at_UC
  use io_units, only: refconf_unit
  implicit none
  integer :: i, j, k
  real(8) :: scale, x_tmp, y_tmp, z_tmp
  character(2), allocatable :: at_tmp(:)
  character(1) :: word

  ! comment
  read(refconf_unit,*) 

  ! scale factor
  read(refconf_unit,*) scale
  if ( abs(scale-1.d0) > tiny ) then
     write(0,*) error_string, 'POSCAR scale factor not implemented yet (it cannot be anything other than 1).'
     stop
  end if
  
  ! read UC box
  do i = 1, 3
     read(refconf_unit,*) (side_UC(i,j),j=1,3) ! Ang
  end do

  ! temporary array
  allocate ( at_tmp(atom_types), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for temporary array AT_TMP.'
     stop
  end if
  
  ! read atomic types
  read(refconf_unit,*) (at_tmp(i),i=1,atom_types) ! atom symbols
  
  ! check if they match those in the input
  do i = 1, atom_types
     if ( at_pertype(i) /= at_tmp(i) ) then
        write(0,*) error_string, 'The symbols of the atomic types in the poscar file do not match what is reported in the ATOMS section.'
        stop
     end if
  end do
  
  ! allocate number of atoms for each type
  allocate ( natoms_UC(atom_types), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array NATOMS_UC.'
     stop
  end if
  
  ! read number of atoms for each type
  read(refconf_unit,*) (natoms_UC(i),i=1,atom_types) 
  
  ! calculate total number of atoms
  atoms_UC = 0
  do i = 1, atom_types
     atoms_UC = atoms_UC + natoms_UC(i) 
  end do
  
  ! allocate equilibrium position array
  allocate ( pos_eq_UC(atoms_UC,3), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array POS_EQ_UC.'
     stop
  end if

  ! read keyword
  read(refconf_unit,*) word
  if ( word == 'S' .or. word == 's' ) read(refconf_unit,*) word
  
  ! read atomic positions and store them in cartesian
  do i = 1, atoms_UC
     read(refconf_unit,*) (pos_eq_UC(i,j),j=1,3) ! Ang or adim
     if ( word /= 'C' .and. word /= 'c' .and. word /= 'K' .and. word /= 'k' ) then
        x_tmp = side_UC(1,1)*pos_eq_UC(i,1) + side_UC(2,1)*pos_eq_UC(i,2) + side_UC(3,1)*pos_eq_UC(i,3) 
        y_tmp = side_UC(1,2)*pos_eq_UC(i,1) + side_UC(2,2)*pos_eq_UC(i,2) + side_UC(3,2)*pos_eq_UC(i,3) 
        z_tmp = side_UC(1,3)*pos_eq_UC(i,1) + side_UC(2,3)*pos_eq_UC(i,2) + side_UC(3,3)*pos_eq_UC(i,3) 
        pos_eq_UC(i,1) = x_tmp
        pos_eq_UC(i,2) = y_tmp
        pos_eq_UC(i,3) = z_tmp
     end if
  end do
     
  ! allocate pertype arrays
  allocate ( mass_UC(atoms_UC), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array MASS_UC.'
     stop
  end if
  allocate ( at_UC(atoms_UC), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array AT_UC.'
     stop
  end if
  
  ! set pertype arrays
  k = 0
  do i = 1, atom_types
     do j = 1, natoms_UC(i)
        k = k + 1
        mass_UC(k) = mass_pertype(i) ! amu
        at_UC(k) = at_pertype(i)
     end do
  end do

  ! free temporary array
  deallocate ( at_tmp, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for temporary array AT_TMP.'
     stop
  end if

  return
end subroutine read_refconf_poscar
