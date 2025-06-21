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
! This subroutine writes the final configuration
! currently LAMMPS data (atomic style only) and POSCAR formats are supported
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine write_finalconf
  use input, only: finalconf_filetype, finalconf_filename
  use io_units, only: finalconf_unit
  implicit none
  
  write(*,'(a)',advance='no') 'Writing the final configuration...'

  ! open the finalconf file
  open( unit=finalconf_unit, file=finalconf_filename, action='write', status='unknown' )

  if ( finalconf_filetype == 'poscar' ) then

     call write_finalconf_poscar
     
  else if ( finalconf_filetype == 'lammpsdata' ) then

     call write_finalconf_lammpsdata

  end if

  close( finalconf_unit )

  write(*,'(a)') 'DONE.'
  
  return
end subroutine write_finalconf

! This subroutine writes the final configuration in LAMMPS data format (atomic style only)
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine write_finalconf_lammpsdata
  ! Note on units:
  !   [side_EC,pos_EC] = Ang
  !   [mass_pertype] = amu
  !   [vel_EC] = Ang fs^-1
  use pars, only: version, error_string, tiny
  use supercell, only: pos_EC, vel_EC, atoms_EC, ncells_tot, side_EC
  use input, only: atom_types, mass_pertype, units
  use refconf, only: natoms_UC
  use io_units, only: finalconf_unit
  implicit none
  integer :: i, j, n, k, l, kk
  
  ! comment
  write(finalconf_unit,'(a)') 'LAMMPS data file written by PINDOL'//version
  write(finalconf_unit,*)
  
  ! number of atoms
  write(finalconf_unit,*) atoms_EC, ' atoms'

  ! number of atom types
  write(finalconf_unit,*) atom_types, ' atom types'

  ! write box
  write(finalconf_unit,*) 
  write(finalconf_unit,*) 0.d0, side_EC(1,1), ' xlo xhi' ! Ang
  write(finalconf_unit,*) 0.d0, side_EC(2,2), ' ylo yhi' ! Ang
  write(finalconf_unit,*) 0.d0, side_EC(3,3), ' zlo zhi' ! Ang
  if ( abs(side_EC(2,1)) > tiny .or. abs(side_EC(3,1)) > tiny .or. abs(side_EC(3,2)) > tiny ) &
       write(finalconf_unit,*) side_EC(2,1), side_EC(3,1), side_EC(3,2), ' xy xz yz'
  if ( abs(side_EC(1,2)) > tiny .or. abs(side_EC(1,3)) > tiny .or. abs(side_EC(2,3)) > tiny ) then
     write(0,*) error_string, 'Settings of the supercell not consistent with the LAMMPS data format.'
     stop
  end if

  ! write masses
  write(finalconf_unit,*) 
  write(finalconf_unit,'(a)') 'Masses'
  write(finalconf_unit,*) 
  do i = 1, atom_types
     write(finalconf_unit,*) i, mass_pertype(i) ! amu
  end do

  ! write atomic positions
  write(finalconf_unit,*) 
  write(finalconf_unit,'(a)') 'Atoms # atomic'
  write(finalconf_unit,*)
  k = 0
  kk = 0
  do i = 1, atom_types
     do j = 1, natoms_UC(i)
        k = k + 1
        do n = 1, ncells_tot 
           kk = kk + 1
           write(finalconf_unit,*) kk, i, (pos_EC(k,l,n),l=1,3), (0,l=1,3)
        end do
     end do
  end do
 
  ! write atomic velocities
  write(finalconf_unit,*) 
  write(finalconf_unit,'(a)') 'Velocities'
  write(finalconf_unit,*)
  k = 0
  kk = 0
  do i = 1, atom_types
     do j = 1, natoms_UC(i)
        k = k + 1
        do n = 1, ncells_tot 
           kk = kk + 1
           ! convert into metal units (if needed)
           if ( units == 'metal' ) then
              write(finalconf_unit,*) kk, (vel_EC(k,l,n)*1.d3,l=1,3)
           else 
              write(finalconf_unit,*) kk, (vel_EC(k,l,n),l=1,3)
           end if
        end do
     end do
  end do
  
  return
end subroutine write_finalconf_lammpsdata

! This subroutine writes the final configuration in POSCAR format
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine write_finalconf_poscar
  ! Note on units:
  !   [side_EC] = Ang
  !   [mass_pertype] = amu
  !   [pos_EC] = Ang
  !   [vel_EC] = Ang fs^-1
  use pars, only: version
  use supercell, only: pos_EC, vel_EC, ncells_tot, side_EC
  use input, only: atom_types, at_pertype
  use refconf, only: natoms_UC
  use io_units, only: finalconf_unit
  use int2str
  implicit none
  integer :: i, j, n, k, l
  
  ! comment
  write(finalconf_unit,'(a)') 'POSCAR file written by PINDOL'//version

  ! scale factor
  write(finalconf_unit,'(a)') '1.0'

  ! write box
  write(finalconf_unit,*) side_EC(1,:) ! Ang
  write(finalconf_unit,*) side_EC(2,:) ! Ang
  write(finalconf_unit,*) side_EC(3,:) ! Ang

  ! write atomic types
  write(finalconf_unit,'('//i2a(atom_types)//'(a,1x))') (at_pertype(i),i=1,atom_types)

  ! write number of atoms for each type
  write(finalconf_unit,*) (natoms_UC(i)*ncells_tot,i=1,atom_types)

  ! write atomic positions 
  write(finalconf_unit,'(a)') 'Cartesian'
  k = 0
  do i = 1, atom_types
     do j = 1, natoms_UC(i)
        k = k + 1
        do n = 1, ncells_tot
           write(finalconf_unit,*) (pos_EC(k,l,n),l=1,3) ! Ang
        end do
     end do
  end do
  
  ! write atomic velocities
  write(finalconf_unit,'(a)') 'Cartesian'
  k = 0
  do i = 1, atom_types
     do j = 1, natoms_UC(i)
        k = k + 1
        do n = 1, ncells_tot
           write(finalconf_unit,*) (vel_EC(k,l,n),l=1,3) ! Ang/fs
        end do
     end do
  end do
  
  return
end subroutine write_finalconf_poscar
