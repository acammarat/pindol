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
! This subroutine reads a restart file in order to continue a previous simulation with the same conditions
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_restart
  ! Note on units:
  !   [restart_Y(1:nptotq)] = amu^1/2 Ang
  !   [restart_Y(nptotq+1:2*nptotq)] = amu^1/2 Ang fs^-1
  !   [restart_time] = fs
  !   [restart_s] = .
  !   [restart_sd] = fs^-1
  !   [restart_thermomass] = amu Ang^2 K^-1
  use pars, only: error_string, tiny
  use eigen, only: nptot, nptotq
  use input, only: atom_types, at_pertype, mass_pertype
  use supercell, only: side_EC, atoms_EC
  use nd, only:  readrestart_filename, restart_s, restart_sd, restart_thermomass, restart_time, restart_Y
  use io_units, only: readrestart_unit
  implicit none
  integer :: i, j, atom_types_tmp, atoms_EC_tmp, nptot_tmp
  real(8) :: mass_pertype_tmp, side_EC_tmp(3,3)
  character(2) :: at_pertype_tmp
  character(3) :: ensemble
  
  write(*,'(a)',advance='no') 'Reading the restart file...'

  ! allocate restart array
  allocate ( restart_Y(2*nptotq), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array RESTART_Y.'
     stop
  end if
  
  ! open the restart file
  open( unit=readrestart_unit, file=readrestart_filename, action='read', status='old' )

  ! comment
  read(readrestart_unit,*) 

  ! read the number of atomic types
  read(readrestart_unit,*) atom_types_tmp
  if ( atom_types_tmp /= atom_types ) then
     write(0,*) error_string, 'The number of atomic types in the restart file does not match what reported in the input.'
     stop
  end if

  ! read the atomic types
  do i = 1, atom_types
     
     read(readrestart_unit,*) at_pertype_tmp, mass_pertype_tmp ! atom symbol, amu

     ! check if the atomic symbol coincides
     if ( at_pertype(i) /= at_pertype_tmp ) then
        write(0,*) error_string, 'The symbols of atomic types in the restart file do not match what reported in the input.'
        stop
     end if

     ! check if the mass coincides
     if ( abs(mass_pertype_tmp-mass_pertype(i)) > tiny ) then
        write(0,*) error_string, 'The atomic masses in the restart file do not match what reported in the input.'
        stop
     end if

  end do

  ! read the total number of atoms and q-points
  read(readrestart_unit,*) atoms_EC_tmp, nptot_tmp

  ! check if the total number of atoms coincides
  if ( atoms_EC_tmp /= atoms_EC ) then
     write(0,*) error_string, 'The total number of atoms in the supercell does not match the current settings.'
     stop
  end if

  ! check if the number of q-points coincides
  if ( nptot_tmp /= nptot ) then
     write(0,*) error_string, 'The number of q-points does not match the current settings.'
     stop
  end if

  ! read the box
  do i = 1, 3
     read(readrestart_unit,*) (side_EC_tmp(i,j),j=1,3)
  end do

  ! check if the box matches
  do i = 1, 3
     do j = 1, 3
        if ( abs(side_EC_tmp(i,j)-side_EC(i,j)) > tiny ) then
           write(0,*) error_string, 'Sizes of the supercell not consistent with the reference cell provided.'
           stop
        end if
     end do
  end do
  
  ! read the current time (in internal units)
  read(readrestart_unit,*) restart_time

  ! read the flag about the kind of dynamics performed (NVE, NVT)
  read(readrestart_unit,*) ensemble

  ! read thermostat variables in case of NVT
  if ( ensemble == 'nvt' ) read(readrestart_unit,*) restart_s, restart_sd, restart_thermomass

  ! read normal coordinates
  read(readrestart_unit,*) (restart_Y(i),i=1,nptotq)
  
  ! read normal velocities
  read(readrestart_unit,*) (restart_Y(i),i=nptotq+1,2*nptotq)
  
  close( readrestart_unit )

  write(*,'(a)') 'DONE.'

  return
end subroutine read_restart
