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
! This subroutine writes a restart file in order to continue a previous simulation with the same conditions
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine write_restart
  ! Note on units:
  !   [mass_pertype] = amu
  !   [side_EC,pos_EC] = Ang
  !   [vel_EC] = Ang fs^-1
  !   [restart_time] = fs
  !   [Y(1:nptotq)] = amu^1/2 Ang
  !   [Y(nptotq+1:2*nptotq)] = amu^1/2 Ang fs^-1
  !   [Y(2*nptotq+1)] = .
  !   [Y(2*nptotq+2)] = fs^-1
  !   [thermomass] = amu Ang^2 K^-1
  use pars, only: error_string, version
  use flags, only: nve_flag, nvt_flag, nd_flag
  use supercell, only: pos_EC, vel_EC, atoms_EC, side_EC
  use eigen, only: nptot, nptotq
  use input, only: atom_types, at_pertype, mass_pertype
  use nd, only: writerestart_filename, restart_time, Y, thermomass
  use io_units, only: writerestart_unit
  implicit none
  integer :: i
  
  write(*,'(a)',advance='no') 'Writing the restart file...'

  ! open the restart file
  open( unit=writerestart_unit, file=writerestart_filename, action='write', status='unknown' )

  ! write comment
  write(writerestart_unit,'(a)') '# restart file written by PINDOL'//version

  ! write the number of atomic types
  write(writerestart_unit,*) atom_types, ' # atom types'

  ! write the atomic types
  do i = 1, atom_types
     write(writerestart_unit,*) at_pertype(i), mass_pertype(i), ' # symbol, mass(amu)' ! atom symbol, amu
  end do

  ! write the total number of atoms and q-points
  write(writerestart_unit,*) atoms_EC, nptot, ' # atoms_EC, npoints'

  ! write the box
  write(writerestart_unit,*) (side_EC(1,i),i=1,3), ' # a(Ang)'
  write(writerestart_unit,*) (side_EC(2,i),i=1,3), ' # b(Ang)'
  write(writerestart_unit,*) (side_EC(3,i),i=1,3), ' # c(Ang)'

  ! write the current time (in internal units)
  if ( nd_flag ) then
     write(writerestart_unit,*) restart_time, ' # time(fs)'
  else
     write(writerestart_unit,*) 0.d0, ' # time(fs)'
  end if

  ! write the flag about the kind of dynamics performed (NVE, NVT...)
  if ( nd_flag ) then
     if ( nve_flag ) write(writerestart_unit,*) 'nve dynamics performed'
     ! write thermostat variables in case of NVT
     if ( nvt_flag ) then
        write(writerestart_unit,*) 'nvt dynamics performed'
        write(writerestart_unit,*) Y(2*nptotq+1), Y(2*nptotq+2), thermomass, ' # eta(.), etadot(fs^-1), thermomass(amu*Ang^2/K)'
     end if
  else
     write(writerestart_unit,*) 'no dynamics performed'
  end if

  ! if nd was not run, generate the vector Y
  if ( .not. nd_flag ) then
     allocate ( Y(2*nptotq), stat = i )
     if ( i /= 0 ) then
        write(0,*) error_string, 'Allocation failed for temporary array Y.'
        stop
     end if
     call cartesian_to_normal ( pos_EC, Y(1:nptotq), 'pos' ) ! project the initial configuration   
     call cartesian_to_normal ( vel_EC, Y(nptotq+1:2*nptotq), 'vel' ) ! project the initial velocities
  end if
  
  ! write normal coordinates
  write(writerestart_unit,*) (Y(i),i=1,nptotq), ' # amu^1/2*Ang'

  
  ! write normal velocities
  write(writerestart_unit,*) (Y(i),i=nptotq+1,2*nptotq), ' # amu^1/2*Ang/fs'

  close( writerestart_unit )

  write(*,'(a)') 'DONE.'
  

  return
end subroutine write_restart
