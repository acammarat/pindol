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
! This is the main of the code, it calls the different subroutine according to the input provided by the user
! it also compute and print the execution walltime
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! v 1.0.1
! - corrected a bug when treating q-points at the border of the Brillouin zone
!
! v 1.0
! - first production release
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
  use pars, only: version, screen_format
  use flags, only: nd_flag, readrestart_flag, writerestart_flag, finalconf_flag, initconf_flag, initvel_flag, distort_flag
  use int2str
  use omp_lib
  implicit none
  integer :: threads, hours, minutes, seconds, milliseconds
  real(8) :: time1, time2
  
  ! call the intrinsic CPU_TIME at the start of simulation
  call cpu_time ( time1 )   

  ! write ascii art
  write(*,'(a)')
  write(*,'(a)')  "        _           _       _  "
  write(*,'(a)')  "  _ __ (_)_ __   __| | ___ | | "
  write(*,'(a)')  " | '_ \| | '_ \ / _` |/ _ \| | "
  write(*,'(a)')  " | |_) | | | | | (_| | (_) | | "
  write(*,'(a)')  " | .__/|_|_| |_|\__,_|\___/|_| "
  write(*,'(2a)') " |_|                           ", version
  write(*,'(a)')

  ! write number of threads
  threads = omp_get_num_threads()
  write(*,'(3a)') 'Running on ', i2a(threads), ' threads'
  
  ! read 'pindol.inp' file
  call read_input

  ! read reference configuration
  call read_refconf
     
  ! read eigenvectors
  call read_eigen

  ! read phi (no need to do it, if no normal dynamics has to be performed)
  if ( nd_flag ) call read_phi
  
  ! setup supercell
  call setup_supercell

  ! read initial configuration
  if ( initconf_flag ) call read_initconf
  
  ! read restart file
  if ( readrestart_flag ) call read_restart

  ! generate a distorted configuration with random displacements
  if ( distort_flag ) call generate_distortions
  
  ! generate a Gaussian distribution of velocities scaled according to the provided temperature
  if ( initvel_flag ) call initialize_velocities
           
  ! perform normal dynamics
  if ( nd_flag ) call run_nd

  ! write the final configuration
  if ( finalconf_flag ) call write_finalconf
  
  ! write restart file
  if ( writerestart_flag ) call write_restart

  ! deallocation
  call deallocate_all

  ! write credits
  call credits
  
  ! call the intrinsic CPU_TIME at the end of simulation
  call cpu_time ( time2 )

  ! print the total simulation wall time
  hours = floor( (time2-time1) / 3600.d0 )
  minutes = floor( ( (time2-time1) - dble(hours)*3600.d0 ) / 60.d0 )
  seconds = floor(time2-time1) - hours*3600 - minutes*60
  milliseconds = nint( ( (time2-time1) - dble( hours*3600 + minutes*60 + seconds ) ) * 1000.d0 )
  write(*,'(a,'//screen_format//',9a)') 'Total simulation walltime is ', time2-time1, ' s (', &
       i2a(hours), ' hours, ', i2a(minutes), ' minutes, ', i2a(seconds), ' seconds and ', i2a(milliseconds), ' milliseconds).'
  
  stop
end program main
