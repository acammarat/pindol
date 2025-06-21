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
! This subroutine performs the normal dynamics
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
!   - Miljan Dašić (Czech Technical University in Prague, Czech Republic & Institute of Physics Belgrade, Serbia), miljan.dasic@scl.rs
subroutine run_nd
  ! Note on units:
  !   [Y(1:nptotq)] = amu^1/2 Ang
  !   [Y(nptotq+1:2*nptotq)] = amu^1/2 Ang fs^-1
  !   [Y(2*nptotq+1)] = .
  !   [Y(2*nptotq+2)] = fs^-1
  !   [*time*] = fs
  !   [Ekin,Epot,Etot,CoM] = amu Ang^2 fs^-2
  !   [Epot_harm,Epot_anh] = amu Ang^2 fs^-2
  !   [kB_ieu] = amu Ang^2 fs^2 K^-1
  !   [temperature,current_temperature] = K  
  !   [thermomass] = amu Ang^2 K^-1
  use dvodemod, only: LRW, LIW, ITOL, ITASK, IOPT, dvode_method, RTOL, ATOL, IWORK, RWORK, ISTATE
  use pars, only: error_string, output_format, kB_ieu
  use input, only: ieu_to_energy_units, itu_to_time_units
  use flags, only: nvt_flag, nve_flag, writerestart_flag, readrestart_flag
  use supercell, only: pos_EC, vel_EC
  use eigen, only: ndof_sub, nptotq, ndof, npGq, npHq, npSq, ndof_temperature
  use nd, only: Y, runsteps, temperature, timeunit, temperature_begin, temperature_end, thermotime, &
       thermomass, time_init, time_total, coordsteps, coord_filename, velsteps, vel_filename, accsteps, &
       acc_filename, restart_Y, restart_s, restart_sd, restart_time, temperature_ramp, restart_thermomass
  use int2str
  use io_units, only: coord_unit, vel_unit, acc_unit
  !!!! test
  use refconf, only: atoms_UC
  use supercell, only: ncells_tot
  implicit none
  integer :: i, step, IPAR
  real(8) :: time, time_out, RPAR
  real(8) :: Qdd(nptotq), Ekin, Epot_harm, Epot_anh, Epot, Etot, CoM, current_temperature
  real(8) :: Epot_thermo, Ekin_thermo
  external :: F, JAC
  !!! test
  integer :: n
  
  write(*,'(a)') 'Performing normal dynamics...'

  ! setting the total number of degrees of freedom of the system
  ndof = 2 * nptotq
  if ( nvt_flag ) ndof = ndof + 2
  ndof_temperature = npGq + npHq + 2*npSq

  ! initializing the dvode solver
  call initialize_dvode

  ! allocate vector Y with dynamical variables
  ! Y(1:npGq)                                  -> real part of normal coordinates at Gamma-point, if present
  ! Y(npGq+1:npGq+npHq)                        -> real part of normal coordinates at points in the set H, if present
  ! Y(npGq+npHq+1:npuniqueq)                   -> real part of normal coordinates at points in the set S, if present
  ! Y(npunique+1:npuniqueq+npHq)               -> imaginary part of normal coordinates at points in the set H, if present
  ! Y(npuniqueq+npHq+1:nptotq)                 -> imaginary part of normal coordinates at points in the set S, if present
  ! Y(nptotq+1:nptotq+npGq)                    -> real part of normal velocities at Gamma-point, if present
  ! Y(nptotq+npGq+1:nptotq+npGq+npHq)          -> real part of normal velocities at points in the set H, if present
  ! Y(nptotq+npGq+npHq+1:nptotq+npuniqueq)     -> real part of normal velocities at points in the set S, if present
  ! Y(nptotq+npunique+1:nptotq+npuniqueq+npHq) -> imaginary part of normal velocities at points in the set H, if present
  ! Y(nptotq+npuniqueq+npHq+1:2*nptotq)        -> imaginary part of normal velocities at points in the set S, if present
  ! if nvt_flag:
  ! Y(2*nptotq+1)                              -> thermostat variable
  ! Y(2*nptotq+2)                              -> thermostat variable time-derivative
  allocate ( Y(ndof), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array Y.'
     stop
  end if

  ! initialize the dynamical variables, time and the thermostat(s) variables
  if ( readrestart_flag ) then

     Y(1:2*nptotq) = restart_Y(1:2*nptotq)
     time_init = restart_time
     if ( nvt_flag ) then
        Y(2*nptotq+1) = restart_s 
        Y(2*nptotq+2) = restart_sd
        thermomass = restart_thermomass
     end if
     
  else
     
     ! project the initial configuration   
     call cartesian_to_normal ( pos_EC, Y(1:nptotq), 'pos' )
     
     ! project the initial velocities
     call cartesian_to_normal ( vel_EC, Y(nptotq+1:2*nptotq), 'vel' )

!!!!!!!!!!!
     open(1111, file='init_pos0.txt', action='write', status='unknown')
     open(1112, file='init_vel0.txt', action='write', status='unknown')
     do i = 1, atoms_UC
     do n = 1, ncells_tot
       write(1111,'(3(g22.14,1x))') pos_EC(i,:,n)
       write(1112,'(3(g22.14,1x))') vel_EC(i,:,n)
     end do
     end do
     close(1111)
     close(1112)

     open(1113, file='init_q.txt', action='write', status='unknown')
     write(1113,'('//i2a(nptotq+1)//output_format//')') time*itu_to_time_units, (Y(i),i=1,nptotq)
     close(1113)
     open(1113, file='init_qd.txt', action='write', status='unknown')
     write(1113,'('//i2a(nptotq+1)//output_format//')') time*itu_to_time_units, (Y(i),i=nptotq+1,2*nptotq)
     close(1113)
!!!!!!!!!!!

!     time_init = 0.d0
!     if ( nvt_flag ) then
!        Y(2*nptotq+1) = 0.d0 
!        Y(2*nptotq+2) = 0.d0
!        ! convert the characteristic timescale of the thermostat into "mass"
!        thermomass = dble(ndof_temperature-ndof_sub) * kB_ieu * timeunit*timeunit * thermotime*thermotime
!     end if

  end if
!  time = time_init
!  time_out = time ! used only if the loop over steps is skipped
  Epot_thermo = 0.d0
  Ekin_thermo = 0.d0

  ! used for temperature ramp (here for printing and in f.f90)
  temperature = temperature_begin
  time_total = dble(runsteps) * timeunit 
  
  ! calculate kinetic energy, harmonic, anharmonic and total potential energies, total energy, constant of motion, current temperature
  call calc_kinetic_energy ( Y(nptotq+1:2*nptotq), Ekin )
  call calc_harmonic_potential_energy ( Y(1:nptotq), Epot_harm )
  call calc_anharmonic_potential_energy ( Y(1:nptotq), Epot_anh )
  Epot = Epot_harm + Epot_anh
  Etot = Epot + Ekin
  CoM = Etot
  if ( nvt_flag ) then
     Ekin_thermo = 0.5d0 * thermomass * temperature * Y(2*nptotq+2)*Y(2*nptotq+2)
     Epot_thermo = dble(ndof_temperature-ndof_sub) * kB_ieu * temperature * Y(2*nptotq+1)
     CoM = CoM + Ekin_thermo + Epot_thermo
  end if
  current_temperature = 2.d0 * Ekin / ( kB_ieu * dble(ndof_temperature-ndof_sub) ) 

  ! header printout
  write(*,'(2a)') '#             time    harm. potential energy anh. potential energy      potential_energy        ', &
       'kinetic_energy          total_energy    conserved_quantity   current_temperature'
  
  ! printout zero-step
  write(*,'(8'//output_format//')') time*itu_to_time_units, Epot_harm*ieu_to_energy_units, &
       Epot_anh*ieu_to_energy_units, Epot*ieu_to_energy_units, Ekin*ieu_to_energy_units, &
       Etot*ieu_to_energy_units, CoM*ieu_to_energy_units, current_temperature

  ! print normal coordinates
  if ( coordsteps > 0 ) then
     open( unit=coord_unit, file=coord_filename, action='write', status='unknown' )
     write(coord_unit,'('//i2a(nptotq+1)//output_format//')') time*itu_to_time_units, (Y(i),i=1,nptotq)
  end if
     
  ! print normal velocities
  if ( velsteps > 0 ) then
     open( unit=vel_unit, file=vel_filename, action='write', status='unknown' )
     write(vel_unit,'('//i2a(nptotq+1)//output_format//')') time*itu_to_time_units, (Y(i),i=nptotq+1,2*nptotq)
  end if
     
  ! print normal accelerations
  if ( accsteps > 0 ) then
     call calc_normal_accelerations ( Y, Qdd )
     open( unit=acc_unit, file=acc_filename, action='write', status='unknown' )
     write(acc_unit,'('//i2a(nptotq+1)//output_format//')') time*itu_to_time_units, (Qdd(i),i=1,nptotq)
  end if

  
!     pos_EC(:,:,:) = 0.d0
!     vel_EC(:,:,:) = 0.d0
     ! project back the initial configuration   
     call normal_to_cartesian (Y(1:nptotq), pos_EC, 'pos' )
     call normal_to_cartesian (Y(nptotq+1:2*nptotq), vel_EC, 'vel' )
     open(1111, file='init_pos1.txt', action='write', status='unknown')
     open(1112, file='init_vel1.txt', action='write', status='unknown')
     do i = 1, atoms_UC
     do n = 1, ncells_tot
       write(1111,'(3(g22.14,1x))') pos_EC(i,:,n)
       write(1112,'(3(g22.14,1x))') vel_EC(i,:,n)
     end do
     end do
     close(1111)
     close(1112)

     call cartesian_to_normal ( pos_EC, Y(1:nptotq), 'pos' )
     call cartesian_to_normal ( vel_EC, Y(nptotq+1:2*nptotq), 'vel' )
     call normal_to_cartesian (Y(1:nptotq), pos_EC, 'pos' )
     call normal_to_cartesian (Y(nptotq+1:2*nptotq), vel_EC, 'vel' )
     open(1111, file='init_pos2.txt', action='write', status='unknown')
     open(1112, file='init_vel2.txt', action='write', status='unknown')
     do i = 1, atoms_UC
     do n = 1, ncells_tot
       write(1111,'(3(g22.14,1x))') pos_EC(i,:,n)
       write(1112,'(3(g22.14,1x))') vel_EC(i,:,n)
     end do
     end do
     close(1111)
     close(1112)

     call cartesian_to_normal ( pos_EC, Y(1:nptotq), 'pos' )
     call cartesian_to_normal ( vel_EC, Y(nptotq+1:2*nptotq), 'vel' )
     call normal_to_cartesian (Y(1:nptotq), pos_EC, 'pos' )
     call normal_to_cartesian (Y(nptotq+1:2*nptotq), vel_EC, 'vel' )
     open(1111, file='init_pos3.txt', action='write', status='unknown')
     open(1112, file='init_vel3.txt', action='write', status='unknown')
     do i = 1, atoms_UC
     do n = 1, ncells_tot
       write(1111,'(3(g22.14,1x))') pos_EC(i,:,n)
       write(1112,'(3(g22.14,1x))') vel_EC(i,:,n)
     end do
     end do
     close(1111)
     close(1112)



stop













!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ND RUN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! loop over steps
  do step = 1, runsteps
     
     time = time_init + dble((step-1)) * timeunit
     time_out = time_init + dble(step) * timeunit
     
     ! call the dvode subroutine
     call dvode (F, ndof, Y, time, time_out, ITOL, RTOL, ATOL, ITASK, ISTATE, &
          IOPT, RWORK, LRW, IWORK, LIW, JAC, dvode_method, RPAR, IPAR)
     
     ! check for errors in dvode
     if ( ISTATE < 0 ) then
        write(0,*) error_string, 'DVODE exited with an error (ISTATE=): ', ISTATE, '.'
        write(0,*) error_string, 'This usually means bad settings or bad input files.'
        write(0,*) error_string, 'Check the dvode.f subroutine for more information.'
        stop
     end if

     ! calculate kinetic energy, harmonic, anharmonic and total potential energies, total energy, constant of motion
     call calc_kinetic_energy ( Y(nptotq+1:2*nptotq), Ekin )
     call calc_harmonic_potential_energy ( Y(1:nptotq), Epot_harm )
     call calc_anharmonic_potential_energy ( Y(1:nptotq), Epot_anh )
     Epot = Epot_harm + Epot_anh
     Etot = Epot + Ekin
     CoM = Etot
     if ( nvt_flag ) then
        if ( temperature_ramp ) temperature = temperature_begin + ( temperature_end - temperature_begin ) / time_total * ( time_out - time_init )
        Ekin_thermo = 0.5d0 * thermomass * temperature * Y(2*nptotq+2)*Y(2*nptotq+2)
        Epot_thermo = dble(ndof_temperature-ndof_sub) * kB_ieu * temperature * Y(2*nptotq+1)
        CoM = CoM + Ekin_thermo + Epot_thermo
     end if
     current_temperature = 2.d0 * Ekin / ( kB_ieu * dble(ndof_temperature-ndof_sub) )
           
     ! printout
     write(*,'(8'//output_format//')') time_out*itu_to_time_units, Epot_harm*ieu_to_energy_units, &
          Epot_anh*ieu_to_energy_units, Epot*ieu_to_energy_units, Ekin*ieu_to_energy_units, &
          Etot*ieu_to_energy_units, CoM*ieu_to_energy_units, current_temperature
     
     ! print normal coordinates
     if ( coordsteps > 0 ) then
        if ( mod(step,coordsteps) == 0 ) write(coord_unit,'('//i2a(nptotq+1)//output_format//')') time_out*itu_to_time_units, (Y(i),i=1,nptotq)
     end if
     
     ! print normal velocities
     if ( velsteps > 0 ) then
        if ( mod(step,velsteps) == 0 ) write(vel_unit,'('//i2a(nptotq+1)//output_format//')') time_out*itu_to_time_units, (Y(i),i=nptotq+1,2*nptotq)
     end if
        
     ! print normal accelerations
     if ( accsteps > 0 ) then
        call calc_normal_accelerations ( Y, Qdd )
        if ( mod(step,accsteps) == 0 ) write(acc_unit,'('//i2a(nptotq+1)//output_format//')') time_out*itu_to_time_units, (Qdd(i),i=1,nptotq)
     end if

  end do

  ! project current values of Q and Qdot into direct coordinates
  call normal_to_cartesian ( Y(1:nptotq), pos_EC, 'pos' )
  call normal_to_cartesian ( Y(nptotq+1:2*nptotq), vel_EC, 'vel' )
  
  ! close units
  if ( coordsteps > 0 ) close( coord_unit )
  if ( velsteps > 0 ) close( vel_unit )
  if ( accsteps > 0 ) close( acc_unit )

  ! store the final time in case a restart file will be written
  if ( writerestart_flag ) restart_time = time_out
  
  write(*,'(a)') 'Performing normal dynamics...DONE.'

  return
end subroutine run_nd
