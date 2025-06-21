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
! This subroutine reads the 'pindol.inp' file
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input
  use pars, only: error_string, warning_string
  use flags, only: nd_flag, finalconf_flag, initconf_flag, readrestart_flag, writerestart_flag, initvel_flag, distort_flag
  use io_units, only: inp_unit
  implicit none
  logical :: atoms_check, refconf_check, units_check, finalconf_check
  logical :: initconf_check, nd_check, readrestart_check
  logical :: writerestart_check, initvel_check, distort_check
  character(100) :: keyword

  write(*,'(a)') 'Reading "pindol.inp"...'
  
  open( unit=inp_unit, file='pindol.inp', action='read', status='old' )

  ! init checks
  atoms_check = .false.
  refconf_check = .false.
  initconf_check = .false.
  nd_check = .false.
  units_check = .false.
  finalconf_check = .false.
  readrestart_check = .false.
  writerestart_check = .false.
  initvel_check = .false.
  distort_check = .false.
  
  ! init flags
  initconf_flag = .false.
  readrestart_flag = .false.
  nd_flag = .false.
  finalconf_flag = .false.
  writerestart_flag = .false.
  initvel_flag = .false.
  distort_flag = .false.
  
  ! loop over file
  do

     ! read keyword
     read(inp_unit,*) keyword
     
     ! ignore comments
     if ( keyword(1:1) == '#' ) cycle
     
     select case ( keyword )
        
        ! file end
     case ( 'END' )
        exit
        
        ! read atoms section - MANDATORY SECTION
     case ( 'ATOMS' )
        call read_input_atoms ( atoms_check )
        
        ! read specifications of the reference configuration (primitive cell) - MANDATORY SECTION
     case ( 'REFCONF' )
        call read_input_refconf ( refconf_check )
        
        ! read units - MANDATORY SECTION
     case ( 'UNITS' )
        call read_input_units ( units_check )

        ! read specifications of the initial configuration (supercell)
     case ( 'INITCONF' )
        call read_input_initconf ( initconf_check )
        initconf_flag = .true.
        
        ! read specifications about the inizialization of atomic velocities
     case ( 'INITVEL' )
        call read_input_initvel ( initvel_check )
        initvel_flag = .true.

     ! read specifications about distortions 
     case ( 'DISTORT' )
        call read_input_distort ( distort_check )
        distort_flag = .true.
        
        ! read specification of the restart file to be read
     case ( 'READRESTART' )
        call read_input_readrestart ( readrestart_check )
        readrestart_flag = .true.

        ! read specifications about normal dynamics to be performed
     case ( 'ND' )
        call read_input_nd ( nd_check )
        nd_flag = .true.

        ! read specifications about the final configuration to be printed
     case ( 'FINALCONF' )
        call read_input_finalconf ( finalconf_check )
        finalconf_flag = .true.

        ! read specification of the restart file to be written
     case ( 'WRITERESTART' )
        call read_input_writerestart ( writerestart_check )
        writerestart_flag = .true.

     ! wrong keyword
     case default
        write(0,*) error_string, 'Unknown keyword in "pindol.inp".'
        stop

     end select

  end do
        
  close( inp_unit )

  ! check if atom types  have been provided
  if ( .not. atoms_check ) then
     write(0,*) error_string, 'ATOMS section not provided.'
     stop
  end if

  ! check if specifications  of the reference configuration have been provided
  if ( .not. refconf_check ) then
     write(0,*) error_string, 'REFCONF not provided.'
     stop
  end if
  
  ! check if the units keyword has been specified
  if ( .not. units_check ) then
     write(0,*) error_string, 'UNITS are not set.'
     stop
  end if

  ! check if none of readrestart and initconf keywords have been provided
  if ( .not. readrestart_flag .and. .not. initconf_flag ) then
     write(0,*) warning_string, 'No starting state has been provided!'
     write(0,*) warning_string, 'The simulation will start with the internally generated supercell.'

  ! check if restart is provided...
  else if ( readrestart_flag .and. .not. initconf_flag ) then
     ! ...together with initvel
     if ( initvel_flag ) then
        write(0,*) warning_string, 'Velocities from restart will be discarded and initialized according to the INITVEL keyword.'
     end if
     ! ...together with distort
     if ( distort_flag ) then
        write(0,*) warning_string, 'Random displacements will be generated according to the DISTORT keyword and added to the restart configuration.'
     end if

  ! check if initconf is provided...
  else if ( .not. readrestart_flag .and. initconf_flag ) then
     ! ...together with initvel
     if ( initvel_flag ) then
        write(0,*) warning_string, 'Velocities from initconf will be discarded and initialized according to the INITVEL keyword.'
     end if
     ! ...together with distort
     if ( distort_flag ) then
        write(0,*) warning_string, 'Random displacements will be generated according to the DISTORT keyword and added to the initconf configuration.'
     end if

  ! check if both are provided...
  else if ( readrestart_flag .and. initconf_flag ) then
     write(0,*) warning_string, 'Both INITCONF and READRESTART keywords have been provided!'
     write(0,*) warning_string, 'The keyword INITCONF will be ignored  (the initial configuration will be taken from restart).'
     initconf_flag = .false.
     ! ...together with initvel
     if ( initvel_flag ) then
        write(0,*) warning_string, 'Velocities from restart will be discarded and initialized according to the INITVEL keyword.'
     end if
     ! ...together with distort
     if ( distort_flag ) then
        write(0,*) warning_string, 'Random displacements will be generated according to the DISTORT keyword and added to the restart configuration.'
     end if

  end if
     
  write(*,'(a)') 'Reading "pindol.inp"...DONE.'
  
  return
end subroutine read_input

! This subroutine reads the ATOMS section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_atoms ( check )
  ! Note on units:
  !   [mass_pertype] = amu
  use pars, only: error_string, tiny
  use input, only: atom_types, mass_pertype, at_pertype
  use int2str
  use masses
  use io_units, only: inp_unit
  implicit none
  logical, intent(inout) :: check
  character(100) :: keyword
  integer :: i, id

  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'There can be only one ATOMS section.'
     stop
  end if
  
  ! init
  id = 0
  atom_types = -1
  
  ! loop over section
  do

     ! read keyword
     read(inp_unit,*) keyword

     ! ignore comments
     if ( keyword(1:1) == '#' ) cycle
     
     select case ( keyword )

        ! section end
     case ( 'END_ATOMS' )
        exit

        ! read number of atomic types and allocate
     case ( 'types' )
        backspace (inp_unit)
        read (inp_unit,*) keyword, atom_types
        ! check if the number is meaningful
        if ( atom_types < 0 ) then
           write(0,*) error_string, 'The number of atomic types must be greater than or equal to 1.'
           stop
        end if
        ! allocate array for atom symbols and masses
        allocate ( mass_pertype(atom_types), stat = i )
        if ( i /= 0 ) then
           write(0,*) error_string, 'Allocation failed for array MASS_PERTYPE.'
           stop
        end if
        allocate ( at_pertype(atom_types), stat = i )
        if ( i /= 0 ) then
           write(0,*) error_string, 'Allocation failed for array AT_PERTYPE.'
           stop
        end if
        
        ! read specific types
     case ( 'atom' )
        id = id + 1
        backspace (inp_unit)
        read (inp_unit,*) keyword, at_pertype(id), mass_pertype(id) ! atom symbol, amu
        ! check if the user-provided mass is the same as hard-coded in phonopy
        do i = 1, n_phonopy
           if ( at_pertype(id) == at_phonopy(i) ) then
              if ( abs(mass_pertype(id)-mass_phonopy(i)) > tiny ) then
                 write(0,*) error_string, 'Non-standard value for the mass.'
                 write(0,*) error_string, 'User-provided mass for species ', id, ' : ', mass_pertype(id), '.'
                 write(0,*) error_string, 'Mass from phonopy: ', mass_phonopy(i), '.'
                 stop
              end if
              exit
           end if
        end do
        
        ! wrong keyword
     case default
        write(0,*) error_string, 'Unknown keyword in ATOMS section.'
        stop
        
     end select
     
  end do

  ! check if all atom types have been provided
  if ( id == 0 .or. id /= atom_types ) then
     write(0,*) error_string, 'Not all atom types have been provided.'
     stop
  end if

  ! check passed
  check = .true.

  write(*,'(2x,2a,'//i2a(atom_types)//'(a,1x))') i2a(atom_types), ' atomic type(s) selected: ', (at_pertype(i),i=1,atom_types)

  return
end subroutine read_input_atoms

! This subroutine reads the REFCONF section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_refconf ( check )
  use pars, only: error_string
  use input, only: refconf_filetype, refconf_filename
  use io_units, only: inp_unit
  implicit none
  logical, intent(inout) :: check
  character(100) :: keyword

  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'The keyword REFCONF cannot be specified more than once.'
     stop
  end if
  
  ! single-line keyword
  backspace(inp_unit)
  
  ! read specifications of the reference configuration
  read(inp_unit,*) keyword, refconf_filetype, refconf_filename

  ! check if filetype is supported
  if ( refconf_filetype /= 'poscar' .and. refconf_filetype /= 'lammpsdata' ) then
     write(0,*) error_string, 'Wrong filetype specified. Supported formats are: "poscar", "lammpsdata".'
     stop
  end if
  
  ! check passed
  check = .true.
  
  write(*,'(5a)') "  '", trim(adjustl(refconf_filename)),"' ", trim(adjustl(refconf_filetype)), ' reference configuration file will be read'

  return
end subroutine read_input_refconf

! This subroutine reads the UNITS section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_units ( check )
  use pars, only: error_string, eV_to_J, eV_to_kcal, amu_to_kg
  use input, only: units, ieu_to_energy_units, itu_to_time_units
  use io_units, only: inp_unit
  use int2str
  implicit none
  logical, intent(inout) :: check
  character(100) :: keyword

  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'The keyword UNITS cannot be specified more than once.'
     stop
  end if
  
  ! single-line keyword
  backspace(inp_unit)

  ! read number of q-points
  read(inp_unit,*) keyword, units

  ! check if the units are one of the styles implemented
  if ( units /= 'real' .and. units /= 'metal' ) then
     write(0,*) error_string, 'Invalid UNITS keyword. Currently implemented: "real" and "metal".'
     stop
  end if

  ! convert constants into the style-dependent units
  if ( units == 'metal' ) then
     itu_to_time_units = 1.d-3
     ieu_to_energy_units = 1.d10 * amu_to_kg / eV_to_J
  else if ( units == 'real' ) then
     itu_to_time_units = 1.d0
     ieu_to_energy_units = 1.d10 * amu_to_kg / eV_to_J * eV_to_kcal
  end if
  
  ! check passed
  check = .true.

  write(*,'(2x,2a)') trim(adjustl(units)), ' units selected'
  
  return
end subroutine read_input_units

! This subroutine reads the INITCONF section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_initconf ( check )
  use pars, only: error_string
  use input, only: initconf_filetype, initconf_filename
  use io_units, only: inp_unit
  implicit none
  logical, intent(inout) :: check
  character(100) :: keyword

  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'The keyword INITCONF cannot be specified more than once.'
     stop
  end if
  
  ! single-line keyword
  backspace(inp_unit)
  
  ! read specifications of the initial configuration
  read(inp_unit,*) keyword, initconf_filetype, initconf_filename

  ! check if filetype is supported
  if ( initconf_filetype /= 'poscar' .and. initconf_filetype /= 'lammpsdata' ) then
     write(0,*) error_string, 'Wrong filetype specified. Supported formats are: "poscar", "lammpsdata".'
     stop
  end if
  
  ! check passed
  check = .true.
  
  write(*,'(5a)') "  '", trim(adjustl(initconf_filename)),"' ", trim(adjustl(initconf_filetype)), ' initial configuration file will be read'

  return
end subroutine read_input_initconf

! This subroutine reads the INITVEL section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_initvel ( check )
  use pars, only: error_string, screen_format
  use input, only: init_temperature, seed_initvel
  use io_units, only: inp_unit
  implicit none
  logical, intent(inout) :: check
  character(100) :: keyword

  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'The keyword INITVEL cannot be specified more than once.'
     stop
  end if
  
  ! single-line keyword
  backspace(inp_unit)
  
  ! read specifications of the initialization of velocities
  read (inp_unit,*) keyword, init_temperature, seed_initvel

  ! check on meaningful values of the parameters
  if ( init_temperature < 0.d0 ) then
     write(0,*) error_string, 'Initial temperature must be non-negative.'
     stop
  end if
  if ( seed_initvel <= 0 ) then
     write(0,*) error_string, 'Seed for random number generator must be positive.'
     stop
  end if
  
  ! check passed
  check = .true.
  
  write(*,'(a,'//screen_format//',a)') '  atomic velocities will be initialized at ', init_temperature, ' K'

  return
end subroutine read_input_initvel

! This subroutine reads the DISTORT section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_distort ( check )
  ! Note on units:
  !   [distort_amp] = Ang
  use pars, only: error_string, screen_format
  use input, only: distort_mode, ndistort, distort, distort_amp, seed_distort
  use io_units, only: inp_unit
  use int2str
  implicit none
  logical, intent(inout) :: check
  integer :: i, j
  character(100) :: keyword

  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'The keyword DISTORT cannot be specified more than once.'
     stop
  end if
  
  ! single-line keyword
  backspace(inp_unit)

  ! read the amplitude of distortions and the distortion mode
  read(inp_unit,*) keyword, seed_distort, distort_amp, distort_mode

  ! check if a meaningful value for the seed has been provided
  if ( seed_distort <= 0 ) then
     write(0,*) error_string, 'Seed for random number generator must be positive.'
     stop
  end if
  
  ! distort all atoms
  if ( distort_mode == 'all' ) then
     write(*,'(2x,a,f12.4,a)') 'global (all atoms) distortions up to a maximum of ', distort_amp, ' Ang will be applied'

  ! distort only along certain modes
  else
     distort_mode = 'mode'
     backspace(inp_unit)
     read(inp_unit,*) keyword, i, distort_amp, ndistort

     ! allocate distortion array
     allocate ( distort(ndistort,2), stat = i )
     if ( i /= 0 ) then
        write(0,*) error_string, 'Allocation failed for array DISTORT.'
        stop
     end if

     ! read modes
     backspace(inp_unit)
     read(inp_unit,*) keyword, i, distort_amp, i, (distort(j,1),distort(j,2),j=1,ndistort)
     write(*,'(2x,a,'//screen_format//','//i2a(5*ndistort+2)//'a)') 'distortions up to a maximum of ', distort_amp, ' Ang projected on ', &
          ('(',i2a(distort(i,1)),',',i2a(distort(i,2)),') ',i=1,ndistort), 'mode(s) will be applied'
  end if

  ! check passed
  check = .true.

  
  return
end subroutine read_input_distort

! This subroutine reads the READRESTART section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_readrestart ( check )
  use pars, only: error_string
  use nd, only: readrestart_filename
  use io_units, only: inp_unit
  implicit none
  logical, intent(inout) :: check
  character(100) :: keyword

  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'The keyword READRESTART cannot be specified more than once.'
     stop
  end if
  
  ! single-line keyword
  backspace(inp_unit)
  
  ! read filename of restart
  read(inp_unit,*) keyword, readrestart_filename

  ! check passed
  check = .true.
  
  write(*,'(4a)') "  '", trim(adjustl(readrestart_filename)),"' ", ' restart file will be read'

  return
end subroutine read_input_readrestart

! This subroutine reads the ND section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_nd ( check )
  ! Note on units:
  !   [thermotime] = fs
  !   [timeunit] = fs
  use pars, only: error_string, screen_format, tiny
  use flags, only: nve_flag, nvt_flag
  use dvodemod, only: dvode_method, dvode_tol_flag, dvode_tol
  use nd, only: timeunit, thermotime, temperature_begin, temperature_end, runsteps, coordsteps, &
       velsteps, accsteps, temperature_ramp, coord_filename, vel_filename, acc_filename
  use input, only: units
  use io_units, only: inp_unit
  implicit none
  logical, intent(inout) :: check
  character(100) :: keyword
  
  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'There can be only one ND section.'
     stop
  end if
  
  ! init
  nve_flag = .false.
  nvt_flag = .false.
  coordsteps = 0
  velsteps = 0
  accsteps = 0
  temperature_ramp = .false.
  
  ! loop over section
  do

     ! read keyword
     read(inp_unit,*) keyword
     
     ! ignore comments
     if ( keyword(1:1) == '#' ) cycle
     
     select case ( keyword )
        
        ! section end
     case ( 'END_ND' )
        exit

        ! type of dynamics: NVE, NVT
     case ( 'nve' )
        nve_flag = .true.

     case ( 'nvt' )
        backspace (inp_unit)
        read (inp_unit,*) keyword, temperature_begin, temperature_end, thermotime
        if ( thermotime < 0.d0 .or. abs(thermotime) < tiny ) then
           write(0,*) error_string, 'Characteristic timescale of the thermostat must be positive.'
           stop
        end if
        if ( temperature_begin < 0.d0 .or. abs(temperature_begin) < tiny .or. temperature_end < 0.d0 .or. abs(temperature_end) < tiny ) then
           write(0,*) error_string, 'Initial and final temperatures must be positive.'
           stop
        end if
        nvt_flag = .true.
        
        ! integrator
     case ( 'integrator' )
        backspace (inp_unit)
        read (inp_unit,*) keyword, timeunit, dvode_method, dvode_tol_flag, dvode_tol
        if ( dvode_method /= 10 .and. dvode_method /= 21 .and. dvode_method /= 22 ) then
           write(0,*) error_string, 'DVODE method not implemented. Possibile choices are:'
           write(0,*) error_string, ' - 10 for nonstiff (Adams) method, no Jacobian used;'
           write(0,*) error_string, ' - 21 for stiff (BDF) method, user-supplied full Jacobian;'
           write(0,*) error_string, ' - 22 for stiff method, internally generated full Jacobian.'
           stop
        end if
        if ( dvode_tol_flag /= 'abs' .and. dvode_tol_flag /= 'rel' .and. dvode_tol_flag /= 'both' ) then
           write(0,*) error_string, 'DVODE tolerance flag invalid. Possible choices are: "abs", "rel" or "both".'
           stop
        end if
        if ( dvode_tol < 0.d0 ) then
           write(0,*) error_string, 'DVODE tolerance must be positive.'
           stop
        end if
        if ( timeunit < 0.d0 ) then
           write(0,*) error_string, 'TIMEUNIT must be positive.'
           stop
        end if
        ! convert into internal units (if needed)
        if ( units == 'metal' ) then
           timeunit = timeunit * 1.d3
        end if
     
        ! runtime in steps
     case ( 'runsteps' )
        backspace (inp_unit)
        read(inp_unit,*) keyword, runsteps
        if ( runsteps < 0 ) then
           write(0,*) error_string, 'RUNSTEPS must be non-negative.'
           stop
        end if

        ! normal coordinates output
     case ( 'printcoord' )
        backspace (inp_unit)
        read(inp_unit,*) keyword, coordsteps, coord_filename
        if ( coordsteps < 1 ) then
           write(0,*) error_string, 'PRINTCOORD must be greater than or equal to 1.'
           stop
        end if
        
        ! normal velocities output
     case ( 'printvel' )
        backspace (inp_unit)
        read(inp_unit,*) keyword, velsteps, vel_filename
        if ( velsteps < 1 ) then
           write(0,*) error_string, 'PRINTVEL must be greater than or equal to 1.'
           stop
        end if
        
        ! normal accelerations output
     case ( 'printacc' )
        backspace (inp_unit)
        read(inp_unit,*) keyword, accsteps, acc_filename
        if ( accsteps < 1 ) then
           write(0,*) error_string, 'PRINTACC must be greater than or equal to 1.'
           stop
        end if
        
        ! wrong keyword
     case default
        write(0,*) error_string, 'Unknown keyword in ND section.'
        stop
        
     end select
     
  end do

  ! check if the type of dynamics is properly specified
  if ( ( nve_flag .and. nvt_flag ) .or. ( .not. nve_flag .and. .not. nvt_flag ) ) then
     write(0,*) error_string, 'Ensemble not specified, or wrong choice of ensemble. Possible mutually exclusive options are: "nve", "nvt".'
     stop
  end if

  ! check if it is a run with a temperature run
  if ( abs(temperature_begin-temperature_end) > tiny ) temperature_ramp = .true.
  
  ! check passed
  check = .true.
  
  ! printout
  if ( nve_flag ) then
     write(*,'(a)') '  NVE dynamics selected'
  else if ( nvt_flag ) then
     if ( temperature_ramp ) then
        write(*,'(3(a,'//screen_format//'),a)') '  temperature-ramp dynamics selected from ', temperature_begin, &
             ' K to ', temperature_end, ' K with a Nose-Hoover thermostat and a characteristic timescale of ', thermotime*timeunit, ' fs'
     else
        write(*,'(2(a,'//screen_format//'),a)') '  NVT dynamics selected at ', temperature_begin, &
             ' K with a Nose-Hoover thermostat and a characteristic timescale of ', thermotime*timeunit, ' fs'
     end if
  end if

  return
end subroutine read_input_nd

! This subroutine reads the FINALCONF section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_finalconf ( check )
  use pars, only: error_string
  use input, only: finalconf_filetype, finalconf_filename
  use io_units, only: inp_unit
  implicit none
  logical, intent(inout) :: check
  character(100) :: keyword

  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'The keyword FINALCONF cannot be specified more than once.'
     stop
  end if
  
  ! single-line keyword
  backspace(inp_unit)
  
  ! read specifications of the reference configuration
  read(inp_unit,*) keyword, finalconf_filetype, finalconf_filename

  ! check if filetype is supported
  if ( finalconf_filetype /= 'poscar' .and. finalconf_filetype /= 'lammpsdata' ) then
     write(0,*) error_string, 'Wrong filetype specified. Supported formats are: "poscar", "lammpsdata".'
     stop
  end if
  
  ! check passed
  check = .true.
  
  write(*,'(5a)') "  '", trim(adjustl(finalconf_filename)),"' ", trim(adjustl(finalconf_filetype)), ' final configuration file will be written'

  return
end subroutine read_input_finalconf

! This subroutine reads the WRITERESTART section of the input
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine read_input_writerestart ( check )
  use pars, only: error_string
  use nd, only: writerestart_filename
  use io_units, only: inp_unit
  implicit none
  logical, intent(inout) :: check
  character(100) :: keyword

  ! check if command has been already provided
  if ( check ) then
     write(0,*) error_string, 'The keyword WRITERESTART cannot be specified more than once.'
     stop
  end if
  
  ! single-line keyword
  backspace(inp_unit)
  
  ! read filename of restart
  read(inp_unit,*) keyword, writerestart_filename

  ! check passed
  check = .true.
  
  write(*,'(4a)') "  '", trim(adjustl(writerestart_filename)),"' ", ' restart file will be written'

  return
end subroutine read_input_writerestart
