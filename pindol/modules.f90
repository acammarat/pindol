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
! Modules containing variables
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
module pars
  ! strings
  character(6), parameter :: version = 'v1.0.1', output_format = 'g22.14', screen_format = 'f12.4'
  character(7), parameter :: error_string = 'ERROR: '
  character(9), parameter :: warning_string = 'WARNING: '
  ! a small number
  real(8), parameter :: tiny = 1.d-13
  ! fundamental physical constants
  real(8), parameter :: kcal_to_J = 4184.d0               ! taken from LAMMPS 1 cal = 4.184 J
  real(8), parameter :: eV_to_J = 1.602176634d-19         ! taken from NIST, 1 e = 1.602176634 x 10-19 C 
  real(8), parameter :: avogadro_constant = 6.02214076d23 ! taken from NIST, 1 mol = 6.02214076 x 10^23 part.
  real(8), parameter :: amu_to_kg = 1.66053906660d-27     ! taken from NIST, 1 amu = 1.66053906660 x 10-27 kg
  real(8), parameter :: kB = 1.380649d-23                 ! taken from NIST, kB = 1.380649 x 10-23 J K^-1
  !real(8), parameter :: kB = 1.3806504d-23 ! taken from LAMMPS, for comparison only
  ! derived physical constants
  real(8), parameter :: eV_to_kcal = eV_to_J * avogadro_constant / kcal_to_J
  real(8), parameter :: eV_to_ieu = 1.d-10 * eV_to_J / amu_to_kg
  real(8), parameter :: kB_ieu = kB / eV_to_J * eV_to_ieu
end module pars

module flags
  logical, save :: nd_flag, nve_flag, nvt_flag
  logical, save :: writerestart_flag, readrestart_flag, initvel_flag
  logical, save :: initconf_flag, finalconf_flag, distort_flag
end module flags

module refconf
  integer, save :: atoms_UC
  real(8), save :: side_UC(3,3)
  integer, save, allocatable :: natoms_UC(:)
  real(8), save, allocatable :: mass_UC(:), pos_eq_UC(:,:)
  character(2), save, allocatable :: at_UC(:)
end module refconf

module eigen
  integer, save :: ncells(3), nq, ndof, ndof_sub, ndof_temperature
  integer, save :: npG, npH, npS, npunique, nptot
  integer, save :: npGq, npHq, npSq, npuniqueq, nptotq, npdiffq
  real(8), save, allocatable :: vec(:,:), vec_red(:,:), freq2(:), freq2_nd(:)
  real(8), save, allocatable :: eigen_real(:,:), eigen_imag(:,:)
end module eigen

module supercell
  integer, save :: ncells_tot, atoms_EC
  real(8), save :: side_EC(3,3)
  real(8), save, allocatable :: mass_EC(:,:), pos_eq_EC(:,:,:)
  real(8), save, allocatable :: pos_EC(:,:,:), vel_EC(:,:,:)
  character(2), save, allocatable :: at_EC(:,:)
end module supercell

module dvodemod
  integer, save :: dvode_method, LRW, LIW, ITOL, ITASK, IOPT, ISTATE
  real(8), save :: dvode_tol, RTOL, ATOL
  character(4), save :: dvode_tol_flag
  integer, save, allocatable :: IWORK(:)
  real(8), save, allocatable :: RWORK(:)
end module dvodemod

module input
  ! units
  real(8), save :: ieu_to_energy_units, itu_to_time_units
  character(10), save :: units
  ! atoms
  integer, save :: atom_types
  real(8), save, allocatable :: mass_pertype(:)
  character(2), save, allocatable :: at_pertype(:)
  ! refconf
  character(10), save :: refconf_filetype
  character(200), save :: refconf_filename
  ! initconf
  character(10), save :: initconf_filetype
  character(200), save :: initconf_filename
  ! finalconf
  character(10), save :: finalconf_filetype
  character(200), save :: finalconf_filename
  ! initvel
  integer, save :: seed_initvel
  real(8), save :: init_temperature
  ! distort
  integer, save :: ndistort, seed_distort
  real(8), save :: distort_amp
  character(4), save :: distort_mode
  integer, save, allocatable :: distort(:,:)
end module input

module nd
  integer, save :: runsteps, coordsteps, velsteps, accsteps
  integer, save :: nphi, jl1(24), jl2(24)
  real(8), save :: timeunit, thermotime, thermomass, temperature_begin
  real(8), save :: temperature_end, time_init, time_total, temperature
  real(8), save :: restart_time, restart_s, restart_sd, restart_thermomass
  character(200), save :: coord_filename, vel_filename, acc_filename
  character(200), save :: writerestart_filename, readrestart_filename
  logical, save :: temperature_ramp
  integer, save, allocatable :: phi_map(:,:)
  real(8), save, allocatable :: Y(:), restart_Y(:), phi(:,:)
end module nd

module int2str
contains
  function i2a(i) result(out)
    character(:), allocatable :: out
    integer, intent(in) :: i
    character(range(i)+2) :: x
    write(x,'(i0)') i
    out = trim(x)
  end function i2a
end module int2str

module masses
  integer, parameter :: n_phonopy = 92
  character(2), parameter, dimension(n_phonopy) :: at_phonopy = (/"H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U "/)
  real(8), parameter, dimension(n_phonopy) :: mass_phonopy =(/1.00794d0, 4.002602d0, 6.941d0, 9.012182d0, 10.811d0, 12.0107d0, 14.0067d0, 15.9994d0, 18.9984032d0, 20.1797d0, 22.98976928d0, 24.3050d0, 26.9815386d0, 28.0855d0, 30.973762d0, 32.065d0, 35.453d0, 39.948d0, 39.0983d0, 40.078d0, 44.955912d0, 47.867d0, 50.9415d0, 51.9961d0, 54.938045d0, 55.845d0, 58.933195d0, 58.6934d0, 63.546d0, 65.38d0, 69.723d0, 72.64d0, 74.92160d0, 78.96d0, 79.904d0, 83.798d0, 85.4678d0, 87.62d0, 88.90585d0, 91.224d0, 92.90638d0, 95.96d0, 98.d0, 101.07d0, 102.90550d0, 106.42d0, 107.8682d0, 112.411d0, 114.818d0, 118.710d0, 121.760d0, 127.60d0, 126.90447d0, 131.293d0, 132.9054519d0, 137.327d0, 138.90547d0, 140.116d0, 140.90765d0, 144.242d0, 145.d0, 150.36d0, 151.964d0, 157.25d0, 158.92535d0, 162.500d0, 164.93032d0, 167.259d0, 168.93421d0, 173.054d0, 174.9668d0, 178.49d0, 180.94788d0, 183.84d0, 186.207d0, 190.23d0, 192.217d0, 195.084d0, 196.966569d0, 200.59d0, 204.3833d0, 207.2d0, 208.98040d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 227.d0, 232.03806d0, 231.03588d0, 238.02891d0/)
end module masses

module io_units
  integer, parameter :: inp_unit          = 10
  integer, parameter :: qmatrix_unit      = 11
  integer, parameter :: freq_unit         = 12
  integer, parameter :: phi_unit          = 13
  integer, parameter :: refconf_unit      = 14
  integer, parameter :: initconf_unit     = 15
  integer, parameter :: readrestart_unit  = 16
  integer, parameter :: finalconf_unit    = 17
  integer, parameter :: writerestart_unit = 18
  integer, parameter :: coord_unit        = 20
  integer, parameter :: vel_unit          = 21
  integer, parameter :: acc_unit          = 22
  integer, parameter :: tmp_unit          = 99
end module io_units
  
