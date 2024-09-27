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

subroutine write_xhiphive
  use functions, only: i2a, f2a, inv
  use var, only: atoms_UC, ncells_tot, atoms_EC, internal_to_eVAng, &
                 infile_q, infile_qdd, nptotq, t1, t2, ncells, &
                 side_eq_UC, mass_pertype, at_pertype, pos_eq_UC, side_eq_EC, &
                 pos_eq_EC, skip, natom_types, natoms_UC, version, atoms_EC
  implicit none
  character(21), parameter :: side_FMT='(a9,9(f20.16,1x),a44)'
  integer :: i, j, k, iconf, cnt
  real(8) :: T, T_next, r, x_tmp, y_tmp, z_tmp, side_eq_UC_inv(3,3)
  real(8) :: Q(nptotq), Qdd(nptotq), for_EC(atoms_UC,3,ncells_tot), pos_EC(atoms_UC,3,ncells_tot)
  character(256) :: command

  if ( trim(infile_qdd(:))/='0' ) then

    open(11,file=infile_q,action='READ')
    open(12,file=infile_qdd,action='READ')
    write(*,'(a)') ' Converting ND trajectory -> hiPhive ...'

    open(21,file='ndhiPhive_prim.xyz',action='WRITE') 

    write(21,'(a)') trim(i2a(atoms_UC))
    write(21,side_FMT) 'Lattice="', side_eq_UC(1,:), side_eq_UC(2,:), side_eq_UC(3,:), &
      '" Properties=species:S:1:pos:R:3 pbc="T T T"'
    cnt = 0
    do i = 1, natom_types
      do j = 1, natoms_UC(i)
        cnt = cnt + 1
        write(21,'(a3,1x,3(f20.16,1x))') at_pertype(i), pos_eq_UC(cnt,:)
      end do
    end do
    close(21)

    open(21,file='ndhiPhive_prim_dir.vasp',action='WRITE') 
    write(21,'(*(a))') ' nd2hiphive v. ',version
    write(21,'(a)') '1.0'
    do i = 1, 3
      write(21,'(3(f24.16,1x))') side_eq_UC(i,:)
    end do
    write(21,'(*(a2,1x))') at_pertype(:)
    write(21,'(*(i4,1x))') natoms_UC(:)
    write(21,'(a)') 'Direct'
    side_eq_UC_inv = inv(side_eq_UC)
    do i = 1, atoms_UC
        x_tmp = side_eq_UC_inv(1,1)*pos_eq_UC(i,1) + side_eq_UC_inv(2,1)*pos_eq_UC(i,2) + side_eq_UC_inv(3,1)*pos_eq_UC(i,3) 
        y_tmp = side_eq_UC_inv(1,2)*pos_eq_UC(i,1) + side_eq_UC_inv(2,2)*pos_eq_UC(i,2) + side_eq_UC_inv(3,2)*pos_eq_UC(i,3) 
        z_tmp = side_eq_UC_inv(1,3)*pos_eq_UC(i,1) + side_eq_UC_inv(2,3)*pos_eq_UC(i,2) + side_eq_UC_inv(3,3)*pos_eq_UC(i,3) 
        write(21,'(3(f20.16,1x))') x_tmp, y_tmp, z_tmp
    end do
    close(21)

    open(21,file='ndhiPhive_phonopy.py',action='WRITE')
    write(21,'(a)',advance='no') 'atoms_phonopy = PhonopyAtoms(cell=['
    write(21,'(a,3(f20.16,a),a2)') '[',side_eq_UC(1,1),',',side_eq_UC(1,2),',',side_eq_UC(1,3),'],'
    write(21,'(a,3(f20.16,a),a2)') '                                   [',side_eq_UC(2,1),',',side_eq_UC(2,2),',',side_eq_UC(2,3),'],'
    write(21,'(a,3(f20.16,a),a2)') '                                   [',side_eq_UC(3,1),',',side_eq_UC(3,2),',',side_eq_UC(3,3),']],'
    if ( atoms_UC /= 1 ) then
      write(21,'(a)',advance='no') '                scaled_positions=['
      x_tmp = side_eq_UC_inv(1,1)*pos_eq_UC(1,1) + side_eq_UC_inv(2,1)*pos_eq_UC(1,2) + side_eq_UC_inv(3,1)*pos_eq_UC(1,3) 
      y_tmp = side_eq_UC_inv(1,2)*pos_eq_UC(1,1) + side_eq_UC_inv(2,2)*pos_eq_UC(1,2) + side_eq_UC_inv(3,2)*pos_eq_UC(1,3) 
      z_tmp = side_eq_UC_inv(1,3)*pos_eq_UC(1,1) + side_eq_UC_inv(2,3)*pos_eq_UC(1,2) + side_eq_UC_inv(3,3)*pos_eq_UC(1,3) 
      write(21,'(a,3(f20.16,a),a2)') '[',x_tmp,',',y_tmp,',',z_tmp,'],'
      do i = 2, atoms_UC-1
        write(21,'(a)',advance='no') '                                   '
        x_tmp = side_eq_UC_inv(1,1)*pos_eq_UC(i,1) + side_eq_UC_inv(2,1)*pos_eq_UC(i,2) + side_eq_UC_inv(3,1)*pos_eq_UC(i,3) 
        y_tmp = side_eq_UC_inv(1,2)*pos_eq_UC(i,1) + side_eq_UC_inv(2,2)*pos_eq_UC(i,2) + side_eq_UC_inv(3,2)*pos_eq_UC(i,3) 
        z_tmp = side_eq_UC_inv(1,3)*pos_eq_UC(i,1) + side_eq_UC_inv(2,3)*pos_eq_UC(i,2) + side_eq_UC_inv(3,3)*pos_eq_UC(i,3) 
        write(21,'(a,3(f20.16,a),a2)') '[',x_tmp,',',y_tmp,',',z_tmp,'],'
      end do
      write(21,'(a)',advance='no') '                                   '
      x_tmp = side_eq_UC_inv(1,1)*pos_eq_UC(atoms_UC,1) + side_eq_UC_inv(2,1)*pos_eq_UC(atoms_UC,2) + side_eq_UC_inv(3,1)*pos_eq_UC(atoms_UC,3) 
      y_tmp = side_eq_UC_inv(1,2)*pos_eq_UC(atoms_UC,1) + side_eq_UC_inv(2,2)*pos_eq_UC(atoms_UC,2) + side_eq_UC_inv(3,2)*pos_eq_UC(atoms_UC,3) 
      z_tmp = side_eq_UC_inv(1,3)*pos_eq_UC(atoms_UC,1) + side_eq_UC_inv(2,3)*pos_eq_UC(atoms_UC,2) + side_eq_UC_inv(3,3)*pos_eq_UC(atoms_UC,3) 
      write(21,'(a,3(f20.16,a),a3)') '[',x_tmp,',',y_tmp,',',z_tmp,']],'
    else
      write(21,'(a)',advance='no') '                scaled_positions=['
      x_tmp = side_eq_UC_inv(1,1)*pos_eq_UC(1,1) + side_eq_UC_inv(2,1)*pos_eq_UC(1,2) + side_eq_UC_inv(3,1)*pos_eq_UC(1,3) 
      y_tmp = side_eq_UC_inv(1,2)*pos_eq_UC(1,1) + side_eq_UC_inv(2,2)*pos_eq_UC(1,2) + side_eq_UC_inv(3,2)*pos_eq_UC(1,3) 
      z_tmp = side_eq_UC_inv(1,3)*pos_eq_UC(1,1) + side_eq_UC_inv(2,3)*pos_eq_UC(1,2) + side_eq_UC_inv(3,3)*pos_eq_UC(1,3) 
      write(21,'(a,3(f20.16,a),a3)') '[',x_tmp,',',y_tmp,',',z_tmp,']],'
    end if

    write(21,'(a)',advance='no') '                symbols='
    do i = 1, natom_types-1
      write(21,'(*(a))',advance='no')  "['",trim(at_pertype(i)),"']*",i2a(natoms_UC(i)),'+'
    end do
    write(21,'(*(a))')  "['",trim(at_pertype(i)),"']*",i2a(natoms_UC(i)),')'
    write(21,*)
    write(21,'(a)',advance='no') 'phonopy = Phonopy(atoms_phonopy, supercell_matrix=['
    write(21,'(*(a))') '[',i2a(ncells(1)),',0,0],'
    write(21,'(*(a))') '                                                   [0,',i2a(ncells(2)),',0],'
    write(21,'(*(a))') '                                                   [0,0,',i2a(ncells(3)),']],'
    write(21,'(a)') '                             primitive_matrix=None)'
    close(21)

    open(21,file='ndhiPhive_superc.xyz',action='WRITE') 

    write(21,'(a)') trim(i2a(atoms_EC))
    write(21,side_FMT) 'Lattice="', side_eq_EC(1,:), side_eq_EC(2,:), side_eq_EC(3,:), &
      '" Properties=species:S:1:pos:R:3 pbc="T T T"'
    cnt = 0
    do i = 1, natom_types
      do j = 1, natoms_UC(i)
        cnt= cnt + 1
        do k = 1, ncells_tot
          write(21,'(a3,1x,3(f20.16,1x))') at_pertype(i), pos_eq_EC(cnt,:,k)
        end do 
      end do
    end do
    close(21)

    open(21,file='ndhiPhive_dispfor.tmp',action='WRITE') 

    T_next = t1
    iconf = 0
    reading: do

      read(11,*,iostat=i) T
      if ( i /= 0 ) exit
      read(12,*,iostat=i) T
      if ( i /= 0 ) exit

      if ( T>=T_next ) then
        if ( T>t2 ) exit reading
        backspace(11)
        read(11,*,iostat=i) T, Q(:)
        if ( i /= 0 ) exit
        backspace(12)
        read(12,*,iostat=i) T, Qdd(:)
        if ( i /= 0 ) exit
        iconf = iconf + 1
        call normal_to_cartesian ( Q, pos_EC, 'dis' ) ! Q -> Ang
        call normal_to_cartesian ( Qdd, for_EC, 'for' )

        cnt = 0
        do i = 1, natom_types
          do j = 1, natoms_UC(i)
            cnt= cnt + 1
            do k = 1, ncells_tot
              for_EC(cnt,:,k) = mass_pertype(i)*for_EC(cnt,:,k)*internal_to_eVAng ! [eV/Ang]
              write(21,'(6(f20.16,1x))') pos_EC(cnt,:,k), for_EC(cnt,:,k)
            end do
          end do
        end do
        call random_number(r)
        if ( (r-1.d-1<tiny(1.d0)) ) r=1.d-1 ! we don't want to sample too close in time
        T_next = T_next + skip*r
      end if

    end do reading

    close(21)
    close(11)
    close(12)

    open(21,file='ndhiPhive_dispfor.xyz',action='WRITE')
    write(21,*) iconf, atoms_EC
    close(21)

    write(command,'(a)') "cat ndhiPhive_dispfor.tmp >> ndhiPhive_dispfor.xyz"
    call system(command, i)
    if ( i /= 0 ) then
      write(*,'(a)') ' ERROR: cat exited with non-zero status.'
      write(*,*)
      stop
    end if

    write(command,'(a)') "rm ndhiPhive_dispfor.tmp"
    call system(command, i)
    if ( i /= 0 ) then
      write(*,'(a)') ' ERROR: rm exited with non-zero status.'
      write(*,*)
      stop
    end if

    write(*,'(*(a))') ' Number of extracted configurations: ', trim(i2a(iconf))

  end if

  return

end subroutine write_xhiphive
