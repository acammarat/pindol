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
! This subroutine calculates the double sum:
! -1/2 ( sum_l1 sum_l2 phi_{l,l1,l2} Q_l1 Q_l2 )*
! which is used to integrate the equations of motion
! and to calculate the anharmonic part of the potential energy
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
!   - Antonio Cammarata (Czech Technical University in Prague), cammaant@fel.cvut.cz
subroutine calc_double_sum_phi ( Q, sum_phi )
  ! Note on units:
  !   [Q] = amu^1/2 Ang
  !   [phi] = amu^-1/2 Ang^-1 fs^-2
  !   [phi_sum] = amu^1/2 Ang fs^-2
  use eigen, only: nptotq, npdiffq, npuniqueq, npGq
  use nd, only: phi, phi_map, jl1, jl2
  implicit none
  real(8), intent(in) :: Q(nptotq)
  real(8), intent(out) :: sum_phi(nptotq)
  integer :: iphi, l, l1, l2
  real(8) :: phi_tmp(2), tmp1r, tmp1i, tmp2r, tmp2i

  ! 64 cases in general can be distinguished, depending if the modes involved
  ! belong to Gamma, the subset H(q) of unique points (no complex-conjugate to be considered),
  ! the subset S(q) of unique points or its complement S*(q) (see the table at the end of the file)
  !
  ! the number of cases immediately shrinks to 48 (S* are not dynamical variables),
  ! it reduces to 13 cases by considering permutation and complex-conjugation relations,
  ! and it goes further down to 9 cases by excluding the cases that a priori can be evaluated null
  !
  ! notice that there are finally 24 sub-cases, as permutation relations
  ! can be exploited within a case as well (e.g. by default, l <= l1 <= l2 )
  !
  ! --------------------------------------------------------------------------------------------------------------------------------------
  ! case    old case lambda lambda’ lambda’’ NEED TO BE CALCULATED?  Notes
  ! --------------------------------------------------------------------------------------------------------------------------------------
  ! I       i        Gamma  Gamma   Gamma    YES                     it has 4 sub-cases
  ! II               Gamma  Gamma   H        IN PRINCIPLE YES        but it should be zero
  ! III     ii       Gamma  Gamma   S        IN PRINCIPLE YES        but it should be zero
  ! IV      iii      Gamma  Gamma   S*       NO                      can be obtained via complex-conjugation of III
  ! V                Gamma  H       Gamma    NO                      can be obtained via permutation of II
  ! VI               Gamma  H       H        YES                     it has 2 sub-cases
  ! VII              Gamma  H       S        IN PRINCIPLE YES        but it should be zero
  ! VIII             Gamma  H       S*       NO                      can be obtained via complex-conjugation of VII
  ! IX      iv       Gamma  S       Gamma    NO                      can be obtained via permutation of III
  ! X                Gamma  S       H        NO                      can be obtained via permutation of VII
  ! XI      v        Gamma  S       S        IN PRINCIPLE YES        but it should be zero
  ! XII     vi       Gamma  S       S*       YES                     it has 2 sub-cases
  ! XIII    vii      Gamma  S*      Gamma    NO                      can be obtained via complex-conjugation and permutation of III
  ! XIV              Gamma  S*      H        NO                      can be obtained via complex-conjugation and permutation of VII
  ! XV      viii     Gamma  S*      S        NO                      can be obtained via complex-conjugation (and/or permutation) of XII
  ! XVI     ix       Gamma  S*      S*       NO                      can be obtained via complex-conjugation of XI
  ! XVII             H      Gamma   Gamma    NO                      can be obtained via permutation of II
  ! XVIII            H      Gamma   H        NO                      can be obtained via permutation of VI
  ! XIX              H      Gamma   S        NO                      can be obtained via permutation of VII
  ! XX               H      Gamma   S*       NO                      can be obtained via complex-conjugation and permutation of VII
  ! XXI              H      H       Gamma    NO                      can be obtained via permutation of VI
  ! XXII             H      H       H        YES                     it has 4 sub-cases
  ! XXIII            H      H       S        YES                     it has 2 sub-cases
  ! XXIV             H      H       S*       NO                      can be obtained via complex-conjugation of XXIII
  ! XXV              H      S       Gamma    NO                      can be obtained via permutation of VII
  ! XXVI             H      S       H        NO                      can be obtained via permutation of XXIII
  ! XXVII            H      S       S        YES                     it has 2 sub-cases
  ! XXVIII           H      S       S*       YES                     it has 2 sub-cases
  ! XXIX             H      S*      Gamma    NO                      can be obtained via complex-conjugation and permutation of VII
  ! XXX              H      S*      H        NO                      can be obtained via complex-conjugation and permutation of XXIII
  ! XXXI             H      S*      S        NO                      can be obtained via complex-conjugation (and/or permutation) of XXVIII
  ! XXXII            H      S*      S*       NO                      can be obtained via complex-conjugation of XXVII
  ! XXXIII  x        S      Gamma   Gamma    NO                      can be obtained via permutation of III
  ! XXXIV            S      Gamma   H        NO                      can be obtained via permutation of VII
  ! XXXV    xi       S      Gamma   S        NO                      can be obtained via permutation of XI
  ! XXXVI   xii      S      Gamma   S*       NO                      can be obtained via complex-conjugation (and/or permutation) of XII
  ! XXXVII           S      H       Gamma    NO                      can be obtained via permutation of VII
  ! XXXVIII          S      H       H        NO                      can be obtained via permutation of XXIII
  ! XXXIX            S      H       S        NO                      can be obtained via permutation of XXVII
  ! XL               S      H       S*       NO                      can be obtained via complex-conjugation (and/or permutation) of XXVIII
  ! XLI     xiii     S      S       Gamma    NO                      can be obtained via permutation of XI
  ! XLII             S      S       H        NO                      can be obtained via permutation of XXVII
  ! XLIII   xiv      S      S       S        YES                     it has 4 sub-cases
  ! XLIV    xv       S      S       S*       YES                     it has 2 sub-cases
  ! XLV     xvi      S      S*      Gamma    NO                      can be obtained via complex-conjugation (and/or permutation) of XII
  ! XLVI             S      S*      H        NO                      can be obtained via complex-conjugation (and/or permutation) of XXVIII
  ! XLVII   xvii     S      S*      S        NO                      can be obtained via permutation of XLIV
  ! XLVIII  xviii    S      S*      S*       NO                      can be obtained via complex-conjugation and permutation of XLIV
  ! XLIX    xix      S*     Gamma   Gamma    |-------------------------------------------------------------|
  ! L                S*     Gamma   H        |                                                             |
  ! LI      xx       S*     Gamma   S        |                                                             |
  ! LII     xxi      S*     Gamma   S*       |                                                             |
  ! LIII             S*     H       Gamma    |                                                             |
  ! LIV              S*     H       H        |                                                             |
  ! LV               S*     H       S        } THIS BLOCK SIMPLY DOES NOT NEED TO BE CALCULATED            {     
  ! LVI              S*     H       S*       |                                                             |
  ! LVII    xxii     S*     S       Gamma    } AS MODES BELONGING TO S* ARE NOT DYNAMICAL VARIABLES FOR US {          
  ! LVIII            S*     S       H        |                                                             |
  ! LIX     xxiii    S*     S       S        |                                                             |
  ! LX      xxiv     S*     S       S*       |                                                             |
  ! LXI     xxv      S*     S*      Gamma    |                                                             |
  ! LXII             S*     S*      H        |                                                             |
  ! LXIII   xxvi     S*     S*      S        |                                                             |
  ! LXIV    xxvii    S*     S*      S*       |-------------------------------------------------------------|
  ! --------------------------------------------------------------------------------------------------------------------------------------
  
  sum_phi(:) = 0.d0

  ! loop over non-null elements of phi
  !
  ! case I) -> l in Gamma, l1 in Gamma, l2 in Gamma
  ! --- sub-case A) -> l == l1 == l2 (this sub-case is unique, i.e. the permutation property does not apply)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q) private(iphi,l,l1,l2,phi_tmp) reduction(+:sum_phi) 
  do iphi = jl1(1), jl2(1)                                                             
     l = phi_map(iphi,1)                                                               
     l1 = phi_map(iphi,2)                                                              
     l2 = phi_map(iphi,3)                                                              
     phi_tmp(1) = phi(iphi,1)
     sum_phi(l) = sum_phi(l) - 0.5d0 * phi_tmp(1) * Q(l1) * Q(l2)                      
  end do                                                                               
  !$omp end parallel do 
  ! --- sub-case B) -> l == l1 /= l2 (in this sub-case we need to update sum_phi for l and l2)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q) private(iphi,l,l1,l2,phi_tmp) reduction(+:sum_phi) 
  do iphi = jl1(2), jl2(2)                                                             
     l = phi_map(iphi,1)                                                               
     l1 = phi_map(iphi,2)                                                              
     l2 = phi_map(iphi,3)                                                              
     phi_tmp(1) = phi(iphi,1)
     sum_phi(l)  = sum_phi(l)  -         phi_tmp(1) * Q(l1) * Q(l2)                              
     sum_phi(l2) = sum_phi(l2) - 0.5d0 * phi_tmp(1) * Q(l)  * Q(l1)                     
  end do
  !$omp end parallel do 
  ! --- sub-case C) -> l /= l1 == l2 (need to update sum_phi for l and l1; in the latter, one permutation is further possible)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q) private(iphi,l,l1,l2,phi_tmp) reduction(+:sum_phi) 
  do iphi = jl1(3), jl2(3) 
     l = phi_map(iphi,1)
     l1 = phi_map(iphi,2)
     l2 = phi_map(iphi,3)
     phi_tmp(1) = phi(iphi,1)
     sum_phi(l)  = sum_phi(l)  - 0.5d0 * phi_tmp(1) * Q(l1) * Q(l2)
     sum_phi(l1) = sum_phi(l1) -         phi_tmp(1) * Q(l)  * Q(l2)
  end do
  !$omp end parallel do 
  ! --- sub-case D) -> l /= l1 /= l2 (need to update sum_phi for all l, l1 and l2)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q) private(iphi,l,l1,l2,phi_tmp) reduction(+:sum_phi) 
  do iphi = jl1(4), jl2(4)
     l = phi_map(iphi,1)
     l1 = phi_map(iphi,2)
     l2 = phi_map(iphi,3)
     phi_tmp(1) = phi(iphi,1)
     sum_phi(l)  = sum_phi(l)  - phi_tmp(1) * Q(l1) * Q(l2)
     sum_phi(l1) = sum_phi(l1) - phi_tmp(1) * Q(l)  * Q(l2)
     sum_phi(l2) = sum_phi(l2) - phi_tmp(1) * Q(l)  * Q(l1)
  end do
  !$omp end parallel do 
  ! end case I)

  ! case VI) -> l in Gamma, l1 in H(q), l2 in H(q)
  ! this covers also:
  ! case XVIII) -> l in H(q), l1 in Gamma, l2 in H(q) via permutation
  ! case XXI) -> l in H(q), l1 in H(q), l2 in Gamma via permutation
  ! --- sub-case A) -> l1 == l2
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(5), jl2(5)
     l = phi_map(iphi,1)                                                                                                            
     l1 = phi_map(iphi,2)                                                                                                           
     l2 = phi_map(iphi,3)                                                                                                           
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - 0.5d0 * tmp2r        ! VI)
     sum_phi(l1)         = sum_phi(l1)         -         Q(l) * tmp1r ! XVIII)+XXI)   
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) -         Q(l) * tmp1i ! XVIII)+XXI)   
  end do                                                                                                                            
  !$omp end parallel do 
  ! --- sub-case B) -> l1 /= l2
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(6), jl2(6)
     l = phi_map(iphi,1)                                                                                                            
     l1 = phi_map(iphi,2)                                                                                                           
     l2 = phi_map(iphi,3)                                                                                                           
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - tmp2r        ! VI)           
     sum_phi(l1)         = sum_phi(l1)         - Q(l) * tmp1r ! XVIII)+XXI)
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - Q(l) * tmp1i ! XVIII)+XXI)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l1), Q(l1+npdiffq), tmp1r, tmp1i )
     sum_phi(l2)         = sum_phi(l2)         - Q(l) * tmp1r ! XVIII)+XXI)         
     sum_phi(l2+npdiffq) = sum_phi(l2+npdiffq) - Q(l) * tmp1i ! XVIII)+XXI) 
  end do
  !$omp end parallel do 
  ! end case VI)
  
  ! case XII) -> l in Gamma, l1 in S(q), l2 in S*(q)
  ! this covers also:
  ! case XV) -> l in Gamma, l1 in S*(q), l2 in S(q) via permutation AND (possibly) complex-conjugation
  ! case XXXVI) -> l in S(q), l1 in Gamma, l2 in S*(q) via permutation AND (possibly) complex-conjugation
  ! case XLV) -> l in S(q), l1 in S*(q), l2 in Gamma via permutation AND (possibly) complex-conjugation
  ! --- sub-case A) -> l1 == l2-npdiffq (i.e. l1 and l2 are complex-conjugates)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(7), jl2(7)
     l = phi_map(iphi,1)                                                                                                    
     l1 = phi_map(iphi,2)                                                                                                   
     l2 = phi_map(iphi,3)                                                                                                   
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2-npdiffq), -Q(l2), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - tmp2r        ! XII)+XV)
     sum_phi(l1)         = sum_phi(l1)         - Q(l) * tmp1r ! XXXVI)+XLV)
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - Q(l) * tmp1i ! XXXVI)+XLV)
  end do                                                                                                                    
  !$omp end parallel do 
  ! --- sub-case B) -> l1 /= l2-npdiffq
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(8), jl2(8)
     l = phi_map(iphi,1)                                                                                                    
     l1 = phi_map(iphi,2)                                                                                                   
     l2 = phi_map(iphi,3)                                                                                                   
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2-npdiffq), -Q(l2), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - 2.d0 * tmp2r ! XII)+XV) + XII)+XV) considering complex-conjugation
     sum_phi(l1)         = sum_phi(l1)         - Q(l) * tmp1r ! XXXVI)+XLV)
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - Q(l) * tmp1i ! XXXVI)+XLV)
     call complex_product ( phi_tmp(1), -phi_tmp(2), Q(l1), -Q(l1+npdiffq), tmp1r, tmp1i )
     sum_phi(l2-npdiffq) = sum_phi(l2-npdiffq) - Q(l) * tmp1r ! XXXVI)+XLV) considering complex-conjugation
     sum_phi(l2)         = sum_phi(l2)         - Q(l) * tmp1i ! XXXVI)+XLV) considering complex-conjugation
  end do
  !$omp end parallel do 
  ! end case XII)
  
  ! case XXII) -> l in H(q), l1 in H(q), l2 in H(q)
  ! --- sub-case A) -> l == l1 == l2 (this sub-case is unique, i.e. the permutation property does not apply)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(9), jl2(9)                                                             
     l = phi_map(iphi,1)                                                               
     l1 = phi_map(iphi,2)                                                              
     l2 = phi_map(iphi,3)                                                              
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l1), Q(l1+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l2), Q(l2+npdiffq), tmp2r, tmp2i )
     sum_phi(l)         = sum_phi(l)         - 0.5d0 * tmp2r
     sum_phi(l+npdiffq) = sum_phi(l+npdiffq) - 0.5d0 * tmp2i
  end do                                                                               
  !$omp end parallel do 
  ! --- sub-case B) -> l == l1 /= l2 (in this sub-case we need to update sum_phi for l and l2)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(10), jl2(10)                                                             
     l = phi_map(iphi,1)                                                               
     l1 = phi_map(iphi,2)                                                              
     l2 = phi_map(iphi,3)                                                              
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l1), Q(l1+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l2), Q(l2+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          -         tmp2r
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  -         tmp2i
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l2)         = sum_phi(l2)         - 0.5d0 * tmp2r
     sum_phi(l2+npdiffq) = sum_phi(l2+npdiffq) - 0.5d0 * tmp2i
  end do
  !$omp end parallel do 
  ! --- sub-case C) -> l /= l1 == l2 (need to update sum_phi for l and l1; in the latter, one permutation is further possible)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(11), jl2(11) 
     l = phi_map(iphi,1)
     l1 = phi_map(iphi,2)
     l2 = phi_map(iphi,3)
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - 0.5d0 * tmp2r
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  - 0.5d0 * tmp2i
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         -         tmp2r
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) -         tmp2i
  end do
  !$omp end parallel do 
  ! --- sub-case D) -> l /= l1 /= l2 (need to update sum_phi for all l, l1 and l2)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(12), jl2(12)
     l = phi_map(iphi,1)
     l1 = phi_map(iphi,2)
     l2 = phi_map(iphi,3)
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - tmp2r
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  - tmp2i
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         - tmp2r
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - tmp2i
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l), Q(l+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l2)         = sum_phi(l2)         - tmp2r
     sum_phi(l2+npdiffq) = sum_phi(l2+npdiffq) - tmp2i
  end do
  !$omp end parallel do 
  ! end case XXII)

  ! case XXIII) -> l in H(q), l1 in H(q), l2 in S(q)
  ! this covers also:
  ! case XXIV) -> l in H(q), l1 in H(q), l2 in S*(q) via complex-conjugation 
  ! case XXVI) -> l in H(q), l1 in S(q), l2 in H(q) via permutation 
  ! case XXX) -> l in H(q), l1 in S*(q), l2 in H(q) via permutation + complex-conjugation 
  ! case XXXVIII) -> l in S(q), l1 in H(q), l2 in H(q) via permutation 
  ! --- sub-case A) -> l == l1
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(13), jl2(13)
     l = phi_map(iphi,1)                                                                                                            
     l1 = phi_map(iphi,2)                                                                                                           
     l2 = phi_map(iphi,3)                                                                                                           
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          -         tmp2r ! XXIII)+XXVI)
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  -         tmp2i ! XXIII)+XXVI)
     call complex_product ( tmp1r, -tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          -         tmp2r ! XXIV)+XXX)
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  -         tmp2i ! XXIV)+XXX)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l), Q(l+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l2)         = sum_phi(l2)         - 0.5d0 * tmp2r ! XXIII)+XXVI)
     sum_phi(l2+npdiffq) = sum_phi(l2+npdiffq) - 0.5d0 * tmp2i ! XXIII)+XXVI)
  end do                                                                                                                            
  !$omp end parallel do 
  ! --- sub-case B) -> l /= l1
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(14), jl2(14)
     l = phi_map(iphi,1)                                                                                                            
     l1 = phi_map(iphi,2)                                                                                                           
     l2 = phi_map(iphi,3)                                                                                                           
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - tmp2r ! XXIII)+XXVI)
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  - tmp2i ! XXIII)+XXVI)
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         - tmp2r ! XXIII)+XXVI)
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - tmp2i ! XXIII)+XXVI)
     call complex_product ( tmp1r, -tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - tmp2r ! XXIV)+XXX)
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  - tmp2i ! XXIV)+XXX)
     call complex_product ( tmp1r, -tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         - tmp2r ! XXIV)+XXX)
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - tmp2i ! XXIV)+XXX)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l), Q(l+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l2)         = sum_phi(l2)         - tmp2r ! XXIII)+XXVI)
     sum_phi(l2+npdiffq) = sum_phi(l2+npdiffq) - tmp2i ! XXIII)+XXVI)
  end do
  !$omp end parallel do 
  ! end case XXIII)

  ! case XXVII) -> l in H(q), l1 in S(q), l2 in S(q)
  ! this covers also:
  ! case XXXII) -> l in H(q), l1 in S*(q), l2 in S*(q) via complex-conjugation
  ! case XXXIX) -> l in S(q), l1 in H(q), l2 in S(q) via permutation
  ! case XLII) -> l in S(q), l1 in S(q), l2 in H(q) via permutation
  ! --- sub-case A) -> l1 == l2
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(15), jl2(15)
     l = phi_map(iphi,1)                                                                                                            
     l1 = phi_map(iphi,2)                                                                                                           
     l2 = phi_map(iphi,3)                                                                                                           
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          -         tmp2r ! XXVII) + XXXII)
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         -         tmp2r ! XXXIX)+XLII)
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) -         tmp2i ! XXXIX)+XLII)
  end do
  !$omp end parallel do 
  ! --- sub-case B) -> l1 /= l2
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(16), jl2(16)
     l = phi_map(iphi,1)                                                                                                            
     l1 = phi_map(iphi,2)                                                                                                           
     l2 = phi_map(iphi,3)                                                                                                           
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - 2.d0 * tmp2r ! XXVII) + XXXII)
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         - tmp2r ! XXXIX)+XLII)
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - tmp2i ! XXXIX)+XLII)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l), Q(l+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l2)         = sum_phi(l2)         - tmp2r ! XXXIX)+XLII)
     sum_phi(l2+npdiffq) = sum_phi(l2+npdiffq) - tmp2i ! XXXIX)+XLII)
  end do
  !$omp end parallel do 
  ! end case XXVII)

  ! case XXVIII) -> l in H(q), l1 in S(q), l2 in S*(q)
  ! this covers also:
  ! case XXXI) -> l in H(q), l1 in S*(q), l2 in S(q) via permutation AND (possibly) complex-conjugation
  ! case XL) -> l in S(q), l1 in H(q), l2 in S*(q) via permutation AND (possibly) complex-conjugation
  ! case XLVI) -> l in S(q), l1 in S*(q), l2 in H(q) via permutation AND (possibly) complex-conjugation
  ! --- sub-case A) -> l1 == l2-npdiffq (i.e. l1 and l2 are complex-conjugates)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(17), jl2(17)
     l = phi_map(iphi,1)                                                                                                    
     l1 = phi_map(iphi,2)                                                                                                   
     l2 = phi_map(iphi,3)                                                                                                   
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2-npdiffq), -Q(l2), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - tmp2r ! XXVIII)+XXXI)
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  - tmp2i ! XXVIII)+XXXI)
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         - tmp2r ! XL)+XLVI)
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - tmp2i ! XL)+XLVI)
  end do
  !$omp end parallel do 
  ! --- sub-case B) -> l1 /= l2-npdiffq
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(18), jl2(18)
     l = phi_map(iphi,1)                                                                                                    
     l1 = phi_map(iphi,2)                                                                                                   
     l2 = phi_map(iphi,3)                                                                                                   
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2-npdiffq), -Q(l2), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - 2.d0 * tmp2r ! XXVIII)+XXXI) + XXVIII)+XXXI) considering complex-conjugation
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         - tmp2r ! XL)+XLVI)
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - tmp2i ! XL)+XLVI)
     call complex_product ( phi_tmp(1), -phi_tmp(2), Q(l), Q(l+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), -Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l2-npdiffq) = sum_phi(l2-npdiffq) - tmp2r ! XL)+XLVI) considering complex-conjugation
     sum_phi(l2)         = sum_phi(l2)         - tmp2i ! XL)+XLVI) considering complex-conjugation 
  end do
  !$omp end parallel do 
  ! end case XXVIII) 

  ! case XLIII) -> l in S(q), l1 in S(q), l2 in S(q)
  ! --- sub-case A) -> l == l1 == l2 (this sub-case is unique, i.e. the permutation property does not apply)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(19), jl2(19)                                                             
     l = phi_map(iphi,1)                                                               
     l1 = phi_map(iphi,2)                                                              
     l2 = phi_map(iphi,3)                                                              
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l1), Q(l1+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l2), Q(l2+npdiffq), tmp2r, tmp2i )
     sum_phi(l)         = sum_phi(l)         - 0.5d0 * tmp2r
     sum_phi(l+npdiffq) = sum_phi(l+npdiffq) - 0.5d0 * tmp2i
  end do                                                                               
  !$omp end parallel do 
  ! --- sub-case B) -> l == l1 /= l2 (in this sub-case we need to update sum_phi for l and l2)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(20), jl2(20)                                                             
     l = phi_map(iphi,1)                                                               
     l1 = phi_map(iphi,2)                                                              
     l2 = phi_map(iphi,3)                                                              
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l1), Q(l1+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l2), Q(l2+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          -         tmp2r
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  -         tmp2i
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l2)         = sum_phi(l2)         - 0.5d0 * tmp2r
     sum_phi(l2+npdiffq) = sum_phi(l2+npdiffq) - 0.5d0 * tmp2i
  end do
  !$omp end parallel do 
  ! --- sub-case C) -> l /= l1 == l2 (need to update sum_phi for l and l1; in the latter, one permutation is further possible)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(21), jl2(21) 
     l = phi_map(iphi,1)
     l1 = phi_map(iphi,2)
     l2 = phi_map(iphi,3)
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - 0.5d0 * tmp2r
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  - 0.5d0 * tmp2i
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         -         tmp2r
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) -         tmp2i
  end do
  !$omp end parallel do 
  ! --- sub-case D) -> l /= l1 /= l2 (need to update sum_phi for all l, l1 and l2)
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(22), jl2(22)
     l = phi_map(iphi,1)
     l1 = phi_map(iphi,2)
     l2 = phi_map(iphi,3)
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2), Q(l2+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - tmp2r
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  - tmp2i
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         - tmp2r
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - tmp2i
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l), Q(l+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l2)         = sum_phi(l2)         - tmp2r
     sum_phi(l2+npdiffq) = sum_phi(l2+npdiffq) - tmp2i
  end do
  !$omp end parallel do 
  ! end case XLIII)
  
  ! case XLIV) -> l in S(q), l1 in S(q), l2 in S*(q)
  ! this covers also:
  ! case XLVII) -> l in S(q), l1 in S*(q), l2 in S(q) via permutation
  ! case XLVIII) -> l in S(q), l1 in S*(q), l2 in S*(q) via permutation + complex-conjugation 
  ! --- sub-case A) -> l == l1
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(23), jl2(23)                                                                                                        
     l = phi_map(iphi,1)                                                                                                            
     l1 = phi_map(iphi,2)                                                                                                           
     l2 = phi_map(iphi,3)
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l1), Q(l1+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l2-npdiffq), -Q(l2), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          -         tmp2r ! XLIV)+XLVII)
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  -         tmp2i ! XLIV)+XLVII)
     call complex_product ( tmp1r, -tmp1i, Q(l), -Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l2-npdiffq) = sum_phi(l2-npdiffq) - 0.5d0 * tmp2r ! XLVIII)
     sum_phi(l2)         = sum_phi(l2)         - 0.5d0 * tmp2i ! XLVIII)
  end do
  !$omp end parallel do 
  ! --- sub-case B) -> l /= l1
  !$omp parallel do default(none) shared(jl1,jl2,phi_map,phi,Q,npdiffq) private(iphi,l,l1,l2,phi_tmp,tmp1r,tmp1i,tmp2r,tmp2i) reduction(+:sum_phi) 
  do iphi = jl1(24), jl2(24)                                                                                                        
     l = phi_map(iphi,1)                                                                                                            
     l1 = phi_map(iphi,2)                                                                                                           
     l2 = phi_map(iphi,3)                                                                                                           
     phi_tmp(1:2) = phi(iphi,1:2)
     call complex_product ( phi_tmp(1), phi_tmp(2), Q(l2-npdiffq), -Q(l2), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l)          = sum_phi(l)          - tmp2r ! XLIV)+XLVII)          
     sum_phi(l+npdiffq)  = sum_phi(l+npdiffq)  - tmp2i ! XLIV)+XLVII)  
     call complex_product ( tmp1r, tmp1i, Q(l), Q(l+npdiffq), tmp2r, tmp2i )
     sum_phi(l1)         = sum_phi(l1)         - tmp2r ! XLIV)+XLVII)           
     sum_phi(l1+npdiffq) = sum_phi(l1+npdiffq) - tmp2i ! XLIV)+XLVII) 
     call complex_product ( phi_tmp(1), -phi_tmp(2), Q(l), -Q(l+npdiffq), tmp1r, tmp1i )
     call complex_product ( tmp1r, tmp1i, Q(l1), -Q(l1+npdiffq), tmp2r, tmp2i )
     sum_phi(l2-npdiffq) = sum_phi(l2-npdiffq) - tmp2r ! XLVIII) 
     sum_phi(l2)         = sum_phi(l2)         - tmp2i ! XLVIII)        
  end do                                                                                                                            
  !$omp end parallel do 
  
  ! take the complex-conjugate of the final sum
  do l = npGq+1, npuniqueq
     sum_phi(l+npdiffq) = -sum_phi(l+npdiffq)
  end do
  
  return
end subroutine calc_double_sum_phi

! This subroutine computes the product between two complex numbers
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine complex_product ( a_real, a_imag, b_real, b_imag, prod_real, prod_imag )
  implicit none
  real(8), intent(in) :: a_real, a_imag, b_real, b_imag
  real(8), intent(out) :: prod_real, prod_imag

  prod_real = a_real * b_real - a_imag * b_imag

  prod_imag = a_real * b_imag + a_imag * b_real
  
  return
end subroutine complex_product

