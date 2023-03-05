MODULE PDECHEB_COMMON
IMPLICIT NONE
INTEGER, PARAMETER :: DP = KIND(1.0D0)

! COMMON /DISCHK/
CHARACTER(LEN=6) :: PDCODE

! COMMON /SCHSZ/
INTEGER :: I2, I3, I4, I5, I6, I7, I8, I9, I10, &
           I10A, I10B, I11, I11A, I11B, I12, I13, I14, I15, &
           I16, I17, I18, I19

! Aliases
INTEGER :: IA(3), IB(3), IC(2), ID(9)

EQUIVALENCE (IA(1),I2), (IA(2),I3), (IA(3),I4)
EQUIVALENCE (IB(1),I6), (IB(2),I7), (IB(3),I8)
EQUIVALENCE (IC(1),I10A), (IC(2),I10B)

EQUIVALENCE (ID(1),I11A)
EQUIVALENCE (ID(2),I11B)
EQUIVALENCE (ID(3),I12)
EQUIVALENCE (ID(4),I13)
EQUIVALENCE (ID(5),I14)
EQUIVALENCE (ID(6),I15)
EQUIVALENCE (ID(7),I16)
EQUIVALENCE (ID(8),I17)
EQUIVALENCE (ID(9),I18)

! COMMON /SCHSZ1/
INTEGER :: NEL, NPTL, NPDE, NPTS, M, NV, NXI, NVST
! NEL ... number of elements
! NPTL = NPOLY + 1
! NPDE ... The number of P.D.E.S
! M = 0,1,2 ... Cartesian, cylindrical or spherical polar
!               co-ordinates in use
! NV ... number of auxillary ODEs
! NXI ... the number of coupling points
! nvst ... the starting point of the ODE components
!          in the solution vector

! COMMON /SCHSZ2/
INTEGER :: IDEV
   !> Output channel for error messages

! COMMON /SCHSZ3/
REAL(DP) :: TWOU
   !> Estimate of unit round-off error
   !     See INICHB for details how this is calculated
   !     TODO: Replace with intrinsic function epsilon()

! COMMON /SCHSZ4
REAL(DP) :: TO
INTEGER :: K1, K2, K3, K4, JTIMES, ILOC

! COMMON /SCHSZ5/
INTEGER :: NNNPTL

! COMMON /SCHSZ6/
REAL(DP) :: CCRULE(50)
   !> Clenshaw-Curtis Quadrature Rule

END MODULE

!
! Eliminate common blocks in PDECHEB:
! - [x] chintr (SCHSZ3)
! - [x] cres (SCHSZ3)
! - [x] cset (SCHSZ6)
! - [x] dres (SCHSZ3)
! - [x] error (SCHSZ2)
! - [x] inichb (DISCHK,SCHSZ,SCHSZ1,SCHSZ2,SCHSZ3,SCHSZ4,SCHSZ5)
! - [x] interc (DISCHK,SCHSZ,SCHSZ1,SCHSZ3)
! - [x] intrch (SCHSZ2)
! - [x] pdecheb (DISCHK,SCHSZ,SCHSZ1)
! - [x] scherr (SCHSZ2)