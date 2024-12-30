      SUBROUTINE DDAWTS(NEQ,IWT,RTOL,ATOL,Y,WT,RPAR,IPAR)
C
C***BEGIN PROLOGUE  DDAWTS
C***REFER TO  DDASSL
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   830315   (YYMMDD)
C***REVISION DATE  830315   (YYMMDD)
C***END PROLOGUE  DDAWTS
C-----------------------------------------------------------------------
C     this subroutine sets the error weight vector
C     wt according to wt(i)=rtol(i)*abs(y(i))+atol(i),
C     i=1,-,n.
C     rtol and atol are scalars if iwt = 0,
C     and vectors if iwt = 1.
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
      INTEGER NEQ, IWT
      REAL(DP) RTOL(*), ATOL(*), Y(*), WT(*), RPAR(*)
      INTEGER IPAR(*)
C     .. Local Scalars ..
      REAL(DP) RTOLI, ATOLI
      INTEGER I
C     .. Executable Statements ..
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO I=1,NEQ
         IF (IWT == 1) THEN
           RTOLI=RTOL(I)
           ATOLI=ATOL(I)
         END IF
         WT(I)=RTOLI*DABS(Y(I))+ATOLI
      END DO
      RETURN
C-----------end of subroutine ddawts------------------------------------
      END