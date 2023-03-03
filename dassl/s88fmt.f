      SUBROUTINE S88FMT(N,IVALUE,IFMT)
C***begin prologue  s88fmt
C***refer to  xerror
C     abstract
C        s88fmt replaces ifmt(1), ... ,ifmt(n) with the
C        characters corresponding to the n least significant
C        digits of ivalue.
C
C     taken from the bell laboratories port library error handler
C     latest revision ---  7 june 1978
C
C***references
C   jones r.e., *slatec common mathematical library error handling
C    package*, sand78-1189, sandia laboratories, 1978.
C***routines called  (none)
C***end prologue  s88fmt
C
      DIMENSION IFMT(N),IDIGIT(10)
      DATA IDIGIT(1),IDIGIT(2),IDIGIT(3),IDIGIT(4),IDIGIT(5),
     1     IDIGIT(6),IDIGIT(7),IDIGIT(8),IDIGIT(9),IDIGIT(10)
     2     /1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/
C***first executable statement  s88fmt
      NT = N
      IT = IVALUE
   10    IF (NT .EQ. 0) RETURN
         INDEX = MOD(IT,10)
         IFMT(NT) = IDIGIT(INDEX+1)
         IT = IT/10
         NT = NT - 1
         GO TO 10
      END
      DOUBLE PRECISION FUNCTION D1MACH(IDUM)
      INTEGER IDUM
      DOUBLE PRECISION U,COMP
      U = 1.0D0
 10   U = U * 0.5D0
      COMP  = 1.0D0 + U
      IF(COMP .NE. 1.0D0) GOTO 10
      D1MACH = U * 2.0D0
      RETURN
      END