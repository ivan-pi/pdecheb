      SUBROUTINE SCHERR(MSG,IERT,NI,I1,I2,NR,R1,R2)
C-----------------------------------------------------------------------
C  ERROR HANDLING ROUTINE FOR THE DASSL  INTEGRATION PACKAGE. THIS
C  ROUTINE IS A FORTRAN77 IMPROVED VERSION OF THE ROUTINE USED IN LSODI
C  AND MAKES USE OF CHARACTER HANDLING FACILITIES.
C-----------------------------------------------------------------------
      USE PDECHEB_COMMON, ONLY: NERR=>IDEV
C     .. Scalar Arguments ..
      DOUBLE PRECISION  R1, R2
      INTEGER           I1, I2, IERT, NI, NR
      CHARACTER*(*)     MSG
C     .. Local Scalars ..
      INTEGER           I, IL, IT, J, K, KP1, LWORD
      CHARACTER(240)    MSG1
C     .. Local Arrays ..
      CHARACTER(60)     MSGOUT(5)
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C ALL ARGUMENTS ARE INPUT ARGUMENTS.
C
C MSG    = THE MESSAGE IN CHARACTER FORMAT
C IERT   = THE ERROR TYPE..
C          1 MEANS RECOVERABLE (CONTROL RETURNS TO CALLER).
C          2 MEANS FATAL (RUN IS ABORTED--SEE NOTE BELOW).
C NI     = NUMBER OF INTEGERS (0, 1, OR 2) TO BE PRINTED WITH MESSAGE.
C I1,I2  = INTEGERS TO BE PRINTED, DEPENDING ON NI.
C NR     = NUMBER OF REALS (0, 1, OR 2) TO BE PRINTED WITH MESSAGE.
C R1,R2  = REALS TO BE PRINTED, DEPENDING ON NR.
C-----------------------------------------------------------------------
      IL = LEN(MSG)
C
C    SET MSG1 BLANK AND GET RID OF UNNECESSARY SPACES IN ERROR MESSAGE
C
      J = 1
      IT = MIN(IL,240)
      DO 20 I = 1, 10
         MSG1(J:) = '                        '
         J = J + 24
   20 CONTINUE
      K = 0
      J = 0
      DO 40 I = 1, IT
         IF (MSG(I:I).EQ.' ') THEN
            K = K + 1
            IF (K.GT.2) GO TO 40
         ELSE
            K = 0
         END IF
         J = J + 1
         MSG1(J:J) = MSG(I:I)
   40 CONTINUE
      IL = J
C
C     FORMAT THE MESSAGE NOW STORED IN MSG1
C
      I = 1
      LWORD = 60
      J = 0
   60 J = J + 1
      IF (J.GT.1) LWORD = 51
      K = I + LWORD - 1
      KP1 = K + 1
   80 IF (MSG1(K:K).NE.' ' .AND. MSG1(KP1:KP1).NE.' ') THEN
         K = K - 1
         IF (K.EQ.I) THEN
            K = I + LWORD - 1
            GO TO 100
         END IF
         GO TO 80
      END IF
  100 IF (J.EQ.1) THEN
         MSGOUT(J) = MSG1(I:K)
      ELSE
         MSGOUT(J) = '         '//MSG1(I:K)
      END IF
      I = K + 1
      IF (K.LT.IL .AND. J.LT.5) GO TO 60
C
C  OUTPUT THE ERROR MESSAGE
C
      WRITE (NERR,FMT=99999) (MSGOUT(I),I=1,J)
C
C  PRINT THE INTEGERS AND REALS IN THE ERROR MESSAGE (IF ANY)
C
      IF (NI.EQ.1) WRITE (NERR,FMT=99998) I1
      IF (NI.EQ.2) WRITE (NERR,FMT=99997) I1, I2
      IF (NR.EQ.1) WRITE (NERR,FMT=99996) R1
      IF (NR.EQ.2) WRITE (NERR,FMT=99995) R1, R2
C ABORT THE RUN IF IERT = 2. -------------------------------------------
      IF (IERT.NE.2) RETURN
      STOP
C----------------------- END OF SUBROUTINE SCHERR ----------------------
99999 FORMAT (1X,A60)
99998 FORMAT (9X,' IN ABOVE MESSAGE I1 =',I10)
99997 FORMAT (9X,' IN ABOVE MESSAGE I1 =',I10,'   I2 =',I10)
99996 FORMAT (9X,' IN ABOVE MESSAGE R1 =',D21.13)
99995 FORMAT (9X,'IN ABOVE,  R1 =',D21.13,3X,'R2 =',D21.13)
      END
C**********************************************************************
C