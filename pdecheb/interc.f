      SUBROUTINE INTERC(XP,UP,NP,U,NEQ,NPDE,IFLAG,ITYPE,WK,IWK)
C********************************************************************
C
C        SPACE INTERPOLATION ROUTINE FOR POST-PROCESSING OF SOLUTION
C        PRODUCED BY DASSL.
C        THIS ROUTINES PROVIDES VALUES OF THE SOLUTION AND POSSIBLY THE
C        FIRST DERIV IN SPACE AND THE FLUX ON THE MESH XP(NP).
C
C        PARAMETERS
C       --------------
C        NPDE     ON ENTRY MUST CONTAIN NO OF PARABOLIC EQUATIONS
C        NPTS     ON ENTRY MUST CONTAIN THE NUMBER OF SPATIAL
C                 MESH POINTS USED IN TIME INTEGRATION.
C        NP       ON ENTRY MUST CONTAIN THE NUMBER OF SPATIAL
C                 INTERPOLATION POINTS
C        XP(NP)   ARRAY WHICH ON ENTRY
C                 CONTAINS THE SPATIAL INTERPOLATION POINTS
C                 WE ASSUME THAT
C                    XP(I) <  XP(I+1)  ,  I = 1,...,NP-1
C        UP(NPDE,NP,ITYPE)  EMPTY ARRAY FOR THE INTERPOLATED VALUES AT
C                           THE CURRENT TIME LEVEL. THE VALUES OF THIS
C                           ARRAY ON EXIT DEPEND ON THE PARAMETER ITYPE.
C        U(NPDE,NPTS) THE CURRENT SOLUTION VECTOR COMPUTED BY THE ODE
C                  TIME INTEGRATOR MUST BE SUPPLIED IN THIS VECTOR.
C        IFLAG          ERROR FLAG       = 0 ON SUCCESSFUL RETURN
C                                        = 1 IF EXTRAPOLATION TRIED.
C                                        = 2 IF WORKSPACE NOT INITIAL
C                                               ISED ON ENTRY BY INICHB.
C                                        = 3 ILLEGAL VALUE OF ITYPE.
C        ITYPE = 1  ONLY THE SOLUTION IS OUTPUT IN THE ARRAY UP
C                   UP(J,K,1) HOLDS U(XP(K),T) FOR PDE J
C                2  AS FOR 1 BUT THE FIRST DERIV IS ALSO OUTPUT.
C                   UP(J,K,2) HOLDS D/DX U(XP(K),T).
C
C        WK(IWK) THE WORKSPACE USED BY THE CHEBYSHEV METHOD. THIS
C                MUST BE THE WORKSPACE INITIALISED BY INICHB.
C**********************************************************************
C     .. Scalar Arguments ..
      INTEGER           IFLAG, ITYPE, IWK, NEQ, NP, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NEQ), UP(NPDE,NP,ITYPE), WK(IWK), XP(NP)
C     .. Scalars in Common ..
      DOUBLE PRECISION  TWOU
      INTEGER           I10, I11, I19, I5, I9, MM, NEL, NNPDE, NNPTS,
     *                  NPTL, NV, NVST, NXI
      CHARACTER*6       PDCODE
C     .. Arrays in Common ..
      INTEGER           IA(3), IB(3), IC(2), ID(9)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, IBK, J, K, NPTS
      CHARACTER*240     ERRMSG
C     .. External Subroutines ..
      EXTERNAL          INTRCH, SCHERR
C     .. Intrinsic Functions ..
      INTRINSIC         DABS
C     .. Common blocks ..
      COMMON            /DISCHK/PDCODE
      COMMON            /SCHSZ/IA, I5, IB, I9, I10, IC, I11, ID, I19
      COMMON            /SCHSZ1/NEL, NPTL, NNPDE, NNPTS, MM, NV, NXI,
     *                  NVST
      COMMON            /SCHSZ3/TWOU
C     .. Save statement ..
      SAVE              /SCHSZ1/, /SCHSZ/, /DISCHK/, /SCHSZ3/
C     .. Executable Statements ..
      IF (PDCODE.NE.'C0CHEB') THEN
         IFLAG = 1
         GO TO 80
      END IF
      IFLAG = 0
      IBK   = NEL + 1
      IF (ITYPE.NE.1 .AND. ITYPE.NE.2) THEN
         ERRMSG =
     *' ILLEGAL VALUE OF ITYPE IN CALL TO SUBROUTINE INTERC
     *  THE VALUE IS (=I1), BUT SHOULD BE 1 OR 2 '
         CALL SCHERR(ERRMSG,1,1,ITYPE,0,0,0.0D0,0.0D0)
         IFLAG = 3
         GO TO 80
      END IF
C
C   TEST THE INTERPOLATION POINTS XP(NP) TO ENSURE THAT THEY ARE IN
C   INCREASING ORDER AND THAT IF ITYPE = 2 (DERIVATIVES REQUIRED) THE
C   POINTS DO NOT CONFLICT WITH THE BREAK-POINTS.
C
      DO 20 I = 2, NP
         TEMP = XP(I) - XP(I-1)
         IF (TEMP.LE.TWOU) THEN
            ERRMSG =
     *' INTERC ROUTINE CALLED WITH INTERP.POINTS NOT            IN STRIC
     *TLY INCREASING ORDER I.E. COMPONENT NO (=I1)              WITH VAL
     *UE (=R1) IS GREATER THAN COMPONENT( =I2)                  WITH VAL
     *UE (=R2).'
            CALL SCHERR(ERRMSG,1,2,J,I,2,XP(J),XP(I))
         END IF
   20 CONTINUE
      IF (ITYPE.GE.2 .AND. IBK.GT.2) THEN
         DO 60 I = 1, NP
            DO 40 J = 2, NEL
               TEMP = DABS(XP(I)-WK(I5-1+J))
               IF (TEMP.LE.TWOU) THEN
                  K = I5 + J - 1
                  ERRMSG =
     *' INTERC ROUTINE CALLED WITH ITYPE = 2                    AND INTE
     *RP. POINTS EQUAL TO BREAK-POINTS I.E.                     COMPONEN
     *T NO (=I1) WITH VALUE (=R1)                               IS CLOSE
     * TO BREAK POINT(=I2) WITH VALUE (=R2).'
                  CALL SCHERR(ERRMSG,1,2,I,J,2,XP(I),WK(K))
               END IF
   40       CONTINUE
   60    CONTINUE
      END IF
C
C    CALL THE INTERPOLATION ROUTINE.
C
      NPTS = NNPTS
      CALL INTRCH(NP,XP,UP,ITYPE,U,NPTS,NPDE,NEL,NPTL,WK,WK(I10),WK(I5),
     *            IBK,IFLAG)
   80 CONTINUE
      RETURN
      END
C***********************************************************************
C