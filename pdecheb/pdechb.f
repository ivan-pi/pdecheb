      SUBROUTINE PDECHB(T,U,UDOT,RESD,IRES,WK,IWK)
C***********************************************************************
C
C THIS IS THE CHEBYSHEV GLOBAL ELEMENT ROUTINE TO EVALUATE THE
C RESIDUAL OF THE IMPLICIT SET OF O.D.E.'S DEFINED BY
C
C        RESIDUAL  =  A(U,T)*DU/DT  -  F(U,T)
C
C PARAMETER LIST
C----------------
C  T               CURRENT TIME INTEGRATION LEVEL , > 0.0
C  U(N)            CURRENT SOLUTION VECTOR
C  RESD(N)          VECTOR WHICH WILL CONTAIN THE RESIDUAL ON EXIT
C  UDOT(N)         CURRENT ESTIMATE OF DU/DT
C  WK(1)           REAL WORKSPACE - DEFINED IN INICHB
C  IWK(1)          INTEGER WORKSPACE - NOT USED HERE.
C  IRES            INDICATOR FOR DASSL FROM RESIDUAL ROUTINE.
C                  ON EXIT = -1 THEN ILLEGAL SOLUTION VALUES HAVE BEEN
C                               FOUND .
C                           =-2 DASSL SHOULD HALT THE INTEGRATION.
C
C  ONLY  RESD(N) IS ALTERED ON EXIT : IT CONTAINS THE CURRENT RESIDUAL
C***********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IRES
C     .. Array Arguments ..
      DOUBLE PRECISION  RESD(*), U(*), UDOT(*), WK(*)
      INTEGER           IWK(*)
C     .. Scalars in Common ..
      INTEGER           I10, I10A, I10B, I11, I11A, I11B, I12, I13, I14,
     *                  I15, I16, I17, I18, I19, I2, I3, I4, I5, I6, I7,
     *                  I8, I9, M, NEL, NPDE, NPTL, NPTS, NV, NVST, NXI
      CHARACTER*6       PDCODE
C     .. Local Scalars ..
      INTEGER           I, IBK, IFL, IR, ITYPE, IV, J, N
      CHARACTER*240     ERRMSG
C     .. External Subroutines ..
      EXTERNAL          CHINTR, CRES, DRES, SCHERR, SODEFN
C     .. Common blocks ..
      COMMON            /DISCHK/PDCODE
      COMMON            /SCHSZ/I2, I3, I4, I5, I6, I7, I8, I9, I10,
     *                  I10A, I10B, I11, I11A, I11B, I12, I13, I14, I15,
     *                  I16, I17, I18, I19
      COMMON            /SCHSZ1/NEL, NPTL, NPDE, NPTS, M, NV, NXI, NVST
C     .. Save statement ..
      SAVE              /SCHSZ1/, /SCHSZ/, /DISCHK/
C     .. Executable Statements ..
C
      IF (PDCODE.NE.'C0CHEB') THEN
         ERRMSG =
     *' C0CHEB-RES ROUTINE ERROR-THE SETUP ROUTINE INICHB       WAS NOT
     *CALLED BEFORE DASSL  WAS ENTERED'
         CALL SCHERR(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
         IRES = -2
         RETURN
      END IF
C
      IR = 1
      IRES = 1
      N = NPDE*NPTS + NV
      DO 20 J = 1, N
         RESD(J) = 0.0D0
   20 CONTINUE
      IBK = NEL + 1
      IV = NPTS*NPDE
      IF (NV.GT.0) THEN
         IV = NVST
C        GENERATE THE SOLUTION VALUES SPACE DERIVS AND FLUXES AT THE
C        COUPLING POINTS
         ITYPE = 3
         IFL = 0
         CALL CHINTR(NXI,WK(I18),WK(I13),ITYPE,U,NPTS,NPDE,NEL,NPTL,WK,
     *               WK(I10),WK(I5),IBK,IFL,NV,U(IV),UDOT(IV),WK(I11),T,
     *               IR)
         IF (IR.NE.1 .OR. IFL.EQ.1) GO TO 60
C        GENERATE TIME DERIV VALUES AND THEIR SPACE DERIVS AT THE
C        COUPLING POINTS.
         ITYPE = 2
         CALL CHINTR(NXI,WK(I18),WK(I16),ITYPE,UDOT,NPTS,NPDE,NEL,NPTL,
     *               WK,WK(I10),WK(I5),IBK,IFL,NV,U(IV),UDOT(IV),WK(I11)
     *               ,T,IR)
         IF (IR.NE.1 .OR. IFL.EQ.1) GO TO 60
C        CALL THE ROUTINE TO DEFINE THE AUXILLARY ODE RESIDUAL.
         CALL SODEFN(T,NV,U(IV),UDOT(IV),NPDE,NXI,WK(I18),WK(I13),
     *               WK(I14),WK(I15),WK(I16),WK(I17),RESD(IV),IRES)
         IF (IRES.NE.1) GO TO 60
      END IF
C       CALL THE CO COLLOCATION DISCRETISATION ROUTINE
      IR = 1
      IF (NPTL.GT.2) THEN
C        GENERAL POLYNOMIAL VERSION.
         CALL CRES(NPDE,NPTS,T,U,RESD,UDOT,M,WK(I6),WK,WK(I2),WK(I5),
     *             WK(I7),WK(I8),WK(I9),WK(I10),WK(I11),NEL,NPTL,WK(I4),
     *             WK(I12),IRES,WK(I10A),WK(I11A),WK(I11B),WK(I10B),NV,
     *             U(IV),UDOT(IV),WK(I19))
      ELSE
C        LINEAR BASIS FUNCTION VERSION.
         CALL DRES(NPDE,NPTS,T,U,RESD,UDOT,M,WK(I6),WK,WK(I2),WK(I5),
     *             WK(I7),WK(I8),WK(I9),WK(I10),WK(I11),NEL,NPTL,WK(I4),
     *             WK(I12),IRES,WK(I10A),WK(I11A),WK(I11B),WK(I10B),NV,
     *             U(IV),UDOT(IV),WK(I19))
      END IF
      DO 40 I = 1, N
         RESD(I) = -RESD(I)
   40 CONTINUE
      IF (IRES.NE.1) THEN
         IR = IRES
         GO TO 60
      END IF
      RETURN
   60 IRES = IR
      IF (IR.EQ.-2) THEN
         ERRMSG =
     *' ROUTINE PDECHB AT TIME T (=R1). THE VALUE OF IRES
     * HAS BEEN SET TO -2 TO TERMINATE INTEGRATION.'
         CALL SCHERR(ERRMSG,1,0,0,0,1,T,0.0D0)
      ELSE IF (IR.NE.-1) THEN
         ERRMSG =
     *' ROUTINE PDECHB AT TIME T (=R1). THE
     *  VALUE OF IRES HAS BEEN SET TO AN ILLEGAL VALUE (=I1).
     *PDECHB HAS RESET IRES TO -1 AND INTEGRATION CONTINUES.'
         CALL SCHERR(ERRMSG,1,0,0,0,1,T,0.0D0)
         IRES = -1
      END IF
      RETURN
C
C---------------------------END OF PDECHB-----------------------------
C
      END