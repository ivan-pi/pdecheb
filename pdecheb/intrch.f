      SUBROUTINE INTRCH(NP,XP,UP,ITYPE,U,NPTS,NPDE,NEL,NPTL,OMEGA,COEFF,
     *                  XBK,IBK,IFLAG)
C
C***********************************************************************
C         PARAMETER LIST
C         **************
C         XP(NP)         THE MESH POINTS AT WHICH INTERPOLATED VALUES
C                        ARE REQUIRED. THESE POINTS SUCH BE IN
C                        INCREASING ORDER.
C         UP(NPDE,NP,ITYPE)  ARRAY THAT HOLDS THE VALUES FOUND BY
C                            INTERPOLATION.
C         IF ITYPE >= 1  UP(J,K,1) HOLDS THE SOLUTION VALUE AT MESH
C                                  POINT XP(K) FOR JTH PDE
C         IF ITYPE >= 2  UP(J,K,2) HOLDS THE SPACE DERIV OF THE SOLUTION
C                                  AT POINT XP(K) FOR JTH PDE.
C
C         U(1..NEQN)     ORIGINAL SOLUTION VECTOR FROM THE ODE CODE.
C
C         NPTS           THE NUMBER OF MESH POINTS USED IN COMPUTING U.
C         NPDE           THE NUMBER OF PDES IN THE PROBLEM.
C         NEL            THE NUMBER OF SPATIAL ELEMENTS IN THE MESH.
C         NPTL           THE NUMBER OF MESH POINTS PER ELEMENT.
C                        THEREFORE NPTS = NEL*(NPTL-1) + 1
C         OMEGA          MATRIX USED IN MAPPING FROM THE SOLUTION ON A
C                        SPATIAL INTERVAL TO ITS CHEBYSHEV COEFFS.
C         COEFFS         WORKSPACE USED TO HOLD THESE COEFFS.
C         XBK(IBK)       ARRAY USED TO HOLD THE BREAKPOINTS BETWEEN THE
C                        SPATIAL ELEMENTS.
C         IFLAG          ERROR FLAG SET TO 0 UNLESS EXTRAPOLATION IS
C                        TRIED AND THEN SET TO 1.
C
C                THE METHOD USED IS DECOMPOSITION OF THE SOLUTION
C         PER ELEMENT INTO CHEBYSHEV COEFFICIENTS. THIS IS DONE BY
C         MATRIX MULTIPLICATION USING THE OMEGA MATRIX .  F.F.T.
C         COULD ALSO BE USED. INTERPOLATION IS USED TO PROVIDE THOSE
C         SOLUTION VALUES IN THE ELEMENT (USING CLENSHAWS ALGORITHM).
C
C***********************************************************************
C COMMON /SCHSZ2/
      USE PDECHEB_COMMON, ONLY: TWOU
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           IBK, IFLAG, ITYPE, NEL, NP, NPDE, NPTL, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  COEFF(NPDE,NPTL,2), OMEGA(NPTL,NPTL),
     *                  U(NPDE,NPTS), UP(NPDE,NP,*), XBK(IBK), XP(NP)
C     .. Local Scalars ..
      DOUBLE PRECISION  AL, BR, BR1, BR2, CU, TEM, TEM1
      INTEGER           I, II, IP, IP1, IX, IY, IZ, J, K, NM1
      CHARACTER(240)    ERRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  XCON(2)
C     .. External Subroutines ..
      EXTERNAL          SCHERR
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C  TREAT EACH ELEMENT SEPARATELY
C
      CU = TWOU
      TEM1 = 1 - CU
      TEM = 1 + CU
      IP = 0
      NM1 = NPTL - 1
      IZ = 0
      OUTER: DO I = 1, NEL
         IP1 = I + 1
         IF (XBK(I).GT.(XBK(I+1)*TEM1-CU)) THEN
            ERRMSG =
     *' INTERC ROUTINE BREAKPOINT NUMBER (=I1)                  WITH VAL
     *UE (=R1) IS TOO CLOSE OR LARGER THAN BREAKPOINT NO        (=I2) WI
     *TH VALUE (=R2). INCORRECT CALL TO INTERC ASSUMED OR       WORKSPAC
     *E CORRUPTED'
            CALL SCHERR(ERRMSG,1,2,I,IP1,2,XBK(I),XBK(IP1))
            GO TO 260
         END IF
   20    IP = IP + 1
         IF (IP.EQ.(NP+1)) GO TO 260
         IF (XP(IP).LT.(XBK(I)*TEM1-CU)) GO TO 20
         IF (XP(IP).GT.(XBK(I+1)*TEM+CU)) THEN
             IP = IP - 1
             CYCLE OUTER
         END IF
         IF (XP(IP).GE.(XBK(I+1)*TEM1-CU)) THEN
            IF (I.LT.NEL .AND. ITYPE.GE.2) IZ = 1
C           IZ = 1 MEANS THAT WEIGHTED AVERAGE MUST BE USED FOR
C           DERIVATIVE VALUES THAT ARE REQUESTED AT XBK(I+1)
         END IF
C        ***************************************************************
C         PROCESS A SEQUENCE OF XP(J) VALUES IN ELEMENT I
C         IX = START OF CORRECT PART OF SOLUTION VECTOR U
C         FORM THE CHEBYSHEV COEFFS IN THE ARRAY COEFF.
C        ***************************************************************
         IX = NM1*(I-1)
         DO K = 1, NPDE
            DO J = 1, NPTL
               COEFF(K,J,1) = 0.0D0
               DO II = 1, NPTL
                  COEFF(K,J,1) = COEFF(K,J,1) + OMEGA(J,II)*U(K,IX+II)
               END DO
            END DO
         END DO
C        FORM THE CHEBYSHEV COEFFS OF THE SPACE DERIV.
         IF (ITYPE.GE.2) THEN
            DO K = 1, NPDE
               COEFF(K,NPTL,2) = 0.0D0
               COEFF(K,NPTL-1,2) = 2.0D0*NM1*COEFF(K,NPTL,1)
               DO J = 2, NM1
                  COEFF(K,NPTL-J,2) = COEFF(K,NPTL-J+2,2) + COEFF(K,
     *                                NPTL-J+1,1)*2*(NPTL-J)
               END DO
               COEFF(K,1,2) = COEFF(K,1,2)*0.5D0
            END DO
         END IF
         XCON(1) = 2.0D0/(XBK(I+1)-XBK(I))
         XCON(2) = -0.5D0*XCON(1)*(XBK(I+1)+XBK(I))
         IY = MIN(2,ITYPE)
  140    DO II = 1, IY
            DO K = 1, NPDE
               BR1 = 0.0D0
               BR2 = 0.0D0
C               COEFF(K,NPTL) IS THE NPTL-TH  COEFF OF SOLUTION OF PDE K
               AL = (XP(IP)*XCON(1)+XCON(2))*2.0D0
               BR = COEFF(K,NPTL,II)
               DO J = 1, NM1
                  BR2 = COEFF(K,NPTL-J,II) + AL*BR - BR1
                  BR1 = BR
                  BR = BR2
               END DO
               IF (II.EQ.1) THEN
                  UP(K,IP,II) = BR - BR1*AL*0.5D0
               ELSE IF (IZ.LT.2) THEN
                  UP(K,IP,II) = (BR-BR1*AL*0.5)*XCON(1)
               ELSE
                  UP(K,IP,II) = 1.D0/(XBK(I+1)-XBK(I-1))*(UP(K,IP,II)
     *                          *(XBK(I)-XBK(I-1))+(BR-BR1*AL*0.5)
     *                          *XCON(1)*(XBK(I+1)-XBK(I)))
               END IF
            END DO
         END DO
         IF (IP.EQ.NP) CYCLE OUTER
         IP = IP + 1
         IF (IZ.EQ.1) THEN
            IZ = 2
            GO TO 220
C           TO CALCULATE THE OTHER ELEMENTS CONTRIBUTION TO DERIV.
         END IF
         IF (IZ.EQ.2) IZ = 0
         IF (XP(IP).LE.(XBK(I+1)*TEM1-CU)) THEN
            GO TO 140
         ELSE IF (XP(IP).LE.(XBK(I+1)*TEM+CU)) THEN
            IF (I.LT.NEL .AND. ITYPE.GE.2) IZ = 1
C           IZ = 1 MEANS THAT WEIGHTED AVERAGE MUST BE USED FOR
C           DERIVATIVE VALUES THAT ARE REQUESTED AT XBK(I+1)
            GO TO 140
         END IF
  220    IP = IP - 2
      END DO OUTER
      RETURN
  260 IFLAG = 1
      RETURN
C---------END OF INTRCH-------------------------------------------
C
      END