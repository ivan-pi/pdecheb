      SUBROUTINE DRES(NPDE,NPTS,T,U,RES,UDOT,M,X,OMEGA,DU,XBK,BETA,
     *                GAMMA,DUDX,R,Q,NEL,NPTL,XC,CCR,IRES,RT,QT,UDT,
     *                UTDX,NV,V,VDOT,VDUM)
C**********************************************************************
C       CHEBYSHEV C0 COLLOCATION ROUTINE
C       THIS VERSION FOR USE WITH LINEAR BASIS FUNCTIONS ONLY
C**********************************************************************
      USE PDECHEB_COMMON, ONLY: TWOU
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION T
      INTEGER         IRES, M, NEL, NPDE, NPTL, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION BETA(NPDE,4), CCR(NPTL), DU(NPTL,NPTL),
     *                DUDX(NPDE,NPTL), GAMMA(NPDE,4), OMEGA(NPTL,NPTL),
     *                Q(NPDE,NPTL), QT(NPDE,NPTL), R(NPDE,NPTL),
     *                RES(NPDE,NPTS), RT(NPDE,NPTL), U(NPDE,NPTS),
     *                UDOT(NPDE,NPTS), UDT(NPDE,NPTL), UTDX(NPDE,NPTL),
     *                V(*), VDOT(*), VDUM(*), X(NPTS), XBK(*), XC(NPTL)
C     .. Local Scalars ..
      DOUBLE PRECISION H, MP1, SAVEL, SAVER, SFIRST, TEM
      INTEGER         I, II, IJ, IK, IV, J, JJ, JK, K, KJ, NM1
C     .. Local Arrays ..
      INTEGER         IZ(3)
C     .. External Subroutines ..
      EXTERNAL        SBNDR, SPDEFN
C     .. Intrinsic Functions ..
      INTRINSIC       MAX, MIN
C     .. Executable Statements ..
      NM1 = NPTL - 1
      IV = MAX(1,NV)
      MP1 = 1.0D0
      DO I = 1, NEL
         JJ = (I-1)*NM1
         IJ = JJ + 1
         H = 2.0D0/(XBK(I+1)-XBK(I))
         DO IK = 1, 3
            IZ(IK) = 1
         END DO
C        ***************************************************************
C        MAIN LOOP OVER ALL THE SPATIAL ELEMENTS START BY FORMING THE
C        SPACE DERIVS OF U AND UDOT IN DUDX AND UTDX RESPECTIVELY.
C        **************************************************************
         DO K = 1, NPDE
            DO II = 1, NPTL
               DUDX(K,II) = 0.0D0
               UTDX(K,II) = 0.0D0
               DO J = 1, NPTL
                  UTDX(K,II) = UTDX(K,II) + DU(II,J)*UDOT(K,JJ+J)*H
                  DUDX(K,II) = DUDX(K,II) + DU(II,J)*U(K,JJ+J)*H
               END DO
            END DO
         END DO
C
         IF (I.EQ.1) THEN
C           SAVE THE VALUES NEEDED FOR LEFT BOUNDARY CONDITIONS
            DO J = 1, NPDE
               BETA(J,3) = DUDX(J,1)
               BETA(J,4) = UTDX(J,1)
            END DO
         END IF
         IF (I.EQ.NEL) THEN
C           SAVE THE VALUES NEEDED FOR RIGHT BOUNDARY CONDITIONS
            DO J = 1, NPDE
               GAMMA(J,3) = DUDX(J,NPTL)
               GAMMA(J,4) = UTDX(J,NPTL)
            END DO
         END IF
C        ---------------------------------------------------------------
C         EVALUATE THE FUNCTIONS Q AND R IN THIS ELEMENT
C        --------------------------------------------------------------
         CALL SPDEFN(T,X(IJ),NPTL,NPDE,U(1,IJ),DUDX,UDOT(1,IJ),UTDX,Q,R,
     *               IV,V,VDOT,IZ(1))
         IF (M.GT.0) THEN
C           MODIFY Q FUNCTION IF POLAR CO-ORDINATES
            KJ = 1
            IF (X(IJ).LE.TWOU) THEN
               MP1 = 1.0D0 + M
               KJ = 2
               DO K = 1, NPDE
C                 R(K,1) = 0.0D0
                  Q(K,1) = Q(K,1)/(M+1)
               END DO
            END IF
            DO J = KJ, NPTL
               DO K = 1, NPDE
                  Q(K,J) = Q(K,J) - R(K,J)*M/X(JJ+J)
               END DO
            END DO
         END IF
C        ---------------------------------------------------------------
C        SET UP SAVEL AND SAVER FOR BOUNDARY AND INTERFACE CONDITIONS
C        --------------------------------------------------------------
         KJ = MAX(2,I)
         JK = MIN(NEL,I+1) + 1
         SAVEL = 1.0D0/(XBK(KJ)+XBK(I+1)-XBK(KJ-1)-XBK(I))
         SAVER = 1.0D0/(XBK(JK)+XBK(I+1)-XBK(JK-1)-XBK(I))
         IF (I.EQ.1) SFIRST = SAVEL
C        ---------------------------------------------------------------
C         FORM THE RESIDUAL AND THE INTERFACE CONDITIONS
C        --------------------------------------------------------------
         DO J = 1, NPDE
            JK = IJ + NM1
            TEM = R(J,1) + R(J,NPTL)
            RES(J,IJ) = RES(J,IJ) + (Q(J,1)*2.0/H-TEM)*SAVEL
            RES(J,JK) = (Q(J,NPTL)*2.0/H+TEM)*SAVER
         END DO
C        TEST TO SEE IF ILLEGAL SOLUTION VALUES HAVE BEEN FOUND.
         IF (IZ(1).NE.1) THEN
            IRES = IZ(1)
            RETURN
         END IF
      END DO
C**********************************************************************
C    EVALUATE THE FUNCTIONS BETA AND GAMMA AT THE BOUNDARY CONDITIONS
C**********************************************************************
C
      CALL SBNDR(T,BETA(1,1),GAMMA(1,1),U(1,1),BETA(1,3),UDOT(1,1),
     *           BETA(1,4),NPDE,.TRUE.,IV,V,VDOT,IZ(2))
      CALL SBNDR(T,BETA(1,2),GAMMA(1,2),U(1,NPTS),GAMMA(1,3),UDOT(1,
     *           NPTS),GAMMA(1,4),NPDE,.FALSE.,IV,V,VDOT,IZ(3))
C
C                             PROCESS THE BOUNDARY CONDITIONS
      DO J = 1, NPDE
C        L.H.--BOUNDARY CONDITION IS PROCESSED
         RES(J,1) = MP1*(RES(J,1)*BETA(J,1)*2.0D0+GAMMA(J,1)
     *              *4.0D0/CCR(1)*SFIRST)
C        R.H.---BOUNDARY CONDITION IS PROCESSED
         RES(J,NPTS) = RES(J,NPTS)*BETA(J,2)*2.0D0 - GAMMA(J,2)
     *                 *4.0D0/CCR(1)*SAVER
      END DO
      DO IK = 2, 3
         IF (IZ(IK).NE.1) IRES = IZ(IK)
      END DO
      RETURN
C-------END OF DRES----------------------------------------------------
C
      END
C***********************************************************************
C