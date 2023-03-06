      SUBROUTINE CRES(NPDE,NPTS,T,U,RES,UDOT,M,X,OMEGA,DU,XBK,BETA,
     *                GAMMA,DUDX,R,Q,NEL,NPTL,XC,CCR,IRES,RT,QT,UDT,
     *                UTDX,NV,V,VDOT,VDUM)
C**********************************************************************
C       CHEBYSHEV C0 COLLOCATION SPATIAL DISCRETISATION ROUTINE
C       FOR POLYNOMIALS OF DEGREE 2 AND ABOVE.
C**********************************************************************
      USE PDECHEB_COMMON, ONLY: TWOU
      IMPLICIT NONE
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
      DOUBLE PRECISION H, MP1, SAVEL, SAVER, SFIRST
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
C        MAIN LOOP OVER ALL THE SPATIAL ELEMENTS START BY
C        FORMING THE SPACE DERIVS OF U AND UDOT IN DUDX AND UTDX
C        RESPECTIVELY.
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
C        ---------------------------------------------------------------
C        EVALUATE THE FUNCTIONS Q AND R IN THIS ELEMENT
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
C        **************************************************************
C        FORM THE FUNCTIONS BETA AND GAMMA IN THE BOUNDARY CONDITIONS
C        **************************************************************
         IF (I.EQ.1) THEN
C           LEFT HAND BOUNDARY CONDITIONS
            CALL SBNDR(T,BETA(1,1),GAMMA(1,1),U(1,1),DUDX,UDOT(1,1),
     *                 UTDX,NPDE,.TRUE.,IV,V,VDOT,IZ(2))
            IF (IZ(2).NE.1) IRES = IZ(2)
         END IF
         IF (I.EQ.NEL) THEN
C           RIGHT HAND BOUNDARY CONDITIONS
            CALL SBNDR(T,BETA(1,2),GAMMA(1,2),U(1,NPTS),DUDX(1,NPTL),
     *                 UDOT(1,NPTS),UTDX(1,NPTL),NPDE,.FALSE.,IV,V,VDOT,
     *                 IZ(3))
            IF (IZ(3).NE.1) IRES = IZ(3)
         END IF
C        ---------------------------------------------------------------
C        SET UP SAVEL AND SAVER  FOR THE BOUNDARY AND INTERFACE
C        CONDITIONS AND FORM DRDX  BY OVERWRITING DUDX
C        --------------------------------------------------------------
         KJ = MAX(2,I)
         JK = MIN(NEL,I+1) + 1
         SAVEL = 1.0D0/(XBK(KJ)+XBK(I+1)-XBK(KJ-1)-XBK(I))
         SAVER = 1.0D0/(XBK(JK)+XBK(I+1)-XBK(JK-1)-XBK(I))
         IF (I.EQ.1) SFIRST = SAVEL
         DO K = 1, NPDE
            DO II = 1, NPTL
               DUDX(K,II) = 0.0D0
               DO J = 1, NPTL
                  DUDX(K,II) = DUDX(K,II) + DU(II,J)*R(K,J)
               END DO
            END DO
         END DO
C        ---------------------------------------------------------------
C         FORM THE RESIDUAL AND THE INTERFACE CONDITIONS
C        --------------------------------------------------------------
         DO J = 1, NPDE
            DO K = 2, NM1
C              COLLOCATION AT INTERIOR POINT
               RES(J,JJ+K) = Q(J,K) - DUDX(J,K)*H
            END DO
            JK = IJ + NM1
            RES(J,IJ) = RES(J,IJ) + ((Q(J,1)/H-DUDX(J,1)-R(J,1)/CCR(1))
     *                  *2.0)*SAVEL
            RES(J,JK) = ((Q(J,NPTL)/H-DUDX(J,NPTL)+R(J,NPTL)/CCR(1))
     *                  *2.0)*SAVER
         END DO
C        TEST TO SEE IF ILLEGAL SOLUTION VALUES HAVE BEEN FOUND.
         IF (IZ(1).NE.1) THEN
            IRES = IZ(1)
            RETURN
         END IF
      END DO
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
      RETURN
C-------END OF CRES----------------------------------------------------
C
      END