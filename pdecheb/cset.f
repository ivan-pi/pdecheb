      SUBROUTINE CSET(NPDE,NPTS,U,X,OMEGA,DU,XBK,NEL,NPTL,XC,CCR,XBH,
     *                IBK,DUTEM,V,NV)
      use pdecheb_common, only: ccrule
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER         IBK, NEL, NPDE, NPTL, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION CCR(NPTL), DU(NPTL,NPTL), DUTEM(NPTL,NPTL),
     *                OMEGA(NPTL,NPTL), U(NPDE,NPTS), V(*), X(NPTS),
     *                XBH(IBK), XBK(IBK), XC(NPTL)
C     .. Local Scalars ..
      DOUBLE PRECISION H1, H2, PI, SINT, SUM, TEMP
      INTEGER         I, IJ, ITEM, J, K, NM1, NT, NTP1
      PARAMETER(PI=3.1415926535897930D0)
C     .. External Subroutines ..
      EXTERNAL        UVINIT
C     .. Intrinsic Functions ..
      INTRINSIC       DBLE, COS, SIN
C     .. Executable Statements ..
C
C  FORM  CONSTANTS FOR WKSPACE INITIALISATION
C
      NM1 = NPTL - 1
C
C  FORMATION OF GRID AND INITIAL VALUES OF U
C
      DO I = 1, NEL
         H1 = XBH(I+1) - XBH(I)
         H2 = XBH(I+1) + XBH(I)
         XBK(I) = XBH(I)
         DO J = 1, NPTL
            IJ = (I-1)*NM1 + J
            IF (I.EQ.1) XC(J) = COS(PI*DBLE(J-NPTL)/NM1)
            X(IJ) = (XC(J)*H1+H2)*0.5D0
            IF (J.EQ.1) X(IJ) = XBH(I)
            IF (J.EQ.NPTL) X(IJ) = XBH(I+1)
         END DO
      END DO
      XBK(IBK) = XBH(IBK)
      XC(1) = -1.0D0
      XC(NPTL) = 1.0D0
C
C  FORM THE MATRIX OMEGA
C
      DO J = 1, NPTL
         DO I = 1, NPTL
            OMEGA(I,J) = 2.D0*COS(PI*(I-1)*(NPTL-J)/NM1)/NM1
         END DO
      END DO
C
C   MODIFY EDGES OF OMEGA AND FORM EDGES OF INTERMEDIATE DU MATRIX
C
      ITEM = 1
      DO I = 1, NPTL
         OMEGA(I,1) = OMEGA(I,1)*0.5D0
         OMEGA(1,I) = OMEGA(1,I)*0.5D0
         OMEGA(NPTL,I) = OMEGA(NPTL,I)*0.5D0
         OMEGA(I,NPTL) = OMEGA(I,NPTL)*0.5D0
         DUTEM(I,1) = 0.0D0
         DUTEM(1,I) = -DBLE((I-1)**2*ITEM)
         DUTEM(NPTL,I) = DBLE((I-1)**2)
         ITEM = -ITEM
      END DO
C
C FINISH FORMING REST OF INTERMEDIATE DU MATRIX THAT IS HELD IN DUTEM.
C
      IF (NPTL.GT.2) THEN
         DO I = 2, NM1
            TEMP = PI*(I-NPTL)/NM1
            SINT = SIN(TEMP)
            DO J = 2, NM1
               DUTEM(I,J) = SIN(TEMP*(J-1))/SINT*(J-1)
            END DO
            DUTEM(I,NPTL) = 0.0D0
         END DO
      END IF
C
C  FORM FULL DU BY MATRIX MULTIPLICATION
C
      DU = MATMUL(DUTEM,OMEGA)
C
C        CALCULATE THE COEFFS OF THE CLENSHAW CURTIS RULE
C
      NT = NM1/2
      IF ((2*NT).NE.NM1) NT = (NM1-1)/2
      NTP1 = NT + 1
      SUM = 0.0D0
      DO I = 1, NPTL
         TEMP = 0.5D0
         CCR(I) = 0.0D0
         DO K = 1, NTP1
            IF (K.EQ.NTP1 .AND. ((2*NT).EQ.NM1)) TEMP = 0.5D0
            CCR(I) = CCR(I) + COS(2.0D0*(I-1)*(K-1)*PI/NM1)
     *               *TEMP/(4.0D0*(K-1)**2-1.0D0)
            TEMP = 1.0D0
         END DO
         IF (I.EQ.1 .OR. I.EQ.NPTL) TEMP = 0.5D0
         CCR(I) = CCR(I)*(-4.0D0)*TEMP/NM1
         SUM = SUM + CCR(I)
      END DO
      DO I = 1, NPTL
         CCRULE(I) = CCR(I)
      END DO
      DO I = 2, NM1
         CCR(I) = CCR(I)/CCR(1)
      END DO
C  FIND THE INITIAL VALUES OF THE O.D.E. AND P.D.E. COMPONENTS.
      CALL UVINIT(NPDE,NPTS,X,U,NV,V)
      RETURN
C
C-----------END  OF  CSET ROUTINE---------------------------------------
C
      END