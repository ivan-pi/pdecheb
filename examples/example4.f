c  problem 4
C     ***********************************************************
C     BP PROBLEM - VAPOUR EVAPORATION OVER POOL
C     REGION OF INTEGRATION CONSISTS OF 2 AREAS A VISCOUS SUB-LAYER AND
C     A TURBULENT REGION, (THE DIVISION OCCURS AT X=0.508D-03).
C     THE PDE IS DIFFERENT IN EACH REGION.
C     ***********************************************************
C     C0 COLLOCATION PARAMETERS
        PARAMETER ( IBK   = 8, NEL  = IBK-1 , NPDE = 1, NV = 3,
     1              NPOLY = 03, NPTS = NEL*NPOLY+1,     NXI = NPTS,
     2              NEQ   = NPTS * NPDE + NV,
     3              NWKRES= (NPOLY+1) * (5*NXI + 3*NPOLY+NEL+5+7*NPDE) +
     4                       NPDE * 8 + 3 + NV + NXI,
C     DDASSL TIME INTEGRATION PARAMETERS
     5              MAXORD = 5, LRW = 40 + (MAXORD+4) * NEQ + NEQ**2,
     6              LIW = 20 + NEQ )
C
        INTEGER IWORK(LIW), INFO(15), IBAND, M, ITIME, I, IDID,
     1          IDEV, ITRACE, GRNPTS, IFL, NOUT, KTIME, ITYPE
        DOUBLE PRECISION XBK(IBK), X(NPTS), Y(NEQ), YDOT(NEQ), TINC(11),
     1          WKRES(NWKRES), RWORK(LRW), XI(NXI), T, TOUT, RTOL, ATOL,
     2          U0,VM,DTX1,DTX2,DM1,DM2,K,SCM,PE,MW,RHO,RT,Q3
     3         ,TEND,W Q1,Q2,TEMP, XOUT(100), UOUT(100,1), CPU, XBAR
        REAL  GRX(800), GRY(800), GRZ(800)
        EXTERNAL PDECHB, DGEJAC
C
C       COMMON BLOCKS  TO PASS ACROSS PROBLEM DEPENDENT CONSTANTS.
        COMMON /C0/     PE,MW,RHO,RT,W
        COMMON /PDES/   U0,VM,DTX1,DTX2,DM1,DM2,K
        COMMON /SDEV2/ ITRACE, IDEV
C IBM CALL TO SWITCH OFF UNDERFLOW COMMENTED OUT
C       CALL ERRSET(208, 256, -1, -1, 0)
C CPU TIMER COMMENTED OUT FOR  PORTABILITY
C       CALL TIMER (CPU, 1)
        PE = 0.39005D+4
        MW = 0.92142D+2
        RHO = 0.3767D+1
        RT = 0.8317D+4*0.29815D+3
        U0 = 0.3164D+0
        VM = 0.147D-04
        DTX1 = 0.0D+0
        SCM = 1.7D+0
        DM1 = VM/SCM
        K = 0.41D+0
        DM2 = 0.0D+0
        DTX2 = U0*K
        W = 0.25D0
        GRNPTS = 1
        WRITE(IDEV,9)NPOLY, NEL
 9      FORMAT(' TEST PROBLEM 4'/' ***********'/' POLY OF DEGREE =',I4,
     1         ' NO OF ELEMENTS = ',I4)
         RTOL = 0.1D-4
         ATOL = 0.1D-4
         ITRACE = 0
         IDEV = 4
         WRITE(IDEV,104)RTOL, ATOL, ITRACE, IDEV
104      FORMAT(//' RTOL=',D12.3,' ATOL=',D12.3,' ITRACE AND IDEV=',2I4)
C
         WRITE(4,55)ATOL, RTOL, NPTS
55       FORMAT(//' SOLUTION TO B.P. POOL EVAPORATION PROBLEM USING
     1   DASSL INTEGRATOR WITH FULL MATRIX ROUTINES '/
     2   '   ATOL = ',D11.3,'  RTOL = ',D11.3,'  NPTS = ',I5/)
         NOUT = 20
         XOUT(1) = 0.0D0
         XOUT(2) = 0.127D-3
         XOUT(3) = 0.254D-3
         XOUT(4) = 0.381D-3
         XOUT(5) = 0.508D-3
         XOUT(6) = 0.635D-3
         XOUT(7) = 0.762D-3
         XOUT(8) = 0.889D-3
         XOUT(9) = 0.1D-2
         XOUT(10)= 0.3D-2
         XOUT(11)= 0.5D-2
         XOUT(12)= 0.75D-2
         XOUT(13)= 0.1D-1
         XOUT(14)= 0.3D-1
         XOUT(15)= 0.5D-1
         XOUT(16)= 0.75D-1
         XOUT(17)= 0.1D0
         XOUT(18)= 0.15D0
         XOUT(19)= 0.2D0
         XOUT(20)= 0.22D0
         XBAR = XOUT(5)
         DO 1000 I = 1,NOUT
            TEMP = DLOG10( 1.0D0 + XOUT(I)/XBAR *2.0D0)
            WRITE(IDEV,999)I,XOUT(I),TEMP
999         FORMAT('   I=',I3,' XOUT=',D13.5,'    LOG10=',D13.5)
1000     CONTINUE
C
C      TEMPORARY VALUES OF XI FOR FIRST CALL TO INICHB
         DO 291 I = 1,NPTS
 291       XI(I) =(I-1.0D0) /(NPTS-1.0D0)
C
          XBK(1) = 0.0D0
          XBK(2) = XBAR* 0.5D0
          XBK(3) = XBAR
          XBK(4) = XBAR * 1.5D0
          XBK(5) = XBAR * 2.0D0
          XBK(6) = XBAR*11.0
          XBK(7) = XBAR * 121
          XBK(8) = 1.0D0
          ITIME = 1
C           INITIALISE THE P.D.E. WORKSPACE
        CALL INICHB(NEQ, NPDE, NPTS, X, Y, WKRES, NWKRES, M, T, IBAND,
     1              ITIME, XBK, IBK, NEL, NPOLY, NV, NXI, XI, IDEV)
        DO 292 I = 1,NPTS
C            FINAL VALUES OF XI
 292       XI(I) = X(I)
        CALL INICHB(NEQ, NPDE, NPTS, X, Y, WKRES, NWKRES, M, T, IBAND,
     1              ITIME, XBK, IBK, NEL, NPOLY, NV, NXI, XI, IDEV)
        IF(ITIME .EQ. -1)THEN
           WRITE(IDEV, 15)
 15        FORMAT(' INITCC ROUTINE RETURNED ITIME = -1 - RUN HALTED ')
           GOTO 900
         ELSE
           WRITE(IDEV,16)(Y(I), I = 1,NPTS)
  16       FORMAT(' INITIAL VALUES ARE =',5D11.3)
        END IF
C         SETUP DASSL PARAMETERS
       DO 20 I = 1,11
 20      INFO(I) = 0
C      INFO(11)= 1
C         BANDED MATRIX OPTION WHEN INFO(6) = 1
       IF(INFO(6) .EQ. 1)THEN
          IWORK(1) = IBAND
          IWORK(2) = IBAND
       END IF
      T    = 0.0D0
      TINC(1) = 0.0001D0
      TINC(2) = 0.0010
      TINC(3) = 0.01D0
      TINC(4) = 0.050D0
      TINC(5) = 0.1D0
      TINC(6) = 0.15D0
      TINC(7) = 0.25D0
      TINC(8) = 0.50D0
      TINC(9) = 0.65D0
      TINC(10)= 0.80D0
      TINC(11)= 1.00D0
      TEND = 1.0D0
      KTIME = 1
      ITYPE = 1
      CALL INTERC(XOUT,UOUT,NOUT,Y,NEQ,NPDE,IFL,ITYPE,WKRES,NWKRES)
      WRITE (IDEV,82) T, (UOUT(I,1),I=1,NOUT,3)
      GRNPTS = 0
      DO 800 I = 1,NOUT
         GRNPTS = GRNPTS + 1
         GRX(GRNPTS) = T
         GRZ(GRNPTS) = UOUT(I,1)/UOUT(1,1)
         GRY(GRNPTS) = DLOG10( 1.D0+XOUT(I)/XBAR * 2.0D0)
         IF(ITRACE .GE.0)WRITE(IDEV,899)GRY(GRNPTS),GRZ(GRNPTS)
 800  CONTINUE
      WRITE(IDEV,81)(XOUT(I), I = 1,NOUT,2)
 81   FORMAT (/'  T/X', 4X,9D11.3)
C TIME LOOP:
100   TOUT = TINC(KTIME)
       CALL DDASSL( PDECHB, NEQ, T, Y, YDOT, TOUT, INFO, RTOL, ATOL,
     1              IDID, RWORK, LRW, IWORK, LIW, WKRES, NWKRES, DGEJAC)
C          DASSL FAILED TO FINISH INTEGRATION.
       WRITE(IDEV,40)T,IDID
 40    FORMAT(' AT TIME T = ',D11.3,' DASSL RETURNED IDID =',I3)
       IF( IDID .LT. 0 )GOTO 900
C        DASSL INTEGRATED TO T = TOUT
C        CALL TO POST PROCESSING HERE E.G. SPACE INTERPOLATION.
       CALL INTERC(XOUT,UOUT,NOUT,Y,NEQ,NPDE,IFL,ITYPE,WKRES,NWKRES)
 82    FORMAT(1X,F3.1,' U   ',9D11.3/)
       WRITE (IDEV,82) TOUT, (UOUT(I,1),I=1,NOUT,3)
       WRITE (6,82) TOUT, (UOUT(I,1),I=1,NOUT,3)
C
C COMPARE RATE OF EVAPORATION Q1 AT SURFACE OF POOL WITH QUANTITY OF
C VAPOUR Q2 WHICH PASSES ABOVE END OF POOL
       Q1 = Y(NEQ-2)
       Q2 = Y(NEQ-1)
       Q3 = Y(NEQ)
       WRITE(IDEV,83) Q1,Q2,Q3
 83    FORMAT(' Q1 , Q2 AND Q3 ARE ',3D13.5)
C
C        PUT INTERPOLATED RESULTS IN ARRAY.
C
      I =(KTIME/2) * 2
      IF (I .EQ. KTIME)GOTO 91
      DO 90 I = 1,NOUT
         GRNPTS = GRNPTS + 1
         GRX(GRNPTS) = TOUT
         GRZ(GRNPTS) = UOUT(I,1)/UOUT(1,1)
         GRY(GRNPTS) = DLOG10( 1.D0+XOUT(I)/XBAR * 2.0D0)
         IF(ITRACE .GE.0)WRITE(IDEV,899)GRY(GRNPTS),GRZ(GRNPTS)
 899     FORMAT(' X AND Y VALUES ARE ',2E12.4)
 90   CONTINUE
 91   KTIME = KTIME + 1
C
C     CHECK IF INTEGRATION WAS SUCCESSFUL AND WHETHER FURTHER TIME
C     STEPS NEEDED
      IF(TOUT.LT.TEND.AND.(IDID.EQ.2 .OR. IDID .EQ. 3)) GO TO 100
      WRITE(IDEV,2112)Q1,Q2,DABS(Q3)
 2112 FORMAT(' RATE OF EVAPORATION AT SURFACE OF POOL Q1 = ',D14.7,/
     -  ' QUANTITY OF VAPOUR ABOVE END OF POOL   Q2 = ',D14.7,/
     -  ' ABSOLUTE DIFFERENCE Q3 = ',D11.4,/
     -  '********************************************************',/)
80    CONTINUE
C     CALL TIMER(CPU,2)
900    WRITE(IDEV,110)IWORK(11),IWORK(12),IWORK(13), CPU
110    FORMAT(' NSTEPS =',I5,' NRESID =',I5,' JAC = ',I4,' CPU=',D11.3)
       STOP
       END
C EXAMPLE PROBLEM FOUR
C *********************
C     THIS PROBLEM IS DEFINED BY
C
C            C X  U  = ( C   U  )            ,  X IN (0 , XBAR)
C             1    T      2   X  X
C  AND
C
C  (C  LOG(X) +C ) U = ( C   X U  )          ,  X IN (XBAR , 1)
C    3          4   T      5    X  X
C
C  WHERE                   -6
C   C = 6810.0  C = 8.65 10    C  =0.7717   C = 9.313  C = 0.1297
C    1           2              3            4          5
C
C     THE LEFT BOUNDARY CONDITION AT X =-1 (LEFT = .TRUE. ) IS GIVEN BY
C         U(0,T) = 0.038475
C
C     THE RIGHT BOUNDARY CONDITION IS  (LEFT = .FALSE.)
C         U (1,T) = 0
C          X
C
C      THE INITIAL CONDITION IS GIVEN BY
C         U(X,0)  = 0
C
C      THE ALGEBRAIC VARIABLES Q (T)  , Q (T)  AND Q (T)  ARE DEFINED BY
C                               1        2          3
C
C      .                  -7
C      Q     =   -7.983 10    U (0 , T)
C       1                      X
C
C                         -2   1
C      Q     =   9.4175 10    I  P(X) U(X,T) DX
C       2                    0
C
C            WHERE     P(X) = C  X               FOR  X IN (0, XBAR)
C                              1
C
C                      P(X) = C  LOG(X) + C      FOR X IN (XBAR, 1)
C                              3           4
C            AND THE VALUES OF THE CONSTANTS ARE GIVEN ABOVE.
C
C       Q (T)    =   Q (T)  - Q (T)
C        3            2        1
C
C**********************************************************************
       SUBROUTINE UVINIT( NPDE, NPTS, X, U, NV,V)
C      ROUTINE FOR P.D.E. INITIAL VALUES.
       INTEGER NPDE, NPTS, NV, I
       DOUBLE PRECISION X(NPTS), U(NPDE,NPTS),PE,MW,RHO,RT,W,V(3)
       COMMON/C0/PE,MW,RHO,RT,W
       DO 10 I= 2,NPTS
 10       U(1,I) = 0.0D+0
       U(1,1) = (PE*MW)/(RHO*RT)
       V(1) = 0.0D0
       V(2) = 0.0D0
       V(3) = 0.0D0
       RETURN
       END
C
       SUBROUTINE SPDEFN( T, X, NPTL, NPDE, U, DUDX, UDOT, UTDX, Q, R,
     1                    NV, V, VDOT, IRES)
C**********************************************************************
C      ROUTINE TO DESCRIBE THE BODY OF THE P.D.E.
C      THE P.D.E. IS WRITEN AS       -M   M
C         Q(X,T,U, U  , U  , U  ) = X   (X  R(X,T,U,U , U , U  ))
C                   X    T    TX                     X   T   TX  X
C      THE FUNCTIONS Q AND R MUST BE DEFINED IN THIS ROUTINE.
C**********************************************************************
       INTEGER NPDE, NPTL, I, NV, IRES
       DOUBLE PRECISION T, X(NPTL), U(NPDE,NPTL), DUDX(NPDE,NPTL),
     1    DM2, UDOT(NPDE,NPTL), Q(NPDE,NPTL), R(NPDE,NPTL), V(3),
     2    K, UTDX(NPDE,NPTL), U0, VM, DTX1, DTX2, DM1, VDOT(3)
      COMMON /PDES/  U0,VM,DTX1,DTX2,DM1,DM2,K
      DO 100 I = 1,NPTL
        IF(X(1) .LT. 0.506D-3 .AND. X(NPTL) .LT. 0.600D-3)THEN
C         ELEMENT TO LEFT OF THE INTERFACE AT 0.508D-3
          Q(1,I) = (X(I)*U0**2)/VM * UDOT(1,I)
          R(1,I) = (DTX1 + DM1)*DUDX(1,I)
        ELSE
          Q(1,I) = ((U0/K)*DLOG(U0*X(I)/VM) + 5.1*U0) * UDOT(1,I)
          R(1,I) = ((DTX2*X(I)) + DM2)*DUDX(1,I)
        ENDIF
 100  CONTINUE
      RETURN
      END
C
       SUBROUTINE SBNDR( T, BETA, GAMMA, U, UX, UDOT, UTDX, NPDE, LEFT,
     1                   NV, V, VDOT, IRES)
C      BOUNDARY CONDITIONS ROUTINE
       INTEGER NPDE, NV, IRES
       DOUBLE PRECISION T, BETA(NPDE), GAMMA(NPDE), U(NPDE), PE,MW,RHO,
     1         UX(NPDE), V(3), VDOT(3), UDOT(NPDE), UTDX(NPDE), RT, W
      LOGICAL LEFT
      COMMON/C0/PE,MW,RHO,RT,W
      IF(LEFT) THEN
          GAMMA(1) = U(1)- (PE*MW)/(RHO*RT)
          BETA(1) = 0.0D+0
      ELSE
          GAMMA(1) = 0.0D0
          BETA(1) = 1.0D0
      END IF
      RETURN
      END
      SUBROUTINE SODEFN(T, NV, V, VDOT, NPDE, NXI, X, Y, UXI, RI,
     1                  UTI, UTXI, VRES, IRES)
C ROUTINE FOR AUXILIARY O.D.E.S (IF ANY) IN MASTER EQN. FORM (4.3)
      INTEGER NPDE, NXI, NV, IRES, NPTL, L, J, I
      DOUBLE PRECISION T, X(NXI), Y(NXI), UXI(NPDE,NXI),
     1        RI(NPDE,NXI), UTI(NPDE,NXI), UTXI(NPDE,NXI), VRES(NV),
     2        V(3), VDOT(3), PE,MW,RHO,RT,W,U0,VM,DTX1,DTX2,DM2,K
     3       ,DM1, Q2, H, CCRULE
      COMMON  /C0/     PE,MW,RHO,RT,W
      COMMON  /PDES/   U0,VM,DTX1,DTX2,DM1,DM2,K
      COMMON  /SCHSZ5/ NPTL
      COMMON  /SCHSZ6/ CCRULE(50)
C
      VRES(1) = VDOT(1) + W*RHO*DM1*UXI(1,1)
      Q2 = 0.0D0
      DO 3 I = 1,2
        J = (NPTL-1) * (I-1) + 1
        L = (NPTL-1) *  I    + 1
        H = ( X(L) - X(J)) * 0.5D0
      DO 3 II = 1,NPTL
         IK = J + II - 1
C        CLENSHAW - CURTIS QUADRATURE UP TO INTERFACE POINT.
         Q2 = Q2 + (W*RHO*U0**2)/VM * X(IK) * Y(IK) * CCRULE(II) * H
 3    CONTINUE
C
C        CLENSHAW - CURTIS QUADRATURE BEYOND THE INTERFACE POINT.
      DO 5 I = 3,7
        J = (NPTL-1) * (I-1) + 1
        L = (NPTL-1) *  I    + 1
        H = ( X(L) - X(J)) * 0.5D0
        DO 5 II = 1, NPTL
         IK = J + II - 1
         Q2=Q2 + H* ((U0/K)*DLOG(U0*X(IK)/VM)+5.1*U0) * Y(IK)*CCRULE(II)
     1            * W * RHO
  5   CONTINUE
      VRES(2) = V(2) - Q2
      VRES(3) = V(3) - (V(2)-V(1))
      RETURN
      END