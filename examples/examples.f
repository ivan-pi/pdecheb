C      ALGORITHM 690, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 17, NO. 2, PP. 178-206.  JUNE, 1991.
c
c  file containg example programs and PDECHEB software.
c
c  This file contains
c       example problem 1
c       example problem 2
c       example problem 3
c       example problem 4
c       example problem 5
c       dassl time integration routine
c       pdecheb spatial discretisation routine
c       interface to AND the linpack full and banded routines
c
c
C     C0 COLLOCATION PARAMETERS
        PARAMETER ( IBK   = 21, NEL  = IBK-1 , NPDE = 1, NV = 1,
     1              NPOLY =  2,  NPTS = NEL*NPOLY+1,     NXI = 1,
     2              NEQ   = NPTS * NPDE + NV,
     3              NWKRES= (NPOLY+1) * (5*NXI + 3*NPOLY+NEL+5+7*NPDE) +
     4                       NPDE * 8 + 3 + NV + NXI,
C     DDASSL TIME INTEGRATION PARAMETERS
     5              MAXORD = 5, LRW = 40 + (MAXORD+4) * NEQ + NEQ**2,
     6              LIW = 20 + NEQ )
C
        INTEGER IWORK(LIW), INFO(15), IBAND, M, ITIME, I, IDID, IRESWK,
     1          IDEV, ITRACE
        DOUBLE PRECISION XBK(IBK), X(NPTS), Y(NEQ), YDOT(NEQ),
     1          WKRES(NWKRES), RWORK(LRW), XI(1), T, TOUT, RTOL, ATOL,
     2          ENORM, GERR, VERROR, CTIME, TOL
        EXTERNAL PDECHB, DGEJAC
        COMMON /SDEV2/ ITRACE, IDEV
        COMMON /PROB1/ TOL
        TOL  = 0.1D-5/50.D0
C N.B. CPU TIMER COMMENTED OUT FOR PORTABILITY
C       CALL TIMER( CTIME, 1)
        M    = 0
        T    = TOL
        IDEV = 4
        ITRACE = 1
        WRITE(IDEV,9)NPOLY, NEL
 9      FORMAT(' TEST PROBLEM 1'/' ***********'/' POLY OF DEGREE =',I4,
     1         ' NO OF ELEMENTS = ',I4)
        XI(1)  = 1.0D0
        DO 10 I = 1,IBK
 10       XBK(I) = (I-1.0D0)/(IBK-1.0D0)
C           INITIALISE THE P.D.E. WORKSPACE
        ITIME  = 1
        CALL INICHB(NEQ, NPDE, NPTS, X, Y, WKRES, NWKRES, M, T, IBAND,
     1              ITIME, XBK, IBK, NEL, NPOLY, NV, NXI, XI, IDEV)
        IF(ITIME .EQ. -1)THEN
           WRITE(IDEV, 15)
 15        FORMAT(' INITCC ROUTINE RETURNED ITIME = -1 - RUN HALTED ')
           GOTO 100
        END IF
C         SETUP DASSL PARAMETERS
       RTOL = TOL
       ATOL = TOL
       DO 20 I = 1,11
 20      INFO(I) = 0
C
C         BANDED MATRIX OPTION WHEN INFO(6) = 1
       IF(INFO(6) .EQ. 1)THEN
          IWORK(1) = IBAND
          IWORK(2) = IBAND
       END IF
 30    TOUT = T * 10.0D0
       IF(TOUT .GE. 2.D0)TOUT =2.0D0
       CALL DDASSL( PDECHB, NEQ, T, Y, YDOT, TOUT, INFO, RTOL, ATOL,
     1              IDID, RWORK, LRW, IWORK, LIW, WKRES, IRESWK, DGEJAC)
       IF( IDID .LT. 0 )THEN
C          DASSL FAILED TO FINISH INTEGRATION.
           WRITE(IDEV,40)T,IDID
 40        FORMAT(' AT TIME T = ',D11.3,' DASSL RETURNED IDID =',I3)
           GOTO 100
       ELSE
C        DASSL INTEGRATED TO T = TOUT
C        CALL TO POST PROCESSING HERE E.G. SPACE INTERPOLATION.
         ITRACE = 1
         CALL ERROR( Y, NPDE, NPTS, X, M, ENORM, GERR, T, RTOL, ATOL,
     1               ITRACE, WKRES, NWKRES)
         ITRACE = 0
         VERROR  = Y(NEQ) - T
         WRITE(IDEV,50)Y(NEQ),VERROR
 50      FORMAT(' MOVING BOUNDARY IS AT ',D12.4,' WITH ERROR=',D12.4)
         IF(TOUT .LT. 1.99D0 ) GOTO 30
       END IF
100    CONTINUE
C      CALL TIMER(CTIME, 2)
       WRITE(IDEV,110)IWORK(11),IWORK(12),IWORK(13), CTIME
110    FORMAT(' NSTEPS =',I5,' NRESID =',I5,' JAC = ',I4,' CPU=',D11.3)
       STOP
       END
C**********************************************************************
C EXAMPLE  PROBLEM 1
C SOLUTION OF MOVING BOUNDARY  PROBLEM BY CO-ORDINATE TRANSFORMATION.
C********************************************************************
C  THIS PROBLEM IS THE ONE PHASE STEFAN PROBLEM (HOFFMAN (1977) ) SEE
C  FURZELAND R.M. A COMPARATIVE STUDY OF NUMERICAL METHODS FOR MOVING
C  BOUNDARY PROBLEMS. J.I.M.A. (1977) ,26, PP 411 - 429.
C  THE PROBLEM HAS  MELTING DUE TO HEAT INPUT AT THE FIXED
C  BOUNDARY . THE P.D.E. IS DEFINED BY THE EQUATIONS
C         U  =  U        0 < Y < S(T) , 0.1 < T < 1
C          T     YY
C            U  = - EXP(T) , Y = 0
C             Y            .
C            U  =  0  AND  S(T) = - U   ON THE MOVING BOUNDARY Y = S(T).
C                                    Y
C  AND THE INITIAL SOLUTION VALUES AT T = 0.1 ARE GIVEN BY THE ANALYTIC
C  SOLUTION
C            U = EXP(T-Y) - 1 , S(T) = T.
C  THE PROBLEM IS REWRITTEN BY USING THE CO-ORDINATE TRANSFORMATION
C  GIVEN BY  X(T)  =  Y / S(T)  . THE EQUATIONS THEN READ
C                      .
C     S * S * U  - S * S  * X * U   =  U     , X IN (0,1).
C              T                 X      XX
C  WITH THE NEUMANN TYPE BOUNDARY CONDITIONS
C                                               .
C     U  = - EXP(T)  AT X=0  AND  U  = - S(T) * S(T) AT X = 1
C      X                           X
C  AND THE O.D.E. COUPLING POINT EQUATION AT X = 1 WHICH IMPLICITLY
C  DEFINES S(T) IS GIVEN  BY
C     U(1,T) = 0
C  THE EXACT SOLUTION IS NOW DEFINED BY
C     U(X,T) = EXP((T - X*S(T))  , S(T) = T
C
C WE SHALL NOW DETAIL THE ROUTINES NEEDED TO DESCRIBE THIS PROBLEM.
C        PROBLEM DESCRIPTION ROUTINES
C       ******************************
C EXACT SOLUTION
       SUBROUTINE EXACT( TIME, NPDE, NPTS, X, U)
C      ROUTINE FOR P.D.E. EXACT VALUES  (IF KNOWN)
       INTEGER NPDE, NPTS
       DOUBLE PRECISION X(NPTS), U(NPDE,NPTS), TIME
       DO 10 I = 1,NPTS
 10       U(1,I) = DEXP( TIME * (1 - X(I))) - 1.0D0
       RETURN
       END
       SUBROUTINE UVINIT( NPDE, NPTS, X, U, NV, V)
C      ROUTINE FOR O.D.E. AND P.D.E. INITIAL VALUES.
       INTEGER NPDE, NPTS, NV
       DOUBLE PRECISION X(NPTS), U(NPDE,NPTS), TIME, V(NV)
       COMMON /PROB1/ TOL
       TIME= TOL
       V(1)= TOL
       CALL EXACT(TIME,NPDE,NPTS,X,U)
       RETURN
       END
C
       SUBROUTINE SPDEFN(T, X, NPTL, NPDE, U, DUDX, UDOT, UTDX, Q, R,
     1                   NV, V, VDOT, IRES)
C      PROBLEM INTERFACE  FOR THE MOVING BOUNDARY PROBLEM.
       INTEGER NPTL, NPDE, NV, I, IRES
       DOUBLE PRECISION X(NPTL), U(NPDE,NPTL), DUDX(NPDE,NPTL), T,
     1         V(1), VDOT(1), Q(NPDE,NPTL) ,R(NPDE,NPTL),
     2         UDOT(NPDE,NPTL), UTDX(NPDE,NPTL)
       DO 10 I = 1,NPTL
          R(1,I) = DUDX(1,I)
          Q(1,I) = V(1)*V(1)*UDOT(1,I) -X(I)*VDOT(1)*DUDX(1,I) * V(1)
 10    CONTINUE
       RETURN
       END
       SUBROUTINE SBNDR( T, BETA, GAMMA, U, UX, UDOT, UTDX, NPDE,
     1                   LEFT, NV, V, VDOT, IRES)
C  THIS SUBROUTINE PROVIDES THE LEFT AND RIGHT BOUNDARY VALUES
C  FOR THE MOVING BOUNDARY PROBLEM IN THE FORM.
C           BETA(I) * DU/DX(I) = GAMMA(I)
C  WHERE I = 1,NPDE AND GAMMA IS A FUNCTION OF U,X AND T
C
       INTEGER NPDE, NV, IRES
       LOGICAL LEFT
       DOUBLE PRECISION BETA(NPDE), GAMMA(NPDE), U(NPDE), UX(NPDE)
     -                  ,T, V(1), VDOT(1), UDOT(NPDE), UTDX(NPDE)
       BETA(1) = 1.0D0
       IF(LEFT)THEN
          GAMMA(1) = -V(1)*DEXP(T)
       ELSE
          GAMMA(1) = -V(1)*VDOT(1)
       END IF
       RETURN
       END
C
       SUBROUTINE SODEFN(T, NV, V, VDOT, NPDE, NXI, XI, UI, UXI, RI,
     1                    UTI, UTXI, VRES, IRES)
C      ROUTINE TO PROVIDE RESIDUAL OF COUPLED ODE SYSTEM FOR THE
C      MOVING BOUNDARY PROBLEM.
C      NOTE HOW IRES CAN BE RESET TO COPE WIH ILLEGAL VALUES OF THE
C           MOVING BOUNDARY POSITION V(1).
       INTEGER NPDE, NXI, NV, IRES
       DOUBLE PRECISION T, XI(NXI), UI(NPDE,NXI), UXI(NPDE,NXI),
     1         RI(NPDE,NXI), UTI(NPDE,NXI), UTXI(NPDE,NXI), VRES(NV),
     2         V(NV), VDOT(NV)
       VRES(1) = UI(1,1)
       IF(V(1) .LT. 0.0D0)IRES = -1
       RETURN
       END
C     C0 COLLOCATION PARAMETERS
        PARAMETER ( IBK   =  2, NEL  = IBK-1 , NPDE = 1, NV = 0,
     1              NPOLY = 10, NPTS = NEL*NPOLY+1,     NXI = 0,
     2              NEQ   = NPTS * NPDE + NV,
C    C              NWKRES= 2*(NPOLY+1)*(NPOLY+NEL+2) + 2 + NV +
     3              NWKRES= (NPOLY+1) * (5*NXI + 3*NPOLY+NEL+5+7*NPDE) +
     4                       NPDE * 8 + 3 + NV + NXI,
C    C                       NPDE * (7 * (NPOLY+1+NXI) + 8),
C     DDASSL TIME INTEGRATION PARAMETERS
     5              MAXORD = 5, LRW = 40 + (MAXORD+4) * NEQ + NEQ**2,
     6              LIW = 20 + NEQ )
C
        INTEGER IWORK(LIW), INFO(15), IBAND, M, ITIME, I, IDID, IRESWK,
     1          IDEV, ITRACE
        DOUBLE PRECISION XBK(IBK), X(NPTS), Y(NEQ), YDOT(NEQ),
     1          WKRES(NWKRES), RWORK(LRW), XI(1), T, TOUT, RTOL, ATOL,
     2          ENORM, GERR, CTIME
        EXTERNAL PDECHB, DGEJAC
        COMMON /SDEV2/ ITRACE, IDEV
C  CPU TIMER COMMENTED OUT FOR PORTABILITY
C       CALL TIMER(CTIME ,1)
        M    = 2
        T    = 0.0D0
        IDEV = 4
        ITRACE = 1
        WRITE(IDEV,9)NPOLY, NEL
 9      FORMAT(' TEST PROBLEM 1'/' ***********'/' POLY OF DEGREE =',I4,
     1         ' NO OF ELEMENTS = ',I4)
        DO 10 I = 1,IBK
 10       XBK(I) =          (I-1.0D0) / (IBK - 1.0D0)
C           INITIALISE THE P.D.E. WORKSPACE
        ITIME = 1
        CALL INICHB(NEQ, NPDE, NPTS, X, Y, WKRES, NWKRES, M, T, IBAND,
     1              ITIME, XBK, IBK, NEL, NPOLY, NV, NXI, XI, IDEV)
        IF(ITIME .EQ. -1)THEN
           WRITE(IDEV, 15)
 15        FORMAT(' INITCC ROUTINE RETURNED ITIME = -1 - RUN HALTED ')
           GOTO 100
        END IF
C         SETUP DASSL PARAMETERS
       RTOL = 1.0D-8
       ATOL = 1.0D-8
       DO 20 I = 1,11
 20      INFO(I) = 0
C      INFO(11)= 1
C         BANDED MATRIX OPTION WHEN INFO(6) = 1
       IF(INFO(6) .EQ. 1)THEN
          IWORK(1) = IBAND
          IWORK(2) = IBAND
       END IF
 30    TOUT = T + 0.1D0
       CALL DDASSL( PDECHB, NEQ, T, Y, YDOT, TOUT, INFO, RTOL, ATOL,
     1              IDID, RWORK, LRW, IWORK, LIW, WKRES, IRESWK, DGEJAC)
       IF( IDID .LT. 0 )THEN
C          DASSL FAILED TO FINISH INTEGRATION.
           WRITE(IDEV,40)T,IDID
 40        FORMAT(' AT TIME T = ',D11.3,' DASSL RETURNED IDID =',I3)
           GOTO 100
       ELSE
C        DASSL INTEGRATED TO T = TOUT
C        CALL TO POST PROCESSING HERE E.G. SPACE INTERPOLATION.
         CALL ERROR( Y, NPDE, NPTS, X, M, ENORM, GERR, T, RTOL, ATOL,
     1               ITRACE, WKRES, NWKRES)
         IF(TOUT .LT. 0.99D0 ) GOTO 30
       END IF
100    CONTINUE
C      CALL TIMER(CTIME, 2)
       WRITE(IDEV,110)IWORK(11),IWORK(12),IWORK(13), CTIME
110    FORMAT(' NSTEPS =',I5,' NRESID =',I5,' JAC = ',I4,' CPU=',D11.3)
       STOP
       END
C EXAMPLE PROBLEM TWO
C ********************
C     THIS PROBLEM IS DEFINED BY
C             -2    2               2
C     U U  = X   ( X  U U  )   + 5 U  + 4 X U U     ,  X IN (0,1)
C        T                X  X                 X
C
C     THE LEFT BOUNDARY CONDITION AT X = 0 (LEFT = .TRUE. ) IS GIVEN BY
C        U (0,T)  = 0.0
C         X
C     THE RIGHT BOUNDARY CONDITION IS  (LEFT = .FALSE.)
C         U( 1,T) = EXP ( -T )
C
C      THE INITIAL CONDITION IS GIVEN BY THE EXACT SOLUTION ;
C        U( X, T )  = EXP ( 1 - X*X - T )  , X IN ( 0,1)
C                            2
C**********************************************************************
       SUBROUTINE UVINIT( NPDE, NPTS, X, U, NV,V)
C      ROUTINE FOR P.D.E. INITIAL VALUES.
       INTEGER NPDE, NPTS, NV
       DOUBLE PRECISION X(NPTS), U(NPDE,NPTS), V(1), T
         T = 0.0D0
C        V(1) IS A DUMMY VARIABLE AS THERE ARE NO COUPLED O.D.E.S
         CALL EXACT( T, NPDE, NPTS, X, U )
       RETURN
       END
C
       SUBROUTINE SPDEFN( T, X, NPTL, NPDE, U, DUDX, UDOT, UTDX, Q, R,
     1                    NV, V, VDOT, IRES)
C      ROUTINE TO DESCRIBE THE BODY OF THE P.D.E.
C      THE P.D.E. IS WRITEN AS       -M   M
C         Q(X,T,U, U  , U  , U  ) = X   (X  R(X,T,U,U , U , U  ))
C                   X    T    TX                     X   T   TX  X
C      THE FUNCTIONS Q AND R MUST BE DEFINED IN THIS ROUTINE.
C      DEFINITIONS FOR THE MODEL PROBLEM ARE GIVEN
C      NOTE NV = 0 : THERE IS NO O.D.E PART.
       INTEGER NPDE, NPTL, I, J, NV, IRES
       DOUBLE PRECISION T, X(NPTL), U(NPDE,NPTL), DUDX(NPDE,NPTL),
     1         UDOT(NPDE,NPTL), Q(NPDE,NPTL), R(NPDE,NPTL), V, VDOT,
     2         UTDX(NPDE,NPTL)
       DO 10 I = 1,NPTL
          R(1,I) = U(1,I) * DUDX(1,I)
          Q(1,I) = U(1,I) * UDOT(1,I) - 5.0D0 * U(1,I)**2
     1                                - 4.0D0 * U(1,I)*DUDX(1,I)*X(I)
 10    CONTINUE
       RETURN
       END
C
       SUBROUTINE SBNDR( T, BETA, GAMMA, U, UX, UDOT, UTDX, NPDE, LEFT,
     1                   NV, V, VDOT, IRES)
C      BOUNDARY CONDITIONS ROUTINE
       INTEGER NPDE, NV, IRES
       DOUBLE PRECISION T, BETA(NPDE), GAMMA(NPDE), U(NPDE), C2,
     1                  UX(NPDE), V, VDOT, UDOT(NPDE), UTDX(NPDE)
       LOGICAL LEFT
       IF(LEFT) THEN
          BETA (1) = 1.0D0
          GAMMA(1) = 0.0D0
       ELSE
C         BETA (1) = 0.0D0
C         GAMMA(1) = U(1) - DEXP( -T )
          BETA (1) = 1.0D0
          GAMMA(1) = - 2.D0 *U(1)**2
       END IF
       RETURN
       END
C
C      DUMMY O.D.E. ROUTINE AS NV IS ZERO
       SUBROUTINE SODEFN
       RETURN
       END
C EXACT SOLUTION
       SUBROUTINE EXACT( TIME, NPDE, NPTS, X, U)
C      ROUTINE FOR P.D.E. EXACT VALUES  (IF KNOWN)
       INTEGER NPDE, NPTS, I
       DOUBLE PRECISION X(NPTS), U(NPDE,NPTS), TIME
       DO 10 I = 1,NPTS
 10       U(1,I) = DEXP( 1.0D0 - X(I)**2 - TIME)
       RETURN
       END
c  problem 3
C     C0 COLLOCATION PARAMETERS
        PARAMETER ( IBK   = 3, NEL  = IBK-1 , NPDE = 1, NV = 0,
     1              NPOLY = 6, NPTS = NEL*NPOLY+1,     NXI = 0,
     2              NEQ   = NPTS * NPDE + NV,
     3              NWKRES= (NPOLY+1) * (5*NXI + 3*NPOLY+NEL+5+7*NPDE) +
     4                       NPDE * 8 + 3 + NV + NXI,
C    3              NWKRES= 2*(NPOLY+1)*(NPOLY+NEL+2) + 2 + NV +
C    4                       NPDE * (7 * (NPOLY+1+NXI) + 8),
C     DDASSL TIME INTEGRATION PARAMETERS
     5              MAXORD = 5, LRW = 40 + (MAXORD+4) * NEQ + NEQ**2,
     6              LIW = 20 + NEQ )
C
        INTEGER IWORK(LIW), INFO(15), IBAND, M, ITIME, I, IDID, IRESWK,
     1          IDEV, ITRACE, IDERIV, IFL, ITYPE
        DOUBLE PRECISION XBK(IBK), X(NPTS), Y(NEQ), YDOT(NEQ), Z(NPTS),
     1          WKRES(NWKRES), RWORK(LRW), XI(1), T, TOUT, RTOL, ATOL,
     2          ENORM, GERR, CTIME, DYDX(NEQ), DYCALC(NPDE,NPTS,2)
        EXTERNAL PDECHB, DGEJAC
        COMMON /SDEV2/ ITRACE, IDEV
        COMMON /PROB3/IDERIV
C CPU TIMER COMMENTED OUT FOR PORTABILITY.
C       CALL TIMER ( CTIME, 1)
        M    = 0
        T    = 0.0D0
        IDEV = 4
        ITRACE = 1
        WRITE(IDEV,9)NPOLY, NEL
 9      FORMAT(' TEST PROBLEM 3'/' ***********'/' POLY OF DEGREE =',I4,
     1         ' NO OF ELEMENTS = ',I4)
        DO 10 I = 1,IBK
 10       XBK(I) = -1.0D0 + 2.0D0 * (I-1.0D0)/(IBK -1.0D0)
C           INITIALISE THE P.D.E. WORKSPACE
        ITIME = 1
        CALL INICHB(NEQ, NPDE, NPTS, X, Y, WKRES, NWKRES, M, T, IBAND,
     1              ITIME, XBK, IBK, NEL, NPOLY, NV, NXI, XI, IDEV)
        IF(ITIME .EQ. -1)THEN
           WRITE(IDEV, 15)
 15        FORMAT(' INITCC ROUTINE RETURNED ITIME = -1 - RUN HALTED ')
           GOTO 100
        END IF
C         SETUP DASSL PARAMETERS
       RTOL = 1.0D-5
       ATOL = 1.0D-5
       DO 20 I = 1,11
 20      INFO(I) = 0
C      INFO(11)= 1
C         BANDED MATRIX OPTION WHEN INFO(6) = 1
       IF(INFO(6) .EQ. 1)THEN
          IWORK(1) = IBAND
          IWORK(2) = IBAND
       END IF
       T = 0.0D0
 30    TOUT = T + 0.1D0
       CALL DDASSL( PDECHB, NEQ, T, Y, YDOT, TOUT, INFO, RTOL, ATOL,
     1              IDID, RWORK, LRW, IWORK, LIW, WKRES, IRESWK, DGEJAC)
       IF( IDID .LT. 0 )THEN
C          DASSL FAILED TO FINISH INTEGRATION.
           WRITE(IDEV,40)T,IDID
 40        FORMAT(' AT TIME T = ',D11.3,' DASSL RETURNED IDID =',I3)
           GOTO 100
       ELSE
C        DASSL INTEGRATED TO T = TOUT
C        CALL TO POST PROCESSING HERE E.G. SPACE INTERPOLATION.
         IDERIV = 0
         CALL ERROR( Y, NPDE, NPTS, X, M, ENORM, GERR, T, RTOL, ATOL,
     1               ITRACE, WKRES, NWKRES)
         IFL   = 0
         ITYPE = 2
         DO 45 I = 1,NPTS
45         Z(I) = X(I)
         CALL INTERC(Z,DYCALC,NPTS,Y,NEQ,NPDE,IFL,ITYPE,WKRES,NWKRES)
         IDERIV = 1
         CALL EXACT(T, NPDE, NPTS, X, DYDX)
         DO 50 I = 1,NPTS
          GERRDX = ABS( DYDX(I) - DYCALC(1,I,2))
          WRITE(IDEV,49)X(I),DYDX(I),DYCALC(1,I,2),GERRDX
 49       FORMAT(' X =',D11.3,' TRUE = ',D11.3,' CALC= ',D11.3,' ERR=',
     1           D11.3)
 50      CONTINUE
         IF(TOUT .LT. 0.99D0 ) GOTO 30
       END IF
100    CONTINUE
C      CALL TIMER(CTIME, 2)
       WRITE(IDEV,110)IWORK(11),IWORK(12),IWORK(13), CTIME
110    FORMAT(' NSTEPS =',I5,' NRESID =',I5,' JAC = ',I4,' CPU=',D11.3)
       STOP
       END
C EXAMPLE PROBLEM THREE
C *********************
C     THIS PROBLEM IS DEFINED BY
C               -1
C       U  = ( C   U  )  - C * EXP(-2U) + EXP(-U)  ,  X IN (-1,0)
C        T      1   X  X    1
C  AND
C               -1
C       U  = ( C   U  )  - C * EXP(-2U) + EXP(-U)  ,  X IN (0,1)
C        T      2   X  X    2
C  WHERE
C       C  = 0.1     AND    C   = 1.0
C        1                   2
C
C     THE LEFT BOUNDARY CONDITION AT X =-1 (LEFT = .TRUE. ) IS GIVEN BY
C         U(-1,T)  = LOG ( - C  + T + P)
C                       1
C     THE RIGHT BOUNDARY CONDITION IS  (LEFT = .FALSE.)
C         U( 1,T) + (C + T + P ) U  = LOG ( - C  + T + P) + 1.0D0
C                                 X
C
C      THE INITIAL CONDITION IS GIVEN BY THE EXACT SOLUTION ;
C        U( X, T )  = LOG ( C X + T + P )  , X IN ( -1, 0)
C                            1
C        U( X, T )  = LOG ( C X + T + P )  , X IN (  0, 1)
C                            2
C**********************************************************************
       SUBROUTINE UVINIT( NPDE, NPTS, X, U, NV,V)
C      ROUTINE FOR P.D.E. INITIAL VALUES.
       INTEGER NPDE, NPTS, NV
       DOUBLE PRECISION X(NPTS), U(NPDE,NPTS), V(1), T
         T = 0.0D0
C        V(1) IS A DUMMY VARIABLE AS THERE ARE NO COUPLED O.D.E.S
         CALL EXACT( T, NPDE, NPTS, X, U )
       RETURN
       END
C
       SUBROUTINE SPDEFN( T, X, NPTL, NPDE, U, DUDX, UDOT, UTDX, Q, R,
     1                    NV, V, VDOT, IRES)
C      ROUTINE TO DESCRIBE THE BODY OF THE P.D.E.
C      THE P.D.E. IS WRITEN AS       -M   M
C         Q(X,T,U, U  , U  , U  ) = X   (X  R(X,T,U,U , U , U  ))
C                   X    T    TX                     X   T   TX  X
C      THE FUNCTIONS Q AND R MUST BE DEFINED IN THIS ROUTINE.
C      DEFINITIONS FOR THE MODEL PROBLEM ARE GIVEN
C      NOTE NV = 0 : THERE IS NO O.D.E PART.
       INTEGER NPDE, NPTL, I, J, NV, IRES
       DOUBLE PRECISION T, X(NPTL), U(NPDE,NPTL), DUDX(NPDE,NPTL),
     1         UDOT(NPDE,NPTL), Q(NPDE,NPTL), R(NPDE,NPTL), V, VDOT,
     2         UTDX(NPDE,NPTL), C
       IF(X(1) .LT. 0.0D0 .AND. X(NPTL) .LE. 0.0D0)THEN
C        ELEMENT TO LEFT OF THE INTERFACE AT 0.0
         C  =  0.1D0
       ELSE
         C =   1.0D0
       END IF
       DO 10 I = 1,NPTL
          R(1,I) = DUDX(1,I) /C
          Q(1,I) = UDOT(1,I) - DEXP(-U(1,I))- DEXP(-2.0D0*U(1,I))* C
 10    CONTINUE
       RETURN
       END
C
       SUBROUTINE SBNDR( T, BETA, GAMMA, U, UX, UDOT, UTDX, NPDE, LEFT,
     1                   NV, V, VDOT, IRES)
C      BOUNDARY CONDITIONS ROUTINE
       INTEGER NPDE, NV, IRES
       DOUBLE PRECISION T, BETA(NPDE), GAMMA(NPDE), U(NPDE), C2,
     1                  UX(NPDE), V, VDOT, UDOT(NPDE), UTDX(NPDE)
       LOGICAL LEFT
       IF(LEFT) THEN
          BETA (1) = 0.0D0
          GAMMA(1) = U(1) - DLOG( -0.1 + T + 1.0D0)
       ELSE
          C2 = 1.0D0
          BETA (1) = C2 * ( C2 + T + 1.0D0)
          GAMMA(1) = U(1) - DLOG( C2 + T + 1.0D0) + 1.0D0
       END IF
       RETURN
       END
C
C      DUMMY O.D.E. ROUTINE AS NV IS ZERO
       SUBROUTINE SODEFN
       RETURN
       END
C EXACT SOLUTION
       SUBROUTINE EXACT( TIME, NPDE, NPTS, X, U)
C      ROUTINE FOR P.D.E. EXACT VALUES  (IF KNOWN)
       INTEGER NPDE, NPTS, I, IDERIV
       DOUBLE PRECISION X(NPTS), U(NPDE,NPTS), TIME, C
       COMMON /PROB3/ IDERIV
       IF(IDERIV .EQ. 0)THEN
          DO 10 I = 1,NPTS
             C = 1.0D0
             IF(X(I) .LT. 0.0D0)C = 0.1D0
 10          U(1,I) = DLOG( TIME + 1.0D0 + C * X(I))
       ELSE
          DO 20 I = 1,NPTS
             C = 1.0D0
             IF(X(I) .LT. 0.0D0)C = 0.1D0
             U(1,I) = C / ( TIME + 1.0D0 + C * X(I))
             IF(X(I) .EQ. 0.0D0) U(1,I) = 0.55D0 / ( TIME + 1.0D0 )
 20       CONTINUE
       END IF
       RETURN
       END
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
c  problem 5
C     ***********************************************************
C     FOURTH ORDER P.D.E. PROBLEM WRITTEN AS ELLIPTIC-PARABOLIC SYSTEM.
C
C     U     =  K U     + UU    - U U
C      XXT        XXXX     XXX    X XX
C
C     ***********************************************************
C
C     C0 COLLOCATION PARAMETERS
        PARAMETER ( IBK   = 21, NEL  = IBK-1 , NPDE = 2, NV = 0,
     1              NPOLY = 02, NPTS = NEL*NPOLY+1,     NXI = 0,
     2              NEQ   = NPTS * NPDE + NV,
     3              NWKRES= (NPOLY+1) * (5*NXI + 3*NPOLY+NEL+5+7*NPDE) +
     4                       NPDE * 8 + 3 + NV + NXI,
C     DDASSL TIME INTEGRATION PARAMETERS
     5              MAXORD = 5, LRW = 40 + (MAXORD+4) * NEQ + NEQ**2,
     6              LIW = 20 + NEQ )
C
        INTEGER IWORK(LIW), INFO(15), IBAND, M, ITIME, I, IDID,
     1          IDEV, ITRACE, GRNPTS, IFL, NOUT, KTIME, ITYPE, NP
        DOUBLE PRECISION XBK(IBK), X(NPTS), Y(NEQ), YDOT(NEQ), TINC(15),
     1          WKRES(NWKRES), RWORK(LRW), XI, T, TOUT, RTOL, ATOL,
     3          TEND, K, XOUT(6), UOUT(2,6)
        EXTERNAL PDECHB, DGEJAC
C
C       COMMON BLOCKS  TO PASS ACROSS PROBLEM DEPENDENT CONSTANTS.
        COMMON /PDES/   K
        COMMON /SDEV2/ ITRACE, IDEV
        DATA XOUT(1)/-1.0D+0/, XOUT(2)/-0.6D+0/, XOUT(3)/-0.2D+0/,
     *  XOUT(4)/0.2D+0/, XOUT(5)/0.6D+0/, XOUT(6)/1.0D+0/
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
55       FORMAT(//' SOLUTION TO FOURTH ORDER P.D.E. PROBLEM USING
     1   DASSL INTEGRATOR WITH BANDED MATRIX ROUTINES '/
     2   '   ATOL = ',D11.3,'  RTOL = ',D11.3,'  NPTS = ',I5/)
C
C    EQUALLY SPACED BREAKPOINTS.
C
         DO 105 I = 1,IBK
          XBK(I) = -1.0D0 + (I -1.0D0)* 2.D0 / (IBK-1.D0)
 105     CONTINUE
          K      =  1.00D0
          ITIME = 1
          T     = 0.0D0
C           INITIALISE THE P.D.E. WORKSPACE
        CALL INICHB(NEQ, NPDE, NPTS, X, Y, WKRES, NWKRES, M, T, IBAND,
     1              ITIME, XBK, IBK, NEL, NPOLY, NV, NXI, XI, IDEV)
        IF(ITIME .EQ. -1)THEN
           WRITE(IDEV, 15)
 15        FORMAT(' INICHB ROUTINE RETURNED ITIME = -1 - RUN HALTED ')
           GOTO 900
        END IF
C         SETUP DASSL PARAMETERS
       DO 20 I = 1,11
 20      INFO(I) = 0
       INFO(6)= 1
       INFO(9)= 1
       INFO(7)= 1
       IWORK(3)= 4
C         BANDED MATRIX OPTION WHEN INFO(6) = 1
       IF(INFO(6) .EQ. 1)THEN
          IWORK(1) = IBAND
          IWORK(2) = IBAND
       END IF
      T    = 0.0D0
      TINC(1) = 0.0001D0
      RWORK(2)= TINC(1) * 0.1D0
      TINC(2) = 0.0010
      TINC(3) = 0.01D0
      TINC(4) = 0.1D0
      TINC(5) = 1.0D0
      TINC(6) = 1.00D1
      TINC(7) = 2.00D1
      TINC(8) = 4.00D1
      TINC(9 )= 6.00D1
      TINC(10)= 8.00D1
      TINC(11)= 1.00D2
      TINC(12)= 1.00D3
      TEND = 1.0D3
      KTIME = 1
      WRITE(IDEV,83)(XOUT(I),I = 1,6)
C TIME LOOP:
100   TOUT = TINC(KTIME)
      IF(KTIME.GT.1)RWORK(2) = 0.05D0 *(TOUT- TINC(KTIME-1))
      IF(KTIME .EQ.12)THEN
         INFO(4) = 1
         RWORK(1) = TEND
       END IF
C
       CALL DDASSL( PDECHB, NEQ, T, Y, YDOT, TOUT, INFO, RTOL, ATOL,
     1              IDID, RWORK, LRW, IWORK, LIW, WKRES, NWKRES, DGEJAC)
C          DASSL FAILED TO FINISH INTEGRATION.
       WRITE(IDEV,40)T,IDID,Y(1),Y(2),Y(NEQ-1), Y(NEQ)
 40    FORMAT(' AT TIME T = ',D11.3,' DASSL RETURNED IDID =',I3/
     1        ' LEFT SOL=',2D11.3,' RIGHT SOL=',2D11.3)
       IF( IDID .LT. 0 )GOTO 900
C        DASSL INTEGRATED TO T = TOUT
C        CALL TO POST PROCESSING HERE E.G. SPACE INTERPOLATION.
         ITYPE = 1
         NP    = 6
         CALL INTERC( XOUT, UOUT, NP, Y, NEQ, NPDE, IFLAG,
     1                      ITYPE, WKRES, NWKRES)
       WRITE(IDEV,82)(UOUT(1,I),I = 1,6)
       WRITE(IDEV,84)(UOUT(2,I),I = 1,6)
 83    FORMAT(1X,'X',9D11.3/)
 82    FORMAT(1X,'U',9D11.3/)
 84    FORMAT(1X,'V',9D11.3/)
C
 91   KTIME = KTIME + 1
C
C     CHECK IF INTEGRATION WAS SUCCESSFUL AND WHETHER FURTHER TIME
C     STEPS NEEDED
      IF(TOUT.LT.TEND.AND.(IDID.EQ.2 .OR. IDID .EQ. 3)) GO TO 100
80    CONTINUE
900    WRITE(IDEV,110)IWORK(11),IWORK(12),IWORK(13)
110    FORMAT(' NSTEPS =',I5,' NRESID =',I5,' JAC = ',I4)
       STOP
       END
C EXAMPLE PROBLEM FIVE
C *********************
C     THIS PROBLEM IS DEFINED BY
C
C       V  =   U
C               XX
C  AND
C
C       V  = ( K   V  )  + U V  - U V              ,  X IN (-1 , 1)
C        T          X  X      X    X
C  WHERE
C       K  = 0.15
C
C     THE LEFT BOUNDARY CONDITION AT X =-1 (LEFT = .TRUE. ) ARE GIVEN BY
C
C       U = 1    U  = 0.0
C                 X
C     THE RIGHT BOUNDARY CONDITION ARE (LEFT = .FALSE.)
C
C       U = -1   U  = 0.0D0
C                 X
C      THE INITIAL CONDITION IS GIVEN BY
C
C         U(X,0) = -SIN ( PI /2  X )
C**********************************************************************
       SUBROUTINE UVINIT( NPDE, NPTS, X, U, NV,V)
C      ROUTINE FOR P.D.E. INITIAL VALUES.
       PARAMETER (PIBY2 = 1.5707963D0)
       INTEGER NPDE, NPTS, NV, I
       DOUBLE PRECISION X(NPTS), U(NPDE,NPTS),V
       DO 10 I= 1,NPTS
          U(1,I) = -SIN( PIBY2 * X(I) )
 10       U(2,I) = - PIBY2**2 * U(1,I)
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
     1    UDOT(NPDE,NPTL), Q(NPDE,NPTL), R(NPDE,NPTL), V,
     2    K, UTDX(NPDE,NPTL), VDOT
      COMMON /PDES/ K
      DO 100 I = 1,NPTL
          Q(1,I) = U(2,I)
          R(1,I) = DUDX(1,I)
          Q(2,I) = UDOT(2,I) - U(1,I)*DUDX(2,I) + DUDX(1,I)*U(2,I)
          R(2,I) = K*DUDX(2,I)
 100  CONTINUE
      RETURN
      END
C
       SUBROUTINE SBNDR( T, BETA, GAMMA, U, UX, UDOT, UTDX, NPDE, LEFT,
     1                   NV, V, VDOT, IRES)
C      BOUNDARY CONDITIONS ROUTINE
       INTEGER NPDE, NV, IRES
       DOUBLE PRECISION T, BETA(NPDE), GAMMA(NPDE), U(NPDE),
     1         UX(NPDE), V, VDOT, UDOT(NPDE), UTDX(NPDE)
      LOGICAL LEFT
      IF(LEFT) THEN
          GAMMA(1) = 0.0D0
          BETA(1) = 1.0D+0
          GAMMA(2) = U(1) - 1.0D0
          BETA(2) = 0.0D+0
      ELSE
          GAMMA(1) = 0.0D0
          BETA(1) = 1.0D+0
          GAMMA(2) = U(1) + 1.0D0
          BETA(2) = 0.0D+0
      END IF
      RETURN
      END
      SUBROUTINE SODEFN(T, NV, V, VDOT, NPDE, NXI, X, Y, UXI, RI,
     1                  UTI, UTXI, VRES, IRES)
C ROUTINE FOR AUXILIARY O.D.E.S (IF ANY) IN MASTER EQN. FORM (4.3)
      INTEGER NPDE, NXI, NV, IRES, NPTL, L, J, I
      DOUBLE PRECISION T, X(NXI), Y(NXI), UXI(NPDE,NXI),
     1        RI(NPDE,NXI), UTI(NPDE,NXI), UTXI(NPDE,NXI), VRES(NV),
     2        V, VDOT
C                     DUMMY ROUTINE AS THERE ARE NO O.D.E.S
      RETURN
      END