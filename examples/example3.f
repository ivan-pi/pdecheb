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
C
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
     1          IDEV, ITRACE, IDERIV, IFL, ITYPE, NEQN
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
        IDEV = 6
        ITRACE = 1
        WRITE(IDEV,9)NPOLY, NEL
 9      FORMAT(' TEST PROBLEM 3'/' ***********'/' POLY OF DEGREE =',I4,
     1         ' NO OF ELEMENTS = ',I4)
        DO 10 I = 1,IBK
 10       XBK(I) = -1.0D0 + 2.0D0 * (I-1.0D0)/(IBK -1.0D0)
C           INITIALISE THE P.D.E. WORKSPACE
        ITIME = 1
        CALL INICHB(NEQN, NPDE, NPTS, X, Y, WKRES, NWKRES, M, T, IBAND,
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
       
