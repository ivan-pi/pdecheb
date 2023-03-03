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
