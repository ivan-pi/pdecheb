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
     1          IDEV, ITRACE, GRNPTS, IFL, NOUT, KTIME, ITYPE, NP, NEQN
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
         IDEV = 6
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
        CALL INICHB(NEQN, NPDE, NPTS, X, Y, WKRES, NWKRES, M, T, IBAND,
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