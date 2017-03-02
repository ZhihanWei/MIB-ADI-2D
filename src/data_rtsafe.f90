      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! RTNEWT --
      !    NEWTON-RAPHSON METHOD
      ! ARGUMENTS:
      !    FUNCD   IN   EXTERNAL FUNCTION DEFINING F AND DF
      !    Y       IN   A GIVEN AXIS
      !    X1      IN   LEFT END OF THE INTERVAL
      !    X2      IN   RIGHT END OF THE INTERVAL
      !    XACC    IN   TOLERANCE TO STOP THE ITERATION
      !    RTNEWT  OUT  RESULT OBTAINED BY NEWTON ITERATION
      ! NOTES:
      !    RTNEWT = RTNEWT - F/DF
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION RTNEWT(FUNCD,Y,X1,X2,XACC)

      REAL Y,X1,X2,XACC,RTNEWT
      EXTERNAL FUNCD
      INTEGER,PARAMETER :: JMAX=20 !MAX NUMBER OF ITERATIONS

      REAL F,DF
      INTEGER J
      REAL DX !LOCAL DECLAIRATIONS TO AVOID CONFLICT WITH GLOBAL DEFINITIONS

      WRITE (*,*) Y,X1/ACOS(-1.D0),X2/ACOS(-1.D0)

      RTNEWT = (X1+X2)*0.5D0

      DO J = 1,JMAX
         CALL FUNCD(RTNEWT,Y,F,DF)
         DX = F/DF
         RTNEWT = RTNEWT-DX
         IF ((X1-RTNEWT)*(RTNEWT-X2) .LT. 0) THEN
            WRITE (*,*) "RTNEWT JUMPED OUT OF BRACKETS"
            STOP
         END IF
         IF (ABS(DX) .LT. XACC) RETURN!CONVERGENCE
      END DO
      WRITE (*,*) "RTNEWT EXCEEDED MAX ITERATIONS WITH ERROR", ABS(DX)
      STOP

      END FUNCTION RTNEWT

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! RTSAFE --
      !    A COMBINATION OF NEWTON-RAPHSON METHOD AND BISECTION
      ! ARGUMENTS:
      !    FUNCD   IN   EXTERNAL FUNCTION DEFINING F AND DF
      !    Y       IN   A GIVEN AXIS
      !    X1      IN   LEFT END OF THE INTERVAL
      !    X2      IN   RIGHT END OF THE INTERVAL
      !    XACC    IN   TOLERANCE TO STOP THE ITERATION
      !    RTSAFE  OUT  RETURNED RESULT
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION RTSAFE(FUNCD,Y,X1,X2,XACC)

      REAL RTSAFE,Y,X1,X2,XACC

      EXTERNAL FUNCD
      INTEGER,PARAMETER :: MAXIT=20!MAX NUMBER OF ITERATIONS

      REAL FL,DF,FH,F
      REAL DXOLD,TEMP,XH
      REAL DX,XL !LOCAL DECLAIRATIONS TO AVOID CONFLICT WITH GLOBAL DEFINITIONS
      INTEGER J

      CALL FUNCD(X1,Y,FL,DF)
      CALL FUNCD(X2,Y,FH,DF)

      IF (((FL .GT. 0.D0) .AND. (FH .GT. 0.D0)) .OR. ((FL .LT. 0.D0) .AND. (FH .LT. 0))) THEN
         WRITE (*,*) "ERROR: ROOT MUST BE BRACKETED IN RTSAFE",Y,X1,X2
         STOP
      END IF
      IF (FL .EQ. 0.D0) THEN
         RTSAFE = X1
         RETURN
      ELSE IF (FH .EQ. 0.D0) THEN
         RTSAFE = X2
         RETURN
      ELSE IF (FL .LT. 0.0D0) THEN !ORIENT THE SEARCH SO THAT F(XL)<0
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
      END IF

      RTSAFE = 0.5D0*(X1+X2)!INITIALIZE THE GUESS FOR ROOT
      DXOLD  = ABS(X2-X1)!THE "STEPSIZE BEFORE LAST,"
      DX    = DXOLD !AND THE LAST STEP
      CALL FUNCD(RTSAFE,Y,F,DF)

      DO J = 1,MAXIT
         !BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING FAST ENOUGH
         IF ((((RTSAFE-XH)*DF-F)*((RTSAFE-XL)*DF-F) .GT. 0.D0) .OR. (ABS(2.D0*F) .GT. ABS(DXOLD*DF))) THEN
            DXOLD  = DX
            DX     = 0.5D0*(XH-XL)
            RTSAFE = XL+DX
            IF (XL .EQ. RTSAFE) THEN
               RETURN !CHANGE IN ROOT IS NEGLIGIBLE
            END IF
         ELSE
            DXOLD  = DX
            DX     = F/DF
            TEMP   = RTSAFE
            RTSAFE = RTSAFE-DX
            IF (TEMP .EQ. RTSAFE) THEN
               RETURN !CHANGE IN ROOT IS NEGLIGIBLE
            END IF
         END IF

         IF (ABS(DX) .LT. XACC) THEN
            RETURN!CONVERGENCE CRITERION
         END IF

         CALL FUNCD (RTSAFE,Y,F,DF)!THE ONE NEW FUNCTION EVALUATION PER ITERATION
         IF (F .LT. 0.0D0) THEN !MAINTIN THE BRACKET ON THE ROOT
            XL = RTSAFE
         ELSE
            XH = RTSAFE
         END IF
      END DO

      WRITE (*,*) "ERROR: RTSAFE EXCEEDED MAX ITERATIONS WITH ERROR",ABS(DX)
      STOP

      END FUNCTION RTSAFE
