      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2D HEAT EQUATION WITH MULTI-LEAVES INTERFACE
      ! ============================================
      ! PDE: U_T = (BETA U_X)_X + (BETA U_Y)_Y + F
      !    BETA1=1     INSIDE
      !    BETA2=10    OUTSIDE
      !    U(X,Y)=SIN(KX)COS(KY)COS(T) INSIDE
      !          =COS(KX)SIN(KY)COS(T) OUTSIDE
      !    F=(2*K^2*BETA1*COS(T)-SIN(T))*SIN(KX)COS(KY)     INSIDE
      !     =(2*K^2*BETA2*COS(T)-SIN(T))*COS(KX)SIN(KY)     OUTSIDE
      !
      ! GAMMA: R = A + B*SIN(K*VARPHI)
      !    VARPHI   THE ARC-LENGTH PARAMETER,
      !    THETA    THE ANGLE OF OUTWARD NORMAL OF AN INTERFACE POINT (XINF,YINF)
      !    VARPHI   [0,2*PI]
      !    (2-LEAVE: A=1/2, B=1/4, K=2; 4-LEAVE: A=1/2, B=1/10, K=4;)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !---------- TEST CASE IDENTIFIER ----------
      !INTEGER :: TESTCASE = 5

      !---------- INTERFACE DEFINITION ----------
      REAL :: STAR_A = 8.0D-1
      REAL :: STAR_B = 3.0D-1
      REAL :: STAR_K = 5.0D0

      !---------- EXACT SOLUTION ----------
      REAL :: VK = 2.0D0 !WAVE NUMBER K

      CONTAINS

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                                                                        !
      !                                INTERFACE                               !
      !                                                                        !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SETS --
      !    LEVEL SET FUNCTION
      ! ARGUMENTS:
      !    X     IN   X-COORDINATE OF A GIVEN POINT
      !    Y     IN   Y-COORDINATE OF A GIVEN POINT
      !    SETS  OUT  SINGED DISTANCE OF A POINT (X,Y) TO THE INTERFACE
      ! NOTES:
      !    RETURNS A POSITIVE VALUE IF THE POINT (X,Y) STAYS OUTSIDE THE GIVEN
      !    INTERFACE, NEGATIVE IF INSIDE
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION SETS(X,Y)

      REAL :: X,Y,SETS

      REAL :: RXY,VARPHI

      RXY    = SQRT(X**2+Y**2)
      VARPHI = GETANGLE(X,Y,RXY)
      SETS   = RXY-(STAR_A+STAR_B*SIN(STAR_K*VARPHI))+TOL_SETUP

      RETURN

      END FUNCTION SETS

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                                                                        !
      !                                EQUATION                                !
      !                                                                        !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ANALYTICAL --
      !    GENERATING ANALYTICAL SOLUTION
      ! ARGUMENTS:
      !    U  OUT  MATRIX TO STORE EXACT SOLUTION ON ALL GRIDS
      !    T  IN   CURRENT TIME
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ANALYTICAL(U,T)

      REAL :: U(NY,NX),T

      REAL :: CT
      INTEGER :: IX,IY

      CT = COS(T)
      DO IY = 1,NY
         DO IX = 1,NX
            IF (SETS(XI(IX),YI(IY)) .GT. TOL_SETUP) THEN
               U(IY,IX) = COS(VK*XI(IX))*SIN(VK*YI(IY))*CT
            ELSE
               U(IY,IX) = SIN(VK*XI(IX))*COS(VK*YI(IY))*CT
            END IF
         END DO
      END DO

      RETURN

      END SUBROUTINE ANALYTICAL

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SETBETA --
      !    ON-GRID BETA VALUES INSIDE AND OUTSIDE THE INTERFACE
      ! ARGUMENTS:
      !
      ! NOTES:
      !    BETA IS A FUNCTION OF X AND Y, I.E., BETA(X,Y).
      !    THREE BETA VALUES ARE NEEDED ON EACH GRID: BETA,(dBETA/dX)/BETA,
      !    AND (dBETA/dY)/BETA.
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SETBETA

      INTEGER :: IX,IY

      DO IY = 1,NY
         DO IX = 1,NX
            IF (SETS(XI(IX),YI(IY)) .GT. TOL_SETUP) THEN
               BETA(IY,IX)  = 10.0D0 !ON-GRID BETA OUTSIDE
               BETAX(IY,IX) =  0.0D0 !(dBETA/dX)/BETA
               BETAY(IY,IX) =  0.0D0 !(dBETA/dY)/BETA
            ELSE
               BETA(IY,IX) =   1.0D0 !ON-GRID BETA INSIDE
               BETAX(IY,IX) =  0.0D0 !(dBETA/dX)/BETA
               BETAY(IY,IX) =  0.0D0 !(dBETA/dY)/BETA
            END IF
         END DO
      END DO

      RETURN

      END SUBROUTINE SETBETA

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SETIFPBETA --
      !    SET INTERFACE BETA VALUES
      ! ARGUMENTS:
      !    DATA     OUT   DATA OF A LISTED INTERFACE POINT
      ! NOTES:
      !    BETA IS A FUNCTION OF X AND Y IN OMEGA^- AND OMEGA^+
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SETIFPBETA(DATA)

      TYPE(LIST_DATA) :: DATA

      DATA%BETA1 = 1.0D0
      DATA%BETA2 = 1.0D1

      RETURN

      END SUBROUTINE SETIFPBETA

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! IFBETA1 --
      !    INTERFACE BETA VALUE IN OMEGA^-
      ! ARGUMENTS:
      !    XINF     IN   X-COORDINATE OF THE INTERFACE POINT
      !    YINF     IN   Y-COORDINATE OF THE INTERFACE POINT
      !    IFBETA1  OUT  BETA^- VALUE ON THIS INTERFACE POINT
      ! NOTES:
      !    BETA IS A FUNCTION OF X AND Y
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      FUNCTION IFBETA1(XINF,YINF)
!
!      REAL :: XINF,YINF,IFBETA1
!
!      IFBETA1 = 1.0D0
!
!      RETURN
!
!      END FUNCTION IFBETA1

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! IFBETA2 --
      !    INTERFACE BETA VALUE IN OMEGA^+
      ! ARGUMENTS:
      !    XINF     IN   X-COORDINATE OF THE INTERFACE POINT
      !    YINF     IN   Y-COORDINATE OF THE INTERFACE POINT
      !    IFBETA2  OUT  BETA^+ VALUE ON THIS INTERFACE POINT
      ! NOTES:
      !    BETA IS A FUNCTION OF X AND Y
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      FUNCTION IFBETA2(XINF,YINF)
!
!      REAL :: XINF,YINF,IFBETA2
!
!      IFBETA2 = 10.0D0
!
!      RETURN
!
!      END FUNCTION IFBETA2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SETSRC --
      !    ON-GRID SOURCE TERM AT TIME T
      ! ARGUMENTS:
      !    T  IN  CURRENT TIME
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SETSRC(T)

      REAL :: T

      INTEGER :: IX,IY
      REAL :: CT

      DO IY=1,NY
         DO IX=1,NX
            CT = 2.0D0*VK**2*BETA(IY,IX)*COS(T)-SIN(T)
            IF (SETS(XI(IX),YI(IY)) .GT. TOL_SETUP) THEN  !OUTSIDE
               SRC(IY,IX)=CT*COS(VK*XI(IX))*SIN(VK*YI(IY))
            ELSE                                          !INSIDE
               SRC(IY,IX)=CT*SIN(VK*XI(IX))*COS(VK*YI(IY))
            END IF
         END DO
      END DO

      RETURN

      END SUBROUTINE SETSRC

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SETBC --
      !    SET AVERAGE BOUNDARY CONDITIONS BETWEEN TIME STEP T AND T+DT
      ! ARGUMENTS:
      !    T    IN   CURRENT TIME
      !    DT   IN   TIME STEP
      !    UHS  OUT  ON-GRID NUMERICAL SOLUTIONS
      ! NOTES:
      !    RESET BOUNDARY CONDITIONS AT T_N      : T = T_N,     DT = 0
      !                                 T_{N+1/2}: T = T_N,     DT = DT
      !                                 T_{N+1}  : T = T_{N+1}, DT = 0
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SETBC(T,DT,UHS)

      REAL :: T,DT,UHS(NY,NX)

      REAL :: CT
      INTEGER :: IX,IY

      CT=(COS(T+DT)+COS(T))/2.0D0

      DO IX=1,NX,NX-1    !FIRST ORDER BOUNDARY CONDITION
         DO IY=1,NY      !U*=(U^N+ U^{N+1})/2
            UHS(IY,IX)=COS(VK*XI(IX))*SIN(VK*YI(IY))*CT
         END DO
      END DO

      DO IY=1,NY,NY-1
         DO IX=1,NX
            UHS(IY,IX)=COS(VK*XI(IX))*SIN(VK*YI(IY))*CT
         END DO
      END DO

      END SUBROUTINE SETBC

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                                                                        !
      !                   SETUPS FOR INTERFACE GAMMA IN MIB                    !
      !                                                                        !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! XDX --
      !    COMPUTER F AND DF FOR GAMMAY
      ! ARGUMENTS:
      !    VARPHI  IN   ARC-LENGTH PARAMETER
      !    X       IN   X-COORDINATE OF THE STARTING POINT
      !    F       OUT  FUNCTION F
      !    DF      OUT  DERIVATIVE OF THE FUNCTION F
      ! NOTES:
      !    CALCULATE F AND DF USED IN NEWTON ITERFATION (X DIRECTION)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE XDX(VARPHI,X,F,DF)

      REAL :: VARPHI,X,F,DF

      F  =  (STAR_A+STAR_B*SIN(STAR_K*VARPHI))*COS(VARPHI)-X
      DF = -(STAR_A+STAR_B*SIN(STAR_K*VARPHI))*SIN(VARPHI)+STAR_B*STAR_K*COS(STAR_K*VARPHI)*COS(VARPHI)

      RETURN

      END SUBROUTINE XDX

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! GAMMAX --
      !    LOCATE ON-GRID INTERFACE POINT X COORDINATE
      ! ARGUMENTS:
      !    Y       IN   Y-COORDINATE OF A IY-AXIS
      !    X1      IN   X-COORDINATE OF LEFT POINT OF THE INTERVAL
      !    X2      IN   X-COORDINATE OF RIGHT POINT OF THE INTERVAL
      !    X       OUT  X-COORDINATE OF THE INTERSECTION
      !    VARPHI  OUT  ARC-LENGTH PARAMETER
      ! NOTES:
      !    NEWTON ITERATION TO LOCATE THE INTERSECTION OF THE INTERFACE AND
      !    IY-AXIS IN THE INTERVAL (X1,X2)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE GAMMAX(Y,X1,X2,X,VARPHI)

      REAL :: Y,X1,X2,X,VARPHI

      REAL :: R1,R2,VARPHI1,VARPHI2

      R1   = SQRT(X1**2+Y**2)
      R2   = SQRT(X2**2+Y**2)
      VARPHI1 = GETANGLE(X1,Y,R1)
      VARPHI2 = GETANGLE(X2,Y,R2)

      IF (ABS(VARPHI1-VARPHI2) .LT. TOL_RTSAFE) THEN
         VARPHI = (VARPHI1+VARPHI2)/2.0D0
      ELSE IF (VARPHI1 .LT. VARPHI2) THEN
         VARPHI = RTSAFE(YDY,Y,VARPHI1,VARPHI2,TOL_RTSAFE)
      ELSE
         VARPHI = RTSAFE(YDY,Y,VARPHI2,VARPHI1,TOL_RTSAFE)
      END IF

      X = (STAR_A+STAR_B*SIN(STAR_K*VARPHI))*COS(VARPHI)
      IF (ABS(X-MIN(X1,X2)) .LT. 1.0D-9) THEN
         X = MIN(X1,X2)+ABS(X2-X1)*1.0D-9
         WRITE (*,*) "WARNING: TRANSLATION IN X",X,Y
      END IF

      IF ((X2-X)*(X-X1) .LT. 0.0D0) THEN
         WRITE (*,*) "ERROR: GAMMAX IS NOT WITHIN THE INTERVAL"
         WRITE (*,*) X1,X,X2
         STOP
      END IF

      RETURN

      END SUBROUTINE GAMMAX

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! YDY --
      !    COMPUTER F AND DF FOR GAMMAX
      ! ARGUMENTS:
      !    VARPHI  IN   ARC-LENGTH PARAMETER
      !    Y       IN   Y-COORDINATE OF THE STARTING POINT
      !    F       OUT  FUNCTION F
      !    DF      OUT  DERIVATIVE OF THE FUNCTION F
      ! NOTES:
      !    CALCULATE F AND DF USED IN NEWTON ITERFATION (Y DIRECTION)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE YDY(VARPHI,Y,F,DF)

      REAL :: VARPHI,Y,F,DF

      F  = (STAR_A+STAR_B*SIN(STAR_K*VARPHI))*SIN(VARPHI)-Y
      DF = (STAR_A+STAR_B*SIN(STAR_K*VARPHI))*COS(VARPHI)+STAR_B*STAR_K*COS(STAR_K*VARPHI)*SIN(VARPHI)

      RETURN

      END SUBROUTINE YDY

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! GAMMAY --
      !   LOCATE ON-GRID INTERFACE POINT Y COORDINATE
      ! ARGUMENTS:
      !    X       IN   X-COORDINATE OF A IX-AXIS
      !    Y1      IN   Y-COORDINATE OF BOTTOM POINT OF THE INTERVAL
      !    Y2      IN   Y-COORDINATE OF TOP POINT OF THE INTERVAL
      !    Y       OUT  Y-COORDINATE OF THE INTERSECTION
      !    VARPHI  OUT   ARC-LENGTH PARAMETER
      ! NOTES:
      !    NEWTON ITERATION TO LOCATE THE INTERSECTION OF THE INTERFACE AND
      !    IX-AXIS IN THE INTERVAL (Y1,Y2)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE GAMMAY(X,Y1,Y2,Y,VARPHI)

      REAL :: X,Y1,Y2,Y,VARPHI

      REAL :: R1,R2,VARPHI1,VARPHI2

      R1   = SQRT(X**2+Y1**2)
      R2   = SQRT(X**2+Y2**2)
      VARPHI1 = GETANGLE(X,Y1,R1)
      VARPHI2 = GETANGLE(X,Y2,R2)

      IF (ABS(VARPHI1-VARPHI2) .GT. PI) THEN !ERROR COULD HAPPEN AROUND VARPHI=0
         WRITE (*,*) "GAMMAY: Y CHANGED SIGN",X,Y1,Y2
         IF (VARPHI1 .GT. PI) THEN !USE PERIODIC VARPHI VALUE TO OVERIDE
            VARPHI1 = VARPHI1-2.0D0*PI
         ELSE
            VARPHI2 = VARPHI2-2.0D0*PI
         END IF
      END IF

      IF (ABS(VARPHI1-VARPHI2) .LT. TOL_RTSAFE) THEN
         VARPHI = (VARPHI1+VARPHI2)/2.0D0
      ELSE IF (VARPHI1 .LT. VARPHI2) THEN
         VARPHI = RTSAFE(XDX,X,VARPHI1,VARPHI2,TOL_RTSAFE)
      ELSE
         VARPHI = RTSAFE(XDX,X,VARPHI2,VARPHI1,TOL_RTSAFE)
      END IF

      Y = (STAR_A+STAR_B*SIN(STAR_K*VARPHI))*SIN(VARPHI)
      IF (ABS(Y-MIN(Y1,Y2)) .LT. 1.0D-9) THEN
         Y = MIN(Y1,Y2)+ABS(Y2-Y1)*1.0D-9
         WRITE (*,*) "WARNING: TRANSLATION IN Y",X,Y
      END IF

      IF ((Y2-Y)*(Y-Y1) .LT. 0.D0) THEN
         WRITE (*,*) "ERROR: GAMMAY IS NOT WITHIN THE INTERVAL"
         WRITE (*,*) Y1,Y,Y2
         STOP
      END IF

      RETURN

      END SUBROUTINE GAMMAY

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! GETNORMAL --
      !    DETERMINE THE ANGLE FOR OUTWARD NORMAL DIRECTION
      ! ARGUMENTS:
      !    XO         IN   X-COORIDNATE OF A GIVEN POINT
      !    YO         IN   Y-COORDINATE OF A GIVEN POINT
      !    VARPHI     IN   ARC-LENGTH PARAMETER
      !    GETNORMAL  OUT  THE ANGLE FORMED BY THE OUTWARD NORMAL DIRECTION AND
      !                    X-AXIS
      ! NOTES:
      !    THE FORMULAR IS GIVEN BY ???
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION GETNORMAL(XO,YO,VARPHI)

      REAL :: XO,YO,VARPHI,GETNORMAL

      REAL :: DEN,XNOR,YNOR,R

      DEN  = SQRT(STAR_A**2+2.0D0*STAR_A*STAR_B*SIN(STAR_K*VARPHI)+&
             STAR_B**2-STAR_B**2*COS(STAR_K*VARPHI)**2+STAR_B**2*STAR_K**2*COS(STAR_K*VARPHI)**2)
      XNOR = (STAR_A*COS(VARPHI)+STAR_B*COS(VARPHI)*SIN(STAR_K*VARPHI)+&
              STAR_B*STAR_K*SIN(VARPHI)*COS(STAR_K*VARPHI))/DEN
      YNOR = (STAR_A*SIN(VARPHI)+STAR_B*SIN(VARPHI)*SIN(STAR_K*VARPHI)-&
              STAR_B*STAR_K*COS(VARPHI)*COS(STAR_K*VARPHI))/DEN

      R = SQRT(XNOR**2+YNOR**2)
      GETNORMAL = GETANGLE(XNOR,YNOR,R)

      IF (WARNINGS) THEN
         WRITE(*,*) "XO = ",XO,"IS UNUSED IN GETNORMAL"
         WRITE(*,*) "YO = ",YO,"IS UNUSED IN GETNORMAL"
      END IF

      RETURN

      END FUNCTION GETNORMAL

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                                                                        !
      !                     INTERFACE JUMP CONDITIONS                          !
      !                                                                        !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! PHI --
      !    PHI = [U] = U^+ - U^-
      ! ARGUMENTS:
      !    T      IN       CURRENT TIME
      !    DATA   IN/OUT   A LISTED INTERFACE POINT, DATA%JUMP(1) IS UPDATED
      ! NOTES:
      !    ZEROTH ORDER JUMP CONDITION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE PHI(T,DATA)

      REAL :: T
      TYPE(LIST_DATA) :: DATA

      REAL :: XINF,YINF

      XINF  = DATA%X
      YINF  = DATA%Y

      DATA%JUMP(1) = (COS(VK*XINF)*SIN(VK*YINF)-SIN(VK*XINF)*COS(VK*YINF))*COS(T)

      RETURN

      END SUBROUTINE PHI

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! PSI --
      !    PSI = [BETA U_N] = BETA^+ U_N^+ - BEAT^- U_N^-
      ! ARGUMENTS:
      !    T      IN       CURRENT TIME
      !    DATA   IN/OUT   A LISTED INTERFACE POINT, DATA%JUMP(2) IS UPDATED
      ! NOTES:
      !    FISRT ORDER JUMP CONDITION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE PSI(T,DATA)

      REAL :: T
      TYPE(LIST_DATA) :: DATA

      REAL :: XINF,YINF,THETA,BETA1,BETA2

      XINF  = DATA%X
      YINF  = DATA%Y
      THETA = DATA%THETA
      BETA1 = DATA%BETA1
      BETA2 = DATA%BETA2

      DATA%JUMP(2) = (SIN(VK*XINF)*SIN(VK*YINF)*VK*(BETA1*SIN(THETA)-BETA2*COS(THETA))+&
                    &COS(VK*XINF)*COS(VK*YINF)*VK*(BETA2*SIN(THETA)-BETA1*COS(THETA)))*COS(T)

      RETURN

      END SUBROUTINE PSI

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! PHI_TAU --
      !    PHI_TAU = [U_TAU] = U_TAU^+ - U_TAU^-
      ! ARGUMENTS:
      !    T      IN       CURRENT TIME
      !    DATA   IN/OUT   A LISTED INTERFACE POINT, DATA%JUMP(3) IS UPDATED
      ! NOTES:
      !    FISRT ORDER JUMP CONDITION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE PHI_TAU(T,DATA)

      REAL :: T
      TYPE(LIST_DATA) :: DATA

      REAL :: XINF,YINF,THETA

      XINF  = DATA%X
      YINF  = DATA%Y
      THETA = DATA%THETA

      DATA%JUMP(3) = (COS(VK*XINF)*COS(VK*YINF)+SIN(VK*XINF)*SIN(VK*YINF))*VK*(COS(THETA)+SIN(THETA))*COS(T)

      RETURN

      END SUBROUTINE PHI_TAU

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! PSI_BAR --
      !    PSI_BAR = [BETA U_X] = BETA^+ U^+_X - BETA^- U^-_X
      ! ARGUMENTS:
      !    T      IN       CURRENT TIME
      !    DT     IN       TIME STEP
      !    DATA   IN/OUT   A LISTED INTERFACE POINT, DATA%JUMP(4) IS UPDATED
      ! NOTES:
      !    FISRT ORDER JUMP CONDITION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE PSI_BAR(T,DT,DATA)

      REAL :: T,DT
      TYPE(LIST_DATA) :: DATA

      REAL :: XINF,YINF,THETA,BETA1,BETA2

      XINF  = DATA%X
      YINF  = DATA%Y
      THETA = DATA%THETA
      BETA1 = DATA%BETA1
      BETA2 = DATA%BETA2

      DATA%JUMP(4) = (-VK*BETA2*SIN(VK*XINF)*SIN(VK*YINF)-VK*BETA1*COS(VK*XINF)*COS(VK*YINF))*COS(T+DT)
      DATA%JUMP(5) = (-VK*BETA2*SIN(VK*XINF)*SIN(VK*YINF)-VK*BETA1*COS(VK*XINF)*COS(VK*YINF))*COS(T)

      RETURN

      END SUBROUTINE PSI_BAR

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! PSI_HAT --
      !    PSI_HAT = [BETA U_Y] = BETA^+ U^+_Y - BEAT^- U^-_Y
      ! ARGUMENTS:
      !    T      IN       CURRENT TIME
      !    DT     IN       TIME STEP
      !    DATA   IN/OUT   A LISTED INTERFACE POINT, DATA%JUMP(4) IS UPDATED
      ! NOTES:
      !    FISRT ORDER JUMP CONDITION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE PSI_HAT(T,DT,DATA)

      REAL :: T,DT
      TYPE(LIST_DATA) :: DATA

      REAL :: XINF,YINF,THETA,BETA1,BETA2

      XINF  = DATA%X
      YINF  = DATA%Y
      THETA = DATA%THETA
      BETA1 = DATA%BETA1
      BETA2 = DATA%BETA2

      DATA%JUMP(4) = (VK*BETA2*COS(VK*XINF)*COS(VK*YINF)+VK*BETA1*SIN(VK*XINF)*SIN(VK*YINF))*COS(T+DT)
      DATA%JUMP(5) = (VK*BETA2*COS(VK*XINF)*COS(VK*YINF)+VK*BETA1*SIN(VK*XINF)*SIN(VK*YINF))*COS(T)

      RETURN

      END SUBROUTINE PSI_HAT

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                                                                        !
      !                              FOR TEST USING                            !
      !                                                                        !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! UTAU_MINUS --
      !    ANALYTICAL U^-_TAU AT A GIVEN INTERFACE PT
      ! ARGUMENTS:
      !    X      IN       X COORDINATE OF THE INTERFACE PT
      !    Y      IN       Y COORDINATE OF THE INTERFACE PT
      !    T      IN       CURRENT TIME
      !    THETA  IN       ANGLE FORMED BY THE NORMAL DIRECTION AND X-AXIS
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION UTAU_MINUS(X,Y,T,THETA)

         REAL :: X,Y,T,THETA,UTAU_MINUS

         UTAU_MINUS = -VK*COS(T)*( SIN(THETA)*COS(VK*X)*COS(VK*Y) + COS(THETA)*SIN(VK*X)*SIN(VK*Y) )

         RETURN

      END FUNCTION UTAU_MINUS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! UTAU_PLUS --
      !    ANALYTICAL U^+_TAU AT A GIVEN INTERFACE PT
      ! ARGUMENTS:
      !    X      IN       X COORDINATE OF THE INTERFACE PT
      !    Y      IN       Y COORDINATE OF THE INTERFACE PT
      !    T      IN       CURRENT TIME
      !    THETA  IN       ANGLE FORMED BY THE NORMAL DIRECTION AND X-AXIS
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION UTAU_PLUS(X,Y,T,THETA)

         REAL :: X,Y,T,THETA,UTAU_PLUS

         UTAU_PLUS = VK*COS(T)*( SIN(THETA)*SIN(VK*X)*SIN(VK*Y) + COS(THETA)*COS(VK*X)*COS(VK*Y) )

         RETURN

      END FUNCTION UTAU_PLUS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! UVAL --
      !    ANALYTICAL U VALUE AT A GIVEN PT
      ! ARGUMENTS:
      !    X      IN       X COORDINATE OF THE INTERFACE PT
      !    Y      IN       Y COORDINATE OF THE INTERFACE PT
      !    T      IN       CURRENT TIME
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION UVAL(X,Y,T)

         REAL :: X,Y,T,UVAL

         IF (SETS(X,Y) .GT. TOL_SETUP) THEN
            UVAL = COS(VK*X)*SIN(VK*Y)*COS(T)
         ELSE
            UVAL = SIN(VK*X)*COS(VK*Y)*COS(T)
         END IF

         RETURN

      END FUNCTION UVAL

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! UVAL_MINUS --
      !    ANALYTICAL U^- VALUE AT A GIVEN PT
      ! ARGUMENTS:
      !    X      IN       X COORDINATE OF THE INTERFACE PT
      !    Y      IN       Y COORDINATE OF THE INTERFACE PT
      !    T      IN       CURRENT TIME
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION UVAL_MINUS(X,Y,T)

         REAL :: X,Y,T,UVAL_MINUS

         UVAL_MINUS = SIN(VK*X)*COS(VK*Y)*COS(T)

         RETURN

      END FUNCTION UVAL_MINUS


      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! UVAL_PLUS --
      !    ANALYTICAL U^+ VALUE AT A GIVEN PT
      ! ARGUMENTS:
      !    X      IN       X COORDINATE OF THE INTERFACE PT
      !    Y      IN       Y COORDINATE OF THE INTERFACE PT
      !    T      IN       CURRENT TIME
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION UVAL_PLUS(X,Y,T)

         REAL :: X,Y,T,UVAL_PLUS

         UVAL_PLUS = COS(VK*X)*SIN(VK*Y)*COS(T)

         RETURN

      END FUNCTION UVAL_PLUS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! FPANALYTICAL --
      !    GENERATING ANALYTICAL SOLUTION FOR FICTITIOUS POINTS
      ! ARGUMENTS:
      !    U  OUT  MATRIX TO STORE EXACT SOLUTION ON ALL GRIDS
      !    T  IN   CURRENT TIME
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE FPANALYTICAL(U,T)

      REAL :: U(NY,NX),T

      REAL :: CT
      INTEGER :: IX,IY

      CT = COS(T)
      DO IY = 1,NY
         DO IX = 1,NX
            IF (SETS(XI(IX),YI(IY)) .GT. TOL_SETUP) THEN
               U(IY,IX) = SIN(VK*XI(IX))*COS(VK*YI(IY))*CT
            ELSE
               U(IY,IX) = COS(VK*XI(IX))*SIN(VK*YI(IY))*CT
            END IF
         END DO
      END DO

      RETURN

      END SUBROUTINE FPANALYTICAL
