      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2D HEAT EQUATION WITH DISC INTERFACE
      ! ============================================
      ! PDE: U_T = (BETA U_X)_X + (BETA U_Y)_Y + F
      !    BETA1=1     INSIDE
      !    BETA2=10    OUTSIDE
      !    U(X,Y)=COS(T)+EXP(X**2+Y**2) INSIDE
      !          =COS(T)+SIN(KX)COS(KY) OUTSIDE
      !    F=-SIN(T)-4*BETA1*EXP(X**2+Y**2)(X**2+Y**2+1)    INSIDE
      !     =-SIN(T)+2*BETA2*K**2*SIN(KX)*COS(KY)           OUTSIDE
      !
      ! GAMMA: R**2 = X**2 + Y**2
      !        R = (1/2)**2
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      !---------- TEST CASE IDENTIFIER ----------
      !INTEGER :: TESTCASE = 3

      !---------- INTERFACE DEFINITION ----------
      REAL :: RADIUS = 0.5D0

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

      SETS = SQRT(X**2+Y**2) - RADIUS + TOL_SETUP

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

      REAL :: CT,RAD
      INTEGER :: IX,IY

      CT = COS(T)
      DO IY = 1,NY
         DO IX = 1,NX
            RAD = SQRT( XI(IX)**2 + YI(IY)**2 )
            IF (RAD .LT. RADIUS) THEN      !INSIDE
               U(IY,IX) = CT + EXP( XI(IX)**2 + YI(IY)**2 )
            ELSE                           !OUTSIDE
               U(IY,IX) = CT + SIN( VK*XI(IX) ) * COS( VK*YI(IY) )
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
               BETA(IY,IX)  =  1.0D0 !ON-GRID BETA INSIDE
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
      ! SETSRC --
      !    ON-GRID SOURCE TERM AT TIME T
      ! ARGUMENTS:
      !    T  IN  CURRENT TIME
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SETSRC(T)

      REAL :: T

      REAL :: RAD
      INTEGER :: IX,IY

      DO IY = 1,NY
         DO IX = 1,NX
            RAD = SQRT( XI(IX)**2 + YI(IY)**2 )
            IF (RAD .LT. RADIUS) THEN  !INSIDE
               SRC(IY,IX) = -SIN(T) - 4.0D0 * BETA(IY,IX) * EXP(RAD**2) * (RAD**2 + 1.0D0)
            ELSE                       !OUTSIDE
               SRC(IY,IX) = -SIN(T) + 2.0D0 * BETA(IY,IX) * VK**2 * SIN( VK*XI(IX) ) * COS( VK*YI(IY) )
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

      CT = ( COS(T+DT) + COS(T) )/2.0D0

      DO IX = 1,NX,NX-1    !FIRST ORDER BOUNDARY CONDITION
         DO IY = 1,NY      !U*=(U^N+ U^{N+1})/2
            UHS(IY,IX) = CT + SIN( VK*XI(IX) ) * COS( VK*YI(IY) )
         END DO
      END DO

      DO IY = 1,NY,NY-1
         DO IX = 1,NX
            UHS(IY,IX) = CT + SIN( VK*XI(IX) ) * COS( VK*YI(IY) )
         END DO
      END DO

      END SUBROUTINE SETBC

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SETBC_ADIPR --
      !    SET BOUNDARY CONDITION FOR PEACEMAN RACHFORD ADI METHOD
      ! ARGUMENTS:
      !    T    IN   CURRENT TIME
      !    DT   IN   TIME STEP
      !    UHS  OUT  ON-GRID NUMERICAL SOLUTIONS
      ! NOTES:
      !    U* = (U^N + U^{N+1})/2 + BETA2 * DT/4 * (U_YY^N - U_YY^{N+1})
      !       = SIN(KX)*COS(KY) + ( COS(T) + COS(T+DT) )/2
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SETBC_ADIPR(T,DT,UHS)

      REAL :: T,DT,UHS(NY,NX)

      REAL :: CT
      INTEGER :: IX,IY

      CT = ( COS(T+DT) + COS(T) )/2.0D0

      DO IX = 1,NX,NX-1    !FIRST ORDER BOUNDARY CONDITION
         DO IY = 1,NY      !U*=(U^N+ U^{N+1})/2
            UHS(IY,IX) = SIN( VK*XI(IX) ) * COS( VK*YI(IY) ) + CT
         END DO
      END DO

      DO IY = 1,NY,NY-1
         DO IX = 1,NX
            UHS(IY,IX) = SIN( VK*XI(IX) ) * COS( VK*YI(IY) ) + CT
         END DO
      END DO

      END SUBROUTINE SETBC_ADIPR

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                                                                        !
      !                   SETUPS FOR INTERFACE GAMMA IN MIB                    !
      !                                                                        !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

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

      REAL :: TMP

      TMP = SQRT( RADIUS**2 - Y**2 )
      IF (ABS( TMP - ABS(X1) ) .LT. 1.0D-9) THEN
         X = X1 + (X2-X1) * 1.0D-9
         WRITE(*,*) "WARNING: TRANSLATION IN X",X,Y
      ELSEIF ((TMP.GT.X1) .AND. (TMP.LT.X2)) THEN
         X =  TMP
      ELSEIF ((-TMP.GT.X2) .AND. (-TMP.LT.X1)) THEN
         X = -TMP
      ELSE
         WRITE(*,*) "ERROR in GAMMAX"
         STOP
      END IF

      VARPHI = GETANGLE(X,Y,RADIUS)

      RETURN

      END SUBROUTINE GAMMAX

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

      REAL :: TMP

      TMP = SQRT( RADIUS**2 - X**2 )
      IF ( ABS( TMP-ABS(Y1) ) .LT. 1.0D-9 ) THEN
         Y = Y1 + (Y2 - Y1) * 1.0D-9
         WRITE(*,*) "WARNING: TRANSLATION IN Y",Y,X
      ELSE IF ( ( TMP.GT.Y1) .AND. ( TMP.LT.Y2) ) THEN
         Y = TMP
      ELSE IF ( (-TMP.GT.Y2) .AND. (-TMP.LT.Y1) ) THEN
         Y = -TMP
      ELSE
         WRITE(*,*) "ERROR in GAMMAY"
         STOP
      END IF

      VARPHI = GETANGLE(X,Y,RADIUS)

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

      GETNORMAL = VARPHI

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
      TYPE (LIST_DATA) :: DATA

      REAL :: XINF,YINF

      XINF = DATA%X
      YINF = DATA%Y

      DATA%JUMP(1) = SIN( VK*XINF ) * COS(VK*YINF) - EXP(RADIUS**2)

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
      TYPE (LIST_DATA) :: DATA

      REAL :: XINF,YINF,THETA,BETA1,BETA2

      XINF  = DATA%X
      YINF  = DATA%Y
      THETA = DATA%THETA
      BETA1 = DATA%BETA1
      BETA2 = DATA%BETA2

      DATA%JUMP(2) = BETA2 * VK * COS(THETA) * COS(VK*XINF) * COS(VK*YINF) - &
                     BETA2 * VK * SIN(THETA) * SIN(VK*XINF) * SIN(VK*YINF) - &
                     BETA1 * EXP( RADIUS**2 )

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

      DATA%JUMP(3) = -VK * ( SIN(THETA) * COS(VK*XINF) * COS(VK*YINF) - COS(THETA) * SIN(VK*XINF) * SIN(VK*YINF) )
      
      RETURN

      END SUBROUTINE PHI_TAU

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! PSI_BAR --
      !    PSI_BAR = [BETA U_X] = BETA^+ U^+_X - BETA^- U^-_X
      ! ARGUMENTS:
      !    T      IN       CURRENT TIME
      !    DATA   IN/OUT   A LISTED INTERFACE POINT, DATA%JUMP(4) IS UPDATED
      ! NOTES:
      !    FISRT ORDER JUMP CONDITION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE PSI_BAR(T,DT,DATA)

      REAL :: T,DT
      TYPE(LIST_DATA) :: DATA

      REAL :: XINF,YINF,BETA1,BETA2

      XINF  = DATA%X
      YINF  = DATA%Y
      BETA1 = DATA%BETA1
      BETA2 = DATA%BETA2

      DATA%JUMP(4) = BETA2 * VK * COS(VK*XINF) * COS(VK*YINF) - BETA1 * 2.0D0 * XINF * EXP(RADIUS**2)
      DATA%JUMP(5) = BETA2 * VK * COS(VK*XINF) * COS(VK*YINF) - BETA1 * 2.0D0 * XINF * EXP(RADIUS**2)

      RETURN

      END SUBROUTINE PSI_BAR

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! PSI_HAT --
      !    PSI_HAT = [BETA U_Y] = BETA^+ U^+_Y - BEAT^- U^-_Y
      ! ARGUMENTS:
      !    T      IN       CURRENT TIME
      !    DATA   IN/OUT   A LISTED INTERFACE POINT, DATA%JUMP(4) IS UPDATED
      ! NOTES:
      !    FISRT ORDER JUMP CONDITION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE PSI_HAT(T,DT,DATA)

      REAL :: T,DT
      TYPE(LIST_DATA) :: DATA

      REAL :: XINF,YINF,BETA1,BETA2

      XINF  = DATA%X
      YINF  = DATA%Y
      BETA1 = DATA%BETA1
      BETA2 = DATA%BETA2

      DATA%JUMP(4) = -BETA2 * VK * SIN(VK*XINF) * SIN(VK*YINF) - BETA1 * 2.0D0 * YINF * EXP(RADIUS**2)
      DATA%JUMP(5) = -BETA2 * VK * SIN(VK*XINF) * SIN(VK*YINF) - BETA1 * 2.0D0 * YINF * EXP(RADIUS**2)

      RETURN

      END SUBROUTINE PSI_HAT

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                                                                        !
      !                              FOR TEST USED                             !
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

         UTAU_MINUS = 2.0D0 * EXP( X**2 + Y**2 ) * ( -X*SIN(THETA) + Y*COS(THETA) )

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

         UTAU_PLUS = -VK * ( SIN(THETA)*COS(VK*X)*COS(VK*Y) + COS(THETA)*SIN(VK*X)*SIN(VK*Y) )

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
            UVAL = COS(T) + SIN(VK*X) * COS(VK*Y)
         ELSE
            UVAL = COS(T) + EXP( X**2 + Y**2 )
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

         UVAL_MINUS = COS(T) + EXP( X**2 + Y**2 )

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

         UVAL_PLUS = COS(T) + SIN(VK*X) * COS(VK*Y)

         RETURN

      END FUNCTION UVAL_PLUS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! FPANALYTICAL --
      !
      ! ARGUMENTS:
      !
      !
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE FPANALYTICAL(U,T)

      REAL :: U(NY,NX),T

      REAL :: CT,RAD
      INTEGER :: IX,IY

      CT = COS(T)
      DO IY = 1,NY
         DO IX = 1,NX
            RAD = SQRT( XI(IX)**2 + YI(IY)**2 )
            IF (RAD .LT. RADIUS) THEN      !INSIDE
               U(IY,IX) = CT + SIN( VK*XI(IX) ) * COS( VK*YI(IY) )
            ELSE                           !OUTSIDE
               U(IY,IX) = CT + EXP( XI(IX)**2 + YI(IY)**2 )
            END IF
         END DO
      END DO

      RETURN

      END SUBROUTINE FPANALYTICAL
