      !-----------------------------------------------------------------!
      !           MODULE FOR INTERFACE POINTS                           !
      !-----------------------------------------------------------------!
      MODULE MOD_INTERFACE

      TYPE INTERFACE
         !-----------------------------------------------------------------!
         !           BASIC INFORMATION ABOUT THE INTERFACE PT              !
         !-----------------------------------------------------------------!

         !INDICATOR OF IX (= 0) OR IY (= 1) AXIS ON WHICH THE INTERFACE POINT IS
         INTEGER :: AXTP

         !INDEX OF THE IX- OR IY- AXIS ON WHICH THIS INTERFACE POINT IS
         INTEGER :: AXID

         !IDENTITY OF AN INTERFACE POINT. IT POSSESES 3 LAYERS OF MEANINGS:
         !   1. ABSOULTE VALUE GIVES THE INDEX OF THIS POINT IN THE LIST
         !   2. +/- INDICATES IRREGULAR/CORNER POINT
         !   3. ODD/EVEN INDICATES IT IS FROM OMEGA^+ TO OMEGA^- OR VICE VERSA
         INTEGER :: ID

         !LOCATION OF THE INTERFACE POINT (INDEX OF THE LEFT GRID OF ONE PAIR OF NEIGHBORING NODES)
         INTEGER :: LGRD

         !COORDINATES OF THE INTERFACE POINT
         REAL :: X,Y

         !ANGLE FORMED BY THE NORMAL DIRECTION AND X-AXIS
         REAL :: THETA

         !DIFFUSION COEFFICIENT ON THE INTERFACE
         REAL :: BETA1,BETA2 !BEAT1 = BETA^-,BETA2 = BETA^+ AT THE INTERFACE POINT

         !0 < GAMMA < DX, DISTANCE BETWEEN THE INTERFACE POINT AND ITS LEFT END
         REAL :: GAMMA

         !TYPE OF THE INTERFACE POINT
         !                                                 NL-1      NL         NR         NR+1
         !   -1: INTERFACE IS CLOSE TO THE LEFT GRID       o----------|--x-------o----------|
         !                                                 F             I       F
         !
         !    0: INTERFACE STAYS IN BETWEEN THE TWO GRIDS  |----------o-----x----o----------|
         !                                                            F     I    F
         !
         !    1: INTERFACE IS CLOSE TO THE RIGHT GRID      |----------o-------x--|----------o
         !                                                            F       I             F
         !
         !(FOR CORNER POINT ONLY)                          NL-1       NL        NM         NR         NR+1
         !   -2: THE LEFT INTERFFACE POINT OF A CORNER     o----------o--x-------o-----x----o----------|
         !       IS CLOSER TO THE NODE                     F          F  IL      F     IR   F
         !
         !                                                 o----------o--x-------o--------x-|----------o
         !                                                 F          F  IL      F       IR            F
         !
         !    2: THE RIGHT INTERFFACE POINT OF A CORNER    |----------o------x---o-------x--o----------o
         !       IS CLOSER TO THE NODE                                F      IL  F       IR F          F
         !
         !                                                 o----------|-x--------o----x-----o----------o
         !                                                 F            IL       F    IR    F          F
         !
         !   "CLOSENESS" IS DETERMINED BY < 0.1*DX RIGHT NOW, CAN BE CHANGED...
         INTEGER :: ITYPE

         !MIB FP REPRESENTATION COEFFS FOR THE LEFT AND RIGHT NODES OF AN INTERFACE POINT
         REAL :: WIJ(2,6)  !2 - TWO FICTITOUS PTS,
                           !6 - WEIGHTS FOR ( IX-1,IX,IX+1,IX+2,
                           !                  [U]_IP,[BEAT U_X]_IP (OR [BETA U_Y]_IP) )
         REAL :: WIJ2(4,9) !4 - FOUR FICTIOUS PTS,
                           !9 - WEIGHTS FOR ( IX-1,IX,IX+1,IX+2,IX+3,
                           !                  [U]_IPL,[BETA U_X]_IPL (OR [BETA U_Y]_IPL),
                           !                  [U]_IPR,[BETA U_X]_IPR (OR [BETA U_Y]_IPR) )

         !-----------------------------------------------------------------!
         !     INFORMATION FOR APPROXMATING U_ATU AT THE INTERFACE PT      !
         !-----------------------------------------------------------------!

         !TWO AUXILIARY POINTS, CALLED AUXL AND AUXR ARE NEEDED TO INTERPLOATE EITHER
         !U_TAU^+ OR U_TAU^- AT EACH INTERFACE POINTS. DEPENDING ON WHERE THE TWO AUXILIARY
         !PTS ARE, THE FOLLOWING CASES MAY OCCUR:
         !CASE 1. BOTH AUXILIARY PTS ARE IN OMEGA^+ SO THAT U_TAU^+ IS INTERPOLATED;
         !CASE 2. BOTH AUXILIARY PTS ARE IN OMEGA^- SO THAT U_TAU^- IS INTERPOLATED;
         !CASE 3. AUXL IS IN OMEGA^- AND AUXR IS IN OMEGA^+;
         !CASE 4. AUXL IS IN OMEGA^+ AND AUXR IS IN OMEGA^-;
         !FOR CASE 3 AND 4, THE FOLLOWING STEPS TO BE CARRIED OUT:
         !STEP 1. INTERPLOATING ONE U VALUE (U^+ OR U^-) AT THE AUXILIARY PTS BASED ON
         !        THEIR LOCATION;
         !STEP 2. SOLVE THE FOLLOWING EQUATIONS FOR THE OTHER VALUE ON THESE TWO
         !        AUXILIARY PTS
         !        (U^+_AUXL + U^+_AUXR)/2  - (U^-_AUXL + U^-_AUXR)/2  = [U]_IP     = PHI_IP
         !        (U^+_AUXR - U^+_AUXL)/2H - (U^-_AUXR - U^-_AUXL)/2H = [U_TAU]_IP = (PHI_TAU)_IP
         !AFTER OBTAINING U^+ AND U^- AT BOTH AUXILIARY PTS, U_TAU^+ OR U_TAU^- CAN BE
         !INTERPOLATED FREELY.

         !TANGENTIAL TAU DIRECTION, EITHER +1 IF POSITVE DIRECTION OR -1 IF NEGATIVE DIRECTION
         INTEGER :: ITAU

         !---------- THE LEFT/LOWER AUXILIARY POINT

         !INDICATOR OF IX (= 0) OR IY (= 1) AXIS ON WHICH THE LEFT/LOWER AUXILIARY POINT IS
         INTEGER :: AUXL_AXTP

         !INDEX OF THE IX- OR IY- GRID LINE ON WHICH THE LEFT/UPPER AUXILIARY POINT IS
         INTEGER :: AUXL_AXID

         !COORDINATES OF THE LEFT AUXILIARY POINT
         REAL :: AUXL_X,AUXL_Y

         !TYPE OF THE LEFT AUXILIARY POINT FOR INTERPOLATING U_TAU,
         != 1 IF U^+ IS INTERPOLATED, = -1 IF U^- IS INTERPOLATED
         INTEGER :: AUXL

         !DISTANCE OF THE LEFT AUXILIARY PT TO THE INTERFACE PT
         REAL :: DAUXL

         !INDICES OF THREE NODES USED FOR INTERPOLATING U^+ OR U^- AT THE LEFT AUXILIARY PT
         INTEGER :: IAUXL(3)

         !WEIGHTS FOR INTERPOLATING U^+ OR U^- AT THE LEFT AUXILIARY PT
         REAL :: WAUXL(3)

         !---------- THE RIGHT/UPPER AUXILIARY POINT

         !INDICATOR OF IX (= 0) OR IY (= 1) AXIS ON WHICH THE RIGHT AUXILIARY POINT IS
         INTEGER :: AUXR_AXTP

         !INDEX OF THE IX- OR IY- AXIS ON WHICH THE RIGHT AUXILIARY POINT IS
         INTEGER :: AUXR_AXID

         !COORDINATES OF THE RIGHT AUXILIARY POINT
         REAL :: AUXR_X,AUXR_Y

         !TYPE OF THE RIGHT AUXILIARY POINT FOR INTERPOLATING U_TAU
         != 1 IF U^+ IS INTERPOLATED, = -1 IF U^- IS INTERPOLATED
         INTEGER :: AUXR

         !DISTANCE OF THE RIGHT AUXILIARY PT TO THE INTERFACE PT
         REAL :: DAUXR

         !INDICES OF 3 NODES USED FOR INTERPOLATING U AT THE RIGHT AUXILIARY POINT
         INTEGER :: IAUXR(3)

         !WEIGHTS FOR INTERPOLATING U^+ OR U^- AT THE RIGHT AUXILIARY POINT
         REAL :: WAUXR(3)

         !---------- U VALUES TO BE EXAMINED

         !INTERFACE JUMP CONDITIONS
         REAL :: JUMP(6) !JUMP(1) = [U]                             AT TIME STEP T^{N+1}
                         !JUMP(2) = [BETA U_N]                      AT TIME STEP T^{N+1}
                         !JUMP(3) = [U_TAU]                         AT TIME STEP T^{N+1}
                         !JUMP(4) = ANA. [BETA U_X] (OR [BETA U_Y]) AT TIME STEP T^{N+1}
                         !JUMP(5) = ANA. [BETA U_X] (OR [BETA U_Y]) AT TIME STEP T^{N}
                         !JUMP(6) = NUM. [BETA U_X] (OR [BETA U_Y]) AT TIME STEP T^{N}

         ! U_TAU AT THE INTERFACE PT
         !   NOTE: IN CURRENT TIME, NOT IN THE FUTURE TIME!
         REAL :: UTAU(4) ! (ANA. U^+_TAU, NUM. U^+_TAU, ANA. U^-_TAU, NUM. U^-_TAU)

         !ANA. U AND NUM. U AT THE LEFT AUXILIARY PT
         REAL :: UAUXL(4) ! UAUXL(1) = INDICATOR, 1 IF IN OMEGA^+, -1 IF IN OMEGA^-
                          ! UAUXL(2) = INDICATOR, 1 IF U^+ IS CALCULATED, -1 IF U^- IS CALCULATED
                          ! UAUXL(3) = ANA. U VALUE AT LEFT AUXILIARY PT
                          ! UAUXL(4) = NUM. U VALUE AT LEFT AUXILIARY PT

         !ANA. U AND NUM. U AT THE RIGHT AUXILIARY PT
         REAL :: UAUXR(4) ! UAUXR(1) = INDICATOR, 1 IF IN OMEGA^+, -1 IF IN OMEGA^-
                          ! UAUXR(2) = INDICATOR, 1 IF U^+ IS CALCULATED, -1 IF U^- IS CALCULATED
                          ! UAUXR(3) = ANA. U VALUE AT RIGHT AUXILIARY PT
                          ! UAUXR(4) = NUM. U VALUE AT RIGHT AUXILIARY PT

         !--------------------------------------------------------------------------

         !FOR MIB PART TEST USED
         REAL :: ERR(2)
         REAL :: ERR2(4)

      END TYPE INTERFACE

      END MODULE MOD_INTERFACE

      !-----------------------------------------------------------------!
      !           MODULE FOR LINKED LIST OF INTERFACE POINTS            !
      !-----------------------------------------------------------------!
      MODULE MOD_LISTS

      USE MOD_INTERFACE,LIST_DATA => INTERFACE

      INCLUDE "data_linkedlist.f90"

      END MODULE MOD_LISTS

      !-----------------------------------------------------------------!
      !           MODULE FOR GLOBAL VARIABLES                           !
      !-----------------------------------------------------------------!
      MODULE MOD_DATA

      USE MOD_LISTS

      IMPLICIT NONE

      !---------- CONSTANT PARAMETERS ----------
      INTEGER,PARAMETER :: NX = 71,NY = NX
      INTEGER,PARAMETER :: NPRINT = 10
      REAL,   PARAMETER :: CD =1.99D0 !PI*0.32D0
      REAL,   PARAMETER :: TSTART = 0.0D0,TFINAL = 1.0D1,TSTEP = 1.0D-3
      REAL,   PARAMETER :: PI = 3.14159265358979323846264338328D0

      !---------- FINITE DIFFERENCE FORMULATION ----------
      INTEGER :: INODE(NY,NX)           !NODE INDEX INDICATING INSIDE OR OUTSIDE GAMMA

      !---------- SOURCE TERM ----------
      REAL :: SRC(NY,NX)                !ON-GRID SOURCE

      !---------- DIFFUSION AND CONVECTION COEFFICIENTS ----------
      REAL :: BETA(NY,NX)               !ON-GRID BETA
      REAL :: BETAX(NY,NX)              !(dBETA/dX)/BETA
      REAL :: BETAY(NY,NX)              !(dBETA/dY)/BETA

      !---------- SPATIAL DOMAIN DISCRETIZATION ----------
      REAL :: XI(NX),YI(NX)             !X AND Y CARTESIAN NODES
      REAL :: XL,YL,DX                  !XLEFT, YLEFT AND DX

      !---------- LINKED LISTS OF INTERFACE POINTS ----------
      !GET ARRAYS OF POINTERS BY DEFINING A NEW TYPE
      TYPE LIST_ARRAY
         TYPE(LINKED_LIST),POINTER :: HEAD
      END TYPE LIST_ARRAY

      TYPE(LIST_ARRAY),DIMENSION(NY) :: IFPY !LIST HEAD FOR EACH Y-AXIS IN X- DIRECTION
      TYPE(LIST_ARRAY),DIMENSION(NX) :: IFPX !LIST HEAD FOR EACH X-AXIS IN Y- DIRECTION

      !---------- TOLERANCES ----------
      REAL :: TOL_SETUP  = 1.0E-12 !TOLERANCE USED IN SETS
      REAL :: TOL_RTSAFE = 1.0E-14 !TOLERANCE USED IN RTSAFE
      REAL :: TOL_ITYPE !LEFT AND RIGHT TOLERANCES TO DETERMINE DATA%ITYPE

      !---------- DEBUGGING FLAGS AND RELATED ----------
      LOGICAL :: DEBUGFLAG = .FALSE.    !debug flag
      LOGICAL :: WARNINGS  = .FALSE.    !FLAG TO COMPRESS WARNINGS GENERATED BY COMPILER

      INTEGER :: NIPXS    = 0 !# OF INTERFACE PTS ON ALL IX GRID LINES
      INTEGER :: NIPYS    = 0 !# OF INTERFACE PTS ON ALL IY GRID LINES
      INTEGER :: NCORNERS = 0 !# OF CORNER PTS
      INTEGER :: NALTAUXS = 0 !# OF INTERFACE PTS THAT U^+ AND U^- ARE INTERPLOATED AT TWO AUXILIARY PTS

      !------ ATTENTION : INCLUDE ONE OF THE FOLLOWING EXAMPLE!!!
!-----------------------------!
! EXAMPLE 1 IN ZHAO JSC 2014  !
!-----------------------------!
#if   example == 1
      INTEGER :: TESTCASE = 1
      INCLUDE "data_testcase_disc1.f90"

!-----------------------------!
! EXAMPLE 2 IN ZHAO JSC 2014  !
!-----------------------------!
#elif   example == 2
      INTEGER :: TESTCASE = 2
      INCLUDE "data_testcase_disc2.f90"

!-----------------------------!
! EXAMPLE 3 IN ZHAO JSC 2014  !
!-----------------------------!
#elif   example == 3
      INTEGER :: TESTCASE = 3
      INCLUDE "data_testcase_disc3.f90"

!-----------------------------!
! EXAMPLE 4 IN ZHAO JSC 2014  !
!-----------------------------!
#elif   example == 4
      INTEGER :: TESTCASE = 4
      INCLUDE "data_testcase_star45.f90"

!-----------------------------!
! EXAMPLE 5 IN ZHAO JSC 2014  !
!-----------------------------!
#elif example == 5
      INTEGER :: TESTCASE = 5
      INCLUDE "data_testcase_star45.f90"
      INCLUDE "data_rtsafe.f90"

!-----------------------------!
! NEW EXAMPLE 6               !
!-----------------------------!
#elif example == 6
      INTEGER :: TESTCASE = 6
      INCLUDE "data_testcase_star6.f90"
      INCLUDE "data_rtsafe.f90"

!-----------------------------!
! NEW EXAMPLE 7               !
!-----------------------------!
#elif example == 7
      INTEGER :: TESTCASE = 7
      INCLUDE "data_testcase_star7.f90"
      INCLUDE "data_rtsafe.f90"

#endif

      !------ SUBROUTINES USED BY ALL EXAMPLES
      INCLUDE "data_test.f90"

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! GETANGLE --
      !    FIND VARPHI BY USING ACOS(VARPHI)=X/R
      ! ARGUMENTS:
      !    XO        IN   X-COORDINATE OF THE POINT
      !    YO        IN   Y-COORDINATE OF THE POINT
      !    RADI      IN   RADIUS/DISTANCE OF THE POINT FROM ORIGIN
      !    GETANGLE  OUT  NORMAL ANGLE OF THE GIVEN POINT
      ! NOTES:
      !    0<=ACOS<=PI, ALL ANGLE GET THIS SUBROUTINE IS BETWEEN [0,2*PI]
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      FUNCTION GETANGLE(XO,YO,RADI)

      REAL :: XO,YO,RADI,GETANGLE

      IF (ABS(RADI) .LT. 1.0D-14) THEN
         GETANGLE = 0.0D0
      ELSE IF (YO .EQ. 0) THEN
         IF (XO .GT. 0) THEN
            GETANGLE = 0.0D0
         ELSE
            GETANGLE = ACOS(-1.0D0)
         END IF
      ELSE IF (YO .GT. 0) THEN
         GETANGLE = ACOS(XO/RADI)
      ELSE
         GETANGLE = 2.0D0 * PI - ACOS(XO/RADI)
      END IF

      RETURN

      END FUNCTION GETANGLE

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! OUTERROR --
      !    ERROR ANALYSIS
      ! ARGUMENTS:
      !    T   IN  CURRENT TIME
      !    U   IN  EXACT SOLUTION
      !    UH  IN  NUMERIAL SOLUTION
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE OUTERROR(T,U,UH)

      REAL :: T,U(NY,NX),UH(NY,NX)

      REAL :: VM,EL2,TEMP
      INTEGER :: IX,IY,IT

      VM  = 0.D0
      EL2 = 0.D0

      DO IY = 1,NY
         DO IX = 1,NX
            TEMP = ABS( U(IY,IX) - UH(IY,IX) )
            IF (TEMP .GT. VM) VM = TEMP
            EL2 = EL2 + TEMP * TEMP
         END DO
      END DO

      EL2 = SQRT( EL2/(NY * 1.0D0 * NX) )

      IT = 10 !CHUAN: ADD IT=10 TO GET INTO THE FOLLOWING IF STATEMENT

      IF (MOD(IT,10) .EQ. 0) THEN
         WRITE (*, 888) T,VM,EL2
      END IF

888   FORMAT ("T=",F12.4,5 X,"LMAX=",E12.6,5 X,"L2=",E12.6)

      RETURN

      END SUBROUTINE OUTERROR

      END MODULE MOD_DATA

