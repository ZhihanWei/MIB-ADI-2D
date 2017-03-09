      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                                                                        !
      !            API SUBROUTINES PROVIDED BY THE MODULE MOD_MIB              !
      !                                                                        !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SETUP --
      !    INITIAL SETUP (SET VALUES ON ALL NODES)
      ! ARGUMENTS:
      !
      ! NOTES:
      !    SETUP SPATIAL DISCRETIZATION AND PDE COEFFICIENTS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SETUP

      USE MOD_DATA

      INTEGER :: IX,IY

      !----- SETUP XI AND YI
      DO IX=1,NX
         XI(IX)=XL+(IX-1.0D0)*DX
         YI(IX)=YL+(IX-1.0D0)*DX
      END DO

      !----- SETUP BETA VALUES ON ALL NODES
      CALL SETBETA

      !----- DETERMINE THE GRIDS INSIDE OR OUTSIDE THE INTERFACE
      DO IY = 1,NY
         DO IX = 1,NX
            IF (SETS(XI(IX),YI(IY)) .GT. TOL_SETUP) THEN
               INODE(IY,IX) =  1 !OUTSIDE
            ELSE
               INODE(IY,IX) = -1 !INSIDE (INCLUDE ON THE INTERFACE)
            END IF
         END DO
      END DO

      !----- CHECK WHETHER THE INTERFACE IS TOO CLOSE TO THE BOUNDARY
      DO IY = 1,NY
         IF (INODE(IY,1)*INODE(IY,2) .LT. 0 .OR. INODE(IY,NX-1)*INODE(IY,NX) .LT. 0) THEN
            WRITE(*,*) "SETUP: INTERFACE IS TOO CLOSE TO THE BOUNDARY. SPATIAL REFINEMENT IS REQUIRED..."
            STOP
         END IF
      END DO

      DO IX = 1,NX
         IF (INODE(1,IX)*INODE(2,IX) .LT. 0 .OR. INODE(NY-1,IX)*INODE(NY,IX) .LT. 0) THEN
            WRITE(*,*) "SETUP: INTERFACE IS TOO CLOSE TO THE BOUNDARY. SPATIAL REFINEMENT IS REQUIRED..."
            STOP
         END IF
      END DO

      RETURN

      END SUBROUTINE SETUP

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! FPSETUP --
      !    FICITIOUS VALUE REPRESENTATION
      ! ARGUMENTS:
      !
      ! NOTES:
      !    MUST BE CALLED ONCE AFTER CALLING SUBROUTINE SETUP
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE FPSETUP

      USE MOD_DATA

      INTEGER :: IX,IY,IFP
      TYPE(LIST_DATA) :: DATA !DATA OF AN IRREGULAR POINT
      TYPE(LIST_DATA) :: DATAL,DATAR !DATA OF A CORNER POINT
      TYPE(LINKED_LIST), POINTER  :: ELEM

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! STEP I: CREATE LINKED LIST FOR INTERFACE POINTS ON EACH IX- OR IY- AXIS
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !---------- FOR EACH IY-AXIS IN X DIRECTION ----------!
      DO IY = 2,NY-2
         IFP  =  0 !TRACK INTERFACE POINT
         IX   =  2
         DO WHILE (IX .LT. (NX-1))
            IF ((INODE(IY,IX)*INODE(IY,IX+1)) .LT. 0) THEN !CROSS INTERFACE
               IF ((INODE(IY,IX+1)*INODE(IY,IX+2)) .GT. 0) THEN !ONCE - IRREGULAR POINT
                  IFP = IFP+1 !ONE IRREGULAR POINT
                  CALL NEWIPY(IY,IX,IFP,DATA)
                  CALL LIST_APPEND(IFPY(IY)%HEAD,DATA) !APPEND A NEW ELEMENT TO THE END OF THE LIST
                  IX    = IX + 2
                  NIPYS = NIPYS + 1
               ELSE !TWICE - CORNER POINT
                  !----- CHECK WHETHER TWO SUCCESSIVE NODE ON ONE SIDE OF THE INTERFACE ARE IN THE SAME SUBDOMAIN
                  IF ((INODE(IY,IX-1)*INODE(IY,IX)) .LT. 0 .OR. (INODE(IY,IX+2)*INODE(IY,IX+3)) .LT. 0) THEN
                     WRITE(*,*) "FPSETUP: TWO INTERFACE POINTS ON AXIS IY = ",IY," IS TOO CLOSE. &
                                &SPATIAL REFINEMENT IS REQUIRED..."
                     STOP
                  END IF

                  IFP = IFP+1 !LEFT INTERFACE POINT
                  CALL NEWIPY2(IY,IX,IFP,DATAL,DATAR)
                  CALL LIST_APPEND(IFPY(IY)%HEAD,DATAL) !APPEND THE LEFT  INTERFACE POINT TO THE END OF THE LIST

                  IFP = IFP+1 !RIGHT INTERFACE POINT
                  CALL LIST_APPEND(IFPY(IY)%HEAD,DATAR) !APPEND THE RIGHT INTERFACE POINT TO THE END OF THE LIST

                  IX       = IX+3
                  NCORNERS = NCORNERS + 1
                  NIPYS    = NIPYS + 2
               END IF
            ELSE !NO INTERFACE
               IX = IX+1
            END IF
         END DO
         IF (MOD(IFP,2) .EQ. 1) THEN      ! IT'S AN OPEN CURVE ON THE DOMAIN 
            WRITE(*,*) "ERROR IN Y-DIRECTION: CHANGE THE SIZE OF DOMAIN! IT'S AN OPEN CURVE ON THE DOMAIN!"
            STOP
         END IF          
      END DO

      !---------- FOR EACH IX-AXIS IN Y DIRECTION ----------!
      DO IX = 2,NX-2
         IFP  =  0 !TRACK INTERFACE POINT
         IY   =  2
         DO WHILE (IY .LT. (NY-1))
            IF ((INODE(IY,IX)*INODE(IY+1,IX)) .LT. 0) THEN !CROSS INTERFACE
               IF ((INODE(IY+1,IX)*INODE(IY+2,IX)) .GT. 0) THEN !ONCE - IRREGULAR POINT
                  IFP = IFP+1 !ONE IRREGULAR POINT
                  CALL NEWIPX(IX,IY,IFP,DATA)
                  CALL LIST_APPEND(IFPX(IX)%HEAD,DATA) !APPEND A NEW ELEMENT TO THE END OF THE LIST
                  IY    = IY + 2
                  NIPXS = NIPXS + 1
               ELSE !TWICE - CORNER POINT
                  !----- CHECK WHETHER TWO SUCCESSIVE NODE ON ONE SIDE OF THE INTERFACE ARE IN THE SAME SUBDOMAIN
                  IF ((INODE(IY-1,IX)*INODE(IY,IX)) .LT. 0 .OR. (INODE(IY+2,IX)*INODE(IY+3,IX)) .LT. 0) THEN
                     WRITE(*,*) "FPSETUP: TWO INTERFACE POINTS ON AXIS IX = ",IX," IS TOO CLOSE. &
                                &SPATIAL REFINEMENT IS REQUIRED..."
                     STOP
                  END IF

                  IFP = IFP+1 !LOWER INTERFACE POINT
                  CALL NEWIPX2(IX,IY,IFP,DATAL,DATAR)
                  CALL LIST_APPEND(IFPX(IX)%HEAD,DATAL) !APPEND THE LOWER INTERFACE POINT TO THE END OF THE LIST

                  IFP = IFP+1 !UPPER INTERFACE POINT
                  CALL LIST_APPEND(IFPX(IX)%HEAD,DATAR) !APPEND THE UPPER INTERFACE POINT TO THE END OF THE LIST

                  IY       = IY + 3
                  NCORNERS = NCORNERS + 1
                  NIPXS    = NIPXS + 2
               END IF
            ELSE !NO INTERFACE
               IY = IY+1
            END IF 
         END DO
         IF (MOD(IFP,2) .EQ. 1) THEN      ! IT'S AN OPEN CURVE ON THE DOMAIN 
            WRITE(*,*) "ERROR IN X-DIRECTION: CHANGE THE SIZE OF DOMAIN! IT'S AN OPEN CURVE ON THE DOMAIN!"
            STOP
         END IF
      END DO

      ! NOTICE: APPROXIMATING U VALUES AT THE AUXILIARY POINTS REQUIRES KNOWING THE TYPE OF ALL INTERFACE POINTS,
      ! I.E., -1, 0, OR 1, SO THAT ANOTHER ITERATION IS NEEDED NEXT TO SET WEIGHTS FOR APPROXIMATING JUMP
      ! CONDITIONS.

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! STEP II: SET JUMP CONDITIONS FOR EACH LISTED INTERFACE POINT
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !---------- FOR EACH IY-AXIS IN X DIRECTION ----------!
      DO IY = 2,NY-2

         IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN !INTERFACE CROSSES THIS AXIS
            ELEM => IFPY(IY)%HEAD !THE 1ST INTERFACE POINT
            DO WHILE ( ASSOCIATED(ELEM) )
               !----- AN IRREGULAR INTERFACE POINT
               IF (ELEM%DATA%ID .GT. 0) THEN
                  DATA = LIST_GET_DATA(ELEM) !GET THE DATA FOR ONE LISTED INTERFACE POINT
                  CALL SETUTAU(DATA)
                  CALL MIB1D(DATA)
                  CALL LIST_PUT_DATA(ELEM,DATA) !PUT BACK THE UPDATED DATA FOR ONE INTERFACE POINT
                  ELEM => ELEM%NEXT !THE NEXT INTERFACE POINT
               !----- A CORNER POINT
               ELSE
                  DATAL = LIST_GET_DATA(ELEM) !GET THE DATA OF THE LEFT INTERFACE POINT
                  CALL SETUTAU(DATAL)
                  DATAR = LIST_GET_DATA(ELEM%NEXT) !GET THE DATA OF THE LEFT INTERFACE POINT
                  CALL SETUTAU(DATAR)
                  CALL CMIB1D(DATAL,DATAR)
                  CALL LIST_PUT_DATA(ELEM,DATAL)      !PUT BACK THE UPDATED LEFT  INTERFACE POINT
                  CALL LIST_PUT_DATA(ELEM%NEXT,DATAR) !PUT BACK THE UPDATED RIGHT INTERFACE POINT
                  ELEM => ELEM%NEXT%NEXT !SKIP THE RIGHT INTERFACE POINT OF A CORNER POINT
               END IF
            END DO
         END IF 

      END DO

      !---------- FOR EACH IX-AXIS IN Y DIRECTION ----------!
      DO IX = 2,NX-2

         IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN !INTERFACE CROSSES THIS AXIS
            ELEM => IFPX(IX)%HEAD !THE 1ST INTERFACE POINT
            DO WHILE ( ASSOCIATED(ELEM) )
               !----- AN IRREGULAR INTERFACE POINT
               IF (ELEM%DATA%ID .GT. 0) THEN
                  DATA = LIST_GET_DATA(ELEM) !GET THE DATA FOR ONE LISTED INTERFACE POINT
                  CALL SETUTAU(DATA)
                  CALL MIB1D(DATA)
                  CALL LIST_PUT_DATA(ELEM,DATA) !PUT BACK THE UPDATED DATA FOR ONE INTERFACE POINT
                  ELEM => ELEM%NEXT !THE NEXT INTERFACE POINT
               !----- A CORNER POINT
               ELSE
                  DATAL = LIST_GET_DATA(ELEM) !GET THE DATA OF THE LEFT INTERFACE POINT
                  CALL SETUTAU(DATAL)
                  DATAR = LIST_GET_DATA(ELEM%NEXT) !GET THE DATA OF THE LEFT INTERFACE POINT
                  CALL SETUTAU(DATAR)
                  CALL CMIB1D(DATAL,DATAR)
                  CALL LIST_PUT_DATA(ELEM,DATAL)      !PUT BACK THE UPDATED LEFT  INTERFACE POINT
                  CALL LIST_PUT_DATA(ELEM%NEXT,DATAR) !PUT BACK THE UPDATED RIGHT INTERFACE POINT
                  ELEM => ELEM%NEXT%NEXT !SKIP THE RIGHT INTERFACE POINT OF A CORNER POINT
               END IF
            END DO
         END IF

      END DO

      RETURN

      END SUBROUTINE FPSETUP

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SETJUMPS --
      !    SET JUMPS AT ALL INTERFACE PTS
      ! ARGUMENTS:
      !    T    IN   CURRENT TIME STEP
      !    DT   IN   TIME STEP
      !    UH   IN   NUMERICAL U AT CURRENT TIME STEP
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SETJUMPS(T,DT,UH)

      USE MOD_DATA

      REAL :: T,DT,UH(NY,NX)

      INTEGER :: IX,IY,I
      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LIST_DATA) :: DATA
      REAL :: SUMJ !NUM. U^+_TAU OR U^-_TAU
      REAL :: DTAU

      !----- RUNNING OVER EACH IX GRID LINE
      DO IX = 1,NX

         IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN

            ELEM => IFPX(IX)%HEAD

            DO WHILE ( ASSOCIATED(ELEM) )

               DATA = LIST_GET_DATA(ELEM)
               CALL PHI(T+DT,DATA)     !DATA%JUMP(1) = [U]              AT TIME STEP T^{N+1}
               CALL PSI(T+DT,DATA)     !DATA%JUMP(2) = [BETA U_N]       AT TIME STEP T^{N+1}
               CALL PHI_TAU(T+DT,DATA) !DATA%JUMP(3) = [U_TAU]          AT TIME STEP T^{N+1}
               CALL PSI_HAT(T,DT,DATA) !DATA%JUMP(4) = ANAL. [BEAT U_Y] AT TIME STEP T^{N+1}
                                       !DATA%JUMP(5) = ANAL. [BEAT U_Y] AT TIME STEP T^{N}

               !----- LOWER AUXILIARY VALUES
               IF (SETS(DATA%AUXL_X,DATA%AUXL_Y) .GT. TOL_SETUP) THEN !INDICATOR, PT IS IN OMEGA^+ OR OMEGA^-
                 DATA%UAUXL(1) =   1.0D0
               ELSE
                  DATA%UAUXL(1) = -1.0D0
               END IF
               DATA%UAUXL(2) = REAL(DATA%AUXL)                        !INDICATOR, U^+ OR U^-
               IF (DATA%AUXL .EQ. 1) THEN                             !ANA. VALUE
                  DATA%UAUXL(3) = UVAL_PLUS(DATA%AUXL_X,DATA%AUXL_Y,T)
               ELSE
                  DATA%UAUXL(3) = UVAL_MINUS(DATA%AUXL_X,DATA%AUXL_Y,T)
               END IF
               !DATA%UAUXL(3) = UVAL(DATA%AUXL_X,DATA%AUXL_Y,T)
               DATA%UAUXL(4) = 0.0D0                                  !NUM. VALUE
               DO I = 1,3
                  DATA%UAUXL(4) = DATA%UAUXL(4) + DATA%WAUXL(I)*UH(DATA%IAUXL(I),DATA%AUXL_AXID)
               END DO
               DATA%UAUXL(4) = -2.0D0*DATA%DAUXL*DATA%UAUXL(4)        !TIMES -2*DAUXL

               !----- UPPER AUXILIARY VALUES
               IF (SETS(DATA%AUXR_X,DATA%AUXR_Y) .GT. TOL_SETUP) THEN !INDICATOR, PT IS IN OMEGA^+ OR OMEGA^-
                  DATA%UAUXR(1) =  1.0D0
               ELSE
                  DATA%UAUXR(1) = -1.0D0
               END IF
               DATA%UAUXR(2) = REAL(DATA%AUXR)                        !INDICATOR, U^+ OR U^-
               IF (DATA%AUXR .EQ. 1) THEN                             !ANA. VALUE
                  DATA%UAUXR(3) = UVAL_PLUS(DATA%AUXR_X,DATA%AUXR_Y,T)
               ELSE
                  DATA%UAUXR(3) = UVAL_MINUS(DATA%AUXR_X,DATA%AUXR_Y,T)
               END IF
               !DATA%UAUXR(3) = UVAL(DATA%AUXR_X,DATA%AUXR_Y,T)
               DATA%UAUXR(4) = 0.0D0                                  !NUM. VALUE
               DO I = 1,3
                  DATA%UAUXR(4) = DATA%UAUXR(4) + DATA%WAUXR(I)*UH(DATA%IAUXR(I),DATA%AUXR_AXID)
               END DO
               DATA%UAUXR(4) = 2.0D0*DATA%DAUXR*DATA%UAUXR(4)         !TIMES 2*DAUXR

               !----- CENTRAL DIFFERENCE FOR CALCULATING U^+_TAU OR U^-_TAU USING 6 SUPPORTING GRIDS
               DTAU = DATA%DAUXL !CENTRAL DIFFERENCE, DATA%DAUXL = DATA%DAUXR

               SUMJ = 0.0D0
               DO I = 1,3
                  SUMJ = SUMJ + DATA%WAUXL(I)*UH(DATA%IAUXL(I),DATA%AUXL_AXID)
               END DO
               DO I = 1,3
                  SUMJ = SUMJ + DATA%WAUXR(I)*UH(DATA%IAUXR(I),DATA%AUXR_AXID)
               END DO

               !----- U^+_TAU IS USED FOR CALCULATING [BETA U_Y]
               IF (DATA%AUXL .EQ. 1) THEN
                  !----- ADD TWO MORE EXCESS TERMS WHEN DATA%AUXR HAS OPPOSITE SIGN
                  IF (DATA%AUXR .EQ. -1) THEN
                     SUMJ = SUMJ + 0.5D0/DTAU*DATA%JUMP(1) + 0.5D0*DATA%JUMP(3) !+ 1/(2L)PHI_IP + 1/2*(PHI_TAU)_IP
                  END IF

                  !----- SET ANALYTICAL U^+_TAU AND SAVE SUMJ FOR COMPARISON
                  DATA%UTAU(1) = UTAU_PLUS(DATA%X,DATA%Y,T,DATA%THETA) !ANA. U^+_TAU
                  DATA%UTAU(2) = SUMJ                                  !NUM. U^+_TAU

                  !+++++ TEST: REPLACE NUM. U^+_TAU W/ ANA. U^+_TAU
                  !SUMJ = DATA%UTAU(1)

                  !----- PSI^HAT := [BETA U_Y]
                  !               = SIN(THETA) PSI + COS(THETA) (BETA^+ - BETA^-) U^+_TAU + COS(THETA) BETA^- PHI_TAU
                  DATA%JUMP(6) = SIN(DATA%THETA)*DATA%JUMP(2) + COS(DATA%THETA)*(DATA%BETA2 - DATA%BETA1)*SUMJ + &
                                &COS(DATA%THETA)*DATA%BETA1*DATA%JUMP(3)

               !----- U^-_TAU IS USED FOR CALCULATING [BETA U_Y]
               ELSE
                  !----- ADD TWO MORE EXCESS TERMS WHEN DATA%AUXR HAS OPPOSITE SIGN
                  IF (DATA%AUXR .EQ. 1) THEN
                     SUMJ = SUMJ - 0.5D0/DTAU*DATA%JUMP(1) - 0.5D0*DATA%JUMP(3) !- 1/(2L)PHI_IP - 1/2*(PHI_TAU)_IP
                  END IF

                  !----- SET ANALYTICAL U^-_TAU AND SAVE SUMJ FOR COMPARISON
                  DATA%UTAU(3) = UTAU_MINUS(DATA%X,DATA%Y,T,DATA%THETA) !ANA. U^-_TAU
                  DATA%UTAU(4) = SUMJ                                   !NUM. U^-_TAU

                  !+++++ TEST - REPLACE NUM. U^-_TAU W/ ANA. U^-_TAU
                  !SUMJ = DATA%UTAU(3)

                  !----- PSI^HAT := [BETA U_Y]
                  !               = SIN(THETA) PSI + COS(THETA) (BETA^+ - BETA^-) U^-_TAU + COS(THETA) BETA^+ PHI_TAU
                  DATA%JUMP(6) = SIN(DATA%THETA)*DATA%JUMP(2) + COS(DATA%THETA)*(DATA%BETA2 - DATA%BETA1)*SUMJ + &
                                &COS(DATA%THETA)*DATA%BETA2*DATA%JUMP(3)

               END IF

               !----- SAVE DATA TO THE LIST AND MOVE ON TO THE NEXT INTERFACE PT
               CALL LIST_PUT_DATA(ELEM,DATA)
               ELEM => ELEM%NEXT ! NEXT INTERFACE

            END DO
         END IF
      END DO

      !----- RUNNING OVER EACH IY GRID LINE
      DO IY = 1,NY

         IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN

            ELEM => IFPY(IY)%HEAD

            DO WHILE ( ASSOCIATED(ELEM) )

               DATA=LIST_GET_DATA(ELEM)
               CALL PHI(T+DT,DATA)     !DATA%JUMP(1) = [U]              AT TIME STEP T^{N+1}
               CALL PSI(T+DT,DATA)     !DATA%JUMP(2) = [BETA U_N]       AT TIME STEP T^{N+1}
               CALL PHI_TAU(T+DT,DATA) !DATA%JUMP(3) = [U_TAU]          AT TIME STEP T^{N+1}
               CALL PSI_BAR(T,DT,DATA) !DATA%JUMP(4) = ANAL. [BEAT U_X] AT TIME STEP T^{N}

               !----- LEFT AUXILIARY VALUES
               IF (SETS(DATA%AUXL_X,DATA%AUXL_Y) .GT. TOL_SETUP) THEN !INDICATOR, PT IS IN OMEGA^+ OR OMEGA^-
                  DATA%UAUXL(1) =  1.0D0
               ELSE
                  DATA%UAUXL(1) = -1.0D0
               END IF
               DATA%UAUXL(2) = REAL(DATA%AUXL)                        !INDICATOR, U^+ OR U^-
               IF (DATA%AUXL .EQ. 1) THEN                             !ANA. VALUE
                  DATA%UAUXL(3) = UVAL_PLUS(DATA%AUXL_X,DATA%AUXL_Y,T)
               ELSE
                  DATA%UAUXL(3) = UVAL_MINUS(DATA%AUXL_X,DATA%AUXL_Y,T)
               END IF
               !DATA%UAUXL(3) = UVAL(DATA%AUXL_X,DATA%AUXL_Y,T)
               DATA%UAUXL(4) = 0.0D0                                  !NUM. VALUE
               DO I = 1,3
                  DATA%UAUXL(4) = DATA%UAUXL(4) + DATA%WAUXL(I)*UH(DATA%AUXL_AXID,DATA%IAUXL(I))
               END DO
               DATA%UAUXL(4) = -2.0D0*DATA%DAUXL*DATA%UAUXL(4)        !TIMES -2*DAUXL

               !----- RIGHT AUXILIARY VALUES
               IF (SETS(DATA%AUXR_X,DATA%AUXR_Y) .GT. TOL_SETUP) THEN !INDICATOR, PT IS IN OMEGA^+ OR OMEGA^-
                  DATA%UAUXR(1) =  1.0D0
               ELSE
                  DATA%UAUXR(1) = -1.0D0
               END IF
               DATA%UAUXR(2) = REAL(DATA%AUXR)                        !INDICATOR, U^+ OR U^-
               IF (DATA%AUXR .EQ. 1) THEN                             !ANA. VALUE
                  DATA%UAUXR(3) = UVAL_PLUS(DATA%AUXR_X,DATA%AUXR_Y,T)
               ELSE
                  DATA%UAUXR(3) = UVAL_MINUS(DATA%AUXR_X,DATA%AUXR_Y,T)
               END IF
               !DATA%UAUXR(3) = UVAL(DATA%AUXR_X,DATA%AUXR_Y,T)
               DATA%UAUXR(4) = 0.0D0                                  !NUM. VALUE
               DO I = 1,3
                  DATA%UAUXR(4) = DATA%UAUXR(4) + DATA%WAUXR(I)*UH(DATA%AUXR_AXID,DATA%IAUXR(I))
               END DO
               DATA%UAUXR(4) = 2.0D0*DATA%DAUXR*DATA%UAUXR(4)         !TIMES 2*DAUXR

               !----- CENTRAL DIFFERENCE FOR CALCULATING U^+_TAU OR U^-_TAU USING 6 SUPPORTING GRIDS
               DTAU = DATA%DAUXL !CENTRAL DIFFERENCE, DATA%DAUXL = DATA%DAUXR

               SUMJ = 0.0D0
               DO I = 1,3
                  SUMJ = SUMJ + DATA%WAUXL(I)*UH(DATA%AUXL_AXID,DATA%IAUXL(I))
               END DO
               DO I = 1,3
                  SUMJ = SUMJ + DATA%WAUXR(I)*UH(DATA%AUXR_AXID,DATA%IAUXR(I))
               END DO

               !----- U^+_TAU IS USED FOR CALCULATING [BETA U_X]
               IF (DATA%AUXL .EQ. 1) THEN
                  !----- ADD TWO MORE EXCESS TERMS WHEN DATA%AUXR HAS OPPOSITE SIGN
                  IF (DATA%AUXR .EQ. -1) THEN
                     SUMJ = SUMJ + 0.5D0/DTAU*DATA%JUMP(1) + 0.5D0*DATA%JUMP(3) !+ 1/(2L)PHI_IP + 1/2*(PHI_TAU)_IP
                  END IF

                  !----- SET ANALYTICAL U^+_TAU AND SAVE SUMJ FOR COMPARISON
                  DATA%UTAU(1) = UTAU_PLUS(DATA%X,DATA%Y,T,DATA%THETA) !ANA. U^+_TAU
                  DATA%UTAU(2) = SUMJ                                  !NUM. U^+_TAU

                  !+++++ TEST: REPLACE NUM. U^+_TAU W/ ANA. U^+_TAU
                  !SUMJ = DATA%UTAU(1)

                  !----- PSI^BAR := [BETA U_X]
                  !               = COS(THETA) PSI - SIN(THETA) (BETA^+ - BETA^-) U^+_TAU - SIN(THETA) BETA^- PHI_TAU
                  DATA%JUMP(6) = COS(DATA%THETA)*DATA%JUMP(2) - SIN(DATA%THETA)*(DATA%BETA2 - DATA%BETA1)*SUMJ - &
                                &SIN(DATA%THETA)*DATA%BETA1*DATA%JUMP(3)

               !----- U^-_TAU IS USED FOR CALCULATING [BETA U_X]
               ELSE
                  !----- ADD TWO MORE EXCESS TERMS WHEN DATA%AUXR HAS OPPOSITE SIGN
                  IF (DATA%AUXR .EQ. 1) THEN
                     SUMJ = SUMJ - 0.5D0/DTAU*DATA%JUMP(1) - 0.5D0*DATA%JUMP(3) !- 1/(2L)PHI_IP - 1/2*(PHI_TAU)_IP
                  END IF

                  !----- SET ANALYTICAL U^-_TAU AND SAVE SUMJ FOR COMPARISON
                  DATA%UTAU(3) = UTAU_MINUS(DATA%X,DATA%Y,T,DATA%THETA) !ANA. U^-_TAU
                  DATA%UTAU(4) = SUMJ                                   !NUM. U^+_TAU

                  !+++++ TEST - REPLACE NUM. U^-_TAU W/ ANA. U^-_TAU
                  !SUMJ = DATA%UTAU(3)

                  !----- PSI^BAR := [BETA U_X]
                  !               = COS(THETA) PSI - SIN(THETA) (BETA^+ - BETA^-) U^-_TAU - SIN(THETA) BETA^+ PHI_TAU
                  DATA%JUMP(6) = COS(DATA%THETA)*DATA%JUMP(2) - SIN(DATA%THETA)*(DATA%BETA2 - DATA%BETA1)*SUMJ - &
                                &SIN(DATA%THETA)*DATA%BETA2*DATA%JUMP(3)

               END IF

               !----- SAVE DATA TO THE LIST AND MOVE ON TO THE NEXT INTERFACE PT
               CALL LIST_PUT_DATA(ELEM,DATA)
               ELEM => ELEM%NEXT ! NEXT INTERFACE

            END DO
         END IF
      END DO

      RETURN

      END SUBROUTINE SETJUMPS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! WEIGHTS --
      !    WEIGHTS OF THE LAGRANGE POLYNOMIAL INTERPOLATION
      ! ARGUMENTS:
      !    Z   IN   LOCATION WHERE APPROXIMATIONS ARE TO BE ACCURATE
      !    X   IN   X(0:ND) GRID POINT LOCATIONS, FOUND IN X(0:N)
      !    N   IN   ONE LESS THAN TOTAL NUMBER OF GRID POINTS; N MUST NOT EXCEED
      !             THE PARAMETER ND BELOW
      !    ND  IN   DIMENSION OF X- AND C-ARRAYS IN CALLING PROGRAM X(0:ND) AND
      !             C(0:ND,0:M), RESPECTIVELY,
      !    M   IN   HIGHEST DERIVATIVE FOR WHICH WEIGHTS ARE SOUGHT,
      !    C   OUT  C(0:ND,0:M) WEIGHTS AT GRID LOCATIONS X(0:N) FOR DERIVATIVES
      !             OF ORDER 0:M, FOUND IN C(0:N,0:M)
      ! NOTES:
      !    FORNBERG, BENGT., "CLASSROOM NOTE: CALCULATION OF WEIGHTS IN FINITE
      !    DIFFERENCE FORMULAS.", SIAM REVIEW 40.3 (1998): 685-691.
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE WEIGHTS(Z,X,N,ND,M,C)

      REAL :: Z,X(0:ND),C(0:ND,0:M)
      INTEGER :: N,ND,M

      REAL :: C1,C2,C3,C4,C5
      INTEGER :: I,J,K,MN

      C1 = 1.0D0
      C4 = X(0)-Z

      DO K = 0,M
         DO J = 0,N
            C(J,K) = 0.0D0
         END DO
      END DO

      C(0,0) = 1.0D0

      DO I = 1,N
         MN = MIN(I,M)
         C2 = 1.0D0
         C5 = C4
         C4 = X(I)-Z
         DO J = 0,I-1
            C3 = X(I)-X(J)
            C2 = C2*C3
            IF (J .EQ. I-1) THEN
               DO K = MN,1,-1
                  C(I,K) = C1*(K*C(I-1,K-1)-C5*C(I-1,K))/C2
               END DO
               C(I,0) = -C1*C5*C(I-1,0)/C2
            END IF
            DO K = MN,1,-1
               C(J,K) = (C4*C(J,K)-K*C(J,K-1))/C3
            END DO
            C(J,0) = C4*C(J,0)/C3
         END DO
         C1 = C2
      END DO

      RETURN

      END SUBROUTINE WEIGHTS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! GETWMIB --
      !    DIFFERENT WEIGHTS GENERATED BASED ON ITYPE ONLY UES IN MIB PART
      ! ARGUMENTS:
      !    Z         IN    INTERFACE LOCATION FOR TWO FPS
      !    ORDER     IN    DERIVATIVE ORDER OF EXPECTED WEIGHTS
      !    N         IN    ONE LESS THAN # OF GRID POINTS INVOLVED TO GET WEIGHTS
      !    IYTPE     IN    CURRENT INTERFACE ITYPE (CAN BE LEFT OR RIGHT INTERFACE)
      !    ITYPE2    IN    %%OPTIONAL FOR CORNER POINT%% SECOND INTERFACE ITYPE
      !    WEI       OUT   MODIFIED WEIGHTS FOR NONUNIFORM GRID
      ! NOTES:
      !    USED IN SUBROUTINE MIB TO GENERATE WEIGHTS FOR DIFFERENT TYPES OF INTERFACES
      !    LOCAL ORIGINAL ALWAYS LEFT GRID POINT OF CURRENT INTERFACE
      ! LIMIT:
      !    ONLY FOR 3 OR 4 POINTS FOR EACH INTERFACE
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE GETWMIB(Z,ORDER,WEI,N,ITYPE,ITYPE2)

      USE MOD_DATA

      INTEGER :: ITYPE,ORDER,N
      INTEGER,OPTIONAL :: ITYPE2
      REAL :: WEI(2,0:N,0:ORDER)
      REAL :: Z

      !FOR CORNER INTERFACE, 1: INTERPOLATION INSIDE, FP1,FP2,FP4; 2: INTERPOLATION OUTSIDE, FP3
      !FOR IRREGULAR INTERFACE, 1: INTERPOLATION FP1; 2: INTERPOLATION FP2
      REAL :: X1(0:N),X2(0:N)
      REAL :: C1(0:N,0:ORDER),C2(0:N,0:ORDER)
      INTEGER :: KEY,I,J

      IF (PRESENT(ITYPE2)) THEN           !CORNER POINT
         IF (ITYPE .LT. 0) THEN           !(-2,0)&(-2,1)&(-1,2)
            KEY=1
         ELSE IF (ITYPE .GT. 0) THEN      !(2,0)&(2,-1)&(1,-2)
            KEY=2
         ELSE
            IF (ITYPE2 .LT. 0) THEN       !(0,-2)
               KEY=2
            ELSE                          !(0,2)
               KEY=1
            END IF
         END IF

         SELECT CASE(KEY)
         CASE (1)      !IT'S LEFT INTERFACE
            !OUTSIDE INTERPOLATION NEVER CHANGE WITH ITYPE
            DO I=0,N
               X2(I)=-DX+I*DX
               X1(I)=X2(I)+DX                !(0,2)
            END DO
            !INSIDE INTERPOLATION MODIFICATION
            IF (ITYPE .EQ. -2) THEN          !CLOSE TO LEFT INTERFACE, 4TH FP CHANGES, (-2,0)
               X1(N)=-DX
               IF (ITYPE2 .EQ. 1) THEN
                  X1(N-1)=X1(N-1)+DX         !FP3 MOVE RIGHT, (-2,1)
               END IF
            ELSE                             !CLOSE TO RIGHT INTERFACE, 4TH FP SAME
               IF (ITYPE .EQ. -1) THEN       !FP1 MOVE LEFT, (-1,2)
                  X1(0)=X1(0)-DX
               END IF
            END IF
        CASE (2)       !IT'S RIGHT INTERFACE
            !OUTSIDE INTERPOLATION NEVER CHANGE WITH ITYPE
            DO I=0,N
               X2(I)=-DX+I*DX
               X1(I)=X2(I)                   !(2,0)
            END DO
            !INSIDE INTERPOLATION MODIFICATION
            IF (ITYPE .EQ. 2) THEN           !CLOSE TO RIGHT INTERFACE, 4TH FP SAME
               IF (ITYPE2 .EQ. -1) THEN      !FP1 MOVE LEFT, (2,-1)
                  X1(0)=X1(0)-DX
               END IF
            ELSE                             !CLOSE TO RIGHT INTERFACE, 4TH FP CHANGE, (0,-2)
               X1(N)=-2.0D0*DX
               IF (ITYPE .EQ. 1) THEN        !FP1 MOVE RIGHT, (1,-2)
                  X1(N-1)=X1(N-1)+DX
               END IF
            END IF
         END SELECT
      ELSE                                 !IRREGULAR POINT
         DO I=0,N
            X2(I)=-DX+I*DX
            X1(I)=X2(I)+DX
         END DO
         IF (ITYPE .EQ. -1) THEN           !FP1 MOVE LEFT ONE UNIT
            X1(0)=X1(0)-DX
         ELSE IF (ITYPE .EQ. 1) THEN       !FP2 MOVE RIGHT ONE UNIT
            X2(N)=X2(N)+DX
         END IF
      END IF
      CALL WEIGHTS(Z,X1,N,N,ORDER,C1)
      CALL WEIGHTS(Z,X2,N,N,ORDER,C2)
      DO J=0,ORDER
         DO I=0,N
            WEI(1,I,J)=C1(I,J)
            WEI(2,I,J)=C2(I,J)
         END DO
      END DO

      RETURN
      END SUBROUTINE GETWMIB

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !  GRAPH EXPLAINATION FOR DIFFERENT SITUATION USED IN MIB
      !  -------- : MESH LINES
      !      *    : REAL POINTS
      !      O    : FP
      !      x    : INTERFACE
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  CORNER LEFT INTERFACE  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     (ITYPE,ITYPE2)
      !  ---------*---------------*-------x--------*--------x-------*---------------*----------
      !  OUTSIDE INTERPOLATION, X2
      !  ---------*---------------*-------x--------O--------x-------*--------------------------
      !  INSIDE INTERPOLATION, X1
      !  ---------O---------------O---x------------*--------x-------O--------------------------     ( -2 , 0 )
      !  ---------O---------------O---x------------*--------------x-----------------O----------     ( -2 , 1 )
      !  -------------------------O-------x--------*------------x---O---------------O----------     (  0 , 2 )
      !  ---------O-----------------x--------------*------------x---O---------------O----------     ( -1 , 2 )

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  CORNER RIGHT INTERFACE  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     (ITYPE,ITYPE2)
      !  ---------*---------------*-------x--------*--------x-------*---------------*----------
      !  OUTSIDE INTERPOLATION, X2
      !  -------------------------*-------x--------O--------x-------*---------------*----------
      !  INSIDE INTERPOLATION, X1
      !  ---------O---------------O---x------------*--------x-------O--------------------------     ( 0 , -2 )
      !  ---------O---------------O---x------------*--------------x-----------------O----------     ( 1 , -2 )
      !  -------------------------O-------x--------*------------x---O---------------O----------     ( 2 ,  0 )
      !  ---------O-----------------x--------------*------------x---O---------------O----------     ( 2 , -1 )

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  IRREGULAR INTERFACE  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !  ITYPE == 0
      !  -----------*---------------*--------x-------*---------------*-------------
      !  ---------------------------O--------x-------*---------------*-------------                 FP1 INTERPOLATION
      !  -----------*---------------*--------x-------O-----------------------------                 FP2 INTERPOLATION
      !  ITYPE == -1
      !  -----------*---------------*-x--------------*---------------*-------------
      !  -----------O-----------------x--------------*---------------*-------------                 FP1 INTERPOLATION
      !  -----------*---------------*-x--------------O-----------------------------                 FP2 INTERPOLATION
      !  ITYPE == 1
      !  -----------*---------------*--------------x-*---------------*-------------
      !  ---------------------------O--------------x-*---------------*-------------                 FP1 INTERPOLATION
      !  -----------*---------------*--------------x-----------------O-------------                 FP2 INTERPOLATION

