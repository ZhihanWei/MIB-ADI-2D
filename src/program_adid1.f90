      PROGRAM HEAT2D_ADID1_MIB2D

      USE MOD_DATA
      USE MOD_MIB2D, ONLY : SETUP,FPSETUP,SETJUMPS,TEST_MIB1D,TEST_CMIB1D
      USE MOD_TIME_STEPPER, ONLY : ADID1

      IMPLICIT NONE

      !---------- NUMERICAL AND EXACT SOLUTIONS ----------
      REAL, ALLOCATABLE :: U(:,:)  !Analytical solution
      REAL, ALLOCATABLE :: UH(:,:) !Numerical solution

      !---------- STABILITY ANALYSIS ----------
!      REAL,    ALLOCATABLE :: A(:)          !nonzero elements of matrix A
!      INTEGER, ALLOCATABLE :: IA(:), JA(:)  !(i,j) index of nonzero elements
!      REAL,    ALLOCATABLE :: B(:)          !nonzero elements of matrix B
!      INTEGER, ALLOCATABLE :: IB(:), JB(:)  !(i,j) index of nonzero elements

      !---------- LOCAL VARIABLES ----------
      REAL :: T,TEND,DT,XR,YR
      INTEGER :: IX,IY,IERR,NUMT,NLOOP,IPRINT,ILOOP

      TYPE(LIST_DATA) :: DATA !DATA OF AN IRREGULAR POINT
      TYPE(LIST_DATA) :: DATAL,DATAR !DATA OF A CORNER POINT
      TYPE(LINKED_LIST), POINTER  :: ELEM

      IF (TESTCASE .EQ. 1) THEN
         WRITE(*,*) " "
         WRITE(*,"(A,I2,A)") "----- SOLVING EXAMPLE <DISC_EX1> WITH METHOD <ADID1> -----"
         WRITE(*,*) " "
      ELSE IF (TESTCASE .EQ. 2) THEN
         WRITE(*,*) " "
         WRITE(*,"(A,I2,A)") "----- SOLVING EXAMPLE <DISC_EX2> WITH METHOD <ADID1> -----"
         WRITE(*,*) " "
      ELSE IF (TESTCASE .EQ. 3) THEN
         WRITE(*,*) " "
         WRITE(*,"(A,I2,A)") "----- SOLVING EXAMPLE <DISC_EX3> WITH METHOD <ADID1> -----"
         WRITE(*,*) " "
      ELSE IF (TESTCASE .EQ. 4) THEN
         WRITE(*,*) " "
         WRITE(*,"(A,I2,A)") "----- SOLVING EXAMPLE <STAR_EX4> WITH METHOD <ADID1> -----"
         WRITE(*,*) " "
      ELSE IF (TESTCASE .EQ. 5) THEN
         WRITE(*,*) " "
         WRITE(*,"(A,I2,A)") "----- SOLVING EXAMPLE <STAR_EX5> WITH METHOD <ADID1> -----"
         WRITE(*,*) " "
      ELSE IF (TESTCASE .EQ. 6) THEN
         WRITE(*,*) " "
         WRITE(*,"(A,I2,A)") "----- SOLVING EXAMPLE <STAR_EX6> WITH METHOD <ADID1> -----"
         WRITE(*,*) " "
      ELSE IF (TESTCASE .EQ. 7) THEN
         WRITE(*,*) " "
         WRITE(*,"(A,I2,A)") "----- SOLVING EXAMPLE <STAR_EX7> WITH METHOD <ADID1> -----"
         WRITE(*,*) " "
      ELSE
         WRITE(*,*) " "
         WRITE(*,"(2(A,I2),A)") "   PROBLEM OCCURS WHEN SOLVING EXAMPLE ",TESTCASE," WITH METHOD <ADIPR>"
         WRITE(*,*) " "
         STOP
      END IF

      !---------- SPATIAL DOMAIN DISCRETIZATION ----------
      XL = -CD; XR = CD !X dimension
      YL = -CD; YR = CD !Y dimension
      DX =  (XR-XL)/(NX-1.0D0)

      !---------- TIME DOMAIN DISCRETIZATION ----------
      T = TSTART; TEND = TFINAL; DT = TSTEP !INITIAL, STOPPING TIME AND TIME STEP
      NUMT = ANINT( (TEND-T)/DT ); NLOOP = NUMT/NPRINT

      !---------- VARIABLES DEFINED IN MOD_MIB ----------
      INODE = 0
      SRC   = 0.0D0
      BETA  = 0.0D0; BETAX = 0.0D0; BETAY = 0.0D0
      TOL_ITYPE = 0.0D0
      !TOL_ITYPE = .0D-1*DX

      CALL SETUP !Initial setup

      WRITE(*,*) "---------------------------------------------------------------------------------"
      WRITE(*,"(2(A,E12.6))") "TOL_PERCENT = ",TOL_ITYPE/DX, ", TOL_ITYPE: ",TOL_ITYPE
      WRITE(*,"(2(A,E12.6),(A,I5))") "XL   = ",XL,  ", DX = ",DX,  ", NX   = ",NX
      WRITE(*,"(2(A,E12.6),(A,I5))") "YL   = ",YL,  ", DY = ",DX,  ", NY   = ",NY
      WRITE(*,"(2(A,E12.6),(A,I5))") "TEND = ",TEND,", DT = ",DT,  ", NUMT = ",NUMT

      CALL FPSETUP !FICTITIOUS POINTS SETUP AT INITIAL TIME

      !-------------------- STABILITY ANALYSIS STARTS--------------------!
!      N2D=NX*NY
!      MA=N2D*6  !Max # of nonzero elements
!      NA=0      !actual # of nonzero elements
!      ALLOCATE(A(MA),IA(MA),JA(MA),STAT=ierr)
!      A=0.0D0; IA=0; JA=0 !INITIALIZATION
!
!      MB=N2D*7+40*NX !Max # of nonzero elements
!      NB=0           !actual # of nonzero elements
!      ALLOCATE(B(MB),IB(MB),JB(MB),STAT=ierr)
!      B=0.0D0; IB=0; JB=0 !INITIALIZATION
!
!      CALL ADIPRSETUP(DT,N2D,MA,NA,A,IA,JA,MB,NB,B,IB,JB) !set up matrix
!
!      TOL=1.D-14   !tolerance
!      MAXTER=10000 !max steps
!
!      CALL EIGENSTAB(N2D,TOL,MAXTER,NA,A,IA,JA,NB,B,IB,JB)
!
!      STOP
      !-------------------- STABILITY ANALYSIS ENDS--------------------!

      ALLOCATE(U(NY,NX),UH(NY,NX),STAT=ierr)
      CALL ANALYTICAL(UH,T) !INITIAL SOLUTION

      !+++++ TEST: TIME STEPS
!      CALL TEST_ALLIFPS(T)
!      STOP

      DO IPRINT = 1,NPRINT
         DO ILOOP = 1,NLOOP

            CALL SETJUMPS(T,DT,UH)

            !+++++ TEST: INTERFACE PTS
            IF (IPRINT .EQ. 1 .AND. ILOOP .EQ. 1) THEN

               !ITERATION TO INITIALIZE ALL ERRORS OF FPS
               DO IY = 2,NY-2
                  IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN
                     ELEM => IFPY(IY)%HEAD
                     DO WHILE ( ASSOCIATED(ELEM) )
                        IF (ELEM%DATA%ID .GT. 0) THEN
                           DATA = LIST_GET_DATA(ELEM)
                           CALL TEST_MIB1D(DATA)
                           CALL LIST_PUT_DATA(ELEM,DATA)
                           ELEM => ELEM%NEXT
                        ELSE
                           DATAL = LIST_GET_DATA(ELEM)
                           DATAR = LIST_GET_DATA(ELEM%NEXT)
                           CALL TEST_CMIB1D(DATAL,DATAR)
                           CALL LIST_PUT_DATA(ELEM,DATAL)
                           CALL LIST_PUT_DATA(ELEM%NEXT,DATAR)
                           ELEM => ELEM%NEXT%NEXT
                        END IF
                     END DO
                  END IF
               END DO

               DO IX = 2,NX-2
                  IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN
                     ELEM => IFPX(IX)%HEAD
                     DO WHILE ( ASSOCIATED(ELEM) )
                        IF ( ELEM%DATA%ID .GT. 0) THEN
                           DATA = LIST_GET_DATA(ELEM)
                           CALL TEST_MIB1D(DATA)
                           CALL LIST_PUT_DATA(ELEM,DATA)
                           ELEM => ELEM%NEXT
                        ELSE
                           DATAL = LIST_GET_DATA(ELEM)
                           DATAR = LIST_GET_DATA(ELEM%NEXT)
                           CALL TEST_CMIB1D(DATAL,DATAL)
                           CALL LIST_PUT_DATA(ELEM%NEXT,DATAR)
                           ELEM => ELEM%NEXT%NEXT
                        END IF
                     END DO
                  END IF
               END DO

               CALL TEST_ALLIFPS(T)
            END IF

            CALL ADID1(T,DT,UH)
            T = T + DT
         END DO !----- END OF ILOOP = 1,NLOOP

         !+++++ TEST: MAX. JUMP ERRORS
         !CALL TEST_JUMPS

         CALL ANALYTICAL(U,T)
         CALL OUTERROR(T,U,UH)
      END DO !----- END OF IPRINT = 1,NPRINT

      DEALLOCATE(U,UH,STAT=ierr)

      END PROGRAM HEAT2D_ADID1_MIB2D
