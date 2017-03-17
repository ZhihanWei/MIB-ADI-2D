      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! MIB1D --
      !    SECOND ORDER MIB IN 1D
      ! ARGUMENTS:
      !    DATA IN/OUT  AN ELEMENT OF THE LIST FOR ONE INTERFACE POINT
      ! NOTES:
      !    WORKS FOR BOTH X- AND Y- DIRECTION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE MIB1D(DATA)

      USE MOD_DATA

      TYPE(LIST_DATA) :: DATA

      REAL,ALLOCATABLE :: AE(:,:),BE(:),VVE(:)           !AE*U=BE, U->EMAT
      REAL,ALLOCATABLE :: WEI(:,:,:),WEI1(:,:),WEI2(:,:) !1: FP1; 2: FP2
      INTEGER,ALLOCATABLE :: INDEXE(:)

      REAL :: EMAT(2,6)
      REAL :: DE,COEF1,COEF2
      INTEGER :: M,L,ME,NE,ND,NS,I,J,IM,IW,IERR

      M = 2          !# OF FICTITIOUS POINTS
      L = 2          !# OF GRID POINTS IS L*2
      ! TWO EXTRA ELEMENTS ARE FOR [U] AND [BETA U_X]
      ME = L*2+2     !UNKNOWN WEIGHTS FOR SINGLE FP
      NE = M*ME      !TOTAL UNKNOWN WEIGHTS
      ND = 1         !CONSIDER MAXIMUM 1ST ORDER JUMP CONDITION
      NS = L+M/2     !LENGTH OF EACH FD STENCIL

      ALLOCATE(AE(NE,NE),BE(NE),INDEXE(NE),VVE(NE),WEI1(NS,0:ND),WEI2(NS,0:ND),WEI(2,0:NS-1,0:ND),STAT=IERR)

      !INITIALIZATION
      AE = 0.0D0
      BE = 0.0D0

      CALL GETWMIB(DATA%GAMMA,ND,WEI,NS-1,DATA%ITYPE)

      DO I = 0,ND
         DO J = 1,NS
            WEI1(J,I) = WEI(1,J-1,I)
            WEI2(J,I) = WEI(2,J-1,I)
         END DO
      END DO

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !   ZERO ORDER:   U^- = U^+ - [U]
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IM = 1       !DISCRETIZE THE FIRST LINEAR EQUATION
      IW = 0       !ZERO ORDER DERIVATIVE

      !FP2
      DO I = 1,NS-1
         BE( (IM-1)*ME+I ) = -WEI2(I,IW)
      END DO

      DO I = 1,ME
         AE( (IM-1)*ME+I,ME+I ) = WEI2(NS,IW)
      END DO

      !FP1
      DO I = 2,NS
         BE( (IM-1)*ME+I+1 ) = WEI1(I,IW)
      END DO

      DO I = 1,ME
         AE( (IM-1)*ME+I,I ) = -WEI1(1,IW)
      END DO

      IF (MOD(DATA%ID,2) .EQ. 0) THEN         !"-" => "+", U^{-}=U^{+}-[U], F=G-[U]
         BE( (IM-1)*ME+L*2+1 ) = -1.0D0
      ELSE                                    !"+" => "-", U^{-}=U^{+}-[U], G+[U]=F
         BE( (IM-1)*ME+L*2+1 ) =  1.0D0
      END IF

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !   FIRST ORDER:  BETA^- U^-_X = BETA^+ U^+_X - [BETA U_X]
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IM = 2       !DISCRETIZE THE SECOND LINEAR EQUATION
      IW = 1       !FIRST ORDER DERIVATIVE

      IF ( MOD(DATA%ID,2) .EQ. 0 ) THEN     !"-" => "+", BETA^{-}U_X^{-}=BETA^{+}U_X^{+}-[BETA U_X], F=G-[BETA U_X]
         COEF1 = DATA%BETA1
         COEF2 = DATA%BETA2
      ELSE                                  !"+" => "-", BETA^{-}U_X^{-}=BETA^{+}U_X^{+}-[BETA U_X], G+[BETA U_X]=F
         COEF1 = DATA%BETA2
         COEF2 = DATA%BETA1
      END IF

      !FP2
      DO I = 1,NS-1
         BE( (IM-1)*ME+I ) = -WEI2(I,IW)*COEF1
      END DO
      DO I = 1,ME
         AE( (IM-1)*ME+I,ME+I ) = WEI2(NS,IW)*COEF1
      END DO

      !FP1
      DO I = 2,NS
         BE( (IM-1)*ME+I+1 ) = WEI1(I,IW)*COEF2
      END DO
      DO I = 1,ME
         AE( (IM-1)*ME+I,I ) = -WEI1(1,IW)*COEF2
      END DO

      IF (MOD(DATA%ID,2) .EQ. 0) THEN         !"-" => "+", BETA^{-}U_X^{-}=BETA^{+}U_X^{+}-[BETA U_X], F=G-[BETA U_X]
         BE( (IM-1)*ME+L*2+2 ) = -1.0D0
      ELSE                                    !"+" => "-", BETA^{-}U_X^{-}=BETA^{+}U_X^{+}-[BETA U_X], G+[BETA U_X]=F
         BE( (IM-1)*ME+L*2+2 ) =  1.0D0
      END IF

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !   LU DECOMPOSITION
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL LUDCMP(AE,NE,NE,INDEXE,DE,VVE)
      CALL LUBKSB(AE,NE,NE,INDEXE,BE)

      DO I = 1,M
         DO J = 1,ME
            EMAT(I,J) = BE( (I-1)*ME+J )
            DATA%WIJ(I,J) = EMAT(I,J)
         END DO
      END DO

      DEALLOCATE(AE,BE,VVE,INDEXE,WEI,WEI1,WEI2,STAT=IERR)

      RETURN

      END SUBROUTINE MIB1D

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! CMIB1D --
      !    SECOND ORDER MIB IN 1D FOR A CORNER POINT
      ! ARGUMENTS:
      !    DATAL IN/OUT  LEFT  INTERFACE POINT OF A CORNER POINT
      !    DATAR IN/OUT  RIGHT INTERFACE POINT OF A CORNER POINT
      ! NOTES:
      !    WORKS FOR BOTH X- AND Y- DIRECTION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE CMIB1D(DATAL,DATAR)

      USE MOD_DATA

      TYPE(LIST_DATA) :: DATAL,DATAR

      REAL,ALLOCATABLE :: AE(:,:),BE(:),VVE(:)            !AE*U=BE, U->EMAT
      REAL,ALLOCATABLE :: WEI(:,:,:),WEI1(:,:),WEI2(:,:)  !1: OUTSIDE OF CORNER; 2: INSIDE OF CORNER
      INTEGER,ALLOCATABLE :: INDEXE(:)

      REAL :: EMAT(4,9)
      REAL :: DE,BETAOUT,BETAIN
      INTEGER :: M,ME,NE,ND,NS,I,J,IM,IW,IERR

      M  = 4           !# OF FICTITIOUS POINTS
      ME = 9           !UNKNOWN WEIGHTS FOR SINGLE FP
      NE = M*ME        !TOTAL UNKNOWN WEIGHTS
      ND = 1           !CONSIDER MAXIMUM 1ST ORDER JUMP CONDITION
      NS = 4           !LENGTH OF EACH FD STENCIL

      ALLOCATE(AE(NE,NE),BE(NE),INDEXE(NE),VVE(NE),WEI1(NS,0:ND),WEI2(NS,0:ND),WEI(2,0:NS-1,0:ND),STAT=IERR)

      !INITIALIZATION
      AE = 0.0D0
      BE = 0.0D0

      !###################### MIB FOR LEFT INTERFACE ############################
      CALL GETWMIB(DATAL%GAMMA,ND,WEI,NS-1,DATAL%ITYPE,DATAR%ITYPE)

      DO I = 0,ND
         DO J = 1,NS
            WEI1(J,I) = WEI(1,J-1,I)
            WEI2(J,I) = WEI(2,J-1,I)
         END DO
      END DO

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !   ZERO ORDER:   U^- = U^+ - [U]
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IM=1       !DISCRETIZE THE FIRST LINEAR EQUATION
      IW=0       !ZERO ORDER DERIVATIVE

      !OUTSIDE INTERPOLATION
      BE((IM-1)*ME+1) = -WEI2(1,IW)                         !GP AT IX-1
      BE((IM-1)*ME+2) = -WEI2(2,IW)                         !GP AT IX
      DO I=1,ME
         AE( (IM-1)*ME+I,(2-1)*ME+I ) = WEI2(3,IW)          !FP AT IX+1, SECOND FP
      END DO
      BE( (IM-1)*ME+4 ) = -WEI2(4,IW)                       !GP AT IX+2

      !INSIDE INTERPOLATION
      DO I = 1,ME
         AE( (IM-1)*ME+I,(1-1)*ME+I ) = -WEI1(1,IW)         !FP AT IX, FIRST FP
      END DO
      BE( (IM-1)*ME+3 ) = WEI1(2,IW)                        !GP AT IX+1
      DO I = 1,ME
         AE( (IM-1)*ME+I,(3-1)*ME+I ) = -WEI1(3,IW)         !FP AT IX+2, THIRD FP
      END DO
      DO I=1,ME
         AE( (IM-1)*ME+I,(4-1)*ME+I ) = -WEI1(4,IW)         !FP AT IX-1 OR IX+3, FORTH FP
      END DO

      IF ( MOD(DATAL%ID,2) .EQ. 0 ) THEN       !"-" => "+", U^{-}=U^{+}-[U], F=G-[U]
         BE( (IM-1)*ME+6 ) = -1.0D0
      ELSE                                     !"+" => "-", U^{-}=U^{+}-[U], G+[U]=F
         BE( (IM-1)*ME+6 ) =  1.0D0
      END IF

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !   FIRST ORDER:  BETA^- U^-_X = BETA^+ U^+_X - [BETA U_X]
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IM = 2       !DISCRETIZE THE SECOND LINEAR EQUATION
      IW = 1       !FIRST ORDER DERIVATIVE

      IF ( MOD(DATAL%ID,2) .EQ. 0 ) THEN       !"-" => "+"
         BETAOUT = DATAL%BETA1
         BETAIN  = DATAL%BETA2
      ELSE                                     !"+" => "-"
         BETAOUT = DATAL%BETA2
         BETAIN  = DATAL%BETA1
      END IF

      !OUTSIDE INTERPOLATION
      BE( (IM-1)*ME+1 ) = -WEI2(1,IW) * BETAOUT                       !GP AT IX-1
      BE((IM-1)*ME+2)   = -WEI2(2,IW) * BETAOUT                       !GP AT IX
      DO I = 1,ME
         AE( (IM-1)*ME+I,(2-1)*ME+I ) = WEI2(3,IW)*BETAOUT            !FP AT IX+1, SECOND FP
      END DO
      BE( (IM-1)*ME+4 ) = -WEI2(4,IW)*BETAOUT                         !GP AT IX+2

      !INSIDE INTERPOLATION
      DO I = 1,ME
         AE( (IM-1)*ME+I,(1-1)*ME+I ) = -WEI1(1,IW)*BETAIN            !FP AT IX, FIRST FP
      END DO
      BE( (IM-1)*ME+3 ) = WEI1(2,IW)*BETAIN                           !GP AT IX+1
      DO I = 1,ME
         AE( (IM-1)*ME+I,(3-1)*ME+I ) = -WEI1(3,IW)*BETAIN            !FP AT IX+2, THIRD FP
      END DO
      DO I = 1,ME
         AE( (IM-1)*ME+I,(4-1)*ME+I ) = -WEI1(4,IW)*BETAIN            !FP AT IX-1 OR IX+3, FORTH FP
      END DO

      IF ( MOD(DATAL%ID,2) .EQ. 0 ) THEN       !"-" => "+", BETA^{-}U_X^{-}=BETA^{+}U_X^{+}-[BETA U_X], F=G-[BETA U_X]
         BE( (IM-1)*ME+7 ) = -1.0D0
      ELSE                                     !"+" => "-", BETA^{-}U_X^{-}=BETA^{+}U_X^{+}-[BETA U_X], G+[BETA U_X]=F
         BE( (IM-1)*ME+7 ) =  1.0D0
      END IF

      !###################### MIB FOR RIGHT INTERFACE  #########################
      CALL GETWMIB(DATAR%GAMMA,ND,WEI,NS-1,DATAR%ITYPE,DATAL%ITYPE)

      DO I = 0,ND
         DO J = 1,NS
            WEI1(J,I) = WEI(1,J-1,I)
            WEI2(J,I) = WEI(2,J-1,I)
         END DO
      END DO

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !   ZERO ORDER:   U^- = U^+ - [U]
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IM = 3       !DISCRETIZE THE THIRD LINEAR EQUATION
      IW = 0       !ZERO ORDER DERIVATIVE

      !OUTSIDE INTERPOLATION
      BE( (IM-1)*ME+2 ) = -WEI2(1,IW)                         !GP AT IX
      DO I = 1,ME
         AE( (IM-1)*ME+I,(2-1)*ME+I ) = WEI2(2,IW)            !FP AT IX+1, SECOND FP
      END DO
      BE( (IM-1)*ME+4 ) = -WEI2(3,IW)                         !GP AT IX+2
      BE( (IM-1)*ME+5 ) = -WEI2(4,IW)                         !GP AT IX+3

      !INSIDE INTERPOLATION
      DO I = 1,ME
         AE( (IM-1)*ME+I,(1-1)*ME+I ) = -WEI1(1,IW)           !FP AT IX, FIRST FP
      END DO
      BE( (IM-1)*ME+3 ) = WEI1(2,IW)                          !GP AT IX+1
      DO I = 1,ME
         AE( (IM-1)*ME+I,(3-1)*ME+I ) = -WEI1(3,IW)           !FP AT IX+2, THIRD FP
      END DO
      DO I = 1,ME
         AE( (IM-1)*ME+I,(4-1)*ME+I ) = -WEI1(4,IW)           !FP AT IX-1 OR IX+3, FORTH FP
      END DO

      IF ( MOD(DATAR%ID,2) .EQ. 0 ) THEN         !"-" => "+", U^{-}=U^{+}-[U], F=G+[U]
         BE( (IM-1)*ME+8 ) =  1.0D0
      ELSE                                       !"+" => "-", U^{-}=U^{+}-[U], G-[U]=F
         BE( (IM-1)*ME+8 ) = -1.0D0
      END IF

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !   FIRST ORDER:  BETA^- U^-_X = BETA^+ U^+_X - [BETA U_X]
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IM = 4       !DISCRETIZE THE FORTH LINEAR EQUATION
      IW = 1       !FIRST ORDER DERIVATIVE

      IF ( MOD(DATAR%ID,2) .EQ. 0 ) THEN         !"-" => "+"
         BETAOUT = DATAR%BETA2
         BETAIN  = DATAR%BETA1
      ELSE                                     !"+" => "-"
         BETAOUT = DATAR%BETA1
         BETAIN  = DATAR%BETA2
      END IF

      !OUTSIDE INTERPOLATION
      BE( (IM-1)*ME+2 ) = -WEI2(1,IW)*BETAOUT                       !GP AT IX
      DO I = 1,ME
         AE( (IM-1)*ME+I,(2-1)*ME+I ) = WEI2(2,IW)*BETAOUT          !FP AT IX+1, SECOND FP
      END DO
      BE( (IM-1)*ME+4 ) = -WEI2(3,IW)*BETAOUT                       !GP AT IX+2
      BE( (IM-1)*ME+5 ) = -WEI2(4,IW)*BETAOUT                       !GP AT IX+3

      !INSIDE INTERPOLATION
      DO I = 1,ME
         AE( (IM-1)*ME+I,(1-1)*ME+I ) = -WEI1(1,IW)*BETAIN          !FP AT IX, FIRST FP
      END DO
      BE( (IM-1)*ME+3 ) = WEI1(2,IW)*BETAIN                         !GP AT IX+1
      DO I = 1,ME
         AE( (IM-1)*ME+I,(3-1)*ME+I ) = -WEI1(3,IW)*BETAIN          !FP AT IX+2, THIRD FP
      END DO
      DO I = 1,ME
         AE( (IM-1)*ME+I,(4-1)*ME+I ) = -WEI1(4,IW)*BETAIN          !FP AT IX-1 OR IX+3, FORTH FP
      END DO

      IF ( MOD(DATAR%ID,2) .EQ. 0 ) THEN         !"-" => "+", BETA^{-}U_X^{-}=BETA^{+}U_X^{+}-[BETA U_X], F=G+[BETA U_X]
         BE( (IM-1)*ME+9 ) =  1.0D0
      ELSE                                       !"+" => "-", BETA^{-}U_X^{-}=BETA^{+}U_X^{+}-[BETA U_X], G-[BETA U_X]=F
         BE( (IM-1)*ME+9 ) = -1.0D0
      END IF

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !   LU DECOMPOSITION
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL LUDCMP(AE,NE,NE,INDEXE,DE,VVE)
      CALL LUBKSB(AE,NE,NE,INDEXE,BE)

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !          ABNORMAL ORDER TO STORE FP
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       DO I=1,M
!          DO J=1,ME
!             EMAT(I,J)=BE((I-1)*ME+J)
!             DATAL%WIJ2(I,J)=EMAT(I,J)
!             DATAR%WIJ2(I,J)=EMAT(I,J)
!          END DO
!       END DO

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !          LEFT TO RIGHT ORDER TO STORE FP
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO I = 1,M
         DO J = 1,ME
            EMAT(I,J) = BE( (I-1)*ME+J )
         END DO
      END DO
      IF ( DATAR%ITYPE .EQ. 2 ) THEN         !CLOSER TO RIGHT INTERFACE
         DO I = 1,M
            DO J = 1,ME
               DATAL%WIJ2(I,J) = EMAT(I,J)
               DATAR%WIJ2(I,J) = EMAT(I,J)
            END DO
         END DO
      ELSE                                 !CLOSER TO LEFT INTERFACE
         DATAL%WIJ2(1,J) = EMAT(4,J)
         DATAR%WIJ2(1,J) = EMAT(4,J)
         DO I = 1,M-1
            DO J = 1,ME
               DATAL%WIJ2(I+1,J) = EMAT(I,J)
               DATAR%WIJ2(I+1,J) = EMAT(I,J)
            END DO
         END DO
      END IF

      DEALLOCATE(AE,BE,VVE,INDEXE,WEI,WEI1,WEI2,STAT=IERR)

      RETURN

      END SUBROUTINE CMIB1D

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TEST_MIB1D --
      !
      ! ARGUMENTS:
      !
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TEST_MIB1D(DATA)

      USE MOD_DATA

      TYPE(LIST_DATA) :: DATA

      REAL :: UH(NY,NX),FPUH(NY,NX),TSUM(2)
      REAL :: T
      INTEGER :: I,J

      T=0.0D0
      CALL ANALYTICAL(UH,T)              !INITIALIZE ANALITICAL SOLUTIONS
      CALL FPANALYTICAL(FPUH,T)          !INITIALIZE FP ANALITYCAL SOLUTIONS

      IF (DATA%AXTP .EQ. 0) THEN         !IT'S ON X-AXIS
         TSUM=0.0D0
         DO J=1,2
            DO I=-1,2
               TSUM(J)=TSUM(J)+UH(DATA%LGRD+I,DATA%AXID)*DATA%WIJ(J,I+2)
            END DO
            TSUM(J)=TSUM(J)+DATA%WIJ(J,5)*DATA%JUMP(1)
            TSUM(J)=TSUM(J)+DATA%WIJ(J,6)*DATA%JUMP(4)
         END DO

         DATA%ERR(1)=ABS((FPUH(DATA%LGRD,DATA%AXID)-TSUM(1))/FPUH(DATA%LGRD,DATA%AXID))              !THE CASE ITYPE==0
         DATA%ERR(2)=ABS((FPUH(DATA%LGRD+1,DATA%AXID)-TSUM(2))/FPUH(DATA%LGRD+1,DATA%AXID))
         IF (DATA%ITYPE .EQ. -1) THEN
            DATA%ERR(1)=ABS((FPUH(DATA%LGRD-1,DATA%AXID)-TSUM(1))/FPUH(DATA%LGRD-1,DATA%AXID))
         ELSE IF (DATA%ITYPE .EQ. 1) THEN
            DATA%ERR(2)=ABS((FPUH(DATA%LGRD+2,DATA%AXID)-TSUM(2))/FPUH(DATA%LGRD+2,DATA%AXID))
         END IF

      ELSE IF (DATA%AXTP .EQ. 1) THEN     !IT'S ON Y-AXIS
         TSUM=0.0D0
         DO J=1,2
            DO I=-1,2
               TSUM(J)=TSUM(J)+UH(DATA%AXID,DATA%LGRD+I)*DATA%WIJ(J,I+2)
            END DO
            TSUM(J)=TSUM(J)+DATA%WIJ(J,5)*DATA%JUMP(1)
            TSUM(J)=TSUM(J)+DATA%WIJ(J,6)*DATA%JUMP(4)
         END DO

         DATA%ERR(1)=ABS((FPUH(DATA%AXID,DATA%LGRD)-TSUM(1))/FPUH(DATA%AXID,DATA%LGRD))               !THE CASE ITYPE==0
         DATA%ERR(2)=ABS((FPUH(DATA%AXID,DATA%LGRD+1)-TSUM(2))/FPUH(DATA%AXID,DATA%LGRD+1))
         IF (DATA%ITYPE .EQ. -1) THEN
            DATA%ERR(1)=ABS((FPUH(DATA%AXID,DATA%LGRD-1)-TSUM(1))/FPUH(DATA%AXID,DATA%LGRD-1))
         ELSE IF (DATA%ITYPE .EQ. 1) THEN
            DATA%ERR(2)=ABS((FPUH(DATA%AXID,DATA%LGRD+2)-TSUM(2))/FPUH(DATA%AXID,DATA%LGRD+2))
         END IF

      ELSE
         WRITE(*,*) "ERROR IN MIB1D_TEST: NO AXTP MATCHED"
         STOP
      END IF

      RETURN

      END SUBROUTINE TEST_MIB1D

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TEST_CMIB1D --
      !
      ! ARGUMENTS:
      !
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TEST_CMIB1D(DATAL,DATAR)

      USE MOD_DATA

      TYPE(LIST_DATA) :: DATAL,DATAR

      REAL :: UH(NY,NX),FPUH(NY,NX),TSUM(4)
      REAL :: T
      INTEGER :: I,J

      T = 0.0D0
      CALL ANALYTICAL(UH,T)              !INITIALIZE ANALITICAL SOLUTIONS
      CALL FPANALYTICAL(FPUH,T)          !INITIALIZE FP ANALITYCAL SOLUTIONS

      IF (DATAL%AXTP .EQ. 0) THEN        !IT'S ON X-AXIS
         TSUM = 0.0D0
         DO J = 1,4
            DO I = -1,3
               TSUM(J) = TSUM(J) + UH(DATAL%LGRD+I,DATAL%AXID) * DATAL%WIJ2(J,I+2)
            END DO
            TSUM(J) = TSUM(J) + DATAL%WIJ2(J,6) * DATAL%JUMP(1) + DATAL%WIJ2(J,7) * DATAL%JUMP(4) + &
                      DATAL%WIJ2(J,8) * DATAR%JUMP(1) + DATAL%WIJ2(J,9) * DATAR%JUMP(4)
         END DO

         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !          ABNORMAL ORDER TO STORE FP
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         DO I = 1,4                                 !(0,2)
!            DATAL%ERR2(I) = ABS( ( FPUH(DATAL%LGRD-1+I,DATAL%AXID) - TSUM(I) ) / FPUH(DATAL%LGRD-1+I,DATAL%AXID) )
!         END DO
!         IF (DATAL%ITYPE .EQ. -2) THEN             !(-2,0), FP4 CHANGES
!             DATAL%ERR2(4) = ABS( ( FPUH(DATAL%LGRD-1,DATAL%AXID) - TSUM(4) ) / FPUH(DATAL%LGRD-1,DATAL%AXID) )
!             IF (DATAR%ITYPE .EQ. 1) THEN          !(-2,1), FP3&FP4 BOTH CHANGES
!                DATAL%ERR2(3) = ABS( ( FPUH(DATAL%LGRD+3,DATAL%AXID) - TSUM(3) ) / FPUH(DATAL%LGRD+3,DATAL%AXID) )
!             END IF
!          ELSE
!             IF (DATAR%ITYPE .EQ. -1) THEN         !(-1,2), FP1 CHANGES
!                DATAL%ERR2(1) = ABS( ( FPUH(DATAL%LGRD-1,DATAL%AXID) - TSUM(1) ) / FPUH(DATAL%LGRD-1,DATAL%AXID) )
!             END IF
!          END IF

         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !          LEFT TO RIGHT ORDER TO STORE FP
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF (DATAL%ITYPE .EQ. -2) THEN      !LEFT INTERFACE TYPE
            DO I = 1,4
               DATAL%ERR2(I) = ABS( ( FPUH(DATAL%LGRD-2+I,DATAL%AXID)-TSUM(I) ) / FPUH(DATAL%LGRD-2+I,DATAL%AXID) )
            END DO
            IF (DATAR%ITYPE .EQ. 1) THEN    !FP4 GOES ONE GRID NODE LEFT
               DATAL%ERR2(4) = ABS( ( FPUH(DATAL%LGRD+3,DATAL%AXID) - TSUM(4) ) / FPUH(DATAL%LGRD+3,DATAL%AXID) )
            END IF
         ELSE                               !RIGHT ITERFACE TYPE
            DO I = 1,4
               DATAL%ERR2(I) = ABS( ( FPUH(DATAL%LGRD-1+I,DATAL%AXID) - TSUM(I) ) / FPUH(DATAL%LGRD-1+I,DATAL%AXID) )
            END DO
            IF (DATAL%ITYPE .EQ. -1) THEN   !FP4 GOES ONE GRID NODE LEFT
               DATAL%ERR2(1) = ABS( ( FPUH(DATAL%LGRD-1,DATAL%AXID) - TSUM(1) ) / FPUH(DATAL%LGRD-1,DATAL%AXID) )
            END IF
         END IF


      ELSE IF (DATAL%AXTP .EQ. 1) THEN    !IT'S ON Y-AXIS
         TSUM = 0.0D0
         DO J = 1,4
            DO I = -1,3
               TSUM(J) = TSUM(J) + UH(DATAL%AXID,DATAL%LGRD+I) * DATAL%WIJ2(J,I+2)
            END DO
            TSUM(J) = TSUM(J) + DATAL%WIJ2(J,6) * DATAL%JUMP(1) + DATAL%WIJ2(J,7) * DATAL%JUMP(4) + &
                      DATAL%WIJ2(J,8) * DATAR%JUMP(1) + DATAL%WIJ2(J,9) * DATAR%JUMP(4)
         END DO

         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !          ABNORMAL ORDER TO STORE FP
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!          DO I = 1,4                                 !(0,2)
!             DATAL%ERR2(I) = ABS( ( FPUH(DATAL%AXID,DATAL%LGRD-1+I)-TSUM(I) ) / FPUH(DATAL%AXID,DATAL%LGRD-1+I) )
!          END DO
!          IF (DATAL%ITYPE .EQ. -2) THEN              !(-2,0), FP4 CHANGES
!             DATAL%ERR2(4) = ABS( ( FPUH(DATAL%AXID,DATAL%LGRD-1) - TSUM(4) ) / FPUH(DATAL%AXID,DATAL%LGRD-1) )
!             IF (DATAR%ITYPE .EQ. 1) THEN            !(-2,1), FP3&FP4 BOTH CHANGES
!                DATAL%ERR2(3) = ABS( ( FPUH(DATAL%AXID,DATAL%LGRD+3) - TSUM(3) ) / FPUH(DATAL%AXID,DATAL%LGRD+3) )
!             END IF
!          ELSE
!             IF (DATAR%ITYPE .EQ. -1) THEN           !(-1,2), FP1 CHANGES
!                DATAL%ERR2(1) = ABS( ( FPUH(DATAL%AXID,DATAL%LGRD-1)-TSUM(1) ) / FPUH(DATAL%AXID,DATAL%LGRD-1) )
!             END IF
!          END IF

         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !          LEFT TO RIGHT ORDER TO STORE FP
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF (DATAL%ITYPE .EQ. -2) THEN       !LEFT INTERFACE TYPE
            DO I = 1,4
               DATAL%ERR2(I) = ABS( ( FPUH(DATAL%AXID,DATAL%LGRD-2+I) - TSUM(I) ) / FPUH(DATAL%AXID,DATAL%LGRD-2+I) )
            END DO
            IF (DATAR%ITYPE .EQ. 1) THEN     !FP4 GOES ONE GRID NODE LEFT
               DATAL%ERR2(4) = ABS( ( FPUH(DATAL%AXID,DATAL%LGRD+3) - TSUM(4) ) / FPUH(DATAL%AXID,DATAL%LGRD+3) )
            END IF
         ELSE                                !RIGHT ITERFACE TYPE
            DO I = 1,4
               DATAL%ERR2(I) = ABS( ( FPUH(DATAL%AXID,DATAL%LGRD-1+I) - TSUM(I) ) / FPUH(DATAL%AXID,DATAL%LGRD-1+I) )
            END DO
            IF (DATAL%ITYPE .EQ. -1) THEN    !FP4 GOES ONE GRID NODE LEFT
               DATAL%ERR2(1) = ABS( ( FPUH(DATAL%AXID,DATAL%LGRD-1) - TSUM(1) ) / FPUH(DATAL%AXID,DATAL%LGRD-1) )
            END IF
         END IF

      ELSE
         WRITE(*,*) "ERROR IN MIB1D_TEST: NO AXTP MATCHED"
         STOP
      END IF

      DO I = 1,4
         DATAR%ERR2(I) = DATAL%ERR2(I)
      END DO

      RETURN

      END SUBROUTINE TEST_CMIB1D
