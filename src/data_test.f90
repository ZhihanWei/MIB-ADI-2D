      !
      !------------- SUBROUTINES USED FOR DEBUGGING PURPOSE ONLY ---------------
      !

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TEST_ALL --
      !    PRINT OUT ALL VARIABLES CONTAINED IN THE DATA MODULE
      ! ARGUMENTS:
      !    TNOW   IN   CURRENT TIME
      ! NOTES:
      !    TEST FUNCTION ONLY USED FOR DEBUGGING PURPOSE
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TEST_ALL(TNOW)

      REAL :: TNOW

      INTEGER :: IX,IY
      TYPE(LINKED_LIST), POINTER  :: LIST

      OPEN(100,FILE="test_all.txt")

      WRITE(100,"(2(A,I4),3(A,F12.6))") "NX = ",NX,", NY = ",NY,", DX = ",DX,", XL = ",XL,", YL = ",YL
      WRITE(100,"(4(A,F12.6))") "TSTART = ",TSTART,", TFINAL = ",TFINAL,", TSTEP = ",TSTEP,", TNOW = ",TNOW

      WRITE(100,"(A)") "INODE(NY,NX):"
      DO IY=1,NY
         DO IX=1,NX
            WRITE(100,"(3(A,I4))"), "   INODE(",IY,",",IX,") = ",INODE(IY,IX)
         END DO
      END DO

      WRITE(100,"(A)") "BETA(NY,NX):"
      DO IY=1,NY
         DO IX=1,NX
            WRITE(100,"(2(A,I4),A,F12.6)"), "   BETA(",IY,",",IX,") = ",BETA(IY,IX)
         END DO
      END DO

      WRITE(100,"(A)") "BETAX(NY,NX):"
      DO IY=1,NY
         DO IX=1,NX
            WRITE(100,"(2(A,I4),A,F12.6)"), "   BETAX(",IY,",",IX,") = ",BETAX(IY,IX)
         END DO
      END DO

      WRITE(100,"(A)") "BETAY(NY,NX):"
      DO IY=1,NY
         DO IX=1,NX
            WRITE(100,"(2(A,I4),A,F12.6)"), "   BETAY(",IY,",",IX,") = ",BETAY(IY,IX)
         END DO
      END DO

      !-------------------------------------------------------------------------

      WRITE(100,"(A)") "IFPY(NY):"
      DO IY=1,NY
         IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN
            WRITE(100,"(A,I4,A)") " -------------------- IY = ",IY," --------------------"

            LIST => IFPY(IY)%HEAD
            DO WHILE ( ASSOCIATED(LIST) )
               WRITE(100,"(4(A,I4),5(A,F12.6),(A,F12.10),(A,I4))") " AXTP = ",LIST%DATA%AXTP,", AXID = ",LIST%DATA%AXID,&
               &", ID = ",LIST%DATA%ID,", LOC = ",LIST%DATA%LGRD,", THETA = ",LIST%DATA%THETA,&
               &", BETA1 = ",LIST%DATA%BETA1,", BETA2 = ",LIST%DATA%BETA2,", X = ",LIST%DATA%X,", Y = ",LIST%DATA%Y,&
               &", GAMMA = ",LIST%DATA%GAMMA,", ITYPE = ",LIST%DATA%ITYPE
               !WRITE(100,"(A,3(I4,A),5F12.6,A)") " NEAREST 3 NODES: (",&
               !&LIST%DATA%ITAUM(1,1),",",LIST%DATA%ITAUM(2,1),",",LIST%DATA%ITAUM(3,1),"), JUMPS = (",LIST%DATA%JUMP,")"

               WRITE(100,*) " "

               LIST => LIST%NEXT
            END DO
         END IF
      END DO

      WRITE(100,"(A)") "IFPX(NX):"
      DO IX=1,NX
         IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN
            WRITE(100,"(A,I4,A)") " -------------------- IX = ",IX," --------------------"

            LIST => IFPX(IX)%HEAD
            DO WHILE ( ASSOCIATED(LIST) )
               WRITE(100,"(4(A,I4),5(A,F12.6),(A,F12.10),(A,I4))") " AXTP = ",LIST%DATA%AXTP,", AXID = ",LIST%DATA%AXID,&
               &", ID = ",LIST%DATA%ID,", LOC = ",LIST%DATA%LGRD,", THETA = ",LIST%DATA%THETA,&
               &", BETA1 = ",LIST%DATA%BETA1,", BETA2 = ",LIST%DATA%BETA2,", X = ",LIST%DATA%X,", Y = ",LIST%DATA%Y,&
               &", GAMMA = ",LIST%DATA%GAMMA,", ITYPE = ",LIST%DATA%ITYPE
               !WRITE(100,"(A,3(I4,A),5F12.6,A)") " NEAREST 3 NODES: (",&
               !&LIST%DATA%ITAUM(1,2),",",LIST%DATA%ITAUM(2,2),",",LIST%DATA%ITAUM(3,2),"), JUMPS = (",LIST%DATA%JUMP,")"

               WRITE(100,*) " "

               LIST => LIST%NEXT
            END DO
         END IF
      END DO

      CLOSE(100)

      RETURN

      END SUBROUTINE TEST_ALL

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TEST_ALLIFPS --
      !    PRINT OUT ALL INTERFACE POINTS
      ! ARGUMENTS:
      !    TNOW   IN   CURRENT TIME
      ! NOTES:
      !    TEST FUNCTION ONLY USED FOR DEBUGGING PURPOSE
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TEST_ALLIFPS(TNOW)

      REAL :: TNOW

      INTEGER :: IX,IY,I
      TYPE(LINKED_LIST), POINTER  :: LIST

      OPEN(100,FILE="test_allifps.txt")

      WRITE(100,"(A,F12.6)") " TNOW = ",TNOW
      WRITE(100,*) "  "

      WRITE(100,"(A,I4)") "TOTAL # OF IPTS ON IX GRID LINES = ", NIPXS
      WRITE(100,"(A,I4)") "TOTAL # OF IPTS ON IY GRID LINES = ", NIPYS
      WRITE(100,"(A,I4)") "TOTAL # OF CORNER PTS            = ", NCORNERS
      WRITE(100,"(A,I4)") "TOTAL # OF ALT. AUX. PTS         = ", NALTAUXS
      WRITE(100,*) "  "

      !----- INTERFACE PTS ON EACH IX GRID LINE
      WRITE(100,"(A)") "IFPX(NX):"
      DO IX=1,NX
         IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN

            LIST => IFPX(IX)%HEAD

            WRITE(100,"(2(A,I4),A)") " -------------------- IX = ",IX,", NUMBER OF INTERFACE PTS = ",LIST_COUNT(LIST),&
            &" --------------------"

            DO WHILE ( ASSOCIATED(LIST) )
               WRITE(100,"(4(A13,I4),2(A14,F12.6))") "AXTP = ",LIST%DATA%AXTP,", AXID = ",LIST%DATA%AXID,&
               &", ID = ",LIST%DATA%ID,", LEFT GRID = ",LIST%DATA%LGRD,", X = ",LIST%DATA%X,", Y = ",LIST%DATA%Y

               WRITE(100,"(4(A13,F12.6))") "THETA = ",LIST%DATA%THETA,&
               &", BETA1 = ",LIST%DATA%BETA1,", BETA2 = ",LIST%DATA%BETA2,", GAMMA = ",LIST%DATA%GAMMA

               WRITE(100,"(A13,I4)") "ITYPE = ",LIST%DATA%ITYPE
               DO I=1,2
                  WRITE(100,"(A14,6F12.6,A2)") "WIJ = (", LIST%DATA%WIJ(I,:),")"
               END DO
               DO I=1,4
                  WRITE(100,"(A14,9F12.6,A2)") "WIJ2 = (", LIST%DATA%WIJ2(I,:),")"
               END DO
               WRITE(100,"(A14,2F12.6,A2)") "ERR = (", LIST%DATA%ERR,")"
               WRITE(100,"(A14,4F12.6,A2)") "ERR2 = (", LIST%DATA%ERR2,")"

               WRITE(100,"(A13,I4)") "ITAU = ",LIST%DATA%ITAU
               WRITE(100,"(2(A13,I4),2(A,F12.6),(A,I4),(A,F12.6))") "AUXL_AXTP = ",LIST%DATA%AUXL_AXTP,&
               &", AUXL_AXID = ",LIST%DATA%AUXL_AXID,", AUXL_X = ",LIST%DATA%AUXL_X,", AUXL_Y = ",LIST%DATA%AUXL_Y,&
               &", AUXL = ",LIST%DATA%AUXL,", DAUXL = ",LIST%DATA%DAUXL
               WRITE(100,"(A14,3I4,A13,3F12.6,A2)") "IAUXL = (", LIST%DATA%IAUXL(:),&
               &"), WAUXL = (",LIST%DATA%WAUXL(:),")"

               WRITE(100,"(2(A13,I4),2(A,F12.6),(A,I4),(A,F12.6))") "AUXR_AXTP = ",LIST%DATA%AUXR_AXTP,&
               &", AUXR_AXID = ",LIST%DATA%AUXR_AXID,", AUXR_X = ",LIST%DATA%AUXR_X,", AUXR_Y = ",LIST%DATA%AUXR_Y,&
               &", AUXR = ",LIST%DATA%AUXR,", DAUXR = ",LIST%DATA%DAUXR
               WRITE(100,"(A14,3I4,A13,3F12.6,A2)") "IAUXR = (", LIST%DATA%IAUXR(:),&
               &"), WAUXR = (",LIST%DATA%WAUXR(:),")"

                WRITE(100,"(A14,6F12.6,A2)") "JUMP = (", LIST%DATA%JUMP(:),")"
               WRITE(100,"(A14,4F12.6,A2)") "UTAU = (", LIST%DATA%UTAU(:),")"
               WRITE(100,"(A14,4F12.6,A2)") "UAUXL = (", LIST%DATA%UAUXL(:),")"
               WRITE(100,"(A14,4F12.6,A2)") "UAUXR = (", LIST%DATA%UAUXR(:),")"

               WRITE(100,*) "  "
               WRITE(100,*) "  "

               LIST => LIST%NEXT
            END DO
         END IF
      END DO

      !----- INTERFACE PTS ON EACH IY GRID LINE
      WRITE(100,"(A)") "IFPY(NY):"
      DO IY=1,NY
         IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN

            LIST => IFPY(IY)%HEAD

            WRITE(100,"(2(A,I4),A)") " -------------------- IY = ",IY,", NUMBER OF INTERFACE PTS = ",LIST_COUNT(LIST),&
            &" --------------------"

            DO WHILE ( ASSOCIATED(LIST) )
               WRITE(100,"(4(A13,I4),2(A14,F12.6))") "AXTP = ",LIST%DATA%AXTP,", AXID = ",LIST%DATA%AXID,&
               &", ID = ",LIST%DATA%ID,", LEFT GRID = ",LIST%DATA%LGRD,", X = ",LIST%DATA%X,", Y = ",LIST%DATA%Y

               WRITE(100,"(4(A13,F12.6))") "THETA = ",LIST%DATA%THETA,&
               &", BETA1 = ",LIST%DATA%BETA1,", BETA2 = ",LIST%DATA%BETA2,", GAMMA = ",LIST%DATA%GAMMA

               WRITE(100,"(A13,I4)") "ITYPE = ",LIST%DATA%ITYPE
               DO I=1,2
                  WRITE(100,"(A14,6F12.6,A2)") "WIJ = (", LIST%DATA%WIJ(I,:),")"
               END DO
               DO I=1,4
                  WRITE(100,"(A14,9F12.6,A2)") "WIJ2 = (", LIST%DATA%WIJ2(I,:),")"
               END DO
               WRITE(100,"(A14,2F12.6,A2)") "ERR = (", LIST%DATA%ERR,")"
               WRITE(100,"(A14,4F12.6,A2)") "ERR2 = (", LIST%DATA%ERR2,")"

               WRITE(100,"(A13,I4)") "ITAU = ",LIST%DATA%ITAU
               WRITE(100,"(2(A13,I4),2(A,F12.6),(A,I4),(A,F12.6))") "AUXL_AXTP = ",LIST%DATA%AUXL_AXTP,&
               &", AUXL_AXID = ",LIST%DATA%AUXL_AXID,", AUXL_X = ",LIST%DATA%AUXL_X,", AUXL_Y = ",LIST%DATA%AUXL_Y,&
               &", AUXL = ",LIST%DATA%AUXL,", DAUXL = ",LIST%DATA%DAUXL
               WRITE(100,"(A14,3I4,A13,3F12.6,A2)") "IAUXL = (", LIST%DATA%IAUXL(:),&
               &"), WAUXL = (",LIST%DATA%WAUXL(:),")"

               WRITE(100,"(2(A13,I4),2(A,F12.6),(A,I4),(A,F12.6))") "AUXR_AXTP = ",LIST%DATA%AUXR_AXTP,&
               &", AUXR_AXID = ",LIST%DATA%AUXR_AXID,", AUXR_X = ",LIST%DATA%AUXR_X,", AUXR_Y = ",LIST%DATA%AUXR_Y,&
               &", AUXR = ",LIST%DATA%AUXR,", DAUXR = ",LIST%DATA%DAUXR
               WRITE(100,"(A14,3I4,A13,3F12.6,A2)") "IAUXR = (", LIST%DATA%IAUXR(:),&
               &"), WAUXR = (",LIST%DATA%WAUXR(:),")"

               WRITE(100,"(A14,6F12.6,A2)") "JUMP = (", LIST%DATA%JUMP(:),")"
               WRITE(100,"(A14,4F12.6,A2)") "UTAU = (", LIST%DATA%UTAU(:),")"
               WRITE(100,"(A14,4F12.6,A2)") "UAUXL = (", LIST%DATA%UAUXL(:),")"
               WRITE(100,"(A14,4F12.6,A2)") "UAUXR = (", LIST%DATA%UAUXR(:),")"

               WRITE(100,*) "  "
               WRITE(100,*) "  "

               LIST => LIST%NEXT
            END DO
         END IF
      END DO

      CLOSE(100)

      RETURN

      END SUBROUTINE TEST_ALLIFPS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TEST_RHS --
      !
      ! ARGUMENTS:
      !
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TEST_U(U)

      REAL :: U(NY,NX)
      INTEGER :: IX,IY

      OPEN(100,FILE="test_rhs.txt")

      DO IY = 1,NY
         DO IX = 1,NX
            WRITE(100,"(2(A,I4),A,F12.6)") " RHS ( ",IY,", ",IX," ) = ",U(IY,IX)
         END DO
      END DO

      RETURN

      END SUBROUTINE TEST_U

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TEST_JUMPS --
      !
      ! ARGUMENTS:
      !
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TEST_JUMPS

      INTEGER :: IX,IY
      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LIST_DATA) :: DATA
      REAL :: MAXERR_JUMP = 0.0D0
      TYPE(LIST_DATA) :: MAXERR_DATA

      !----- RUNNING OVER EACH IX GRID LINE
      DO IX=1,NX
         IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN
            ELEM => IFPX(IX)%HEAD
            DO WHILE ( ASSOCIATED(ELEM) )
               DATA = LIST_GET_DATA(ELEM)

               IF ( ABS( DATA%JUMP(5) - DATA%JUMP(6) ) .GT. MAXERR_JUMP) THEN
                  MAXERR_JUMP = ABS( DATA%JUMP(5) - DATA%JUMP(6) )
                  MAXERR_DATA = DATA
               END IF

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

               !+++++ TEST THE DIFF. BETWEEN ANA. AND NUM. CALCULATED U_TAU
               IF ( ABS( DATA%JUMP(5) - DATA%JUMP(6) ) .GT. MAXERR_JUMP) THEN
                  MAXERR_JUMP = ABS( DATA%JUMP(5) - DATA%JUMP(6) )
                  MAXERR_DATA = DATA
               END IF

               ELEM => ELEM%NEXT ! NEXT INTERFACE
            END DO
         END IF
      END DO

      !+++++ TEST THE DIFF. BETWEEN ANA. AND NUM. CALCULATED U_TAU
      WRITE(*,"(A,F12.6,F12.6,A,F12.6)") "MAX JUMP ERROR OCCURS AT IFP ( ",MAXERR_DATA%X,MAXERR_DATA%Y,&
                                       & " ) W/ ERROR = ",MAXERR_JUMP
      CALL TEST_IFP(IFPX,NX,22,1)

      RETURN

      END SUBROUTINE TEST_JUMPS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TEST_IFP --
      !
      ! ARGUMENTS:
      !
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TEST_IFP(LIST,DIM,AXID,PID)

      TYPE(LIST_ARRAY),DIMENSION(DIM) :: LIST
      INTEGER :: DIM,AXID,PID

      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LIST_DATA) :: DATA
      INTEGER :: I

      ELEM => LIST(AXID)%HEAD
      I = 1
      DO WHILE (I .LT. PID)
         ELEM => ELEM%NEXT
         I = I + 1
      END DO
      DATA=LIST_GET_DATA(ELEM)

      WRITE(*,"(4(A13,I4),2(A14,F12.6))") "AXTP = ",DATA%AXTP,", AXID = ",DATA%AXID,&
      &", ID = ",DATA%ID,", LEFT GRID = ",DATA%LGRD,", X = ",DATA%X,", Y = ",DATA%Y

      WRITE(*,"(4(A13,F12.6))") "THETA = ",DATA%THETA,&
      &", BETA1 = ",DATA%BETA1,", BETA2 = ",DATA%BETA2,", GAMMA = ",DATA%GAMMA

      WRITE(*,"(A13,I4)") "ITAU = ",DATA%ITAU
      WRITE(*,"(2(A13,I4),2(A,F12.6),(A,I4),(A,F12.6))") "AUXL_AXTP = ",DATA%AUXL_AXTP,&
      &", AUXL_AXID = ",DATA%AUXL_AXID,", AUXL_X = ",DATA%AUXL_X,", AUXL_Y = ",DATA%AUXL_Y,&
      &", AUXL = ",DATA%AUXL,", DAUXL = ",DATA%DAUXL
      WRITE(*,"(A14,3I4,A13,3F12.6,A2)") "IAUXL = (", DATA%IAUXL(:),&
      &"), WAUXL = (",DATA%WAUXL(:),")"

      WRITE(*,"(2(A13,I4),2(A,F12.6),(A,I4),(A,F12.6))") "AUXR_AXTP = ",DATA%AUXR_AXTP,&
      &", AUXR_AXID = ",DATA%AUXR_AXID,", AUXR_X = ",DATA%AUXR_X,", AUXR_Y = ",DATA%AUXR_Y,&
      &", AUXR = ",DATA%AUXR,", DAUXR = ",DATA%DAUXR
      WRITE(*,"(A14,3I4,A13,3F12.6,A2)") "IAUXR = (",DATA%IAUXR(:),&
      &"), WAUXR = (",DATA%WAUXR(:),")"

      WRITE(*,"(A14,6F12.6,A2)") "JUMP = (", DATA%JUMP(:),")"
      WRITE(*,"(A14,4F12.6,A2)") "UTAU = (", DATA%UTAU(:),")"
      WRITE(*,"(A14,4F12.6,A2)") "UAUXL = (",DATA%UAUXL(:),")"
      WRITE(*,"(A14,4F12.6,A2)") "UAUXR = (",DATA%UAUXR(:),")"

      WRITE(*,*) "  "
      WRITE(*,*) "  "


      RETURN

      END SUBROUTINE TEST_IFP
