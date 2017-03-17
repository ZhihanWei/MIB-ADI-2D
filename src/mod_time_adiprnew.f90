      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !       PEACEMAN RACHAFORD ADI
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ADIPR(T,DT,UH)

      USE MOD_DATA

      REAL :: T,DT,UH(NY,NX)

      REAL :: UHS(NY,NX)
      REAL, ALLOCATABLE :: A(:),B(:),C(:),R(:),UT(:)
      INTEGER :: IX,IY,ierr

      !---------- UPDATE JUMPS TO PRESENT TIME STEP T^N
      CALL SETJUMPS(T,DT,UH)
      !+++++ TEST: INTERFACE PTS
      !IF (IPRINT .EQ. 1 .AND. ILOOP .EQ. 1) THEN
      !   CALL TEST_ALLIFPS(T)
      !END IF

      !---------- SET SOURCE SRC(IY,IX) TO TIME T_{N+1/2}
      CALL SETSRC(T+DT/2.0D0)

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !
      !PEACEMAN RACHFORD ADI FIRST STEP:
      !   (1/beta - dt/2 D_xx) u* = (1/beta + dt/2 D_yy)u^n + dt/2 f^{n+1/2}/beta
      !                           = dt/2 D_yy u^n + (u^n + dt/2 f^{n+1/2})/beta
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !----- STEP I: CALCULATE THE RHS -----!
      CALL D_YY(UH,UHS) ! D_YY U^N

      DO IX=2,NX-1
         DO IY=2,NY-1
            UHS(IY,IX) = DT/2.0D0 * UHS(IY,IX) ! DT/2 * D_YY U^N
            UHS(IY,IX) = UHS(IY,IX) + ( UH(IY,IX) + DT/2.0D0 * SRC(IY,IX) ) / BETA(IY,IX) ! + (U^N + DT/2 SRC^{N+1/2})/BETA
         END DO
      END DO

      !..... NOW UHS IS THE RHS OF THE FIRST STEP FORMULA .....!

      !----- STEP II: RESTORE FIRST ORDER BOUNDARY CONDITION FOR TIME T^{N+1/2}-----!
      CALL SETBC_ADIPR(T,DT,UHS)

      !----- STEP III: CALCULATE THE LHS MATRIX -----!
      ALLOCATE(A(NX),B(NX),C(NX),R(NX),UT(NX),STAT=ierr)

      DO IY = 2,NY-1
         CALL D_XX(IY,DT/2.0D0,UHS,A,B,C,R)  ! CONVERT LHS OPERATOR (1/BETA - DT/2 D_XX) TO A TRI-DIAGONAL MATRIX
                                             ! AND STORE IT IN A, B, AND C WITH CORRESPONDING RHS R
         
         CALL TRIDAG(A,B,C,R,UT,NX) ! THOMAS ALGORITHM

         DO IX = 1,NX
            UHS(IY,IX) = UT(IX)
         END DO       

      END DO

      !+++++ TEST: ADI LHS
      !CALL TEST_U(UHS); STOP

      DEALLOCATE(A,B,C,R,UT,STAT=ierr)

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !
      !PEACEMAN RACHFORD ADI SECOND STEP:
      !   (1/beta - dt/2 D_yy) u^(n+1) = (1/beta + dt/2 D_xx) u* + dt/2 f^{n+1/2}/beta
      !                                = dt/2 D_xx u* + (u* + dt/2 f^{n+1/2})/beta
      !
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !CALL SETJUMPS(T,DT,UHS)

      !----- STEP I: CALCULATE THE RHS -----!
      CALL D_XX(UHS,UH) ! D_XX U*

      DO IY = 2,NY-1
         DO IX = 2,NX-1
            UH(IY,IX) = DT/2.0D0 * UH(IY,IX) ! DT/2 * D_XX U*
            UH(IY,IX) = UH(IY,IX) + ( UHS(IY,IX) + DT/2.0D0 * SRC(IY,IX) ) / BETA(IY,IX) ! + (U* + DT/2 SRC^{N+1/2})/BETA
         END DO
      END DO

      !..... NOW UH IS THE RHS OF THE FIRST STEP FORMULA .....!

      !----- STEP II: RESTORE FIRST ORDER BOUNDARY CONDITION FOR U^{N+1}-----!
      CALL SETBC_ADIPR(T+DT,0.0D0,UH)

      !----- STEP III: CALCULATE THE LHS MATRIX -----!
      ALLOCATE(A(NY),B(NY),C(NY),R(NY),UT(NY),STAT=ierr)

      DO IX = 2,NX-1

         CALL D_YY(IX,DT/2.0D0,UH,A,B,C,R)   ! CONVERT LHS OPERATOR (1/BETA - DT/2 D_YY) TO A TRI-DIAGONAL MATRIX
                                             ! AND STORE IT IN A, B, AND C WITH CORRESPONDING RHS R
               
         CALL TRIDAG(A,B,C,R,UT,NY) ! THOMAS ALGORITHM

         DO IY = 1,NY
            UH(IY,IX) = UT(IY)
         END DO      

      END DO

      DEALLOCATE(A,B,C,R,UT,STAT=ierr)

      RETURN

      END SUBROUTINE ADIPR

