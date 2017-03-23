      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DXX_ADIL --
      !    SECOND ORDER DERIVATIVE D_XX WHEN APPEARING ON THE LHS OF
      !    ADI FORMULATION
      ! ARGUMENTS:
      !    IY     IN     A Y GRID LINE
      !    COEF   IN     COEFFICIENT IN FRONT OF D_XX
      !    RHS    IN     RHS OF ADI FORMULATION AT CURRENT TIME STEP
      !    A      OUT    LOWER DIAGONAL ELEMENTS
      !    B      OUT    DIAGONAL ELEMENTS
      !    C      OUT    UPPER DIAGONAL ELEMENTS
      !    R      OUT    CORRESPONDING RHS VECTOR
      ! NOTES:
      !    THE SUBROUTINE CONVERTS THE MATRIX RESULTING FROM DISCRETIZING
      !    THE OPERATOR
      !         1/BETA - COEF * D_XX
      !    ON THE LEFT HAND SIDE TO A TRI-DIAGONAL MATRIX
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE DXX_ADIL(IY,COEF,RHS,A,B,C,R)

      USE MOD_DATA

      INTEGER :: IY ! A Y GRID LINE
      REAL :: COEF,RHS(NY,NX),A(NX),B(NX),C(NX),R(NX)

      INTEGER,PARAMETER :: ORDER = 2 ! 2ND DERIVATIVE D_XX
      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LIST_DATA) :: DATA,DATAL,DATAR
      REAL :: VDE2(-1:2,-1:1) ! WEIGHTS FOR FD FORMULA TO APPROX. D_XX
                              ! VDE2(-1,:) - LEFT BIASED 3-STENCIL DIFF. FORMULA
                              ! VDE2( 0,:) - CENTRAL DIFF. FORMULA WITH DX
                              ! VDE2( 1,:) - RIGHT BIASED 3-STENCIL DIFF. FORMULA
                              ! VDE2( 2,:) - CENTRAL DIFF. FORMULA WITH 2*DX
      INTEGER :: I,J,K,IX
      REAL :: UJP,UJPL,UJPR,BUXJP,BUXJPR,BUXJPL
      REAL :: ROW(2,5),ROW2(3,6)

      !---------- GET WEIGHTS ASSOICATED WITH EACH TYPE OF 3-STENCIL FINITE DIFF. FORMULAS
      CALL SETVDE2(ORDER,DX,VDE2)

      !----- STEP I: GENERATE LHS USING CENTRAL FD MATRIX FOR D_XX WITHOUT MIB -----!
      I = 0 ! CENTRAL DIFF. FORMULA
      A(1) = 0.0D0; B(1) = 1.0D0; C(1) = 0.0D0
      DO IX = 2,NX-1
         A(IX) = - COEF * VDE2(I,-1)
         B(IX) =   1.0D0/BETA(IY,IX) - COEF * VDE2(I, 0)
         C(IX) = - COEF * VDE2(I, 1)
      END DO
      A(NX) = 0.0D0; B(NX) = 1.0D0; C(NX) = 0.0D0

      DO IX = 1,NX  !R IS FOR RHS, THE IY GRID LINE
         R(IX) = RHS(IY,IX)
      END DO

      !----- STEP II: MODIFY LHS MATRIX ACCORDING TO MIB -----!
      IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN

         ELEM => IFPY(IY)%HEAD
         DO WHILE ( ASSOCIATED(ELEM) )

            !----- AN IRREGULAR INTERFACE POINT
            IF (ELEM%DATA%ID .GT. 0) THEN

               DATA  = LIST_GET_DATA(ELEM)

               UJP   = DATA%JUMP(1) ! [U]
               BUXJP = DATA%JUMP(6) ! [BETA U_X]

               ROW   = 0.0D0

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE LEFT GRID TO THE INTERFACE PT
               IX = DATA%LGRD
               IF ( DATA%ITYPE .EQ. 1) THEN
                  I = DATA%ITYPE ! RIGHT BIASED 3 STENCIL
               ELSE
                  I = 0 ! CENTRAL DIFF.
               END IF
               ROW(1,1) = -COEF * VDE2(I,-1)                      !1ST COFF. IN THE CENTRAL DIFF.
               ROW(1,2) =  1.0D0/BETA(IY,IX) - COEF * VDE2(I, 0)  !2ND COFF. IN THE CENTRAL DIFF.
               DO K=1,4
                  ROW(1,K) = ROW(1,K) - COEF * VDE2(I,1) * DATA%WIJ(2,K)
               END DO
               ROW(1,5) = RHS(IY,IX) + COEF * VDE2(I,1) * ( UJP * DATA%WIJ(2,5) + BUXJP * DATA%WIJ(2,6) )

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE RIGHT GRID TO THE INTERFACE PT
               IX = IX + 1
               IF ( DATA%ITYPE .EQ. -1) THEN
                  I = DATA%ITYPE ! LEFT BIASED 3 STENCIL
               ELSE
                  I = 0 ! CENTRAL DIFF.
               END IF
               ROW(2,3) =   1.0D0/BETA(IY,IX) - COEF * VDE2(I,0)  !2ND COFF. IN THE CENTRAL DIFF.
               ROW(2,4) = - COEF * VDE2(I,1)                      !3RD COFF. IN THE CENTRAL DIFF.
               DO K = 1,4
                  ROW(2,K) = ROW(2,K) - COEF * VDE2(I,-1) * DATA%WIJ(1,K)
               END DO
               ROW(2,5) = RHS(IY,IX) + COEF * VDE2(I,-1) * ( UJP * DATA%WIJ(1,5) + BUXJP * DATA%WIJ(1,6) )

               CALL GE1(NX,A,B,C,R,DATA%LGRD,ROW)

               ELEM => ELEM%NEXT ! NEXT INTERFACE

            !----- A CORNER POINT
            ELSE

               DATAL  = LIST_GET_DATA(ELEM)      ! LEFT  INTERFACE PT
               DATAR  = LIST_GET_DATA(ELEM%NEXT) ! RIGHT INTERFACE PT

               UJPL   = DATAL%JUMP(1) ! [U]        AT  LEFT INTERFACE PT
               BUXJPL = DATAL%JUMP(6) ! [BETA U_X] AT  LEFT INTERFACE PT
               UJPR   = DATAR%JUMP(1) ! [U]        AT RIGHT INTERFACE PT
               BUXJPR = DATAR%JUMP(6) ! [BETA U_X] AT RIGHT INTERFACE PT

               ROW2   = 0.0D0

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE LEFT GRID TO THE CORNER
               IX = DATAL%LGRD
               I  = 0 ! CENTRAL DIFF. FORMULA
               ROW2(1,1) = - COEF * VDE2(I,-1)                      !1ST COFF. IN THE DIFF. FORMULA
               ROW2(1,2) =   1.0D0/BETA(IY,IX) - COEF * VDE2(I, 0)  !2ND COFF. IN THE DIFF. FORMULA
               IF ( DATAL%ITYPE .EQ. -2 ) THEN
                  K = 3 ! CORNER IS THE 3RD ROW IN WIJ2(4,9)
               ELSE
                  K = 2 ! CORNER IS THE 2ND ROW IN WIJ2(4,9)
               END IF
               DO J = 1,5
                  ROW2(1,J) = ROW2(1,J) - COEF * VDE2(I,1) * DATAL%WIJ2(K,J)
               END DO
               ROW2(1,6) = RHS(IY,IX) + COEF * VDE2(I,1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                                           + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE CORNER
               IX = IX + 1
               IF ( DATAL%ITYPE .EQ. 0 .OR. DATAR%ITYPE .EQ. 0 ) THEN
                  I = 0 ! CENTRAL DIFF. FORMULA WITH DX
               ELSE
                  I = 2 ! CENTRAL DIFF. FORMULA WITH 2*DX
               END IF
               !------ THE CORNER PT
               ROW2(2,3) = 1.0D0/BETA(IY,IX) - COEF * VDE2(I,0)
               !------ THE GRID TO THE LEFT OF THE CORNER PT
               IF      ( DATAL%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 0 ) THEN
                  K = 2 ! CENTRAL DIFF. W/ DX,    THE LEFT GRID TO THE CORNER IS THE 2ND ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 1 ) THEN
                  K = 1 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ.  0 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                  K = 1 ! CENTRAL DIFF. W/ DX,   MOST LEFT GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ. -1 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                  K = 1 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
               ELSE
                  STOP
               END IF
               DO J = 1,5
                  ROW2(2,J) = ROW2(2,J) - COEF * VDE2(I,-1) * DATAL%WIJ2(K,J)
               END DO
               ROW2(2,6) = RHS(IY,IX) + COEF * VDE2(I,-1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                                            + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
               !------ THE GRID TO THE RIGHT OF THE CORNER PT
               IF      ( DATAL%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 0 ) THEN
                  K = 4 ! CENTRAL DIFF. W/ DX,    THE RIGHT GRID TO THE CORNER IS THE 4TH ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 1 ) THEN
                  K = 4 ! CENTRAL DIFF. W/ 2*DX, MOST RIGHT GRID TO THE CORNER IS THE 4TH ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ.  0 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                  K = 3 ! CENTRAL DIFF. W/ DX,   MOST RIGHT GRID TO THE CORNER IS THE 3RD ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ. -1 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                  K = 3 ! CENTRAL DIFF. W/ 2*DX, MOST RIGHT GRID TO THE CORNER IS THE 3RD ROW IN WIJ2(4,9)
               ELSE
                  STOP
               END IF
               DO J = 1,5
                  ROW2(2,J) = ROW2(2,J) - COEF * VDE2(I,1) * DATAL%WIJ2(K,J)
               END DO
               ROW2(2,6) = ROW2(2,6) + COEF * VDE2(I,1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                                          + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE RIGHT GRID TO THE CORNER
               IX = IX + 1
               I  = 0 ! CENTRAL DIFF. FORMULA
               ROW2(3,4) =   1.0D0/BETA(IY,IX) - COEF * VDE2(I,0)  !2ND COFF. IN THE DIFF. FORMULA
               ROW2(3,5) = - COEF * VDE2(I,1)                      !3RD COFF. IN THE DIFF. FORMULA
               IF ( DATAR%ITYPE .EQ. 2 ) THEN
                  K = 2 ! CORNER IS THE 2ND ROW IN WIJ2(4,9)
               ELSE
                  K = 3 ! CORNER IS THE 3RD ROW IN WIJ2(4,9)
               END IF
               DO J = 1,5
                  ROW2(3,J) = ROW2(3,J) - COEF * VDE2(I,-1) * DATAL%WIJ2(K,J)
               END DO
               ROW2(3,6) = RHS(IY,IX) + COEF * VDE2(I,-1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                                            + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )

               CALL GE2(NX,A,B,C,R,DATAL%LGRD,ROW2)

               ELEM => ELEM%NEXT%NEXT ! SKIP THE COUPLE OF CORNER INTERFACES

            END IF !---------- END OF IF (ELEM%DATA%ID .GT. 0)

         END DO !---------- END OF DO WHILE ( ASSOCIATED(ELEM) )

      END IF !---------- END OF IF ( ASSOCIATED(IFPY(IY)%HEAD) )

      RETURN

      END SUBROUTINE DXX_ADIL

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DYY_ADIL --
      !    SECOND ORDER DERIVATIVE D_YY WHEN APPEARING ON THE LEFT HAND SIDE OF
      !    ADI FORMULATION
      ! ARGUMENTS:
      !    IX     IN     A X GRID LINE
      !    COEF   IN     COEFFICIENT IN FRONT OF D_YY
      !    RHS    IN     U(NY,NX) AT CURRENT TIME STEP
      !    A      OUT    LOWER DIAGONAL ELEMENTS
      !    B      OUT    DIAGONAL ELEMENTS
      !    C      OUT    UPPER DIAGONAL ELEMENTS
      !    R      OUT    CORRESPONDING RHS VECTOR
      ! NOTES:
      !    THE SUBROUTINE CONVERTS THE MATRIX RESULTING FROM DISCRETIZING
      !    THE OPERATOR
      !         1/BETA - COEF * D_YY
      !    ON THE LEFT HAND SIDE TO A TRI-DIAGONAL MATRIX
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE DYY_ADIL(IX,COEF,RHS,A,B,C,R)

      USE MOD_DATA

      INTEGER :: IX ! A X GRID LINE
      REAL :: COEF,RHS(NY,NX),A(NY),B(NY),C(NY),R(NY)

      INTEGER,PARAMETER :: ORDER = 2 ! 2ND DERIVATIVE D_YY
      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LIST_DATA) :: DATA,DATAL,DATAR
      REAL :: VDE2(-1:2,-1:1) ! WEIGHTS FOR FD FORMULA TO APPROX. D_XX
                              ! VDE2(-1,:) - LOWER BIASED 3-STENCIL DIFF. FORMULA
                              ! VDE2( 0,:) - CENTRAL DIFF. FORMULA WITH DX
                              ! VDE2( 1,:) - UPPER BIASED 3-STENCIL DIFF. FORMULA
                              ! VDE2( 2,:) - CENTRAL DIFF. FORMULA WITH 2*DX
      INTEGER :: I,J,K,IY
      REAL :: UJP,UJPL,UJPR,BUXJP,BUXJPR,BUXJPL
      REAL :: ROW(2,5),ROW2(3,6)

      !---------- GET WEIGHTS ASSOICATED WITH EACH TYPE OF 3-STENCIL FINITE DIFF. FORMULAS
      CALL SETVDE2(ORDER,DX,VDE2)

      !----- STEP I: GENERATE LHS USING CENTRAL FD MATRIX FOR D_YY WITHOUT MIB -----!
      I = 0 ! CENTRAL DIFF. FORMULA
      A(1) = 0.0D0; B(1) = 1.0D0; C(1) = 0.0D0
      DO IY = 2,NY-1
         A(IY) = - COEF * VDE2(I,-1)
         B(IY) =   1.0D0/BETA(IY,IX) - COEF * VDE2(I, 0)
         C(IY) = - COEF * VDE2(I, 1)
      END DO
      A(NY) = 0.0D0; B(NY) = 1.0D0; C(NY) = 0.0D0

      DO IY = 1,NY  !R IS FOR RHS, THE IX GRID LINE
         R(IY) = RHS(IY,IX)
      END DO

      !----- STEP II: MODIFY LHS MATRIX ACCORDING TO MIB -----!
      IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN
         ELEM => IFPX(IX)%HEAD
         DO WHILE ( ASSOCIATED(ELEM) )

            !----- AN IRREGULAR INTERFACE POINT
            IF (ELEM%DATA%ID .GT. 0) THEN

               DATA  = LIST_GET_DATA(ELEM)

               UJP   = DATA%JUMP(1) ! [U]
               BUXJP = DATA%JUMP(6) ! [BEAT U_Y]

               ROW   = 0.0D0

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE LOWER GRID TO THE INTERFACE PT
               IY = DATA%LGRD
               IF ( DATA%ITYPE .EQ. 1) THEN
                  I = DATA%ITYPE ! UPPER BIASED 3 STENCIL
               ELSE
                  I = 0 ! CENTRAL DIFF.
               END IF
               ROW(1,1) = - COEF * VDE2(I,-1)                      !1ST COFF. IN THE CENTRAL DIFF.
               ROW(1,2) =   1.0D0/BETA(IY,IX) - COEF * VDE2(I, 0)  !2ND COFF. IN THE CENTRAL DIFF.
               DO K = 1,4
                  ROW(1,K) = ROW(1,K) - COEF * VDE2(I,1) * DATA%WIJ(2,K)
               END DO
               ROW(1,5) = RHS(IY,IX) + COEF * VDE2(I,1) * ( UJP * DATA%WIJ(2,5) + BUXJP * DATA%WIJ(2,6) )

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE UPPER GRID TO THE INTERFACE PT
               IY = IY + 1
               IF ( DATA%ITYPE .EQ. -1) THEN
                  I = DATA%ITYPE ! LOWER BIASED 3 STENCIL
               ELSE
                  I = 0 ! CENTRAL DIFF.
               END IF
               ROW(2,3) =   1.0D0/BETA(IY,IX) - COEF * VDE2(I,0)  !2ND COFF. IN THE CENTRAL DIFF.
               ROW(2,4) = - COEF * VDE2(I,1)                      !3RD COFF. IN THE CENTRAL DIFF.
               DO K = 1,4
                  ROW(2,K) = ROW(2,K) - COEF * VDE2(I,-1) * DATA%WIJ(1,K)
               END DO
               ROW(2,5) = RHS(IY,IX) + COEF * VDE2(I,-1) * ( UJP * DATA%WIJ(1,5) + BUXJP * DATA%WIJ(1,6) )

               CALL GE1(NY,A,B,C,R,DATA%LGRD,ROW)

               ELEM => ELEM%NEXT ! NEXT INTERFACE

            !----- A CORNER POINT
            ELSE

               DATAL  = LIST_GET_DATA(ELEM)      ! LOWER INTERFACE PT
               DATAR  = LIST_GET_DATA(ELEM%NEXT) ! UPPER INTERFACE PT

               UJPL   = DATAL%JUMP(1) ! [U]        AT LOWER INTERFACE PT
               UJPR   = DATAR%JUMP(1) ! [BETA U_Y] AT LOWER INTERFACE PT
               BUXJPL = DATAL%JUMP(6) ! [U]        AT UPPER INTERFACE PT
               BUXJPR = DATAR%JUMP(6) ! [BETA U_Y] AT UPPER INTERFACE PT

               ROW2   = 0.0D0

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE LOWER GRID TO THE CORNER
               IY = DATAL%LGRD
               I  = 0 ! CENTRAL DIFF. FORMULA
               ROW2(1,1) = - COEF * VDE2(I,-1)                      !1ST COFF. IN THE DIFF. FORMULA
               ROW2(1,2) =   1.0D0/BETA(IY,IX) - COEF * VDE2(I, 0)  !2ND COFF. IN THE DIFF. FORMULA
               IF ( DATAL%ITYPE .EQ. -2 ) THEN
                  K = 3 ! CORNER IS THE 3RD ROW IN WIJ2(4,9)
               ELSE
                  K = 2 ! CORNER IS THE 2ND ROW IN WIJ2(4,9)
               END IF
               DO J = 1,5
                  ROW2(1,J) = ROW2(1,J) - COEF * VDE2(I,1) * DATAL%WIJ2(K,J)
               END DO
               ROW2(1,6) = RHS(IY,IX) + COEF * VDE2(I,1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                                           + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE CORNER
               IY = IY + 1
               IF ( DATAL%ITYPE .EQ. 0 .OR. DATAR%ITYPE .EQ. 0 ) THEN
                  I = 0 ! CENTRAL DIFF. FORMULA WITH DX
               ELSE
                  I = 2 ! CENTRAL DIFF. FORMULA WITH 2*DX
               END IF
               !------ THE CORNER PT
               ROW2(2,3) = 1.0D0/BETA(IY,IX) - COEF * VDE2(I,0)
               !------ THE GRID TO THE LOWER OF THE CORNER PT
               IF      ( DATAL%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 0 ) THEN
                  K = 2 ! CENTRAL DIFF. W/ DX,    THE LOWER GRID TO THE CORNER IS THE 2ND ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 1 ) THEN
                  K = 1 ! CENTRAL DIFF. W/ 2*DX, MOST LOWER GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ.  0 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                  K = 1 ! CENTRAL DIFF. W/ DX,   MOST LOWER GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ. -1 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                  K = 1 ! CENTRAL DIFF. W/ 2*DX, MOST LOWER GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
               ELSE
                  STOP
               END IF
               DO J = 1,5
                  ROW2(2,J) = ROW2(2,J) - COEF * VDE2(I,-1) * DATAL%WIJ2(K,J)
               END DO
               ROW2(2,6) = RHS(IY,IX) + COEF * VDE2(I,-1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                                            + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
               !------ THE GRID TO THE UPPER OF THE CORNER PT
               IF      ( DATAL%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 0 ) THEN
                  K = 4 ! CENTRAL DIFF. W/ DX,    THE UPPER GRID TO THE CORNER IS THE 4TH ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 1 ) THEN
                  K = 4 ! CENTRAL DIFF. W/ 2*DX, MOST UPPER GRID TO THE CORNER IS THE 4TH ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ.  0 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                  K = 3 ! CENTRAL DIFF. W/ DX,   MOST UPPER GRID TO THE CORNER IS THE 3RD ROW IN WIJ2(4,9)
               ELSE IF ( DATAR%ITYPE .EQ. -1 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                  K = 3 ! CENTRAL DIFF. W/ 2*DX, MOST UPPER GRID TO THE CORNER IS THE 3RD ROW IN WIJ2(4,9)
               ELSE
                  STOP
               END IF
               DO J = 1,5
                  ROW2(2,J) = ROW2(2,J) - COEF * VDE2(I,1) * DATAL%WIJ2(K,J)
               END DO
               ROW2(2,6) = ROW2(2,6) + COEF * VDE2(I,1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                                          + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )

               !+++++ MODIFY FINITE DIFF. FORMULA AT THE UPPER GRID TO THE CORNER
               IY = IY + 1
               I  = 0 ! CENTRAL DIFF. FORMULA
               ROW2(3,4) =   1.0D0/BETA(IY,IX) - COEF * VDE2(I,0)  !2ND COFF. IN THE DIFF. FORMULA
               ROW2(3,5) = - COEF * VDE2(I,1)                      !3RD COFF. IN THE DIFF. FORMULA
               IF ( DATAR%ITYPE .EQ. 2 ) THEN
                  K = 2 ! CORNER IS THE 2ND ROW IN WIJ2(4,9)
               ELSE
                  K = 3 ! CORNER IS THE 3RD ROW IN WIJ2(4,9)
               END IF
               DO J = 1,5
                  ROW2(3,J) = ROW2(3,J) - COEF * VDE2(I,-1) * DATAL%WIJ2(K,J)
               END DO
               ROW2(3,6) = RHS(IY,IX) + COEF * VDE2(I,-1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                                            + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )

               CALL GE2(NY,A,B,C,R,DATAL%LGRD,ROW2)

               ELEM => ELEM%NEXT%NEXT ! SKIP THE COUPLE OF CORNER INTERFACES

            END IF !---------- END OF IF (ELEM%DATA%ID .GT. 0)

         END DO !---------- END OF DO WHILE ( ASSOCIATED(ELEM) )

      END IF !---------- END OF IF ( ASSOCIATED(IFPX(IX)%HEAD) )

      RETURN

      END SUBROUTINE DYY_ADIL

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DXX_RHS --
      !    SECOND ORDER DERIVATIVE D_XX WHEN APPEARING ON THE RIGHT HAND SIDE OF
      !    ADI FORMULATION
      ! ARGUMENTS:
      !    U      IN     U(NY,NX) AT CURRENT TIME STEP
      !    UXX    OUT    DXX APPLIED TO U AT CURRENT TIME STEP
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE DXX_RHS(U,UXX)

      USE MOD_DATA

      REAL :: U(NY,NX),UXX(NY,NX)

      INTEGER,PARAMETER :: ORDER = 2 !2ND DERIVATIVE D_XX
      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LIST_DATA) :: DATA,DATAL,DATAR
      REAL :: VDE2(-1:2,-1:1) ! WEIGHTS FOR FD FORMULA TO APPROX. D_XX
                              ! VDE2(-1,:) - LEFT BIASED 3-STENCIL DIFF. FORMULA
                              ! VDE2( 0,:) - CENTRAL DIFF. FORMULA WITH DX
                              ! VDE2( 1,:) - RIGHT BIASED 3-STENCIL DIFF. FORMULA
                              ! VDE2( 2,:) - CENTRAL DIFF. FORMULA WITH 2*DX
      INTEGER :: I,J,K,IX,IY
      REAL :: SUM,UJP,UJPL,UJPR,BUXJP,BUXJPR,BUXJPL

      !---------- GET WEIGHTS ASSOICATED WITH EACH TYPE OF 3-STENCIL FINITE DIFF. FORMULAS
      CALL SETVDE2(ORDER,DX,VDE2)

      DO IY = 2,NY-1 !---------- FOR EACH IY GRID LINE

         !---------- STEP I: SET UXX = D_XX U WITHOUT MIB
         I = 0 ! CENTRAL DIFF. FORMULA
         DO IX = 2,NX-1  !INITLIAL WITHOUT MIB
            SUM = 0.0D0
            DO J = -1,1
               SUM = SUM + U(IY,IX+J) * VDE2(I,J)
            END DO
            UXX(IY,IX) = SUM
         END DO

         !---------- STEP II: MODIFY 3 STENCIL DIFF FORMULA AT GRIDS NEAR INTERFACE ACCORDING TO MIB
         IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN
            ELEM => IFPY(IY)%HEAD

            DO WHILE ( ASSOCIATED(ELEM) )

               !----- AN IRREGULAR INTERFACE POINT
               IF (ELEM%DATA%ID .GT. 0) THEN
                  DATA  = LIST_GET_DATA(ELEM)

                  UJP   = DATA%JUMP(1) ! [U]
                  BUXJP = DATA%JUMP(6) ! [BETA U_X]

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE LEFT GRID TO THE INTERFACE PT
                  IX  = DATA%LGRD
                  IF ( DATA%ITYPE .EQ. 1 ) THEN
                     I = DATA%ITYPE ! RIGHT BIASED 3 STENCIL
                  ELSE
                     I = 0 ! CENTRAL DIFF.
                  END IF
                  SUM = 0.0D0
                  DO J = -1,0
                     SUM = SUM + VDE2(I,J) * U(IY,IX+J)
                  END DO
                  DO K = 1,4
                     SUM = SUM + VDE2(I,1) * U(IY,DATA%LGRD-2+K) * DATA%WIJ(2,K)
                  END DO
                  SUM  = SUM + VDE2(I,1) * ( UJP * DATA%WIJ(2,5) + BUXJP * DATA%WIJ(2,6) )
                  UXX(IY,IX) = SUM

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE RIGHT GRID TO THE INTERFACE PT
                  IX = IX + 1
                  IF ( DATA%ITYPE .EQ. -1) THEN
                     I = DATA%ITYPE ! LEFT BIASED 3 STENCIL
                  ELSE
                     I = 0 ! CENTRAL DIFF.
                  END IF
                  SUM = 0.0D0
                  DO J = 0,1
                     SUM = SUM + VDE2(I,J) * U(IY,IX+J)
                  END DO
                  DO K = 1,4
                     SUM = SUM + VDE2(I,-1) * U(IY,DATA%LGRD-2+K) * DATA%WIJ(1,K)
                  END DO
                  SUM = SUM + VDE2(I,-1) * ( UJP * DATA%WIJ(1,5) + BUXJP * DATA%WIJ(1,6) )
                  UXX(IY,IX) = SUM

                  ELEM => ELEM%NEXT  ! NEXT INTERFACE PT

               !----- A CORNER POINT
               ELSE
                  DATAL  = LIST_GET_DATA(ELEM)      ! LEFT  INTERFACE PT
                  DATAR  = LIST_GET_DATA(ELEM%NEXT) ! RIGHT INTERFACE PT

                  UJPL   = DATAL%JUMP(1) ! [U]        AT THE  LEFT INTERFACE PT
                  BUXJPL = DATAL%JUMP(6) ! [BETA U_X] AT THE  LEFT INTERFACE PT
                  UJPR   = DATAR%JUMP(1) ! [U]        AT THE RIGHT INTERFACE PT
                  BUXJPR = DATAR%JUMP(6) ! [BETA U_X] AT THE RIGHT INTERFACE PT

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE LEFT GRID TO THE CORNER
                  IX  = DATAL%LGRD
                  I   = 0 ! CENTRAL DIFF. FORMULA
                  SUM  = 0.0D0
                  DO J = -1,0
                     SUM = SUM + VDE2(I,J) * U(IY,IX+J)
                  END DO
                  IF ( DATAL%ITYPE .EQ. -2 ) THEN
                     K = 3 ! CORNER IS THE 3RD ROW IN WIJ2(4,9)
                  ELSE
                     K = 2 ! CORNER IS THE 2ND ROW IN WIJ2(4,9)
                  END IF
                  DO J = 1,5
                     SUM = SUM + VDE2(I,1) * U(IY,DATAL%LGRD-2+J) * DATAL%WIJ2(K,J)
                  END DO
                  SUM = SUM + VDE2(I,1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                          + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
                  UXX(IY,IX) = SUM

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE CORNER
                  IX = IX + 1
                  IF ( DATAL%ITYPE .EQ. 0 .OR. DATAR%ITYPE .EQ. 0 ) THEN
                     I = 0 ! CENTRAL DIFF. FORMULA WITH DX
                  ELSE
                     I = 2 ! CENTRAL DIFF. FORMULA WITH 2*DX
                  END IF
                  SUM = 0.0D0
                  !------ THE GRID TO THE LEFT OF THE CORNER PT
                  IF      ( DATAL%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 0 ) THEN
                     K = 2 ! CENTRAL DIFF. W/ DX,   THE LEFT GRID TO THE CORNER IS THE 2ND ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 1 ) THEN
                     K = 1 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ.  0 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                     K = 1 ! CENTRAL DIFF. W/ DX,   MOST LEFT GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ. -1 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                     K = 1 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
                  ELSE
                     STOP
                  END IF
                  DO J = 1,5
                     SUM = SUM + VDE2(I,-1) * U(IY,DATAL%LGRD-2+J) * DATAL%WIJ2(K,J)
                  END DO
                  SUM = SUM + VDE2(I,-1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                           + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
                  !------ THE CORNER PT
                  SUM = SUM + VDE2(I,0) * U(IY,IX)
                  !------ THE GRID TO THE RIGHT OF THE CORNER PT
                  IF      ( DATAL%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 0 ) THEN
                     K = 4 ! CENTRAL DIFF. W/ DX,    THE LEFT GRID TO THE CORNER IS THE 4TH ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 1 ) THEN
                     K = 4 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 4TH ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ.  0 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                     K = 3 ! CENTRAL DIFF. W/ DX,   MOST LEFT GRID TO THE CORNER IS THE 3RD ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ. -1 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                     K = 3 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 3RD ROW IN WIJ2(4,9)
                  ELSE
                     STOP
                  END IF
                  DO J = 1,5
                     SUM = SUM + VDE2(I,1) * U(IY,DATAL%LGRD-2+J) * DATAL%WIJ2(K,J)
                  END DO
                  SUM = SUM + VDE2(I,1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                          + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
                  UXX(IY,IX) = SUM

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE RIGHT GRID TO THE CORNER
                  IX  = IX + 1
                  I   = 0 ! CENTRAL DIFF. FORMULA
                  SUM = 0.0D0
                  DO J = 0,1
                     SUM = SUM + VDE2(I,J) * U(IY,IX+J)
                  END DO
                  IF ( DATAR%ITYPE .EQ. 2 ) THEN
                     K = 2 ! CORNER IS THE 2ND ROW IN WIJ2(4,9)
                  ELSE
                     K = 3 ! CORNER IS THE 3RD ROW IN WIJ2(4,9)
                  END IF
                  DO J = 1,5
                     SUM = SUM + VDE2(I,-1) * U(IY,DATAL%LGRD-2+J) * DATAL%WIJ2(K,J)
                  END DO
                  SUM = SUM + VDE2(I,-1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                           + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
                  UXX(IY,IX) = SUM

                  ELEM => ELEM%NEXT%NEXT  ! SKIP THE PAIR OF CORNER INTERFACE PTS

               END IF !---------- END OF IF (ELEM%DATA%ID .GT. 0)

            END DO !---------- END OF DO WHILE ( ASSOCIATED(ELEM) )

         END IF !---------- END OF IF ( ASSOCIATED(IFPY(IY)%HEAD) )

      END DO !---------- END OF FOR EACH IY GRID LINE

      RETURN

      END SUBROUTINE DXX_RHS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! DYY_RHS --
      !    SECOND ORDER DERIVATIVE D_YY WHEN APPEARING ON THE RIGHT HAND SIDE OF
      !    ADI FORMULATION
      ! ARGUMENTS:
      !    U      IN     U(NY,NX) AT CURRENT TIME STEP
      !    UYY    OUT    DYY APPLIED TO U AT CURRENT TIME STEP
      ! NOTES:
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE DYY_RHS(U,UYY)

      USE MOD_DATA

      REAL :: U(NY,NX),UYY(NY,NX)

      INTEGER,PARAMETER :: ORDER = 2 !2ND DERIVATIVE D_YY
      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LIST_DATA) :: DATA,DATAL,DATAR
      REAL :: VDE2(-1:2,-1:1) ! WEIGHTS FOR FD FORMULA TO APPROX. D_XX
                              ! VDE2(-1,:) - LEFT BIASED 3-STENCIL DIFF. FORMULA
                              ! VDE2( 0,:) - CENTRAL DIFF. FORMULA WITH DX
                              ! VDE2( 1,:) - RIGHT BIASED 3-STENCIL DIFF. FORMULA
                              ! VDE2( 2,:) - CENTRAL DIFF. FORMULA WITH 2*DX
      INTEGER :: I,J,K,IX,IY
      REAL :: SUM,UJP,UJPL,UJPR,BUXJP,BUXJPR,BUXJPL

      !---------- GET WEIGHTS ASSOICATED WITH EACH TYPE OF 3-STENCIL FINITE DIFF. FORMULAS
      CALL SETVDE2(ORDER,DX,VDE2)

      DO IX=2,NX-1 !---------- FOR EACH IX GRID LINE

         !---------- STEP I: SET UYY = D_YY U WITHOUT MIB
         I = 0 ! CENTRAL DIFF. FORMULA
         DO IY=2,NY-1 !INITIAL WITHOUT MIB
            SUM = 0.0D0
            DO J = -1,1
               SUM = SUM + U(IY+J,IX)*VDE2(I,J)
            END DO
            UYY(IY,IX) = SUM
         END DO

         !---------- STEP II: MODIFY 3 STENCIL DIFF FORMULA AT GRIDS NEAR INTERFACE ACCORDING TO MIB
         IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN
            ELEM => IFPX(IX)%HEAD
            DO WHILE ( ASSOCIATED(ELEM) )

               !----- AN IRREGULAR INTERFACE POINT
               IF ( ELEM%DATA%ID .GT. 0 ) THEN
                  DATA = LIST_GET_DATA(ELEM)

                  UJP   = DATA%JUMP(1) ! [U]
                  BUXJP = DATA%JUMP(6) ! [BETA U_Y]

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE LOWER GRID TO THE INTERFACE PT
                  IY = DATA%LGRD
                  IF ( DATA%ITYPE .EQ. 1) THEN
                     I = DATA%ITYPE ! RIGHT BIASED 3 STENCIL
                  ELSE
                     I = 0 ! CENTRAL DIFF.
                  END IF
                  SUM = 0.0D0
                  DO J = -1,0
                     SUM = SUM + VDE2(I,J) * U(IY+J,IX)
                  END DO
                  DO K=1,4
                     SUM = SUM + VDE2(I,1) * U(DATA%LGRD-2+K,IX) * DATA%WIJ(2,K)
                  END DO
                  SUM = SUM + VDE2(I,1) * ( UJP * DATA%WIJ(2,5) + BUXJP * DATA%WIJ(2,6) )
                  UYY(IY,IX) = SUM

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE UPPER GRID TO THE INTERFACE PT
                  IY = IY + 1
                  IF ( DATA%ITYPE .EQ. -1) THEN
                     I = DATA%ITYPE ! LEFT BIASED 3 STENCIL
                  ELSE
                     I = 0 ! CENTRAL DIFF.
                  END IF
                  SUM = 0.0D0
                  DO J = 0,1
                     SUM = SUM + VDE2(I,J) * U(IY+J,IX)
                  END DO
                  DO K = 1,4
                     SUM = SUM + VDE2(I,-1) * U(DATA%LGRD-2+K,IX) * DATA%WIJ(1,K)
                  END DO
                  SUM = SUM + VDE2(I,-1) * ( UJP * DATA%WIJ(1,5) + BUXJP * DATA%WIJ(1,6) )
                  UYY(IY,IX) = SUM

                  ELEM => ELEM%NEXT ! NEXT INTERFACE

               !----- A CORNER POINT
               ELSE
                  DATAL  = LIST_GET_DATA(ELEM)      ! LEFT  INTERFACE PT
                  DATAR  = LIST_GET_DATA(ELEM%NEXT) ! RIGHT INTERFACE PT

                  UJPL   = DATAL%JUMP(1) ! [U]        AT THE  LEFT INTERFACE PT
                  BUXJPL = DATAL%JUMP(6) ! [BETA U_X] AT THE  LEFT INTERFACE PT
                  UJPR   = DATAR%JUMP(1) ! [U]        AT THE RIGHT INTERFACE PT
                  BUXJPR = DATAR%JUMP(6) ! [BETA U_X] AT THE RIGHT INTERFACE PT

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE LOWER GRID TO THE CORNER
                  IY  = DATAL%LGRD
                  I   = 0 ! CENTRAL DIFF. FORMULA
                  SUM = 0.0D0
                  DO J = -1,0
                     SUM = SUM + VDE2(I,J)*U(IY+J,IX)
                  END DO
                  IF ( DATAL%ITYPE .EQ. -2 ) THEN
                     K = 3 ! CORNER IS THE 3RD ROW IN WIJ2(4,9)
                  ELSE
                     K = 2 ! CORNER IS THE 2ND ROW IN WIJ2(4,9)
                  END IF
                  DO J = 1,5
                     SUM = SUM + VDE2(I,1) * U(DATAL%LGRD-2+J,IX) * DATAL%WIJ2(K,J)
                  END DO
                  SUM = SUM + VDE2(I,1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                          + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
                  UYY(IY,IX) = SUM

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE CORNER
                  IY = IY + 1
                  IF ( DATAL%ITYPE .EQ. 0 .OR. DATAR%ITYPE .EQ. 0 ) THEN
                     I = 0 ! CENTRAL DIFF. FORMULA WITH DX
                  ELSE
                     I = 2 ! CENTRAL DIFF. FORMULA WITH 2*DX
                  END IF
                  SUM = 0.0D0
                  !------ THE GRID TO THE LOWER OF THE CORNER PT
                  IF      ( DATAL%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 0 ) THEN
                     K = 2 ! CENTRAL DIFF. W/ DX,   THE LEFT GRID TO THE CORNER IS THE 2ND ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 1 ) THEN
                     K = 1 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ.  0 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                     K = 1 ! CENTRAL DIFF. W/ DX,   MOST LEFT GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ. -1 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                     K = 1 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 1ST ROW IN WIJ2(4,9)
                  ELSE
                     STOP
                  END IF
                  DO J = 1,5
                     SUM = SUM + VDE2(I,-1) * U(DATAL%LGRD-2+J,IX) * DATAL%WIJ2(K,J)
                  END DO
                  SUM = SUM + VDE2(I,-1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                           + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
                  !------ THE CORNER PT
                  SUM = SUM + VDE2(I,0)*U(IY,IX)
                  !------ THE GRID TO THE UPPER OF THE CORNER PT
                  IF      ( DATAL%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 0 ) THEN
                     K = 4 ! CENTRAL DIFF. W/ DX,    THE LEFT GRID TO THE CORNER IS THE 4TH ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ. -2 .AND. DATAR%ITYPE .EQ. 1 ) THEN
                     K = 4 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 4TH ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ.  0 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                     K = 3 ! CENTRAL DIFF. W/ DX,   MOST LEFT GRID TO THE CORNER IS THE 3RD ROW IN WIJ2(4,9)
                  ELSE IF ( DATAR%ITYPE .EQ. -1 .AND. DATAR%ITYPE .EQ. 2 ) THEN
                     K = 3 ! CENTRAL DIFF. W/ 2*DX, MOST LEFT GRID TO THE CORNER IS THE 3RD ROW IN WIJ2(4,9)
                  ELSE
                     STOP
                  END IF
                  DO J = 1,5
                     SUM = SUM + VDE2(I,1) * U(DATAL%LGRD-2+J,IX) * DATAL%WIJ2(K,J)
                  END DO
                  SUM = SUM + VDE2(I,1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                          + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
                  UYY(IY,IX) = SUM

                  !+++++ MODIFY FINITE DIFF. FORMULA AT THE UPPER GRID TO THE CORNER
                  IY  = IY + 1
                  I   = 0 ! CENTRAL DIFF. FORMULA
                  SUM = 0.0D0
                  DO J = 0,1
                     SUM = SUM + VDE2(I,J)*U(IY+J,IX)
                  END DO
                  IF ( DATAR%ITYPE .EQ. 2 ) THEN
                     K = 2 ! CORNER IS THE 2ND ROW IN WIJ2(4,9)
                  ELSE
                     K = 3 ! CORNER IS THE 3RD ROW IN WIJ2(4,9)
                  END IF
                  DO J = 1,5
                     SUM = SUM + VDE2(I,-1) * U(DATAL%LGRD-2+J,IX) * DATAL%WIJ2(K,J)
                  END DO
                  SUM = SUM + VDE2(I,-1) * ( UJPL * DATAL%WIJ2(K,6) + BUXJPL * DATAL%WIJ2(K,7) &
                                           + UJPR * DATAR%WIJ2(K,8) + BUXJPR * DATAR%WIJ2(K,9) )
                  UYY(IY,IX) = SUM

                  ELEM => ELEM%NEXT%NEXT ! SKIP THE COUPLE OF CORNER INTERFACES

               END IF !---------- END OF IF (ELEM%DATA%ID .GT. 0)

            END DO !---------- END OF DO WHILE ( ASSOCIATED(ELEM) )

         END IF !---------- END OF IF ( ASSOCIATED(IFPY(IX)%HEAD) )

      END DO !---------- END OF FOR EACH IX GRID LINE

      RETURN

      END SUBROUTINE DYY_RHS

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! SETVDE2 --
      !    SET FINITE DIFF. WEIGHTS FOR APPROX. 1ST OR 2ND ORDER DERIVATIVES
      !    USING 3-POINT STENCIL
      ! ARGUMENTS:
      !    ORDE   IN     ORDER OF DERIVATIVE
      !    DX     IN     SPACIAL STEP
      !    W      OUT    CALCULATED WEIGHTS
      ! NOTES:
      !    THE TARGET POSITION IS FIXED AT THE MIDDLE PT
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SETVDE2(ORDER,DX,W)

      INTEGER :: ORDER ! DERIVATIVE ORDER
      REAL :: DX
      REAL :: W(-1:2,-1:1) ! WEIGHTS FOR FD FORMULA TO APPROX. D_XX
                           ! W(-1,:) - LEFT BIASED 3-STENCIL DIFF. FORMULA
                           ! W( 0,:) - CENTRAL DIFF. FORMULA WITH DX
                           ! W( 1,:) - RIGHT BIASED 3-STENCIL DIFF. FORMULA
                           ! W( 2,:) - CENTRAL DIFF. FORMULA WITH 2*DX

      REAL :: XP(-1:2,0:2),CP(0:2,0:2),Z ! NODES, WEIGHTS AND TARGET PT USED IN SUBROUTINE WEIGHTS
      INTEGER :: I,J

      Z = 0.0D0 ! ALWAYS ASSUME THE TARGET POSITIONED AT 0
      XP(-1,0) = -2.0D0*DX; XP(-1,1) = 0.0D0; XP(-1,2) =       DX; ! LEFT BIASED 3 STENCIL
      XP( 0,0) =       -DX; XP( 0,1) = 0.0D0; XP( 0,2) =       DX; ! CENTRAL DIFF. W/ DX
      XP( 1,0) =       -DX; XP( 1,1) = 0.0D0; XP( 1,2) = 2.0D0*DX; ! RIGHT BIASED 3 STENCIL
      XP( 2,0) = -2.0D0*DX; XP( 1,1) = 0.0D0; XP( 1,2) = 2.0D0*DX; ! CENTRAL DIFF. W/ 2*DX

      !---------- APPROX. 2ND ORDER DERIVATIVE
      DO I = -1,2
         CP = 0.0D0
         CALL WEIGHTS(Z,XP(I,:),2,2,ORDER,CP)
         DO J = 0,2
            W(I,J-1) = CP(J,ORDER)
         END DO
      END DO

      RETURN

      END SUBROUTINE SETVDE2

