      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !   ELEMENTARY ROW OPERATIONS, GUASSIAN ELIMINATION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE GE1(N,A,B,C,R,IX,ROW)

      REAL :: A(N),B(N),C(N),R(N),ROW(2,5)
      INTEGER :: N,IX

      REAL :: RATE
      INTEGER :: J

      RATE = ROW(2,1)/ROW(1,1)!ELIMINATE ROW(2,1)
      DO J = 1,5
         ROW(2,J) = ROW(2,J)-RATE*ROW(1,J)
      END DO

      RATE = ROW(1,4)/ROW(2,4)!ELIMINATE ROW(1,4)
      DO J = 1,5
         ROW(1,J) = ROW(1,J)-RATE*ROW(2,J)
      END DO

      DO J = IX,IX+1 !AFTER ELIMINATION, IT IS TRIDIAGONAL NOW
         A(J) = ROW(J-IX+1,J-IX+1)
         B(J) = ROW(J-IX+1,J-IX+2)
         C(J) = ROW(J-IX+1,J-IX+3)
         R(J) = ROW(J-IX+1,5)
      END DO

      RETURN

      END SUBROUTINE GE1

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !   ELEMENTARY ROW OPERATIONS, GUASSIAN ELIMINATION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE GE2(N,A,B,C,R,IX,ROW2)

      REAL :: A(N),B(N),C(N),R(N),ROW2(3,6)
      INTEGER :: N,IX

      REAL :: RATE1,RATE2,RATE3
      INTEGER :: J

      RATE2 = ROW2(2,1)/ROW2(1,1)!ELIMINATE ROW2(2,1)
      RATE3 = ROW2(3,1)/ROW2(1,1)!ELIMINATE ROW2(3,1)
      DO J = 1,6
         ROW2(2,J) = ROW2(2,J)-RATE2*ROW2(1,J)
         ROW2(3,J) = ROW2(3,J)-RATE3*ROW2(1,J)
      END DO

      RATE1 = ROW2(1,5)/ROW2(3,5)!ELIMINATE ROW2(1,5)
      RATE2 = ROW2(2,5)/ROW2(3,5)!ELIMINATE ROW2(2,5)
      DO J = 1,6
         ROW2(1,J) = ROW2(1,J)-RATE1*ROW2(3,J)
         ROW2(2,J) = ROW2(2,J)-RATE2*ROW2(3,J)
      END DO

      RATE3 = ROW2(3,2)/ROW2(2,2)!ELIMINATE ROW2(3,2)
      DO J = 1,6
         ROW2(3,J) = ROW2(3,J)-RATE3*ROW2(2,J)
      END DO

      RATE1= ROW2(1,4)/ROW2(2,4)!ELIMINATE ROW2(1,4)
      DO J = 1,6
         ROW2(1,J) = ROW2(1,J) - RATE1*ROW2(2,J)
      END DO

      DO J = IX,IX+2 !AFTER ELIMINATION, IT IS TRIDIAGONAL NOW
         A(J) = ROW2(J-IX+1,J-IX+1)
         B(J) = ROW2(J-IX+1,J-IX+2)
         C(J) = ROW2(J-IX+1,J-IX+3)
         R(J) = ROW2(J-IX+1,6)
      END DO

      RETURN

      END SUBROUTINE GE2

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !   FROM BOOK: NUMERICAL RECEIPT.
      !   SOLVES FOR A VECTOR U(1:N) OF LENGTH N THE TRIDIAGONAL LINEAR SET GIVEN BY EQUATION(2.4.1).
      !   A(1:N), B(1:N), C(1:N), AND R(1:N) ARE INPUT VECTORS AND ARE NOT MODIFIED.
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TRIDAG(A,B,C,R,U,N)

      INTEGER :: N
      REAL :: A(N),B(N),C(N),R(N),U(N)

      INTEGER :: J
      REAL :: BET,GAM(N)!ONE VECTOR OF WORKSPACE, GAM IS NEEDED.

      !IF THIS HAPPENS THEN YOU SHOULD REWRITE YOUR EQUATIONS AS A SET OF ORDER N-1, WITH U2 TRIVIALLY ELIMINATED.
      IF (ABS(B(1)) .LT. 1.0D-15) THEN
         WRITE (*,*) 'TRIDAG: REWRITE EQUATIONS'
         STOP
      END IF

      BET  = B(1)
      U(1) = R(1) / BET
      DO J = 2,N
         !DECOMPOSITION AND FORWARD SUBSTITUTION.
         GAM(J) = C(J-1)/BET
         BET    = B(J)-A(J)*GAM(J)

         !ALGORITHM FAILS; SEE BOOK
         IF (ABS(BET) .LT. 1.0D-15) THEN
            WRITE (*,*) 'TRIDAG FAILED'
            STOP
         END IF

         U(J) = (R(J)-A(J)*U(J-1))/BET
      END DO

      !BACKSUBSTITUTION.
      DO J = N-1,1,-1
         U(J) = U(J)-GAM(J+1)*U(J+1)
      END DO

      RETURN

      END SUBROUTINE TRIDAG

