      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !   DIMENSION INDEX1(N),VV(N)
      !   CALL LUDCMP(U,N,N,INDEX1,D,VV)
      !   CALL LUBKSB(U,N,N,INDEX1,B)
      !
      !   GUASSIAN ELELIMINATION METHOD TO FIND THE SOLUTION
      !   OF LINEAR ALGEBRA SYSTEM, COLUMN MAIN ENTRY SKILL
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE LUDCMP(A,N,NP,INDX,D,VV)

      REAL,PARAMETER :: TINY = 1.D-20
      INTEGER :: N,NP
      REAL :: D
      REAL :: A(NP,NP),VV(N)
      INTEGER :: INDX(N)

      INTEGER :: I,J,K,IMAX
      REAL :: AAMAX,SUM,DUM

      D = 1.D0
      DO I = 1,N
         VV(I) = 0.D0
      END DO

      DO 12 I = 1,N
         AAMAX = 0.D0
         DO 11 J = 1,N
            IF ( ABS(A(I,J)) .GT. AAMAX) AAMAX = ABS(A(I,J) )
11       CONTINUE

         !IF (AAMAX.EQ.0.D0) PAUSE 'SINGULAR MATRIX'
         VV(I) = 1.D0/AAMAX
12    CONTINUE

      DO 19 J = 1,N
         IF (J .GT. 1) THEN
            DO 14 I = 1,J-1
               SUM = A(I,J)
               IF (I .GT. 1) THEN
                  DO 13 K = 1,I-1
                     SUM = SUM - A(I,K) * A(K,J)
13                CONTINUE
                  A(I,J) = SUM
               END IF
14          CONTINUE
         END IF
         AAMAX = 0.D0
         DO 16 I = J,N
            SUM = A (I,J)
            IF (J .GT. 1) THEN
               DO 15 K = 1,J-1
                  SUM = SUM - A(I,K) * A(K,J)
15             CONTINUE
               A (I,J) = SUM
            END IF
            DUM = VV(I) * ABS(SUM)
            IF (DUM .GE. AAMAX) THEN
               IMAX  = I
               AAMAX = DUM
            END IF
16       CONTINUE
         IF (J .NE. IMAX) THEN
            DO 17 K = 1,N
               DUM       = A(IMAX,K)
               A(IMAX,K) = A(J,K)
               A(J,K)    = DUM
17          CONTINUE
            D = -D
            VV(IMAX) = VV(J)
         END IF
         INDX(J) = IMAX
         IF (J .NE. N) THEN
            IF (A(J,J) .EQ. 0.D0) A (J,J) = TINY
            DUM = 1.D0/A(J,J)
            DO 18 I = J+1, N
               A (I,J) = A (I,J) * DUM
18          CONTINUE
         END IF
19    CONTINUE
      IF (A(N,N) .EQ. 0.D0) A(N,N) = TINY
      RETURN

      END SUBROUTINE LUDCMP

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !   BACKWARD SUBSTITUTION
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

      INTEGER :: N,NP      
      REAL :: A(NP,NP),B(N)
      INTEGER :: INDX(N)

      INTEGER :: II,I,J,LL
      REAL :: SUM

      II = 0
      DO 12 I = 1,N
         LL    = INDX(I)
         SUM   = B(LL)
         B(LL) = B(I)
         IF (II .NE. 0) THEN
            DO 11 J = II,I-1
               SUM = SUM - A(I,J) * B(J)
11          CONTINUE
         ELSE IF (SUM .NE. 0.D0) THEN
            II = I
         END IF
         B(I) = SUM
12    CONTINUE

      DO 14 I = N,1,-1
         SUM = B(I)
         IF (I .LT. N) THEN
            DO 13 J = I+1, N
               SUM = SUM - A(I,J) * B(J)
13          CONTINUE
         END IF
         B(I) = SUM / A(I,I)
14    CONTINUE

      RETURN

      END SUBROUTINE LUBKSB
