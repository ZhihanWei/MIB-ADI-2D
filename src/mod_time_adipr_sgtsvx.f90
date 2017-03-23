      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !       PEACEMAN RACHAFORD ADI
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ADIPR(T,DT,UH)

      USE MOD_DATA

      REAL :: T,DT,UH(NY,NX)

      REAL :: UHS(NY,NX)
      INTEGER :: IX,IY,ierr

      !----- USED TO GENERATE D_XX AND D_YY ON THE LHS
      REAL,ALLOCATABLE :: AA(:),BB(:),CC(:),RR(:),UT(:)


      !----- USED BY LAPACK SUBROUTINE DGTSVX
      ! Purpose:
      !
      !   DGTSVX uses the LU factorization to compute the solution to a real system of linear equations A * X = B
      !   or A**T * X = B, where A is a tridiagonal matrix of order N and X and B are N-by-NRHS matrices.
      !
      !   Error bounds on the solution and a condition estimate are also provided.
      !
      ! Description:
      !
      !   The following steps are performed:
      !
      !   1. If FACT = 'N', the LU decomposition is used to factor the matrix A as A = L * U, where L is a product
      !      of permutation and unit lower bidiagonal matrices and U is upper triangular with nonzeros in only the
      !      main diagonal and first two superdiagonals.
      !
      !   2. If some U(i,i)=0, so that U is exactly singular, then the routine returns with INFO = i. Otherwise, the
      !      factored form of A is used to estimate the condition number of the matrix A. If the reciprocal of the
      !      condition number is less than machine precision, INFO = N+1 is returned as a warning, but the routine
      !      still goes on to solve for X and compute error bounds as described below.
      !
      !   3. The system of equations is solved for X using the factored form of A.
      !
      !   4. Iterative refinement is applied to improve the computed solution matrix and calculate error bounds and
      !      backward error estimates for it.
      !
      !   http://www.netlib.org/lapack/explore-html-3.4.2/d4/d62/group__double_g_tsolve.html#gaa9e938f737eedf395c4429393c769d07

      CHARACTER*1 :: FACT  = 'N'      ! [in] Specifies whether or not the factored form of A has been supplied on entry.

                                      ! = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored form of A;
                                      ! DL, D, DU, DLF, DF, DUF, DU2 and IPIV will not be modified.
                                      !
                                      ! = 'N':  The matrix will be copied to DLF, DF, and DUF and factored.

      CHARACTER*1 :: TRANS = 'N'      ! [in] Specifies the form of the system of equations:
                                      ! = 'N':  A * X = B     (No transpose)
                                      ! = 'T':  A**T * X = B  (Transpose)
                                      ! = 'C':  A**H * X = B  (Conjugate transpose = Transpose)

      INTEGER     :: N                ! [in] The order of the matrix A.  N >= 0.

      INTEGER     :: NRHS = 1         ! [in] The number of right hand sides, i.e., the number of columns of the matrix B.
                                      ! NRHS >= 0.

      REAL,ALLOCATABLE :: DL(:)       ! [in] DL is REAL array, dimension (N-1).
                                      ! The (n-1) subdiagonal elements of A.

      REAL,ALLOCATABLE :: D(:)        ! [in] D is REAL array, dimension (N)
                                      ! The n diagonal elements of A.

      REAL,ALLOCATABLE :: DU(:)       ! [in] DU is REAL array, dimension (N-1)
                                      ! The (n-1) superdiagonal elements of A.

      REAL,ALLOCATABLE :: DLF(:)      ! [in,out] DLF is REAL array, dimension (N-1)
                                      ! If FACT = 'F', then DLF is an input argument and on entry contains the (n-1) multipliers
                                      ! that define the matrix L from the LU factorization of A as computed by DGTTRF.
                                      !
                                      ! If FACT = 'N', then DLF is an output argument and on exit contains the (n-1) multipliers
                                      ! that define the matrix L from the LU factorization of A.

      REAL,ALLOCATABLE :: DF(:)       ! [in,out] DF is REAL array, dimension (N)
                                      ! If FACT = 'F', then DF is an input argument and on entry contains the n diagonal elements
                                      ! of the upper triangular matrix U from the LU factorization of A.
                                      !
                                      ! If FACT = 'N', then DF is an output argument and on exit contains the n diagonal elements
                                      ! of the upper triangular matrix U from the LU factorization of A.

      REAL,ALLOCATABLE :: DUF(:)      ! [in,out] DUF is REAL array, dimension (N-1)
                                      ! If FACT = 'F', then DUF is an input argument and on entry contains the (n-1) elements of
                                      ! the first superdiagonal of U.
                                      !
                                      ! If FACT = 'N', then DUF is an output argument and on exit contains the (n-1) elements of
                                      ! the first superdiagonal of U.

      REAL,ALLOCATABLE :: DU2(:)      ! [in,out] DU2 is REAL array, dimension (N-2)
                                      ! If FACT = 'F', then DU2 is an input argument and on entry contains the (n-2) elements of
                                      ! the second superdiagonal of U.
                                      !
                                      ! If FACT = 'N', then DU2 is an output argument and on exit contains the (n-2) elements of
                                      ! the second superdiagonal of U.

      INTEGER,ALLOCATABLE :: IPIV(:)  ! [in,out] IPIV is INTEGER array, dimension (N)
                                      ! If FACT = 'F', then IPIV is an input argument and on entry contains the pivot indices
                                      ! from the LU factorization of A as computed by DGTTRF.

                                      ! If FACT = 'N', then IPIV is an output argument and on exit contains the pivot indices
                                      ! from the LU factorization of A; row i of the matrix was interchanged with row IPIV(i).
                                      ! IPIV(i) will always be either i or i+1; IPIV(i) = i indicates a row interchange was
                                      ! not required.

      REAL,ALLOCATABLE :: B(:,:)      ! [in] B is REAL array, dimension (LDB,NRHS)
                                      ! The N-by-NRHS right hand side matrix B.

      INTEGER :: LDB                  ! [in] LDB is INTEGER
                                      ! The leading dimension of the array B.  LDB >= max(1,N).

      REAL,ALLOCATABLE :: X(:,:)      ! [out] X is REAL array, dimension (LDX,NRHS)
                                      ! If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X.

      INTEGER :: LDX                  ! [in] LDX is INTEGER
                                      ! The leading dimension of the array X.  LDX >= max(1,N).

      REAL :: RCOND                   ! [out] RCOND is REAL
                                      ! The estimate of the reciprocal condition number of the matrix A. If RCOND is less than
                                      ! the machine precision (in particular, if RCOND = 0), the matrix is singular to working
                                      ! precision. This condition is indicated by a return code of INFO > 0.

      REAL,ALLOCATABLE :: FERR(:)     ! [out] FERR is REAL array, dimension (NRHS)
                                      ! The estimated forward error bound for each solution vector X(j) (the j-th column of the
                                      ! solution matrix X). If XTRUE is the true solution corresponding to X(j), FERR(j) is an
                                      ! estimated upper bound for the magnitude of the largest element in (X(j) - XTRUE) divided
                                      ! by the magnitude of the largest element in X(j). The estimate is as reliable as the
                                      ! estimate for RCOND, and is almost always a slight overestimate of the true error.

      REAL,ALLOCATABLE :: BERR(:)     ! [out] BERR is REAL array, dimension (NRHS)
                                      ! The componentwise relative backward error of each solution vector X(j) (i.e., the
                                      ! smallest relative change in any element of A or B that makes X(j) an exact solution).

      REAL,ALLOCATABLE :: WORK(:)     ! [out] WORK is REAL array, dimension (3*N)

      INTEGER,ALLOCATABLE :: IWORK(:) ! [out] IWORK is INTEGER array, dimension (N)

      INTEGER :: INFO                 ! [out] INFO is INTEGER
                                      ! = 0:  successful exit
                                      ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                      ! > 0:  if INFO =  i, and i is
                                      !       <= N:  U(i,i) is exactly zero.  The factorization has not been completed unless
                                      !              i = N, but the factor U is exactly singular, so the solution and error
                                      !              bounds could not be computed. RCOND = 0 is returned.
                                      !       = N+1: U is nonsingular, but RCOND is less than machine precision, meaning that
                                      !              the matrix is singular to working precision. Nevertheless, the solution
                                      !              and error bounds are computed because there are a number of situations
                                      !              where the computed solution can be more accurate than the value of RCOND
                                      !              would suggest.


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
      ALLOCATE(AA(NX),BB(NX),CC(NX),RR(NX),UT(NX),STAT=ierr)

      N     = NX
      LDB   = NX
      LDX   = NX

      ALLOCATE(DL(N-1),D(N),DU(N-1),DLF(N-1),DF(N),DUF(N-1),DU2(N-2),IPIV(N),B(LDB,NRHS),X(LDX,NRHS),STAT=ierr)
      ALLOCATE(FERR(NRHS),BERR(NRHS),WORK(3*N),IWORK(N),STAT=ierr)

      DO IY = 2,NY-1
         CALL D_XX(IY,DT/2.0D0,UHS,AA,BB,CC,RR)  ! CONVERT LHS OPERATOR (1/BETA - DT/2 D_XX) TO A TRI-DIAGONAL MATRIX
                                                 ! AND STORE IT IN AA, BB, AND CC WITH CORRESPONDING RHS RR

         DL    = AA(2:NX)
         D     = BB
         DU    = CC(1:(NX-1))
         DLF   = 0.0D0
         DF    = 0.0D0
         DUF   = 0.0D0
         DU2   = 0.0D0
         IPIV  = 0
         DO IX = 1,NX
            B(IX,NRHS) = RR(IX)
         END DO
         X     = 0.0D0
         RCOND = 0.0D0
         FERR  = 0.0D0
         BERR  = 0.0D0
         WORK  = 0.0D0
         IWORK = 0
         INFO  = 1

         CALL SGTSVX(FACT,TRANS,N,NRHS,DL,D,DU,DLF,DF,DUF,DU2,IPIV,B,LDB,X,LDX,RCOND,FERR,BERR,WORK,IWORK,INFO)

         IF (INFO .NE. 0) STOP

         DO IX = 1,NX
            UHS(IY,IX) = X(IX,NRHS)
         END DO

!         WRITE(*,*) "IY = ",IY,", RCOND = ",1.0D0/RCOND

      END DO

      DEALLOCATE(AA,BB,CC,RR,UT,STAT=ierr)

      DEALLOCATE(DL,D,DU,DLF,DF,DUF,DU2,IPIV,B,X,FERR,BERR,WORK,IWORK,STAT=ierr)

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
      ALLOCATE(AA(NY),BB(NY),CC(NY),RR(NY),UT(NY),STAT=ierr)

      N     = NY
      LDB   = NY
      LDX   = NY

      ALLOCATE(DL(N-1),D(N),DU(N-1),DLF(N-1),DF(N),DUF(N-1),DU2(N-2),IPIV(N),B(LDB,NRHS),X(LDX,NRHS),STAT=ierr)
      ALLOCATE(FERR(NRHS),BERR(NRHS),WORK(3*N),IWORK(N),STAT=ierr)

      DO IX = 2,NX-1

         CALL D_YY(IX,DT/2.0D0,UH,AA,BB,CC,RR)   ! CONVERT LHS OPERATOR (1/BETA - DT/2 D_YY) TO A TRI-DIAGONAL MATRIX
                                                 ! AND STORE IT IN AA, BB, AND CC WITH CORRESPONDING RHS RR

         DL    = AA(2:NY)
         D     = BB
         DU    = CC(1:(NY-1))
         DLF   = 0.0D0
         DF    = 0.0D0
         DUF   = 0.0D0
         DU2   = 0.0D0
         IPIV  = 0
         DO IY = 1,NY
            B(IY,NRHS) = RR(IY)
         END DO
         X     = 0.0D0
         RCOND = 0.0D0
         FERR  = 0.0D0
         BERR  = 0.0D0
         WORK  = 0.0D0
         IWORK = 0
         INFO  = 1

         CALL SGTSVX(FACT,TRANS,N,NRHS,DL,D,DU,DLF,DF,DUF,DU2,IPIV,B,LDB,X,LDX,RCOND,FERR,BERR,WORK,IWORK,INFO)

         IF (INFO .NE. 0) STOP

         DO IY = 1,NY
            UH(IY,IX) = X(IY,NRHS)
         END DO

!         WRITE(*,*) "IX = ",IX,", RCOND = ",1.0D0/RCOND

!         IF (IX .EQ. 23) THEN
!            WRITE(*,*) "IPIV = ",IPIV
!            STOP
!         END IF
      END DO

      DEALLOCATE(AA,BB,CC,RR,UT,STAT=ierr)

      DEALLOCATE(DL,D,DU,DLF,DF,DUF,DU2,IPIV,B,X,FERR,BERR,WORK,IWORK,STAT=ierr)

      RETURN

      END SUBROUTINE ADIPR

