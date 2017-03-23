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


      !----- USED BY LAPACK SUBROUTINE SGTSV
      ! Purpose:
      !
      !   SGTSV  solves the equation
      !
      !      A*X = B,
      !
      !   where A is an n by n tridiagonal matrix, by Gaussian elimination with partial pivoting.
      !
      !   Note that the equation  A**T*X = B  may be solved by interchanging the order of the arguments DU and DL.
      !
      !   http://www.netlib.org/lapack/explore-html-3.4.2/d1/d88/group__real_g_tsolve.html#gae1cbb7cd9c376c9cc72575d472eba346

      INTEGER     :: N                ! [in] The order of the matrix A.  N >= 0.

      INTEGER     :: NRHS = 1         ! [in] The number of right hand sides, i.e., the number of columns of the matrix B.
                                      ! NRHS >= 0.

      REAL,ALLOCATABLE :: DL(:)       ! [in,out] DL is REAL array, dimension (N-1).
                                      ! On entry, DL must contain the (n-1) sub-diagonal elements of A.
                                      !
                                      ! On exit, DL is overwritten by the (n-2) elements of the second super-diagonal of
                                      ! the upper triangular matrix U from the LU factorization of A, in DL(1), ..., DL(n-2).

      REAL,ALLOCATABLE :: D(:)        ! [in,out] D is REAL array, dimension (N)
                                      ! On entry, D must contain the diagonal elements of A.
                                      !
                                      ! On exit, D is overwritten by the n diagonal elements of U.

      REAL,ALLOCATABLE :: DU(:)       ! [in,out] DU is REAL array, dimension (N-1)
                                      ! On entry, DU must contain the (n-1) super-diagonal elements of A.
                                      !
                                      ! On exit, DU is overwritten by the (n-1) elements of the first super-diagonal of U.

      REAL,ALLOCATABLE :: B(:,:)      ! [in,out] B is REAL array, dimension (LDB,NRHS)
                                      ! On entry, the N by NRHS matrix of right hand side matrix B.
                                      !
                                      ! On exit, if INFO = 0, the N by NRHS solution matrix X.

      INTEGER :: LDB                  ! [in] LDB is INTEGER
                                      ! The leading dimension of the array B.  LDB >= max(1,N).

      INTEGER :: INFO                 ! [out] INFO is INTEGER
                                      ! = 0:  successful exit
                                      ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                      ! > 0:  if INFO = i, U(i,i) is exactly zero, and the solution has not been computed.
                                      !       The factorization has not been completed unless i = N.


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
      ALLOCATE(DL(N-1),D(N),DU(N-1),B(LDB,NRHS),STAT=ierr)

      DO IY = 2,NY-1

         CALL D_XX(IY,DT/2.0D0,UHS,AA,BB,CC,RR)  ! CONVERT LHS OPERATOR (1/BETA - DT/2 D_XX) TO A TRI-DIAGONAL MATRIX
                                                 ! AND STORE IT IN AA, BB, AND CC WITH CORRESPONDING RHS RR

         DL    = AA(2:NX)
         D     = BB
         DU    = CC(1:(NX-1))
         DO IX = 1,LDB
            B(IX,NRHS) = RR(IX)
         END DO
         INFO  = 0

!         WRITE(*,*) "I AM HERE ..."
!         WRITE(*,*) "DL = ", DL
!         WRITE(*,*) "D  = ", D
!         WRITE(*,*) "DU = ", DU

         CALL SGTSV(N,NRHS,DL,D,DU,B,LDB,INFO)

         IF (INFO .NE. 0) STOP

         DO IX = 1,NX
            UHS(IY,IX) = B(IX,NRHS)
         END DO

      END DO

      DEALLOCATE(AA,BB,CC,RR,UT,STAT=ierr)

      DEALLOCATE(DL,D,DU,B,STAT=ierr)

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
      ALLOCATE(DL(N-1),D(N),DU(N-1),B(LDB,NRHS),STAT=ierr)

      DO IX = 2,NX-1

         CALL D_YY(IX,DT/2.0D0,UH,AA,BB,CC,RR)   ! CONVERT LHS OPERATOR (1/BETA - DT/2 D_YY) TO A TRI-DIAGONAL MATRIX
                                                 ! AND STORE IT IN AA, BB, AND CC WITH CORRESPONDING RHS RR

         DL    = AA(2:NY)
         D     = BB
         DU    = CC(1:(NY-1))
         DO IY = 1,LDB
            B(IY,NRHS) = RR(IY)
         END DO
         INFO  = 0

         CALL SGTSV(N,NRHS,DL,D,DU,B,LDB,INFO)

         IF (INFO .NE. 0) STOP

         DO IY = 1,NY
            UH(IY,IX) = B(IY,NRHS)
         END DO

      END DO

      DEALLOCATE(AA,BB,CC,RR,UT,STAT=ierr)

      DEALLOCATE(DL,D,DU,B,STAT=ierr)

      RETURN

      END SUBROUTINE ADIPR


