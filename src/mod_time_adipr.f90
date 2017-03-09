      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !       PEACEMAN RACHAFORD ADI
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ADIPR(T,DT,UH)

      USE MOD_DATA

      REAL :: UH(NY,NX),UHS(NY,NX),UHS2(NY,NX)
      REAL :: ROW(2,5),ROW2(3,6)
      REAL, ALLOCATABLE :: A(:),B(:),C(:),R(:),UT(:)
      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LIST_DATA) :: DATA,DATAL,DATAR

      REAL :: SUM,T,DT,UJP,UJPL,UJPR,BUXJP,BUXJPR,BUXJPL
      INTEGER :: IX,IY,J,K,I,ierr
      REAL :: VDE1(-1:1),VDE2(-1:1)
      REAL,ALLOCATABLE :: WTYPE1(:,:),WTYPE2(:,:)       ! FOR IRREGULAR INTERFACE
      REAL :: WTYPE(-1:1),WTYPEC1(-1:1),WTYPEC2(-1:1)                                     

      !INITIALIZE GEBERAL WEIGHTS OF 1ST&2ND ORDER  
      CALL GETWTYPE(0,1,VDE1,1)
      CALL GETWTYPE(0,2,VDE2,1)

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !
      !FIRST STEP:
      !(1/BETA - DT/2 D_XX - DT/2 BETA_X/BETA D_X) U*
      !                    = (1/BETA + DT/2 D_YY + DT/2 BETA_Y/BETA D_Y)U^N + DT/2 SRC^{N+1/2}/BETA
      !                    = DT/2 D_YY U^N + (U^N + DT/2 SRC^{N+1/2})/BETA + DT/2 BETA_Y/BETA D_Y U^N
      !
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SETSRC(T+DT/2.0D0)  !SET SOURCE SRC(IY,IX) AT TIME T_{N+1/2}

      ALLOCATE(A(NX),B(NX),C(NX),R(NX),UT(NX),STAT=ierr)

      DO IX=2,NX-1 !GENERATE RHS

         !----- STEP I:   FIRST PART OF RHS: DT/2 D_YY U^N -----!
         !----- STEP I-1: FIRST PART OF RHS: DT/2 D_YY U^N WITHOUT MIB -----!
         DO IY=2,NY-1 !INITIAL WITHOUT MIB
            SUM=0.0D0
            DO J=-1,1
               SUM=SUM+UH(IY+J,IX)*VDE2(J)
            END DO            
            UHS(IY,IX)=SUM*DT/2.0D0
         END DO
         
         !----- STEP I-2: ADD MIB TO THE FRIST PART OF RHS -----!
         IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN
            ELEM => IFPX(IX)%HEAD
            DO WHILE ( ASSOCIATED(ELEM) )

               !----- AN IRREGULAR INTERFACE POINT                
               IF (ELEM%DATA%ID .GT. 0) THEN
                  DATA=LIST_GET_DATA(ELEM)

                  UJP=DATA%JUMP(1)
                  BUXJP=DATA%JUMP(6)

                  ALLOCATE(WTYPE2(2,-1:1),STAT=ierr)
                  CALL GETWTYPE(DATA%ITYPE,2,WTYPE2,2)
                  
                  IY=DATA%LGRD       !APPROXIMATE IY
                  
                  !IY+1 CROSS INTERFACE, USE FP2   
                  SUM=0.0D0                
                  DO J=-1,0
                     SUM=SUM+WTYPE2(2,J)*UH(IY+J,IX)
                  END DO
                  DO K=1,4          
                     SUM=SUM+WTYPE2(2,1)*UH(DATA%LGRD-2+K,IX)*DATA%WIJ(2,K)
                  END DO
                  SUM=SUM+WTYPE2(2,1)*(UJP*DATA%WIJ(2,5)+BUXJP*DATA%WIJ(2,6))
                  UHS(IY,IX)=SUM*DT/2.0D0                

                  IY=IY+1           !APPROXIMATE IY+1
                  !IY-1 CROSS INTERFACE, USE FP1
                  SUM=0.0D0
                  DO J=0,1
                     SUM=SUM+WTYPE2(1,J)*UH(IY+J,IX)
                  END DO
                  DO K=1,4          
                     SUM=SUM+WTYPE2(1,-1)*UH(DATA%LGRD-2+K,IX)*DATA%WIJ(1,K)
                  END DO
                  SUM=SUM+WTYPE2(1,-1)*(UJP*DATA%WIJ(1,5)+BUXJP*DATA%WIJ(1,6))
                  UHS(IY,IX)=SUM*DT/2.0D0  

                  ELEM => ELEM%NEXT     ! NEXT INTERFACE 

               !----- A CORNER POINT
               ELSE 
                  DATAL=LIST_GET_DATA(ELEM)
                  DATAR=LIST_GET_DATA(ELEM%NEXT)

                  UJPL=DATAL%JUMP(1)
                  UJPR=DATAR%JUMP(1)
                  BUXJPL=DATAL%JUMP(6)
                  BUXJPR=DATAR%JUMP(6)

                  ALLOCATE(WTYPE2(2,-1:1),STAT=ierr)

                  IY=DATAL%LGRD  !APPROXIMATE IY
                  CALL GETWTYPE(DATAL%ITYPE,2,WTYPE2,2)                 
                  DO J=-1,1
                     WTYPE(J)=WTYPE2(2,J)
                  END DO
                  !IY+1 CROSS INTERFACE, USE FP2 
                  SUM=0.0D0
                  DO J=-1,0
                     SUM=SUM+WTYPE(J)*UH(IY+J,IX)
                  END DO
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     SUM=SUM+WTYPE(1)*UH(DATAL%LGRD-2+K,IX)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                   +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UHS(IY,IX)=SUM*DT/2.0D0
                  
                  IY=IY+1       !APPROXIMATE IY+1
                  !DETERMINE USE WHICH INTERFACE DATA APPROXIMATE THE MIDDLE POINT, IY+1
                  IF (ABS(DATAL%ITYPE) .LT. 2) THEN              !USE LEFT INTERFACE WEIGHTS
                     CALL GETWTYPE(DATAL%ITYPE,2,WTYPE2,2)    
                     DO J=-1,1
                        WTYPE(J)=WTYPE2(1,J)
                     END DO            
                  ELSE                                           !USE RIGHT INTERFACE WEIGHTS
                     CALL GETWTYPE(DATAR%ITYPE,2,WTYPE2,2)    
                     DO J=-1,1
                        WTYPE(J)=WTYPE2(2,J)
                     END DO        
                  END IF                                
                  !IY-1 CROSS INTERFACE, USE FP1 
                  SUM=0.0D0                  
                  I=1                              !USE FP1, INDEX USED FOR FP WIJ2(:,:)                 
                  DO K=1,5
                     SUM=SUM+WTYPE(-1)*UH(DATAL%LGRD-2+K,IX)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(-1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                    +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  !IY, NO FP
                  SUM=SUM+WTYPE(0)*UH(IY,IX)
                  !IY+1 CROSS INTERFACE, USE FP3                  
                  I=I+2                 
                  DO K=1,5
                     SUM=SUM+WTYPE(1)*UH(DATAL%LGRD-2+K,IX)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                   +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UHS(IY,IX)=SUM*DT/2.0D0

                  IY=IY+1        !APPROXIMATE IY+2
                  CALL GETWTYPE(DATAR%ITYPE,2,WTYPE2,2) 
                  DO J=-1,1
                     WTYPE(J)=WTYPE2(1,J)
                  END DO
                  !IY-1 CROSS INTERFACE, USE FP2
                  SUM=0.0D0
                  DO J=0,1
                     SUM=SUM+WTYPE(J)*UH(IY+J,IX)
                  END DO
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     SUM=SUM+WTYPE(-1)*UH(DATAL%LGRD-2+K,IX)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(-1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                    +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UHS(IY,IX)=SUM*DT/2.0D0                  

                  ELEM => ELEM%NEXT%NEXT    ! SKIP THE COUPLE OF CORNER INTERFACES
               
               END IF 

               DEALLOCATE(WTYPE2,STAT=ierr)
               
            END DO

         END IF

         !----- STEP II: ADD SECOND PART OF RHS: +(U^N+DT/2 SRC^{N+1/2})/BETA -----!
         DO IY=2,NY-1
            UHS(IY,IX)=UHS(IY,IX)+(UH(IY,IX)+DT/2.0D0*SRC(IY,IX))/BETA(IY,IX)
         END DO

         !----- STEP III: ADD THIRD PART OF RHS: + DT/2 BETA_Y/BETA D_Y U^N -----!
         UHS2(:,IX) = 0.0D0 !SAVE RHS VALUE OBTAINED IN PREVIOUS TWO STEPS

         !----- STEP III-1: THIRD PART OF RHS: + DT/2 BETA_Y/BETA D_Y U^N WITHOUT MIB -----!
         DO IY=2,NY-1 !INITIALLY WITHOUT MIB
            SUM=0.0D0
            DO J=-1,1
               SUM=SUM+UH(IY+J,IX)*VDE1(J)
            END DO
            UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAY(IY,IX)

         END DO

         !----- STEP III-2: ADD MIB TO THE THIRD PART OF RHS -----!
         IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN
            ELEM => IFPX(IX)%HEAD
            DO WHILE ( ASSOCIATED(ELEM) )

               !----- AN IRREGULAR INTERFACE POINT                
               IF (ELEM%DATA%ID .GT. 0) THEN
                  DATA=LIST_GET_DATA(ELEM)

                  UJP=DATA%JUMP(1)
                  BUXJP=DATA%JUMP(6)

                  ALLOCATE(WTYPE1(2,-1:1),STAT=ierr)
                  CALL GETWTYPE(DATA%ITYPE,1,WTYPE1,2)

                  IY=DATA%LGRD  !APPROXIMATE IY
                  !IY+1 CROSS INTERFACE, USE FP2
                  SUM=0.0D0                   
                  DO J=-1,0
                     SUM=SUM+WTYPE1(2,J)*UH(IY+J,IX)
                  END DO
                  DO K=1,4          
                     SUM=SUM+WTYPE1(2,1)*UH(DATA%LGRD-2+K,IX)*DATA%WIJ(2,K)
                  END DO
                  SUM=SUM+WTYPE1(2,1)*(UJP*DATA%WIJ(2,5)+BUXJP*DATA%WIJ(2,6))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAY(IY,IX)
                  
                  IY=IY+1           !APPROXIMATE IY+1
                  !IY-1 CROSS INTERFACE, USE FP1
                  SUM=0.0D0
                  DO J=0,1
                     SUM=SUM+UH(IY+J,IX)*WTYPE1(1,J)
                  END DO
                  DO K=1,4          
                     SUM=SUM+WTYPE1(1,-1)*UH(DATA%LGRD-2+K,IX)*DATA%WIJ(1,K)
                  END DO
                  SUM=SUM+WTYPE1(1,-1)*(UJP*DATA%WIJ(1,5)+BUXJP*DATA%WIJ(1,6))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAY(IY,IX)

                  ELEM => ELEM%NEXT         ! NEXT INTERFACE

               !----- A CORNER POINT
               ELSE
                  DATAL=LIST_GET_DATA(ELEM)
                  DATAR=LIST_GET_DATA(ELEM%NEXT)

                  UJPL=DATAL%JUMP(1)
                  UJPR=DATAR%JUMP(1)
                  BUXJPL=DATAL%JUMP(6)
                  BUXJPR=DATAR%JUMP(6)

                  ALLOCATE(WTYPE1(2,-1:1),STAT=ierr)

                  IY=DATAL%LGRD  !APPROXIMATE IY
                  CALL GETWTYPE(DATAL%ITYPE,2,WTYPE1,2)                 
                  DO J=-1,1
                     WTYPE(J)=WTYPE1(2,J)
                  END DO
                  !IY+1 CROSS INTERFACE, USE FP2 
                  SUM=0.0D0
                  DO J=-1,0
                     SUM=SUM+WTYPE(J)*UH(IY+J,IX)
                  END DO
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     SUM=SUM+WTYPE(1)*UH(DATAL%LGRD-2+K,IX)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                   +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAY(IY,IX)
                  
                  IY=IY+1       !APPROXIMATE IY+1
                  !DETERMINE USE WHICH INTERFACE DATA APPROXIMATE THE MIDDLE POINT, IY+1
                  IF (ABS(DATAL%ITYPE) .LT. 2) THEN             !USE LEFT INTERFACE WEIGHTS 
                     CALL GETWTYPE(DATAL%ITYPE,2,WTYPE1,2)    
                     DO J=-1,1
                        WTYPE(J)=WTYPE1(1,J)
                     END DO            
                  ELSE                                          !USE RIGHT INTERFACE WEIGHTS 
                     CALL GETWTYPE(DATAR%ITYPE,2,WTYPE1,2)    
                     DO J=-1,1
                        WTYPE(J)=WTYPE1(2,J)
                     END DO        
                  END IF                                
                  !IY-1 CROSS INTERFACE, USE FP1 
                  SUM=0.0D0                  
                  I=1                              !USE FP1, INDEX USED FOR FP WIJ2(:,:)                  
                  DO K=1,5
                     SUM=SUM+WTYPE(-1)*UH(DATAL%LGRD-2+K,IX)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(-1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                    +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  !IY, NO FP
                  SUM=SUM+WTYPE(0)*UH(IY,IX)
                  !IY+1 CROSS INTERFACE, USE FP3                  
                  I=I+2                 
                  DO K=1,5
                     SUM=SUM+WTYPE(1)*UH(DATAL%LGRD-2+K,IX)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                   +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAY(IY,IX)

                  IY=IY+1        !APPROXIMATE IY+2
                  CALL GETWTYPE(DATAR%ITYPE,2,WTYPE1,2) 
                  DO J=-1,1
                     WTYPE(J)=WTYPE1(1,J)
                  END DO
                  !IY-1 CROSS INTERFACE, USE FP2 
                  SUM=0.0D0
                  DO J=0,1
                     SUM=SUM+WTYPE(J)*UH(IY+J,IX)
                  END DO
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     SUM=SUM+WTYPE(-1)*UH(DATAL%LGRD-2+K,IX)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(-1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                    +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAY(IY,IX)               

                  ELEM => ELEM%NEXT%NEXT    ! SKIP THE COUPLE OF CORNER INTERFACES

               END IF

               DEALLOCATE(WTYPE1,STAT=ierr)               

            END DO
         END IF 
   
         UHS(:,IX) = UHS(:,IX)+UHS2(:,IX) !ADD THE THIRD PART OF RHS

      END DO !----- END OF DO IX=2,NX-1 TO GENERATE RHS   

      !----- STEP IV: SET FIRST ORDER BOUNDARY CONDITION FOR N+1/2-----!
      CALL SETBC(T,DT,UHS)

      DO IY=2,NY-1

         !----- STEP V:   GENERATE LHS MATRIX -----!
         !----- STEP V-1: GENERATE LHS MATRIX: (1/BETA-DT/2*D_XX-DT/2*BETA_X/BETA*D_X) U* WITHOUT MIB -----!
         A(1)=0.0D0; B(1)=1.0D0; C(1)=0.0D0 
         DO IX=2,NX-1
            A(IX) = -(VDE2(-1)*DT/2.0D0+VDE1(-1)*DT/2.0D0*BETAX(IY,IX))
            B(IX) = 1.0D0/BETA(IY,IX)-(VDE2(0)*DT/2.0D0+VDE1(0)*DT/2.0D0*BETAX(IY,IX))
            C(IX) = -(VDE2(1)*DT/2.0D0+VDE1(1)*DT/2.0D0*BETAX(IY,IX))
         END DO
         A(NX)=0.0D0; B(NX)=1.0D0; C(NX)=0.0D0

         DO IX=1,NX  !R IS FOR RHS
            R(IX)=UHS(IY,IX)
         END DO
         
         !----- STEP V-2: GENERATE LHS MATRIX WITH MIB -----!
         IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN
            ELEM => IFPY(IY)%HEAD
            DO WHILE ( ASSOCIATED(ELEM) )

               !----- AN IRREGULAR INTERFACE POINT                
               IF (ELEM%DATA%ID .GT. 0) THEN
                  DATA=LIST_GET_DATA(ELEM)

                  UJP=DATA%JUMP(1)
                  BUXJP=DATA%JUMP(6)

                  ALLOCATE(WTYPE1(2,-1:1),WTYPE2(2,-1:1),STAT=ierr)
                  CALL GETWTYPE(DATA%ITYPE,1,WTYPE1,2)
                  CALL GETWTYPE(DATA%ITYPE,2,WTYPE2,2)
                  
                  ROW=0.0D0
                  IX=DATA%LGRD !APPROXIMATE IX
                  !USE FP2
                  ROW(1,1)=-(WTYPE2(2,-1)*DT/2.0D0+WTYPE1(2,-1)*DT/2.0D0*BETAX(IY,IX))
                  ROW(1,2)=1.0D0/BETA(IY,IX)-(WTYPE2(2,0)*DT/2.0D0+WTYPE1(2,0)*DT/2.0D0*BETAX(IY,IX))
                  DO K=1,4
                     ROW(1,K)=ROW(1,K)-(WTYPE2(2,1)*DT/2.0D0+WTYPE1(2,1)*DT/2.0D0*BETAX(IY,IX))*DATA%WIJ(2,K)
                  END DO
                  ROW(1,5)=UHS(IY,IX)+(WTYPE2(2,1)*DT/2.0D0+WTYPE1(2,1)*DT/2.0D0*BETAX(IY,IX))*&
                           (UJP*DATA%WIJ(2,5)+BUXJP*DATA%WIJ(2,6))

                  IX=IX+1 !APPROXIMATE IX+1
                  !USE FP1
                  ROW(2,3)=1.0D0/BETA(IY,IX)-(WTYPE2(1,0)*DT/2.0D0+WTYPE1(1,0)*DT/2.0D0*BETAX(IY,IX))
                  ROW(2,4)=-(WTYPE2(1,1)*DT/2.0D0+WTYPE1(1,1)*DT/2.0D0*BETAX(IY,IX))
                  DO K=1,4
                     ROW(2,K)=ROW(2,K)-(WTYPE2(1,-1)*DT/2.0D0+WTYPE1(1,-1)*DT/2.0D0*BETAX(IY,IX))*DATA%WIJ(1,K)
                  END DO
                  ROW(2,5)=UHS(IY,IX)+(WTYPE2(1,-1)*DT/2.0D0+WTYPE1(1,-1)*DT/2.0D0*BETAX(IY,IX))*&
                           (UJP*DATA%WIJ(1,5)+BUXJP*DATA%WIJ(1,6))

                  CALL GE1(NX,A,B,C,R,DATA%LGRD,ROW)

                  ELEM => ELEM%NEXT         ! NEXT INTERFACE

               !----- A CORNER POINT
               ELSE
                  DATAL=LIST_GET_DATA(ELEM)
                  DATAR=LIST_GET_DATA(ELEM%NEXT)

                  UJPL=DATAL%JUMP(1)
                  UJPR=DATAR%JUMP(1)
                  BUXJPL=DATAL%JUMP(6)
                  BUXJPR=DATAR%JUMP(6)

                  ALLOCATE(WTYPE1(2,-1:1),WTYPE2(2,-1:1),STAT=ierr)

                  ROW2=0.0D0

                  IX=DATAL%LGRD   !APPROXIMATE IX
                  CALL GETWTYPE(DATAL%ITYPE,1,WTYPE1,2)
                  CALL GETWTYPE(DATAL%ITYPE,2,WTYPE2,2)
                  DO J=-1,1
                     WTYPEC1(J)=WTYPE1(2,J)
                     WTYPEC2(J)=WTYPE2(2,J)
                  END DO
                  !IX+1 CROSS INTERFACE, USE FP2 
                  ROW2(1,1)=-(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAX(IY,IX))
                  ROW2(1,2)=1.0D0/BETA(IY,IX)-(WTYPEC2(0)*DT/2.0D0+WTYPEC1(0)*DT/2.0D0*BETAX(IY,IX))
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     ROW2(1,K)=ROW2(1,K)-(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAX(IY,IX))*DATAL%WIJ2(I,K)
                  END DO
                  ROW2(1,6)=UHS(IY,IX)+(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAX(IY,IX))*&
                                       (UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)+&
                                        UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))

                  IX=IX+1        !APPROXIMATE IX+1
                  !DETERMINE USE WHICH INTERFACE DATA APPROXIMATE THE MIDDLE POINT, IY+1
                  IF (ABS(DATAL%ITYPE) .LT. 2) THEN            !USE LEFT INTERFACE WEIGHTS
                     CALL GETWTYPE(DATAL%ITYPE,1,WTYPE1,2)  
                     CALL GETWTYPE(DATAL%ITYPE,2,WTYPE2,2)
                     DO J=-1,1
                        WTYPEC1(J)=WTYPE1(1,J)
                        WTYPEC2(J)=WTYPE2(1,J)
                     END DO
                  ELSE                                         !USE RIGHT INTERFACE WEIGHTS
                     CALL GETWTYPE(DATAR%ITYPE,1,WTYPE1,2)    
                     CALL GETWTYPE(DATAR%ITYPE,2,WTYPE2,2)
                     DO J=-1,1
                        WTYPEC1(J)=WTYPE1(2,J)
                        WTYPEC2(J)=WTYPE2(2,J)
                     END DO
                  END IF
                  ROW2(2,3)=1.0D0/BETA(IY,IX)-(WTYPEC2(0)*DT/2.0D0+WTYPEC1(0)*DT/2.0D0*BETAX(IY,IX))
                  !IX-1 CROSS INTERFACE, USE FP1
                  I=1                              !USE FP1, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     ROW2(2,K)=ROW2(2,K)-(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAX(IY,IX))*DATAL%WIJ2(I,K)
                  END DO
                  ROW2(2,6)=UHS(IY,IX)+(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAX(IY,IX))*&
                                       (UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)+&
                                        UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))              

                  !IX+1 CROSS INTERFACE, USE FP3 
                  I=I+2
                  DO K=1,5
                     ROW2(2,K)=ROW2(2,K)-(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAX(IY,IX))*DATAL%WIJ2(I,K)
                  END DO
                  ROW2(2,6)=ROW2(2,6)+(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAX(IY,IX))*&
                                     (UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)+&
                                      UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9)) 

                  IX=IX+1        !APPROXIMATE IX+2
                  CALL GETWTYPE(DATAR%ITYPE,1,WTYPE1,2)
                  CALL GETWTYPE(DATAR%ITYPE,2,WTYPE2,2)
                  DO J=-1,1
                     WTYPEC1(J)=WTYPE1(1,J)
                     WTYPEC2(J)=WTYPE2(1,J)
                  END DO
                  !IX-1 CROSS INTERFACE, USE FP2 
                  ROW2(3,4)=1.0D0/BETA(IY,IX)-(WTYPEC2(0)*DT/2.0D0+WTYPEC1(0)*DT/2.0D0*BETAX(IY,IX))
                  ROW2(3,5)=-(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAX(IY,IX))                  
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     ROW2(3,K)=ROW2(3,K)-(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAX(IY,IX))*DATAL%WIJ2(I,K)
                  END DO
                  ROW2(3,6)=UHS(IY,IX)+(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAX(IY,IX))*&
                                       (UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)+&
                                        UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))

                  CALL GE2(NX,A,B,C,R,DATAL%LGRD,ROW2)

                  ELEM => ELEM%NEXT%NEXT    ! SKIP THE COUPLE OF CORNER INTERFACES

               END IF 
               
               DEALLOCATE(WTYPE1,WTYPE2,STAT=ierr)

            END DO
         END IF
         
         !----- STEP VI: THOMAS ALGORITHM -----!
         CALL TRIDAG(A,B,C,R,UT,NX)

         DO IX=1,NX
            UHS(IY,IX)=UT(IX)
         END DO       

      END DO         !END OF DO IY=2,NY-1 TO GENERATE RHS MATRIX

      DEALLOCATE(A,B,C,R,UT,STAT=ierr)

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !
      !SECOND STEP:
      !(1/BETA - DT/2 D_Y^2 - DT/2 BETA_Y/BETA D_Y) U^(N+1)
      !                     = (1/BETA + DT/2 D_X^2 + DT/2 BETA_X/BETA D_X) U* + DT/2 SRC^{N+1/2}/BETA
      !                     = DT/2 D_XX U* + (U* + DT/2 SRC^{N+1/2})/BETA + DT/2 BETA_X/BETA D_X U*
      !
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ALLOCATE(A(NY),B(NY),C(NY),R(NY),UT(NY),STAT=ierr) 

      DO IY=2,NY-1

         !----- STEP I:   FIRST PART OF RHS: DT/2 D_XX U*  -----!
         !----- STEP I-1: FIRST PART OF RHS: DT/2 D_XX U* WITHOUT MIB -----!
         DO IX=2,NX-1  !INITLIAL WITHOUT MIB 
            SUM=0.0D0
            DO J=-1,1
               SUM=SUM+UHS(IY,IX+J)*VDE2(J)
            END DO
            UH(IY,IX)=SUM*DT/2.0D0
         END DO

         !----- STEP I-2: ADD MIB TO THE FRIST PART OF RHS -----!         
         IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN
            ELEM =>IFPY(IY)%HEAD
            DO WHILE ( ASSOCIATED(ELEM) )

               !----- AN IRREGULAR INTERFACE POINT                
               IF (ELEM%DATA%ID .GT. 0) THEN
                  DATA=LIST_GET_DATA(ELEM)

                  UJP=DATA%JUMP(1)
                  BUXJP=DATA%JUMP(6)

                  ALLOCATE(WTYPE2(2,-1:1),STAT=ierr)
                  CALL GETWTYPE(DATA%ITYPE,2,WTYPE2,2)
                  
                  IX=DATA%LGRD      !APPROXIMATE IX
                  !IX+1 CROSS INTERFACE, USE FP2   
                  SUM=0.0D0                
                  DO J=-1,0
                     SUM=SUM+WTYPE2(2,J)*UHS(IY,IX+J)
                  END DO
                  DO K=1,4          
                     SUM=SUM+WTYPE2(2,1)*UHS(IY,DATA%LGRD-2+K)*DATA%WIJ(2,K)
                  END DO
                  SUM=SUM+WTYPE2(2,1)*(UJP*DATA%WIJ(2,5)+BUXJP*DATA%WIJ(2,6))
                  UH(IY,IX)=SUM*DT/2.0D0
   
                  IX=IX+1           !APPROXIMATE IX+1
                  !IX-1 CROSS INTERFACE, USE FP1
                  SUM=0.0D0
                  DO J=0,1
                     SUM=SUM+WTYPE2(1,J)*UHS(IY,IX+J)
                  END DO
                  DO K=1,4          
                     SUM=SUM+WTYPE2(1,-1)*UHS(IY,DATA%LGRD-2+K)*DATA%WIJ(1,K)
                  END DO
                  SUM=SUM+WTYPE2(1,-1)*(UJP*DATA%WIJ(1,5)+BUXJP*DATA%WIJ(1,6))
                  UH(IY,IX)=SUM*DT/2.0D0
                  
                  ELEM => ELEM%NEXT     ! NEXT INTERFACE

               !----- A CORNER POINT
               ELSE 
                  DATAL=LIST_GET_DATA(ELEM)
                  DATAR=LIST_GET_DATA(ELEM%NEXT)

                  UJPL=DATAL%JUMP(1)
                  UJPR=DATAR%JUMP(1)
                  BUXJPL=DATAL%JUMP(6)
                  BUXJPR=DATAR%JUMP(6)

                  ALLOCATE(WTYPE2(2,-1:1),STAT=ierr)

                  IX=DATAL%LGRD  !APPROXIMATE IX
                  CALL GETWTYPE(DATAL%ITYPE,2,WTYPE2,2)                 
                  DO J=-1,1
                     WTYPE(J)=WTYPE2(2,J)
                  END DO
                  !IX+1 CROSS INTERFACE, USE FP2 
                  SUM=0.0D0
                  DO J=-1,0
                     SUM=SUM+WTYPE(J)*UHS(IY,IX+J)
                  END DO
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     SUM=SUM+WTYPE(1)*UHS(IY,DATAL%LGRD-2+K)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                   +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UH(IY,IX)=SUM*DT/2.0D0
                  
                  IX=IX+1       !APPROXIMATE IX+1
                  !DETERMINE USE WHICH INTERFACE DATA APPROXIMATE THE MIDDLE POINT, IX+1
                  IF (ABS(DATAL%ITYPE) .LT. 2) THEN             !USE LEFT INTERFACE WEIGHTS 
                     CALL GETWTYPE(DATAL%ITYPE,2,WTYPE2,2)    
                     DO J=-1,1
                        WTYPE(J)=WTYPE2(1,J)
                     END DO            
                  ELSE                                          !USE RIGHT INTERFACE WEIGHTS
                     CALL GETWTYPE(DATAR%ITYPE,2,WTYPE2,2)    
                     DO J=-1,1
                        WTYPE(J)=WTYPE2(2,J)
                     END DO        
                  END IF                                
                  !IX-1 CROSS INTERFACE, USE FP1 
                  SUM=0.0D0                  
                  I=1                              !USE FP1, INDEX USED FOR FP WIJ2(:,:)                  
                  DO K=1,5
                     SUM=SUM+WTYPE(-1)*UHS(IY,DATAL%LGRD-2+K)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(-1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                    +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  !IX, NO FP
                  SUM=SUM+WTYPE(0)*UHS(IY,IX)
                  !IX+1 CROSS INTERFACE, USE FP3                  
                  I=I+2                 
                  DO K=1,5
                     SUM=SUM+WTYPE(1)*UHS(IY,DATAL%LGRD-2+K)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                   +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UH(IY,IX)=SUM*DT/2.0D0

                  IX=IX+1        !APPROXIMATE IX+2
                  CALL GETWTYPE(DATAR%ITYPE,2,WTYPE2,2) 
                  DO J=-1,1
                     WTYPE(J)=WTYPE2(1,J)
                  END DO
                  !IX-1 CROSS INTERFACE, USE FP2 
                  SUM=0.0D0
                  DO J=0,1
                     SUM=SUM+WTYPE(J)*UHS(IY,IX+J)
                  END DO
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     SUM=SUM+WTYPE(-1)*UHS(IY,DATAL%LGRD-2+K)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(-1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                    +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UH(IY,IX)=SUM*DT/2.0D0                  

                  ELEM => ELEM%NEXT%NEXT    ! SKIP THE COUPLE OF CORNER INTERFACES
               
               END IF 

               DEALLOCATE(WTYPE2,STAT=ierr)
               
            END DO 

         END IF

         !----- STEP II: ADD SECOND PART OF RHS: +(U*+DT/2 SRC^{N+1/2})/BETA -----!
         DO IX=2,NX-1
            UH(IY,IX)=UH(IY,IX)+(UHS(IY,IX)+DT/2.0D0*SRC(IY,IX))*1.0D0/BETA(IY,IX)
         END DO

         !----- STEP III: ADD THIRD PART OF RHS: + DT/2 BETA_X/BETA D_X U* -----!
         UHS2(IY,:) = 0.0D0 !SAVE RHS VALUE OBTAINED IN PREVIOUS TWO STEPS

         !----- STEP III-1: THIRD PART OF RHS: + DT/2 BETA_X/BETA D_X U* WITHOUT MIB -----!
         DO IX=2,NX-1 !INITIALLY WITHOUT MIB
            SUM=0.0D0
            DO J=-1,1
               SUM=SUM+UHS(IY,IX+J)*VDE1(J)
            END DO
            UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAX(IY,IX)
         END DO

         !----- STEP III-2: ADD MIB TO THE THIRD PART OF RHS -----!
         IF ( ASSOCIATED(IFPY(IY)%HEAD) ) THEN
            ELEM => IFPY(IY)%HEAD
            DO WHILE ( ASSOCIATED(ELEM) )

               !----- AN IRREGULAR INTERFACE POINT                
               IF (ELEM%DATA%ID .GT. 0) THEN
                  DATA=LIST_GET_DATA(ELEM)

                  UJP=DATA%JUMP(1)
                  BUXJP=DATA%JUMP(6)

                  ALLOCATE(WTYPE1(2,-1:1),STAT=ierr)
                  CALL GETWTYPE(DATA%ITYPE,1,WTYPE1,2)

                  IX=DATA%LGRD  !APPROXIMATE IX
                  !IX+1 CROSS INTERFACE, USE FP2
                  SUM=0.0D0                   
                  DO J=-1,0
                     SUM=SUM+WTYPE1(2,J)*UHS(IY,IX+J)
                  END DO
                  DO K=1,4          
                     SUM=SUM+WTYPE1(2,1)*UHS(IY,DATA%LGRD-2+K)*DATA%WIJ(2,K)
                  END DO
                  SUM=SUM+WTYPE1(2,1)*(UJP*DATA%WIJ(2,5)+BUXJP*DATA%WIJ(2,6))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAX(IY,IX)
                  
                  IX=IX+1           !APPROXIMATE IX+1
                  !IX-1 CROSS INTERFACE, USE FP1
                  SUM=0.0D0
                  DO J=0,1
                     SUM=SUM+UHS(IY,IX+J)*WTYPE1(1,J)
                  END DO
                  DO K=1,4          
                     SUM=SUM+WTYPE1(1,-1)*UHS(IY,DATA%LGRD-2+K)*DATA%WIJ(1,K)
                  END DO
                  SUM=SUM+WTYPE1(1,-1)*(UJP*DATA%WIJ(1,5)+BUXJP*DATA%WIJ(1,6))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAX(IY,IX)

                  ELEM => ELEM%NEXT         ! NEXT INTERFACE

               !----- A CORNER POINT
               ELSE
                  DATAL=LIST_GET_DATA(ELEM)
                  DATAR=LIST_GET_DATA(ELEM%NEXT)

                  UJPL=DATAL%JUMP(1)
                  UJPR=DATAR%JUMP(1)
                  BUXJPL=DATAL%JUMP(6)
                  BUXJPR=DATAR%JUMP(6)

                  ALLOCATE(WTYPE1(2,-1:1),STAT=ierr)

                  IX=DATAL%LGRD  !APPROXIMATE IX
                  CALL GETWTYPE(DATAL%ITYPE,2,WTYPE1,2)                 
                  DO J=-1,1
                     WTYPE(J)=WTYPE1(2,J)
                  END DO
                  !IX+1 CROSS INTERFACE, USE FP2
                  SUM=0.0D0
                  DO J=-1,0
                     SUM=SUM+WTYPE(J)*UHS(IY,IX+J)
                  END DO
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     SUM=SUM+WTYPE(1)*UHS(IY,DATAL%LGRD-2+K)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                   +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAX(IY,IX)
                  
                  IX=IX+1       !APPROXIMATE IX+1
                  !DETERMINE USE WHICH INTERFACE DATA APPROXIMATE THE MIDDLE POINT, IX+1
                  IF (ABS(DATAL%ITYPE) .LT. 2) THEN             !USE LEFT INTERFACE WEIGHTS  
                     CALL GETWTYPE(DATAL%ITYPE,2,WTYPE1,2)    
                     DO J=-1,1
                        WTYPE(J)=WTYPE1(1,J)
                     END DO            
                  ELSE                                          !USE RIGHT INTERFACE WEIGHTS 
                     CALL GETWTYPE(DATAR%ITYPE,2,WTYPE1,2)    
                     DO J=-1,1
                        WTYPE(J)=WTYPE1(2,J)
                     END DO        
                  END IF                                
                  !IX-1 CROSS INTERFACE, USE FP1 
                  SUM=0.0D0                  
                  I=1                              !USE FP1, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     SUM=SUM+WTYPE(-1)*UHS(IY,DATAL%LGRD-2+K)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(-1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                    +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  !IX, NO FP
                  SUM=SUM+WTYPE(0)*UHS(IY,IX)
                  !IX+1 CROSS INTERFACE, USE FP3                  
                  I=I+2                 
                  DO K=1,5
                     SUM=SUM+WTYPE(1)*UHS(IY,DATAL%LGRD-2+K)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                   +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAX(IY,IX)

                  IX=IX+1        !APPROXIMATE IX+2
                  CALL GETWTYPE(DATAR%ITYPE,2,WTYPE1,2) 
                  DO J=-1,1
                     WTYPE(J)=WTYPE1(1,J)
                  END DO
                  !IX-1 CROSS INTERFACE, USE FP2
                  SUM=0.0D0
                  DO J=0,1
                     SUM=SUM+WTYPE(J)*UHS(IY,IX+J)
                  END DO
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     SUM=SUM+WTYPE(-1)*UHS(IY,DATAL%LGRD-2+K)*DATAL%WIJ2(I,K)
                  END DO
                  SUM=SUM+WTYPE(-1)*(UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)&
                                    +UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  UHS2(IY,IX)=SUM*(DT/2.0D0)*BETAX(IY,IX)               

                  ELEM => ELEM%NEXT%NEXT    ! SKIP THE COUPLE OF CORNER INTERFACES

               END IF

               DEALLOCATE(WTYPE1,STAT=ierr)               

            END DO
         END IF 
   
         UH(IY,:) = UH(IY,:)+UHS2(IY,:) !ADD THE THIRD PART OF RHS

      END DO !----- END OF DO IX=2,NX-1 TO GENERATE RHS 

      !----- STEP IV: SET FIRST ORDER BOUNDARY CONDITION FOR U^{N+1}-----!
      CALL SETBC(T+DT,0.0D0,UH)  

      DO IX=2,NX-1

         !----- STEP V:   GENERATE LHS MATRIX -----!
         !----- STEP V-1: GENERATE LHS MATRIX: (1/BETA-DT/2 D_YY -DT/2 BETA_Y/BETA D_Y) U^(N+1) WITHOUT MIB -----!
         A(1)=0.0D0; B(1)=1.0D0; C(1)=0.0D0 
         DO IY=2,NY-1
            A(IY)=-(VDE2(-1)*DT/2.0D0+VDE1(-1)*DT/2.0D0*BETAY(IY,IX))
            B(IY)=1.0D0/BETA(IY,IX)-(VDE2(0)*DT/2.0D0+VDE1(0)*DT/2.0D0*BETAY(IY,IX))
            C(IY)=-(VDE2(1)*DT/2.0D0+VDE1(1)*DT/2.0D0*BETAY(IY,IX))
         END DO
         A(NY)=0.0D0; B(NY)=1.0D0; C(NY)=0.0D0
         
         DO IY=1,NY  !R IS FOR RHS
            R(IY)=UH(IY,IX)
         END DO
         
         !----- STEP V-2: GENERATE LHS MATRIX WITH MIB -----!
         IF ( ASSOCIATED(IFPX(IX)%HEAD) ) THEN
            ELEM => IFPX(IX)%HEAD
            DO WHILE ( ASSOCIATED(ELEM) )

               !----- AN IRREGULAR INTERFACE POINT                
               IF (ELEM%DATA%ID .GT. 0) THEN
                  DATA=LIST_GET_DATA(ELEM)

                  UJP=DATA%JUMP(1)
                  BUXJP=DATA%JUMP(6)

                  ALLOCATE(WTYPE1(2,-1:1),WTYPE2(2,-1:1),STAT=ierr)
                  CALL GETWTYPE(DATA%ITYPE,1,WTYPE1,2)
                  CALL GETWTYPE(DATA%ITYPE,2,WTYPE2,2)
                  
                  ROW=0.0D0
                  IY=DATA%LGRD  !APPROXIMATE IY
                  !USE FP2
                  ROW(1,1)=-(WTYPE2(2,-1)*DT/2.0D0+WTYPE1(2,-1)*DT/2.0D0*BETAY(IY,IX))
                  ROW(1,2)=1.0D0/BETA(IY,IX)-(WTYPE2(2,0)*DT/2.0D0+WTYPE1(2,0)*DT/2.0D0*BETAY(IY,IX))
                  DO K=1,4
                     ROW(1,K)=ROW(1,K)-(WTYPE2(2,1)*DT/2.0D0+WTYPE1(2,1)*DT/2.0D0*BETAY(IY,IX))*DATA%WIJ(2,K)
                  END DO
                  ROW(1,5)=UH(IY,IX)+(WTYPE2(2,1)*DT/2.0D0+WTYPE1(2,1)*DT/2.0D0*BETAY(IY,IX))*&
                           (UJP*DATA%WIJ(2,5)+BUXJP*DATA%WIJ(2,6))

                  IY=IY+1      !APPROXIMATE IY+1
                  !USE FP1
                  ROW(2,3)=1.0D0/BETA(IY,IX)-(WTYPE2(1,0)*DT/2.0D0+WTYPE1(1,0)*DT/2.0D0*BETAY(IY,IX))
                  ROW(2,4)=-(WTYPE2(1,1)*DT/2.0D0+WTYPE1(1,1)*DT/2.0D0*BETAY(IY,IX))
                  DO K=1,4
                     ROW(2,K)=ROW(2,K)-(WTYPE2(1,-1)*DT/2.0D0+WTYPE1(1,-1)*DT/2.0D0*BETAY(IY,IX))*DATA%WIJ(1,K)
                  END DO
                  ROW(2,5)=UH(IY,IX)+(WTYPE2(1,-1)*DT/2.0D0+WTYPE1(1,-1)*DT/2.0D0*BETAY(IY,IX))*&
                           (UJP*DATA%WIJ(1,5)+BUXJP*DATA%WIJ(1,6))

                  CALL GE1(NY,A,B,C,R,DATA%LGRD,ROW)

                  ELEM => ELEM%NEXT         ! NEXT INTERFACE

               !----- A CORNER POINT
               ELSE
                  DATAL=LIST_GET_DATA(ELEM)
                  DATAR=LIST_GET_DATA(ELEM%NEXT)

                  UJPL=DATAL%JUMP(1)
                  UJPR=DATAR%JUMP(1)
                  BUXJPL=DATAL%JUMP(6)
                  BUXJPR=DATAR%JUMP(6)

                  ALLOCATE(WTYPE1(2,-1:1),WTYPE2(2,-1:1),STAT=ierr)

                  ROW2=0.0D0

                  IY=DATAL%LGRD   !APPROXIMATE IY
                  CALL GETWTYPE(DATAL%ITYPE,1,WTYPE1,2)
                  CALL GETWTYPE(DATAL%ITYPE,2,WTYPE2,2)
                  DO J=-1,1
                     WTYPEC1(J)=WTYPE1(2,J)
                     WTYPEC2(J)=WTYPE2(2,J)
                  END DO
                  !IY+1 CROSS INTERFACE, USE FP2 
                  ROW2(1,1)=-(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAY(IY,IX))
                  ROW2(1,2)=1.0D0/BETA(IY,IX)-(WTYPEC2(0)*DT/2.0D0+WTYPEC1(0)*DT/2.0D0*BETAY(IY,IX))
                  I=2                              !USE FP2, INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     ROW2(1,K)=ROW2(1,K)-(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAY(IY,IX))*DATAL%WIJ2(I,K)
                  END DO
                  ROW2(1,6)=UH(IY,IX)+(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAY(IY,IX))*&
                                      (UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)+&
                                       UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))

                  IY=IY+1        !APPROXIMATE IY+1
                  !DETERMINE USE WHICH INTERFACE DATA APPROXIMATE THE MIDDLE POINT, IX+1
                  IF (ABS(DATAL%ITYPE) .LT. 2) THEN            !USE LEFT INTERFACE WEIGHTS
                     CALL GETWTYPE(DATAL%ITYPE,1,WTYPE1,2)  
                     CALL GETWTYPE(DATAL%ITYPE,2,WTYPE2,2)
                     DO J=-1,1
                        WTYPEC1(J)=WTYPE1(1,J)
                        WTYPEC2(J)=WTYPE2(1,J)
                     END DO
                  ELSE                                         !USE RIGHT INTERFACE WEIGHTS  
                     CALL GETWTYPE(DATAR%ITYPE,1,WTYPE1,2)  
                     CALL GETWTYPE(DATAR%ITYPE,2,WTYPE2,2)
                     DO J=-1,1
                        WTYPEC1(J)=WTYPE1(2,J)
                        WTYPEC2(J)=WTYPE2(2,J)
                     END DO
                  END IF    
                  ROW2(2,3)=1.0D0/BETA(IY,IX)-(WTYPEC2(0)*DT/2.0D0+WTYPEC1(0)*DT/2.0D0*BETAY(IY,IX))
                  I=1                              !INDEX USED FOR FP WIJ2(:,:)  
                  !IY-1 CROSS INTERFACE, USE FP1 
                  DO K=1,5
                     ROW2(2,K)=ROW2(2,K)-(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAY(IY,IX))*DATAL%WIJ2(I,K)
                  END DO
                  ROW2(2,6)=UH(IY,IX)+(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAY(IY,IX))*&
                                       (UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)+&
                                        UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))
                  !IY+1 CROSS INTERFACE, USE FP3 
                  I=I+2
                  DO K=1,5
                     ROW2(2,K)=ROW2(2,K)-(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAY(IY,IX))*DATAL%WIJ2(I,K)
                  END DO
                  ROW2(2,6)=ROW2(2,6)+(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAY(IY,IX))*&
                                      (UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)+&
                                       UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))

                  IY=IY+1        !APPROXIMATE IY+2
                  CALL GETWTYPE(DATAR%ITYPE,1,WTYPE1,2)
                  CALL GETWTYPE(DATAR%ITYPE,2,WTYPE2,2)
                  DO J=-1,1
                     WTYPEC1(J)=WTYPE1(1,J)
                     WTYPEC2(J)=WTYPE2(1,J)
                  END DO
                  !IY-1 CROSS INTERFACE, USE FP2 
                  ROW2(3,4)=1.0D0/BETA(IY,IX)-(WTYPEC2(0)*DT/2.0D0+WTYPEC1(0)*DT/2.0D0*BETAY(IY,IX))
                  ROW2(3,5)=-(WTYPEC2(1)*DT/2.0D0+WTYPEC1(1)*DT/2.0D0*BETAY(IY,IX))                  
                  I=2                              !INDEX USED FOR FP WIJ2(:,:)
                  DO K=1,5
                     ROW2(3,K)=ROW2(3,K)-(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAY(IY,IX))*DATAL%WIJ2(I,K)
                  END DO
                  ROW2(3,6)=UH(IY,IX)+(WTYPEC2(-1)*DT/2.0D0+WTYPEC1(-1)*DT/2.0D0*BETAY(IY,IX))*&
                                       (UJPL*DATAL%WIJ2(I,6)+BUXJPL*DATAL%WIJ2(I,7)+&
                                        UJPR*DATAR%WIJ2(I,8)+BUXJPR*DATAR%WIJ2(I,9))

                  CALL GE2(NY,A,B,C,R,DATAL%LGRD,ROW2)

                  ELEM => ELEM%NEXT%NEXT    !SKIP THE COUPLE OF CORNER INTERFACES

               END IF 
               
               DEALLOCATE(WTYPE1,WTYPE2,STAT=ierr)

            END DO
         END IF
         
         !----- STEP VI: THOMAS ALGORITHM -----!
         CALL TRIDAG(A,B,C,R,UT,NY)

         DO IY=1,NY
            UH(IY,IX)=UT(IY)
         END DO      

      END DO         !END OF DO IY=2,NY-1 TO GENERATE RHS MATRIX

      DEALLOCATE(A,B,C,R,UT,STAT=ierr)

      RETURN

      END SUBROUTINE ADIPR
