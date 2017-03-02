      ! LINKEDLIST.F90 --
      !     INCLUDE FILE FOR DEFINING LINKED LISTS WHERE EACH ELEMENT HOLDS
      !     THE SAME KIND OF DATA
      !
      !     SEE THE EXAMPLE/TEST PROGRAM FOR THE WAY TO USE THIS
      !
      !     NOTE:
      !     YOU SHOULD ONLY USE POINTER VARIABLES OF THIS TYPE, NO
      !     ORDINARY VARIABLES, AS SOMETIMES THE MEMORY POINTED TO
      !     WILL BE DEALLOCATED. THE SUBROUTINES AND FUNCTIONS
      !     ARE DESIGNED TO MINIMIZE MISTAKES (FOR INSTANCE: USING
      !     = INSTEAD OF =>)
      !
      !     $ID: LINKEDLIST.F90,V 1.3 2007/01/26 09:56:43 ARJENMARKUS EXP $
      !
      ! DEFINE THE LINKED-LIST DATA TYPE
      TYPE LINKED_LIST

      TYPE(LINKED_LIST), POINTER :: NEXT
      TYPE(LIST_DATA)            :: DATA

      END TYPE LINKED_LIST

      ! DEFINE A PRIVATE (!) INTERFACE TO PREVENT
      ! MISTAKES WITH ORDINARY ASSIGNMENT
      !
      !INTERFACE ASSIGNMENT(=)
      !    MODULE PROCEDURE LIST_ASSIGN
      !END INTERFACE
      !PRIVATE :: LIST_ASSIGN

      ! DEFINE THE SUBROUTINES AND FUNCTIONS
      !
      CONTAINS

      ! LIST_ASSIGN
      !     SUBROUTINE TO PREVENT ERRORS WITH ASSIGNMENT
      ! ARGUMENTS:
      !     LIST_LEFT   LIST ON THE LEFT-HAND SIDE
      !     LIST_RIGHT  LIST ON THE RIGHT-HAND SIDE
      !
      ! NOTE:
      !     THIS DOES NOT WORK BECAUSE OF A PRIVATE/PUBLIC
      !     CONFLICT
      !
      !SUBROUTINE LIST_ASSIGN( LIST_LEFT, LIST_RIGHT )
      !    TYPE(LINKED_LIST), INTENT(OUT)  :: LIST_LEFT
      !    TYPE(LINKED_LIST), INTENT(IN)   :: LIST_RIGHT
      !   !TYPE(LINKED_LIST), POINTER      :: LIST_LEFT
      !   !TYPE(LINKED_LIST), POINTER      :: LIST_RIGHT
      !
      !    !
      !    ! NOTE THE ORDER!
      !    !
      !    STOP 'ERROR: ORDINARY ASSIGNMENT FOR LISTS'
      !    LIST_LEFT%NEXT => NULL()
      !END SUBROUTINE LIST_ASSIGN

      !-------------------------------------------------------------------------
      ! LIST_CREATE --
      !     CREATE AND INITIALISE A LIST
      ! ARGUMENTS:
      !     LIST       POINTER TO NEW LINKED LIST
      !     DATA       THE DATA FOR THE FIRST ELEMENT
      ! NOTE:
      !     THIS VERSION ASSUMES A SHALLOW COPY IS ENOUGH (THAT IS, THERE ARE
      !     NO POINTERS WITHIN THE DATA TO BE STORED)
      !     IT ALSO ASSUMES THE ARGUMENT LIST DOES NOT ALREADY REFER TO A LIST.
      !     USE LIST_DESTROY FIRST TO DESTROY UP AN OLD LIST.
      SUBROUTINE LIST_CREATE( LIST, DATA )

      TYPE(LINKED_LIST), POINTER  :: LIST
      TYPE(LIST_DATA), INTENT(IN) :: DATA

      ALLOCATE( LIST )
      LIST%NEXT => NULL()
      LIST%DATA =  DATA

      END SUBROUTINE LIST_CREATE

      !-------------------------------------------------------------------------
      ! LIST_DESTROY --
      !     DESTROY AN ENTIRE LIST
      ! ARGUMENTS:
      !     LIST       POINTER TO THE LIST TO BE DESTROYED
      ! NOTE:
      !     THIS VERSION ASSUMES THAT THERE ARE NO
      !     POINTERS WITHIN THE DATA THAT NEED DEALLOCATION
      !
      SUBROUTINE LIST_DESTROY( LIST )

      TYPE(LINKED_LIST), POINTER  :: LIST

      TYPE(LINKED_LIST), POINTER  :: CURRENT
      TYPE(LINKED_LIST), POINTER  :: NEXT

      CURRENT => LIST
      DO WHILE ( ASSOCIATED(CURRENT%NEXT) )
         NEXT => CURRENT%NEXT
         DEALLOCATE( CURRENT )
         CURRENT => NEXT
      END DO

      END SUBROUTINE LIST_DESTROY

      !-------------------------------------------------------------------------
      ! LIST_COUNT --
      !     COUNT THE NUMBER OF ITEMS IN THE LIST
      ! ARGUMENTS:
      !     LIST       POINTER TO THE LIST
      !
      INTEGER FUNCTION LIST_COUNT( LIST )

      TYPE(LINKED_LIST), POINTER  :: LIST

      TYPE(LINKED_LIST), POINTER  :: CURRENT
      !TYPE(LINKED_LIST), POINTER  :: NEXT

      IF ( ASSOCIATED(LIST) ) THEN
         LIST_COUNT = 1
         CURRENT => LIST
         DO WHILE ( ASSOCIATED(CURRENT%NEXT) )
            CURRENT => CURRENT%NEXT
            LIST_COUNT = LIST_COUNT + 1
         END DO
      ELSE
         LIST_COUNT = 0
      END IF

      END FUNCTION LIST_COUNT

      !-------------------------------------------------------------------------
      ! LIST_NEXT
      !     RETURN THE NEXT ELEMENT (IF ANY)
      ! ARGUMENTS:
      !     ELEM       ELEMENT IN THE LINKED LIST
      ! RESULT:
      !
      FUNCTION LIST_NEXT( ELEM ) RESULT(NEXT)

      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LINKED_LIST), POINTER :: NEXT

      NEXT => ELEM%NEXT

      END FUNCTION LIST_NEXT

      !-------------------------------------------------------------------------
      ! LIST_INSERT
      !     INSERT A NEW ELEMENT
      ! ARGUMENTS:
      !     ELEM       ELEMENT IN THE LINKED LIST AFTER WHICH TO INSERT THE NEW
      !                ELEMENT
      !     DATA       THE DATA FOR THE NEW ELEMENT
      SUBROUTINE LIST_INSERT( ELEM, DATA )

      TYPE(LINKED_LIST), POINTER  :: ELEM
      TYPE(LIST_DATA), INTENT(IN) :: DATA
      TYPE(LINKED_LIST), POINTER :: NEXT

      ALLOCATE(NEXT)

      NEXT%NEXT => ELEM%NEXT
      ELEM%NEXT => NEXT
      NEXT%DATA =  DATA

      END SUBROUTINE LIST_INSERT

      !-------------------------------------------------------------------------
      ! LIST_APPEND
      !     APPEND A NEW ELEMENT AT THE END OF THE LIST
      ! ARGUMENTS:
      !     LIST       START OF THE LIST
      !     DATA       THE DATA TO BE APPENDED AT THE END OF THE LINKED LIST
      SUBROUTINE LIST_APPEND( LIST, DATA )

      TYPE(LINKED_LIST), POINTER  :: LIST
      TYPE(LIST_DATA), INTENT(IN) :: DATA

      TYPE(LINKED_LIST), POINTER  :: ELEM
      TYPE(LINKED_LIST), POINTER  :: NEXT

      IF (.NOT. ASSOCIATED(LIST)) THEN
         CALL LIST_CREATE(LIST,DATA)
      ELSE
         ELEM => LIST

         DO WHILE ( ASSOCIATED(ELEM%NEXT) )
            ELEM => LIST_NEXT(ELEM)
         END DO

         ALLOCATE(NEXT)

         ELEM%NEXT => NEXT
         NEXT%NEXT => NULL()
         NEXT%DATA =  DATA
      END IF

      RETURN

      END SUBROUTINE LIST_APPEND

      !-------------------------------------------------------------------------
      ! LIST_INSERT_HEAD
      !     INSERT A NEW ELEMENT BEFORE THE FIRST ELEMENT
      ! ARGUMENTS:
      !     LIST       START OF THE LIST
      !     DATA       THE DATA FOR THE NEW ELEMENT
      !
      SUBROUTINE LIST_INSERT_HEAD( LIST, DATA )

      TYPE(LINKED_LIST), POINTER  :: LIST
      TYPE(LIST_DATA), INTENT(IN) :: DATA
      TYPE(LINKED_LIST), POINTER :: ELEM

      ALLOCATE(ELEM)
      ELEM%DATA =  DATA

      ELEM%NEXT => LIST
      LIST      => ELEM

      END SUBROUTINE LIST_INSERT_HEAD

      !-------------------------------------------------------------------------
      ! LIST_DELETE_ELEMENT
      !     DELETE AN ELEMENT FROM THE LIST
      ! ARGUMENTS:
      !     LIST       HEADER OF THE LIST
      !     ELEM       ELEMENT IN THE LINKED LIST TO BE
      !                REMOVED
      !
      SUBROUTINE LIST_DELETE_ELEMENT( LIST, ELEM )

      TYPE(LINKED_LIST), POINTER  :: LIST
      TYPE(LINKED_LIST), POINTER  :: ELEM
      TYPE(LINKED_LIST), POINTER  :: CURRENT
      TYPE(LINKED_LIST), POINTER  :: PREV

      IF ( ASSOCIATED(LIST,ELEM) ) THEN
         LIST => ELEM%NEXT
         DEALLOCATE( ELEM )
      ELSE
         CURRENT => LIST
         PREV    => LIST
         DO WHILE ( ASSOCIATED(CURRENT) )
            IF ( ASSOCIATED(CURRENT,ELEM) ) THEN
               PREV%NEXT => CURRENT%NEXT
               DEALLOCATE( CURRENT ) ! IS ALSO "ELEM"
               EXIT
            END IF
            PREV    => CURRENT
            CURRENT => CURRENT%NEXT
         END DO
      END IF

      !ALLOCATE(NEXT)
      !
      !NEXT%NEXT => ELEM%NEXT
      !ELEM%NEXT => NEXT
      !NEXT%DATA =  DATA

      END SUBROUTINE LIST_DELETE_ELEMENT

      !-------------------------------------------------------------------------
      ! LIST_GET_DATA
      !     GET THE DATA STORED WITH A LIST ELEMENT
      ! ARGUMENTS:
      !     ELEM       ELEMENT IN THE LINKED LIST
      !
      FUNCTION LIST_GET_DATA( ELEM ) RESULT(DATA)

      TYPE(LINKED_LIST), POINTER :: ELEM
      TYPE(LIST_DATA)            :: DATA

      DATA = ELEM%DATA

      END FUNCTION LIST_GET_DATA

      !-------------------------------------------------------------------------
      ! LIST_PUT_DATA
      !     STORE NEW DATA WITH A LIST ELEMENT
      ! ARGUMENTS:
      !     ELEM       ELEMENT IN THE LINKED LIST
      !     DATA       THE DATA TO BE STORED
      !
      SUBROUTINE LIST_PUT_DATA( ELEM, DATA )

      TYPE(LINKED_LIST), POINTER  :: ELEM
      TYPE(LIST_DATA), INTENT(IN) :: DATA

      ELEM%DATA = DATA

      END SUBROUTINE LIST_PUT_DATA

