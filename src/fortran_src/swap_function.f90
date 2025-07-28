!=======================================================================
!
!  Functions to swap the values of two scalar numbers or two arrays
!  (either 1- or 2-D) of numbers. The function is overloaded, so it
!  can be called regardless of the type (integer or double) or size
!  (scalar or array) of the input variables.
!
!-----------------------------------------------------------------------
MODULE SWAP_FUNCTION

   INTERFACE SWAP
      MODULE PROCEDURE SWAP_INTEGER,SWAP_REAL,SWAP_INTEGER_1D_ARRAY,SWAP_REAL_1D_ARRAY,SWAP_INTEGER_2D_ARRAY,SWAP_REAL_2D_ARRAY
   END INTERFACE SWAP

CONTAINS

   SUBROUTINE SWAP_INTEGER(A,B)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: A,B
      INTEGER :: D
      D = A
      A = B
      B = D
      RETURN
   END SUBROUTINE SWAP_INTEGER

   SUBROUTINE SWAP_REAL(A,B)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: A,B
      DOUBLE PRECISION :: D
      D = A
      A = B
      B = D
      RETURN
   END SUBROUTINE SWAP_REAL

   SUBROUTINE SWAP_INTEGER_1D_ARRAY(A,B)
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: A,B
      INTEGER, DIMENSION(:), ALLOCATABLE   :: D
      IF(SIZE(A).NE.SIZE(B)) THEN
         WRITE(6,*) 'ERROR! Cannot swap values between arrays of different dimensions'
         WRITE(6,*) 'SHAPE(A) =',SHAPE(A),'; SHAPE(B) =',SHAPE(B)
         STOP
      END IF
      ALLOCATE(D(SIZE(A)))
      D = A
      A = B
      B = D
      DEALLOCATE(D)
      RETURN
   END SUBROUTINE SWAP_INTEGER_1D_ARRAY

   SUBROUTINE SWAP_REAL_1D_ARRAY(A,B)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: A,B
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: D
      IF(SIZE(A).NE.SIZE(B)) THEN
         WRITE(6,*) 'ERROR! Cannot swap values between arrays of different dimensions'
         WRITE(6,*) 'SHAPE(A) =',SHAPE(A),'; SHAPE(B) =',SHAPE(B)
         STOP
      END IF
      ALLOCATE(D(SIZE(A)))
      D = A
      A = B
      B = D
      DEALLOCATE(D)
      RETURN
   END SUBROUTINE SWAP_REAL_1D_ARRAY


   SUBROUTINE SWAP_INTEGER_2D_ARRAY(A,B)
      IMPLICIT NONE
      INTEGER, DIMENSION(:,:), INTENT(INOUT) :: A,B
      INTEGER, DIMENSION(:,:), ALLOCATABLE   :: D
      IF(SIZE(A,1).NE.SIZE(B,1) .OR. SIZE(A,2).NE.SIZE(B,2)) THEN
         WRITE(6,*) 'ERROR! Cannot swap values between arrays of different dimensions'
         WRITE(6,*) 'SHAPE(A) =',SHAPE(A),'; SHAPE(B) =',SHAPE(B)
         STOP
      END IF
      ALLOCATE(D(SIZE(A,1),SIZE(A,2)))
      D = A
      A = B
      B = D
      DEALLOCATE(D)
      RETURN
   END SUBROUTINE SWAP_INTEGER_2D_ARRAY

   SUBROUTINE SWAP_REAL_2D_ARRAY(A,B)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: A,B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: D
      IF(SIZE(A,1).NE.SIZE(B,1) .OR. SIZE(A,2).NE.SIZE(B,2)) THEN
         WRITE(6,*) 'ERROR! Cannot swap values between arrays of different dimensions'
         WRITE(6,*) 'SHAPE(A) =',SHAPE(A),'; SHAPE(B) =',SHAPE(B)
         STOP
      END IF
      ALLOCATE(D(SIZE(A,1),SIZE(A,2)))
      D = A
      A = B
      B = D
      DEALLOCATE(D)
      RETURN
   END SUBROUTINE SWAP_REAL_2D_ARRAY

END MODULE SWAP_FUNCTION
!=======================================================================
