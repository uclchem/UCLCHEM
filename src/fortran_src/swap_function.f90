!=======================================================================
!
!  Functions to swap the values of two scalar numbers or two arrays
!  (either 1- or 2-D) of numbers. The function is overloaded, so it
!  can be called regardless of the type (integer or double) or size
!  (scalar or array) of the input variables.
!
!-----------------------------------------------------------------------
module SWAP_FUNCTION

   interface SWAP
      module procedure SWAP_INTEGER,SWAP_REAL,SWAP_INTEGER_1D_ARRAY,SWAP_REAL_1D_ARRAY,SWAP_INTEGER_2D_ARRAY,SWAP_REAL_2D_ARRAY
   end interface SWAP

contains

   subroutine SWAP_INTEGER(A,B)
      implicit none
      integer, intent(inout) :: A,B
      integer :: D
      D = A
      A = B
      B = D
   end subroutine SWAP_INTEGER

   subroutine SWAP_REAL(A,B)
      implicit none
      double precision, intent(inout) :: A,B
      double precision :: D
      D = A
      A = B
      B = D
   end subroutine SWAP_REAL

   subroutine SWAP_INTEGER_1D_ARRAY(A,B)
      implicit none
      integer, dimension(:), intent(inout) :: A,B
      integer, dimension(:), allocatable   :: D
      if(SIZE(A)/=SIZE(B)) then
         write(6,*) "ERROR! Cannot swap values between arrays of different dimensions"
         write(6,*) "SHAPE(A) =",SHAPE(A),"; SHAPE(B) =",SHAPE(B)
         stop
      end if
      allocate(D(SIZE(A)))
      D = A
      A = B
      B = D
      deallocate(D)
   end subroutine SWAP_INTEGER_1D_ARRAY

   subroutine SWAP_REAL_1D_ARRAY(A,B)
      implicit none
      double precision, dimension(:), intent(inout) :: A,B
      double precision, dimension(:), allocatable   :: D
      if(SIZE(A)/=SIZE(B)) then
         write(6,*) "ERROR! Cannot swap values between arrays of different dimensions"
         write(6,*) "SHAPE(A) =",SHAPE(A),"; SHAPE(B) =",SHAPE(B)
         stop
      end if
      allocate(D(SIZE(A)))
      D = A
      A = B
      B = D
      deallocate(D)
   end subroutine SWAP_REAL_1D_ARRAY


   subroutine SWAP_INTEGER_2D_ARRAY(A,B)
      implicit none
      integer, dimension(:,:), intent(inout) :: A,B
      integer, dimension(:,:), allocatable   :: D
      if(SIZE(A,1)/=SIZE(B,1) .OR. SIZE(A,2)/=SIZE(B,2)) then
         write(6,*) "ERROR! Cannot swap values between arrays of different dimensions"
         write(6,*) "SHAPE(A) =",SHAPE(A),"; SHAPE(B) =",SHAPE(B)
         stop
      end if
      allocate(D(SIZE(A,1),SIZE(A,2)))
      D = A
      A = B
      B = D
      deallocate(D)
   end subroutine SWAP_INTEGER_2D_ARRAY

   subroutine SWAP_REAL_2D_ARRAY(A,B)
      implicit none
      double precision, dimension(:,:), intent(inout) :: A,B
      double precision, dimension(:,:), allocatable   :: D
      if(SIZE(A,1)/=SIZE(B,1) .OR. SIZE(A,2)/=SIZE(B,2)) then
         write(6,*) "ERROR! Cannot swap values between arrays of different dimensions"
         write(6,*) "SHAPE(A) =",SHAPE(A),"; SHAPE(B) =",SHAPE(B)
         stop
      end if
      allocate(D(SIZE(A,1),SIZE(A,2)))
      D = A
      A = B
      B = D
      deallocate(D)
   end subroutine SWAP_REAL_2D_ARRAY

end module SWAP_FUNCTION
!=======================================================================
