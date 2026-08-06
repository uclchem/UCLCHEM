module numerics
    use constants, only: dp
    use iso_fortran_env, only: output_unit

    implicit none

    public

   interface SWAP
      module procedure SWAP_INTEGER,SWAP_REAL,SWAP_INTEGER_1D_ARRAY,SWAP_REAL_1D_ARRAY,SWAP_INTEGER_2D_ARRAY,SWAP_REAL_2D_ARRAY
   end interface SWAP

contains

    ! Evaluate a polynomial using Horner's method.
    ! Taken (and corrected) from https://rosettacode.org/wiki/Horner%27s_rule_for_polynomial_evaluation#Fortran
    pure function evaluate_polynomial(coeffs, x) result(value)
      real(dp), dimension(:), intent(in) :: coeffs  ! Polynomial coefficients (starting at x^0)
      real(dp), intent(in) :: x  ! x-value where to evaluate the polynomial
      real(dp) :: value
      integer :: i

      value = coeffs(size(coeffs))
      do i = size(coeffs)-1, 1, -1
        value = value * x + coeffs(i)
      end do
    end function evaluate_polynomial

    ! Taken from https://fortranwiki.org/fortran/show/integration
    pure function integrate_trapezoid(x, y) result(integral)
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(in), dimension(size(x)) :: y
        real(dp) :: integral

        associate(n => size(x))
            integral = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
        end associate
    end function integrate_trapezoid

    pure function logspace(start, stop, num) result(points)
        real(dp), intent(in) :: start, stop
        integer, intent(in) :: num

        real(dp), dimension(num) :: points

        integer :: i
        do i = 1, num
            points(i) = 10.0_dp**(start + (i-1)*(stop-start)/DBLE(num-1))
        end do
    end function logspace

!=======================================================================
!
!  Functions to swap the values of two scalar numbers or two arrays
!  (either 1- or 2-D) of numbers. The function is overloaded, so it
!  can be called regardless of the type (integer or double) or size
!  (scalar or array) of the input variables.
!
!-----------------------------------------------------------------------
   pure subroutine SWAP_INTEGER(A,B)
      integer, intent(inout) :: A,B
      integer :: D
      D = A
      A = B
      B = D
   end subroutine SWAP_INTEGER

   pure subroutine SWAP_REAL(A,B)
      real(dp), intent(inout) :: A,B
      real(dp) :: D
      D = A
      A = B
      B = D
   end subroutine SWAP_REAL

   subroutine SWAP_INTEGER_1D_ARRAY(A,B)
      integer, dimension(:), intent(inout) :: A,B
      integer, dimension(:), allocatable   :: D
      if(SIZE(A)/=SIZE(B)) then
         write(output_unit,*) "ERROR! Cannot swap values between arrays of different dimensions"
         write(output_unit,*) "SHAPE(A) =",SHAPE(A),"; SHAPE(B) =",SHAPE(B)
         stop
      end if
      allocate(D(SIZE(A)))
      D = A
      A = B
      B = D
      deallocate(D)
   end subroutine SWAP_INTEGER_1D_ARRAY

   subroutine SWAP_REAL_1D_ARRAY(A,B)
      real(dp), dimension(:), intent(inout) :: A,B
      real(dp), dimension(:), allocatable   :: D
      if(SIZE(A)/=SIZE(B)) then
         write(output_unit,*) "ERROR! Cannot swap values between arrays of different dimensions"
         write(output_unit,*) "SHAPE(A) =",SHAPE(A),"; SHAPE(B) =",SHAPE(B)
         stop
      end if
      allocate(D(SIZE(A)))
      D = A
      A = B
      B = D
      deallocate(D)
   end subroutine SWAP_REAL_1D_ARRAY


   subroutine SWAP_INTEGER_2D_ARRAY(A,B)
      integer, dimension(:,:), intent(inout) :: A,B
      integer, dimension(:,:), allocatable   :: D
      if(SIZE(A,1)/=SIZE(B,1) .OR. SIZE(A,2)/=SIZE(B,2)) then
         write(output_unit,*) "ERROR! Cannot swap values between arrays of different dimensions"
         write(output_unit,*) "SHAPE(A) =",SHAPE(A),"; SHAPE(B) =",SHAPE(B)
         stop
      end if
      allocate(D(SIZE(A,1),SIZE(A,2)))
      D = A
      A = B
      B = D
      deallocate(D)
   end subroutine SWAP_INTEGER_2D_ARRAY

   subroutine SWAP_REAL_2D_ARRAY(A,B)
      real(dp), dimension(:,:), intent(inout) :: A,B
      real(dp), dimension(:,:), allocatable   :: D
      if(SIZE(A,1)/=SIZE(B,1) .OR. SIZE(A,2)/=SIZE(B,2)) then
         write(output_unit,*) "ERROR! Cannot swap values between arrays of different dimensions"
         write(output_unit,*) "SHAPE(A) =",SHAPE(A),"; SHAPE(B) =",SHAPE(B)
         stop
      end if
      allocate(D(SIZE(A,1),SIZE(A,2)))
      D = A
      A = B
      B = D
      deallocate(D)
   end subroutine SWAP_REAL_2D_ARRAY

!=======================================================================

end module numerics
