module numerics
    use constants, only: dp
    use iso_fortran_env, only: output_unit

    implicit none

    public

   interface SWAP
      module procedure SWAP_INTEGER,SWAP_REAL,SWAP_INTEGER_1D_ARRAY,SWAP_REAL_1D_ARRAY,SWAP_INTEGER_2D_ARRAY,SWAP_REAL_2D_ARRAY
   end interface SWAP

contains

    ! Test that two values are equal, up to machine precision.
    elemental function is_equal(a, b, atol) result(equal)
        real(dp), intent(in) :: a, b
        real(dp), intent(in), optional :: atol

        logical :: equal

        real(dp) :: tol

        if (.not. present(atol)) then
            tol = epsilon(1.0_dp)
        else
            tol = atol
        end if
    
        if (abs(a-b) <= tol) then
            equal = .true.
        else
            equal = .false.
        end if
    end function is_equal
    
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


    pure subroutine pair_insertion_sort(array)
        real(dp), intent(inout) :: array(:)
        integer :: i,j,last
        real(dp) :: t1,t2

        last=size(array)
        do i=2,last-1,2
           t1=min(array(i),array(i+1))
           t2=max(array(i),array(i+1))
           j=i-1
           do while((j>=1).and.(array(j)>t2))
              array(j+2)=array(j)
              j=j-1
           end do
           array(j+2)=t2
           do while((j>=1).and.(array(j)>t1))
              array(j+1)=array(j)
              j=j-1
           end do
           array(j+1)=t1
        end do

        if(mod(last,2)==0)then
           t1=array(last)
           loop: do j=last-1,1,-1
              if (array(j)<=t1) exit loop
              array(j+1)=array(j)
           end do loop
           array(j+1)=t1
        end if
    end subroutine pair_insertion_sort

    !-----------------------------------------------------------------------
    ! pair_insertion_sort_with_perm - Sort a 1D array and track permutation
    !
    ! Sorts a 1D array using pair insertion sort while maintaining a
    ! permutation array that tracks where each element originally came from.
    ! This allows you to reorder other arrays using the same permutation.
    !
    ! Arguments:
    !   array(:)    - 1D array to sort (modified in place)
    !-----------------------------------------------------------------------
    pure subroutine pair_insertion_sort_with_perm(array, perm)
        real(dp), intent(inout) :: array(:)
        integer, intent(out), dimension(size(array)) :: perm
        integer :: i, j, last
        real(dp) :: t1, t2
        integer :: p1, p2

        last = SIZE(array)

        ! Initialize permutation array (caller must allocate with correct size)
        perm = [(i, i=1, last)]

        ! Pair insertion sort - process elements two at a time
        do i = 2, last-1, 2
           ! Get the pair and their permutation indices
           if (array(i) <= array(i+1)) then
              t1 = array(i)
              t2 = array(i+1)
              p1 = perm(i)
              p2 = perm(i+1)
           else
              t1 = array(i+1)
              t2 = array(i)
              p1 = perm(i+1)
              p2 = perm(i)
           end if

           ! Find position for larger element (t2)
           j = i - 1
           do while ((j >= 1) .AND. (array(j) > t2))
              array(j+2) = array(j)
              perm(j+2) = perm(j)
              j = j - 1
           end do
           array(j+2) = t2
           perm(j+2) = p2

           ! Find position for smaller element (t1)
           do while ((j >= 1) .AND. (array(j) > t1))
              array(j+1) = array(j)
              perm(j+1) = perm(j)
              j = j - 1
           end do
           array(j+1) = t1
           perm(j+1) = p1
        end do

        ! Handle last element if array has even number of elements
        if (MOD(last, 2) == 0) then
           t1 = array(last)
           p1 = perm(last)
           loop: do j = last-1, 1, -1
              if (array(j) <= t1) exit loop
              array(j+1) = array(j)
              perm(j+1) = perm(j)
           end do loop
           array(j+1) = t1
           perm(j+1) = p1
        end if
    end subroutine pair_insertion_sort_with_perm

end module numerics
