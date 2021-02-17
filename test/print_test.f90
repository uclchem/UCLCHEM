PROGRAM PRINT_TEST
	IMPLICIT NONE
	DOUBLE PRECISION :: awkward=1.0d-102,ok=1.0d-10
	write(*,11) awkward
	write(*,11) ok
	11 format(1pe15.5e3)
END PROGRAM PRINT_TEST