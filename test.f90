program trapezetest
implicit none
DOUBLE PRECISION :: a, b, int,dum,func
INTEGER :: i
dum=0.0
a=0.0
b=2.0
int=0.0
do i=1,50
CALL trapzd(a,b,int,i,dum)

write(*,*) i,int
end do
end program trapezetest

 function func(x,dummy)
	implicit none
	DOUBLE PRECISION :: func,dummy,x
	func=2.0*x
end function

SUBROUTINE trapzd(a,b,s,n,func_arg)
	implicit none

        INTEGER n
        DOUBLE PRECISION a,b,s,func_arg,func
        INTEGER it,j
        DOUBLE PRECISION del,sum,tnm,x
        IF(n.eq.1) THEN
            s=0.5*(b-a)*(func(a,func_arg)+func(b,func_arg))
        ELSE
            it=2**(n-2)
            tnm=it
            del=(b-a)/tnm
            x=a+0.5*del
            sum=0.
            DO j=1,it
                sum=sum+func(x,func_arg)
                x=x+del
            END DO
            s=0.5*(s+(b-a)*sum/tnm)
        ENDIF
        RETURN
    END