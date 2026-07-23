module photoreactions
    use constants
    use DEFAULTPARAMETERS
    !f2py INTEGER, parameter :: dp
implicit none
    real(dp) :: UV_FAC=3.02

    !Below are arrays for self-shielding of CO and H2
    logical :: start=.true.
    integer :: NUM_LAMBDA=30
    real(dp), dimension(30) :: LAMBDA_GRID=(/ &
      &  910.0D0, 950.0D0,1000.0D0,1050.0D0,1110.0D0, &
      & 1180.0D0,1250.0D0,1390.0D0,1490.0D0,1600.0D0, &
      & 1700.0D0,1800.0D0,1900.0D0,2000.0D0,2100.0D0, &
      & 2190.0D0,2300.0D0,2400.0D0,2500.0D0,2740.0D0, &
      & 3440.0D0,4000.0D0,4400.0D0,5500.0D0,7000.0D0, &
      & 9000.0D0,12500.0D0,22000.0D0,34000.0D0,1.0D9/)
    real(dp), dimension(30) :: XLAMBDA_GRID=(/ &
      & 5.76D0,5.18D0,4.65D0,4.16D0,3.73D0, &
      & 3.40D0,3.11D0,2.74D0,2.63D0,2.62D0, &
      & 2.54D0,2.50D0,2.58D0,2.78D0,3.01D0, &
      & 3.12D0,2.86D0,2.58D0,2.35D0,2.00D0, &
      & 1.58D0,1.42D0,1.32D0,1.00D0,0.75D0, &
      & 0.48D0,0.28D0,0.12D0,0.05D0,0.00D0/)
    real(dp), dimension(30) :: XLAMBDA_DERIV

    logical :: startr=.true.

    !  12CO line shielding data from van Dishoeck & Black (1988, ApJ, 334, 771, Table 5)
    integer, parameter ::  DIMCO=7, DIMH2=6
    real(KIND=DP), dimension(8) :: NCO_GRID=(/12.0D0,13.0D0,14.0D0,15.0D0,16.0D0,17.0D0,18.0D0,19.0D0/)
    real(KIND=DP), dimension(6) :: NH2_GRID=(/18.0D0,19.0D0,20.0D0,21.0D0,22.0D0,23.0D0/)
    real(KIND=DP), dimension(8,6) :: SCO_GRID=RESHAPE((/ &
      &  0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00, &
      & -8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00, &
      & -1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00, &
      & -4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00, &
      & -1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00, &
      & -3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/), (/8,6/))
    real(KIND=DP), dimension(8,6) :: SCO_DERIV

    real(dp) :: ICE_GAS_PHOTO_CROSSSECTION_RATIO = 0.3  ! Kalvans 2018
contains

!photodissociation rate of H2  accounting for self-shielding
real(dp) function H2PhotoDissRate(NH2,radField,av,turbVel)
    !H2 Column density to surface, UV at that surface, visual extinction to surface, and turbulent velocity of cloud
    real(dp), intent(in) :: NH2 ,radField,av ,turbVel
    !unshielded reaction rate, characteristic wavelendth of radiation, radiative linewidth
    real(dp), parameter :: baseRate=5.18d-11,xl=1000.0,radWidth=8.0d7
    real(dp) :: dopplerWidth

    dopplerWidth=turbVel/(xl*1.0d-8)
    H2PhotoDissRate = baseRate * (radField/1.7) * scatter(xl,av) * H2SelfShielding(NH2,dopplerWidth,radWidth)
end function H2PhotoDissRate

real(dp) function COPhotoDissRate(NH2,NCO,radField,av)
    real(dp), intent(in) :: NH2,NCO,radField,av  !column densities of H2 and CO
    real(dp) :: ssf,lba,sca
    !calculate photodissociation rates for co (species # nco; reaction
    !# nrco) according to van dishoeck and black (apj 334, p771 (1988))
    !cocol is the co column density (in cm-2); the scaling of the pdrate
    !by a factor of 1.8 is due to the change from draine's is uv field
    !to the habing field
    ssf = COSelfShielding(NH2,NCO)
    lba = lbar(NCO,NH2)
    sca = scatter(lba,av)

    !The reason why rad is divided by 1.7 is that the alphas are for Draine and the rad is in
    !Habing units
    COPhotoDissRate = (2.0d-10) * (radfield/1.7) * ssf * sca
end function COPhotoDissRate

function cIonizationRate(alpha,gamma,gasTemp,NC,NH2,av,radfield) result(RATE)
   real(DP) :: RATE
   real(DP), intent(in) :: alpha,gamma,gasTemp,NC,NH2,av,radField
   real(DP) :: TAUC

!  Calculate the optical depth in the CI absorption band, accounting
!  for grain extinction and shielding by CI and overlapping H2 lines
!  1.1D-17 seems the value of the ionization cross-section of C
   TAUC=gamma*av+1.1D-17*NC+(0.9D0*gasTemp**0.27D0*(NH2/1.59D21)**0.45D0)
!  Calculate the CI photoionization rate
   RATE=alpha*(radfield/1.7)*EXP(-TAUC)
end function cIonizationRate

!-----------------------------------------------------------------------
!  H2 line self-shielding, adopting the treatment of
!  Federman, Glassgold & Kwan (1979, ApJ, 227, 466)
!-----------------------------------------------------------------------
real(dp) function H2SelfShielding(NH2,dopplerWidth,radWidth)
    real(dp), intent(in) :: NH2,dopplerWidth,radWidth
    real(dp) ::  r, sj, sr, t, u, taud
    real(dp), parameter :: FPARA=0.5,FOSC  = 1.0d-2
    !--------------------------------------------------------------
    !taud = opt. depth at line center (assuming ortho:para h2=1)
    !pi**0.5 * e2 / (m(electr) * c) = 1.5e-2 cm2/s

    taud  = FPARA * NH2 * 1.5e-2 * FOSC / dopplerWidth

    !calculate doppler contribution of self shielding function sj
    if (taud == 0.0d0) then
       sj = 1.0d0
    else if (taud < 2.0d0) then
       sj = exp(-0.6666667d0*taud)
    else if (taud < 10.0d0) then
       sj = 0.638d0*taud**(-1.25d0)
    else if (taud < 100.0d0) then
       sj = 0.505d0*taud**(-1.15d0)
    else
       sj = 0.344d0*taud**(-1.0667d0)
    end if

    !calculate wing contribution of self shielding function sr
    !IF (taud.lt.0.0d0)  taud=0.0d0
    if (radWidth == 0.0d0) then
       sr = 0.0d0
    else
       r  = radWidth/(1.7724539d0*dopplerWidth)
       t  = 3.02d0 * ((r*1.0d+03)**(-0.064d0))
       u  = SQRT(taud*r)/t
       sr = r/(t*SQRT(0.78539816d0+u**2.0))
    end if

    !calculate total self shielding function fgk
    H2SelfShielding = sj + sr
end function H2SelfShielding



!calculate the influence of dust extinction (g=0.8, omega=0.3)
!wagenblast&hartquist, mnras237, 1019 (1989)
real(dp) function scatter(x1,av)


!---------------------------------------------------------------------
!         i/o variables type declaration
!       scat   : factor describing the influence of grain scattering
!                on the uv flux dependent on the total h number
!                density and wavelength of the uv radiation
!       x1      : wavelength (in angstrom)
!       cdntot : total h number density (in cm-2)
!
!         program variables
!       av     : visual extinction in magnitudes (cf. savage et al.,
!                 1977 apj 216, p.291)
!        expo   : exponent
!        i      : loop index
!        tl     : tau(lambda)
!        tv     : tau(visual=5500a)
!        xlambda : function which calculates tl/tv
!        c(0)   : c(0) * exp(-k(0)*tau) : (rel.) intensity
!                 decrease for 0<=tau<=1 caused by grain
!                 scattering with g=0.8, omega=0.3
!                 (approximation)
!        c(i)   : sum (c(i) * exp(-k(i)*tau)) i=1,5  (rel.)
!                 intensity decrease for 1<=tau<=oo caused by
!                 grain scattering with g=0.8, omega=0.3.
!                 (cf. flannery, roberge, and rybicki 1980,
!                 apj 236, p.598).
!        k(0)   : see c0
!        k(i)   : see ci
!---------------------------------------------------------------------

!   i/o variables type declaration
    real(DP), intent(in) :: x1,av
    !
    !program variables type declaration
    real(dp), dimension(6) :: c=(/1.0d0,2.006d0,-1.438d0,7.364d-01,-5.076d-01,-5.920d-02/)
    real(dp), dimension(6) ::  k1=(/7.514d-01,8.490d-01,1.013d0,1.282d0,2.005d0,5.832d0/)
    real(dp)  :: expo, tl, tv
    integer :: i

    !optical depth in the visual
    tv = av/ 1.086d0

    !make correction for the wavelength considered
    tl = tv * xlambda(x1)

    !calculate attuentuation  due to dust scattering
    scatter = 0.0d0
    if (tl<1.0d0) then
        expo = k1(1)*tl
        if (expo<100.0d0) then
            scatter = c(1) * EXP(-expo)
        end if
    else
        do i=2,6
            expo = k1(i)*tl
            if (expo<100.0d0) then
            scatter = scatter + c(i)*EXP(-expo)
            end if
        end do
    end if

end function scatter

!=======================================================================
!
!  Determine the ratio of the optical depth at a given wavelength to
!  that at visual wavelength (lambda=5500 Angstrom) using the extinction curve of
!  Savage & Mathis (1979, ARA&A, 17, 73, Table 2)
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  LAMBDA  = wavelength (in Angstrom)
!
!  Program variables:
!  XLAMBDA       = Ratio of tau(wavelength)/tau(V) at the desired wavelength
!                  by 1D spline interpolation over the grid values
!  XLAMBDA_GRID  = tau(wavelength)/tau(V) values, determined by dividing the
!                  A_wavelength/E(B-V) values from Savage & Mathis (1979) by
!                  an assumed reddening of R=AV/E(B-V)=3.1
!  XLAMBDA_DERIV = 2nd derivative of XLAMBDA_GRID values from SPLINE
!  LAMBDA_GRID   = wavelengths (in Angstrom) listed in Table 2
!  NUM_LAMBDA    = number of wavelengths
!  START         = .TRUE. when XLAMBDA is first called
!
!  Functions called:
!  SPLINE = second derivative of the supplied 1D function (in spline.f90)
!  SPLINT = 1-dimensional cubic spline interpolated value (in spline.f90)
!
!-----------------------------------------------------------------------
real(dp) function xlambda(LAMBDA)
    real(dp), intent(in) :: LAMBDA

    real(dp) :: LAMBDA_VALUE

    if (START) then
      call SPLINE(LAMBDA_GRID,XLAMBDA_GRID,NUM_LAMBDA,1.0D30,1.0D30,XLAMBDA_DERIV)
    end if

    LAMBDA_VALUE=LAMBDA
    if(LAMBDA<LAMBDA_GRID(1)) LAMBDA_VALUE=LAMBDA_GRID(1)
    if(LAMBDA>LAMBDA_GRID(NUM_LAMBDA)) LAMBDA_VALUE=LAMBDA_GRID(NUM_LAMBDA)

    call SPLINT(LAMBDA_GRID,XLAMBDA_GRID,XLAMBDA_DERIV,NUM_LAMBDA,LAMBDA_VALUE,xlambda)
    if(XLAMBDA<0.0D0) XLAMBDA=0.0D0

end function xlambda

!-----------------------------------------------------------------------
!self-shielding of CO due to 12CO self-shielding and H2 screening
!Use Van Dishoeck & Black APJ 334, 1988 for value
!-----------------------------------------------------------------------


real(dp) function COSelfShielding(NH2,NCO)
    real(dp), intent(in) :: NH2, NCO
    real(dp) :: lognco,lognh2

    if (startr)  then
        call splie2(NCO_GRID,NH2_GRID,SCO_GRID,DIMCO,DIMH2,SCO_DERIV)
        startr = .false.
    end if

    lognco = dlog10(NCO+1.0)
    lognh2 = dlog10(NH2+1.0)

    if (lognco<NCO_GRID(1))      lognco = NCO_GRID(1)
    if (lognh2<NH2_GRID(1))      lognh2 = NH2_GRID(1)
    if (lognco>NCO_GRID(DIMCO))  lognco = NCO_GRID(DIMCO)
    if (lognh2>NH2_GRID(DIMH2))  lognh2 = NH2_GRID(DIMH2)

    call splin2(NCO_GRID,NH2_GRID,SCO_GRID,SCO_DERIV,DIMCO,DIMH2,lognco,&
    &               lognh2,COSelfShielding)
    COSelfShielding = 10.0d0**COSelfShielding
end function COSelfShielding



real(dp) function lbar(u,w)
!calculate lambda bar (in a) according to equ. 4 of van dishoeck
!and black, apj 334, p771 (1988)
! --------------------------------------------------------------
!       i/o parameter
!       u : co column density in (cm-2)
!       w : h2 column density in (cm-2)

!        program variables
!        lu : log10(co column density in cm-2)
!        lw : log10(h2 column density in cm-2)

!--------------------------------------------------------------
    !i/o parameter type declaration
    real(dp)  :: u, w, lu, lw

    lu = dlog10(dabs(u)+1.0d0)
    lw = dlog10(dabs(w)+1.0d0)

    lbar = (5675.0d0 - 200.6d0*lw) - (571.6d0 - 24.09d0*lw)*lu +&
    &(18.22d0 - 0.7664d0*lw)*lu**2

    !lbar represents the mean of the wavelengths of the 33
    !dissociating bands weighted by their fractional contribution
    !to the total rate of each depth. lbar cannot be larger than
    !the wavelength of band 33 (1076.1a) and not be smaller than
    !the wavelength of band 1 (913.6a).
    if (lbar>1076.1d0)  lbar = 1076.1d0
    if (lbar< 913.6d0)  lbar =  913.6d0
end function lbar

subroutine splie2(x1a,x2a,ya,m,n,y2a)
!given an m by n tabulated function ya, and tabulated indepen-
!dent variables x1a (m values) and x2a (n values), this routine
!constructs one-dimensional natural cubic splines of the rows
!of ya and returns the second-derivatives in the array y2a.
!(copied from numerical recipes)

!--------------------------------------------------------------
!i/o parameter and program variables
          integer           :: nn
          parameter         (nn=100)
          integer           :: m, n, j, k
          real(dp)  :: x1a(m), x2a(n), ya(m,n), y2a(m,n), ytmp(nn),&
     &                      y2tmp(nn)
!--------------------------------------------------------------
    do j=1,m
        do k=1,n
            ytmp(k) = ya(j,k)
        end do
        !values 1.0d30 signal a natural spline.
        call spline(x2a,ytmp,n,1.0d30,1.0d30,y2tmp)
        do k=1,n
            y2a(j,k) = y2tmp(k)
        end do
    end do
!==============================================================
end subroutine splie2

subroutine spline(x,y,n,yp1,ypn,y2)

!calculate cubic spline for a set of points (x,y)

!(cf. "numerical recipes" 3.3 : routine spline)
!given arrays x and y of length n containing a tabulated
!function, i.e. y(i) = f(x(i)), with x(1) < x(2) < ... < x(n),
!and given values yp1 and ypn for the first derivative of the
!interpolating function at points 1 and n, respectively, this
!routine returns an array y2 of length n which contains the
!second derivatives of the interpolating function at the
!tabulated points x(i). if yp1 and/or ypn are equal to 1.0e+30
!or larger, the routine is signaled to set the corresponding
!boundary condition for a natural spline, with zero second
!derivative on that boundary.

!--------------------------------------------------------------
!i/o parameter
!x   : vector for independent variable x; dimension x(n)
!y   : vector for x-dependent variable y; dimension y(n)
!n   : dimension of vectors containing the tabulated function
!yp1 : 1. derivative of the interpolating function at point 1
!ypn : 1. derivative of the interpolating function at point n
!y2  : 2. derivative of the interpolating function
!--------------------------------------------------------------

!i/o parameter type declaration
    integer           :: n
    real(dp)  :: x(n), y(n), yp1, ypn, y2(n)

!program variables type declaration
    integer           :: i, k
    real(dp)  :: p, qn, sig, u(100), un
!--------------------------------------------------------------

    if (yp1 >= 1.0d30) then
    !the lower boundary condition is set either to be
    !"natural"
        y2(1) =  0.0d0
        u(1)  =  0.0d0
    else
    !or else to have a specified first derivative.
        y2(1) = -0.5d0
        u(1)  = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if

    !this is the decomposition loop of the tridiagonal algorithm.
    !y2 and u are used for temporary storage of decomposed factors.
    do  i=2,n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1) + 2.0d0
        y2(i) = (sig-1.0d0)/p
        u(i)  = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
        &                /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do

    if (ypn >= 1.0d30) then
    !     the upper boundary condition is set either to be
    !     "natural"
        qn = 0.0d0
        un = 0.0d0
    else
    !     or else to have a specified first derivative.
        qn = 0.5d0
        un = (3.0d0/(x(n)-x(n-1))) *&
    &              (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if

    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)

    !this is the backsubstitution loop of the tridiagonal algorithm
    do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    end do
end subroutine SPLINE

subroutine splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
!given x1a, x2a, ya, m, n as described in splie2 and y2a as
!produced by that routine; and given a desired interpolating
!point x1, x2; this routine returns an interpolated function
!value y by bicubic spline interpolation.

!--------------------------------------------------------------
!i/o parameter and program variables type declaration
integer           :: nn
parameter         (nn=100)
integer           :: m, n, j, k
real(dp)  :: x1a(m), x2a(n), ya(m,n), y2a(m,n), ytmp(nn),&
&                      y2tmp(nn), yytmp(nn), x1, x2, y
!--------------------------------------------------------------

! perform m evaluations of the row splines constructed by splie2
!using the one-dimensional spline evaluator splint.
    do j=1,m
        do k=1,n
            ytmp(k)  = ya(j,k)
            y2tmp(k) = y2a(j,k)
        end do
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
    end do
!construct the one-dimensional column spline and evaluate it.
    call spline(x1a,yytmp,m,1.0d30,1.0d30,y2tmp)
    call splint(x1a,yytmp,y2tmp,m,x1,y)
end subroutine splin2

subroutine splint(xa,ya,y2a,n,x,y)
save
!cubic spline interpolation

!(cf. "numerical recipes" 3.3 : routine splint, and 3.4.
!routine hunt)
!given the arrays xa and ya of length n, which tabulate a
!function (with the xa(i)'s in order), and given the array y2a,
!which is the output of routine cubspl, and given a value x,
!this routine returns a cubic-spline interpolated value y.

!--------------------------------------------------------------
!-i/o parameters
!-xa  : vector for independent variable x; dimension xa(n)
!-ya  : vector for x-dependent variable y; dimension ya(n)
!-y2a : 2. derivative of the interpolating function; dim. y2a(n)
!-n   : dimension of input vectors
!-x   : x value for which y is to be interpolated
!-y   : result of interpolation
!--------------------------------------------------------------

!--------------------------------------------------------------
    !i/o parameter type declaration
    integer           :: n,nstore
    real(dp)  :: x, xa(n), y, ya(n), y2a(n)

    !program variables type declaration
    integer           :: inc, jhi, jlo, jm
    real(dp)  :: h, a, b
    logical           :: ascnd

    !find interval xa(jlo) <= x <= xa(jlo+1) = xa(jhi)
    !ascnd is true if ascending order of table, false otherwise
    ascnd = xa(n)>xa(1)
    if (jlo<=0 .or. jlo>n) then
    !input guess not useful. go immediately to bisection.
        jlo = 0
        jhi = n+1
    else
        !set the hunting increment
        inc = 1

        !hunt up if ascending array or down if descending.
        if (x>=xa(jlo) .eqv. ascnd) then
            !hunt up:
            jhi=jlo+inc
            if (jhi > n) then
                !done hunting since off end of table
                jhi=n+1
            else
                !nstore is a work around for old 'go to' logic, if jhi exceeds n, that is fine
                !but the do while loop will break so jhi equals n temporarily and nstore holds
                !real value until we exit loop.
                nstore=1
                do while (((x>=xa(jhi)) .eqv. ascnd) .and. (jhi < n))
                    !not done hunting
                    jlo=jhi
                    !so double increment
                    inc=inc+inc
                    !try again
                    jhi=jlo+inc
                    if (jhi > n) then
                        jhi=n
                        nstore=n+1
                    end if
                end do
                if (nstore == n+1) jhi=nstore
            !done hunting, value bracketed.
            end if
        else
            jhi = jlo
            !hunt down:
            jlo = jhi-inc
            if (jlo < 1) then
            jlo=0
            else
                nstore=1
                do while (((x<xa(jlo)) .eqv. ascnd) .and. (jlo > 1))
                    !not done hunting,
                    jhi = jlo
                    !so double the increment
                    inc = inc+inc
                    !and try again.
                    jlo = jhi-inc
                    if (jlo < 1) then
                        jlo=1
                        nstore=0
                    end if
                end do
                if (nstore == 0) jlo=nstore
            end if
            !done hunting, since off end of table.
        end if
    end if
    do while (jhi-jlo/=1)
    !hunt is done, so begin final bisection phase:
        jm = (jhi+jlo)/2
        if (x>xa(jm) .eqv. ascnd) then
           jlo = jm
        else
           jhi = jm
        end if
    end do

    if (jlo==0)  then
        jlo = 1
        jhi = 2
    end if

    !jlo and jhi now bracket the input value of x.
    !cubic spline polynomial is now evaluated.
    if (jlo==n)  then
        jlo = n-1
        jhi = n
    end if
    h = xa(jhi) - xa(jlo)
    a = (xa(jhi) - x) / h
    b = (x - xa(jlo)) / h
    y = a*ya(jlo) + b*ya(jhi) +&
    &  ((a**3-a)*y2a(jlo) + (b**3-b)*y2a(jhi)) * (h**2)/6.0d0
end subroutine splint

end module photoreactions


