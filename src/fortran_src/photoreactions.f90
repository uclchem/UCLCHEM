module photoreactions
    use constants, only: dp, HABING_TO_DRAINE
    use DEFAULTPARAMETERS
    !f2py INTEGER, parameter :: dp

    implicit none

    public

    real(dp) :: UV_FAC=3.02_dp

    ! Below are arrays for self-shielding of CO and H2
    ! They cannot be parameters, because they can be adjusted
    ! at runtime from within python.
    logical :: start=.true.
    integer, parameter :: MAX_NUM_LAMBDA = 30
    integer :: NUM_LAMBDA = MAX_NUM_LAMBDA  ! Initialize with maximum number of wavelengths
    real(dp), dimension(MAX_NUM_LAMBDA) :: LAMBDA_GRID=(/ &
      &  910.0_dp, 950.0_dp,1000.0_dp,1050.0_dp,1110.0_dp, &
      & 1180.0_dp,1250.0_dp,1390.0_dp,1490.0_dp,1600.0_dp, &
      & 1700.0_dp,1800.0_dp,1900.0_dp,2000.0_dp,2100.0_dp, &
      & 2190.0_dp,2300.0_dp,2400.0_dp,2500.0_dp,2740.0_dp, &
      & 3440.0_dp,4000.0_dp,4400.0_dp,5500.0_dp,7000.0_dp, &
      & 9000.0_dp,12500.0_dp,22000.0_dp,34000.0_dp,1.0e9_dp/)
    real(dp), dimension(MAX_NUM_LAMBDA) :: XLAMBDA_GRID=(/ &
      & 5.76_dp,5.18_dp,4.65_dp,4.16_dp,3.73_dp, &
      & 3.40_dp,3.11_dp,2.74_dp,2.63_dp,2.62_dp, &
      & 2.54_dp,2.50_dp,2.58_dp,2.78_dp,3.01_dp, &
      & 3.12_dp,2.86_dp,2.58_dp,2.35_dp,2.00_dp, &
      & 1.58_dp,1.42_dp,1.32_dp,1.00_dp,0.75_dp, &
      & 0.48_dp,0.28_dp,0.12_dp,0.05_dp,0.00_dp/)
    real(dp), dimension(MAX_NUM_LAMBDA) :: XLAMBDA_DERIV

    logical :: startr=.true.

    !  12CO line shielding data from van Dishoeck & Black (1988, ApJ, 334, 771, Table 5)
    integer, parameter ::  DIMCO=8, DIMH2=6
    real(dp), dimension(DIMCO) :: NCO_GRID=(/12.0_dp,13.0_dp,14.0_dp,15.0_dp,16.0_dp,17.0_dp,18.0_dp,19.0_dp/)
    real(dp), dimension(DIMH2) :: NH2_GRID=(/18.0_dp,19.0_dp,20.0_dp,21.0_dp,22.0_dp,23.0_dp/)
    real(dp), dimension(DIMCO,DIMH2) :: SCO_GRID=RESHAPE((/ &
      &  0.000e+00_dp,-1.408e-02_dp,-1.099e-01_dp,-4.400e-01_dp,-1.154e+00_dp,-1.888e+00_dp,-2.760e+00_dp,-4.001e+00_dp, &
      & -8.539e-02_dp,-1.015e-01_dp,-2.104e-01_dp,-5.608e-01_dp,-1.272e+00_dp,-1.973e+00_dp,-2.818e+00_dp,-4.055e+00_dp, &
      & -1.451e-01_dp,-1.612e-01_dp,-2.708e-01_dp,-6.273e-01_dp,-1.355e+00_dp,-2.057e+00_dp,-2.902e+00_dp,-4.122e+00_dp, &
      & -4.559e-01_dp,-4.666e-01_dp,-5.432e-01_dp,-8.665e-01_dp,-1.602e+00_dp,-2.303e+00_dp,-3.146e+00_dp,-4.421e+00_dp, &
      & -1.303e+00_dp,-1.312e+00_dp,-1.367e+00_dp,-1.676e+00_dp,-2.305e+00_dp,-3.034e+00_dp,-3.758e+00_dp,-5.077e+00_dp, &
      & -3.883e+00_dp,-3.888e+00_dp,-3.936e+00_dp,-4.197e+00_dp,-4.739e+00_dp,-5.165e+00_dp,-5.441e+00_dp,-6.446e+00_dp/),&
        & (/DIMCO,DIMH2/))
    real(dp), dimension(DIMCO,DIMH2) :: SCO_DERIV

    real(dp) :: ICE_GAS_PHOTO_CROSSSECTION_RATIO = 0.3_dp  ! Kalvans 2018
contains

!photodissociation rate of H2  accounting for self-shielding
function getH2PhotoDissRateConstant(NH2,radField,av,turbVel) result(H2PhotoDissRateConstant)
    !H2 Column density to surface, UV at that surface, visual extinction to surface, and turbulent velocity of cloud
    real(dp), intent(in) :: NH2 ,radField,av ,turbVel
    real(dp) :: H2PhotoDissRateConstant
    !unshielded reaction rate, characteristic wavelendth of radiation, radiative linewidth
    real(dp), parameter :: baseRateConstant=5.18e-11_dp,xl=1000.0,radWidth=8.0e7_dp
    real(dp) :: dopplerWidth

    dopplerWidth=turbVel/(xl*1.0e-8_dp)
    H2PhotoDissRateConstant = baseRateConstant * (radField*HABING_TO_DRAINE) * calculate_grain_scatter(xl,av) * &
        calculateH2SelfShielding(NH2,dopplerWidth,radWidth)
end function getH2PhotoDissRateConstant

function getCOPhotoDissRateConstant(NH2,NCO,radField,av) result(COPhotoDissRateConstant)
    real(dp), intent(in) :: NH2,NCO,radField,av  !column densities of H2 and CO
    real(dp) :: COPhotoDissRateConstant

    real(dp) :: ssf,lba,sca
    !calculate photodissociation rates for co (species # nco; reaction
    !# nrco) according to van dishoeck and black (apj 334, p771 (1988))
    !cocol is the co column density (in cm-2); the scaling of the pdrate
    !by a factor of 1.8 is due to the change from draine's is uv field
    !to the habing field
    ssf = calculateCOSelfShielding(NH2,NCO)
    lba = calculate_lbar(NCO,NH2)
    sca = calculate_grain_scatter(lba,av)

    !The reason why rad is divided by 1.7 is that the alphas are for Draine and the rad is in
    !Habing units
    COPhotoDissRateConstant = (2.0e-10_dp) * (radfield*HABING_TO_DRAINE) * ssf * sca
end function getCOPhotoDissRateConstant

pure function getCarbonIonizationRateConstant(alpha,gamma,gasTemp,NC,NH2,av,radfield) &
        result(carbonIonizationRateConstant)
   real(dp), intent(in) :: alpha,gamma,gasTemp,NC,NH2,av,radField
   real(dp) :: carbonIonizationRateConstant
   real(dp) :: TAUC

!  Calculate the optical depth in the CI absorption band, accounting
!  for grain extinction and shielding by CI and overlapping H2 lines
!  1.1e-17_dp seems the value of the ionization cross-section of C
   TAUC=gamma*av+1.1e-17_dp*NC+(0.9_dp*gasTemp**0.27_dp*(NH2/1.59e21_dp)**0.45_dp)
!  Calculate the CI photoionization rate
   carbonIonizationRateConstant=alpha*(radfield*HABING_TO_DRAINE)*EXP(-TAUC)
end function getCarbonIonizationRateConstant

!-----------------------------------------------------------------------
!  H2 line self-shielding, adopting the treatment of
!  Federman, Glassgold & Kwan (1979, ApJ, 227, 466)
!-----------------------------------------------------------------------
pure function calculateH2SelfShielding(NH2,dopplerWidth,radWidth) result(H2SelfShielding)
    real(dp), intent(in) :: NH2,dopplerWidth,radWidth
    real(dp) :: H2SelfShielding

    real(dp) ::  r, sj, sr, t, u, taud
    real(dp), parameter :: FPARA = 0.5_dp
    real(dp), parameter :: FOSC  = 1.0e-2_dp
    !--------------------------------------------------------------
    !taud = opt. depth at line center (assuming ortho:para h2=1)
    !pi**0.5 * e2 / (m(electr) * c) = 1.5e-2 cm2/s

    taud  = FPARA * NH2 * 1.5e-2_dp * FOSC / dopplerWidth

    !calculate doppler contribution of self shielding function sj
    if (taud == 0.0_dp) then
       sj = 1.0_dp
    else if (taud < 2.0_dp) then
       sj = exp(-0.6666667_dp*taud)
    else if (taud < 10.0_dp) then
       sj = 0.638_dp*taud**(-1.25_dp)
    else if (taud < 100.0_dp) then
       sj = 0.505_dp*taud**(-1.15_dp)
    else
       sj = 0.344_dp*taud**(-1.0667_dp)
    end if

    !calculate wing contribution of self shielding function sr
    !IF (taud.lt.0.0_dp)  taud=0.0_dp
    if (radWidth == 0.0_dp) then
       sr = 0.0_dp
    else
       r  = radWidth/(1.7724539_dp*dopplerWidth)
       t  = 3.02_dp * ((r*1.0e+03_dp)**(-0.064_dp))
       u  = SQRT(taud*r)/t
       sr = r/(t*SQRT(0.78539816_dp+u**2.0_dp))
    end if

    !calculate total self shielding function fgk
    H2SelfShielding = sj + sr
end function calculateH2SelfShielding



!calculate the influence of dust extinction (g=0.8, omega=0.3)
!wagenblast&hartquist, mnras237, 1019 (1989)
function calculate_grain_scatter(x1,av) result(grain_scatter)
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
!        calculate_xlambda : function which calculates tl/tv
!        c(0)   : c(0) * exp(-k(0)*tau) : (rel.) intensity
!                 decrease for 0<=tau<=1 caused by grain
!                 calculate_grain_scattering with g=0.8, omega=0.3
!                 (approximation)
!        c(i)   : sum (c(i) * exp(-k(i)*tau)) i=1,5  (rel.)
!                 intensity decrease for 1<=tau<=oo caused by
!                 grain calculate_grain_scattering with g=0.8, omega=0.3.
!                 (cf. flannery, roberge, and rybicki 1980,
!                 apj 236, p.598).
!        k(0)   : see c0
!        k(i)   : see ci
!---------------------------------------------------------------------

!   i/o variables type declaration
    real(dp), intent(in) :: x1, av
    real(dp) :: grain_scatter
    !
    !program variables type declaration
    integer, parameter :: NUM_FIT_CONSTANTS = 6
    real(dp), parameter, dimension(NUM_FIT_CONSTANTS) :: c=(/1.0_dp,2.006_dp,-1.438_dp,7.364e-01_dp,-5.076e-01_dp,-5.920e-02_dp/)
    real(dp), parameter, dimension(NUM_FIT_CONSTANTS) :: k1=(/7.514e-01_dp,8.490e-01_dp,1.013_dp,1.282_dp,2.005_dp,5.832_dp/)
    real(dp)  :: expo, tl, tv
    integer :: i

    !optical depth in the visual
    tv = av / 1.086_dp

    !make correction for the wavelength considered
    tl = tv * calculate_xlambda(x1)

    !calculate attuentuation  due to dust calculate_grain_scattering
    grain_scatter = 0.0_dp
    if (tl<1.0_dp) then
        expo = k1(1)*tl
        if (expo<100.0_dp) then
            grain_scatter = c(1) * EXP(-expo)
        end if
    else
        do i=2,NUM_FIT_CONSTANTS
            expo = k1(i)*tl
            if (expo<100.0_dp) then
            grain_scatter = grain_scatter + c(i)*EXP(-expo)
            end if
        end do
    end if

end function calculate_grain_scatter

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
function calculate_xlambda(LAMBDA) result(xlambda)
    real(dp), intent(in) :: LAMBDA
    real(dp) :: xlambda

    real(dp) :: LAMBDA_VALUE

    if (START) then
      call SPLINE(LAMBDA_GRID,XLAMBDA_GRID,NUM_LAMBDA,1.0e30_dp,1.0e30_dp,XLAMBDA_DERIV)
    end if

    LAMBDA_VALUE=LAMBDA
    if(LAMBDA<LAMBDA_GRID(1)) LAMBDA_VALUE=LAMBDA_GRID(1)
    if(LAMBDA>LAMBDA_GRID(NUM_LAMBDA)) LAMBDA_VALUE=LAMBDA_GRID(NUM_LAMBDA)

    call SPLINT(LAMBDA_GRID,XLAMBDA_GRID,XLAMBDA_DERIV,NUM_LAMBDA,LAMBDA_VALUE,xlambda)
    if(xlambda<0.0_dp) xlambda=0.0_dp

end function calculate_xlambda

!-----------------------------------------------------------------------
!self-shielding of CO due to 12CO self-shielding and H2 screening
!Use Van Dishoeck & Black APJ 334, 1988 for value
!-----------------------------------------------------------------------
function calculateCOSelfShielding(NH2,NCO) result(COSelfShielding)
    real(dp), intent(in) :: NH2, NCO
    real(dp) :: COSelfShielding
    real(dp) :: lognco,lognh2

    if (startr)  then
        call splie2(NCO_GRID,NH2_GRID,SCO_GRID,DIMCO,DIMH2,SCO_DERIV)
        startr = .false.
    end if

    lognco = log10(NCO+1.0_dp)
    lognh2 = log10(NH2+1.0_dp)

    if (lognco<NCO_GRID(1))      lognco = NCO_GRID(1)
    if (lognh2<NH2_GRID(1))      lognh2 = NH2_GRID(1)
    if (lognco>NCO_GRID(DIMCO))  lognco = NCO_GRID(DIMCO)
    if (lognh2>NH2_GRID(DIMH2))  lognh2 = NH2_GRID(DIMH2)

    call splin2(NCO_GRID,NH2_GRID,SCO_GRID,SCO_DERIV,DIMCO,DIMH2,lognco,&
    &               lognh2,COSelfShielding)
    COSelfShielding = 10.0_dp**COSelfShielding
end function calculateCOSelfShielding



pure function calculate_lbar(u,w) result(lbar)
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
    real(dp), intent(in) :: u, w
    real(dp) :: lbar
    !i/o parameter type declaration
    real(dp) :: lu, lw

    lu = log10(abs(u)+1.0_dp)
    lw = log10(abs(w)+1.0_dp)

    lbar = (5675.0_dp - 200.6_dp*lw) - (571.6_dp - 24.09_dp*lw)*lu +&
    &(18.22_dp - 0.7664_dp*lw)*lu**2

    !calculate_lbar represents the mean of the wavelengths of the 33
    !dissociating bands weighted by their fractional contribution
    !to the total rate of each depth. calculate_lbar cannot be larger than
    !the wavelength of band 33 (1076.1a) and not be smaller than
    !the wavelength of band 1 (913.6a).
    if (lbar>1076.1_dp) lbar = 1076.1_dp
    if (lbar< 913.6_dp) lbar =  913.6_dp
end function calculate_lbar

pure subroutine splie2(x1a,x2a,ya,m,n,y2a)
!given an m by n tabulated function ya, and tabulated indepen-
!dent variables x1a (m values) and x2a (n values), this routine
!constructs one-dimensional natural cubic splines of the rows
!of ya and returns the second-derivatives in the array y2a.
!(copied from numerical recipes)

!--------------------------------------------------------------
    integer, intent(in) :: m, n
    real(dp), intent(in), dimension(m) :: x1a
    real(dp), intent(in), dimension(n) :: x2a
    real(dp), intent(in), dimension(m, n) :: ya
    real(dp), intent(out), dimension(m, n) :: y2a

    !i/o parameter and program variables
    integer, parameter :: nn = 100
    real(dp), dimension(nn) :: ytmp, y2tmp
    integer :: j, k
!--------------------------------------------------------------
    do j=1,m
        do k=1,n
            ytmp(k) = ya(j,k)
        end do
        !values 1.0e30_dp signal a natural spline.
        call spline(x2a,ytmp,n,1.0e30_dp,1.0e30_dp,y2tmp)
        do k=1,n
            y2a(j,k) = y2tmp(k)
        end do
    end do
!==============================================================
end subroutine splie2

pure subroutine spline(x,y,n,yp1,ypn,y2)

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
    integer, intent(in)           :: n
    real(dp), intent(in), dimension(n) :: x, y
    real(dp), intent(in) :: yp1, ypn

    real(dp), intent(out), dimension(n) :: y2

!program variables type declaration
    integer, parameter :: nn = 100
    integer           :: i, k
    real(dp)  :: p, qn, sig, un
    real(dp), dimension(nn) :: u
!--------------------------------------------------------------

    if (yp1 >= 1.0e30_dp) then
    !the lower boundary condition is set either to be
    !"natural"
        y2(1) =  0.0_dp
        u(1)  =  0.0_dp
    else
    !or else to have a specified first derivative.
        y2(1) = -0.5_dp
        u(1)  = (3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if

    !this is the decomposition loop of the tridiagonal algorithm.
    !y2 and u are used for temporary storage of decomposed factors.
    do  i=2,n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1) + 2.0_dp
        y2(i) = (sig-1.0_dp)/p
        u(i)  = (6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
        &                /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do

    if (ypn >= 1.0e30_dp) then
    !     the upper boundary condition is set either to be
    !     "natural"
        qn = 0.0_dp
        un = 0.0_dp
    else
    !     or else to have a specified first derivative.
        qn = 0.5_dp
        un = (3.0_dp/(x(n)-x(n-1))) *&
    &              (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if

    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0_dp)

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
    integer, intent(in) :: m, n
    real(dp), intent(in), dimension(m) :: x1a
    real(dp), intent(in), dimension(n) :: x2a
    real(dp), intent(in), dimension(m, n) :: ya, y2a
    real(dp), intent(in) :: x1, x2
    real(dp), intent(out) :: y

    integer, parameter :: nn = 100
    integer :: j, k
    real(dp), dimension(nn) :: ytmp, y2tmp, yytmp
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
    call spline(x1a,yytmp,m,1.0e30_dp,1.0e30_dp,y2tmp)
    call splint(x1a,yytmp,y2tmp,m,x1,y)
end subroutine splin2

subroutine splint(xa,ya,y2a,n,x,y)
save
!cubic spline interpolation

!(cf. "numerical recipes" 3.3 : routine splint, and 3.4.
!routine hunt)
!given the arrays xa and ya of length n, which tabulate a
!function (with the xa(i)'s in order), and given the array y2a,
!which is theZZ output of routine cubspl, and given a value x,
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
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in)  :: xa, ya, y2a

    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    !program variables type declaration
    integer           :: nstore, inc, jhi, jlo, jm
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
    &  ((a**3-a)*y2a(jlo) + (b**3-b)*y2a(jhi)) * (h**2)/6.0_dp
end subroutine splint

end module photoreactions
