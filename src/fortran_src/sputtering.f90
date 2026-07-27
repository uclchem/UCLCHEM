!Sputtering module for the shock parameterizations
!Takes a velocity from the shock module and uses the shock sputtering
!treatment described in  Jimenez-Serra et al. 2008 A&A 482
!http://adsabs.harvard.edu/abs/2008A&A...482..549J
! to calculate the amount of material removed from dust grains
!
! Additionally uses Guillet et al. 2011 (http://www.aanda.org/10.1051/0004-6361/201015973)
! result that for shocks >19km/s you get vaporization which sends dust grain material
! into the gas phase. We make use of this by setting 19km/s as limit above which
! refratory dust grain material is sputtered.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module sputtering
    use constants
    use DEFAULTPARAMETERS
    use f2py_constants
    !f2py INTEGER, parameter :: dp
    use network
    use physicscore, only: timeInYears
    use SurfaceReactions, only: GAS_DUST_DENSITY_RATIO,GRAIN_RADIUS
    implicit none

    integer :: projectiles(6)
    real(dp) :: sConst,eta,epso
    !Speed at which refractory species are also removed from dust grains during sputtering. 19.0 km/s taken from Guillet et al. 2011. (see above)
    real(dp), parameter :: VAPORIZE_SPEED=19.0
    integer, allocatable :: sputters(:),gasSputters(:)
contains
    subroutine sputteringSetup
        integer :: i,j,k,new_size
        logical :: found

        !Need abundance of major projectiles to get sputtering rate
        projectiles=(/nh2,nhe,nc,no,nsi,nco/)


        !check for refractory species and create a sublist of ice species
        !that only includes species not in refractory list
        !This allows us to sputter only volatile species when shockvelocity is low
        if (ALLOCATED(sputters)) deallocate(sputters,gasSputters)

        if (refractoryList(1) > 0) then
            new_size=size(iceList)-size(refractoryList)
            allocate(sputters(new_size),gasSputters(new_size))
            k=1
            do i=1,size(iceList)
                found=.false.
                do j=1,size(refractoryList)
                    if (iceList(i) == refractoryList(j)) found=.true.
                end do
                if (.not. found) then
                    sputters(k)=iceList(i)
                    gasSputters(k)=gasIceList(i)
                    k=k+1
                end if

            end do
        else
            allocate(sputters(size(iceList)))
            allocate(gasSputters(size(iceList)))
            gasSputters=gasIceList
            sputters=iceList
        end if
    end subroutine sputteringSetup

    subroutine sputterIces(abund,shockVel,gasTemp,density,timeDelta)
      ! Sputter ices following Jimenez-Serra 2008 treatment
      ! Args:
      !     abund: abundances of all species
      !     shockVel: relative velocity of dust and colliding gas (shockVel in c-shock)
      real(dp), intent(inout) :: abund(nspec+1)
      real(dp), intent(in) :: shockVel,gasTemp,density,timeDelta
      real(dp) :: sputterRate,abundChangeFrac,totalMantle, grainNumberDensity
      integer :: iSpec

      !Constant relating mass and speed of projectile to energy
      sConst=(shockVel*shockVel*km*km)/(2.0*gasTemp*K_BOLTZ)
      sConst=sqrt(sConst)

      !loop over projectile species and get rates of change of mantle for each, summing them
      sputterRate=0.0
      do iSpec=1,SIZE(projectiles)  !!!! Make projectiles array in initialize
          sputterRate=sputterRate+iceYieldRate(mass(projectiles(iSpec))*MH,density*abund(projectiles(iSpec)),gasTemp)
      end do

      grainNumberDensity=density/GAS_DUST_DENSITY_RATIO
      !Total rate/cm3 (ie released particles /cm3/s) is sputterRate (per grain) multiplied by grain number density
      sputterRate=sputterRate*grainNumberDensity

      !integrate that forward from currentTimeOld to currentTime. to get total number of particles released
      abundChangeFrac=sputterRate*(timeDelta)  !/density
      !I think that commented out dens is required for units. However, sputtering doesn't happen if it is uncommented
      !and sputtering matches Jimenez-Serra et al. 2008 curves when it's commented out.


      !if M particles are released and there are N particles on the grain total
      !then a species with X particles on the grain will release M*(X/N)
      !this is M/N and we'll multiply by X below
      totalMantle=sum(abund(iceList))
      abundChangeFrac=abundChangeFrac/totalMantle
      if (abundChangeFrac > 1.0d0) abundChangeFrac=1.0d0
      if (abundChangeFrac < 0.00d0) abundChangeFrac=0.0d0

      !write(87,*) timeInYears,shockVel,abundChangeFrac,timeDelta/SECONDS_PER_YEAR
      !multiply M/N by x and add to gas phase
      if (shockVel >= VAPORIZE_SPEED) then
        do iSpec = 1, SIZE(iceList)
          abund(gasIceList(iSpec)) = abund(gasIceList(iSpec)) + abundChangeFrac * abund(iceList(iSpec))
          abund(iceList(iSpec))    = abund(iceList(iSpec))    * (1.0d0 - abundChangeFrac)
        end do
      else
        do iSpec = 1, SIZE(sputters)
          abund(gasSputters(iSpec)) = abund(gasSputters(iSpec)) + abundChangeFrac * abund(sputters(iSpec))
          abund(sputters(iSpec))    = abund(sputters(iSpec))    * (1.0d0 - abundChangeFrac)
        end do
      end if
  end subroutine sputterIces

  !Function calculates rate of change of ice mantle abundance of a species!
  !due to the impact of molecules of a given mass. actual rate is         !
  !proportional to projectile abundance                                   !
  function iceYieldRate(projectileMass,projectileDensity,gasTemp)
      real(dp) :: iceYieldRate
      real(dp) :: projectileMass,projectileDensity,gasTemp
      real(dp) :: lowerLimit,upperLimit,s

      real(dp), parameter :: iceBindingEnergy=0.53*1.6d-12,targetMass=18.0*MH
      real(dp), parameter :: iceYieldEfficiency=0.8  !

      !eta is effectively reduced mass of the collision
      eta=4.0*iceYieldEfficiency*projectileMass*targetMass*((projectileMass+targetMass)**(-2.0))
      epso=max(1.0,4.0*eta)
      s=sConst*sqrt(projectileMass)

      !Lower limit is xth in Jimenez-Serra et al. 2008
      lowerLimit=sqrt(epso*iceBindingEnergy/(eta*K_BOLTZ*gasTemp))

      !Upper limit is just where the integrand goes to zero
      upperLimit=iceYieldIntegralLimit(lowerLimit,projectileMass,gasTemp)
      !calculate eq B.1 from Jimenez-Serra et al. 2008
      if ((upperlimit-lowerLimit) > 1d-4) then
          !first get integral from Eq B.1 including 1/s factor
          iceYieldRate=trapezoidIntegrate(lowerLimit,upperLimit,projectileMass,gasTemp)/s
          !multiply through by constants
          iceYieldRate=iceYieldRate*GRAIN_RADIUS*GRAIN_RADIUS*sqrt(8.0*K_BOLTZ*gasTemp*pi/projectileMass)
          !need projectile number density
          iceYieldRate=iceYieldRate*projectileDensity
      else
          iceYieldRate=0.0
      end if
  end function iceYieldRate


  !Function calculates integrand from Eq B.1 of Jimenez-Serra et al. 2008 !
  !                                                                       !
  !Inputs are mass of projectile and x. Returns value of integrand at x   !
  !allowing trapezium rule to integrate from xth to infinity              !
  function iceYieldIntegrand(x,projectileMass,gasTemp)
      real(dp) :: iceYieldIntegrand,x,projectileMass,gasTemp
      real(dp) :: yield,s,eps

      real(dp), parameter :: yieldConst=8.3d-4
      real(dp), parameter :: iceBindingEnergy=0.53*1.6d-12

      !this is s from exp(x+s) in eq B.1, varies only with mass so constant precalculated in initialize
      s=sConst*sqrt(projectileMass)

      !epsilon is calculated from inmpact energy (Ep)
      eps=(x**2)*K_BOLTZ*gasTemp
      !and some other factors
      eps=eta*eps/iceBindingEnergy
      !this yield is for ice averaged over all angles. There's a different one for cores (Appendix B Jimenez-Serra 2008)
      !it's 2 times the normal incidence yield, but there's a factor of 0.5 in integrand so we drop both
      yield=yieldConst*((eps-epso)**2)/(1.0+((eps/30.0)**(1.3333)))
      iceYieldIntegrand=yield*(x**2)*(exp(-((x-s)**2))-exp(-((x+s)**2)))
  end function iceYieldIntegrand

  !Function to calculate the upper limit beyond which there's no point
  !evaluating the ice yield integrand. Ie trapezoids from upper limit to
  !upperlimit+dx will have an area~0
  function iceYieldIntegralLimit(xth,projectileMass,gasTemp)
      real(dp) :: iceYieldIntegralLimit
      real(dp) :: xth,projectileMass,gasTemp
      integer :: i
      i=1
      !Take upperlimit to be half way between lowerLimit and 1000.
      iceYieldIntegralLimit=xth+(1d3-xth)*(0.5**i)
      !decrease upper limit for as long as f(upperlimit) is <1.0e-20 and
      !difference between lower and upper limit is not zero.
      do while (iceYieldIntegrand(iceYieldIntegralLimit,projectileMass,gasTemp) < 1d-200 .and.&
          & (iceYieldIntegralLimit-xth) > 1.0d-3)
              i=i+1
              iceYieldIntegralLimit=xth+(1d3-xth)*(0.5**i)
      end do
  end function iceYieldIntegralLimit



  function trapezoidIntegrate(lowerLimit,upperlimit,projectileMass,gasTemp)
     !Subroutine that calculates an integral using the trapezoidal method. It repeatedly
    !tries smaller intervals until the area is small enough to be accurate.
    ! It used to take a function to integrate as an argument but I removed it
    ! since we just want to integrate iceYieldIntegrand anyway.
      real(dp) :: trapezoidIntegrate
      integer,parameter :: JMAX=25
      real(dp) :: lowerLimit,upperlimit,projectileMass,gasTemp
      real(dp), parameter :: tolerance=1.0e-3
      integer :: j
      real(dp) :: olds
      olds=-1.0e30
      do j=1,JMAX
          call trapzd(lowerLimit,upperlimit,trapezoidIntegrate,j,projectileMass,gasTemp)
          if (abs(trapezoidIntegrate-olds)<=tolerance*abs(olds)) return
          olds=trapezoidIntegrate
      end do
  end function trapezoidIntegrate

  subroutine trapzd(a,b,s,n,projectileMass,gasTemp)
    ! Subroutine to integrate with trapazoidal rule using n intervals.
      integer :: n
      real(dp) :: a,b,s,projectileMass,gasTemp
      integer :: it,j
      real(dp) :: del,sum,tnm,x
      if(n==1) then
          s=0.5*(b-a)*(iceYieldIntegrand(a,projectileMass,gasTemp)+&
          &iceYieldIntegrand(b,projectileMass,gasTemp))
      else
          it=2**(n-2)
          tnm=it
          del=(b-a)/tnm
          x=a+0.5*del
          sum=0.0
          do j=1,it
              sum=sum+iceYieldIntegrand(x,projectileMass,gasTemp)
              x=x+del
          end do
          s=0.5*(s+(b-a)*sum/tnm)
      end if
  end subroutine trapzd
end module sputtering
