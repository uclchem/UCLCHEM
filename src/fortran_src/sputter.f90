!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mock subroutine for cshock that will sputter the ices based on Jimenez-Serra!
! 2008 paper.                                                                 !
!
! TO DO:          
!  - actual loop to get masses and abundances for each projectile and do calculation!
!  - More careful decisions over what should be module variable rather than Function
!  - Tests against the izaskun paper    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sputtering(abund)
      DOUBLE PRECISION :: abund(nspec+1,points)
      INTENT(INOUT) :: abund
      DOUBLE PRECISION :: sputterRate=0.0
      INTEGER :: iSpec

      !loop over projectile species and get rates of change of mantle for each, summing them
      DO iSpec=1,SIZE(projectiles) !!!! Make projectiles array in initialize
            sputterRate+=iceYieldRate(mass(projectiles(i),density(dstep)*abund(projectiles(i,dstep))))
      END DO

      !Total rate is sputterRate (per grain) multiplied by grain number density
      

      !integrate that forward (check currentTime/targetTime) think I want to go
      !from currentTimeOld to currentTime.

      !Should get a portion of current mantle to remove, add that to gas phase.
      !eg abund(gasGrainList)+=sputterFrac*abund(grainList)
            !abund(grainList)-=sputterFrac*abund(grainList)
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function calculates rate of change of ice mantle abundance of a species!
!due to the impact of molecules of a given mass. actual rate is         !
!proportional to projectile abundance                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION iceYieldRate(projectileMass,projectileAbundance)
      DOUBLE PRECISION projectileMass,projectileAbundance,sConst,driftVel
      DOUBLE PRECISION, PARAMETER :: iceBindingEnergy=0.53*1.6d-12
      DOUBLE PRECISION :: lowerLimit,upperLimit

      sConst=(driftVel*driftVel)/(2.0*temp(dstep)*K_BOLTZ_CGS)
      sConst=sqrt(sConst)

      !eta is effectively reduced mass of the collision
      eta=4.*iceYieldEfficiency*projectileMass*targetMass*(projectileMass+targetMass)**(-2)
      epso=max(1.,4.*eta)

      !Lower limit is xth in Jimenez-Serra et al. 2008
      lowerLimit=sqrt(epso*iceBindingEnergy/(eta*K_BOLTZ_CGS*temp(dstep)))
      !Upper limit is just where the integrand goes to zero
      upperLimit=iceYieldIntegralLimit(lowerLimit,projectileMass)

      iceYieldRate=trapezoidIntegrate(iceYieldIntegrand,lowerLimit,upperLimit)
      iceYieldRate=iceYieldRate*grainRadius*grainRadius*sqrt(8.0*K_BOLTZ_CGS*temp(dstep)*pi/projectileMass)
      iceYieldRate=iceYieldRate*projectileAbundance
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function calculates integrand from Eq B.1 of Jimenez-Serra et al. 2008 !
!                                                                       !
!Inputs are mass of projectile and x. Returns value of integrand at x   !
!allowing trapezium rule to integrate from xth to infinity              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
FUNCTION iceYieldIntegrand(x,projectileMass)
      DOUBLE PRECISION :: iceYieldIntegrand,x,projectileMass,
      DOUBLE PRECISION :: yield

      DOUBLE PRECISION, PARAMETER :: yieldConst=8.2d-4,
      DOUBLE PRECISION, PARAMETER :: iceYieldEfficiency=0.8 !
      DOUBLE PRECISION, PARAMETER :: targetMass=18.0 !Mass of H2O in amu
      DOUBLE PRECISION, PARAMETER :: iceBindingEnergy=0.53*1.6d-12

      !this is s from exp(x+s), varies only with mass so should probably precalculate and save constant 
      s=sConst*sqrt(projectileMass)
      E=(x**2)*K_BOLTZ_CGS*temp
      eps=eta*E/iceBindingEnergy
      !this yield is for ice. There's a different one for cores (Appendix B Jimenez-Serra 2008)
      yield=yieldConst*(eps-epso)**2/(1.+(eps/30.)**(4./3.))
      func_iceHe=yield*(x**2)*(exp(-(x-s)**2)-exp(-(x+s)**2))
END FUNCTION iceYieldIntegrand

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function to calculate the upper limit beyond which there's no point   !
!evaluating the ice yield integrand. Ie trapezoids from upper limit to !
!upperlimit+dx will have an area~0                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION iceYieldIntegralLimit(xth,projectileMass)
      DOUBLE PRECISION iceYieldIntegralLimit,xth,projectileMass
      INTEGER :: i=1
      DO WHILE (iceYieldIntegrand(iceYieldIntegralLimit,projectileMass) .lt. 1e-200)
            iceYieldIntegralLimit=xth+(1d3-xth)*(0.5**i)
      END DO
END FUNCTION iceYieldIntegralLimit

