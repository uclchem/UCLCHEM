SUBROUTINE calculateReactionRates
!Assuming the user has temperature changes or uses the desorption features of phase 1, these need working out on a timestep by time step basis
    DO j=1,nreac
        !This case structure looks at the reaction type. species-species happens in default.
        !Other cases are special reactions, particularly desorption events (photons, CRs etc)
        SELECT CASE (re2(j))
        !Cosmic ray reactions            
        CASE ('CRP')
            rate(j) = alpha(j)*zeta
        !UV photons, radfield has (factor of 1.7 conversion from habing to Draine)
        CASE ('PHOTON')
            rate(j) = alpha(j)*dexp(-gama(j)*av(dstep))*radfield/1.7
            !co photodissoction number is stored as nrco
            IF (re1(j).eq.'CO') THEN
                IF(p1(j).eq.'O' .and. p2(j).eq.'C') nrco=j
                IF(p1(j).eq.'C' .and. p2(j).eq.'O') nrco=j
            ENDIF
        !cosmic ray induced photon
        CASE ('CRPHOT')
            rate(j)=alpha(j)*gama(j)*1.0/(1.0-omega)*zeta*(temp(dstep)/300)**beta(j)
        !freeze out only happens if fr>0 and depending on evap choice 
        CASE ('FREEZE')             
            IF (evap .ne. 0 .or. fr .eq. 0.0) then
                rate(j)=0.0
            ELSE
                IF (re1(j).eq."E-") THEN
                    cion=1.0+16.71d-4/(GRAIN_RADIUS*temp(dstep))
                    rate(j)=4.57d4*alpha(j)*GRAIN_AREA*fr*cion
                ELSE
                    DO i=1,nspec-1
                        IF (specname(i).eq.re1(j)) THEN
                            IF (beta(j).eq.0.0 ) THEN
                                !taken from Rawlings et al. 1992
                                rate(j)=4.57d4*alpha(j)*dsqrt(temp(dstep)/mass(i))*GRAIN_AREA*fr
                            ELSE
                                !Make rates sets beta=1 for ion freeze out. this catches that and
                                !freezes differently
                                cion=1.0+16.71d-4/(GRAIN_RADIUS*temp(dstep))
                                rate(j)=4.57d4*alpha(j)*dsqrt(temp(dstep)/mass(i))*GRAIN_AREA*fr*cion
                            ENDIF
                        ENDIF
                    END DO
                END IF
            ENDIF
        !The below desorption mechanisms are from Roberts et al. 2007 MNRAS with
        !the addition of direct UV photodesorption. DESOH2,DESCR1,DEUVCR
        CASE ('DESOH2')
            IF (desorb .eq. 1 .and. h2desorb .eq. 1&
            & .and. gama(j) .le. ebmaxh2 .and.&
            &  mantle(dstep) .ge. 1.0d-30) THEN
                !Epsilon is efficieny of this process, number of molecules removed per event
                !h2form is formation rate of h2, dependent on hydrogen abundance. 
                rate(j) = epsilon*h2form*abund(nh,dstep)*1.0/mantle(dstep)
            ELSE
                rate(j) = 0.0
            ENDIF
        CASE ('DESCR')
            IF (desorb .eq. 1 .and. crdesorb .eq. 1&
            &.and.mantle(dstep).ge. 1d-30&
            &.and. gama(j) .le. ebmaxcr) THEN
                !4*pi*zeta = total CR flux. 1.64d-4 is iron to proton ratio of CR
                !as iron nuclei are main cause of CR heating.
                !GRAIN_AREA is the total area per hydrogen atom. ie total grain area per cubic cm when multiplied by density.
                !phi is efficieny of this reaction, number of molecules removed per event.
                rate(j) = 4.0*pi*zeta*1.64d-4*(GRAIN_AREA)*&
                      &(1.0/mantle(dstep))*phi
            ELSE
                rate(j) = 0.0
            ENDIF
        CASE ('DEUVCR')
            IF (desorb .eq. 1 .and. uvcr .eq. 1 &
             &.and. gama(j) .le. ebmaxuvcr .and. mantle(dstep) .ge. 1.0d-30) THEN
                !4.875d3 = photon flux, Checchi-Pestellini & Aiello (1992) via Roberts et al. (2007)
                !UVY is yield per photon.
                rate(j) = GRAIN_AREA*uvy*4.875d3*zeta*(1.0/mantle(dstep))
                !additional factor accounting for UV desorption from ISRF. UVCREFF is ratio of 
                !CR induced UV to ISRF UV.
                rate(j) = rate(j) * (1+(radfield/uvcreff)*(1.0/zeta)*dexp(-1.8*av(dstep)))
            ELSE
                rate(j) = 0.0
            ENDIF
        CASE DEFAULT
            !Reactions on surface can be treated considering diffusion of reactants
            !See work of David Quenard 2017 Arxiv:1711.05184
            !abstracted to functions below for ease of reading
            IF (re3(j).eq.'DIFF' .or. re3(j) .eq. 'CHEMDES') THEN
                rate(j) = diffusionReactionRate()
            ELSE
            !--------------------------------------------------------------------------------------------------------
            !Basic gas phase reactions 
                rate(j) = alpha(j)*((temp(dstep)/300.)**beta(j))*dexp(-gama(j)/temp(dstep))
            ENDIF

        END SELECT
    END DO

    h2dis=h2d()
    rate(nrco)=knrco()
END SUBROUTINE calculateReactionRates


!----------------------------------------------------------------------------------------------------
!Reactions on the surface treated by evaluating diffusion rates across grains and accounting
!For competition with chemical desorption. Products remain bound ('DIFF') or are desorbed ('CHEMDES')
!Assuming Eb = 0.5 Ed. Units of s-1. 
!David Quenard 2017 Arxiv:1711.05184
!----------------------------------------------------------------------------------------------------
double precision FUNCTION diffusionReactionRate()
    double precision :: diffuseRate,activationBarrier,reducedMass,tunnelProb
    double precision :: diffuseProb,desorbProb
    integer :: ns,index1,index2
    character(len=10) :: reactant1,reactant2


    !H and H2 can also be reactants but are not in grain list
    SELECT CASE (re1(j))
        CASE ('H')
            reactant1="#H"
        CASE ('H2')
            reactant1="#H2"
        CASE DEFAULT
            reactant1=re1(j)
    END SELECT

    SELECT CASE (re2(j))
        CASE ('H')
            reactant2="#H"
        CASE ('H2')
            reactant2="#H2"
        CASE DEFAULT
            reactant2=re2(j)
    END SELECT

    ! now loop through mantle species and match reactants to species
    !index is the index of the species in the grain array. NOT in the full species array.
    DO i=lbound(grainList,1),ubound(grainList,1)
        IF (specname(grainList(i)) .eq. reactant1) index1 = i
        IF (specname(grainList(i)) .eq. reactant2) index2 = i
    END DO

    !Hasegawa 1992 diffusion rate. Rate that two species diffuse and meet on grain surface
    diffuseRate = vdiff(index1)*dexp(-0.5*bindingEnergy(index1)/temp(dstep))
    diffuseRate = (diffuseRate+ (vdiff(index2)*dexp(-0.5*bindingEnergy(index2)/temp(dstep))))/NUM_SITES_PER_GRAIN

    !Calculate classical activation energy barrier exponent
    activationBarrier = gama(j)/temp(dstep)

    !Calculate quantum activation energy barrier exponent
    reducedMass = mass(grainList(index1)) * mass(grainList(index2)) / (mass(grainList(index1)) + mass(grainList(index2)))
    tunnelProb = 2.0d0 *CHEMICAL_BARRIER_THICKNESS/REDUCED_PLANCK * dsqrt(2.0d0*AMU*reducedMass*K_BOLTZ*gama(j))

    !Choose fastest between classical and tunnelling
    IF (activationBarrier.GT.tunnelProb) activationBarrier=tunnelProb
    !set activationBarrier to probability of reaction Ruaud+2016
    activationBarrier=dexp(-activationBarrier)

    ! Keff from Garrod & Pauly 2011 and Ruaud+2016
    ! Actual reaction probability is Preac/(Preac+Pevap+Pdiffuse), accounting for the other possible processes
    IF(DIFFUSE_REACT_COMPETITION) THEN
       activationBarrier = max(vdiff(index1),vdiff(index2)) * activationBarrier       
       !probability a reactant will just desorb
       desorbProb = vdiff(index1)*dexp(-bindingEnergy(index1)/temp(dstep))
       desorbProb = desorbProb + vdiff(index2)*dexp(-bindingEnergy(index2)/temp(dstep)) 
       !Probability reactants wills diffuse onwards when they meet
       diffuseProb = (diffuseRate) * NUM_SITES_PER_GRAIN
       !Therefore, total probability the reactants that meet will react
       activationBarrier = activationBarrier/(activationBarrier + desorbProb + diffuseProb)
    END IF
    
    diffusionReactionRate=alpha(j) * diffuseRate * activationBarrier* GAS_DUST_DENSITY_RATIO / dens(dstep)

    !Now adjust for fraction of this reaction's products that will desorb due to energy released
    IF (re3(j).eq.'DIFF') THEN
        diffusionReactionRate = diffusionReactionRate * (1.0-desorptionFraction(j,index1,index2))
    ELSE IF(re3(j).eq.'CHEMDES') THEN
        diffusionReactionRate = diffusionReactionRate * desorptionFraction(j,index1,index2)
    ENDIF
END FUNCTION diffusionReactionRate

! ---------------------------------------------------------------------
!  Chemical Reactive Desorption (CRD)
! David Quenard 2017 Arxiv:1711.05184
! From Minissalle+ 2016 and Vasyunin+ 2016
! ---------------------------------------------------------------------
double precision FUNCTION desorptionFraction(j,reactIndex1,reactIndex2)
    integer :: j,reactIndex1,reactIndex2,degreesOfFreedom
    integer :: productIndex(4)

    double precision :: deltaEnthalpy,maxBindingEnergy,epsilonCd,productEnthalpy
    double precision, parameter :: EFFECTIVE_SURFACE_MASS = 120.0

    
    !Get indices of grain surface version of products products 
    productIndex=0
    productIndex = 0.0
    maxBindingEnergy=0.0
    productEnthalpy=0.0

    IF (re3(j).eq.'DIFF') THEN
        DO i=lbound(grainList,1),ubound(grainList,1)
            IF (specname(grainList(i)) .eq. p1(j)) productIndex(1)=grainList(i)
            !Go through grain list and try to find species in product list
            !If it has a binding energy larger than largest energy found so far, update maxBindingEnergy
            IF (specname(grainList(i)) .eq. p1(j)) THEN
                productIndex(1) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF

            IF (specname(grainList(i)) .eq. p2(j)) THEN
                productIndex(2) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF

            IF (specname(grainList(i)) .eq. p3(j)) THEN
                productIndex(3) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF

            IF (specname(grainList(i)) .eq. p4(j)) THEN
                productIndex(4) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF

        END DO
    ELSE IF (re3(j).eq.'CHEMDES') THEN
        DO i=lbound(gasGrainList,1),ubound(gasGrainList,1)
            IF (specname(gasGrainList(i)) .eq. p1(j)) THEN
                productIndex(1) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF
            IF (specname(gasGrainList(i)) .eq. p2(j)) THEN
                productIndex(2) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF
            IF (specname(gasGrainList(i)) .eq. p3(j)) THEN
                productIndex(3) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF
            IF (specname(gasGrainList(i)) .eq. p4(j)) THEN
                productIndex(4) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF
        END DO
    ENDIF
    
    !
    !epsilonCd is the fraction of kinetic energy kept my the product when it collides with grain surface
    epsilonCd = mass(productIndex(1)) + mass(productIndex(2)) + mass(productIndex(3)) + mass(productIndex(4))
    epsilonCd = ((epsilonCd - EFFECTIVE_SURFACE_MASS) / (epsilonCd + EFFECTIVE_SURFACE_MASS))**2
    
    !Now calculate the change in enthalpy of the reaction.
    deltaEnthalpy= formationEnthalpy(reactIndex1)+formationEnthalpy(reactIndex2)-productEnthalpy
    
    !Convert from kcal to J, from J to K and from moles-1 to reactions-1
    deltaEnthalpy = deltaEnthalpy*4.184d03/(1.38054D-23*6.02214129d23)
    ! Total energy change includes activation energy of the reaction !
    deltaEnthalpy = deltaEnthalpy + gama(j)


    IF (deltaEnthalpy.eq.0.00) deltaEnthalpy = 1e-30 

    !Degrees of freedom = 3 * number of atoms in the molecule
    degreesOfFreedom = atomCounts(productIndex(1))
    if (productIndex(2).NE.0) degreesOfFreedom = max(degreesOfFreedom,atomCounts(productIndex(2)))
    if (productIndex(3).NE.0) degreesOfFreedom = max(degreesOfFreedom,atomCounts(productIndex(3)))
    if (productIndex(4).NE.0) degreesOfFreedom = max(degreesOfFreedom,atomCounts(productIndex(4)))                    
    degreesOfFreedom = 3 * degreesOfFreedom
        
    desorptionFraction = dexp((-maxBindingEnergy*real(degreesOfFreedom)) / (epsilonCd * deltaEnthalpy))
    
   IF (deltaEnthalpy.lt.0.d0) THEN        !< If reaction is endothermic, no CRD
        desorptionFraction = 0.d0
    END IF
    
    IF (GRAINS_HAVE_ICE) THEN
        desorptionFraction = desorptionFraction/10    !< See Minisalle et al. 2016 for icy grain surface.
        ! Special case of OH+H, O+H, N+N on ices, see same paper
        if (re1(j).eq.'#N'.and.re2(j).eq.'#N') desorptionFraction = 0.5
        if ((re1(j).eq.'#O'.and.re2(j).eq.'H') .or. (re1(j).eq.'H'.and.re2(j).eq.'#O')) desorptionFraction = 0.3
        if ((re1(j).eq.'#OH'.and.re2(j).eq.'H') .or. (re1(j).eq.'H'.and.re2(j).eq.'#OH')) desorptionFraction = 0.25
    ENDIF
END FUNCTION desorptionFraction


!-----------------------------------------------------------------------------------------
!h2d() and knrco() recalculate h2 and CO photodissociation using a self-shielding 
!treatment from van dishoeck and black (apj 334, p771 (1988))
!All other functions are helpers for this purpose
!Code by Serena Viti
!-----------------------------------------------------------------------------------------

double precision FUNCTION h2d()
    double precision ::ch2

    !h2col is h2 column density. Half total column density for uniform sphere.
    !Sum of shells for multidepth point (worked out in chem.f90.updateChemistry)

    !taud = opt. depth at line centre (assum. ortho:parah2=1)
    !pi**0.5 * e2 / (m(electr) * c) = 1.5e-2 cm2/s
    taud  = 0.5 * h2col * 1.5e-2 * fosc / dopw

    !c Here, the constant 5.1e-11 is by assuming rad is in Habing [double check]
    h2d = 5.1d-11 * radfield * scat(xl) * fgk(taud)
END FUNCTION h2d 

double precision FUNCTION knrco()
    double precision ::ch2,ssf,lba,sca,chtot
    !double precision :: ssfco,lbar,scat,chtot
    
    !calculate photodissociation rates for co (species # nco; reaction
    !# nrco) according to van dishoeck and black (apj 334, p771 (1988))
    !cocol is the co column density (in cm-2); the scaling of the pdrate
    !by a factor of 1.8 is due to the change from draine's is uv field
    !to the habing field
    ssf = ssfco(h2col)
    lba = lbar(cocol,h2col)
    sca = scat(lba)

    !The reason why rad is divided by 1.7 is that the alphas are for Draine and the rad is in 
    !Habing units
    
    knrco = (2.d-10) * radfield/1.7 * ssf * sca
END FUNCTION knrco
        
double precision function fgk(taud)
!calculates the line self shielding function
!described in federman et al. apj vol.227 p.466.
   
!--------------------------------------------------------------
!input parameters
!dopw : doppler   line width (in s-1) from sr setpmp
!radw : radiative line width (in s-1) from sr setpmp
!taud : optical depth at line center  from sr setequ
!   
!program variables
!fgk  : total self shielding function containing the radiative
!       and the doppler contribution.
!r    : parameter r  (eq. a2) in federman's paper
!sj   : parameter jd (eq. a8) in federman's paper
!       doppler contribution to self shielding function
!sr   : parameter jr (eq. a9) in federman's paper
!       radiative contribution to self shielding function
!t    : parameter t1 (eq. a6) in federman's paper
!u    : parameter u1 (eq. a6) in federman's paper
!-----------------------------------------------------------

    !program variables type declaration
    double precision taud
    double precision  r, sj, sr, t, u
    !--------------------------------------------------------------


    !calculate wing contribution of self shielding function sr
    IF (taud.lt.0.0d0)  taud=0.0d0
    IF (radw .eq. 0.0d0) then
       sr = 0.0d0
    ELSE
       r  = radw/(1.7724539d0*dopw)
       t  = 3.02d0 * ((r*1.0d+03)**(-0.064d0))
       u  = ((taud*r)**0.5d0)/t
       sr = ((0.78539816d0+u**2)**(-0.5d0))*r/t
    ENDIF
       
    !calculate doppler contribution of self shielding function sj
    IF (taud .eq. 0.0d0) THEN
       sj = 1.0d0
    ELSE IF (taud .lt. 2.0d0) THEN
       sj = exp(-0.6666667d0*taud)
    ELSE IF (taud .lt. 10.0d0) THEN
       sj = 0.638d0*taud**(-1.25d0)
    ELSE IF (taud .lt. 100.0d0) THEN
       sj = 0.505d0*taud**(-1.15d0)
    ELSE
       sj = 0.344d0*taud**(-1.0667d0)
    END IF
    
    !calculate total self shielding function fgk
    fgk = sj + sr
END FUNCTION fgk
     
double precision FUNCTION scat(x1)
!calculate the influence of dust extinction (g=0.8, omega=0.3)
!wagenblast&hartquist, mnras237, 1019 (1989)

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
!        xlamda : function which calculates tl/tv
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
    double precision  cdntot,x1
    !   
    !program variables type declaration
    double precision, dimension(6) :: c=(/1.0d0,2.006d0,-1.438d0,7.364d-01,-5.076d-01,-5.920d-02/)
    double precision, dimension(6) ::  k1=(/7.514d-01,8.490d-01,1.013d0,1.282d0,2.005d0,5.832d0/)
    double precision  expo, tl, tv

    !calculate visual extinction
    !total column density of hydrogen nuclei / e(b-v) = 5.8e+21
    !atoms cm-2 mag-1 (bohlin, savage, and drake 1978; apj 224,132)
    !for lambda**-1 scattering : r = av / e(b-v) = 3.6 .
    
    !optical depth in the visual
    tv = av(dstep)/ 1.086d0
      
    !make correction for the wavelength considered
    tl = tv * xlamda(x1)
       
    !calculate scat
    scat = 0.0d0
    IF (tl.lt.1.0d0) THEN
        expo = k1(1)*tl
        IF (expo.lt.35.0d0) THEN
            scat = c(1) * dexp(-expo)
        ENDIF
    ELSE
        DO i=2,6
            expo = k1(i)*tl
            IF (expo.lt.35.0d0) THEN
            scat = scat + c(i)*dexp(-expo)
            ENDIF
        END DO
    ENDIF

END FUNCTION scat

double precision function xlamda(x)
!calculate  xlamda : =  tau(lambda) / tau(visual)
!   
! --------------------------------------------------------------
!         i/o parameter
!xlamda : =  tau(lambda) / tau(visual);
!         tau(lambda) is the opt. depth for dust extinction at
!         wavelength x (cf. b.d.savage and j.s.mathis, annual
!         review of astronomy and astrophysics vol.17(1979),p.84)
!x      : wavelength in angstrom
! --------------------------------------------------------------

          !i/o parameter type declaration
          double precision  x
  
          if (x.lt.  910.0d0) then
               xlamda = 5.76d0
          else if (x.lt.  950.0d0) then
               xlamda = 5.76d0 - 1.45d-02*(x-  910.0d0)
          else if (x.lt. 1000.0d0) then
               xlamda = 5.18d0 - 1.06d-02*(x-  950.0d0)
          else if (x.lt. 1050.0d0) then
               xlamda = 4.65d0 - 9.68d-03*(x- 1000.0d0)
          else if (x.lt. 1110.0d0) then
               xlamda = 4.16d0 - 7.26d-03*(x- 1050.0d0)
          else if (x.lt. 1180.0d0) then
               xlamda = 3.73d0 - 4.61d-03*(x- 1110.0d0)
          else if (x.lt. 1250.0d0) then
               xlamda = 3.40d0 - 4.15d-03*(x- 1180.0d0)
          else if (x.lt. 1390.0d0) then
               xlamda = 3.11d0 - 2.67d-03*(x- 1250.0d0)
          else if (x.lt. 1490.0d0) then
               xlamda = 2.74d0 - 1.10d-03*(x- 1390.0d0)
          else if (x.lt. 1600.0d0) then
               xlamda = 2.63d0 - 8.80d-05*(x- 1490.0d0)
          else if (x.lt. 1700.0d0) then
               xlamda = 2.62d0 - 8.06d-04*(x- 1600.0d0)
          else if (x.lt. 1800.0d0) then
               xlamda = 2.54d0 - 3.87d-04*(x- 1700.0d0)
          else if (x.lt. 1900.0d0) then
               xlamda = 2.50d0 + 8.07d-04*(x- 1800.0d0)
          else if (x.lt. 2000.0d0) then
               xlamda = 2.58d0 + 2.00d-03*(x- 1900.0d0)
          else if (x.lt. 2100.0d0) then
               xlamda = 2.78d0 + 2.29d-03*(x- 2000.0d0)
          else if (x.lt. 2190.0d0) then
               xlamda = 3.01d0 + 1.22d-03*(x- 2100.0d0)
          else if (x.lt. 2300.0d0) then
               xlamda = 3.12d0 - 2.35d-03*(x- 2190.0d0)
          else if (x.lt. 2400.0d0) then
               xlamda = 2.86d0 - 2.81d-03*(x- 2300.0d0)
          else if (x.lt. 2500.0d0) then
               xlamda = 2.58d0 - 2.29d-03*(x- 2400.0d0)
          else if (x.lt. 2740.0d0) then
               xlamda = 2.35d0 - 1.46d-03*(x- 2500.0d0)
          else if (x.lt. 3440.0d0) then
               xlamda = 2.00d0 - 5.99d-04*(x- 2740.0d0)
          else if (x.lt. 4000.0d0) then
               xlamda = 1.58d0 - 2.88d-04*(x- 3440.0d0)
          else if (x.lt. 4400.0d0) then
               xlamda = 1.42d0 - 2.42d-04*(x- 4000.0d0)
          else if (x.lt. 5500.0d0) then
               xlamda = 1.32d0 - 2.93d-04*(x- 4400.0d0)
          else if (x.lt. 7000.0d0) then
               xlamda = 1.00d0 - 1.68d-04*(x- 5500.0d0)
          else if (x.lt. 9000.0d0) then
               xlamda = 0.75d0 - 1.32d-04*(x- 7000.0d0)
          else if (x.lt.12500.0d0) then
               xlamda = 0.48d0 - 5.81d-05*(x- 9000.0d0)
          else if (x.lt.22000.0d0) then
               xlamda = 0.28d0 - 1.66d-05*(x-12500.0d0)
          else if (x.lt.34000.0d0) then
               xlamda = 0.12d0 - 5.91d-06*(x-22000.0d0)
          else
               xlamda = 0.05d0 - 5.16d-11*(x-34000.0d0)
          ENDIF

END FUNCTION xlamda

double precision FUNCTION ssfco(ch2)
!calculates self-shielding factors for 12co transitions due to
!12co self-shielding and h2 screening. values given in table 5
!of van dishoeck and black, apj 334, p771 (1988) are used for
!2dim spline interpolation. nco and nh2 are the 12co resp. h2
!column densities for which the self-shielding factor is
!interpolated.

!common parameter
!ssfcor : corates(7,6) : self-shielding factors for 12co
!                     considering h2 line overlap according
!                     to van dishoeck and black, apj 334, p771
!                     (1988). logarithmic (base 10) values for
!                     the self-shielding factors are stored.
!                     1st index : variation of co column dens.
!                     2nd index : variation of h2 column dens.
!y2r(7,6) :   2nd derivative of rates from sr splie2
!ncogr(7) :   12co column density grid (log10 values
!             of column densities in cm-2)
!nh2gr(6) :   h2 column density grid (log10 values
!             of column densities in cm-2)
!dimco :      dimension of ncogr
!dimh2 :      dimension of nh2gr
!startr :     is .true. when ssfcor is entered first
!actual values are supplied by block data routine
   
!program variables type declaration
    double precision  ch2, lognco, lognh2
    if (startr)  THEN
        call splie2(ncogr,nh2gr,corates,dimco,dimh2,y2r)
        startr = .false.
    endif
   
    lognco = dlog10(max(cocol,1.0d5))
    lognh2 = dlog10(max(ch2,1.0d10))
       
    if (lognco.lt.ncogr(1))      lognco = ncogr(1)
    if (lognh2.lt.nh2gr(1))      lognh2 = nh2gr(1)
    if (lognco.gt.ncogr(dimco))  lognco = ncogr(dimco)
    if (lognh2.gt.nh2gr(dimh2))  lognh2 = nh2gr(dimh2)
       
    call splin2(ncogr,nh2gr,corates,y2r,dimco,dimh2,lognco,&
    &               lognh2,ssfco)
    ssfco = 10.0d0**ssfco
END FUNCTION ssfco
   
   
   
double precision FUNCTION lbar(u,w)
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
    double precision  u, w, lu, lw

    lu = dlog10(dabs(u)+1.0d0)
    lw = dlog10(dabs(w)+1.0d0)
    
    lbar = (5675.0d0 - 200.6d0*lw) - (571.6d0 - 24.09d0*lw)*lu +&
    &(18.22d0 - 0.7664d0*lw)*lu**2
       
    !lbar represents the mean of the wavelengths of the 33
    !dissociating bands weighted by their fractional contribution
    !to the total rate of each depth. lbar cannot be larger than
    !the wavelength of band 33 (1076.1a) and not be smaller than
    !the wavelength of band 1 (913.6a).
    if (lbar.gt.1076.1d0)  lbar = 1076.1d0
    if (lbar.lt. 913.6d0)  lbar =  913.6d0
END FUNCTION lbar

SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
!given an m by n tabulated function ya, and tabulated indepen-
!dent variables x1a (m values) and x2a (n values), this routine
!constructs one-dimensional natural cubic splines of the rows
!of ya and returns the second-derivatives in the array y2a.
!(copied from numerical recipes)

!--------------------------------------------------------------
!i/o parameter and program variables
          integer           nn
          parameter         (nn=100)
          integer           m, n, j, k
          double precision  x1a(m), x2a(n), ya(m,n), y2a(m,n), ytmp(nn),&
     &                      y2tmp(nn)
!--------------------------------------------------------------
    DO j=1,m
        DO k=1,n
            ytmp(k) = ya(j,k)
        END DO
        !values 1.0d30 signal a natural spline.
        call spline(x2a,ytmp,n,1.0d30,1.0d30,y2tmp)
        DO k=1,n
            y2a(j,k) = y2tmp(k)
        END DO
    END DO
return
!==============================================================
END SUBROUTINE splie2
  
SUBROUTINE spline(x,y,n,yp1,ypn,y2)

!calculate cubic spline for a set of points (x,y)
 
!(cf. "numerical recipes" 3.3 : routine spline)
!given arrays x and y of length n containing a tabulated
!function, i.e. y(i) = f(x(i)), with x(1) < x(2) < ... < x(n),
!and given values yp1 and ypn for the first derivative of the
!interpolating function at points 1 and n, respectively, this
!routine returns an array y2 of length n which contains the
!second derivatives of the interpolating function at the
!tabulated points x(i). if yp1 and/or ypn are equal to 1.0e+30
!or larger, the routine is signalled to set the corresponding
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
    integer           n
    double precision  x(n), y(n), yp1, ypn, y2(n)
   
!program variables type declaration
    integer           i, k
    double precision  p, qn, sig, u(100), un
!--------------------------------------------------------------

    IF (yp1 .ge. 1.0d30) THEN
    !the lower boundary condition is set either to be
    !"natural"
        y2(1) =  0.0d0
        u(1)  =  0.0d0
    ELSE
    !or else to have a specified first derivative.
        y2(1) = -0.5d0
        u(1)  = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    ENDIF
   
    !this is the decomposition loop of the tridiagonal algorithm.
    !y2 and u are used for temporary storage of decomposed factors.
    DO  i=2,n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1) + 2.0d0
        y2(i) = (sig-1.0d0)/p
        u(i)  = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
        &                /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    END DO
       
    IF (ypn .ge. 1.0d30) THEN
    !     the upper boundary condition is set either to be
    !     "natural"
        qn = 0.0d0
        un = 0.0d0
    ELSE
    !     or else to have a specified first derivative.
        qn = 0.5d0
        un = (3.0d0/(x(n)-x(n-1))) *&
    &              (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    ENDIF
       
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
       
    !this is the backsubstitution loop of the tridiagonal algorithm
    DO k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    END DO
END SUBROUTINE SPLINE
 
SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
!given x1a, x2a, ya, m, n as described in splie2 and y2a as
!produced by that routine; and given a desired interpolating
!point x1, x2; this routine returns an interpolated function
!value y by bicubic spline interpolation.
  
!--------------------------------------------------------------
!i/o parameter and program variables type declaration
integer           nn
parameter         (nn=100)
integer           m, n, j, k
double precision  x1a(m), x2a(n), ya(m,n), y2a(m,n), ytmp(nn),&
&                      y2tmp(nn), yytmp(nn), x1, x2, y
!--------------------------------------------------------------

! perform m evaluations of the row splines constructed by splie2
!using the one-dimensional spline evaluator splint.
    DO j=1,m
        DO k=1,n
            ytmp(k)  = ya(j,k)
            y2tmp(k) = y2a(j,k)
        END DO
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
    END DO
!construct the one-dimensional column spline and evaluate it.
    call spline(x1a,yytmp,m,1.0d30,1.0d30,y2tmp)
    call splint(x1a,yytmp,y2tmp,m,x1,y)
END SUBROUTINE splin2   
     
SUBROUTINE splint(xa,ya,y2a,n,x,y)
SAVE
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
    integer           n,nstore
    double precision  x, xa(n), y, ya(n), y2a(n)
       
    !program variables type declaration
    integer           inc, jhi, jlo, jm
    double precision  h, a, b
    logical           ascnd
   
    !find interval xa(jlo) <= x <= xa(jlo+1) = xa(jhi)
    !ascnd is true if ascending order of table, false otherwise
    ascnd = xa(n).gt.xa(1)
    if (jlo.le.0 .or. jlo.gt.n) then
    !input guess not useful. go immediately to bisection.
        jlo = 0
        jhi = n+1
    ELSE           
        !set the hunting increment
        inc = 1
        
        !hunt up if ascending array or down if descending.
        IF (x.ge.xa(jlo) .eqv. ascnd) THEN
            !hunt up:
            jhi=jlo+inc
            IF (jhi .gt. n) THEN
                !done hunting since off end of table
                jhi=n+1
            ELSE
                !nstore is a work around for old 'go to' logic, if jhi exceeds n, that is fine
                !but the do while loop will break so jhi equals n temporarily and nstore holds
                !real value until we exit loop.
                nstore=1
                DO WHILE (((x.ge.xa(jhi)) .eqv. ascnd) .and. (jhi .lt. n))
                    !not done hunting
                    jlo=jhi
                    !so double increment
                    inc=inc+inc
                    !try again
                    jhi=jlo+inc
                    IF (jhi .gt. n) THEN
                        jhi=n
                        nstore=n+1
                    ENDIF
                END DO
                IF (nstore .eq. n+1) jhi=nstore
            !done hunting, value bracketed.
            END IF      
        ELSE
            jhi = jlo
            !hunt down:
            jlo = jhi-inc
            IF (jlo .lt. 1) THEN
            jlo=0
            ELSE
                nstore=1
                DO WHILE (((x.lt.xa(jlo)) .eqv. ascnd) .and. (jlo .gt. 1))
                    !not done hunting,
                    jhi = jlo
                    !so double the increment
                    inc = inc+inc
                    !and try again.
                    jlo = jhi-inc
                    IF (jlo .lt. 1) THEN
                        jlo=1
                        nstore=0
                    END IF
                END DO
                IF (nstore .eq. 0) jlo=nstore
            END IF
            !done hunting, since off end of table.
        ENDIF
    END IF   
    DO WHILE (jhi-jlo.ne.1)  
    !hunt is done, so begin final bisection phase:
        jm = (jhi+jlo)/2
        IF (x.gt.xa(jm) .eqv. ascnd) THEN
           jlo = jm
        ELSE
           jhi = jm
        ENDIF    
    END DO
       
    IF (jlo.eq.0)  THEN
        jlo = 1
        jhi = 2
    ENDIF
       
    !jlo and jhi now bracket the input value of x.
    !cubic spline polynomial is now evaluated.
    IF (jlo.eq.n)  THEN
        jlo = n-1
        jhi = n
    ENDIF
    h = xa(jhi) - xa(jlo)
    a = (xa(jhi) - x) / h
    b = (x - xa(jlo)) / h
    y = a*ya(jlo) + b*ya(jhi) +&
    &  ((a**3-a)*y2a(jlo) + (b**3-b)*y2a(jhi)) * (h**2)/6.0d0
END SUBROUTINE splint
