! Chemistry module of UCL_CHEM. Contains all the core machinery of the code, not really intended to be altered.
! Use physics module to alter temp/density behaviour etc. This module should solve chemistry for a cloud of gas
MODULE chemistry
USE physics
USE dvode_f90_m

IMPLICIT NONE
    !one of three incuded files. Odes.f90 and network.f90 come from makerates. rates.f90 hides a large subroutine
    include "network.f90"
   !These integers store the array index of important species and reactions, x is for ions    
    integer :: nh,nh2,nc,ncx,no,nn,ns,nhe,nco,nmg,nf,nh2o,nsi,nsix,ncl,nclx,nch3oh,np
    integer :: nrco,njunk,evapevents,ngrainco,readAbunds
    !loop counters    
    integer :: i,j,l,writeStep,writeCounter=0

    !Array to store reaction rates
    double precision :: rate(nreac)
    
    !Option column output
    character(LEN=15),allocatable :: outSpecies(:)
    logical :: columnFlag
    integer :: nout
    integer, allocatable :: outIndx(:)


    !DLSODE variables    
    integer :: ITASK,ISTATE,NEQ,MXSTEP
    double precision :: reltol
    double precision, allocatable :: abstol(:)
    TYPE(VODE_OPTS) :: OPTIONS

    !initial fractional elemental abudances and arrays to store abundances
    double precision :: fh,fhe,fc,fo,fn,fs,fmg,fsi,fcl,fp,ff,h2col,cocol,junk1,junk2
    double precision,allocatable :: abund(:,:),mantle(:)
    
    !Variables controlling chemistry
    double precision :: radfield,zeta,fr,omega,grainArea,cion,h2form,h2dis
    double precision :: ebmaxh2,epsilon,ebmaxcrf,ebmaxcr,phi,ebmaxuvcr,uvy,uvcreff
    double precision, allocatable ::vdiff(:)

    !Variables for self-shielding of CO and H2
    !dopw = doppler width (in s-1) of a typical transition
    !(assuming turbulent broadening with beta=3e5cms-1)
    !radw = radiative line width of typ. transition (in s-1)
    !fosc = oscillator strength of a typical transition
    double precision  :: dopw=3.0e10,radw=8.0e07,xl=1000.0,fosc  = 1.0d-2,taud

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !variables for diffusion reactions on the grains, CGS unless otherwise stated.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    LOGICAL, parameter :: DIFFUSE_REACT_COMPETITION=.True., GRAINS_HAVE_ICE=.True.
    double precision, parameter :: GAS_DUST_MASS_RATIO=100.0,REDUCED_PLANCK=1.054571628d-27,AMU=1.66053892d-24
    double precision, parameter :: K_BOLTZ=1.3806588d-16,GRAIN_AREA=2.4d-22,GRAIN_RADIUS=1.d-5 
    double precision, parameter :: CHEMICAL_BARRIER_THICKNESS = 1.40d-8  !gre Parameter used to compute the probability for a surface reaction with 
    !! activation energy to occur through quantum tunneling (Hasegawa et al. Eq 6 (1992).)
    double precision, parameter :: SURFACE_SITE_DENSITY = 1.5d15 ! site density on one grain [cm-2]
    double precision, parameter :: VDIFF_PREFACTOR=2.0*K_BOLTZ*SURFACE_SITE_DENSITY/PI/PI/AMU
    double precision, parameter :: GRAIN_DENSITY = 3.0 ! Mass density of a dust grain
    double precision, parameter :: NUM_SITES_PER_GRAIN = GRAIN_RADIUS*GRAIN_RADIUS*SURFACE_SITE_DENSITY*4.0*PI
    double precision, parameter :: GAS_DUST_DENSITY_RATIO = (4.0*PI*(GRAIN_RADIUS**3)*GRAIN_DENSITY*GAS_DUST_MASS_RATIO)/(3.0 * AMU)
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !CO and H2 self-shielding
    !Used by functions in rates.f90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical :: startr=.True.
    integer,parameter :: dimco=7, dimh2=6
    double precision :: corates(7,6)=reshape((/0.000d+00, -1.408d-02, -1.099d-01, -4.400d-01,&
     &  -1.154d+00, -1.888d+00, -2.760d+00,&
     &  -8.539d-02, -1.015d-01, -2.104d-01, -5.608d-01,&
     &  -1.272d+00, -1.973d+00, -2.818d+00,&
     &  -1.451d-01, -1.612d-01, -2.708d-01, -6.273d-01,&
     &  -1.355d+00, -2.057d+00, -2.902d+00,&
     &  -4.559d-01, -4.666d-01, -5.432d-01, -8.665d-01,&
     &  -1.602d+00, -2.303d+00, -3.146d+00,&
     &  -1.303d+00, -1.312d+00, -1.367d+00, -1.676d+00,&
     &  -2.305d+00, -3.034d+00, -3.758d+00,&
     &  -3.883d+00, -3.888d+00, -3.936d+00, -4.197d+00,&
     &  -4.739d+00, -5.165d+00, -5.441d+00 /),shape(corates))
    double precision :: y2r(7,6)
    double precision :: ncogr(dimco) =(/12.0d+00, 13.0d+00, 14.0d+00, 15.0d+00,&
      &16.0d+00, 17.0d+00, 18.0d+00 /)
    double precision :: nh2gr(dimh2)=(/18.0d+00, 19.0d+00, 20.0d+00, 21.0d+00,&
       &22.0d+00, 23.0d+00 /)
CONTAINS
!This gets called immediately by main so put anything here that you want to happen before the time loop begins, reader is necessary.
    SUBROUTINE initializeChemistry
        NEQ=nspec+1
        IF (ALLOCATED(abund)) DEALLOCATE(abund,vdiff,mantle)
        ALLOCATE(abund(NEQ,points),vdiff(nspec))
        CALL reader
        !if this is the first step of the first phase, set initial abundances
        !otherwise reader will fix it
        IF (readAbunds.eq.0) THEN
            !ensure abund is initially zero
            abund= 0.
            !As default, have half in molecular hydrogen and half in atomic hydrogen
            abund(nh2,:) = 0.5*(0.5*(1.0e0-fh))
            abund(nh,:) = (0.5*(1.0e0-fh))     
            abund(nhe,:) = fhe                       
            abund(no,:) = fo  
            abund(nn,:) = fn               
            abund(ns,:) = fs
            abund(nmg,:) = fmg
            abund(nsix,:) = fsi                
            abund(nclx,:) = fcl 
            !abund(np,:) = fp
            !abund(nf,:) = ff

            !abund(nfe,:) = ffe
            !abund(nna,:) = fna
            abund(nspec+1,:)=dens      

            !Decide how much iron is initiall ionized using parameters.f90
            SELECT CASE (ion)
                CASE(0)
                    abund(nc,:)=fc
                    abund(ncx,:)=1.d-10
                CASE(1)
                    abund(nc,:)=fc*0.5
                    abund(ncx,:)=fc*0.5
                CASE(2)
                    abund(nc,:)=1.d-10
                    abund(ncx,:)=fc
            END SELECT
            abund(nspec,:)=abund(ncx,:)

        ENDIF
        !Initial calculations of diffusion frequency for each species bound to grain
        !and other parameters required for diffusion reactions
        DO  i=lbound(grainList,1),ubound(grainList,1)
            j=grainList(i)
            vdiff(i)=VDIFF_PREFACTOR*bindingEnergy(i)/mass(j)
            vdiff(i)=dsqrt(vdiff(i))
        END DO

        !h2 formation rate initially set
        h2form = 1.0d-17*dsqrt(initialTemp)
        ALLOCATE(mantle(points))
        DO l=1,points
            mantle(l)=sum(abund(grainList,l))
        END DO
        
        !DVODE SETTINGS
        ISTATE=1;;ITASK=1
        reltol=1e-4;MXSTEP=10000

        IF (.NOT. ALLOCATED(abstol)) THEN
            ALLOCATE(abstol(NEQ))
        END IF
        !OPTIONS = SET_OPTS(METHOD_FLAG=22, ABSERR_VECTOR=abstol, RELERR=reltol,USER_SUPPLIED_JACOBIAN=.FALSE.)
        

    END SUBROUTINE initializeChemistry

!Reads input reaction and species files as well as the final step of previous run if this is phase 2
    SUBROUTINE reader
        IMPLICIT NONE
        integer i,j,l,m
        double precision junktemp

        IF (ALLOCATED(outIndx)) DEALLOCATE(outIndx)
        IF (columnFlag) THEN
            nout = SIZE(outSpecies)
            ALLOCATE(outIndx(nout))
        END IF
        !assign array indices for important species to the integers used to store them.
        DO i=1,nspec
            IF (specname(i).eq.'H')   nh  = i
            IF (specname(i).eq.'H2')  nh2 = i
            IF (specname(i).eq.'C')   nc  = i
            IF (specname(i).eq.'C+')  ncx = i
            IF (specname(i).eq.'O')   no  = i
            IF (specname(i).eq.'N')   nn  = i
            IF (specname(i).eq.'S+')  ns  = i
            IF (specname(i).eq.'HE')  nhe = i
            IF (specname(i).eq.'CO')  nco = i
            IF (specname(i).eq.'MG')  nmg = i
            IF (specname(i).eq.'H2O') nh2o = i
            IF (specname(i).eq.'SI')  nsi = i
            IF (specname(i).eq.'SI+') nsix= i
            IF (specname(i).eq.'CL')  ncl = i
            IF (specname(i).eq.'CL+') nclx= i
            IF (specname(i).eq.'CH3OH') nch3oh= i
            IF (specname(i).eq.'#CO') ngrainco = i
            IF (specname(i).eq. 'P') np=i
            IF (specname(i).eq.'F') nf=i
            IF (columnFlag) THEN
                DO j=1,nout
                    IF (specname(i).eq.outSpecies(j)) outIndx(j)=i
                END DO
            END IF
        END DO

        !read start file if choosing to use abundances from previous run 
        !
        IF (readAbunds .eq. 1) THEN
            DO l=1,points
                READ(7,*)
                READ(7,7000) abund(nspec+1,l),junktemp,av(l)
                READ(7,*)
                READ(7,7010) h2form,fc,fo,&
                            &fmg,fhe,dstep
                READ(7,*)
                READ(7,7030) (abund(i,l),i=1,nspec)
                REWIND(7)
                (nspec+1,l)=dens(l)
            END DO
            7000 format(&
            &33x,1pe10.4,5x,/,&
            &33x,0pf8.2,2x,/,&
            &33x,0pf12.4,4x,/)
            7010 format(&
            &33x,1pe8.2,8x,/,&
            &11x,1pe7.1,4x,12x,1pe7.1,/&
            &12x,1pe7.1,13x,1pe7.1,&
            &13x,i3,/)
            7020  format(//)
            7030  format(4(18x,1pe10.3,:))     
        END IF
    END SUBROUTINE reader

!Writes physical variables and fractional abundances to output file, called every time step.
    SUBROUTINE output
        !write out cloud properties
        write(10,8020) timeInYears,dens(dstep),temp(dstep),av(dstep),radfield,zeta,h2form,fc,fo,&
                        &fmg,fhe,dstep
        !and a blank line
        write(10,8000)
        !and then all the abundances for this step
        write(10,8010) (specname(i),abund(i,dstep),i=1,nspec) 
        write(10,8000)
        !If this is the last time step of phase I, write a start file for phase II
        IF (readAbunds .eq. 0) THEN
           IF (switch .eq. 0 .and. timeInYears .ge. finalTime& 
               &.or. switch .eq. 1 .and.dens(dstep) .ge. finalDens) THEN
               write(7,8020) timeInYears,dens(dstep),temp(dstep),av(dstep),radfield,zeta,h2form,fc,fo,&
                       &fmg,fhe,dstep
               write(7,8000)
               write(7,8010) (specname(i),abund(i,dstep),i=1,nspec)
               write(7,8000)
           ENDIF
        ENDIF
        8000  format(/)
        8010  format(4(1x,a15,'=',1x,1pe10.3,:))
        8020 format(&
        &'age of cloud             time  = ',1pe10.3,' years',/,&
        &'total hydrogen density   dens  = ',1pe10.4,' cm-3',/,&
        &'cloud temperature        temp  = ',0pf8.2,' k',/,&
        &'visual extinction        av    = ',0pf12.4,' mags',/,&
        &'radiation field          rad   = ',0pf10.2,' (habing = 1)',/,&
        &'cosmic ray ioniz. rate   zeta  = ',0pf10.2,' (unit = 1.3e-17s-1)',/,&
        &'h2 formation rate coef.        = ',1pe8.2,' cm3 s-1',/,&
        &'c / htot = ',1pe7.1,4x,' o / htot = ',1pe7.1,/&
        &'mg / htot = ',1pe7.1,&
        &' he / htot = ',1pe7.1,&
        &' depth     = ',i3)

        !Every 'writestep' timesteps, write the chosen species out to separate file
        !choose species you're interested in by looking at parameters.f90
        IF (writeCounter==writeStep .and. columnFlag) THEN
            writeCounter=0
            write(11,8030) timeInYears,dens(dstep),temp(dstep),abund(outIndx,dstep)
            8030  format(1pe11.3,1x,1pe11.4,1x,0pf8.2,6(1x,1pe10.3))
        ELSE
            writeCounter=writeCounter+1
        END IF
    END SUBROUTINE output

    SUBROUTINE updateChemistry
    !Called every time/depth step and updates the abundances of all the species
        !allow option for dens to have been changed elsewhere.
        IF (collapse .ne. 1) abund(nspec+1,dstep)=dens(dstep)
        !y is at final value of previous depth iteration so set to initial values of this depth with abund
        !reset other variables for good measure        
        h2form = 1.0d-17*dsqrt(temp(dstep))
    
        !Sum of abundaces of all mantle species. mantleindx stores the indices of mantle species.
        mantle(dstep)=sum(abund(grainList,dstep))
        !evaluate co and h2 column densities for use in rate calculations
        !sum column densities of each point up to dstep. boxlength and dens are pulled out of the sum as common factors  
        IF (dstep.gt.1) THEN
            h2col=(sum(abund(nh2,:dstep-1))+0.5*abund(nh2,dstep))*dens(dstep)*(cloudSize/real(points))
            cocol=(sum(abund(nco,:dstep-1))+0.5*abund(nco,dstep))*dens(dstep)*(cloudSize/real(points))
        ELSE
            h2col=0.5*abund(nh2,dstep)*dens(dstep)*(cloudSize/real(points))
            cocol=0.5*abund(nco,dstep)*dens(dstep)*(cloudSize/real(points))
        ENDIF
        !call the actual ODE integrator
        CALL integrate
        !call evaporation to remove species from grains at certain temperatures

        CALL thermalEvaporation
        !1.d-30 stops numbers getting too small for fortran.
        WHERE(abund<1.0d-30) abund=1.0d-30
        dens(dstep)=abund(nspec+1,dstep)
    END SUBROUTINE updateChemistry

    SUBROUTINE integrate
    !This subroutine calls DVODE (3rd party ODE solver) until it can reach targetTime with acceptable errors (reltol/abstol)
        DO WHILE(currentTime .lt. targetTime)         
            !reset parameters for DVODE
            ITASK=1 !try to integrate to targetTime
            ISTATE=1 !pretend every step is the first
            reltol=1e-4 !relative tolerance effectively sets decimal place accuracy
            abstol=1.0d-14*abund(:,dstep) !absolute tolerances depend on value of abundance
            WHERE(abstol<1d-30) abstol=1d-30 ! to a minimum degree

            !get reaction rates for this iteration
            CALL calculateReactionRates
            !Call the integrator.
            OPTIONS = SET_OPTS(METHOD_FLAG=22, ABSERR_VECTOR=abstol, RELERR=reltol,USER_SUPPLIED_JACOBIAN=.FALSE.,MXSTEP=MXSTEP)
            CALL DVODE_F90(F,NEQ,abund(:,dstep),currentTime,targetTime,ITASK,ISTATE,OPTIONS)
            SELECT CASE(ISTATE)
                CASE(-1)
                    !More steps required for this problem
                    MXSTEP=MXSTEP*2    
                CASE(-2)
                    !Tolerances are too small for machine but succesful to current currentTime
                    abstol=abstol*10.0
                CASE(-3)
                    write(*,*) "DVODE found invalid inputs"
                    write(*,*) "abstol:"
                    write(*,*) abstol
                    STOP
                CASE(-4)
                    !Successful as far as currentTime but many errors.
                    !Make targetTime smaller and just go again
                    targetTime=currentTime+10.0/year
                CASE(-5)
                    targetTime=currentTime*1.01
            END SELECT
        END DO                   
    END SUBROUTINE integrate

    !This is where reacrates subroutine is hidden
    include 'rates.f90'

    SUBROUTINE  F (NEQ, T, Y, YDOT)
        INTEGER, PARAMETER :: WP = KIND(1.0D0)
        INTEGER NEQ
        REAL(WP) T
        REAL(WP), DIMENSION(NEQ) :: Y, YDOT
        INTENT(IN)  :: NEQ, T, Y
        INTENT(OUT) :: YDOT
        DOUBLE PRECISION :: D,loss,prod
        !Set D to the gas density for use in the ODEs
        D=y(NEQ)
        ydot=0.0
        !The ODEs created by MakeRates go here, they are essentially sums of terms that look like k(1,2)*y(1)*y(2)*dens. Each species ODE is made up
        !of the reactions between it and every other species it reacts with.
        INCLUDE 'odes.f90'

        !updated just in case temp changed
        h2form=1.0d-17*dsqrt(temp(dstep))

        !H2 formation should occur at both steps - however note that here there is no 
        !temperature dependence. y(nh) is hydrogen fractional abundance.
        ydot(nh)  = ydot(nh) - 2.0*( h2form*y(nh)*D - h2dis*y(nh2) )
        !                             h2 formation - h2-photodissociation
        ydot(nh2) = ydot(nh2) + h2form*y(nh)*D - h2dis*y(nh2)
        !                       h2 formation  - h2-photodissociation

        ! get density change from physics module to send to DLSODE
        IF (collapse .eq. 1) ydot(NEQ)=densdot(y(NEQ))
    END SUBROUTINE F

!integrate calls reacrates to get the reaction rates at every iteration. reacrates calls further functions.
!This file is already long so I've hidden those subroutines in rates.f90
 
    SUBROUTINE thermalEvaporation
    !Evaporation is based on Viti et al. 2004. A proportion of the frozen species is released into the gas phase
    !in specific events. These events are activated by flags (eg solidflag) which can be set in physics module.
    !The species evaporated are in lists, created by Makerates and based on groupings. see the viti 2004 paper.
    IF (mantle(dstep) .gt. 1d-30) THEN
        !Viti 04 evap

        IF (evap .eq. 1) THEN
            !Solid Evap
            IF (solidflag .eq. 1) THEN
                abund(gasGrainList,dstep)=abund(gasGrainList,dstep)+solidFractions*abund(grainList,dstep)
                abund(grainList,dstep)=(1.0-solidFractions)*abund(grainList,dstep)
                !Set flag to 2 to stop it being recalled
                solidflag=2
            ENDIF
            !monotonic evaporation at binding energy of species
            CALL bindingEnergyEvap
            !Volcanic evap
            IF (volcflag .eq. 1) THEN
                abund(gasGrainList,dstep)=abund(gasGrainList,dstep)+volcanicFractions*abund(grainList,dstep)
                abund(grainList,dstep)=(1.0-volcanicFractions)*abund(grainList,dstep)
                volcflag=2 !Set flag to 2 to stop it being recalled
            ENDIF

            !Co-desorption
            IF (coflag .eq. 1) THEN
                abund(gasGrainList,dstep)=abund(gasGrainList,dstep)+abund(grainList,dstep)
                abund(grainList,dstep)=1d-30
                coflag=2
            ENDIF

        ELSE IF (evap .eq. 2 .and. coflag .ne. 2) THEN
            !Alternative evap. Instaneous evaporation of all grain species
            abund(gasGrainList,dstep)=abund(gasGrainList,dstep)+abund(grainList,dstep)
            abund(grainList,dstep)=1d-30
            coflag = 2
        ENDIF
    ENDIF

    END SUBROUTINE thermalEvaporation

    SUBROUTINE bindingEnergyEvap
        !Subroutine to handle mono-evaporation. See viti 2004
        double precision en,newm,expdust,freq,kevap
        integer speci
        !mono evaporation at the binding energy of each species
        DO i=lbound(grainList,1),ubound(grainList,1)
            speci=grainList(i)
            en=bindingEnergy(i)*kbolt
            expdust=bindingEnergy(i)/temp(dstep)
            newm = mass(speci)*1.66053e-27
            freq = dsqrt((2*(SURFACE_SITE_DENSITY)*en)/((pi**2)*newm))
            kevap=freq*exp(-expdust)
            IF (kevap .ge. 0.99) THEN
                abund(gasGrainList(i),dstep)=abund(gasGrainList(i),dstep)+(monoFractions(i)*abund(speci,dstep))
                abund(speci,dstep)=(1.0-monoFractions(i))*abund(speci,dstep)
                !set to 1d50 so it can't happen again
                bindingEnergy(i)=1d50
            END IF 
        END DO
    END SUBROUTINE bindingEnergyEvap

    SUBROUTINE debugout
        open(79,file='output/debuglog',status='unknown')       !debug file.
        write(79,*) "Integrator failed, printing relevant debugging information"
        write(79,*) "dens",dens(dstep)
        write(79,*) "density in integration array",abund(nspec+1,dstep)
        write(79,*) "Av", av(dstep)
        write(79,*) "Mantle", mantle(dstep)
        write(79,*) "Temp", temp(dstep)
        DO i=1,nreac
            if (rate(i) .ge. huge(i)) write(79,*) "Rate(",i,") is potentially infinite"
        END DO
    END SUBROUTINE debugout
END MODULE chemistry