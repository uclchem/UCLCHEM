!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Chemistry module of UCL_CHEM.                                                               !
! Contains all the core machinery of the code, not really intended to be altered in standard  !
! use. Use a (custom) physics module to alter temp/density behaviour etc.                     !
!                                                                                             !
! chemistry module contains rates.f90, a series of subroutines to calculate all reaction rates!
! when updateChemistry is called from main, these rates are calculated, the ODEs are solved   !
! from currentTime to targetTime to get abundances at targetTime and then all abundances are  !
! written to the fullOutput file.                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE chemistry
USE physics
USE dvode_f90_m
USE network
USE heating
USE photoreactions
USE constants
IMPLICIT NONE
   !These integers store the array index of important species and reactions, x is for ions    
    integer :: nrco,njunk,evapevents,ngrainco,readAbunds
    !loop counters    
    integer :: i,j,l,writeStep,writeCounter=0

    !Array to store reaction rates
    REAL(dp) :: rate(nreac)
    
    !Option column output
    character(LEN=15),allocatable :: outSpecies(:)
    logical :: columnFlag,heatingFlag,fullOutput,heatWriteFlag
    integer :: nout
    integer, allocatable :: outIndx(:)


    !DLSODE variables    
    integer :: ITASK,ISTATE,NEQ,MXSTEP
    REAL(dp) :: reltol
    REAL(dp), allocatable :: abstol(:)
    TYPE(VODE_OPTS) :: OPTIONS

    !initial fractional elemental abudances and arrays to store abundances
    REAL(dp) :: fh,fhe,fc,fo,fn,fs,fmg,fsi,fcl,fp,ff,h2col,cocol,ccol,junk1,junk2,metallicity
    REAL(dp),allocatable :: abund(:,:),mantle(:)
    
    !Variables controlling chemistry
    REAL(dp) :: radfield,zeta,fr,omega,grainArea,cion,h2form,h2dis,tempDot,oldTemp=0.0
    REAL(dp) :: ebmaxh2,epsilon,ebmaxcrf,ebmaxcr,phi,ebmaxuvcr,uv_yield,uvcreff
    REAL(dp), allocatable ::vdiff(:)

    !Turbulent velocity of gas in cm/s for heating functions
    REAL(dp)  :: turbVel=1.0d5 !1 km/s in cm/s


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Grain surface parameters
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(dp) :: gasDustMassRatio=100.0,grainRadius=1.d-5, GRAIN_DENSITY = 3.0 ! Mass density of a dust grain
    REAL(dp)  :: THERMAL_VEL= SQRT(8.0d0*K_BOLTZ/(PI*AMU)) !Thermal velocity without the factor of SQRT(T/m) where m is moelcular mass in amu

    !reciprocal of fractional abundance of dust grains (we only divide by number density so better to store reciprocal)
    REAL(dp) :: GAS_DUST_DENSITY_RATIO
    !Grain area per h nuclei, calculated from average radius.
    REAL(dp) :: GRAIN_AREA_PER_H
    !Below are values for grain surface reactions
    LOGICAL, PARAMETER :: DIFFUSE_REACT_COMPETITION=.True., GRAINS_HAVE_ICE=.True.
    REAL(dp), PARAMETER :: CHEMICAL_BARRIER_THICKNESS = 1.40d-8  !gre Parameter used to compute the probability for a surface reaction with 
    !! activation energy to occur through quantum tunneling (Hasegawa et al. Eq 6 (1992).)
    REAL(dp), PARAMETER :: SURFACE_SITE_DENSITY = 1.5d15 ! site density on one grain [cm-2]
    REAL(dp) :: VDIFF_PREFACTOR=2.0*K_BOLTZ*SURFACE_SITE_DENSITY/PI/PI/AMU
    REAL(dp) :: NUM_SITES_PER_GRAIN 
CONTAINS
!This gets called immediately by main so put anything here that you want to happen before the time loop begins, reader is necessary.
    SUBROUTINE initializeChemistry
        !calculate basic dust grain properties from the two free parameters (grain radius and dust to gas mass ratio)
        GAS_DUST_DENSITY_RATIO = (4.0*PI*(grainRadius**3)*GRAIN_DENSITY*gasDustMassRatio)/(3.0 * AMU)
        GRAIN_AREA_PER_H=4.0*PI*grainRadius*grainRadius/GAS_DUST_DENSITY_RATIO
        NUM_SITES_PER_GRAIN = grainRadius*grainRadius*SURFACE_SITE_DENSITY*4.0*PI

        NEQ=nspec+2 !one ODE per species +1 each for temperature and density
        IF (ALLOCATED(abund)) DEALLOCATE(abund,vdiff,mantle)
        ALLOCATE(abund(NEQ,points),vdiff(nspec))
        CALL reader
        !if this is the first step of the first phase, set initial abundances
        !otherwise reader will fix it
        IF (readAbunds.eq.0) THEN
            !ensure abund is initially zero
            abund= 0.0D0

            !then assign metals
            abund(no,:) = fo  
            abund(nn,:) = fn               
            abund(nsx,:) = fs
            abund(nmg,:) = fmg
            abund(nsix,:) = fsi                
            abund(nclx,:) = fcl 
            !abund(np,:) = fp
            !abund(nf,:) = ff
            !abund(nfe,:) = ffe
            !abund(nna,:) = fna
            
            !Decide how much carbon is initially ionized using parameters.f90
            SELECT CASE (ion)
                CASE(0)
                    abund(nc,:)=fc
                    abund(ncx,:)=0.0
                CASE(1)
                    abund(nc,:)=fc*0.5
                    abund(ncx,:)=fc*0.5
                CASE(2)
                    abund(nc,:)=0.0
                    abund(ncx,:)=fc
            END SELECT
            abund(nspec,:)=abund(ncx,:)!+abund(nsix,:)+abund(nsx,:)+abund(nclx,:)
            !and adjust for metallicity
            abund=abund*metallicity
            !As default, have half in molecular hydrogen and half in atomic hydrogen
            abund(nh,:) = fh 
            abund(nh2,:) = (0.5*(1.0e0-fh))
            abund(nhe,:) = fhe                       
            abund(NEQ-1,:)=gasTemp
            abund(NEQ,:)=density      


        ENDIF
        !Initial calculations of diffusion frequency for each species bound to grain
        !and other parameters required for diffusion reactions
        DO  i=lbound(grainList,1),ubound(grainList,1)
            j=grainList(i)
            vdiff(i)=VDIFF_PREFACTOR*bindingEnergy(i)/mass(j)
            vdiff(i)=dsqrt(vdiff(i))
        END DO
        !h2 formation rate initially set
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
        CALL initializeHeating(initialTemp,initialDens,abund(:,1),colDens(dstep),cloudSize,heatWriteFlag)
        IF (columnFlag) write(11,333) specName(outIndx)
        333 format("Time,Density,gasTemp,dustTemp,av,",(999(A,:,',')))
        

         INQUIRE( UNIT=10, OPENED=fullOutput )
         if (fullOutput) write(10,334) specName
         334 format("Time,Density,gasTemp,dustTemp,av,radfield,zeta,",(999(A,:,',')))

    END SUBROUTINE initializeChemistry

!Reads input reaction and species files as well as the final step of previous run if this is phase 2
    SUBROUTINE reader
        IMPLICIT NONE
        integer i,l
        REAL(dp) junktemp

        !read start file if choosing to use abundances from previous run 
        !
        IF (readAbunds .eq. 1) THEN
            DO l=1,points
                READ(7,*)
                READ(7,7000) abund(NEQ,l),junktemp,av(l)
                READ(7,*)
                READ(7,7010) h2form,fc,fo,&
                            &fmg,fhe,dstep
                READ(7,*)
                READ(7,7030) (abund(i,l),i=1,nspec)
                REWIND(7)
                abund(neq,l)=density(l)
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
            7030  format(4(18x,1pe10.3,:))     
        END IF
    END SUBROUTINE reader

!Writes physical variables and fractional abundances to output file, called every time step.
    SUBROUTINE output

         INQUIRE( UNIT=10, OPENED=fullOutput )
         IF (fullOutput) THEN
            IF ((.not. REGULAR_STEPS) .or. (MOD(timeInYears,1000.0).lt.0.001)) THEN
                write(10,8020) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep),av(dstep),radfield,zeta, abund(:neq-2,dstep)
            END IF
        END IF
        8020 format(1pe11.3,',',1pe11.4,',',0pf8.2,',',0pf8.2,',',1pe11.4,',',1pe11.4,',',0pf8.2,',',(999(1pe15.5,:,',')))

        !Every 'writestep' timesteps, write the chosen species out to separate file
        !choose species you're interested in by looking at parameters.f90
        IF ((writeCounter==writeStep .or. timeInYears .lt. 1000.0) .and. columnFlag) THEN
            writeCounter=1
            write(11,8030) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep),av(dstep),abund(outIndx,dstep)
            8030  format(1pe11.3e3,',',1pe11.4e3,',',0pf8.2,',',0pf8.2,',',1pe11.4e3,',',(999(1pe15.5e3,:,',')))
        ELSE
            writeCounter=writeCounter+1
        END IF
    END SUBROUTINE output

    SUBROUTINE updateChemistry
    !Called every time/depth step and updates the abundances of all the species
        !gasTemp(dstep)=getEquilibriumTemp(gasTemp(dstep),abund(NEQ,dstep),radfield,abund(:,dstep),h2dis,zeta,rate(nR_C_hv)&
        !    &,1.0/GAS_DUST_DENSITY_RATIO,abund(exoReactants1,dstep),abund(exoReactants2,dstep),RATE(exoReacIdxs),exothermicities)
        abund(NEQ-1,dstep)=gasTemp(dstep)
        !allow option for dens to have been changed elsewhere.
        IF (collapse .ne. 1) abund(NEQ,dstep)=density(dstep)
    
        !Sum of abundaces of all mantle species. mantleindx stores the indices of mantle species.
        mantle(dstep)=sum(abund(grainList,dstep))
        !evaluate co and h2 column densities for use in rate calculations
        !sum column densities of each point up to dstep. boxlength and dens are pulled out of the sum as common factors  
        IF (dstep.gt.1) THEN
            !h2col=(sum(abund(nh2,:dstep-1)*density(:dstep-1))+0.5*abund(nh2,dstep)*density(dstep))*(cloudSize/real(points))
            cocol=(sum(abund(nco,:dstep-1)*density(:dstep-1))+0.5*abund(nco,dstep)*density(dstep))*(cloudSize/real(points))
        ELSE
            !h2col=abund(nh2,dstep)*density(dstep)*(cloudSize/real(points))
            cocol=abund(nco,dstep)*density(dstep)*(cloudSize/real(points))
        ENDIF

        !call the actual ODE integrator
        CALL integrate


        !1.d-30 stops numbers getting too small for fortran.
        WHERE(abund<1.0d-30) abund=1.0d-30

        density(dstep)=abund(NEQ,dstep)
        gasTemp(dstep)=abund(NEQ-1,dstep)
    END SUBROUTINE updateChemistry

    SUBROUTINE integrate
    !This subroutine calls DVODE (3rd party ODE solver) until it can reach targetTime with acceptable errors (reltol/abstol)
        DO WHILE(currentTime .lt. targetTime)         
            !reset parameters for DVODE
            ITASK=1 !try to integrate to targetTime
            ISTATE=1 !pretend every step is the first
            reltol=1e-5 !relative tolerance effectively sets decimal place accuracy
            abstol=1.0d-14*abund(:,dstep) !absolute tolerances depend on value of abundance
            abstol(neq-1)=1.0d-6
            WHERE(abstol<1d-20) abstol=1d-20 ! to a minimum degree


            !get reaction rates for this iteration - do it here since we may reloop a few times
            !First get dust temperature - requires local and surface field in Habing
            dustTemp(dstep)=calculateDustTemp(radfield*EXP(-UV_FAC*av(dstep)),radfield)
            !Then we can use that to get H2 formation rate
            !h2form= 3.0d-17*(SQRT(gastemp(dstep)/100.0))*(100.0/gasDustMassRatio)
            h2form= h2FormRate(gasTemp(dstep),dustTemp(dstep))!1.0d-17*dsqrt(gasTemp(dstep))

            !Finally calculate all other rates
            CALL calculateReactionRates
            !Tempdot is only recalculated after change in temperature of 1 degree so we force a calculation here
            IF (heatingFlag) tempDot=getTempDot(abund(NEQ-1,dstep),abund(NEQ,dstep),radfield*EXP(-UV_FAC*av(dstep)),abund(:,dstep),h2dis,h2form,zeta,rate(nR_C_hv),&
                &1.0/GAS_DUST_DENSITY_RATIO,grainRadius,metallicity,abund(exoReactants1,dstep),abund(exoReactants2,dstep),RATE(exoReacIdxs),exothermicities,heatWriteFlag,&
                &dustTemp(dstep),turbVel)

            !Call the integrator.
            OPTIONS = SET_OPTS(METHOD_FLAG=22, ABSERR_VECTOR=abstol, RELERR=reltol,USER_SUPPLIED_JACOBIAN=.FALSE.,MXSTEP=MXSTEP)
            CALL DVODE_F90(F,NEQ,abund(:,dstep),currentTime,targetTime,ITASK,ISTATE,OPTIONS)
            !write(*,*) ISTATE

            SELECT CASE(ISTATE)
                CASE(-1)
                    !More steps required for this problem
                    MXSTEP=MXSTEP*2    
                CASE(-2)
                    !Tolerances are too small for machine but succesful to current currentTime
                    abstol=abstol*10.0
                    write(*,*) "ISTATE -2, abstol too small"
                CASE(-3)
                    write(*,*) "DVODE found invalid inputs"
                    write(*,*) "abstol:"
                    write(*,*) abstol
                    writeCounter=writeStep
                    call output
                    STOP
                CASE(-4)
                    !Successful as far as currentTime but many errors.
                    !Make targetTime smaller and just go again
                    targetTime=currentTime*1.001
                    ISTATE=1
                CASE(-5)
                    write(*,*) "ISTATE -5"
                    write(*,*) gasTemp(dstep),currentTime/SECONDS_PER_YEAR
                    ISTATE=1
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
        REAL(dp) :: D,loss,prod


        ! write(*,*) "F start ",h2form,h2dis
        ! write(*,*) y(neq-1),ydot(nh),y(nh2),y(nh),D
        ! write(*,*) "***"

        !Set D to the gas density for use in the ODEs
        D=y(NEQ)
        ydot=0.0
        !The ODEs created by MakeRates go here, they are essentially sums of terms that look like k(1,2)*y(1)*y(2)*dens. Each species ODE is made up
        !of the reactions between it and every other species it reacts with.
        INCLUDE 'odes.f90'

        !H2 formation should occur at both steps - however note that here there is no 
        !temperature dependence. y(nh) is hydrogen fractional abundance.
       
        !                       h2 formation  - h2-photodissociation
        !get temperature change from heating module


        IF (heatingFlag) THEN
            IF (ABS(y(NEQ-1)-oldTemp).gt.0.1) THEN
                gasTemp(dstep)=y(NEQ-1)
                IF (gasTemp(dstep) .lt. 10) gasTemp(dstep)=10.0
                IF (gasTemp(dstep) .gt. 1.0d4) gasTemp(dstep)=1.0d6
                !CALL calculateReactionRates
                h2form= h2FormRate(gasTemp(dstep),dustTemp(dstep))
                tempDot=getTempDot(Y(NEQ-1),Y(NEQ),radfield*EXP(-UV_FAC*av(dstep)),Y,h2dis,h2form,zeta,rate(nR_C_hv),1.0/GAS_DUST_DENSITY_RATIO&
                    &,grainRadius,metallicity,y(exoReactants1),y(exoReactants2),RATE(exoReacIdxs),exothermicities,heatWriteFlag,dustTemp(dstep),turbVel)
                oldTemp=y(NEQ-1)
            END IF
            ydot(NEQ-1)=tempDot

        ELSE
            ydot(NEQ-1)=0.0
        END IF

        ydot(nh)  = ydot(nh) - 2.0*( h2form*y(nh)*D - h2dis*y(nh2) )
        !                             h2 formation - h2-photodissociation
        ydot(nh2) = ydot(nh2) + h2form*y(nh)*D - h2dis*y(nh2)
        ! write(*,*) "F see Tdot=",ydot(NEQ-1)
        ! get density change from physics module to send to DLSODE
        ydot(NEQ)=densdot(y(NEQ))
    END SUBROUTINE F


    SUBROUTINE debugout
        open(79,file='output/debuglog',status='unknown')       !debug file.
        write(79,*) "Integrator failed, printing relevant debugging information"
        write(79,*) "dens",density(dstep)
        write(79,*) "density in integration array",abund(NEQ,dstep)
        write(79,*) "Av", av(dstep)
        write(79,*) "Mantle", mantle(dstep)
        write(79,*) "Temp", gasTemp(dstep)
        DO i=1,nreac
            if (rate(i) .ge. huge(i)) write(79,*) "Rate(",i,") is potentially infinite"
        END DO
    END SUBROUTINE debugout
END MODULE chemistry