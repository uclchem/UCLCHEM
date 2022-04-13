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
USE photoreactions
USE surfacereactions
USE constants
IMPLICIT NONE
   !These integers store the array index of important species and reactions, x is for ions    
    INTEGER :: njunk
    !loop counters    
    INTEGER :: i,j,l,writeStep,writeCounter=0,loopCounter
    INTEGER, PARAMETER :: maxLoops=10

    !Flags to control desorption processes
    INTEGER :: h2desorb,crdesorb,uvdesorb,desorb,thermdesorb


    !Array to store reaction rates
    REAL(dp) :: rate(nreac)
    
    !IO Options
    character(LEN=15),ALLOCATABLE :: outSpecies(:)
    LOGICAL :: columnOutput=.False.,fullOutput=.False.,readAbunds=.False.,writeAbunds=.False.
    INTEGER :: nout
    INTEGER, ALLOCATABLE :: outIndx(:)


    !DLSODE variables    
    INTEGER :: ITASK,ISTATE,NEQ,MXSTEP
    REAL(dp) :: reltol,abstol_factor,abstol_min
    REAL(dp), ALLOCATABLE :: abstol(:)
    TYPE(VODE_OPTS) :: OPTIONS
    !initial fractional elemental abudances and arrays to store abundances
    REAL(dp) :: fh,fd,fhe,fc,fo,fn,fs,fmg,fsi,fcl,fp,ff,fli,fna,fpah,f15n,f13c,f18O,metallicity
    REAL(dp) :: h2col,cocol,ccol,h2colToCell,cocolToCell,ccolToCell
    REAL(dp),ALLOCATABLE :: abund(:,:)
    
    !Variables controlling chemistry
    LOGICAL :: PARAMETERIZE_H2FORM=.True.
    REAL(dp) :: radfield,fr,omega,grainArea,cion,h2dis,lastTemp=0.0
    REAL(dp) :: ebmaxh2,epsilon,ebmaxcrf,ebmaxcr,phi,ebmaxuvcr,uv_yield,uvcreff
    

    REAL(dp) :: turbVel=1.0


CONTAINS
    SUBROUTINE initializeChemistry
    ! Sets variables at the start of every run.
    ! Since python module persists, it's not enough to set initial
    ! values in module definitions above. Reset here.
        NEQ=nspec+1
        IF (ALLOCATED(abund)) DEALLOCATE(abund,vdiff)
        ALLOCATE(abund(NEQ,points),vdiff(SIZE(iceList)))
        CALL fileSetup
        !Set abundances to initial elemental if not reading them in.
        IF (.NOT. readAbunds) THEN
            !ensure abund is initially zero
            abund= 0.

            !Start by filling all metallicity scaling elements
            !neutral atoms  
            abund(no,:) = fo  
            abund(nn,:) = fn               
            abund(nmg,:) = fmg
            abund(np,:) = fp
            abund(nf,:) = ff
            abund(nna,:) = fna
            abund(nli,:) = fli
            abund(npah,:) = fpah
            !default to ions
            abund(nsx,:) = fs
            abund(nsix,:) = fsi                
            abund(nclx,:) = fcl 
            !Decide how much carbon is initiall ionized using parameters.f90
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

            !isotopes
            abund(n18o,:) = f18o  
            abund(n15n,:) = f15n           
            abund(n13c,:) = f13c    

            abund(nelec,:)=abund(ncx,:)+abund(nsix,:)+abund(nsx,:)+abund(nclx,:)

            abund=abund*metallicity

            !Total H nuclei is always 1 so put fh into H and whatever is left over in H2
            abund(nh,:) = fh
            abund(nh2,:) = 0.5*(1.0e0-fh) 
            abund(nd,:)=fd

            abund(nhe,:) = fhe  
            abund(nspec+1,:)=density  

        ENDIF
        !Initial calculations of diffusion frequency for each species bound to grain
        !and other parameters required for diffusion reactions
        DO  i=lbound(iceList,1),ubound(iceList,1)
            j=iceList(i)
            vdiff(i)=VDIFF_PREFACTOR*bindingEnergy(i)/mass(j)
            vdiff(i)=dsqrt(vdiff(i))
        END DO
        
        !DVODE SETTINGS
        ISTATE=1
        ITASK=1

        IF (.NOT. ALLOCATED(abstol)) THEN
            ALLOCATE(abstol(NEQ))
        END IF
        !OPTIONS = SET_OPTS(METHOD_FLAG=22, ABSERR_VECTOR=abstol, RELERR=reltol,USER_SUPPLIED_JACOBIAN=.FALSE.)
        
        !Set rates to zero to ensure they don't hold previous values or random ones if we don't set them in calculateReactionRates
        rate=0.0
        lastTemp=99.0d99!to save time, we always store last time step's temperature here and don't recalcuate some rates unless there's a change
    END SUBROUTINE initializeChemistry

!Reads input reaction and species files as well as the final step of previous run if this is phase 2
    SUBROUTINE fileSetup
        IMPLICIT NONE
        INQUIRE(UNIT=11, OPENED=columnOutput)
        IF (columnOutput) write(11,333) specName(outIndx)
        333 format("Time,Density,gasTemp,av,zeta,",(999(A,:,',')))

        INQUIRE(UNIT=10, OPENED=fullOutput)
        IF (fullOutput) THEN
            write(10,334) fc,fo,fn,fs
            write(10,*) "Radfield ", radfield, " Zeta ",zeta
            write(10,335) specName
        END IF
        335 format("Time,Density,gasTemp,av,point,zeta,",(999(A,:,',')))
        334 format("Elemental abundances, C:",1pe15.5e3," O:",1pe15.5e3," N:",1pe15.5e3," S:",1pe15.5e3)

        INQUIRE(UNIT=71, OPENED=readAbunds)
        INQUIRE(UNIT=72, OPENED=writeAbunds)

        !read start file if choosing to use abundances from previous run 
        !
        IF (readAbunds) THEN
            DO l=1,points
                READ(71,*) fhe,fc,fo,fn,fs,fmg
                READ(71,*) abund(:nspec,l)
                REWIND(71)
                abund(nspec+1,l)=density(l)
            END DO
        END IF
    END SUBROUTINE fileSetup

!Writes physical variables and fractional abundances to output file, called every time step.
    SUBROUTINE output

        IF (fullOutput) THEN
            write(10,8020) timeInYears,density(dstep),gasTemp(dstep),av(dstep),dstep,zeta,abund(:neq-1,dstep)
            8020 format(1pe11.3,',',1pe11.4,',',0pf8.2,',',1pe11.4,',',I4,',',(999(1pe15.5,:,',')))
        END IF

        !If this is the last time step of phase I, write a start file for phase II
        IF (writeAbunds) THEN
           IF (switch .eq. 0 .and. timeInYears .ge. finalTime& 
               &.or. switch .eq. 1 .and.density(dstep) .ge. finalDens) THEN
               write(72,*) fhe,fc,fo,fn,fs,fmg
               write(72,8010) abund(:neq-1,dstep)
           ENDIF
        ENDIF
        8010  format((999(1pe15.5,:,',')))
        

        !Every 'writestep' timesteps, write the chosen species out to separate file
        !choose species you're interested in by looking at parameters.f90
        IF (writeCounter==writeStep .and. columnOutput) THEN
            writeCounter=1
            write(11,8030) timeInYears,density(dstep),gasTemp(dstep),av(dstep),zeta,abund(outIndx,dstep)
            8030  format(1pe11.3,',',1pe11.4,',',0pf8.2,',',1pe11.4,',',(999(1pe15.5,:,',')))
        ELSE
            writeCounter=writeCounter+1
        END IF
    END SUBROUTINE output

    SUBROUTINE updateChemistry(successFlag)
    !Called every time/depth step and updates the abundances of all the species
        INTEGER, INTENT(OUT) :: successFlag
        loopCounter=0
        successFlag=1
        DO WHILE((currentTime .lt. targetTime) .and. (loopCounter .lt. maxLoops))   
            !allow option for dens to have been ch anged elsewhere.
            IF (collapse .ne. 1) abund(nspec+1,dstep)=density(dstep)

            !First sum the total column density over all points further towards edge of cloud
            IF (dstep.gt.1) THEN
                h2ColToCell=(sum(abund(nh2,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
                coColToCell=(sum(abund(nco,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
                cColToCell=(sum(abund(nc,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
            ELSE
                h2ColToCell=0.0
                coColToCell=0.0
                cColToCell=0.0
            ENDIF
            !then add half the column density of the current point to get average in this "cell"
            h2Col=h2ColToCell+0.5*abund(nh2,dstep)*density(dstep)*(cloudSize/real(points))
            coCol=coColToCell+0.5*abund(nco,dstep)*density(dstep)*(cloudSize/real(points))
            cCol=cColToCell+0.5*abund(nc,dstep)*density(dstep)*(cloudSize/real(points))

            !recalculate coefficients for ice processes
            safeMantle=MAX(1d-30,abund(nSurface,dstep))
            safeBulk=MAX(1d-30,abund(nBulk,dstep))
            bulkLayersReciprocal=MIN(1.0,NUM_SITES_PER_GRAIN/(GAS_DUST_DENSITY_RATIO*safeBulk))
            surfaceCoverage=bulkGainFromMantleBuildUp()

            
            CALL calculateReactionRates

            CALL integrateODESystem(successFlag)
            IF (successFlag .lt. 0) THEN
                write(*,*) "Integration failed, exiting"
                RETURN
            END IF

            !1.d-30 stops numbers getting too small for fortran.
            WHERE(abund<1.0d-30) abund=1.0d-30
            density(dstep)=abund(NEQ,dstep)
            loopCounter=loopCounter+1
        END DO
        IF (loopCounter .eq. maxLoops) successFlag=-1
    END SUBROUTINE updateChemistry

    SUBROUTINE integrateODESystem(successFlag)
        INTEGER, INTENT(OUT) :: successFlag
        successFlag=0
    !This subroutine calls DVODE (3rd party ODE solver) until it can reach targetTime with acceptable errors (reltol/abstol)
        !reset parameters for DVODE
        ITASK=1 !try to integrate to targetTime
        ISTATE=1 !pretend every step is the first
        abstol=abstol_factor*abund(:,dstep) !absolute tolerances depend on value of abundance
        WHERE(abstol<abstol_min) abstol=abstol_min ! to a minimum degree

        !Call the integrator.
        OPTIONS = SET_OPTS(METHOD_FLAG=22, ABSERR_VECTOR=abstol, RELERR=reltol,USER_SUPPLIED_JACOBIAN=.False.,MXSTEP=MXSTEP)
        CALL DVODE_F90(F,NEQ,abund(:,dstep),currentTime,targetTime,ITASK,ISTATE,OPTIONS)

        SELECT CASE(ISTATE)
            CASE(-1)
                write(*,*) "ISTATE -1: MAXSTEPS will be increased"
                !More steps required for this problem
                MXSTEP=MXSTEP*2   
                targetTime=currentTime*1.01 
            CASE(-2)
                write(*,*) "ISTATE -2: Tolerances too small"
                !Tolerances are too small for machine but succesful to current currentTime
                abstol_factor=abstol_factor*10.0
            CASE(-3)
                write(*,*) "DVODE found invalid inputs"
                write(*,*) "abstol:"
                write(*,*) abstol
                successFlag=-1
                RETURN
            CASE(-4)
                !Successful as far as currentTime but many errors.
                !Make targetTime smaller and just go again
                write(*,*) "ISTATE -4 - shortening step"
                targetTime=currentTime*1.01
            CASE(-5)
                timeInYears=currentTime/SECONDS_PER_YEAR
                write(*,*) "ISTATE -5 - shortening step at time", timeInYears,"years"
                !WHERE(abund<1.0d-30) abund=1.0d-30
                targetTime=currentTime*1.01
            CASE default
                MXSTEP=10000    
        END SELECT
    END SUBROUTINE integrateODESystem

    !This is where reacrates subroutine is hidden
    include 'rates.f90'

    SUBROUTINE F (NEQUATIONS, T, Y, YDOT)
        INTEGER, PARAMETER :: WP = KIND(1.0D0)
        INTEGER NEQUATIONS,looper
        REAL(WP) T
        REAL(WP), DIMENSION(NEQUATIONS) :: Y, YDOT
        INTENT(IN)  :: NEQUATIONS, T, Y
        INTENT(OUT) :: YDOT
        REAL(dp) :: D,loss,prod
        !Set D to the gas density for use in the ODEs
        D=y(NEQ)
        ydot=0.0
    
        !changing abundances of H2 and CO can causes oscillation since their rates depend on their abundances
        !recalculating rates as abundances are updated prevents that.
        !thus these are the only rates calculated each time the ODE system is called.
        cocol=coColToCell+0.5*Y(nco)*D*(cloudSize/real(points))
        h2col=h2ColToCell+0.5*Y(nh2)*D*(cloudSize/real(points))
        rate(nR_H2_hv)=H2PhotoDissRate(h2Col,radField,av(dstep),turbVel) !H2 photodissociation
        rate(nR_CO_hv)=COPhotoDissRate(h2Col,coCol,radField,av(dstep)) !CO photodissociation

        !recalculate coefficients for ice processes
        safeMantle=MAX(1d-30,Y(nSurface))
        safeBulk=MAX(1d-30,Y(nBulk)-SUM(Y(refractoryList)))
        bulkLayersReciprocal=MIN(1.0,NUM_SITES_PER_GRAIN/(GAS_DUST_DENSITY_RATIO*safeBulk))
        surfaceCoverage=bulkGainFromMantleBuildUp()

        !The ODEs created by MakeRates go here, they are essentially sums of terms that look like k(1,2)*y(1)*y(2)*dens. Each species ODE is made up
        !of the reactions between it and every other species it reacts with.
        INCLUDE 'odes.f90'
        ! get density change from physics module to send to DLSODE

        IF (collapse .eq. 1) ydot(NEQUATIONS)=densdot(y(NEQUATIONS))
    END SUBROUTINE F

    ! SUBROUTINE JAC(NEQ, T, Y, ML, MU, J, NROWPD)
    !     INTEGER NEQ,ML,MU,NROWPD
    !     DOUBLE PRECISION T, Y(NEQ), J(NROWPD,NEQ)
    !     REAL(DP) :: D
    !     INTENT(IN)  :: NEQ, T, Y,ML,MU,NROWPD
    !     INTENT(INOUT) :: J
    !     D=y(NEQ)

    !     J=0.0d0
    !     INCLUDE 'jacobian.f90'
    !     J(nh,nh2)=J(nh,nh2)+2.0*h2dis
    !     J(nh2,nh2)=J(nh,nh2)-h2dis
    ! END SUBROUTINE JAC



    SUBROUTINE debugout
        open(79,file='output/debuglog',status='unknown')       !debug file.
        write(79,*) "Integrator failed, printing relevant debugging information"
        write(79,*) "dens",density(dstep)
        write(79,*) "density in integration array",abund(nspec+1,dstep)
        write(79,*) "Av", av(dstep)
        write(79,*) "Temp", gasTemp(dstep)
        DO i=1,nreac
            if (rate(i) .ge. huge(i)) write(79,*) "Rate(",i,") is potentially infinite"
        END DO
    END SUBROUTINE debugout
END MODULE chemistry