MODULE Wrap
    USE physics
    USE chemistry
    IMPLICIT NONE

    CHARACTER (LEN=100) :: abundFile, outputFile, columnFile, outFile
    REAL(dp) :: colDensIn
CONTAINS
    !Subroutine that just runs UCLCHEM as standard, writing to finalDens
    !Intended to be called from python script in place of running uclchem binary
    SUBROUTINE GENERAL(dictionary, outSpeciesIn,results_dir,model_no)
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn,results_dir,model_no
        !f2py intent(in) dictionary, outSpeciesIn,results_dir,model_no

        CALL dictionary_parser(dictionary,outSpeciesIn)
        if (heatWriteFlag) open(15,file=results_dir//"cooling/"//model_no//".dat",status="unknown")
        if (heatWriteFlag) open(16,file=results_dir//"heating/"//model_no//".dat",status="unknown")

        dstep=1
        currentTime=0.0
        timeInYears=0.0

        CALL initializePhysics

        colDens=colDensIn
        av= coldens*avFactor
        CALL initializeChemistry
        CALL output
        CALL integrateToTarget 
        dstep=1
        IF (heatWriteFlag) THEN
            tempDot=getTempDot(abund(NEQ-1,dstep),abund(NEQ,dstep),radfield*EXP(-UV_FAC*av(dstep))&
                   &,abund(:,dstep),h2dis,h2form,zeta,rate(nR_C_hv),&
                    &1.0/GAS_DUST_DENSITY_RATIO,grainRadius,abund(exoReactants1,dstep),&
                    &abund(exoReactants2,dstep),RATE(exoReacIdxs),exothermicities,heatWriteFlag,&
                    &dustTemp(dstep),turbVel)
            CLOSE(15)
            CLOSE(16)
        END IF
        CLOSE(10)
        CLOSE(11)
        CLOSE(7)
    END SUBROUTINE GENERAL


    SUBROUTINE run_model_for_abundances(dictionary, outSpeciesIn,abundance_out,outTemp)
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(12),outTemp
        !f2py intent(in) dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,outTemp

        ! Set the boolean to True if you want to be able to return a python array
        ! call the dictionary_parser function in order to read the dictionary of parameters
        CALL dictionary_parser(dictionary, outSpeciesIn)
        dstep=1
        currentTime=0.0
        timeInYears=0.0
        CALL initializePhysics

        colDens=colDensIn

        CALL initializeChemistry

        CALL integrateToTarget 
        outTemp=gasTemp(1)
        abundance_out=abund(outIndx,1)

        write(*,*) "SUCCESS",outTemp
    END SUBROUTINE run_model_for_abundances

    !Subroutine that runs UCLCHEM but does not write to file.
    !Instead, it returns the abundance array (including dens and temp) and dust temp to python script
    SUBROUTINE Sampler(dictionary, outSpeciesIn,abundances_in,abundances_out,dust_temp)
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        double precision :: abundances_in(34),abundances_out(36),dust_temp
        !f2py intent(in) dictionary, outSpeciesIn,abundances_in
        !f2py intent(out) abundances_out,dust_temp
        CALL dictionary_parser(dictionary,outSpeciesIn)
        heatWriteFlag=.False.

        dstep=1
        currentTime=0.0
        timeInYears=0.0
        write(*,*) initialDens
        CALL initializePhysics

        colDens=colDensIn

        CALL initializeChemistry
        abund(:nspec,1)=abundances_in
        CALL output
        CALL integrateToTarget
        write(*,*) abund
        abundances_out=abund(:,1)
        dust_temp=dustTemp(1)
    END SUBROUTINE Sampler


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Helper Functions
    !
    ! The subroutines/functions below are not intended to be called by your python script.
    ! Some will work but they are called by the funtions above to do tasks that most subroutines
    ! need such as parsing in the inputs.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !The core UCLCHEM loop. Abstracted out of the above subroutines. 
    SUBROUTINE integrateToTarget()
         !loop until the end condition of the model is reached
        DO WHILE ((switch .eq. 1 .and. density(1) < finalDens) .or. (switch .eq. 0 .and. timeInYears < finalTime))
            !store current time as starting point for each depth step
            IF (points .gt. 1) THEN
                currentTimeold=targetTime
                currentTime=currentTimeold
            END IF
            !Each physics module has a subroutine to set the target time from the current time
            CALL updateTargetTime

            !loop over parcels, counting from centre out to edge of cloud
            DO dstep=1,points
                !if (timeInYears .gt. 1000) heatingFlag=.True.
                !update chemistry
                CALL updateChemistry

                !set time to the final time of integrator rather than target
                targetTime=currentTime
                !reset target for next depth point
                if (points .gt. 1)currentTime=currentTimeold
                !get time in years for output
                timeInYears= currentTime/SECONDS_PER_YEAR
                write(*,*) timeInYears
                !update physics
                CALL updatePhysics
                
                !rely on thermal desorption reactions
                !CALL sublimation(abund)
                !write this depth step
                CALL output
                ! write(*,*) timeInYears
            END DO
        END DO  
    END SUBROUTINE integrateToTarget

    !A subroutine to parse the parameter dictionary sent from python. 
    SUBROUTINE dictionary_parser(dictionary, outSpeciesIn)
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        LOGICAL :: ChemicalDuplicateCheck
        INTEGER :: posStart, posEnd, whileInteger
        CHARACTER(LEN=100) :: inputParameter, inputValue

        INCLUDE "defaultparameters.f90"
        close(10)
        close(11)
        close(7)

        IF (scan(dictionary, 'columnFile') .EQ. 0) THEN
            columnFlag=.False.
        END IF

        whileInteger = 0

        posStart = scan(dictionary, '{')

        DO WHILE (whileInteger .NE. 1)
            posEnd = scan(dictionary, ':')
            inputParameter = dictionary(posStart+2:posEnd-2)
            dictionary = dictionary(posEnd:)
            posStart = scan(dictionary, ' ')
            IF (scan(dictionary, ',') .EQ. 0) THEN
                posEnd = scan(dictionary, '}')
                whileInteger = 1
            ELSE
                posEnd = scan(dictionary, ',')
            END IF
            inputValue = dictionary(posStart+1:posEnd-1)
            dictionary = dictionary(posEnd:)


            SELECT CASE (inputParameter)
                CASE('pahAbund')
                    READ(inputValue,*) pahAbund
                CASE('grainRadius')
                    READ(inputValue,*) grainRadius
                CASE('gasDustMassRatio')
                    READ(inputValue,*) gasDustMassRatio    
                CASE('heatWriteFlag')
                    READ(inputValue,*) heatWriteFlag
                CASE('metallicity') 
                    READ(inputValue,*) metallicity
                CASE('avFactor')
                    READ(inputValue,*) avFactor
                CASE('heatingFlag') 
                    READ(inputValue,*) heatingFlag
                CASE('h2col') 
                    READ(inputValue,*) h2col
                CASE('ccol') 
                    READ(inputValue,*) ccol
                CASE('coldens') 
                    READ(inputValue,*) colDensIn
                CASE('initialTemp')
                    READ(inputValue,*) initialTemp
                CASE('maxTemp')
                    READ(inputValue,*) maxTemp
                CASE('initialDens')
                    READ(inputValue,*) initialDens
                CASE('finalDens')
                    READ(inputValue,*) finalDens
                CASE('currentTime')
                    READ(inputValue,*) currentTime
                CASE('finalTime')
                    READ(inputValue,*) finalTime
                CASE('radfield')
                    READ(inputValue,*) radfield
                CASE('zeta')
                    READ(inputValue,*) zeta
                CASE('fr')
                    READ(inputValue,*) fr
                CASE('rout')
                    READ(inputValue,*) rout
                CASE('rin')
                    READ(inputValue,*) rin
                CASE('baseAv')
                    READ(inputValue,*) baseAv
                CASE('points')
                    READ(inputValue,*) points
                CASE('switch')
                    Read(inputValue,*) switch
                CASE('collapse')
                    READ(inputValue,*) collapse
                CASE('bc')
                    READ(inputValue,*) bc
                CASE('readAbunds')
                    READ(inputValue,*) readAbunds
                CASE('phase')
                    READ(inputValue,*) phase
                CASE('desorb')
                    READ(inputValue,*) desorb
                CASE('h2desorb')
                    READ(inputValue,*) h2desorb
                CASE('crdesorb')
                    READ(inputValue,*) crdesorb
                CASE('uvcr')
                    READ(inputValue,*) uvcr
                CASE('instantSublimation')
                    READ(inputValue,*) instantSublimation
                CASE('ion')
                    READ(inputValue,*) ion
                CASE('tempindx')
                    READ(inputValue,*) tempindx
                CASE('fh')
                    READ(inputValue,*) fh
                CASE('fhe')
                    READ(inputValue,*) fhe
                CASE('fc')
                    READ(inputValue,*) fc
                CASE('fo')
                    READ(inputValue,*) fo
                CASE('fn')
                    READ(inputValue,*) fn
                CASE('fs')
                    READ(inputValue,*) fs
                CASE('fmg')
                    READ(inputValue,*) fmg
                CASE('fsi')
                    READ(inputValue,*) fsi
                CASE('fcl')
                    READ(inputValue,*) fcl
                CASE('fp')
                    READ(inputValue,*) fp
                CASE('ff')
                    READ(inputValue,*) ff
                CASE('outSpecies')
                    IF (ALLOCATED(outIndx)) DEALLOCATE(outIndx)
                    IF (ALLOCATED(outSpecies)) DEALLOCATE(outSpecies)
                    READ(inputValue,*) nout
                    ALLOCATE(outIndx(nout))
                    ALLOCATE(outSpecies(nout))
                    IF (outSpeciesIn .eq. "") THEN
                        write(*,*) "Outspecies parameter set but no outspecies string given"
                        write(*,*) "general(parameter_dict,outSpeciesIn) requires a delimited string of species names"
                        write(*,*) "if outSpecies or columnFlag is set in the parameter dictionary"
                        STOP
                    ELSE
                        READ(outSpeciesIn,*, END=22) outSpecies
                        IF (outSpeciesIn .eq. "") THEN
    22                      write(*,*) "mismatch between outSpeciesIn and number given in dictionary"
                            write(*,*) "Number:",nout
                            write(*,*) "Species list:",outSpeciesIn
                            STOP
                        END IF
                    END IF  
                    !assign array indices for important species to the integers used to store them.
                    DO i=1,nspec
                        DO j=1,nout
                            IF (specname(i).eq.outSpecies(j)) outIndx(j)=i
                        END DO
                    END DO
                CASE('writeStep')
                    READ(inputValue,*) writeStep
                CASE('ebmaxh2')
                    READ(inputValue,*) ebmaxh2
                CASE('epsilon')
                    READ(inputValue,*) epsilon
                CASE('ebmaxcrf')
                    READ(inputValue,*) ebmaxcrf
                CASE('uvcreff')
                    READ(inputValue,*) uvcreff
                CASE('ebmaxcr')
                    READ(inputValue,*) ebmaxcr
                CASE('phi')
                    READ(inputValue,*) phi
                CASE('ebmaxuvcr')
                    READ(inputValue,*) ebmaxuvcr
                CASE('uv_yield')
                    READ(inputValue,*) uv_yield
                CASE('omega')
                    READ(inputValue,*) omega
                CASE('vs')
                    READ(inputValue,*) vs
                CASE('abundFile')
                    READ(inputValue,*) abundFile
                    abundFile = trim(abundFile)
                    open(7,file=abundFile,status='unknown')
                CASE('outputFile')
                    READ(inputValue,*) outFile
                    outputFile = trim(outFile)
                    open(10,file=outputFile,status='unknown')
                CASE('columnFile')
                    IF (trim(outSpeciesIn) .NE. '' ) THEN
                        READ(inputValue,*) columnFile
                        columnFile = trim(columnFile)
                        columnFlag=.True.
                        open(11,file=columnFile,status='unknown')
                    ELSEIF (trim(outSpeciesIn) .NE. '' ) THEN
                        WRITE(*,*) "Error in output species. No species were given but a column file was given."
                        WRITE(*,*) "columnated output requires output species to be chosen."
                        STOP
                    END IF

                CASE DEFAULT
                    WRITE(*,*) "Problem with given parameter: ", trim(inputParameter)
                    WRITE(*,*) "This is either not supported yet, or invalid"
            END SELECT
        END DO
    END SUBROUTINE dictionary_parser

END Module Wrap

! SUBROUTINE equilibrium(dictionary, outSpeciesIn,results_dir,model_no)
!     CHARACTER(LEN=*) :: dictionary, outSpeciesIn,results_dir,model_no
!     REAL(dp) :: colDensIn,newTemp,fixedCooling,fixedHeating
!     LOGICAL :: fully_converged,coolingFlag,heatingFixFlag
!     REAL(dp), parameter :: tempCriterion=0.1    
!     close(10)
!     close(11)
!     close(7)

!     open(15,file=results_dir//"cooling/"//model_no//".dat",status="unknown")
!     open(16,file=results_dir//"heating/"//model_no//".dat",status="unknown")


!     CALL initializePhysics

!     colDens=colDensIn

!     CALL initializeChemistry

!     !loop until the end condition of the model is reached
!     fully_converged=.False.
!     newTemp=0.0
!     loopCount=0
!     dstep=1

!     !send a 0 when we want real cooling so check if that's true
!     coolingFlag=(fixedCooling .lt. 1.0d-29)
!     heatingFixFlag=(fixedHeating .eq. 0)
!     write(*,*) coolingFlag,fixedCooling
!     DO WHILE (.NOT. fully_converged)
!         loopCount=loopCount+1
!         write(*,*) "loop",loopCount,"temp",gasTemp(dstep)
!         currentTime=0.0d0
!         targetTime=finalTime*SECONDS_PER_YEAR
!         heatingFlag=.False.

!         abund(NEQ-1,:)=gasTemp
!         abund(NEQ,:)=density 

!         CALL updatePhysics

!         CALL updateChemistry
!         CALL output

!         newTemp=getEquilibriumTemp(abund(NEQ-1,dstep),abund(NEQ,dstep),radfield*EXP(-UV_FAC*av(dstep)),&
!             &abund(:,dstep),h2dis,h2form,zeta,rate(nR_C_hv),&
!                 &1.0/GAS_DUST_DENSITY_RATIO,abund(exoReactants1,dstep),abund(exoReactants2,dstep),RATE(exoReacIdxs),&
!                 exothermicities,.True.,&
!                 &dustTemp(dstep),turbVel,fixedCooling,coolingFlag,fixedHeating,heatingFixFlag)

!         IF (ABS(newTemp-gasTemp(dstep)) .lt. tempCriterion) fully_converged=.True.

!         gasTemp(dstep)=newTemp
!         abund(NEQ-1,dstep)=gasTemp(dstep)

!         IF (loopCount .gt. 50) EXIT
!     END DO
!     CALL output
!     ! CALL writePopulations(results_dir//model_no//'.pop',model_no)
!     ! CALL writeOpacities(results_dir//model_no//'.tau',model_no)
!     CLOSE(10)
!     CLOSE(11)
!     CLOSE(7)
!     CLOSE(15)
!     CLOSE(16)
! END SUBROUTINE equilibrium