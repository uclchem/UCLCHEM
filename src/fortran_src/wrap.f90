!Python wrapper for uclchem, compiled with "make python"
!Subroutines become functions in the Python module
MODULE uclchemwrap
    USE physicscore
    USE chemistry
    USE io
    USE constants
    IMPLICIT NONE
CONTAINS
    SUBROUTINE cloud(dictionary, outSpeciesIn,abundance_out,successFlag)
        USE cloud_mod

        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(500)
        INTEGER :: successFlag
        !f2py intent(in) dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,successFlag
        successFlag=1

        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,updatePhysics,updateTargetTime,sublimation)

        IF ((ALLOCATED(outIndx)) .and. (successFlag .ge. 0)) THEN 
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF 
    END SUBROUTINE cloud

    SUBROUTINE collapse(collapseIn,collapseFileIn,writeOut,dictionary, outSpeciesIn,abundance_out,successFlag)
        USE collapse_mod

        CHARACTER(LEN=*) :: dictionary, outSpeciesIn,collapseFileIn
        DOUBLE PRECISION :: abundance_out(500)
        INTEGER :: successFlag,collapseIn
        LOGICAL :: writeOut
        !f2py intent(in) collapseIn,dictionary,outSpeciesIn,collapseFileIn,writeOut
        !f2py intent(out) abundance_out,successFlag
        successFlag=1
        collapse_mode=collapseIn
        writePhysics = writeOut
        collapseFile = collapseFileIn
        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,updatePhysics,updateTargetTime,sublimation)
        
        IF ((ALLOCATED(outIndx)) .and. (successFlag .ge. 0)) THEN 
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF 
    END SUBROUTINE collapse

    SUBROUTINE hot_core(temp_indx,max_temp,dictionary, outSpeciesIn,abundance_out,successFlag)
        USE hotcore

        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(500),max_temp
        INTEGER :: temp_indx,successFlag
        !f2py intent(in) temp_indx,max_temp,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,successFlag
        maxTemp=max_temp
        tempIndx=temp_indx
        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,updatePhysics,updateTargetTime,sublimation)
        IF ((ALLOCATED(outIndx)) .and. (successFlag .ge. 0)) THEN 
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF 
    END SUBROUTINE hot_core

    SUBROUTINE cshock(shock_vel,timestep_factor,minimum_temperature,dictionary, outSpeciesIn,&
        &abundance_out,dissipation_time,successFlag)
        USE cshock_mod

        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(500),shock_vel,timestep_factor
        DOUBLE PRECISION :: minimum_temperature,dissipation_time
        INTEGER :: successFlag
        !f2py intent(in) shock_vel,timestep_factor,minimum_temperature,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,dissipation_time,successFlag
        vs=shock_vel
        timestepFactor=timestep_factor
        minimumPostshockTemp=minimum_temperature
        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,updatePhysics,updateTargetTime,sublimation)

        IF (successFlag .ge. 0) THEN 
            IF (ALLOCATED(outIndx)) abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
            dissipation_time=dissipationTime
        END IF
    END SUBROUTINE cshock

    SUBROUTINE jshock(shock_vel,dictionary, outSpeciesIn,abundance_out,successFlag)
        USE jshock_mod

        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(500),shock_vel
        INTEGER :: successFlag
        !f2py intent(in) shock_vel,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,successFlag
        vs=shock_vel
        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,updatePhysics,updateTargetTime,sublimation)

        IF ((ALLOCATED(outIndx)) .and. (successFlag .ge. 0)) THEN 
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF 
    END SUBROUTINE jshock




    SUBROUTINE get_rates(dictionary,abundancesIn,speciesIndx,rateIndxs,&
        &speciesRates,successFlag,transfer,swap,bulk_layers)
        USE cloud_mod
        CHARACTER(LEN=*) :: dictionary
        DOUBLE PRECISION :: abundancesIn(500),speciesRates(500),transfer,swap,bulk_layers
        DOUBLE PRECISION :: ydot(NEQ)
        INTEGER :: rateIndxs(500),speciesIndx,successFlag
        INTEGER :: speci,bulk_version,surface_version
        !f2py intent(in) dictionary,abundancesIn,speciesIndx
        !f2py intent(out) :: speciesRates,successFlag,transfer,swap,bulk_layers

        INCLUDE 'defaultparameters.f90'

        CALL dictionary_parser(dictionary, "",successFlag)
        IF (successFlag .lt. 0) THEN
            WRITE(*,*) 'Error reading parameter dictionary'
            RETURN
        END IF

        CALL coreInitializePhysics(successFlag)
        CALL initializePhysics(successFlag)
        IF (successFlag .lt. 0) then
            WRITE(*,*) 'Error initializing physics'
            RETURN
        END IF

        CALL initializeChemistry(readAbunds)
        
        dstep=1
        successFlag=1
        abund(:nspec,dstep)=abundancesIn(:nspec)
        abund(neq,dstep)=initialDens
        currentTime=0.0
        timeInYears=0.0

        targetTime=1.0d-7
        CALL updateChemistry(successFlag)
        CALL F(NEQ,currentTime,abund(:,dstep),ydot)

        speciesRates=rate(rateIndxs)

        IF ((specname(speciesIndx)(1:1) .eq. "#") .or.&
        & (specname(speciesIndx)(1:1) .eq. "@")) THEN
            DO speci=1,nSpec
                IF (specname(speci) .eq. "@"//specname(speciesIndx)(2:)) bulk_version=speci
                IF (specname(speci) .eq. "#"//specname(speciesIndx)(2:)) surface_version=speci
            END DO
            IF (YDOT(nsurface) .lt. 0) THEN
                transfer=YDOT(nsurface)*surfaceCoverage*abund(bulk_version,1)/safeBulk
            ELSE
                transfer=YDOT(nsurface)*surfaceCoverage*abund(surface_version,1)
            END If
            swap=totalSwap
            bulk_layers=bulkLayersReciprocal
        ELSE
            swap=0.0
            transfer=0.0
            bulk_layers=0.0
        END IF

    END SUBROUTINE get_rates

    SUBROUTINE get_odes(dictionary,abundancesIn,ratesOut)
        USE cloud_mod
        CHARACTER(LEN=*) :: dictionary
        DOUBLE PRECISION :: abundancesIn(500),ratesOut(500)
        INTEGER :: successFlag
        !f2py intent(in) dictionary,abundancesIn
        !f2py intent(out) :: ratesOut
        INCLUDE 'defaultparameters.f90'
        CALL dictionary_parser(dictionary, "",successFlag)

        call coreInitializePhysics(successFlag)
        CALL initializePhysics(successFlag)
        IF (successFlag .lt. 0) then
            WRITE(*,*) 'Error initializing physics'
            RETURN
        END IF

        CALL initializeChemistry(readAbunds)
        
        dstep=1
        successFlag=1
        abund(:nspec,dstep)=abundancesIn(:nspec)
        abund(neq,dstep)=initialDens
        currentTime=0.0
        timeInYears=0.0

        targetTime=1.0d-7
        CALL updateChemistry(successFlag)
        CALL F(NEQ,currentTime,abund(:,dstep),ratesOut(:NEQ))
    END SUBROUTINE get_odes

    SUBROUTINE solveAbundances(dictionary,outSpeciesIn,successFlag,&
        &modelInitializePhysics,modelUpdatePhysics,updateTargetTime,&
        &sublimation)
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        EXTERNAL modelInitializePhysics,updateTargetTime,modelUpdatePhysics,sublimation

        INTEGER, INTENT(OUT) :: successFlag
        successFlag=1

        ! Set the boolean to True if you want to be able to return a python array
        ! call the dictionary_parser function in order to read the dictionary of parameters
        INCLUDE 'defaultparameters.f90'
        !Read input parameters from the dictionary
        CALL dictionary_parser(dictionary, outSpeciesIn,successFlag)
        IF (successFlag .lt. 0) THEN
            successFlag=PARAMETER_READ_ERROR
            WRITE(*,*) 'Error reading parameter dictionary'
            RETURN
        END IF
        
        dstep=1
        currentTime=0.0
        timeInYears=0.0

        CALL fileSetup

        !Initialize the physics, first do core physics
        !Then do model specific. This allows it to overwrite core
        call coreInitializePhysics(successFlag)
        CALL modelInitializePhysics(successFlag)

        IF (successFlag .lt. 0) then
            successFlag=PHYSICS_INIT_ERROR
            WRITE(*,*) 'Error initializing physics'
            RETURN
        END IF

        CALL initializeChemistry(readAbunds)
        CALL readInputAbunds !this won't do anything if no abundLoadFile was in input
        !CALL simpleDebug("Initialized")

    
        dstep=1

        call output

        !loop until the end condition of the model is reached
        DO WHILE (((endAtFinalDensity) .and. (density(1) < finalDens)) .or. &
            &((.not. endAtFinalDensity) .and. (timeInYears < finalTime)))
          
            currentTimeold=currentTime

            !Each physics module has a subroutine to set the target time from the current time
            CALL updateTargetTime
            !loop over parcels, counting from centre out to edge of cloud
            DO dstep=1,points
                !reset time if this isn't first depth point
                currentTime=currentTimeold
                !update chemistry from currentTime to targetTime
                CALL updateChemistry(successFlag)
                IF (successFlag .lt. 0) THEN
                    write(*,*) 'Error updating chemistry'
                    RETURN
                END IF

                !get time in years for output, currentTime is now equal to targetTime
                timeInYears= currentTime/SECONDS_PER_YEAR

                !Update physics so it's correct for new currentTime and start of next time step
                Call coreUpdatePhysics
                CALL modelUpdatePhysics
                !Sublimation checks if Sublimation should happen this time step and does it
                CALL sublimation(abund)

                !write this depth step now time, chemistry and physics are consistent
                CALL output
            END DO
        END DO
        CALL finalOutput
        CALL closeFiles
    END SUBROUTINE solveAbundances

    SUBROUTINE dictionary_parser(dictionary, outSpeciesIn,successFlag)
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        INTEGER, INTENT(OUT) :: successFlag
        INTEGER, ALLOCATABLE, DIMENSION(:) :: locations
        LOGICAL :: ChemicalDuplicateCheck
        INTEGER :: posStart, posEnd, whileInteger,inputindx
        CHARACTER(LEN=100) :: inputParameter, inputValue

        close(10)
        close(11)
        close(7)

        !always deallocate these so that if user didn't specify them,
        ! they don't remain from previous run
        IF (ALLOCATED(outIndx)) DEALLOCATE(outIndx)
        IF (ALLOCATED(outSpecies)) DEALLOCATE(outSpecies)

        !All reads use IOSTAT which will change successFlag from 0 if an error occurs
        !so set zero and check for non-zero each loop.
        successFlag=0

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

            SELECT CASE (inputParameter)
                CASE('alpha')
                    !To provide alphas, set keyword alpha in inputdictionary with a dictionary value
                    !that dictionary should be index:value pairs for the alpha array    
                    posStart=scan(dictionary,'{')
                    posEnd=scan(dictionary,'}')
                    CALL coefficientParser(dictionary(posStart+1:posEnd),alpha)
                CASE('beta')
                    !To provide alphas, set keyword alpha in inputdictionary with a dictionary value
                    !that dictionary should be index:value pairs for the alpha array    
                    posStart=scan(dictionary,'{')
                    posEnd=scan(dictionary,'}')
                    CALL coefficientParser(dictionary(posStart+1:posEnd),beta)
                CASE('gamma')
                    !To provide alphas, set keyword alpha in inputdictionary with a dictionary value
                    !that dictionary should be index:value pairs for the alpha array    
                    posStart=scan(dictionary,'{')
                    posEnd=scan(dictionary,'}')
                    CALL coefficientParser(dictionary(posStart+1:posEnd),gama)
                CASE('initialtemp')
                    READ(inputValue,*,iostat=successFlag) initialTemp
                CASE('initialdens')
                    READ(inputValue,*,iostat=successFlag) initialDens
                CASE('finaldens')
                    READ(inputValue,*,iostat=successFlag) finalDens
                CASE('currenttime')
                    READ(inputValue,*,iostat=successFlag) currentTime
                CASE('finaltime')
                    READ(inputValue,*,iostat=successFlag) finalTime
                CASE('radfield')
                    READ(inputValue,*,iostat=successFlag) radfield
                CASE('zeta')
                    READ(inputValue,*,iostat=successFlag) zeta
                CASE('freezefactor')
                    READ(inputValue,*,iostat=successFlag) freezeFactor
                CASE('rout')
                    READ(inputValue,*,iostat=successFlag) rout
                CASE('rin')
                    READ(inputValue,*,iostat=successFlag) rin
                CASE('baseav')
                    READ(inputValue,*,iostat=successFlag) baseAv
                CASE('points')
                    READ(inputValue,*,iostat=successFlag) points
                CASE('endatfinaldensity')
                    Read(inputValue,*,iostat=successFlag) endAtFinalDensity
                CASE('freefall')
                    READ(inputValue,*,iostat=successFlag) freefall
                CASE('freefallfactor')
                    READ(inputValue,*,iostat=successFlag) freefallFactor
                CASE('desorb')
                    READ(inputValue,*,iostat=successFlag) desorb
                CASE('h2desorb')
                    READ(inputValue,*,iostat=successFlag) h2desorb
                CASE('crdesorb')
                    READ(inputValue,*,iostat=successFlag) crdesorb
                CASE('uvdesorb')
                    READ(inputValue,*,iostat=successFlag) uvdesorb
                CASE('thermdesorb')
                    READ(inputValue,*,iostat=successFlag) uvdesorb
                CASE('instantsublimation')
                    READ(inputValue,*,iostat=successFlag) instantSublimation
                CASE('cosmicrayattenuation')
                    READ(inputValue,*,iostat=successFlag) cosmicRayAttenuation
                CASE('ionmodel')
                    READ(inputValue,*,iostat=successFlag) ionModel
                CASE('improvedh2crpdissociation')
                    READ(inputValue,*,iostat=successFlag) improvedH2CRPDissociation
                CASE('ion')
                    READ(inputValue,*,iostat=successFlag) ion
                CASE('fhe')
                    READ(inputValue,*,iostat=successFlag) fhe
                CASE('fc')
                    READ(inputValue,*,iostat=successFlag) fc
                CASE('fo')
                    READ(inputValue,*,iostat=successFlag) fo
                CASE('fn')
                    READ(inputValue,*,iostat=successFlag) fn
                CASE('fs')
                    READ(inputValue,*,iostat=successFlag) fs
                CASE('fmg')
                    READ(inputValue,*,iostat=successFlag) fmg
                CASE('fsi')
                    READ(inputValue,*,iostat=successFlag) fsi
                CASE('fcl')
                    READ(inputValue,*,iostat=successFlag) fcl
                CASE('fp')
                    READ(inputValue,*,iostat=successFlag) fp
                CASE('ff')
                    READ(inputValue,*,iostat=successFlag) ff
                CASE('outspecies')
                    READ(inputValue,*,iostat=successFlag) nout
                    ALLOCATE(outIndx(nout))
                    ALLOCATE(outSpecies(nout))
                    IF (outSpeciesIn .eq. "") THEN
                        write(*,*) "Outspecies parameter set but no outspecies string given"
                        write(*,*) "general(parameter_dict,outSpeciesIn) requires a delimited string of species names"
                        write(*,*) "if outSpecies or columnFlag is set in the parameter dictionary"
                        successFlag=-1
                        RETURN
                    ELSE
                        READ(outSpeciesIn,*, END=22) outSpecies
                        IF (outSpeciesIn .eq. "") THEN
22                              write(*,*) "mismatch between outSpeciesIn and number given in dictionary"
                            write(*,*) "Number:",nout
                            write(*,*) "Species list:",outSpeciesIn
                            successFlag=-1
                            RETURN
                        END IF
                    END IF
                    !assign array indices for important species to the integers used to store them.
                    DO i=1,nspec
                        DO j=1,nout
                            IF (specname(i).eq.outSpecies(j)) outIndx(j)=i
                        END DO
                    END DO
                CASE('writestep')
                    READ(inputValue,*,iostat=successFlag) writeStep
                CASE('ebmaxh2')
                    READ(inputValue,*,iostat=successFlag) ebmaxh2
                CASE('epsilon')
                    READ(inputValue,*,iostat=successFlag) epsilon
                CASE('uvcreff')
                    READ(inputValue,*,iostat=successFlag) uvcreff
                CASE('ebmaxcr')
                    READ(inputValue,*,iostat=successFlag) ebmaxcr
                CASE('phi')
                    READ(inputValue,*,iostat=successFlag) phi
                CASE('ebmaxuvcr')
                    READ(inputValue,*,iostat=successFlag) ebmaxuvcr
                CASE('uv_yield')
                    READ(inputValue,*,iostat=successFlag) uv_yield
                CASE('metallicity')
                    READ(inputValue,*,iostat=successFlag) metallicity
                CASE('omega')
                    READ(inputValue,*,iostat=successFlag) omega
                CASE('reltol')
                    READ(inputValue,*,iostat=successFlag) reltol
                CASE('abstol_factor')
                    READ(inputValue,*,iostat=successFlag) abstol_factor
                CASE('abstol_min')
                    READ(inputValue,*,iostat=successFlag) abstol_min
                ! CASE('jacobian')
                !     READ(inputValue,*) jacobian
                CASE('abundsavefile')
                    READ(inputValue,*,iostat=successFlag) abundSaveFile
                    abundSaveFile = TRIM(abundSaveFile)
                    open(abundSaveID,file=abundSaveFile,status="unknown")
                CASE('abundloadfile')
                    READ(inputValue,*,iostat=successFlag) abundLoadFile
                    abundLoadFile = TRIM(abundLoadFile)
                    open(abundLoadID,file=abundLoadFile,status='old')
                CASE('outputfile')
                    READ(inputValue,*,iostat=successFlag) outFile
                    outputFile = trim(outFile)
                    fullOutput=.True.
                    open(outputId,file=outputFile,status='unknown',iostat=successFlag)
                    IF (successFlag .lt. 0) THEN
                        write(*,*) "An error occured when opening the output file!"//&
                                        & NEW_LINE('A')//&
                                    &" The failed file was ",outputFile&
                                    &, NEW_LINE('A')//"A common error is that the directory doesn't exist"&
                                    &//NEW_LINE('A')//"************************"
                        successFlag=-1
                    END IF
                CASE('columnfile')
                    IF (trim(outSpeciesIn) .NE. '' ) THEN
                        columnOutput=.True.
                        READ(inputValue,*,iostat=successFlag) columnFile
                        columnFile = trim(columnFile)
                        open(columnId,file=columnFile,status='unknown')
                    ELSE
                        WRITE(*,*) "Error in output species. No species were given but a column file was given."
                        WRITE(*,*) "columnated output requires output species to be chosen."
                        successFlag=-1
                        RETURN
                    END IF

                CASE DEFAULT
                    WRITE(*,*) "Problem with given parameter: '", trim(inputParameter), "'."
                    WRITE(*,*) "This is either not supported yet, or invalid."
                    successFlag=-1
                    RETURN
            END SELECT
            dictionary = dictionary(posEnd:)
            IF (SCAN(dictionary,',') .eq. 0) whileInteger=1

            !check for failure
            IF (successFlag .ne. 0) THEN
                WRITE(*,*) "Error reading ",inputParameter
                write(*,*) "This is usually due to wrong type."
                successFlag=PARAMETER_READ_ERROR
                RETURN
            END IF 
        END DO

    END SUBROUTINE dictionary_parser

    SUBROUTINE coefficientParser(coeffDictString,coeffArray)
        CHARACTER(LEN=*) :: coeffDictString
        REAL(dp), INTENT(INOUT) :: coeffArray(*)
        INTEGER :: inputIndx,posStart,posEnd
        CHARACTER(LEN=100) :: inputValue
        LOGICAL :: continue_flag
        
        continue_flag=.True.
        DO WHILE (continue_flag)
            !substring containing integer key
            posStart=1
            posEnd=SCAN(coeffDictString,':')
            !read it into index integer
            READ(coeffDictString(posStart:posEnd-1),*) inputindx

            !substring including alpha value for the index.
            posStart=posEnd+1
            posEnd=SCAN(coeffDictString,',')
            !last value will have a } instead of , so grab index and tell loop to finish
            IF (posEnd .eq. 0) THEN
                posEnd=SCAN(coeffDictString,"}")
                continue_flag=.False.
            END IF

            !read that substring
            inputValue=coeffDictString(posStart:posEnd-1)
            READ(inputValue,*) coeffArray(inputIndx)
            !update string to remove this entry
            coeffDictString=coeffDictString(posEnd+1:)
        END DO
    END SUBROUTINE coefficientParser

END MODULE uclchemwrap