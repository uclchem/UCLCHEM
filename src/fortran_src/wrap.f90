!Python wrapper for uclchem, compiled with "make python"
!Subroutines become functions in the Python module
MODULE uclchemwrap
    USE physicscore
    USE chemistry
    USE io
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

    SUBROUTINE cshock(shock_vel,dictionary, outSpeciesIn,abundance_out,dissipation_time,successFlag)
        USE cshock_mod

        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(500),shock_vel,dissipation_time
        INTEGER :: successFlag
        !f2py intent(in) shock_vel,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,dissipation_time,successFlag
        vs=shock_vel
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




    SUBROUTINE get_rates(dictionary,abundancesIn,rateIndxs,speciesRates,successFlag)
        USE cloud_mod
        CHARACTER(LEN=*) :: dictionary
        DOUBLE PRECISION :: abundancesIn(500),speciesRates(500)
        INTEGER :: rateIndxs(500),successFlag
        !f2py intent(in) dictionary,abundancesIn
        !f2py intent(out) :: speciesRates,successFlag

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
        abund(:nspec,1)=abundancesIn(:nspec)

        !allow option for dens to have been changed elsewhere.
        IF (.not. freefall) abund(nspec+1,dstep)=density(dstep)

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
        speciesRates=rate(rateIndxs)

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

        CALL updateTargetTime
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
            WRITE(*,*) 'Error reading parameter dictionary'
            RETURN
        END IF
        
        dstep=1
        currentTime=0.0
        timeInYears=0.0


        !Initialize the physics, first do core physics
        !Then do model specific. This allows it to overwrite core
        call coreInitializePhysics(successFlag)
        CALL modelInitializePhysics(successFlag)

        IF (successFlag .lt. 0) then
            WRITE(*,*) 'Error initializing physics'
            RETURN
        END IF
        CALL initializeChemistry(readAbunds)
        CALL fileSetup

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
        !close outputs to attempt to force flush
        close(10)
        close(11)
        close(71)
        close(72)
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
        successFlag=1

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
                    CALL alpha_parser(dictionary(posStart+1:posEnd))
                CASE('initialtemp')
                    READ(inputValue,*,err=666) initialTemp
                CASE('initialdens')
                    READ(inputValue,*,err=666) initialDens
                CASE('finaldens')
                    READ(inputValue,*,err=666) finalDens
                CASE('currenttime')
                    READ(inputValue,*,err=666) currentTime
                CASE('finaltime')
                    READ(inputValue,*,err=666) finalTime
                CASE('radfield')
                    READ(inputValue,*,err=666) radfield
                CASE('zeta')
                    READ(inputValue,*,err=666) zeta
                CASE('freezefactor')
                    READ(inputValue,*,err=666) freezeFactor
                CASE('rout')
                    READ(inputValue,*,err=666) rout
                CASE('rin')
                    READ(inputValue,*,err=666) rin
                CASE('baseav')
                    READ(inputValue,*,err=666) baseAv
                CASE('points')
                    READ(inputValue,*,err=666) points
                CASE('endatfinaldensity')
                    Read(inputValue,*,err=666) endAtFinalDensity
                CASE('freefall')
                    READ(inputValue,*,err=666) freefall
                CASE('freefallfactor')
                    READ(inputValue,*,err=666) freefallFactor
                CASE('desorb')
                    READ(inputValue,*,err=666) desorb
                CASE('h2desorb')
                    READ(inputValue,*,err=666) h2desorb
                CASE('crdesorb')
                    READ(inputValue,*,ERR=666) crdesorb
                CASE('uvdesorb')
                    READ(inputValue,*,ERR=666) uvdesorb
                CASE('thermdesorb')
                    READ(inputValue,*,ERR=666) uvdesorb
                CASE('instantsublimation')
                    READ(inputValue,*,ERR=666) instantSublimation
                CASE('cosmicrayattenuation')
                    READ(inputValue,*,ERR=666) cosmicRayAttenuation
                CASE('ionmodel')
                    READ(inputValue,*,ERR=666) ionModel
                CASE('improvedh2crpdissociation')
                    READ(inputValue,*,ERR=666) improvedH2CRPDissociation
                CASE('ion')
                    READ(inputValue,*,ERR=666) ion
                CASE('fhe')
                    READ(inputValue,*,ERR=666) fhe
                CASE('fc')
                    READ(inputValue,*,ERR=666) fc
                CASE('fo')
                    READ(inputValue,*,ERR=666) fo
                CASE('fn')
                    READ(inputValue,*,ERR=666) fn
                CASE('fs')
                    READ(inputValue,*,ERR=666) fs
                CASE('fmg')
                    READ(inputValue,*,ERR=666) fmg
                CASE('fsi')
                    READ(inputValue,*,ERR=666) fsi
                CASE('fcl')
                    READ(inputValue,*,ERR=666) fcl
                CASE('fp')
                    READ(inputValue,*,ERR=666) fp
                CASE('ff')
                    READ(inputValue,*,ERR=666) ff
                CASE('outspecies')
                    IF (ALLOCATED(outIndx)) DEALLOCATE(outIndx)
                    IF (ALLOCATED(outSpecies)) DEALLOCATE(outSpecies)
                    READ(inputValue,*,ERR=666) nout
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
                    READ(inputValue,*,ERR=666) writeStep
                CASE('ebmaxh2')
                    READ(inputValue,*,ERR=666) ebmaxh2
                CASE('epsilon')
                    READ(inputValue,*,ERR=666) epsilon
                CASE('uvcreff')
                    READ(inputValue,*,ERR=666) uvcreff
                CASE('ebmaxcr')
                    READ(inputValue,*,ERR=666) ebmaxcr
                CASE('phi')
                    READ(inputValue,*,ERR=666) phi
                CASE('ebmaxuvcr')
                    READ(inputValue,*,ERR=666) ebmaxuvcr
                CASE('uv_yield')
                    READ(inputValue,*,ERR=666) uv_yield
                CASE('metallicity')
                    READ(inputValue,*,ERR=666) metallicity
                CASE('omega')
                    READ(inputValue,*,ERR=666) omega
                CASE('reltol')
                    READ(inputValue,*,ERR=666) reltol
                CASE('abstol_factor')
                    READ(inputValue,*,ERR=666) abstol_factor
                CASE('abstol_min')
                    READ(inputValue,*,ERR=666) abstol_min
                ! CASE('jacobian')
                !     READ(inputValue,*) jacobian
                CASE('abundsavefile')
                    writeAbunds=.True.
                    READ(inputValue,*,ERR=666) abundSaveFile
                    abundSaveFile = TRIM(abundSaveFile)
                    open(72,file=abundSaveFile,status="unknown")
                CASE('abundloadfile')
                    READ(inputValue,*,ERR=666) abundLoadFile
                    abundLoadFile = TRIM(abundLoadFile)
                    readAbunds=.True.
                    open(71,file=abundLoadFile,status='old')
                CASE('outputfile')
                    READ(inputValue,*,ERR=666) outFile
                    outputFile = trim(outFile)
                    fullOutput=.True.
                    open(10,file=outputFile,status='unknown')
                CASE('columnfile')
                    IF (trim(outSpeciesIn) .NE. '' ) THEN
                        columnOutput=.True.
                        READ(inputValue,*,ERR=666) columnFile
                        columnFile = trim(columnFile)
                        open(11,file=columnFile,status='unknown')
                    ELSE
                        WRITE(*,*) "Error in output species. No species were given but a column file was given."
                        WRITE(*,*) "columnated output requires output species to be chosen."
                        successFlag=-1
                        RETURN
                    END IF

                CASE DEFAULT
                    WRITE(*,*) "Problem with given parameter: '", trim(inputParameter), "'."
                    WRITE(*,*) "This is either not supported yet, or invalid."
            END SELECT
            dictionary = dictionary(posEnd:)
            IF (SCAN(dictionary,',') .eq. 0) whileInteger=1
        END DO
    IF (successFlag .lt. 0) THEN
        666 successFlag=-1
        write(*,*) "Error reading",inputParameter
    END IF 
    END SUBROUTINE dictionary_parser

    SUBROUTINE alpha_parser(alpha_string)
        INTEGER :: inputIndx,posStart,posEnd
        CHARACTER(LEN=100) :: inputValue
        CHARACTER(LEN=*) :: alpha_string
        LOGICAL :: continue_flag
        
        continue_flag=.True.
        DO WHILE (continue_flag)
            !substring containing integer key
            posStart=1
            posEnd=SCAN(alpha_string,':')
            !read it into index integer
            READ(alpha_string(posStart:posEnd-1),*) inputindx

            !substring including alpha value for the index.
            posStart=posEnd+1
            posEnd=SCAN(alpha_string,',')
            !last value will have a } instead of , so grab index and tell loop to finish
            IF (posEnd .eq. 0) THEN
                posEnd=SCAN(alpha_string,"}")
                continue_flag=.False.
            END IF

            !read that substring
            inputValue=alpha_string(posStart:posEnd-1)
            READ(inputValue,*) alpha(inputIndx)
            !update string to remove this entry
            alpha_string=alpha_string(posEnd+1:)
        END DO
    END SUBROUTINE alpha_parser

END MODULE uclchemwrap