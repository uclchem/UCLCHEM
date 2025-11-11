MODULE IO
    USE constants
    USE DEFAULTPARAMETERS
    USE chemistry
    USE physicscore
    USE network
    USE heating
    
    ! CHARACTER (LEN=100) :: abundSaveFile, abundLoadFile, outputFile, columnFile
    LOGICAL :: columnOutput=.False.,fullOutput=.False.,rateOutput=.False.,fluxOutput=.False.,&
    &readAbunds=.False.,writeAbunds=.False.,heatingOutput=.False.
    CHARACTER (LEN=15),ALLOCATABLE :: outSpecies(:)
    INTEGER :: nout
    INTEGER, ALLOCATABLE :: outIndx(:)

    INTEGER, PARAMETER :: outputId=10,columnId=11,rateConstantId=12,ratesId=13,heatingId=14,abundLoadID=71,abundSaveID=72,outID=74,debugId=79,inputId=21
CONTAINS
    !Reads input reaction and species files as well as the final step of previous run if this is phase 2
    SUBROUTINE fileSetup
        IMPLICIT NONE
        ! TODO: improve the file setup to interact with in memory mode (not opening the file handles if not needed.)
        INQUIRE(UNIT=columnId, OPENED=columnOutput)
        IF (columnOutput) WRITE(columnId,333) specName(outIndx)
        333 FORMAT("Time,Density,gasTemp,dustTemp,av,radfield,zeta,",(999(A,:,',')))

        INQUIRE(UNIT=outputId, OPENED=fullOutput)
        IF (fullOutput) THEN
            WRITE(outputId,335) specName
        END IF
        335 FORMAT("Time,Density,gasTemp,dustTemp,Av,radfield,zeta,point,",(999(A,:,',')))
        
        INQUIRE(UNIT=rateConstantId, OPENED=rateOutput)
        INQUIRE(UNIT=ratesId, OPENED=fluxOutput)
        INQUIRE(UNIT=abundLoadID, OPENED=readAbunds)
        INQUIRE(UNIT=abundSaveID, OPENED=writeAbunds)
        INQUIRE(UNIT=heatingId, OPENED=heatingOutput)


        if (heatingOutput) then
            ! Write cooling/heating rates headers dynamically
            ! Format: Time, [4 cooling values], [NCOOL line cooling], [8 heating values], [1 chem heating]
            WRITE(heatingId,'(A)',advance='no') "Time,"
            WRITE(heatingId,'(A)',advance='no') "AtomicCooling,CollisionallyInducedEmission,ComptonCooling,ContinuumEmission,"
            WRITE(heatingId,'(A)',advance='no') "LineCooling_H,LineCooling_C+,LineCooling_O,LineCooling_C,LineCooling_CO,LineCooling_p-H2,LineCooling_o-H2,"
            WRITE(heatingId,'(A)') "Photoelectric,H2Formation,FUVPumping,Photodissociation,CIonization,CRHeating,TurbHeating,GasGrainColls,ChemicalHeating"
        END IF

    END SUBROUTINE fileSetup

    SUBROUTINE readInputAbunds
        !read start file if choosing to use abundances from previous run 
        !f2py integer, intent(aux) :: points
        IF (readAbunds) THEN
            DO l=1,points
                READ(abundLoadID,*) abund(:nspec,l)
                REWIND(abundLoadID)
            END DO
        END IF
    END SUBROUTINE readInputAbunds

    SUBROUTINE finalOutput
        !f2py integer, intent(aux) :: points
        IF (writeAbunds) THEN
            DO dstep=1,points
                ! WRITE(abundSaveID,*) fhe,fc,fo,fn,fs,fmg
                WRITE(abundSaveID,8010) abund(:nspec+2,dstep)
            8010  FORMAT((999(1pe15.5,:,',')))
            END DO
        END IF
    END SUBROUTINE finalOutput

    SUBROUTINE output(returnArray,writerates,successflag,physicsarray, chemicalabunarray, ratesarray, heatarray, dtime, timepoints)
        DOUBLE PRECISION, DIMENSION(:, :, :), OPTIONAL :: physicsarray
        DOUBLE PRECISION, DIMENSION(:, :, :), OPTIONAL :: chemicalabunarray
        DOUBLE PRECISION, DIMENSION(:, :, :), OPTIONAL :: ratesarray
        DOUBLE PRECISION, DIMENSION(:, :, :), OPTIONAL :: heatarray
        INTEGER, OPTIONAL :: dtime, timepoints
        INTEGER, intent(out) :: successflag
        INTEGER :: i  ! Loop variable for heating array assignment
        LOGICAL :: returnArray, writerates
        successflag = 0
        IF (returnArray) THEN
            ! Try to catch out of bounds errors before they create a segfault
            if (dtime .gt. timepoints+1) then
                write(*,*) "Ran out of timepoints in arrays, trying to stop gracefully"
                successflag=NOT_ENOUGH_TIMEPOINTS_ERROR
                return
            else 
                physicsarray(dtime, dstep, 1) = timeInYears
                physicsarray(dtime, dstep, 2) = density(dstep)
                physicsarray(dtime, dstep, 3) = gasTemp(dstep)
                physicsarray(dtime, dstep, 4) = dustTemp(dstep)
                physicsarray(dtime, dstep, 5) = av(dstep)
                physicsarray(dtime, dstep, 6) = radfield
                physicsarray(dtime, dstep, 7) = zeta
                physicsarray(dtime, dstep, 8) = dstep
                ! chemicalabunarray(dtime, dstep, :) = abund(:neq-1,dstep)
                chemicalabunarray(dtime, dstep, :) = abund(1:nspec,dstep)
            end if 
        ELSE IF (fullOutput .AND. .NOT. returnArray) THEN
            WRITE(outputId,8020) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep),&
                & av(dstep),radfield,zeta,dstep,abund(1:nspec,dstep)
            8020 FORMAT(1pe11.3,',',1pe11.4,',',0pf8.2,',',0pf8.2,',',1pe11.4,',',1pe11.4,&
            &','1pe11.4,',',I4,',',(999(1pe15.5,:,',')))
        END IF
        IF (writerates) THEN
            IF (returnArray) THEN
                ! If returnArray is true, we write the rates to the rates array, we compute the flux in Python.
                ratesarray(dtime, dstep, :) = rate(:nreac)
                ! Only populate the heating array if it is present, properly sized, AND heating is enabled
                IF (SIZE(heatarray, 1) .ge. timePoints .AND. heatingFlag .AND. ALLOCATED(lineCoolingArray)) THEN
                    heatarray(dtime, dstep, 1) = timeInYears
                    heatarray(dtime, dstep, 2:5) = coolingValues(:)
                    ! Write all NCOOL line cooling terms (currently 7: H, C+, O, C, CO, p-H2, o-H2)
                    heatarray(dtime, dstep, 6:(5+NCOOL)) = lineCoolingArray(median_line_index, :)

                    ! Heating terms (8 values)
                    heatarray(dtime, dstep, (6+NCOOL):(13+NCOOL)) = heatingValues(:)
                    ! Chemical heating (1 value)
                    heatarray(dtime, dstep, 14+NCOOL) = chemheating
                END IF
            ELSE 
                ! Else, we write the rate constants and rates to the file.
                IF (rateOutput) THEN
                    WRITE(rateConstantId,8021) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep),av(dstep),radfield,zeta,dstep,RATE
                    8021 FORMAT(1pe11.3,',',1pe11.4,',',0pf8.2,',',0pf8.2,',',1pe11.4,',',1pe11.4,','1pe11.4,',',I4,',',(9999(1pe15.5e3,:,',')))
                END IF
                if (fluxOutput) THEN
                    WRITE(ratesId,8022) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep),av(dstep),radfield,zeta,dstep,REACTIONRATE
                    8022 FORMAT(1pe11.3,',',1pe11.4,',',0pf8.2,',',0pf8.2,',',1pe11.4,',',1pe11.4,','1pe11.4,',',I4,',',(9999(1pe15.5e3,:,',')))
                END IF
                IF (heatingOutput) THEN
                    ! Write: time, cooling values (4), line cooling array (NCOOL), heating values (8), chem heating
                    WRITE(heatingId,8023) timeInYears, coolingValues(:), &
                        lineCoolingArray(median_line_index, :), &
                        heatingValues(:), chemheating
                    8023 FORMAT(1PE16.6E3,:,(999(',',1PE16.6E3)))
                END IF
            END IF
        END IF

        !Every 'writestep' timesteps, write the chosen species out to separate file
        !choose species you're interested in by looking at parameters.f90
        IF (.NOT. PRESENT(dtime)) THEN
            IF (writeCounter==writeStep .and. columnOutput) THEN
                writeCounter=1
                WRITE(columnId,8030) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep),&
                &av(dstep),radfield,zeta,abund(outIndx,dstep)
                8030  FORMAT(1pe11.3,',',1pe11.4,',',0pf8.2,',',0pf8.2,',',1pe11.4,',',1pe11.4,&
                &',',1pe11.4,',',(999(1pe15.5,:,',')))
            ELSE
                writeCounter=writeCounter+1
            END IF
        END IF
    END SUBROUTINE output

    SUBROUTINE closeFiles
        CLOSE(outputId)
        CLOSE(rateConstantId)
        CLOSE(columnId)
        CLOSE(abundSaveID)
        CLOSE(abundLoadID)
        CLOSE(heatingId)  ! heating rates file

    END SUBROUTINE closeFiles

    SUBROUTINE debugout
        OPEN(debugId,file='output/debuglog',status='unknown')       !debug file.
        WRITE(debugId,*) "Integrator failed, printing relevant debugging information"
        WRITE(debugId,*) "dens",density(dstep)
        WRITE(debugId,*) "gas temperature in integration array",abund(nspec+1,dstep)
        WRITE(debugId,*) "density in integration array",abund(nspec+2,dstep)
        WRITE(debugId,*) "Av", av(dstep)
        WRITE(debugId,*) "Temp", gasTemp(dstep)
        DO i=1,nreac
            if (rate(i) .ge. huge(i)) write(debugId,*) "Rate(",i,") is potentially infinite"
        END DO
    END SUBROUTINE debugout

    SUBROUTINE simpleDebug(message)
    !A very simply subroutine for debugging, just write a bunch of variables to screen
    !so we can check they're all as expected. 
    !Argument message is a string to print before the variables
        CHARACTER(LEN=*) :: message
        WRITE(*,*) message
        WRITE(*,*)"endAtFinalDensity",endAtFinalDensity
        WRITE(*,*)"freefall",freefall
        WRITE(*,*)"initialDens",initialDens
        WRITE(*,*)"finalDens",finalDens
        WRITE(*,*)"initialTemp",initialTemp
        WRITE(*,*)"finalTime",finalTime
        WRITE(*,*)"rout",rout
        WRITE(*,*)"baseAv",baseAv
        WRITE(*,*) "freezeFactor",freezeFactor
        WRITE(*,*) "abstol_factor",abstol_factor
        WRITE(*,*) "neq",neq
        WRITE(*,*) "density abund",abund(neq,1)
        WRITE(*,*) "density arr",density
        WRITE(*,*) "gasTemp",gasTemp
        WRITE(*,*) "coldens",coldens
        WRITE(*,*) "av",av
    END SUBROUTINE simpleDebug

END MODULE IO
