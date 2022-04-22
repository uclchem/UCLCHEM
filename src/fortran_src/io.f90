MODULE IO
    USE chemistry
    USE physicscore
    
    CHARACTER (LEN=100) :: abundSaveFile, abundLoadFile, outputFile, columnFile, outFile
    LOGICAL :: columnOutput=.False.,fullOutput=.False.,readAbunds=.False.,writeAbunds=.False.
    character(LEN=15),ALLOCATABLE :: outSpecies(:)
    INTEGER :: nout
    INTEGER, ALLOCATABLE :: outIndx(:)

    INTEGER, PARAMETER :: outputId=10,columnId=11,abundLoadID=71,abundSaveID=72,outID=74,debugId=79,inputId=21
CONTAINS
    !Reads input reaction and species files as well as the final step of previous run if this is phase 2
    SUBROUTINE fileSetup
        IMPLICIT NONE
        INQUIRE(UNIT=columnId, OPENED=columnOutput)
        IF (columnOutput) write(columnId,333) specName(outIndx)
        333 format("Time,Density,gasTemp,av,zeta,",(999(A,:,',')))

        INQUIRE(UNIT=outputId, OPENED=fullOutput)
        IF (fullOutput) THEN
            write(outputId,334) fc,fo,fn,fs
            write(outputId,*) "Radfield ", radfield
            write(outputId,335) specName
        END IF
        335 format("Time,Density,gasTemp,av,zeta,point,",(999(A,:,',')))
        334 format("Elemental abundances, C:",1pe15.5e3," O:",1pe15.5e3," N:",1pe15.5e3," S:",1pe15.5e3)

        INQUIRE(UNIT=abundLoadID, OPENED=readAbunds)
        INQUIRE(UNIT=abundSaveID, OPENED=writeAbunds)

        !read start file if choosing to use abundances from previous run 
        !
        IF (readAbunds) THEN
            DO l=1,points
                READ(abundLoadID,*) fhe,fc,fo,fn,fs,fmg
                READ(abundLoadID,*) abund(:nspec,l)
                REWIND(abundLoadID)
            END DO
        END IF
    END SUBROUTINE fileSetup

    !Writes physical variables and fractional abundances to output file, called every time step.
    SUBROUTINE output

        IF (fullOutput) THEN
            write(outputId,8020) timeInYears,density(dstep),gasTemp(dstep),av(dstep),zeta,dstep,abund(:neq-1,dstep)
            8020 format(1pe11.3,',',1pe11.4,',',0pf8.2,',',1pe11.4,',',1pe11.4,',',I4,',',(999(1pe15.5,:,',')))
        END IF
       
        !Every 'writestep' timesteps, write the chosen species out to separate file
        !choose species you're interested in by looking at parameters.f90
        IF (writeCounter==writeStep .and. columnOutput) THEN
            writeCounter=1
            write(columnId,8030) timeInYears,density(dstep),gasTemp(dstep),av(dstep),zeta,abund(outIndx,dstep)
            8030  format(1pe11.3,',',1pe11.4,',',0pf8.2,',',1pe11.4,',',1pe11.4,',',(999(1pe15.5,:,',')))
        ELSE
            writeCounter=writeCounter+1
        END IF
    END SUBROUTINE output

    SUBROUTINE finalOutput
        IF (writeAbunds) THEN
            DO dstep=1,points
                write(abundSaveID,*) fhe,fc,fo,fn,fs,fmg
                write(abundSaveID,8010) abund(:neq-1,dstep)
            8010  format((999(1pe15.5,:,',')))
            END DO
        END IF
    END SUBROUTINE finalOutput

    SUBROUTINE debugout
        open(debugId,file='output/debuglog',status='unknown')       !debug file.
        write(debugId,*) "Integrator failed, printing relevant debugging information"
        write(debugId,*) "dens",density(dstep)
        write(debugId,*) "density in integration array",abund(nspec+1,dstep)
        write(debugId,*) "Av", av(dstep)
        write(debugId,*) "Temp", gasTemp(dstep)
        DO i=1,nreac
            if (rate(i) .ge. huge(i)) write(debugId,*) "Rate(",i,") is potentially infinite"
        END DO
    END SUBROUTINE debugout

    SUBROUTINE simpleDebug(message)
    !A very simply subroutine for debugging, just write a bunch of variables to screen
    !so we can check they're all as expected. 
    !Argument message is a string to print before the variables
        CHARACTER(LEN=*) :: message
        write(*,*) message
        write(*,*)"endAtFinalDensity",endAtFinalDensity
        write(*,*)"freefall",freefall
        write(*,*)"initialDens",initialDens
        write(*,*)"finalDens",finalDens
        write(*,*)"initialTemp",initialTemp
        write(*,*)"finalTime",finalTime
        write(*,*)"rout",rout
        write(*,*)"baseAv",baseAv
        write(*,*) "freezeFactor",freezeFactor
        write(*,*) "abstol_factor",abstol_factor
        write(*,*) "neq",neq
        write(*,*) "density abund",abund(neq,1)
        write(*,*) "density arr",density
        write(*,*) "gasTemp",gasTemp
        write(*,*) "coldens",coldens
        write(*,*) "av",av
    END SUBROUTINE simpleDebug

END MODULE IO
