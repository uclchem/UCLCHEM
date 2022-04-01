MODULE IO
    USE chemistry
    USE physicscore
    
    CHARACTER (LEN=100) :: abundSaveFile, abundLoadFile, outputFile, columnFile, outFile
    LOGICAL :: columnOutput=.False.,fullOutput=.False.,readAbunds=.False.,writeAbunds=.False.
    character(LEN=15),ALLOCATABLE :: outSpecies(:)
    INTEGER :: nout
    INTEGER, ALLOCATABLE :: outIndx(:)
CONTAINS
    !Reads input reaction and species files as well as the final step of previous run if this is phase 2
    SUBROUTINE fileSetup
        IMPLICIT NONE
        INQUIRE(UNIT=11, OPENED=columnOutput)
        IF (columnOutput) write(11,333) specName(outIndx)
        333 format("Time,Density,gasTemp,av,",(999(A,:,',')))

        INQUIRE(UNIT=10, OPENED=fullOutput)
        IF (fullOutput) THEN
            write(10,334) fc,fo,fn,fs
            write(10,*) "Radfield ", radfield, " Zeta ",zeta
            write(10,335) specName
        END IF
        335 format("Time,Density,gasTemp,av,point,",(999(A,:,',')))
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
            write(10,8020) timeInYears,density(dstep),gasTemp(dstep),av(dstep),dstep,abund(:neq-1,dstep)
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
            write(11,8030) timeInYears,density(dstep),gasTemp(dstep),av(dstep),abund(outIndx,dstep)
            8030  format(1pe11.3,',',1pe11.4,',',0pf8.2,',',1pe11.4,',',(999(1pe15.5,:,',')))
        ELSE
            writeCounter=writeCounter+1
        END IF
    END SUBROUTINE output

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

END MODULE IO
