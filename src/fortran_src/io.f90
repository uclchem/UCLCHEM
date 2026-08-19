module IO
    ! allow(use-all)
    use chemistry
    use constants, only: dp
    use coolant_module, only: NCOOLANTS, coolants, populationsId
    use DEFAULTPARAMETERS
    use f2py_constants, only: coolantNames
    use heating, only: coolingLabels, NCOOLING, coolingValues, heatingValues, lineCoolingArray, &
        chemheating, median_line_index, NHEATING, se_coolant_iterations, se_coolant_max_rel_change, &
        heatingLabels
    use network, only: nSpec, nReac
    use physicscore, only: timeInYears, density, gasTemp, dustTemp, av, &
        radfield, zeta, dstep, parcel_radius, av_internal, radfield_internal

    implicit none

    public

    logical :: columnOutput=.false.,fullOutput=.false.,rateConstantsOutput=.false.,fluxOutput=.false.,&
    &readAbunds=.false.,writeAbunds=.false.,heatingOutput=.false.
    character (LEN=15),allocatable :: outSpecies(:)
    integer :: nout
    integer, allocatable :: outIndx(:)

    integer, parameter :: outputId=10,columnId=11,rateConstantId=12,ratesId=13,heatingId=14,abundLoadID=71, &
        & abundSaveID=72,outID=74,debugId=79,inputId=21
contains
    !Reads input reaction and species files as well as the final step of previous run if this is phase 2
    subroutine fileSetup
        inquire(unit=columnId, OPENED=columnOutput)
        if (columnOutput) write(columnId,333) specName(outIndx)
        333 format("Time,Density,gasTemp,dustTemp,av,radfield,zeta,",(999(A,:,",")))

        inquire(unit=outputId, OPENED=fullOutput)
        if (fullOutput) then
            write(outputId,335) specName
        end if
        335 format("Time,Density,gasTemp,dustTemp,Av,radfield,zeta,point,",(999(A,:,",")))

        inquire(unit=rateConstantId, OPENED=rateConstantsOutput)
        inquire(unit=ratesId, OPENED=fluxOutput)
        inquire(unit=abundLoadID, OPENED=readAbunds)
        inquire(unit=abundSaveID, OPENED=writeAbunds)
        inquire(unit=heatingId, OPENED=heatingOutput)


        if (heatingOutput) then
            ! Write cooling/heating rates headers dynamically
            ! Format: Time, [NCOOLING cooling values], [NCOOLANTS line cooling], [NHEATING heating values], [1 chem heating]
            write(heatingId,"(A)",advance="no") "Time,"

            ! Write cooling mechanism labels
            do i = 1, NCOOLING
                write(heatingId,"(A,A)",advance="no") TRIM(coolingLabels(i)), ","
            end do

            ! Write line cooling labels with species names
            do i = 1, NCOOLANTS
                write(heatingId,"(A,A,A)",advance="no") "LineCooling_", TRIM(coolantNames(i)), ","
            end do

            ! Write heating mechanism labels
            do i = 1, NHEATING
                write(heatingId,"(A,A)",advance="no") TRIM(heatingLabels(i)), ","
            end do

            ! Write chemical heating label (last column, no comma)
            write(heatingId,"(A)") "ChemicalHeating"
        end if

    end subroutine fileSetup

    subroutine readInputAbunds
        !read start file if choosing to use abundances from previous run
        !f2py integer, intent(aux) :: points
        if (readAbunds) then
            if (enable_radiative_transfer .AND. points>1) then
                print*,"we are reading the initial abundances"
                do l=points,1,-1
                    read(abundLoadID,*) abund(:nspec,l)
                end do
                rewind(abundLoadID)
            else
                do l=1,points
                    read(abundLoadID,*) abund(:nspec,l)
                    rewind(abundLoadID)
                end do
            end if
        end if
    end subroutine readInputAbunds

    subroutine finalOutput
        !f2py integer, intent(aux) :: points
        if (writeAbunds) then
            if (enable_radiative_transfer .AND. points>1) then
                do dstep=points,1,-1
                    ! WRITE(abundSaveID,*) fhe,fc,fo,fn,fs,fmg
                    write(abundSaveID,8010) abund(:nspec+2,dstep)
                8010  format((999(1pe15.5,:,",")))
                end do
            else
                do dstep=1,points
                    write(abundSaveID,8010) abund(:nspec+2,dstep)
                end do
            end if
        end if
    end subroutine finalOutput

    subroutine output(returnArray,writerates,successflag,physicsarray, chemicalabunarray, rateConstantsArray, &
            & heatarray, statsarray, levelpopulationsarray, sestatsarray, dtime, timePoints)
        real(dp), intent(out), dimension(:, :, :), optional :: physicsarray
        real(dp), intent(out), dimension(:, :, :), optional :: chemicalabunarray
        real(dp), intent(out), dimension(:, :, :), optional :: rateConstantsArray
        real(dp), intent(out), dimension(:, :, :), optional :: heatarray
        real(dp), intent(out), dimension(:, :, :), optional :: statsarray
        real(dp), intent(out), dimension(:, :, :), optional :: levelpopulationsarray
        real(dp), intent(out), dimension(:, :, :), optional :: sestatsarray
        integer, intent(in), optional :: dtime, timePoints
        logical, intent(in) :: returnArray, writerates

        integer, intent(out) :: successflag
        successflag = 0
        if (returnArray) then
            ! Try to catch out of bounds errors before they create a segfault
            if (dtime > timepoints+1) then
                write(*,*) "Ran out of timepoints in arrays, trying to stop gracefully"
                successflag=NOT_ENOUGH_TIMEPOINTS_ERROR
                return
            end if

            physicsarray(dtime, dstep, 1) = timeInYears
            physicsarray(dtime, dstep, 2) = density(dstep)
            physicsarray(dtime, dstep, 3) = gasTemp(dstep)
            physicsarray(dtime, dstep, 4) = dustTemp(dstep)
            physicsarray(dtime, dstep, 5) = av(dstep)
            physicsarray(dtime, dstep, 6) = radfield
            physicsarray(dtime, dstep, 7) = zeta
            physicsarray(dtime, dstep, 8) = dstep
            physicsarray(dtime, dstep, 9) = parcel_radius(dstep)
            physicsarray(dtime, dstep, 10) = av_internal(dstep)
            physicsarray(dtime, dstep, 11) = radfield_internal(dstep)
            chemicalabunarray(dtime, dstep, :) = abund(1:nspec,dstep)
            ! DVODE solver statistics are now written in chemistry.f90
            ! after each solver attempt (including retries)

            ! Level populations (SIZE-based check - don't use PRESENT)
            if (SIZE(levelpopulationsarray, 1) >= timePoints .AND. heatingFlag) then
                call WRITE_LEVEL_POPULATIONS(levelpopulationsarray, dtime, dstep)
            end if

            ! SE solver statistics (SIZE-based check - don't use PRESENT)
            if (SIZE(sestatsarray, 1) >= timePoints .AND. heatingFlag) then
                call WRITE_SE_STATISTICS(sestatsarray, dtime, dstep)
            end if
        else if (fullOutput .AND. .NOT. returnArray) then
            write(outputId,8020) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep),&
                & av(dstep),radfield,zeta,dstep,abund(1:nspec,dstep)
            8020 format(1pe11.3,",",1pe11.4,",",0pf8.2,",",0pf8.2,",",1pe11.4,",",1pe11.4,&
                ",",1pe11.4,",",I4,",",(999(1pe15.5,:,",")))
        end if
        if (writerates) then
            if (returnArray) then
                ! If returnArray is true, we write the rate constants to the rate constants array
                ! We compute the flux in Python.
                rateConstantsArray(dtime, dstep, :) = rate_constants(:nreac)
                ! Only populate the heating array if it is present, properly sized, AND heating is enabled
                if (SIZE(heatarray, 1) >= timePoints .AND. heatingFlag) then
                    heatarray(dtime, dstep, 1) = timeInYears
                    heatarray(dtime, dstep, 2:(1+NCOOLING)) = coolingValues(:)
                    ! Write all NCOOLANTS line cooling terms
                    heatarray(dtime, dstep, (2+NCOOLING):(1+NCOOLING+NCOOLANTS)) = lineCoolingArray(median_line_index, :NCOOLANTS)

                    ! Heating terms (NHEATING values)
                    heatarray(dtime, dstep, (2+NCOOLING+NCOOLANTS):(1+NCOOLING+NCOOLANTS+NHEATING)) = heatingValues(:)
                    ! Chemical heating (1 value)
                    heatarray(dtime, dstep, 2+NCOOLING+NCOOLANTS+NHEATING) = chemheating
                end if
            else
                ! Else, we write the rate constants and rates to the file.
                if (rateConstantsOutput) then
                    write(rateConstantId,8021) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep), &
                        & av(dstep),radfield,zeta,dstep,rate_constants
                    8021 format(1pe11.3,",",1pe11.4,",",0pf8.2,",",0pf8.2,",",1pe11.4,",",1pe11.4,",",1pe11.4,",", &
                        & I4,",",(9999(1pe15.5e3,:,",")))
                end if
                if (fluxOutput) then
                    write(ratesId,8022) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep), &
                        & av(dstep),radfield,zeta,dstep,REACTIONRATE
                    8022 format(1pe11.3,",",1pe11.4,",",0pf8.2,",",0pf8.2,",",1pe11.4,",",1pe11.4,&
                        & ",",1pe11.4,",",I4,",",(9999(1pe15.5e3,:,",")))
                end if
                if (heatingOutput) then
                    ! Write: time, cooling values (4), line cooling array (NCOOLANTS), heating values (8), chem heating
                    write(heatingId,8023) timeInYears, coolingValues(:), &
                        lineCoolingArray(median_line_index, :NCOOLANTS), &
                        heatingValues(:), chemheating
                    8023 format(1PE16.6E3,:,(999(",",1PE16.6E3)))
                end if
            end if
        end if

        !Every 'writestep' timesteps, write the chosen species out to separate file
        !choose species you're interested in by looking at parameters.f90
        if (.NOT. PRESENT(dtime)) then
            if (writeCounter==writeStep .and. columnOutput) then
                writeCounter=1
                write(columnId,8030) timeInYears,density(dstep),gasTemp(dstep),dustTemp(dstep),&
                &av(dstep),radfield,zeta,abund(outIndx,dstep)
                8030  format(1pe11.3,",",1pe11.4,",",0pf8.2,",",0pf8.2,",",1pe11.4,",",1pe11.4,&
                &",",1pe11.4,",",(999(1pe15.5,:,",")))
            else
                writeCounter=writeCounter+1
            end if
        end if
    end subroutine output

    subroutine closeFiles
        close(outputId)
        close(rateConstantId)
        close(columnId)
        close(abundSaveID)
        close(abundLoadID)
        close(heatingId)  ! heating rates file

    end subroutine closeFiles

    subroutine debugout
        open(debugId,file="output/debuglog",status="unknown", action="write")       !debug file.
        write(debugId,*) "Integrator failed, printing relevant debugging information"
        write(debugId,*) "dens",density(dstep)
        write(debugId,*) "gas temperature in integration array",abund(nspec+1,dstep)
        write(debugId,*) "density in integration array",abund(nspec+2,dstep)
        write(debugId,*) "Av", av(dstep)
        write(debugId,*) "Temp", gasTemp(dstep)
        do i=1,nreac
            if (rate_constants(i) >= huge(i)) write(debugId,*) "rate_constant(",i,") is potentially infinite"
        end do
    end subroutine debugout

    subroutine simpleDebug(message)
    !A very simply subroutine for debugging, just write a bunch of variables to screen
    !so we can check they're all as expected.
    !Argument message is a string to print before the variables
        character(LEN=*), intent(in) :: message
        write(*,*) message
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
    end subroutine simpleDebug

    subroutine WRITE_LEVEL_POPULATIONS(levelpopulationsarray, dtime, dstep)
        real(dp), dimension(:, :, :), intent(inout) :: levelpopulationsarray
        integer, intent(in) :: dtime, dstep
        integer :: N, level_offset, i

        level_offset = 0
        do N = 1, NCOOLANTS
            ! Safety net: prevent out-of-bounds write if levels exceed array size
            if (level_offset + coolants(N)%NLEVEL > SIZE(levelpopulationsarray, 3)) then
                write(*,"(A,I3,A,A,A,I6,A,I6,A,I6)") &
                    "[LEVPOP] ERROR coolant ", N, " (", TRIM(coolants(N)%NAME), &
                    ") would exceed array dim3: offset=", level_offset, &
                    " nlevel=", coolants(N)%NLEVEL, &
                    " dim3=", SIZE(levelpopulationsarray, 3)
                return
            end if
            do i = 1, coolants(N)%NLEVEL
                levelpopulationsarray(dtime, dstep, level_offset + i) = coolants(N)%POPULATION(i)
            end do
            level_offset = level_offset + coolants(N)%NLEVEL
        end do
    end subroutine WRITE_LEVEL_POPULATIONS

    subroutine WRITE_SE_STATISTICS(sestatsarray, dtime, dstep)
        real(dp), dimension(:, :, :), intent(inout) :: sestatsarray
        integer, intent(in) :: dtime, dstep
        integer :: N, idx

        do N = 1, NCOOLANTS
            idx = (N-1) * 3  ! 3 stats per coolant
            sestatsarray(dtime, dstep, idx + 1) = MERGE(1.0_dp, 0.0_dp, coolants(N)%CONVERGED)
            sestatsarray(dtime, dstep, idx + 2) = DBLE(se_coolant_iterations(N))
            sestatsarray(dtime, dstep, idx + 3) = se_coolant_max_rel_change(N)
        end do
    end subroutine WRITE_SE_STATISTICS

end module IO
