!The interface between main fortran code and python.
!wrap.f90 subroutines are all accessible through the python wrap.
! Core algorithm is found in solveAbundances subroutine below, all others call it
module uclchemwrap
    ! allow(use-all)
    use chemistry
    use constants, only: dp
    use f2py_constants, only: nspec, nReac, NCOOLANTS, N_DVODE_STATS, N_PHYSICS_PARAMS, &
        N_TOTAL_LEVELS, N_SE_STATS_PER_COOLANT
    use heating, only: nHeatingTerms
    use io, only: fileSetup, readInputAbunds, finalOutput, output, closeFiles, outIndx, &
        outSpecies, nout, outputId, columnId, rateConstantId, ratesId, heatingId, abundLoadID, &
        abundSaveID, columnOutput, fullOutput, rateConstantsOutput, fluxOutput, &
        readAbunds, writeAbunds, heatingOutput
    !f2py INTEGER, parameter :: dp
    use physicscore, only: coreInitializePhysics
    use postprocess_mod, only: postprocess_error

    implicit none

    private
    public :: cloud, jshock, cshock, postprocess, collapse, hot_core, get_rates, &
        get_odes, get_specname, get_coolant_restart_mode_wrap, set_coolant_restart_mode_wrap

    character(LEN=1) :: dummyString = ""

contains
    subroutine cloud(dictionary, outSpeciesIn,returnArray,returnRateConstants,&
            &givestartabund,timePoints,gridPoints,physicsarray,chemicalabunarray,&
            &rateConstantsArray, heatarray, statsarray, levelpopulationsarray, sestatsarray,&
            &abundanceStart ,abundance_out,specname_out,successFlag)
        !f2py threadsafe
        !Subroutine to call a cloud model, used to interface with python
        ! Loads cloud specific subroutines and send to solveAbundances
        !
        !Args:
        ! dictionary - python parameter dictionary
        ! outSpeciesIn - list of species to output as a space separated string
        ! returnArray - boolean on whether arrays will be returned
        ! returnRateConstants - boolean on whether to write rates to file, or return to memory (if returnArray is True)
        ! givestartabund -  boolean on whether starting abundances were given
        ! gridPoints - number of points uclchem should simulate
        ! physicsarray - array to be filled with physical information for each timestep
        ! chemicalabunarray - array to be filled with chemical abundances for each timestep
        ! abundanceStart - array containing starting chemical conditions
        !Returns:
        ! abundance_out - list of abundances of species in outSpeciesIn
        ! specname_out - array of species that are in the chemicalabunarray
        ! successFlag - integer flag indicating success or fail

        use cloud_mod, only: initializePhysics, updatePhysics, updateTargetTime, sublimation

        !f2py integer, intent(aux) :: nspec, n_physics_params, nHeatingTerms, N_DVODE_STATS, N_TOTAL_LEVELS, N_SE_STATS_PER_COOLANT
        !f2py intent(out) abundance_out, specname_out
        character(LEN=*), intent(inout) :: dictionary, outSpeciesIn
        real(dp), intent(out), dimension(nspec) :: abundance_out
        character(LEN=32), intent(out) :: specname_out(nspec)
        integer, intent(out) :: successFlag
        !f2py intent(in) dictionary,outSpeciesIn
        !f2py intent(out) successFlag
        logical, intent(in) :: returnArray
        !f2py intent(in) returnArray
        logical, intent(in) :: returnRateConstants
        !f2py intent(in) returnRateConstants
        logical, intent(in) :: givestartabund
        !f2py intent(in) givestartabund
        integer, intent(in) :: gridPoints
        !f2py intent(in) gridPoints
        integer, intent(in) :: timePoints
        !f2py intent(in) timePoints

        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, n_physics_params) :: physicsarray
        !f2py intent(in,out) physicsarray
        !f2py depend(timePoints,gridPoints, n_physics_params) physicsarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nspec) :: chemicalabunarray
        !f2py intent(in,out) chemicalabunarray
        !f2py depend(timePoints,gridPoints,nspec) chemicalabunarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nReac) :: rateConstantsArray
        !f2py intent(in,out) rateConstantsArray
        !f2py depend(timePoints,gridPoints,nReac) rateConstantsArray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nHeatingTerms) :: heatarray
        !f2py intent(in,out) heatarray
        !f2py depend(timePoints,gridPoints,nHeatingTerms) heatarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_DVODE_STATS) :: statsarray
        !f2py intent(in,out) statsarray
        !f2py depend(timePoints,gridPoints,N_DVODE_STATS) statsarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_TOTAL_LEVELS) :: levelpopulationsarray
        !f2py intent(in,out) levelpopulationsarray
        !f2py depend(timePoints,gridPoints,N_TOTAL_LEVELS) levelpopulationsarray
        real(dp), intent(inout), optional, &
            & dimension(timePoints+1, gridPoints, NCOOLANTS*N_SE_STATS_PER_COOLANT) :: sestatsarray
        !f2py intent(in,out) sestatsarray
        !f2py depend(timePoints,gridPoints,NCOOLANTS,N_SE_STATS_PER_COOLANT) sestatsarray
        real(dp), intent(in), optional, dimension(gridPoints, nspec) :: abundanceStart
        !f2py intent(in) abundanceStart
        !f2py depend(gridPoints, nspec) abundanceStart


        successFlag=0
        specname_out = specName
        call solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray,returnRateConstants,givestartabund,&
                &timePoints,physicsarray,chemicalabunarray,rateConstantsArray, heatarray, statsarray,&
                &levelpopulationsarray, sestatsarray, abundanceStart)
        if ((ALLOCATED(outIndx)) .and. (successFlag == 0)) then
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        end if

    end subroutine cloud

    subroutine collapse(collapseIn,dictionary,outSpeciesIn,&
            &returnArray,returnRateConstants,givestartabund,timePoints,gridPoints,physicsarray,chemicalabunarray,&
            &rateConstantsArray, heatarray, statsarray, levelpopulationsarray, sestatsarray,&
            &abundanceStart, abundance_out,specname_out,successFlag)
        !f2py threadsafe
        !Subroutine to call a collapse model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! collapseIn - integer indicating which collapse mode to run
        ! dictionary - python parameter dictionary
        ! outSpeciesIn - list of species to output as a space separated string
        ! returnArray - boolean on whether arrays will be returned
        ! givestartabund -  boolean on whether starting abundances were given
        ! gridPoints - number of points uclchem should simulate
        ! physicsarray - array to be filled with physical information for each timestep
        ! chemicalabunarray - array to be filled with chemical abundances for each timestep
        ! abundanceStart - array containing starting chemical conditions
        !Returns:
        ! abundance_out - list of abundances of species in outSpeciesIn
        ! specname_out - array of species that are in the chemicalabunarray
        ! successFlag - integer flag indicating success or fail

        use collapse_mod, only: initializePhysics, updatePhysics, updateTargetTime, sublimation, &
            collapse_mode

        !f2py integer,parameter intent(aux) nspec, n_physics_params, nHeatingTerms, N_DVODE_STATS, N_TOTAL_LEVELS, N_SE_STATS_PER_COOLANT
        character(LEN=*), intent(inout) :: dictionary, outSpeciesIn
        real(dp), intent(out), dimension(nspec) :: abundance_out
        character(LEN=32), intent(out) :: specname_out(nspec)
        integer, intent(out) :: successFlag
        integer, intent(in) :: collapseIn
        !f2py intent(in) collapseIn,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,specname_out,successFlag
        logical, intent(in) :: returnArray
        !f2py intent(in) returnArray
        logical, intent(in) :: returnRateConstants
        !f2py intent(in) returnRateConstants
        logical, intent(in) :: givestartabund
        !f2py intent(in) givestartabund
        integer, intent(in) :: gridPoints
        !f2py intent(in) gridPoints
        integer, intent(in) :: timePoints
        !f2py intent(in) timePoints
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, n_physics_params) :: physicsarray
        !f2py intent(in,out) physicsarray
        !f2py depend(gridPoints) physicsarray
        real(dp), intent(inout), dimension(timePoints+1, gridPoints, nspec), optional :: chemicalabunarray
        !f2py intent(in,out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nReac) :: rateConstantsArray
        !f2py intent(in,out) rateConstantsArray
        !f2py depend(timePoints,gridPoints,nReac) rateConstantsArray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nHeatingTerms) :: heatarray
        !f2py intent(in,out) heatarray
        !f2py depend(timePoints,gridPoints,nHeatingTerms) heatarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_DVODE_STATS) :: statsarray
        !f2py intent(in,out) statsarray
        !f2py depend(timePoints,gridPoints,N_DVODE_STATS) statsarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_TOTAL_LEVELS) :: levelpopulationsarray
        !f2py intent(in,out) levelpopulationsarray
        !f2py depend(timePoints,gridPoints,N_TOTAL_LEVELS) levelpopulationsarray
        real(dp), intent(inout), optional, &
            & dimension(timePoints+1, gridPoints, NCOOLANTS*N_SE_STATS_PER_COOLANT) :: sestatsarray
        !f2py intent(in,out) sestatsarray
        !f2py depend(timePoints,gridPoints,NCOOLANTS,N_SE_STATS_PER_COOLANT) sestatsarray
        real(dp), intent(in), optional, dimension(gridPoints, nspec) :: abundanceStart
        !f2py intent(in) abundanceStart
        !f2py depend(gridPoints, nspec) abundanceStart
        successFlag=0
        specname_out(:nspec) = specName
        collapse_mode=collapseIn

        call solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray, returnRateConstants,givestartabund,&
                &timepoints,physicsarray,chemicalabunarray,rateConstantsArray, heatarray, statsarray,&
                &levelpopulationsarray, sestatsarray, abundanceStart)

        if ((ALLOCATED(outIndx)) .and. (successFlag == 0)) then
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        end if
    end subroutine collapse

    subroutine hot_core(temp_index,max_temp,dictionary,outSpeciesIn,returnArray,&
            &returnRateConstants,givestartabund,timePoints,gridPoints,physicsarray,chemicalabunarray,&
            &rateConstantsArray, heatarray, statsarray, levelpopulationsarray, sestatsarray,&
            &abundanceStart, abundance_out,specname_out,successFlag)
        !f2py threadsafe
        !Subroutine to call a hot core model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! temp_index - integer indicating which mass hot core to run - see hotcore.90
        ! max_temp - maximum temperature before we stop increasing.
        ! dictionary - python parameter dictionary
        ! outSpeciesIn - list of species to output as a space separated string
        ! returnArray - boolean on whether arrays will be returned
        ! givestartabund -  boolean on whether starting abundances were given
        ! gridPoints - number of points uclchem should simulate
        ! physicsarray - array to be filled with physical information for each timestep
        ! chemicalabunarray - array to be filled with chemical abundances for each timestep
        ! abundanceStart - array containing starting chemical conditions
        !Returns:
        ! abundance_out - list of abundances of species in outSpeciesIn
        ! successFlag - integer flag indicating success or fail
        use hotcore, only: initializePhysics, updatePhysics, updateTargetTime, sublimation, &
            maxTemp, tempIndx


        !f2py integer, parameter intent(aux) nspec, n_physics_params, nHeatingTerms, N_DVODE_STATS
        character(LEN=*), intent(inout) :: dictionary, outSpeciesIn
        real(dp), intent(out) :: abundance_out(nspec)
        real(dp), intent(in) :: max_temp
        integer, intent(in) :: temp_index
        integer, intent(out) :: successFlag
        character(LEN=32), intent(out) :: specname_out(nspec)
        !f2py intent(in) temp_index,max_temp,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,specname_out,successFlag
        logical, intent(in) :: returnArray
        !f2py intent(in) returnArray
        logical, intent(in) :: returnRateConstants
        !f2py intent(in) returnRateConstants
        logical, intent(in) :: givestartabund
        !f2py intent(in) givestartabund
        integer, intent(in) :: gridPoints
        !f2py intent(in) gridPoints
        integer, intent(in) :: timePoints
        !f2py intent(in) timePoints
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, n_physics_params) :: physicsarray
        !f2py intent(in, out) physicsarray
        !f2py depend(gridPoints) physicsarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nspec) :: chemicalabunarray
        !f2py intent(in, out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nReac) :: rateConstantsArray
        !f2py intent(in,out) rateConstantsArray
        !f2py depend(timePoints,gridPoints,nReac) rateConstantsArray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nHeatingTerms) :: heatarray
        !f2py intent(in,out) heatarray
        !f2py depend(timePoints,gridPoints,nHeatingTerms) heatarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_DVODE_STATS) :: statsarray
        !f2py intent(in,out) statsarray
        !f2py depend(timePoints,gridPoints,N_DVODE_STATS) statsarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_TOTAL_LEVELS) :: levelpopulationsarray
        !f2py intent(in,out) levelpopulationsarray
        !f2py depend(timePoints,gridPoints,N_TOTAL_LEVELS) levelpopulationsarray
        real(dp), intent(inout), optional, &
            & dimension(timePoints+1, gridPoints, NCOOLANTS*N_SE_STATS_PER_COOLANT) :: sestatsarray
        !f2py intent(in,out) sestatsarray
        !f2py depend(timePoints,gridPoints,NCOOLANTS,N_SE_STATS_PER_COOLANT) sestatsarray
        real(dp), intent(in), optional, dimension(gridPoints, nspec) :: abundanceStart
        !f2py intent(in) abundanceStart
        !f2py depend(gridPoints, nspec) abundanceStart
        specname_out(:nspec) = specName
        maxTemp=max_temp
        tempIndx=temp_index

        call solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
        &updatePhysics,updateTargetTime,sublimation,returnArray, returnRateConstants,givestartabund,&
        &timepoints,physicsarray,chemicalabunarray,rateConstantsArray, heatarray, statsarray,&
        &levelpopulationsarray, sestatsarray, abundanceStart)

        if ((ALLOCATED(outIndx)) .and. (successFlag == 0)) then
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        end if
    end subroutine hot_core

    subroutine cshock(shock_vel,timestep_factor,minimum_temperature,dictionary, outSpeciesIn,&
            &returnArray,returnRateConstants,givestartabund,timePoints,gridPoints,physicsarray,chemicalabunarray,&
            &rateConstantsArray, heatarray, statsarray, levelpopulationsarray, sestatsarray,&
            &abundanceStart, abundance_out,dissipation_time,specname_out,successFlag)
        !f2py threadsafe
        !Subroutine to call a C-shock model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! shock_vel - real(dp) shock velocity in km/s
        ! timestep_factor - Multiply dissipation time of the shock by this real(dp) factor to
        !                  get the timestep for the simulation up until dissipation time is reached
        ! minimum_temperature - float indicating minimum temperature before we stop post shock cooling.
        !                        set to zero to disable.
        ! dictionary - python parameter dictionary
        ! outSpeciesIn - list of species to output as a space separated string
        ! returnArray - boolean on whether arrays will be returned
        ! givestartabund -  boolean on whether starting abundances were given
        ! gridPoints - number of points uclchem should simulate
        ! physicsarray - array to be filled with physical information for each timestep
        ! chemicalabunarray - array to be filled with chemical abundances for each timestep
        ! abundanceStart - array containing starting chemical conditions
        !Returns:
        ! abundance_out - list of abundances of species in outSpeciesIn
        ! dissipation_time - float, dissipation time in years
        ! successFlag - integer flag indicating success or fail
        use cshock_mod, only: initializePhysics, updatePhysics, updateTargetTime, sublimation, &
            minimumPostshockTemp, dissipationTime, timestepFactor, vs

        !f2py integer, parameter intent(aux) nspec, n_physics_params, nHeatingTerms, N_DVODE_STATS
        character(LEN=*), intent(inout) :: dictionary, outSpeciesIn
        real(dp), intent(out) :: abundance_out(nspec)
        real(dp), intent(in) :: shock_vel,timestep_factor
        real(dp), intent(in) :: minimum_temperature
        real(dp), intent(out) :: dissipation_time
        integer, intent(out) :: successFlag
        character(LEN=32), intent(out) :: specname_out(nspec)
        !f2py intent(in) shock_vel,timestep_factor,minimum_temperature,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,dissipation_time,specname_out,successFlag
        logical, intent(in) :: returnArray
        !f2py intent(in) returnArray
        logical, intent(in) :: returnRateConstants
        !f2py intent(in) returnRateConstants
        logical, intent(in) :: givestartabund
        !f2py intent(in) givestartabund
        integer, intent(in) :: gridPoints
        !f2py intent(in) gridPoints
        integer, intent(in) :: timePoints
        !f2py intent(in) timePoints
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, n_physics_params) :: physicsarray
        !f2py intent(in, out) physicsarray
        !f2py depend(gridPoints) physicsarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nspec) :: chemicalabunarray
        !f2py intent(in, out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nReac) :: rateConstantsArray
        !f2py intent(in,out) rateConstantsArray
        !f2py depend(timePoints,gridPoints,nReac) rateConstantsArray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nHeatingTerms) :: heatarray
        !f2py intent(in,out) heatarray
        !f2py depend(timePoints,gridPoints,nHeatingTerms) heatarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_DVODE_STATS) :: statsarray
        !f2py intent(in,out) statsarray
        !f2py depend(timePoints,gridPoints,N_DVODE_STATS) statsarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_TOTAL_LEVELS) :: levelpopulationsarray
        !f2py intent(in,out) levelpopulationsarray
        !f2py depend(timePoints,gridPoints,N_TOTAL_LEVELS) levelpopulationsarray
        real(dp), intent(inout), optional, &
            & dimension(timePoints+1, gridPoints, NCOOLANTS*N_SE_STATS_PER_COOLANT) :: sestatsarray
        !f2py intent(in,out) sestatsarray
        !f2py depend(timePoints,gridPoints,NCOOLANTS,N_SE_STATS_PER_COOLANT) sestatsarray
        real(dp), intent(in), optional, dimension(gridPoints, nspec) :: abundanceStart
        !f2py intent(in) abundanceStart
        !f2py depend(gridPoints, nspec) abundanceStart

        vs=shock_vel
        timestepFactor=timestep_factor
        minimumPostshockTemp=minimum_temperature

        specname_out = specName
        successFlag=0
        call solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray, returnRateConstants,givestartabund,&
                &timepoints,physicsarray,chemicalabunarray,rateConstantsArray,heatarray,statsarray,&
                &levelpopulationsarray, sestatsarray,abundanceStart)

        if (successFlag == 0) then
            if (ALLOCATED(outIndx)) abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
            dissipation_time=dissipationTime
        end if
    end subroutine cshock

    subroutine jshock(shock_vel,dictionary,outSpeciesIn,returnArray,returnRateConstants,givestartabund,&
            &timePoints,gridPoints,physicsarray,chemicalabunarray,&
            &rateConstantsArray, heatarray, statsarray, levelpopulationsarray, sestatsarray, &
            &abundanceStart,abundance_out,specname_out,successFlag)
        !f2py threadsafe
        !Subroutine to call a J-shock model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! shock_vel - real(dp) shock velocity in km/s
        ! dictionary - python parameter dictionary
        ! outSpeciesIn - list of species to output as a space separated string
        ! returnArray - boolean on whether arrays will be returned
        ! givestartabund -  boolean on whether starting abundances were given
        ! gridPoints - number of points uclchem should simulate
        ! physicsarray - array to be filled with physical information for each timestep
        ! chemicalabunarray - array to be filled with chemical abundances for each timestep
        ! abundanceStart - array containing starting chemical conditions
        !Returns:
        ! abundance_out - list of abundances of species in outSpeciesIn
        ! successFlag - integer flag indicating success or fail
        use jshock_mod, only: initializePhysics, updatePhysics, updateTargetTime, sublimation, vs

        !f2py integer, parameter intent(aux) nspec, n_physics_params, nHeatingTerms, N_DVODE_STATS
        character(LEN=*), intent(inout) :: dictionary, outSpeciesIn
        real(dp), intent(in) :: shock_vel
        real(dp), intent(out), dimension(nspec) :: abundance_out
        integer, intent(out) :: successFlag
        character(LEN=32), intent(out) :: specname_out(nspec)
        !f2py intent(in) shock_vel,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out, specname_out, successFlag
        logical, intent(in) :: returnArray
        !f2py intent(in) returnArray
        logical, intent(in) :: returnRateConstants
        !f2py intent(in) returnRateConstants
        logical, intent(in) :: givestartabund
        !f2py intent(in) givestartabund
        integer, intent(in) :: gridPoints
        !f2py intent(in) gridPoints
        integer, intent(in) :: timePoints
        !f2py intent(in) timePoints
        real(dp), intent(inout), optional, dimension(timePoints + 1, gridPoints, n_physics_params) :: physicsarray
        !f2py intent(in, out) physicsarray
        !f2py depend(gridPoints) physicsarray
        real(dp), intent(inout), optional, dimension(timePoints + 1, gridPoints, nspec) :: chemicalabunarray
        !f2py intent(in, out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nReac) :: rateConstantsArray
        !f2py intent(in,out) rateConstantsArray
        !f2py depend(timePoints,gridPoints,nReac) rateConstantsArray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nHeatingTerms) :: heatarray
        !f2py intent(in,out) heatarray
        !f2py depend(timePoints,gridPoints,nHeatingTerms) heatarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_DVODE_STATS) :: statsarray
        !f2py intent(in,out) statsarray
        !f2py depend(timePoints,gridPoints,N_DVODE_STATS) statsarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_TOTAL_LEVELS) :: levelpopulationsarray
        !f2py intent(in,out) levelpopulationsarray
        !f2py depend(timePoints,gridPoints,N_TOTAL_LEVELS) levelpopulationsarray
        real(dp), intent(inout), optional, &
            dimension(timePoints+1, gridPoints, NCOOLANTS*N_SE_STATS_PER_COOLANT) :: sestatsarray
        !f2py intent(in,out) sestatsarray
        !f2py depend(timePoints,gridPoints,NCOOLANTS,N_SE_STATS_PER_COOLANT) sestatsarray
        real(dp), intent(in), optional, dimension(gridPoints, nspec) :: abundanceStart
        !f2py intent(in) abundanceStart
        !f2py depend(gridPoints, nspec) abundanceStart
        vs=shock_vel
        call solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
        &updatePhysics,updateTargetTime,sublimation,returnArray, returnRateConstants,givestartabund,&
        &timepoints,physicsarray,chemicalabunarray,rateConstantsArray,heatarray,statsarray,&
        &levelpopulationsarray, sestatsarray,abundanceStart)

        specname_out(:nspec) = specName
        if ((ALLOCATED(outIndx)) .and. (successFlag == 0)) then
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        end if
    end subroutine jshock

    subroutine postprocess(dictionary,outSpeciesIn,returnArray,returnRateConstants,givestartabund,&
        &timePoints,gridPoints,physicsarray,chemicalabunarray,rateConstantsArray,heatarray,statsarray,&
        &levelpopulationsarray, sestatsarray, abundanceStart,timegrid,densgrid,gastempgrid,&
        &dusttempgrid,radfieldgrid,zetagrid,useav,avgrid,&
        &usecoldens,nhgrid,nh2grid,ncogrid,ncgrid,abundance_out,specname_out,successFlag)
        !f2py threadsafe
        !Subroutine to call a J-shock model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! shock_vel - real(dp) shock velocity in km/s
        ! dictionary - python parameter dictionary
        ! outSpeciesIn - list of species to output as a space separated string
        ! returnArray - boolean on whether arrays will be returned
        ! givestartabund -  boolean on whether starting abundances were given
        ! gridPoints - number of points uclchem should simulate
        ! physicsarray - array to be filled with physical information for each timestep
        ! chemicalabunarray - array to be filled with chemical abundances for each timestep
        ! abundanceStart - array containing starting chemical conditions
        !Returns:
        ! abundance_out - list of abundances of species in outSpeciesIn
        ! successFlag - integer flag indicating success or fail
        use postprocess_mod, only: initializePhysics, updatePhysics, updateTargetTime, sublimation

        !f2py integer, parameter intent(aux) nspec, n_physics_params, nHeatingTerms, N_DVODE_STATS
        character(LEN=*), intent(inout) :: dictionary, outSpeciesIn
        real(dp), intent(out) :: abundance_out(nspec)
        integer, intent(out) :: successFlag
        character(LEN=32), intent(out) :: specname_out(nspec)
        !f2py intent(in) shock_vel,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out, specname_out, successFlag
        logical, intent(in) :: returnArray
        !f2py intent(in) returnArray
        logical, intent(in) :: returnRateConstants
        !f2py intent(in) returnRateConstants
        logical, intent(in) :: givestartabund
        !f2py intent(in) givestartabund
        integer, intent(in) :: gridPoints
        !f2py intent(in) gridPoints
        integer, intent(in) :: timePoints
        !f2py intent(in) timePoints
        logical, intent(in) :: useav
        !f2py intent(in) useav
        logical, intent(in) :: usecoldens
        !f2py intent(in) usecoldens
        real(dp), intent(inout), optional, dimension(timePoints + 1, gridPoints, n_physics_params) :: physicsarray
        !f2py intent(in, out) physicsarray
        !f2py depend(gridPoints) physicsarray
        real(dp), intent(inout), optional, dimension(timePoints + 1, gridPoints, nspec) :: chemicalabunarray
        !f2py intent(in, out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nReac) :: rateConstantsArray
        !f2py intent(in,out) rateConstantsArray
        !f2py depend(timePoints,gridPoints,nReac) rateConstantsArray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, nHeatingTerms) :: heatarray
        !f2py intent(in,out) heatarray
        !f2py depend(timePoints,gridPoints,nHeatingTerms) heatarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_DVODE_STATS) :: statsarray
        !f2py intent(in,out) statsarray
        !f2py depend(timePoints,gridPoints,N_DVODE_STATS) statsarray
        real(dp), intent(inout), optional, dimension(timePoints+1, gridPoints, N_TOTAL_LEVELS) :: levelpopulationsarray
        !f2py intent(in,out) levelpopulationsarray
        !f2py depend(timePoints,gridPoints,N_TOTAL_LEVELS) levelpopulationsarray
        real(dp), intent(inout), optional, &
            & dimension(timePoints+1, gridPoints, NCOOLANTS*N_SE_STATS_PER_COOLANT) :: sestatsarray
        !f2py intent(in,out) sestatsarray
        !f2py depend(timePoints,gridPoints,NCOOLANTS,N_SE_STATS_PER_COOLANT) sestatsarray
        real(dp), intent(in), dimension(gridPoints, nspec), optional :: abundanceStart
        !f2py intent(in) abundanceStart
        !f2py depend(gridPoints, nspec) abundanceStart
        real(dp), intent(in), dimension(timePoints) :: timegrid
        real(dp), intent(in), dimension(timePoints) :: densgrid
        real(dp), intent(in), dimension(timePoints) :: gastempgrid
        real(dp), intent(in), dimension(timePoints) :: dusttempgrid
        real(dp), intent(in), dimension(timePoints) :: radfieldgrid
        real(dp), intent(in), dimension(timePoints) :: zetagrid
        real(dp), intent(in), dimension(timePoints), optional :: avgrid
        real(dp), intent(in), dimension(timePoints), optional :: nhgrid
        real(dp), intent(in), dimension(timePoints), optional :: nh2grid
        real(dp), intent(in), dimension(timePoints), optional :: ncogrid
        real(dp), intent(in), dimension(timePoints), optional :: ncgrid
        !f2py  intent(in) timegrid,densgrid,gastempgrid,dusttempgrid,radfieldgrid,zetagrid,useav,avgrid
        !f2py  intent(in) nhgrid,nh2grid,ncogrid,ncgrid

        successFlag=0

        if (usecoldens) then
            call solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray,returnRateConstants,givestartabund,&
                &timepoints,physicsarray,chemicalabunarray,rateConstantsArray,heatarray,statsarray,&
                &levelpopulationsarray, sestatsarray, abundanceStart,&
                &timegrid,densgrid,gastempgrid,dusttempgrid,radfieldgrid,zetagrid,useav,avgrid,&
                &usecoldens,nhgrid,nh2grid,ncogrid,ncgrid)
        else
            call solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray,returnRateConstants,givestartabund,&
                &timepoints,physicsarray,chemicalabunarray,rateConstantsArray,heatarray,statsarray,&
                &levelpopulationsarray, sestatsarray, abundanceStart,&
                &timegrid,densgrid,gastempgrid,dusttempgrid,radfieldgrid,zetagrid,useav,avgrid,&
                &usecoldens)
        end if
        specname_out(:nspec) = specName
        if ((ALLOCATED(outIndx)) .and. (successFlag == 0)) then
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        end if
    end subroutine postprocess

    subroutine get_rates(dictionary,abundancesIn,speciesIndx,rateIndxs,&
        &speciesRates,successFlag,transfer,swap,bulk_layers)
        !Given a species of interest, some parameters and abundances, this subroutine
        !returns the rate of all reactions that include that species plus some extra variables
        !to allow for the calculation of the rate of bulk/surface ice transfer.
        use cloud_mod, only: initializePhysics
        use f2py_constants, only: nspec, nReac
        !f2py integer, intent(aux) :: nspec
        character(LEN=*), intent(inout) :: dictionary
        real(dp), intent(in) :: abundancesIn(nspec)
        integer, intent(in) :: speciesIndx
        integer, intent(in), dimension(nReac) :: rateIndxs
        real(dp), intent(out), dimension(nReac) :: speciesRates
        real(dp), intent(out) :: transfer, swap, bulk_layers
        integer, intent(out) :: successFlag

        real(dp), dimension(nspec+1) :: ydot
        real(dp) :: surfaceCoverage
        integer :: speci,bulk_version,surface_version
        !f2py intent(in) dictionary,abundancesIn,speciesIndx,rateIndxs
        !f2py intent(out) speciesRates,successFlag,transfer,swap,bulk_layers
        ! INCLUDE 'defaultparameters.f90'

        call dictionaryParser(dictionary, dummyString, successFlag)
        if (successFlag < 0) then
            write(*,*) "Error reading parameter dictionary"
            return
        end if
        call coreInitializePhysics(successFlag)
        call initializePhysics(successFlag)
        if (successFlag < 0) then
            write(*,*) "Error initializing physics"
            return
        end if

        call initializeChemistry(readAbunds, successFlag)
        if (successFlag < 0) return
        dstep=1
        abund(:nspec,dstep)=abundancesIn(:nspec)
        abund(neq,dstep)=initialDens
        currentTime=0.0
        timeInYears=0.0

        targetTime=1.0e-7_dp
        call updateChemistry(successFlag)

        call F(NEQ,currentTime,abund(:,dstep),ydot)

        speciesRates=rate(rateIndxs)

        if ((specname(speciesIndx)(1:1) == "#") .or.&
        & (specname(speciesIndx)(1:1) == "@")) then
            do speci=1,nspec
                if (specname(speci) == "@"//specname(speciesIndx)(2:)) bulk_version=speci
                if (specname(speci) == "#"//specname(speciesIndx)(2:)) surface_version=speci
            end do
            if (SURFGROWTHUNCORRECTED < 0) then
                surfaceCoverage = MIN(1.0_dp, safeBulk/safeMantle)
                transfer=SURFGROWTHUNCORRECTED*surfaceCoverage*abund(bulk_version,1)/safeBulk
            else
                surfaceCoverage = bulkGainFromMantleBuildUp()
                transfer=SURFGROWTHUNCORRECTED*surfaceCoverage*abund(surface_version,1)
            end if
            swap=totalSwap
            bulk_layers=bulkLayersReciprocal
        else
            swap=0.0
            transfer=0.0
            bulk_layers=0.0
        end if

    end subroutine get_rates

    subroutine get_odes(dictionary,abundancesIn,ratesOut)
        !Obtain the ODE values for some given parameters and abundances.
        !Essentially runs one time step of solveAbundances  then calls the ODE subroutine (F)
        use cloud_mod, only: initializePhysics
        use f2py_constants, only: nspec
        !f2py integer, intent(aux) :: nspec
        character(LEN=*), intent(inout) :: dictionary
        real(dp), intent(in) :: abundancesIn(nspec)
        real(dp), intent(out) :: ratesOut(nspec+1)
        integer :: successFlag
        !f2py intent(in) :: dictionary, abundancesIn
        !f2py intent(out) :: ratesOut
        call dictionaryParser(dictionary, dummyString, successFlag)

        call coreInitializePhysics(successFlag)
        call initializePhysics(successFlag)
        if (successFlag < 0) then
            write(*,*) "Error initializing physics"
            return
        end if

        call initializeChemistry(readAbunds, successFlag)
        if (successFlag < 0) return
        dstep=1
        abund(:nspec,dstep)=abundancesIn(:nspec)
        abund(neq,dstep)=initialDens
        currentTime=0.0
        timeInYears=0.0
        targetTime=1.0e-7_dp
        call updateChemistry(successFlag)
        call F(NEQ,currentTime,abund(:,dstep),ratesOut(:NEQ))
    end subroutine get_odes

    subroutine solveAbundances(dictionary,outSpeciesIn,successFlag,modelInitializePhysics,&
            &modelUpdatePhysics,updateTargetTime, sublimation, returnArray, returnRateConstants,givestartabund,&
            &timePoints, physicsarray, chemicalabunarray,rateConstantsArray,heatarray,statsarray,&
            &levelpopulationsarray,sestatsarray,abundanceStart, timegrid,densgrid,gtempgrid,dtempgrid,radgrid,&
            &zetagrid,useav,avgrid,usecoldens,nhgrid,nh2grid,ncogrid,ncgrid)
        use, intrinsic :: iso_c_binding, only: C_NULL_PTR, C_INT
        ! Core UCLCHEM routine. Solves the chemical equations for a given set of parameters through time
        ! for a specified physical model.
        ! Change behavior of physics by sending different subroutine arguments - hence the need for model subroutines above
        ! dictionary - the parameter dictionary string representing a python dictionary
        ! outSpeciesIn - the species to output to separate file
        ! successFlag - Integer to indicate whether code completed successfully
        ! modelInitializePhysics - subroutine to initialize physics from a physics module
        ! modelUpdatePhysics - subroutine to update physics from a physics module
        ! updateTargetTime - subroutine to update the target time from a physics module
        ! sublimation - subroutine allowing physics module to directly modify abundances once per time step.
        ! returnArray - boolean on whether arrays will be returned
        ! givestartabund -  boolean on whether starting abundances were given
        ! gridPoints - number of points uclchem should simulate
        ! physicsarray - array to be filled with physical information for each timestep (optional)
        ! chemicalabunarray - array to be filled with chemical abundances for each timestep (optional)
        ! abundanceStart - array containing starting chemical conditions (optional)
        ! USE constants, only : nspec
        !f2py integer, intent(aux) :: nspec, nHeatingTerms, NCOOLANTS, N_DVODE_STATS, N_TOTAL_LEVELS, N_SE_STATS_PER_COOLANT
        character(LEN=*), intent(inout) :: dictionary, outSpeciesIn
        external modelInitializePhysics,updateTargetTime,modelUpdatePhysics,sublimation
        integer, intent(out) :: successFlag
        logical, intent(in) :: returnArray, givestartabund, returnRateConstants
        integer, intent(in) :: timePoints
        integer :: dtime, dtime_local
        real(dp) :: initialTime, outermost_stop_time
        logical :: parcelDone, outermost_stopped
        ! Arrays needed to work return physics in memory mode
        real(dp), dimension(:, :, :), optional, intent(in) :: physicsarray
        real(dp), dimension(:, :, :), optional, intent(in) :: chemicalabunarray
        real(dp), dimension(:, :, :), optional, intent(in) :: rateConstantsArray
        real(dp), dimension(:, :, :), optional, intent(in) :: heatarray
        real(dp), dimension(:, :, :), optional, intent(inout) :: statsarray
        real(dp), dimension(:, :, :), optional, intent(in) :: levelpopulationsarray
        real(dp), dimension(:, :, :), optional, intent(inout) :: sestatsarray
        real(dp), dimension(:, :), optional, intent(in) :: abundanceStart
        ! Arrays needed to work with custom density/temperature profiles
        !  &timegrid,densgrid,gastempgrid,dusttempgrid,nhgrid,nh2grid,ncogrodi,ncgrid)
        real(dp), dimension(:), optional, intent(in) :: timegrid
        real(dp), dimension(:), optional, intent(in) :: densgrid
        real(dp), dimension(:), optional, intent(in) :: gtempgrid
        real(dp), dimension(:), optional, intent(in) :: dtempgrid
        real(dp), dimension(:), optional, intent(in) :: radgrid
        real(dp), dimension(:), optional, intent(in) :: avgrid
        real(dp), dimension(:), optional, intent(in) :: zetagrid
        logical, optional, intent(in) :: useav
        logical, optional, intent(in) :: usecoldens
        real(dp), dimension(:), optional, intent(in) :: nhgrid
        real(dp), dimension(:), optional, intent(in) :: nh2grid
        real(dp), dimension(:), optional, intent(in) :: ncogrid
        real(dp), dimension(:), optional, intent(in) :: ncgrid
        integer(C_INT) :: flush_rc
        interface
            ! allow(function-missing-result)
            integer(C_INT) function c_fflush(stream) bind(C, name="fflush")
                use, intrinsic :: iso_c_binding, only: C_INT, C_PTR

                implicit none

                type(C_PTR), value :: stream
            end function c_fflush
        end interface
        successFlag=0
        ! Debug: Print array dimensions
        ! Set variables to default values
        ! INCLUDE 'defaultparameters.f90'
        !Read input parameters from the dictionary
        call dictionaryParser(dictionary, outSpeciesIn, successFlag)
        if (successFlag < 0) then
            successFlag=PARAMETER_READ_ERROR
            write(*,*) "Error reading parameter dictionary"
            return
        end if
        dstep=1
        currentTime=0.0
        timeInYears=0.0
        if (returnArray) then
            ! Fix to make sure that running in memory mode after running in file mode works correctly
            readAbunds=.false.
            writeAbunds=.false.
        else
            call fileSetup
        end if
        !Initialize core physics first then model specific
        !This allows model to overrule changes made by core
        call coreInitializePhysics(successFlag)
        if (successFlag < 0) then
            successFlag=PHYSICS_INIT_ERROR
            write(*,*) "Error initializing physics"
            return
        end if

        if (present(timegrid)) then
            if (usecoldens) then
                call modelInitializePhysics(successflag, timegrid,densgrid,radgrid,zetagrid,gtempgrid,&
                &dtempgrid,useav,avgrid,usecoldens,timepoints,nhgrid,nh2grid,ncogrid,ncgrid)
            else
                call modelInitializePhysics(successflag, timegrid,densgrid,radgrid,zetagrid,gtempgrid,&
                &dtempgrid,useav,avgrid,usecoldens,timepoints)
            end if
        else
            call modelInitializePhysics(successFlag)
        end if
        if (successFlag < 0) then
            successFlag=PHYSICS_INIT_ERROR
            write(*,*) "Error initializing physics"
            return
        end if

        ! Initialize the chemistry
        call initializeChemistry(readAbunds, successFlag)
        if (successFlag < 0) return
        if (returnArray .AND. givestartabund) then
            ! In case we have custom abundances, set them here
            if (enable_radiative_transfer .AND. points>1) then
                do l=points,1,-1
                    abund(:nspec,l) = abundanceStart(l, :nspec)
                end do
            else
                do l=1,points
                    abund(:nspec,l) = abundanceStart(l, :nspec)
                end do
            end if
        else
            ! Else just use the default readInputAbunds routine:
            call readInputAbunds  !this won't do anything if no abundLoadFile was in input
        end if
        ! Set conservation baseline after all starting abundances are loaded.
        ! For 1D RT the per-parcel baseline is updated inside the parcel loop below.
        if (.NOT. (enable_radiative_transfer .AND. points>1)) call setConservationBaseline()
        !CALL simpleDebug("Initialized")

        dstep = 1
        dtime = 0  ! Start at 0 so first loop iteration writes to dtime=1
        ! Set timeInYears for consistency
        timeInYears = currentTime / SECONDS_PER_YEAR
        ! Don't call output here - let the first loop iteration handle it
        ! This avoids duplicate output at t=0

        ! One-time setup for the (points, time) loop order used in 1D RT models
        if (enable_radiative_transfer .AND. points>1) then
            if (ALLOCATED(coldens_history)) deallocate(coldens_history)
            allocate(coldens_history(timepoints+1))
            coldens_history = 0.0_dp
            initialTime = currentTime
            outermost_stop_time = finalTime * SECONDS_PER_YEAR
            outermost_stopped = .false.
        end if

        if (enable_radiative_transfer .AND. points>1) then
            ! -- (points outer, time inner) loop -------------------------------------
            ! Processes each spatial parcel through its full time evolution before
            ! moving inward. DVODE retains BDF history (ISTATE=2) across all timesteps
            ! for a given parcel. coldens_history stores coldens(dstep+1, dtime) so the
            ! edge-to-core AV accumulation is reproduced exactly without approximation.
            parcel_loop: do dstep=points,1,-1

                ! -- per-parcel reset ------------------------------------------------
                currentTime  = initialTime
                timeInYears  = initialTime / SECONDS_PER_YEAR
                dtime_local  = 0
                parcelDone   = .false.
                coolant_levpop_force_recompute = .true.

                if (givestartabund) abund(:nspec, dstep) = abundanceStart(dstep, :nspec)
                call resetConservationBaselineForPoint(dstep)
                call resetDVODEForNewPoint()

                ! Outermost shell has no outer neighbor; initialize to zero
                if (dstep == points) outer_coldens_for_current_step = 0.0_dp

                ! -- time integration for this parcel --------------------------------
                time_loop: do while ((successFlag == 0) .AND. (timeInYears < finalTime) .AND. (.NOT. parcelDone))

                    ! Modes 1 & 2: stop inner parcels when the outermost parcel stopped
                    if ((parcelStoppingMode == 1 .OR. parcelStoppingMode == 2) &
                            .AND. outermost_stopped .AND. dstep < points) then
                        if (currentTime >= outermost_stop_time) then
                            exit time_loop
                        end if
                    end if

                    dtime_local    = dtime_local + 1
                    currentTimeold = currentTime
                    call updateTargetTime
                    coolant_levpop_force_recompute = .true.
                    if (targetTime/SECONDS_PER_YEAR > finalTime) then
                        exit time_loop
                    end if

                    ! Freefall stopping logic
                    if (freefall .AND. density(dstep) >= density_max(dstep)) then
                        select case (parcelStoppingMode)
                            case (0)
                                density(dstep) = density_max(dstep)
                            case (1)
                                if (dstep == points) then
                                    outermost_stop_time = currentTime
                                    outermost_stopped   = .true.
                                end if
                                density(dstep) = density_max(dstep)
                            case (2)
                                if (dstep == points) then
                                    outermost_stop_time = currentTime
                                    outermost_stopped   = .true.
                                    parcelDone = .true.    ! outermost exits at its own freefall time
                                else
                                    density(dstep) = density_max(dstep)  ! inner: cap density, keep running
                                end if
                            case default
                                if (dstep == points) then
                                    outermost_stop_time = currentTime
                                    outermost_stopped   = .true.
                                end if
                                density(dstep) = density_max(dstep)
                        end select
                    end if

                    if (writeTimestepInfo) then
                        write(*,"(A,I3,A,ES10.2,A,ES10.2,A,F9.1,A,F9.1,A,ES10.2,A,ES10.2,A,ES10.2)") &
                            "p:", dstep, &
                            "t:", currentTimeold/SECONDS_PER_YEAR, &
                            " dT:", dvode_rstats(11)/SECONDS_PER_YEAR, &
                            " Tg:", gasTemp(dstep), &
                            " Td:", dustTemp(dstep), &
                            " nH:", density(dstep), &
                            " G0:", radfield, &
                            " z:", zeta
                        flush(6)
                        flush_rc = c_fflush(C_NULL_PTR)
                    end if

                    if (.NOT. parcelDone) then
                        ! Load outer neighbor's coldens at this timestep for the AV accumulation
                        if (dstep < points) then
                          outer_coldens_for_current_step = coldens_history(dtime_local)
                        end if

                        call updateChemistry(successFlag, statsarray, timePoints+1, dtime_local)
                        if (successFlag < 0) then
                            write(*,*) "Error updating chemistry"
                            return
                        end if
                        if (coolant_error_flag /= 0) then
                            write(*,*) "Coolant error: ", TRIM(coolant_error_message)
                            successFlag = coolant_error_flag
                            coolant_error_flag = 0
                            return
                        end if
                        timeInYears = targetTime / SECONDS_PER_YEAR

                        call coreUpdatePhysics
                        call modelUpdatePhysics
                        if (postprocess_error /= 0) then
                            successFlag = PHYSICS_UPDATE_ERROR
                            write(*,*) "ERROR: postprocess physics update failed with code=", postprocess_error
                            return
                        end if

                        ! Store this parcel's coldens for the next inner parcel to use
                        coldens_history(dtime_local) = coldens(dstep)

                        call sublimation(abund, points)
                    end if

                    if (returnArray) then
                        call output(returnArray, returnRateConstants, successflag, physicsarray, &
                            chemicalabunarray, rateConstantsArray, heatarray, statsarray, &
                            levelpopulationsarray, sestatsarray, dtime_local, timepoints)
                    else
                        call output(returnArray, returnRateConstants, successflag)
                    end if

                    if (parcelDone) then
                        exit time_loop
                    end if

                end do time_loop  ! inner time loop

            end do parcel_loop  ! outer parcel loop

        else
            ! -- non-RT path: original (time outer, points inner) loop ----------------
            time_loop_non_rt: do while ((successFlag == 0) .and. (timeInYears < finalTime))
                dtime = dtime + 1
                currentTimeold=currentTime
                call updateTargetTime
                coolant_levpop_force_recompute = .true.
                if ((.not. endAtFinalDensity) .and. (targetTime/SECONDS_PER_YEAR > finalTime)) then
                    exit time_loop_non_rt
                end if

                ! Exit loop if density exceeds finalDens (when using density-based stopping)
                if ((parcelStoppingMode /= 0) .and. (density(1) >= finalDens)) then
                    exit time_loop_non_rt
                end if
                !loop over parcels, counting from center out to edge of cloud
                do dstep=1,points
                    !reset time if this isn't first depth point
                    currentTime=currentTimeold

                    if (writeTimestepInfo) then
                        write(*,"(A,ES10.2,A,ES10.2,A,F9.1,A,F9.1,A,ES10.2,A,ES10.2,A,ES10.2)") &
                            "t:", currentTimeold/SECONDS_PER_YEAR, &
                            " dT:", dvode_rstats(11)/SECONDS_PER_YEAR, &
                            " Tg:", gasTemp(1), &
                            " Td:", dustTemp(1), &
                            " nH:", density(1), &
                            " G0:", radfield, &
                            " z:", zeta
                        flush(6)
                        flush_rc = c_fflush(C_NULL_PTR)
                    end if
                    call updateChemistry(successFlag, statsarray, timePoints+1, dtime)
                    if (successFlag < 0) then
                        write(*,*) "Error updating chemistry"
                        return
                    end if
                    if (coolant_error_flag /= 0) then
                        write(*,*) "Coolant error: ", TRIM(coolant_error_message)
                        successFlag = coolant_error_flag
                        coolant_error_flag = 0
                        return
                    end if
                    timeInYears= targetTime/SECONDS_PER_YEAR

                    call coreUpdatePhysics
                    call modelUpdatePhysics
                    if (postprocess_error /= 0) then
                        successFlag = PHYSICS_UPDATE_ERROR
                        write(*,*) "ERROR: postprocess physics update failed with code=", postprocess_error
                        return
                    end if
                    call sublimation(abund, points)
                    if (returnArray) then
                        call output(returnArray, returnRateConstants, successFlag, physicsarray, &
                            & chemicalabunarray, rateConstantsArray, heatarray, statsarray, levelpopulationsarray, &
                            & sestatsarray, dtime, timepoints)
                    else
                        call output(returnArray, returnRateConstants, successFlag)
                    end if
                end do
            end do time_loop_non_rt
        end if
        if (.NOT. returnArray) then
            call finalOutput
            call closeFiles
        end if
        ! IF (ALLOCATED(outIndx)) DEALLOCATE(outIndx)
        ! IF (ALLOCATED(outSpecies)) DEALLOCATE(outSpecies)
    end subroutine solveAbundances

    ! allow(stat-without-message)
    subroutine dictionaryParser(dictionary, outSpeciesIn, successFlag)
        !Reads the input parameters from a string containing a python dictionary/JSON format
        !set of parameter names and values.
        !dictionary - lowercase keys matching the names of the parameters in the select case below
        !OutSpeciesIn - string containing the species to output
        !successFlag - integer flag to indicate success

        !f2py integer, intent(aux) :: nspec
        character(LEN=*), intent(inout) :: dictionary
        character(LEN=*), intent(inout) :: outSpeciesIn
        integer, intent(out) :: successFlag
        integer, allocatable, dimension(:) :: locations
        logical :: ChemicalDuplicateCheck, isOpen
        integer :: posStart, posEnd, whileInteger,inputindx
        character(LEN=100) :: inputParameter, inputValue
        character(256) :: fullPath

        !always deallocate these so that if user didn't specify them,
        ! they don't remain from previous run
        if (ALLOCATED(outIndx)) deallocate(outIndx)
        if (ALLOCATED(outSpecies)) deallocate(outSpecies)

        !All reads use IOSTAT which will change successFlag from 0 if an error occurs
        !so set zero and check for non-zero each loop.
        successFlag=0

        whileInteger = 0

        posStart = scan(dictionary, "{")
        do while (whileInteger /= 1)
            posEnd = scan(dictionary, ":")
            inputParameter = dictionary(posStart+2:posEnd-2)
            dictionary = dictionary(posEnd:)
            posStart = scan(dictionary, " ")
            if (scan(dictionary, ",") == 0) then
                posEnd = scan(dictionary, "}")
                whileInteger = 1
            else
                posEnd = scan(dictionary, ",")
            end if
            inputValue = dictionary(posStart+1:posEnd-1)

            select case (inputParameter)
                case("alpha")
                    !To provide alphas, set keyword alpha in inputdictionary with a dictionary value
                    !that dictionary should be index:value pairs for the alpha array
                    posStart=scan(dictionary,"{")
                    posEnd=scan(dictionary,"}")
                    call coefficientParser(dictionary(posStart+1:posEnd),alpha)
                case("beta")
                    !To provide alphas, set keyword alpha in inputdictionary with a dictionary value
                    !that dictionary should be index:value pairs for the alpha array
                    posStart=scan(dictionary,"{")
                    posEnd=scan(dictionary,"}")
                    call coefficientParser(dictionary(posStart+1:posEnd),beta)
                case("gamma")
                    !To provide alphas, set keyword alpha in inputdictionary with a dictionary value
                    !that dictionary should be index:value pairs for the alpha array
                    posStart=scan(dictionary,"{")
                    posEnd=scan(dictionary,"}")
                    call coefficientParser(dictionary(posStart+1:posEnd),gama)
                case("initialtemp")
                    read(inputValue,*,iostat=successFlag) initialTemp
                case("initialdens")
                    read(inputValue,*,iostat=successFlag) initialDens
                case("finaldens")
                    read(inputValue,*,iostat=successFlag) finalDens
                case("currenttime")
                    read(inputValue,*,iostat=successFlag) currentTime
                case("finaltime")
                    read(inputValue,*,iostat=successFlag) finalTime
                case("enable_radiative_transfer")
                    read(inputValue,*,iostat=successFlag) enable_radiative_transfer
                case("radfield")
                    read(inputValue,*,iostat=successFlag) radfield
                case("zeta")
                    read(inputValue,*,iostat=successFlag) zeta
                case("freezefactor")
                    read(inputValue,*,iostat=successFlag) freezeFactor
                case("rout")
                    read(inputValue,*,iostat=successFlag) rout
                case("rin")
                    read(inputValue,*,iostat=successFlag) rin
                case("baseav")
                    read(inputValue,*,iostat=successFlag) baseAv
                case("points")
                    read(inputValue,*,iostat=successFlag) points
                case("bm0")
                    read(inputValue,*,iostat=successFlag) bm0
                case("density_scale_radius")
                    read(inputValue,*,iostat=successFlag) density_scale_radius
                case("density_power_index")
                    read(inputValue,*,iostat=successFlag) density_power_index
                case("lum_star")
                    read(inputValue,*,iostat=successFlag) lum_star
                case("temp_star")
                    read(inputValue,*,iostat=successFlag) temp_star
                case("parcelstoppingmode")
                    read(inputValue,*,iostat=successFlag) parcelStoppingMode
                case("endatfinaldensity")
                    read(inputValue,*,iostat=successFlag) endAtFinalDensity
                case("freefall")
                    read(inputValue,*,iostat=successFlag) freefall
                case("freefallfactor")
                    read(inputValue,*,iostat=successFlag) freefallFactor
                case("desorb")
                    read(inputValue,*,iostat=successFlag) desorb
                case("h2desorb")
                    read(inputValue,*,iostat=successFlag) h2desorb
                case("crdesorb")
                    read(inputValue,*,iostat=successFlag) crdesorb
                case("uvdesorb")
                    read(inputValue,*,iostat=successFlag) uvdesorb
                case("chemdesorb")
                    read(inputValue,*,iostat=successFlag) chemdesorb
                case("thermdesorb")
                    read(inputValue,*,iostat=successFlag) thermdesorb
                case("instantsublimation")
                    read(inputValue,*,iostat=successFlag) instantSublimation
                case("cosmicrayattenuation")
                    read(inputValue,*,iostat=successFlag) cosmicRayAttenuation
                case("ionmodel")
                    read(inputValue,*,iostat=successFlag) ionModel
                case("improvedh2crpdissociation")
                    read(inputValue,*,iostat=successFlag) improvedH2CRPDissociation
                case("ion")
                    read(inputValue,*,iostat=successFlag) ion
                case("fhe")
                    read(inputValue,*,iostat=successFlag) fhe
                case("fc")
                    read(inputValue,*,iostat=successFlag) fc
                case("fo")
                    read(inputValue,*,iostat=successFlag) fo
                case("fn")
                    read(inputValue,*,iostat=successFlag) fn
                case("fs")
                    read(inputValue,*,iostat=successFlag) fs
                case("fmg")
                    read(inputValue,*,iostat=successFlag) fmg
                case("fsi")
                    read(inputValue,*,iostat=successFlag) fsi
                case("fcl")
                    read(inputValue,*,iostat=successFlag) fcl
                case("fp")
                    read(inputValue,*,iostat=successFlag) fp
                case("ffe")
                    read(inputValue,*,iostat=successFlag) ffe
                case("ff")
                    read(inputValue,*,iostat=successFlag) ff
                case("fd")
                    read(inputValue,*,iostat=successFlag) fd
                case("fli")
                    read(inputValue,*,iostat=successFlag) fli
                case("fna")
                    read(inputValue,*,iostat=successFlag) fna
                case("fpah")
                    read(inputValue,*,iostat=successFlag) fpah
                case("f15n")
                    read(inputValue,*,iostat=successFlag) f15n
                case("f13c")
                    read(inputValue,*,iostat=successFlag) f13c
                case("f18o")
                    read(inputValue,*,iostat=successFlag) f18o
                case("outspecies")
                    ! allow(unchecked-stat)
                    read(inputValue,*,iostat=successFlag) nout
                    allocate(outIndx(nout))
                    allocate(outSpecies(nout))
                    if (outSpeciesIn == "") then
                        write(*,*) "Outspecies parameter set but no outspecies string given"
                        write(*,*) "general(parameter_dict,outSpeciesIn) requires a delimited string of species names"
                        write(*,*) "if outSpecies or columnFlag is set in the parameter dictionary"
                        successFlag=-1
                        return
                    end if

                                            read(outSpeciesIn,*, END=22) outSpecies
                                            if (outSpeciesIn == "") then
                    22                              write(*,*) "mismatch between outSpeciesIn and number given in dictionary"
                                                write(*,*) "Number:",nout
                                                write(*,*) "Species list:",outSpeciesIn
                                                successFlag=-1
                                                return
                                            end if
                    !assign array indices for important species to the integers used to store them.
                    do i=1,nspec
                        do j=1,nout
                            if (specname(i)==outSpecies(j)) outIndx(j)=i
                        end do
                    end do
                case("writestep")
                    read(inputValue,*,iostat=successFlag) writeStep
                case("writetimestepinfo","writetimestep")
                    read(inputValue,*,iostat=successFlag) writeTimestepInfo
                case("heatingFlag", "heatingflag")
                    read(inputValue,*,iostat=successFlag) heatingFlag
                case("ebmaxh2")
                    read(inputValue,*,iostat=successFlag) ebmaxh2
                case("epsilon")
                    read(inputValue,*,iostat=successFlag) epsilon
                case("uvcreff")
                    read(inputValue,*,iostat=successFlag) uvcreff
                case("ebmaxcr")
                    read(inputValue,*,iostat=successFlag) ebmaxcr
                case("phi")
                    read(inputValue,*,iostat=successFlag) phi
                case("ebmaxuvcr")
                    read(inputValue,*,iostat=successFlag) ebmaxuvcr
                case("uv_yield")
                    read(inputValue,*,iostat=successFlag) uv_yield
                case("metallicity")
                    read(inputValue,*,iostat=successFlag) metallicity
                case("omega")
                    read(inputValue,*,iostat=successFlag) omega
                case("difftobindratio")
                    read(inputValue,*,iostat=successFlag) diffToBindRatio
                case("min_desorption_rate")
                    read(inputValue,*,iostat=successFlag) min_desorption_rate
                case("max_desorption_rate_factor")
                    read(inputValue,*,iostat=successFlag) max_desorption_rate_factor
                case("min_desorption_rate_cap")
                    read(inputValue,*,iostat=successFlag) min_desorption_rate_cap
                case("max_desorption_rate_cap")
                    read(inputValue,*,iostat=successFlag) max_desorption_rate_cap
                case("enforcechargeconservation")
                    read(inputValue,*,iostat=successFlag) enforceChargeConservation
                case("reltol")
                    read(inputValue,*,iostat=successFlag) reltol
                case("abstol_factor")
                    read(inputValue,*,iostat=successFlag) abstol_factor
                case("abstol_min")
                    read(inputValue,*,iostat=successFlag) abstol_min
                case("abstol_ice_factor")
                    read(inputValue,*,iostat=successFlag) abstol_ice_factor
                case("abstol_ice_min")
                    read(inputValue,*,iostat=successFlag) abstol_ice_min
                case("negative_abundance_tol")
                    read(inputValue,*,iostat=successFlag) negative_abundance_tol
                case("runtime_conservation_tolerance")
                    read(inputValue,*,iostat=successFlag) runtime_conservation_tolerance
                case("reltol_phys")
                    read(inputValue,*,iostat=successFlag) reltol_phys
                case("abstol_phys_factor")
                    read(inputValue,*,iostat=successFlag) abstol_phys_factor
                case("abstol_t_min")
                    read(inputValue,*,iostat=successFlag) abstol_T_min
                case("abstol_nh_min")
                    read(inputValue,*,iostat=successFlag) abstol_nH_min
                ! CASE('jacobian')
                !     READ(inputValue,*) jacobian
                case("abundsavefile")
                    read(inputValue,*,iostat=successFlag) abundSaveFile
                    abundSaveFile = TRIM(abundSaveFile)
                    if (LEN(abundSaveFile) > 0) then
                        open(abundSaveID,file=abundSaveFile,status="unknown", action="write")
                    end if
                case("abundloadfile")
                    read(inputValue,*,iostat=successFlag) abundLoadFile
                    abundLoadFile = TRIM(abundLoadFile)
                    if (LEN(abundLoadFile) > 0) then
                        open(abundLoadID,file=abundLoadFile,status="old", action="read")
                    end if
                case("outputfile")
                    ! allow(unchecked-stat)
                    read(inputValue,*,iostat=successFlag) outputFile
                    outputFile = trim(outputFile)
                    fullOutput=.true.
                    if (LEN(outputFile) > 0) then
                        open(outputId,file=outputFile,status="unknown",iostat=successFlag, action="write")
                    end if

                    if (successFlag /= 0) then
                        write(*,*) "An error occurred when opening the output file!"//&
                                        & NEW_LINE("A")//&
                                    &" The failed file was ",outputFile&
                                    &, NEW_LINE("A")//"A common error is that the directory doesn't exist"&
                                    &//NEW_LINE("A")//"************************"
                        successFlag=-1
                        return
                    end if
                case("rateConstantFile")
                    ! allow(unchecked-stat)
                    read(inputValue,*,iostat=successFlag) rateConstantFile
                    rateConstantFile = trim(rateConstantFile)
                    if ((LEN(ratesFile) > 0) .and. (.not. storeRatesComputation))then
                        write(*,*) "Error: a ratesFile was specified as output, but fortran was"//&
                        & NEW_LINE("A")//&
                        &"compiled without the rate storage option, please enable `enable_rates_storage` "//&
                        & NEW_LINE("A")//"in the Makerates compilation process."
                    end if
                    open(rateConstantId,file=rateConstantFile,status="unknown",iostat=successFlag, action="write")
                    if (successFlag /= 0) then
                        write(*,*) "An error occurred when opening the rate file!"//&
                                        & NEW_LINE("A")//&
                                    &" The failed file was ",rateConstantFile&
                                    &, NEW_LINE("A")//"A common error is that the directory doesn't exist"&
                                    &//NEW_LINE("A")//"************************"
                        successFlag=-1
                        return
                    end if
                case("ratesFile")
                    ! allow(unchecked-stat)
                    read(inputValue,*,iostat=successFlag) ratesFile
                    ratesFile = trim(ratesFile)
                    if ((LEN(ratesFile) > 0) .and. (.not. storeRatesComputation))then
                        write(*,*) "Error: a ratesFile was specified as output, but fortran was"//&
                        & NEW_LINE("A")//&
                        &"compiled without the rate storage option, please enable `enable_rates_storage` "//&
                        & NEW_LINE("A")//"in the Makerates compilation process."
                    end if
                    open(ratesId,file=ratesFile,status="unknown",iostat=successFlag, action="write")
                    if (successFlag /= 0) then
                        write(*,*) "An error occurred when opening the rate file!"//&
                                        & NEW_LINE("A")//&
                                    &" The failed file was ",ratesFile&
                                    &, NEW_LINE("A")//"A common error is that the directory doesn't exist"&
                                    &//NEW_LINE("A")//"************************"
                        successFlag=-1
                        return
                    end if
                case("heatingfile")
                    ! allow(unchecked-stat)
                    read(inputValue,*,iostat=successFlag) heatingFile
                    heatingFile = trim(heatingFile)
                    open(heatingId,file=heatingFile,status="unknown",iostat=successFlag, action="read")
                    if (successFlag /= 0) then
                        write(*,*) "An error occurred when opening the heating rate file!"//&
                                        & NEW_LINE("A")//&
                                    &" The failed file was ",heatingFile&
                                    &, NEW_LINE("A")//"A common error is that the directory doesn't exist"&
                                    &//NEW_LINE("A")//"************************"
                        successFlag=-1
                        return
                    end if
                case("columnfile")
                    if (trim(outSpeciesIn) /= "") then
                        columnOutput=.true.

                        read(inputValue,*,iostat=successFlag) columnFile
                        columnFile = trim(columnFile)
                        if (LEN(columnFile) > 0) then
                            open(columnId,file=columnFile,status="unknown", action="write")
                        end if
                    else
                        write(*,*) "Error in output species. No species were given but a column file was given."
                        write(*,*) "columnated output requires output species to be chosen."
                        successFlag=-1
                        return
                     end if
                ! Additional parameters for postprocessing mode
                case("fh")
                   read(inputValue,*,iostat=successFlag) fh
                case("mxstep")
                   read(inputValue,*,iostat=successFlag) MXSTEP
                case("h2encounterdesorption")
                   read(inputValue,*,iostat=successFlag) h2EncounterDesorption
                case("hencounterdesorption")
                   read(inputValue,*,iostat=successFlag) hEncounterDesorption
                case("edendothermicityfactor")
                   read(inputValue,*,iostat=successFlag) EDEndothermicityFactor
                case("h2stickingcoeffbyh2coverage")
                   read(inputValue,*,iostat=successFlag) h2StickingCoeffByh2Coverage
                case("hstickingcoeffbyh2coverage")
                   read(inputValue,*,iostat=successFlag) hStickingCoeffByh2Coverage
                case("hdiffusionbarrier")
                   read(inputValue,*,iostat=successFlag) HdiffusionBarrier
                case("usecustomdiffusionbarriers")
                   read(inputValue,*,iostat=successFlag) useCustomDiffusionBarriers
                case("separatediffanddesorbprefactor")
                   read(inputValue,*,iostat=successFlag) separateDiffAndDesorbPrefactor
                case("usetstprefactors")
                   read(inputValue,*,iostat=successFlag) useTSTprefactors
                case("usecustomprefactors")
                   read(inputValue,*,iostat=successFlag) useCustomPrefactors
                case("useminissaleicechemdesefficiency")
                   read(inputValue,*,iostat=successFlag) useMinissaleIceChemdesEfficiency
                case("freq_rel_tol")
                   read(inputValue,*,iostat=successFlag) freq_rel_tol
                case("pop_rel_tol")
                   read(inputValue,*,iostat=successFlag) pop_rel_tol
                case("maxgraintemp")
                   read(inputValue,*,iostat=successFlag) maxGrainTemp
                case("parameterizeh2form")
                   read(inputValue,*,iostat=successFlag) parameterizeH2Form
                case("heating_temp_abstol")
                   read(inputValue,*,iostat=successFlag) heating_temp_abstol
                case("heating_temp_reltol")
                   read(inputValue,*,iostat=successFlag) heating_temp_reltol
                case("solver_mode")
                   read(inputValue,*,iostat=successFlag) solverMode
                case("log_change_threshold")
                   read(inputValue,*,iostat=successFlag) logChangeThreshold
                ! CASE('trajecfile')
                !    READ(inputValue,*,iostat=successFlag) trajecfile
                case default
                    write(*,*) "Problem with given parameter: '", trim(inputParameter), "'."
                    write(*,*) "This is either not supported yet, or invalid."
                    successFlag=-1
                    return
            end select
            dictionary = dictionary(posEnd:)
            if (SCAN(dictionary,",") == 0) whileInteger=1

            !check for failure
            if (successFlag /= 0) then
                write(*,*) "Error reading ",inputParameter
                write(*,*) "This is usually due to wrong type."
                successFlag=PARAMETER_READ_ERROR
                return
            end if
        end do
    end subroutine dictionaryParser

    subroutine coefficientParser(coeffDictString,coeffArray)
        !Similar to dictionaryParser, it reads a python dictionary
        !however, it's intended to read pairs of reaction indices and coefficient values
        !for the alpha, beta, and gama arrays.
        ! No return value, just modifies the coeffArray
        character(LEN=*), intent(inout) :: coeffDictString
        real(dp), intent(inout), dimension(:) :: coeffArray
        integer :: inputIndx, posStart, posEnd
        character(LEN=100) :: inputValue
        logical :: continue_flag

        continue_flag=.true.
        do while (continue_flag)
            !substring containing integer key
            posStart=1
            posEnd=SCAN(coeffDictString,":")
            !read it into index integer
            read(coeffDictString(posStart:posEnd-1),*) inputindx

            !substring including alpha value for the index.
            posStart=posEnd+1
            posEnd=SCAN(coeffDictString,",")
            !last value will have a } instead of , so grab index and tell loop to finish
            if (posEnd == 0) then
                posEnd=SCAN(coeffDictString,"}")
                continue_flag=.false.
            end if

            !read that substring
            inputValue=coeffDictString(posStart:posEnd-1)
            read(inputValue,*) coeffArray(inputIndx)
            !update string to remove this entry
            coeffDictString=coeffDictString(posEnd+1:)
        end do
    end subroutine coefficientParser

    subroutine get_specname(specname_out)
        !Returns:
        ! specname_out - array of species that are in the chemicalabunarray
        !f2py intent(out) specname_out
        character(LEN=32), intent(out) :: specname_out(nspec)
        specname_out(:nspec) = specName
    end subroutine get_specname

    ! TODO: move this to coolant_module, but coolant_module is not being exposed
    ! to f2py currently. So for now, keep it here.
    !=======================================================================
    !
    !  Wrapper function to get the current coolant restart mode
    !  Accessible from Python via f2py
    !
    !-----------------------------------------------------------------------
    integer function get_coolant_restart_mode_wrap() result(coolant_restart_mode_wrap)

        !f2py intent(out) coolant_restart_mode_wrap
        coolant_restart_mode_wrap = GET_COOLANT_RESTART_MODE()
    end function get_coolant_restart_mode_wrap


    !=======================================================================
    !
    !  Wrapper subroutine to set the coolant restart mode
    !  Accessible from Python via f2py
    !
    !  mode values:
    !    0 = WARM (default): Initialize to LTE on first call, rescale on density change
    !    1 = FORCE_LTE: Always reset to LTE before SE iteration
    !    2 = FORCE_GROUND: Always reset to ground state before SE iteration
    !
    !-----------------------------------------------------------------------
    subroutine set_coolant_restart_mode_wrap(mode)
        integer, intent(in) :: mode
        !f2py intent(in) mode
        call SET_COOLANT_RESTART_MODE(mode)
    end subroutine set_coolant_restart_mode_wrap

end module uclchemwrap
