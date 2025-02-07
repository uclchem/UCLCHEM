!The interface between main fortran code and python.
!wrap.f90 subroutines are all accessible through the python wrap.
! Core algorithm is found in solveAbundances subroutine below, all others call it
MODULE uclchemwrap
    USE physicscore
    USE chemistry
    USE io
    USE constants
    USE F2PY_CONSTANTS
    USE postprocess_mod, ONLY: ntime
    IMPLICIT NONE
CONTAINS
    SUBROUTINE cloud(dictionary, outSpeciesIn,&
            &returnArray,givestartabund,timePoints,gridPoints,physicsarray,chemicalabunarray,&
            &abundanceStart,abundance_out,specname_out,successFlag)
        !Subroutine to call a cloud model, used to interface with python
        ! Loads cloud specific subroutines and send to solveAbundances
        !
        !Args:
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
        
        USE cloud_mod
            
        !f2py integer, intent(aux) :: nspec, n_physics_params
        !f2py intent(out) abundance_out, specname_out 
        CHARACTER(LEN=*), INTENT(IN) :: dictionary, outSpeciesIn
        DOUBLE PRECISION, INTENT(OUT) :: abundance_out(nspec)
        CHARACTER(LEN=32), INTENT(OUT) :: specname_out(nspec)
        INTEGER :: successFlag
        !f2py intent(in) dictionary,outSpeciesIn
        !f2py intent(out) successFlag
        LOGICAL, INTENT(IN) :: returnArray
        !f2py intent(in) returnArray
        LOGICAL, INTENT(IN) :: givestartabund
        !f2py intent(in) givestartabund
        INTEGER, INTENT(IN) :: gridPoints
        !f2py intent(in) gridPoints
        INTEGER, INTENT(IN) :: timePoints
        !f2py intent(in) timePoints

        DOUBLE PRECISION, INTENT(INOUT), OPTIONAL, DIMENSION(timePoints+1, gridPoints, n_physics_params) :: physicsarray
        !f2py intent(in,out) physicsarray
        !f2py depend(timePoints,gridPoints, n_physics_params) physicsarray
        DOUBLE PRECISION, INTENT(INOUT), OPTIONAL, DIMENSION(timePoints+1, gridPoints, nspec) :: chemicalabunarray
        !f2py intent(in,out) chemicalabunarray
        !f2py depend(timePoints,gridPoints,nspec) chemicalabunarray
        DOUBLE PRECISION, OPTIONAL, DIMENSION(nspec) :: abundanceStart
        !f2py intent(in) abundanceStart
        !f2py depend(nspec) abundanceStart

        successFlag=0
        specname_out = specName
        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray,givestartabund,&
                &timePoints,physicsarray,chemicalabunarray,abundanceStart)
        IF ((ALLOCATED(outIndx)) .and. (successFlag .eq. 0)) THEN
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF

    END SUBROUTINE cloud


    SUBROUTINE collapse(collapseIn,collapseFileIn,writeOut,dictionary,outSpeciesIn,&
            &returnArray,givestartabund,timePoints,gridPoints,physicsarray,chemicalabunarray,abundanceStart,&
            &abundance_out,specname_out,successFlag)
        !Subroutine to call a collapse model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! collapseIn - integer indicating which collapse mode to run
        ! collapseFileIn - string indicating file to write collapse data to
        ! writeOut - flag indicating whether to write to collapseFileIn
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

        USE collapse_mod

        !f2py integer,parameter intent(aux) nspec, n_physics_params
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn, collapseFileIn
        DOUBLE PRECISION :: abundance_out(nspec)
        CHARACTER(LEN=32) :: specname_out(nspec)
        INTEGER :: successFlag,collapseIn
        LOGICAL :: writeOut
        !f2py intent(in) collapseIn,dictionary,outSpeciesIn,collapseFileIn,writeOut
        !f2py intent(out) abundance_out,specname_out,successFlag
        LOGICAL, INTENT(IN) :: returnArray
        !f2py intent(in) returnArray
        LOGICAL, INTENT(IN) :: givestartabund
        !f2py intent(in) givestartabund
        INTEGER, INTENT(IN) :: gridPoints
        !f2py intent(in) gridPoints
        INTEGER, INTENT(IN) :: timePoints
        !f2py intent(in) timePoints
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints+1, gridPoints, n_physics_params), OPTIONAL :: physicsarray
        !f2py intent(in out) physicsarray
        !f2py depend(gridPoints) physicsarray
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints+1, gridPoints, nspec), OPTIONAL :: chemicalabunarray
        !f2py intent(in out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        DOUBLE PRECISION, DIMENSION(nspec), OPTIONAL :: abundanceStart
        !f2py  intent(in) abundanceStart
        !f2py  depend(gridPoints) abundanceStart
        successFlag=0
        specname_out(:nspec) = specName
        collapse_mode=collapseIn
        writePhysics = writeOut
        collapseFile = collapseFileIn
        
        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray,givestartabund,&
                &timepoints,physicsarray,chemicalabunarray,abundanceStart)
        
        IF ((ALLOCATED(outIndx)) .and. (successFlag .eq. 0)) THEN 
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF
    END SUBROUTINE collapse

    SUBROUTINE hot_core(temp_indx,max_temp,dictionary,outSpeciesIn,&
            &returnArray,givestartabund,timePoints,gridPoints,physicsarray,chemicalabunarray,&
            &abundanceStart,abundance_out,specname_out,successFlag)
        !Subroutine to call a hot core model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! temp_indx - integer indicating which mass hot core to run - see hotcore.90
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
        USE hotcore

        !f2py integer, parameter intent(aux) nspec, n_physics_params
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(nspec),max_temp
        INTEGER :: temp_indx,successFlag
        CHARACTER(LEN=32) :: specname_out(nspec)
        !f2py intent(in) temp_indx,max_temp,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,specname_out,successFlag
        LOGICAL, INTENT(IN) :: returnArray
        !f2py intent(in) returnArray
        LOGICAL, INTENT(IN) :: givestartabund
        !f2py intent(in) givestartabund
        INTEGER, INTENT(IN) :: gridPoints
        !f2py intent(in) gridPoints
        INTEGER, INTENT(IN) :: timePoints
        !f2py intent(in) timePoints
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints+1, gridPoints, n_physics_params), OPTIONAL :: physicsarray
        !f2py intent(in out) physicsarray
        !f2py depend(gridPoints) physicsarray
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints+1, gridPoints, nspec), OPTIONAL :: chemicalabunarray
        !f2py intent(in out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        DOUBLE PRECISION, DIMENSION(nspec), OPTIONAL :: abundanceStart
        !f2py  intent(in) abundanceStart
        !f2py  depend(gridPoints) abundanceStart
        specname_out(:nspec) = specName
        maxTemp=max_temp
        tempIndx=temp_indx

        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
        &updatePhysics,updateTargetTime,sublimation,returnArray,givestartabund,&
        &timepoints,physicsarray,chemicalabunarray,abundanceStart)
    
        IF ((ALLOCATED(outIndx)) .and. (successFlag .eq. 0)) THEN 
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF 
    END SUBROUTINE hot_core

    SUBROUTINE cshock(shock_vel,timestep_factor,minimum_temperature,dictionary, outSpeciesIn,&
            &returnArray,givestartabund,timePoints,gridPoints,physicsarray,chemicalabunarray,abundanceStart,&
            &abundance_out,dissipation_time,specname_out,successFlag)
        !Subroutine to call a C-shock model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! shock_vel - double precision shock velocity in km/s
        ! timestep_factor - Multiply dissipation time of the shock by this double precision factor to
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
        USE cshock_mod

        !f2py integer, parameter intent(aux) nspec, n_physics_params
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(nspec),shock_vel,timestep_factor
        DOUBLE PRECISION :: minimum_temperature,dissipation_time
        INTEGER :: successFlag
        CHARACTER(LEN=32) :: specname_out(nspec)
        !f2py intent(in) shock_vel,timestep_factor,minimum_temperature,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,dissipation_time,specname_out,successFlag
        LOGICAL, INTENT(IN) :: returnArray
        !f2py intent(in) returnArray
        LOGICAL, INTENT(IN) :: givestartabund
        !f2py intent(in) givestartabund
        INTEGER, INTENT(IN) :: gridPoints
        !f2py intent(in) gridPoints
        INTEGER, INTENT(IN) :: timePoints
        !f2py intent(in) timePoints
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints+1, gridPoints, n_physics_params), OPTIONAL :: physicsarray
        !f2py intent(in out) physicsarray
        !f2py depend(gridPoints) physicsarray
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints+1, gridPoints, nspec), OPTIONAL :: chemicalabunarray
        !f2py intent(in out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        DOUBLE PRECISION, DIMENSION(nspec), OPTIONAL :: abundanceStart
        !f2py  intent(in) abundanceStart
        !f2py  depend(gridPoints) abundanceStart

        vs=shock_vel
        timestepFactor=timestep_factor
        minimumPostshockTemp=minimum_temperature

        specname_out = specName
        successFlag=0
        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray,givestartabund,&
                &timepoints,physicsarray,chemicalabunarray,abundanceStart)
        
        IF (successFlag .eq. 0) THEN 
            IF (ALLOCATED(outIndx)) abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
            dissipation_time=dissipationTime
        END IF
    END SUBROUTINE cshock

    SUBROUTINE jshock(shock_vel,dictionary,outSpeciesIn,returnArray,givestartabund,&
            &timePoints,gridPoints,physicsarray,chemicalabunarray,&
            &abundanceStart,abundance_out,specname_out,successFlag)
        !Subroutine to call a J-shock model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! shock_vel - double precision shock velocity in km/s
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
        USE jshock_mod

        !f2py integer, parameter intent(aux) nspec, n_physics_params
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(nspec),shock_vel
        INTEGER :: successFlag
        CHARACTER(LEN=32) :: specname_out(nspec)
        !f2py intent(in) shock_vel,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out, specname_out, successFlag
        LOGICAL, INTENT(IN) :: returnArray
        !f2py intent(in) returnArray
        LOGICAL, INTENT(IN) :: givestartabund
        !f2py intent(in) givestartabund
        INTEGER, INTENT(IN) :: gridPoints
        !f2py intent(in) gridPoints
        INTEGER, INTENT(IN) :: timePoints
        !f2py intent(in) timePoints
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints + 1, gridPoints, n_physics_params), OPTIONAL :: physicsarray
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints + 1, gridPoints, nspec), OPTIONAL :: chemicalabunarray
        !f2py intent(in out) physicsarray
        !f2py depend(gridPoints) physicsarray
        !f2py intent(in out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        DOUBLE PRECISION, DIMENSION(nspec), OPTIONAL :: abundanceStart
        !f2py  intent(in) abundanceStart
        !f2py  depend(gridPoints) abundanceStart
        vs=shock_vel

     
        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
        &updatePhysics,updateTargetTime,sublimation,returnArray,givestartabund,&
        &timepoints,physicsarray,chemicalabunarray,abundanceStart)

        specname_out(:nspec) = specName
        IF ((ALLOCATED(outIndx)) .and. (successFlag .eq. 0)) THEN 
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF 
    END SUBROUTINE jshock

    SUBROUTINE postprocess(dictionary,outSpeciesIn,returnArray,givestartabund,&
        &timePoints,gridPoints,physicsarray,chemicalabunarray,&
        &abundanceStart,timegrid,densgrid,gastempgrid,dusttempgrid,radfieldgrid,zetagrid,&
        &usecoldens,nhgrid,nh2grid,ncogrid,ncgrid,abundance_out,specname_out,successFlag)
        !Subroutine to call a J-shock model, used to interface with python
        ! Loads model specific subroutines and send to solveAbundances
        !
        !Args:
        ! shock_vel - double precision shock velocity in km/s
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
        USE postprocess_mod  

        !f2py integer, parameter intent(aux) nspec, n_physics_params
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(nspec)
        INTEGER :: successFlag
        CHARACTER(LEN=32) :: specname_out(nspec)
        !f2py intent(in) shock_vel,dictionary,outSpeciesIn
        !f2py intent(out) abundance_out, specname_out, successFlag
        LOGICAL, INTENT(IN) :: returnArray
        !f2py intent(in) returnArray
        LOGICAL, INTENT(IN) :: givestartabund
        !f2py intent(in) givestartabund
        INTEGER, INTENT(IN) :: gridPoints
        !f2py intent(in) gridPoints
        INTEGER, INTENT(IN) :: timePoints
        !f2py intent(in) timePoints
        LOGICAL, INTENT(IN) :: usecoldens
        !f2py intent(in) usecoldens
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints+1, gridPoints, n_physics_params), OPTIONAL :: physicsarray
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(timePoints+1, gridPoints, nspec), OPTIONAL :: chemicalabunarray
        !f2py intent(in out) physicsarray
        !f2py depend(gridPoints) physicsarray
        !f2py intent(in out) chemicalabunarray
        !f2py depend(gridPoints) chemicalabunarray
        DOUBLE PRECISION, DIMENSION(nspec), OPTIONAL :: abundanceStart
        !f2py  intent(in) abundanceStart
        !f2py  depend(gridPoints) abundanceStart
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: timegrid
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: densgrid
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: gastempgrid
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: dusttempgrid
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: radfieldgrid
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: zetagrid
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints), OPTIONAL :: nhgrid
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints), OPTIONAL :: nh2grid
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints), OPTIONAL :: ncogrid
        DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints), OPTIONAL :: ncgrid
        !f2py  intent(in) timegrid,densgrid,gastempgrid,dusttempgrid,radfieldgrid,zetagrid
        !f2py  intent(in) nhgrid,nh2grid,ncogrid,ncgrid

        successFlag=0
        if (usecoldens) then
            call solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray,givestartabund,&
                &timepoints,physicsarray,chemicalabunarray,abundanceStart,&
                &timegrid,densgrid,gastempgrid,dusttempgrid,radfieldgrid,zetagrid,&
                &usecoldens,nhgrid,nh2grid,ncogrid,ncgrid)
        else
            call solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,&
                &updatePhysics,updateTargetTime,sublimation,returnArray,givestartabund,&
                &timepoints,physicsarray,chemicalabunarray,abundanceStart,&
                &timegrid,densgrid,gastempgrid,dusttempgrid,radfieldgrid,zetagrid,&
                &usecoldens)
        end if
        
        IF ((ALLOCATED(outIndx)) .and. (successFlag .eq. 0)) THEN 
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF 
    END SUBROUTINE postprocess


    SUBROUTINE get_rates(dictionary,abundancesIn,speciesIndx,rateIndxs,&
        &speciesRates,successFlag,transfer,swap,bulk_layers)
        !Given a species of interest, some parameters and abundances, this subroutine
        !returns the rate of all reactions that include that species plus some extra variables
        !to allow for the calculation of the rate of bulk/surface ice transfer.
        USE cloud_mod
        ! USE constants, only : nspec
        !f2py integer, intent(aux) :: nspec
        CHARACTER(LEN=*):: dictionary
        DOUBLE PRECISION :: abundancesIn(nspec),speciesRates(nspec)
        DOUBLE PRECISION :: transfer,swap,bulk_layers
        INTEGER:: rateIndxs(nspec),speciesIndx, successFlag
        DOUBLE PRECISION :: ydot(nspec+1)
        INTEGER :: speci,bulk_version,surface_version
        !f2py intent(in) dictionary,abundancesIn,speciesIndx,rateIndxs
        !f2py intent(out) speciesRates,successFlag,transfer,swap,bulk_layers
        INCLUDE 'defaultparameters.f90'

        CALL dictionaryParser(dictionary, "",successFlag)
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
        successFlag=0
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
            DO speci=1,nspec
                IF (specname(speci) .eq. "@"//specname(speciesIndx)(2:)) bulk_version=speci
                IF (specname(speci) .eq. "#"//specname(speciesIndx)(2:)) surface_version=speci
            END DO
            IF (SURFGROWTHUNCORRECTED .lt. 0) THEN
                surfaceCoverage = MIN(1.0, safeBulk/safeMantle)
                transfer=SURFGROWTHUNCORRECTED*surfaceCoverage*abund(bulk_version,1)/safeBulk
            ELSE
                surfaceCoverage = bulkGainFromMantleBuildUp()
                transfer=SURFGROWTHUNCORRECTED*surfaceCoverage*abund(surface_version,1)
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
        !Obtain the ODE values for some given parameters and abundances.
        !Essentially runs one time step of solveAbundances  then calls the ODE subroutine (F)
        USE cloud_mod
        use f2py_constants
        !f2py integer, intent(aux) :: nspec
        CHARACTER(LEN=*) :: dictionary
        DOUBLE PRECISION :: abundancesIn(nspec),ratesOut(nspec+1)
        INTEGER :: successFlag
        !f2py intent(in) :: dictionary, abundancesIn
        !f2py intent(out) :: ratesOut
        INCLUDE 'defaultparameters.f90'
        CALL dictionaryParser(dictionary, "",successFlag)

        call coreInitializePhysics(successFlag)
        CALL initializePhysics(successFlag)
        IF (successFlag .lt. 0) then
            WRITE(*,*) 'Error initializing physics'
            RETURN
        END IF

        CALL initializeChemistry(readAbunds)
        dstep=1
        successFlag=0
        abund(:nspec,dstep)=abundancesIn(:nspec)
        abund(neq,dstep)=initialDens
        currentTime=0.0
        timeInYears=0.0
        targetTime=1.0d-7
        CALL updateChemistry(successFlag)
        CALL F(NEQ,currentTime,abund(:,dstep),ratesOut(:NEQ))
    END SUBROUTINE get_odes

    SUBROUTINE solveAbundances(dictionary,outSpeciesIn,successFlag,modelInitializePhysics,&
            &modelUpdatePhysics,updateTargetTime, sublimation, returnArray, givestartabund,&
            &timePoints, physicsarray, chemicalabunarray, abundanceStart,&
            &timegrid,densgrid,gtempgrid,dtempgrid,radgrid,zetagrid,usecoldens,nhgrid,nh2grid,ncogrid,ncgrid)
        ! Core UCLCHEM routine. Solves the chemical equations for a given set of parameters through time
        ! for a specified physical model.
        ! Change behaviour of physics by sending different subroutine arguments - hence the need for model subroutines above
        ! dictionary - the parameter dictionary string reprenting a python dictionary
        ! outSpeciesIn - the species to output to seperate file
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
        !f2py integer, intent(aux) :: nspec
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        EXTERNAL modelInitializePhysics,updateTargetTime,modelUpdatePhysics,sublimation
        INTEGER, INTENT(OUT) :: successFlag
        LOGICAL :: returnArray, givestartabund
        INTEGER :: dtime, timePoints
        ! Arrays needed to work return physics in memory mode
        DOUBLE PRECISION, DIMENSION(:, :, :), OPTIONAL :: physicsarray
        DOUBLE PRECISION, DIMENSION(:, :, :), OPTIONAL :: chemicalabunarray
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: abundanceStart
        ! Arrays neede to work with custom density/temperature profiles
        !  &timegrid,densgrid,gastempgrid,dusttempgrid,nhgrid,nh2grid,ncogrodi,ncgrid)
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: timegrid
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: densgrid
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: gtempgrid
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: dtempgrid
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: radgrid
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: zetagrid
        LOGICAL, OPTIONAL :: usecoldens
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: nhgrid
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: nh2grid
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: ncogrid
        DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: ncgrid
        successFlag=0
        ! Set variables to default values
        INCLUDE 'defaultparameters.f90'
        !Read input parameters from the dictionary
        CALL dictionaryParser(dictionary, outSpeciesIn, successFlag)
        IF (successFlag .lt. 0) THEN
            successFlag=PARAMETER_READ_ERROR
            WRITE(*,*) 'Error reading parameter dictionary'
            RETURN
        END IF
        dstep=1
        currentTime=0.0
        timeInYears=0.0
        IF (returnArray) THEN
            ! Fix to make sure that running in memory mode after running in file mode works correctly
            readAbunds=.False.
            writeAbunds=.False.
        ELSE
            CALL fileSetup
        END IF
        !Initialize core physics first then model specific
        !This allows model to overrule changes made by core
        CALL coreInitializePhysics(successFlag)
        if (present(timegrid)) then
            if (usecoldens) then
                call modelInitializePhysics(successflag, timegrid,densgrid,radgrid,zetagrid,&
                &gtempgrid,dtempgrid,usecoldens,timepoints,nhgrid,nh2grid,ncogrid,ncgrid)
            else 
                call modelInitializePhysics(successflag, timegrid,densgrid,radgrid,zetagrid,&
                &gtempgrid,dtempgrid,usecoldens,timepoints) 
            end if
        else 
            call modelInitializePhysics(successFlag)
        end if 
        IF (successFlag .lt. 0) then
            successFlag=PHYSICS_INIT_ERROR
            WRITE(*,*) 'Error initializing physics'
            RETURN
        END IF
        
        ! Initialize the chemistry
        CALL initializeChemistry(readAbunds)
        IF (returnArray .AND. givestartabund) THEN
            ! In case we have custom abundances, set them here
            DO l=1,points
                abund(:nspec+1,l) = abundanceStart(:nspec+1)
            END DO
        ELSE
            ! Else just use the default readInputAbunds routine:
            CALL readInputAbunds !this won't do anything if no abundLoadFile was in input
        END IF
        !CALL simpleDebug("Initialized")

        dstep = 1
        dtime = 1
        IF (returnArray) THEN
            CALL output(returnArray, successflag, physicsarray, chemicalabunarray, dtime, timepoints)
        ELSE
            CALL output(returnArray, successflag)
        END IF
        

        !loop until the end condition of the model is reached
        DO WHILE ((successFlag .eq. 0) .and. (((endAtFinalDensity) .and. &
            &(density(1) < finalDens)) .or. &
            &((.not. endAtFinalDensity) .and. (timeInYears < finalTime))))
            dtime = dtime + 1
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
                call coreUpdatePhysics
                call modelUpdatePhysics()
                !Sublimation checks if Sublimation should happen this time step and does it
                CALL sublimation(abund)
                !write this depth step now time, chemistry and physics are consistent
                IF (returnArray) THEN
                    CALL output(returnArray, successFlag, physicsarray, chemicalabunarray, dtime, timepoints)
                ELSE
                    CALL output(returnArray, successFlag)
                END IF
            END DO
        END DO
        IF (.NOT. returnArray) THEN
            CALL finalOutput
            CALL closeFiles
        END IF
        ! IF (ALLOCATED(outIndx)) DEALLOCATE(outIndx)
        ! IF (ALLOCATED(outSpecies)) DEALLOCATE(outSpecies)
    END SUBROUTINE solveAbundances

    SUBROUTINE dictionaryParser(dictionary, outSpeciesIn,successFlag)
        !Reads the input parameters from a string containing a python dictionary/JSON format
        !set of parameter names and values.
        !dictionary - lowercase keys matching the names of the parameters in the select case below
        !OutSpeciesIn - string containing the species to output
        !successFlag - integer flag to indicate success
        
        !f2py integer, intent(aux) :: nspec
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
                CASE('bm0')
                    READ(inputValue,*,iostat=successFlag) bm0
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
                    IF (successFlag .ne. 0) THEN
                        write(*,*) "An error occured when opening the output file!"//&
                                        & NEW_LINE('A')//&
                                    &" The failed file was ",outputFile&
                                    &, NEW_LINE('A')//"A common error is that the directory doesn't exist"&
                                    &//NEW_LINE('A')//"************************"
                        successFlag=-1
                        RETURN
                    END IF
                CASE('ratefile')
                    READ(inputValue,*,iostat=successFlag) rateFile
                    rateFile = trim(rateFile)
                    rateOutput=.True.
                    open(rateId,file=rateFile,status='unknown',iostat=successFlag)
                    IF (successFlag .ne. 0) THEN
                        write(*,*) "An error occured when opening the rate file!"//&
                                        & NEW_LINE('A')//&
                                    &" The failed file was ",rateFile&
                                    &, NEW_LINE('A')//"A common error is that the directory doesn't exist"&
                                    &//NEW_LINE('A')//"************************"
                        successFlag=-1
                        RETURN
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
                ! Additional parameters for postprocessing mode
                CASE('fh')
                   READ(inputValue,*,iostat=successFlag) fh
                CASE('ntime')
                   READ(inputValue,*,iostat=successFlag) ntime
                ! CASE('trajecfile')
                !    READ(inputValue,*,iostat=successFlag) trajecfile
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
    END SUBROUTINE dictionaryParser

    SUBROUTINE coefficientParser(coeffDictString,coeffArray)
        !Similar to dictionaryParser, it reads a python dictionary
        !however, it's intended to read pairs of reaction indices and coefficient values
        !for the alpha, beta, and gama arrays.
        ! No return value, just modifies the coeffArray
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
