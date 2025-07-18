'''
    Advanced Physical Modelling
In the previous tutorial, we simply modelled the chemistry of a static cloud for 1 Myr.
This is unlikely to meet everybody's modelling needs and UCLCHEM is capable of modelling much more complex environments
such as prestellar cores and shocks. In this tutorial, we model both a prestellar core and a shock to explore how these models work
and to demonstrate the workflow that the UCLCHEM team normally follow.

'''
import numpy as np
import uclchem
import matplotlib.pyplot as plt
import os
import time

run_p_core = False
run_shockStart = True
run_cshock = False
run_jshock = False


if run_p_core:
    '''
        The Prestellar Core
        Initial Conditions (Phase 1)
    UCLCHEM typically starts with the gas in atomic/ionic form with no molecules. However, this clearly is not appropriate 
    when modelling an object such as a prestellar core. In these objects, the gas is already evolved and there should be molecules
    in the gas phase as well as ice mantles on the dust. To allow for this, one must provide some initial abundances to 
    the model. There are many ways to do this but we typically chose to run a preliminary model to produce our abundances. 
    In many UCLCHEM papers, we refer to the preliminary model as *phase 1* and the science model as *phase 2*. Phase 1 
    simply models a collapsing cloud and phase 2 models the object in question.

    To do this, we will use `uclchem.model.cloud()` to run a model where a cloud of gas collapses from a density of 
    $10^2 cm^{-3}$ to our prestellar core density of $10^6 cm^{-3}$, keeping all other parameters constant. During this collapse, 
    chemistry will occur and we can assume the final abundances of this model will be reasonable starting abundances for 
    the prestellar core.
    '''
    # set a parameter dictionary for cloud collapse model
    param_dict = {
        "endAtFinalDensity": False,  # stop at finalTime
        "freefall": True,  # increase density in freefall
        "initialDens": 1e2,  # starting density
        "finalDens": 1e6,  # final density
        "initialTemp": 10.0,  # temperature of gas
        "finalTime": 6.0e6,  # final time
        "rout": 0.1,  # radius of cloud in pc
        "baseAv": 1.0,  # visual extinction at cloud edge.
        "outputFile": "../examples/test-output/phase1-full.dat",
    }
    if not os.path.exists("../examples/test-output/"):
        os.makedirs("../examples/test-output/")

    cloud = uclchem.trialmodel.Cloud(param_dict=param_dict)
    '''
    With that done, we now have a file containing the final abundances of a cloud of gas after this collapse: 
    `param_dict["abundSaveFile"]` we can pass this to our prestellar core model to use those abundances as our initial abundances.
    
        Running the Science Model (Phase 2)
    We need to change just a few things in `param_dict` to set up the prestellar core model and pass the cloud model as the
    starting point for the Prestellar core model.
    
    We also want to turn off freefall and change how long the model runs for.
    '''
    # change other bits of input to set up phase 2
    param_dict["initialDens"] = 1e6
    param_dict["finalTime"] = 1e6
    param_dict["freefall"] = False

    # freeze out is completely overwhelmed by thermal desorption
    # so turning it off has no effect on abundances but speeds up integrator.
    param_dict["freezeFactor"] = 0.0

    param_dict["abstol_factor"] = 1e-18
    param_dict["reltol"] = 1e-12

    # pop is dangerous, it removes the original key so you can't rerun this cell.
    param_dict["outputFile"] = "../examples/test-output/phase2-full.dat"

    '''
    Note that we've changed made two changes to the parameters here which aren't strictly necessary but can be helpful 
    in certain situations.

    Since the gas temperature increases throughout a prestellar core model, freeze out is much slower than thermal desorption 
    for all but the first few time steps. Turning it off doesn't affect the abundances but will speed up the solution.
    
    We also change abstol and reltol here, largely to demonstrate their use. They control the integrator accuracy and 
    whilst making them smaller does slow down successful runs, it can make runs complete that stall completely otherwise
     or give correct solutions where lower tolerances allow issues like element conservation failure to sneak in. If 
     your code does not complete or element conservation fails, you can change them.
    
        Checking the Result
    With a successful run, we can check the output. We first load the file and check the abundance conservation, then we can plot it up.
    '''
    p_core = uclchem.trialmodel.PrestellarCore(temp_indx=3, max_temperature=300.0, param_dict=param_dict, previous_model=cloud)
    p_core.check_conservation()
    df_p_core = p_core.get_dataframes()

    species = ["CO", "H2O", "CH3OH", "#CO", "#H2O", "#CH3OH", "@H2O", "@CO", "@CH3OH"]
    fig, [ax, ax2] = plt.subplots(1, 2, figsize=(16, 9))
    ax = p_core.plot_species(ax, species)
    settings = ax.set(
        yscale="log",
        xlim=(1e2, 1e6),
        ylim=(1e-10, 1e-2),
        xlabel="Time / years",
        ylabel="Fractional Abundance",
        xscale="log",
    )

    ax2.plot(df_p_core["Time"], df_p_core["Density"], color="black")
    ax2.set(xscale="log")
    ax3 = ax2.twinx()
    ax3.plot(df_p_core["Time"], df_p_core["gasTemp"], color="red")
    ax2.set(xlabel="Time / year", ylabel="Density")
    ax3.set(ylabel="Temperature", facecolor="red", xlim=(1e2, 1e6))
    ax3.tick_params(axis="y", colors="red")
    plt.show()

    '''
    Here, we see the value of running a collapse phase before the science run. Having run a collapse, we start this 
    model with well developed ices and having material in the surface and bulk allows us to properly model the effect 
    of warm up in a prestellar core. For example, the @CO abundance is $\sim10^{-4}$ and #CO is $\sim10^{-6}$. As the gas 
    warms to around 30K, the #CO abundance drops drastically as CO's binding energy is such that it is efficiently 
    desorbed from the surface at this temperature. However, the rest of the CO is trapped in the bulk, surrounded by 
    more strongly bound H2O molecules. Thus, the @CO abundance stays high until the gas reaches around 130K, when the 
    H2O molecules are released along with the entire bulk.
    '''

if run_cshock or run_jshock or run_shockStart:
    '''
        Shocks
    Essentially the same process should be followed for shocks. Let's run a C-type and J-type shock through a gas of
    density. Again, we first run a simple cloud model to obtain some reasonable starting abundances, then we can run 
    the shocks.
    '''

    param_dict = {
        "endAtFinalDensity": False,  # stop at finalTime
        "freefall": True,  # increase density in freefall
        "initialDens": 1e2,  # starting density
        "finalDens": 1e4,  # final density
        "initialTemp": 10.0,  # temperature of gas
        "finalTime": 6.0e6,  # final time
        "rout": 0.1,  # radius of cloud in pc
        "baseAv": 1.0,  # visual extinction at cloud edge.
        "abundSaveFile": "../examples/test-output/shockstart.dat",
        "outputFile": "../examples/test-output/shockstart-full.dat",
    }

    shock_start = uclchem.trialmodel.Cloud(param_dict=param_dict)
    print(len(shock_start.next_starting_chemistry))
    print(shock_start.chemical_abun_array[235,:,:])

    if run_cshock:
        '''
            C-shock
        We'll first run a c-shock. We'll run a 40 km s $^{-1}$ shock through a gas of density $10^4$ cm $^{-3}$, using 
        the abundances we just produced. For this model, the object will have an extra attribute. The dissipation time 
        of the shock. We'll use that time to make our plots look nicer, cutting to a reasonable time. You can obtain it 
        by using `cshock.dissipation_time`.
        
        '''
        # change other bits of input to set up phase 2
        param_dict["initialDens"] = 1e4
        param_dict["finalTime"] = 1e6

        cshock = uclchem.trialmodel.CShock(shock_vel=40,
                                           param_dict=param_dict,
                                           starting_chemistry=shock_start.next_starting_chemistry)
        '''
        The code completes fine. We do get a couple of warnings though. First, we're informed that `freefall` must be 
        set to False for the C-shock model. Then we get a few integrator warnings. These are not important and can be 
        ignored as long as the element conservation looks ok. However, it is an indication that the integrator did 
        struggle with these ODEs under these conditions.
        '''
        df_cshock = cshock.get_dataframes()
        cshock.check_conservation()

        species = ["CO", "H2O", "CH3OH", "NH3", "$CO", "$H2O", "$CH3OH", "$NH3"]

        fig, [ax, ax2] = plt.subplots(1, 2, figsize=(16, 9))
        ax = cshock.plot_species(ax, species)
        settings = ax.set(
            yscale="log",
            xlim=(1, 20 * cshock.dissipation_time),
            ylim=(1e-10, 1e-2),
            xlabel="Time / years",
            ylabel="Fractional Abundance",
            xscale="log",
        )

        ax2.plot(df_cshock["Time"], df_cshock["Density"], color="black")
        ax2.set(xscale="log")
        ax3 = ax2.twinx()
        ax3.plot(df_cshock["Time"], df_cshock["gasTemp"], color="red")
        ax2.set(xlabel="Time / year", ylabel="Density")
        ax3.set(ylabel="Temperature", facecolor="red", xlim=(1, 20 * cshock.dissipation_time))
        ax3.tick_params(axis="y", colors="red")
        plt.show()

    if run_jshock:
        '''
            J-shock
        Running a j-shock is a simple case of changing function. We'll run a 10 km s $^{-1}$ shock through a gas of 
        density $10^3$ cm $^{-3}$ gas this time. Note that nothing stops us using the intial abundances we produced for 
        the c-shock. UCLCHEM will not check that the initial density matches the density of the previous model. It may 
        not always be a good idea to do this but we should remember the intial abundances really are just a rough 
        approximation.
        '''

        param_dict["initialDens"] = 1e3
        param_dict["freefall"] = False  # lets remember to turn it off this time
        param_dict["reltol"] = 1e-12

        shock_vel = 10.0
        jshock = uclchem.trialmodel.JShock(shock_vel=shock_vel,
                                           param_dict=param_dict,
                                           previous_model=shock_start.next_starting_chemistry,
                                           timepoints=1500)

        print(jshock.success_flag)

        df_jshock = jshock.get_dataframes()
        jshock.check_conservation()

        species = ["CO", "H2O", "CH3OH", "NH3", "$CO", "$H2O", "$CH3OH", "$NH3"]

        fig, [ax, ax2] = plt.subplots(1, 2, figsize=(16, 9))
        ax = jshock.plot_species(ax, species)
        settings = ax.set(
            yscale="log",
            xlim=(1e-7, 1e6),
            ylim=(1e-10, 1e-2),
            xlabel="Time / years",
            ylabel="Fractional Abundance",
            xscale="log",
        )

        ax2.plot(df_jshock["Time"], df_jshock["Density"], color="black")
        ax2.set(xscale="log", yscale="log")
        ax3 = ax2.twinx()
        ax3.plot(df_jshock["Time"], df_jshock["gasTemp"], color="red")
        ax2.set(xlabel="Time / year", ylabel="Density")
        ax3.set(ylabel="Temperature", facecolor="red", xlim=(1e-7, 1e6))
        ax3.tick_params(axis="y", colors="red")
        plt.show()

'''
That's everything! We've run various science models using reasonable starting abundances that we produced by running a 
simple UCLCHEM model beforehand. One benefit of this method is that the abundances are consistent with the network. If 
we start with arbitrary, perhaps observationally motivated, abundances, it would be possible to initiate the model in 
a state our network could never produce.

However, one should be aware of the limitations of this method. A freefall collapse from low density to high is not 
really how a molecular cloud forms and so the abundances are only approximately similar to values they'd truly have in 
a real cloud. Testing whether your results are sensitive to things like the time you run the preliminary for or the 
exact density is a good way to make sure these approximations are not problematic.

Bear in mind that you can pass previous models to new models. This lets you chain model runs 
together. For example, you could run a c-shock from a cloud model as we did here and then a j-shock with the c-shock's 
abundances as the initial abundances.
'''