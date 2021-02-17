DOUBLE PRECISION FUNCTION atomicCooling(gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund)
        use, intrinsic :: iso_fortran_env, dp=>real64

        DOUBLE PRECISION, INTENT(IN) :: gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund
        DOUBLE PRECISION :: t5,invT,rootT,collTFactor !temp/10^5, 1/T and a weird factor from the table
        DOUBLE PRECISION :: hDens,elecDens,heDens,hxDens,hexDens,gauntFactor
        hDens=gasDensity*hAbund
        elecDens=gasDensity*electronAbund
        heDens=gasDensity*heAbund
        hxDens=gasDensity*hxAbund
        hexDens=gasDensity*hexAbund
        t5=1.0d-5*gasT
        invT=1.0/gasT
        rootT=SQRT(gasT)
        collTFactor=1.0/(1.0+SQRT(t5))

        !gauntFactor from Neal et al. 1995
        gauntFactor=1.1+(0.34*EXP(-((5.5-LOG10(gasT))**2.0)/3.0))
        !Neal et al. 1995 lists several fits to cooling each in ergs/cm3/s so we'll just sum them
        !see table 1 of that paper
        !These are just numerical fits so there's loads of magic numbers
        !I've shorted variable names to make it easier to write/read (tn is temperature/10^n)

        !collisional excitation and ionization
        atomicCooling=(7.5d-19*collTFactor*EXP(-118348.0*invT)*elecDens*hDens) &
            &+(5.54d-17*(gasT**-0.397)*collTFactor*EXP(-473638.0*invT)*elecDens*hexDens)&
            &+(1.27d-21*rootT*EXP(-157809.1*invT)*elecDens*hDens*collTFactor)&
            &+(9.38d-22*rootT*EXP(-285335.4*invT)*elecDens*heDens*collTFactor)&
            &+(4.95d-22*rootT*EXP(-631515.0*invT)*elecDens*hexDens*collTFactor)&

        !recombination
            &+(8.7d-27*rootT*((1.0d-3*gasT)**-0.2)*elecDens*hxDens/(1.0+((0.1*t5)**0.7)))&
            &+(1.55d-26*(gasT**0.3647)*elecDens*hexDens)&
            !&+(3.48d-26*rootT*((0.001*gasT)**-0.2)*nelec*nhexI/(1+(0.1*t5)**0.7))&
            &+(1.24d-13*(gasT**-1.5)*EXP(-470000.0*invT)*(1.0+0.3*EXP(-94000.0*invT))*elecDens*hexDens)&

        !free-free emission from all interactions
            &+(1.42d-27*rootT*nelec*(nhex+nhx)*gauntFactor)

    END FUNCTION atomicCooling