! Chemistry module of UCL_CHEM. Contains all the core machinery of the code, not really intended to be altered.
! Use physics module to alter temp/density behaviour etc. This module should solve chemistry for a cloud of gas
MODULE chem
USE physics
IMPLICIT NONE
EXTERNAL dvode
   !These integers store the array index of important species and reactions, x is for ions    
    integer :: nh,nh2,nc,ncx,no,nn,ns,nhe,nco,nmg,nf,nh2o,nsi,nsix,ncl,nclx,nch3oh,np
    integer :: nrco,nout,nspec,nreac,njunk,evapevents,ngrainco
    integer, allocatable :: outIndx(:)
    !loop counters    
    integer :: i,j,l,writestep

    !These are variables for reaction rates, alpha/beta/gamas are combined each time step to make rate,the total reaction rate
    double precision,allocatable :: rate(:),alpha(:),beta(:),gama(:),mass(:)
    character(LEN=10),allocatable :: re1(:),re2(:),re3(:),p1(:),p2(:),p3(:),p4(:)
    character(LEN=15),allocatable :: outSpecies(:),specname(:)
    
    !DLSODE variables    
    integer :: ITOL,ITASK,ISTATE,IOPT,MESFLG,NEQ,lrw,liw
    integer :: MXSTEP,MF
    integer,allocatable :: IWORK(:)
    double precision :: reltol,rpar,ipar
    double precision, allocatable :: RWORK(:),abstol(:)

    !initial fractional elemental abudances and arrays to store abundances
    double precision :: fh,fhe,fc,fo,fn,fs,fmg,fsi,fcl,fp,ff,h2col,cocol,junk1,junk2
    double precision,allocatable :: abund(:,:),mantle(:)
    
    !Variables controlling chemistry
    double precision :: radfield,zeta,fr,omega,grainArea,cion,h2form,h2dis
    double precision :: ebmaxh2,epsilon,ebmaxcrf,ebmaxcr,phi,ebmaxuvcr,uvy,uvcreff
    double precision :: taud,dopw,radw,xl,fosc

    !evaporation lists, these are used to evaporation species in specific events
    !See viti 2004 for more information.
    integer, allocatable :: colist(:),mcolist(:),intlist(:),mintlist(:),grainlist(:),mgrainlist(:)
    integer, allocatable :: co2list(:),mco2list(:),int2list(:),mint2list(:)
    double precision, allocatable :: cobindener(:),co2bindener(:),intbindener(:)
    double precision, allocatable :: comono(:),covolc(:),co2mono(:),co2volc(:),intmono(:),intvolc(:)
    double precision, allocatable :: bindener(:),vdiff(:)

    !Variables for selfshielding rates for CO and H2
    logical :: startr
    integer :: dimco, dimh2
    double precision :: corates(7,6), y2r(7,6), ncogr(7), nh2gr(6)    

CONTAINS
!This gets called immediately by main so put anything here that you want to happen before the time loop begins, reader is necessary.
    SUBROUTINE chem_initialise
        CALL reader
        !if this is the first step of the first phase, set initial abundances
        !otherwise reader will fix it
        IF (first.eq.1) THEN
            !ensure abund is initially zero
            abund= 0.
      
            !As default, have half in molecular hydrogen and half in atomic hydrogen
            abund(nh2,:) = 0.5*(0.5*(1.0e0-fh))
            abund(nh,:) = (0.5*(1.0e0-fh))     
            abund(nhe,:) = fhe                       
            abund(no,:) = fo  
            abund(nn,:) = fn               
            abund(ns,:) = fs
            abund(nmg,:) = fmg
            abund(nsix,:) = fsi                
            abund(nclx,:) = fcl 
            !abund(np,:) = fp
            !abund(nf,:) = ff

            !abund(nfe,:) = ffe
            !abund(nna,:) = fna
            abund(nspec+1,:)=dens      

            !Decide how much iron is initiall ionized using parameters.f90
            SELECT CASE (ion)
                CASE(0)
                    abund(nc,:)=fc
                    abund(ncx,:)=1.d-10
                CASE(1)
                    abund(nc,:)=fc*0.5
                    abund(ncx,:)=fc*0.5
                CASE(2)
                    abund(nc,:)=1.d-10
                    abund(ncx,:)=fc
            END SELECT
            abund(nspec,:)=abund(ncx,:)

        ENDIF

        DO  i=lbound(mgrainlist,1),ubound(mgrainlist,1)
            j=mgrainlist(i)
            vdiff(j)=2.5e14*bindener(j)/mass(j)
            vdiff(j)=dsqrt(vdiff(j))
        END DO

        !h2 formation rate initially set
        h2form = 1.0d-17*dsqrt(initialTemp)
        allocate(mantle(points))
        DO l=1,points
            mantle(l)=sum(abund(mgrainlist,l))
        END DO
        
        !DVODE SETTINGS
        ISTATE=1;MF=22;ITOL=1;ITASK=1;IOPT=1;MESFLG=1
        reltol=1e-4;MXSTEP=10000

        NEQ=nspec+1
        LIW=30+NEQ
        LRW=22+(9*NEQ)+(2*NEQ*NEQ)
        allocate(IWORK(LIW),RWORK(LRW),abstol(NEQ))


    END SUBROUTINE chem_initialise

!Reads input reaction and species files as well as the final step of previous run if this is phase 2
    SUBROUTINE reader
        IMPLICIT NONE
        integer i,j,l,m

        !read species file and allocate sufficient space to relevant arrays
        read(21,*)nspec
        allocate(abund(nspec+1,points),specname(nspec),mass(nspec),bindener(nspec),vdiff(nspec))
        read(21,*)(specname(j),mass(j),bindener(j),j=1,nspec-1)
        
        nout = SIZE(outSpecies)
        allocate(outIndx(nout))

        !assign array indices for important species to the integers used to store them.
        specname(nspec)='electr'
        DO i=1,nspec-1
            IF (specname(i).eq.'H')   nh  = i
            IF (specname(i).eq.'H2')  nh2 = i
            IF (specname(i).eq.'C')   nc  = i
            IF (specname(i).eq.'C+')  ncx = i
            IF (specname(i).eq.'O')   no  = i
            IF (specname(i).eq.'N')   nn  = i
            IF (specname(i).eq.'S+')  ns  = i
            IF (specname(i).eq.'HE')  nhe = i
            IF (specname(i).eq.'CO')  nco = i
            IF (specname(i).eq.'MG')  nmg = i
            IF (specname(i).eq.'H2O') nh2o = i
            IF (specname(i).eq.'SI')  nsi = i
            IF (specname(i).eq.'SI+') nsix= i
            IF (specname(i).eq.'CL')  ncl = i
            IF (specname(i).eq.'CL+') nclx= i
            IF (specname(i).eq.'CH3OH') nch3oh= i
            IF (specname(i).eq.'#CO') ngrainco = i
            IF (specname(i).eq. 'P') np=i
            IF (specname(i).eq.'F') nf=i
            DO j=1,nout
                IF (specname(i).eq.outSpecies(j)) outIndx(j)=i
            END DO
        END DO

        !read reac file, assign array space
        !alpha, beta and gama are used for working out reaction rate each time step
        read(22,*) nreac
        allocate(re1(nreac),re2(nreac),re3(nreac),p1(nreac),p2(nreac),p3(nreac),&
            &p4(nreac),alpha(nreac),beta(nreac),gama(nreac),rate(nreac))
        DO j=1,nreac
            read(22,*) re1(j),re2(j),re3(j),p1(j),p2(j),p3(j),&
            &p4(j),alpha(j),beta(j),gama(j),junk1,junk2
        END DO    

        !finally, read in evaporation lists
        !what we do is read in length of each list, then the numbers in the list
        !these are used for evaporation sums.
        read(23,*) l
        allocate(colist(l),mcolist(l), cobindener(l),comono(l),covolc(l))
        read(23,*)colist
        read(23,*)mcolist
        read(23,*) comono
        read(23,*) covolc
        read(23,*) l
        allocate(co2list(l),mco2list(l), co2bindener(l),co2mono(l),co2volc(l))
        read(23,*)co2list
        read(23,*)mco2list
        read(23,*) co2mono
        read(23,*) co2volc
        read(23,*) l
        allocate(intlist(l),mintlist(l),intbindener(l),intmono(l),intvolc(l))
        read(23,*) intlist
        read(23,*) mintlist
        read(23,*) intmono
        read(23,*) intvolc
        read(23,*)l
        allocate(grainlist(l),mgrainlist(l))
        read(23,*) grainlist
        read(23,*) mgrainlist

        !read start file IF not first phase to get finale abundances from previous phase 
        !density, temp and av read but NOT zeta or radfield
        IF (first .eq. 0) THEN
            DO l=1,points
                read(7,*)
                read(7,7000) abund(nspec+1,l),temp(l),av(l)
                read(7,*)
                read(7,7010) h2form,fc,fo,&
                            &fmg,fhe,dstep
                read(7,*)
                read(7,7030) (specname(i),abund(i,l),i=1,nspec)
                rewind(7)
            END DO
            7000 format(&
            &33x,1pe11.4,5x,/,&
            &33x,0pf8.2,2x,/,&
            &33x,0pf12.4,4x,/)
            7010 format(&
            &33x,1pe8.2,8x,/,&
            &11x,1pe7.1,4x,12x,1pe7.1,/&
            &12x,1pe7.1,13x,1pe7.1,&
            &13x,i3,/)
            7020  format(//)
            7030  format(4(1x,a15,2x,1pe10.3,:))     
        END IF
    END SUBROUTINE reader

!Writes physical variables and fractional abundances to output file, called every time step.
    SUBROUTINE output
        !write out cloud properties
        write(10,8020) tage,dens(dstep),temp(dstep),av(dstep),radfield,zeta,h2form,fc,fo,&
                        &fmg,fhe,dstep
        !and a blank line
        write(10,8000)
        !and then all the abundances for this step
        write(10,8010) (specname(i),abund(i,dstep),i=1,nspec) 
        write(10,8000)
        !If this is the last time step of phase I, write a start file for phase II
        IF (first .eq. 1) THEN
           IF (switch .eq. 0 .and. tage .ge. finalTime& 
               &.or. switch .eq. 1 .and.dens(dstep) .ge. finalDens) THEN
               write(7,8020) tage,dens(dstep),temp(dstep),av(dstep),radfield,zeta,h2form,fc,fo,&
                       &fmg,fhe,dstep
               write(7,8000)
               write(7,8010) (specname(i),abund(i,dstep),i=1,nspec)
               write(7,8000)
           ENDIF
        ENDIF
        8000  format(/)
        8010  format(4(1x,a15,'=',1x,1pe10.3,:))
        8020 format(&
        &'age of cloud             time  = ',1pe11.3,' years',/,&
        &'total hydrogen density   dens  = ',1pe11.4,' cm-3',/,&
        &'cloud temperature        temp  = ',0pf8.2,' k',/,&
        &'visual extinction        av    = ',0pf12.4,' mags',/,&
        &'radiation field          rad   = ',0pf10.2,' (habing = 1)',/,&
        &'cosmic ray ioniz. rate   zeta  = ',0pf10.2,' (unit = 1.3e-17s-1)',/,&
        &'h2 formation rate coef.        = ',1pe8.2,' cm3 s-1',/,&
        &'c / htot = ',1pe7.1,4x,' o / htot = ',1pe7.1,/&
        &'mg / htot = ',1pe7.1,&
        &' he / htot = ',1pe7.1,&
        &' depth     = ',i3)

        !Every 'writestep' timesteps, write the chosen species out to separate file
        !choose species you're interested in by looking at parameters.f90
        IF ( mod(tstep,writestep) .eq. 0) THEN
            write(11,8030) tage,dens(dstep),temp(dstep),abund(outIndx,dstep)
            8030  format(1pe11.3,1x,1pe11.4,1x,0pf8.2,6(1x,1pe10.3))
            write(*,*)'Call to LSODE successful at time: ',(TOUT*year),' years'
            write(*,*)'        Steps: ',IWORK(6)
        END IF
    END SUBROUTINE output

    SUBROUTINE chem_update
    !Called every time/depth step and updates the abundances of all the species

        !y is at final value of previous depth iteration so set to initial values of this depth with abund
        !reset other variables for good measure        
        h2form = 1.0d-17*dsqrt(temp(dstep))
    
        !Sum of abundaces of all mantle species. mantleindx stores the indices of mantle species.
        mantle(dstep)=sum(abund(mgrainlist,dstep))

        !evaluate co and h2 column densities for use in rate calculations
        !sum column densities of each point up to dstep. boxlength and dens are pulled out of the sum as common factors  
        IF (dstep.gt.1) THEN
            h2col=(sum(abund(nh2,:dstep-1))+0.5*abund(nh2,dstep))*dens(dstep)*(cloudSize/real(points))
            cocol=(sum(abund(nco,:dstep-1))+0.5*abund(nco,dstep))*dens(dstep)*(cloudSize/real(points))
        ELSE
            h2col=0.5*abund(nh2,dstep)*dens(dstep)*(cloudSize/real(points))
            cocol=0.5*abund(nco,dstep)*dens(dstep)*(cloudSize/real(points))
        ENDIF

        !call the actual ODE integrator
        CALL integrate

        !call evaporation to remove species from grains at certain temperatures
        CALL evaporate

        !1.d-30 stops numbers getting too small for fortran.
        WHERE(abund<1.0d-30) abund=1.0d-30
    END SUBROUTINE chem_update

    SUBROUTINE integrate
    !This subroutine calls DVODE (3rd party ODE solver) until it can reach tout with acceptable errors (reltol/abstol)

        DO WHILE(t0 .lt. tout)            
            !reset parameters for DVODE
            ITOL=2 !abstol is an array
            ITASK=1 !try to integrate to tout

            !first step only, set some stuff up
            IF(ISTATE .EQ. 1) THEN
                IOPT=1
                IWORK(6)=MXSTEP
            ENDIF

            abstol=1.0d-16*abund(:,dstep)
            WHERE(abstol<1d-30) abstol=1d-30
            !get reaction rates for this iteration
            CALL reacrates
            !Call the integrator.
            CALL DVODE(F,NEQ,abund(:,dstep),T0,TOUT,ITOL,reltol,abstol,ITASK,ISTATE,IOPT,&
            &             RWORK,LRW,IWORK,LIW,JAC,MF,RPAR,IPAR)

            SELECT CASE(ISTATE)
                CASE(-1)
                    !More steps required for this problem
                    MXSTEP=MXSTEP*2    
                    write(79,*)'Call to LSODE returned -1 meaning that MXSTEP exceeded'
                    write(79,*)'but the integration was successful'
                    write(79,*)'Doubling MXSTEP from:',MXSTEP,' to:',MXSTEP*2
                    ISTATE=3
                    IOPT=1
                    IWORK(6)=MXSTEP
                CASE(-2)
                    !Tolerances are too small for machine but succesful to current t0
                    abstol=abstol*10.0
                    ISTATE=3
                CASE(-3)
                    write(79,*) "DVODE found invalid inputs"
                    write(79,*) "abstol"
                    write(79,*) abstol
                    STOP
                CASE(-4)
                    !Successful as far as t0 but many errors.
                    !Make tout smaller and just go again
                    tout=(t0+tout)/2.0
                    ISTATE=2
                CASE DEFAULT
                    IOPT=0
                    ISTATE=3
            END SELECT
        END DO                   
    END SUBROUTINE integrate

    !This is where reacrates subroutine is hidden
    include 'rates.f90'

    SUBROUTINE  F (NEQ, T, Y, YDOT)
        !DLSODE calls this subroutine to ask it what the RHS of the equations dy/dt=... are    

        INTEGER :: NEQ
        DOUBLE PRECISION :: T,Y(nspec+1),YDOT(nspec+1)
        DOUBLE PRECISION :: D,loss,prod
        
        !For collapse =1 Dens is updated by DLSODE just like abundances so this ensures dens is at correct value
        !For collapse =0 allow option for dens to have been changed elsewhere.
        IF (collapse .ne. 1) THEN
            y(nspec+1)=dens(dstep)
        ELSE
            dens(dstep)=y(nspec+1)
        END IF

        !Set D to the gas density for use in the ODEs
        D=dens(dstep)
        !The ODEs created by MakeRates go here, they are essentially sums of terms that look like k(1,2)*y(1)*y(2)*dens. Each species ODE is made up
        !of the reactions between it and every other species it reacts with.
        INCLUDE 'odes.f90'

        !updated just in case temp changed
        h2form=1.0d-17*dsqrt(temp(dstep)) 

        !H2 formation should occur at both steps - however note that here there is no 
        !temperature dependence. y(nh) is hydrogen fractional abundance.
        ydot(nh)  = ydot(nh) - 2.0*( h2form*dens(dstep)*y(nh) - h2dis*y(nh2) )
        !                             h2 formation - h2-photodissociation
        ydot(nh2) = ydot(nh2) + h2form*dens(dstep)*y(nh) - h2dis*y(nh2)
        !                       h2 formation  - h2-photodissociation

        ! get density change from physics module to send to DLSODE
        IF (collapse .eq. 1) ydot(nspec+1)=densdot()

    END SUBROUTINE F

!integrate calls reacrates to get the reaction rates at every iteration. reacrates calls further functions.
!This file is already long so I've hidden those subroutines in rates.f90
 
    SUBROUTINE evaporate
    !Evaporation is based on Viti et al. 2004. A proportion of the frozen species is released into the gas phase
    !in specific events. These events are activated by flags (eg solidflag) which can be set in physics module.
    !The species evaporated are in lists, created by Makerates and based on groupings. see the viti 2004 paper.
    IF (tstep .gt. 1) THEN
        !Viti 04 evap
        IF (evap .eq. 1) THEN
            !Solid Evap
            IF (solidflag .eq. 1) THEN
                abund(colist,dstep)=abund(colist,dstep)+0.35*abund(mcolist,dstep)
                abund(mcolist,dstep)=0.65*abund(mcolist,dstep)
                !Set flag to 2 to stop it being recalled
                solidflag=2
            ENDIF

            !mono evap
            CALL temper

            !Volcanic evap
            IF (volcflag .eq. 1) THEN
                DO i=lbound(colist,1),ubound(colist,1)
                    abund(colist,dstep)=abund(colist,dstep)+covolc(i)*abund(mcolist,dstep)
                    abund(mcolist,dstep)=(1.0-covolc(i))*abund(mcolist,dstep)
                END DO
                DO i=lbound(co2list,1),ubound(co2list,1)
                    abund(co2list,dstep)=abund(co2list,dstep)+co2volc(i)*abund(mco2list,dstep)
                    abund(mco2list,dstep)=(1.0-co2volc(i))*abund(mco2list,dstep)
                END DO
                DO i=lbound(intlist,1),ubound(intlist,1)
                    abund(intlist,dstep)=abund(intlist,dstep)+intvolc(i)*abund(mintlist,dstep)
                    abund(mintlist,dstep)=(1.0-intvolc(i))*abund(mintlist,dstep)
                END DO
                !abund(int2list,dstep)=abund(int2list,dstep)+0.5*abund(mint2list,dstep)
                !abund(mint2list,dstep)=0.5*abund(mint2list,dstep)
                !Set flag to 2 to stop it being recalled
                volcflag=2
            ENDIF

            !Co-desorption
            IF (coflag .eq. 1) THEN
                abund(grainlist,dstep)=abund(grainlist,dstep)+abund(mgrainlist,dstep)
                abund(mgrainlist,dstep)=1d-30
                coflag=2
            ENDIF
        ELSE IF (evap .eq. 2 .and. coflag .ne. 2) THEN
            !Alternative evap. Instaneous evaporation of all grain species
            abund(grainlist,dstep)=abund(grainlist,dstep)+abund(mgrainlist,dstep)
            abund(mgrainlist,dstep)=1d-30
            coflag = 2
        ENDIF
    ENDIF

    END SUBROUTINE evaporate

    SUBROUTINE temper
        !Subroutine to handle mono-evaporation. See viti 2004
        double precision en,newm,surden,expdust,freq,kevap
        integer speci
        parameter(surden=1.5e15)

        !co monoevap
        DO i=lbound(colist,1),ubound(colist,1)
            speci=colist(i)
            en=bindener(speci)*kbolt
            expdust=bindener(speci)/temp(dstep)
            newm = mass(speci)*1.66053e-27
            freq = dsqrt((2*(surden)*en)/((pi**2)*newm))
            kevap=freq*exp(-expdust)
            IF (kevap .ge. 0.99) THEN
                write(*,*)i
                abund(speci,dstep)=abund(speci,dstep)+(comono(i)*abund(mcolist(i),dstep))
                abund(mcolist(i),dstep)=(1.0-comono(i))*abund(mcolist(i),dstep)
                !set to 1d50 so it can't happen again
                bindener(speci)=1d50
            END IF 
        END DO

        !co2 monoevap
        DO i=lbound(co2list,1),ubound(co2list,1)
            speci=co2list(i)
            en=bindener(speci)*kbolt
            expdust=bindener(speci)/temp(dstep)
            newm = mass(speci)*1.66053e-27
            freq = dsqrt((2*(surden)*en)/((pi**2)*newm))
            kevap=freq*exp(-expdust)
            IF (kevap .ge. 0.99) THEN
                write(*,*)i
                abund(speci,dstep)=abund(speci,dstep)+(co2mono(i)*abund(mco2list(i),dstep))
                abund(mco2list(i),dstep)=(1.0-co2mono(i))*abund(mco2list(i),dstep)
                bindener(speci)=1d50
            END IF 
        END DO

        !int mono evap
        DO i=lbound(intlist,1),ubound(intlist,1)
            speci=intlist(i)
            en=bindener(speci)*kbolt
            expdust=bindener(speci)/temp(dstep)
            newm = mass(speci)*1.66053e-27
            freq = dsqrt((2*(surden)*en)/((pi**2)*newm))
            kevap=freq*exp(-expdust)
            IF (kevap .ge. 0.99) THEN
                write(*,*)i
                abund(speci,dstep)=abund(speci,dstep)+(intmono(i)*abund(mintlist(i),dstep))
                abund(mintlist(i),dstep)=(1.0-intmono(i))*abund(mintlist(i),dstep)
                bindener(speci)=1d50
            END IF 
        END DO
    END SUBROUTINE temper
!This is a dummy for DLSODE, it has to call it but we do not use it.
    subroutine JAC (NEQ,T,Y,ML,MU,PD,NROWPD)
         INTEGER  NEQ, ML, MU, NROWPD
         DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
    END SUBROUTINE JAC        

    SUBROUTINE debugout
        write(79,*) "Integrator failed, printing relevant debugging information"
        write(79,*) "dens",dens(dstep)
        write(79,*) "density in integration array",abund(nspec+1,dstep)
        write(79,*) "Av", av(dstep)
        write(79,*) "Mantle", mantle(dstep)
        write(79,*) "Temp", temp(dstep)
        DO i=1,nreac
            if (rate(i) .ge. huge(i)) write(79,*) "Rate(",i,") is potentially infinite"
        END DO
    END SUBROUTINE debugout
END MODULE chem
