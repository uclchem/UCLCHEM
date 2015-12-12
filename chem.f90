! Chemistry module of UCL_CHEM. Contains all the core machinery of the code, not really intended to be altered.
! Use physics module to alter temp/density behaviour etc. This module should solve chemistry for a cloud of gas
MODULE chem
USE physics
IMPLICIT NONE
EXTERNAL dlsode
    !makerates gives these numbers, nspec includes electrons
    integer,parameter :: nout=6

    !These integers store the array index of important species and reactions, x is for ions    
    integer :: nh,nh2,nc,ncx,no,nn,ns,nhe,nco,nmg,nh2o,nsi,nsix,ncl,nclx,nch3oh
    integer ::nrco,outindx(nout),nspec,nreac,njunk,evapevents
    
    !loop counters    
    integer :: i,j,l,writestep

    !These are variables for reaction rates, alpha/beta/gamas are combined each time step to make rate,the total reaction rate
    double precision,allocatable :: rate(:),alpha(:),beta(:),gama(:),mass(:)
    character(LEN=10),allocatable :: re1(:),re2(:),re3(:),p1(:),p2(:),p3(:),p4(:)
    character(LEN=10),allocatable :: specname(:)  
    
    !DLSODE variables    
    integer :: ITOL,ITASK,ISTATE,IOPT,MESFLG,LUNIT,NEQ,mf
    integer :: LRW,LIW,MXSTEP,IWORK(500)
    double precision :: RTOL,ATOL,RWORK(100000)

    !initial fractional elemental abudances and arrays to store abundances
    double precision :: fh,fhe,fc,fo,fn,fs,fmg,fsi,fcl,h2col,cocol,mantle,junk1,junk2
    double precision,allocatable :: y(:),abund(:,:)
    
    !Variables controlling chemistry
    double precision :: radfield,zeta,fr,omega,grain,radg,cion,h2form,h2dis
    double precision :: ebmaxh2,epsilon,ebmaxcrf,ebmaxcr,phi,ebmaxuvcr,uvy,uvcreff
    double precision :: taud,dopw,radw,xl,fosc

    !evaporation lists, these are used to evaporation species in specific events
    !See viti 2004 for more information.
    integer, allocatable :: colist(:),mcolist(:),intlist(:),mintlist(:),grainlist(:),mgrainlist(:)
    integer, allocatable :: co2list(:),mco2list(:),int2list(:),mint2list(:)
    double precision, allocatable :: cobindener(:),co2bindener(:),intbindener(:)

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
            abund(nspec+1,:)=dens      
            !abund(nfe,:) = ffe
            !abund(nna,:) = fna

            !Decide how much iron is initiall ionized using parameters.f90
            SELECT CASE (ion)
                CASE(0)
                    abund(nc,:)=fc
                    abund(ncx,:)=1.d-10
                CASE(1)
                    abund(nc,:)=fc/2.0
                    abund(ncx,:)=fc/2.0
                CASE(2)
                    abund(nc,:)=1.d-10
                    abund(ncx,:)=fc
            END SELECT
            abund(nspec,:)=abund(ncx,:)

       ENDIF
       !h2 formation rate initially set
       h2form = 1.0d-17*dsqrt(temp)
        
    END SUBROUTINE chem_initialise

!Reads input reaction and species files as well as the final step of previous run if this is phase 2
    SUBROUTINE reader
        IMPLICIT NONE
        integer i,j,l,m

        !read species file and allocate sufficient space to relevant arrays
        read(3,*)nspec
        allocate(y(nspec+1),abund(nspec+1,points),specname(nspec),mass(nspec))
        read(3,*)(specname(j),junk1,mass(j),j=1,nspec-1)
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
        END DO

        !read reac file, assign array space
        !alpha, beta and gama are used for working out reaction rate each time step
        read(2,*) nreac
        allocate(re1(nreac),re2(nreac),re3(nreac),p1(nreac),p2(nreac),p3(nreac),&
            &p4(nreac),alpha(nreac),beta(nreac),gama(nreac),rate(nreac))
        DO j=1,nreac
            read(2,*) re1(j),re2(j),re3(j),p1(j),p2(j),p3(j),&
            &p4(j),alpha(j),beta(j),gama(j),junk1,junk2
        END DO    

        !finally, read in evaporation lists
        !what we do is read in length of each list, then the numbers in the list
        !these are used for evaporation sums.
        read(8,*) l
        allocate(colist(l))
        allocate(mcolist(l))
        read(8,*)colist
        read(8,*)mcolist
        read(8,*) l
        allocate(intlist(l))
        allocate(mintlist(l))
        read(8,*) intlist
        read(8,*) mintlist
        read(8,*)l
        allocate(grainlist(l))
        allocate(mgrainlist(l))
        read(8,*) grainlist
        read(8,*) mgrainlist

        !read start file IF not first phase to get finale abundances from previous phase 
        !density, temp and av read but NOT zeta or radfield
        IF (first .eq. 0) THEN
            DO l=1,points
                read(7,*)
                read(7,7000) dens,temp,av(l)
                write(*,*)dens,temp,av(l)
                read(7,*)
                read(7,7010) h2form,fc,fo,&
                            &fmg,fhe,dstep
                read(7,*)
                read(7,7030) (specname(i),abund(i,l),i=1,nspec)
            END DO
            7000 format(&
            &33x,0pf15.4,5x,/,&
            &33x,0pf8.2,2x,/,&
            &33x,0pf12.4,4x,/)
            7010 format(&
            &33x,1pe8.2,8x,/,&
            &11x,1pe7.1,4x,12x,1pe7.1,/&
            &12x,1pe7.1,13x,1pe7.1,&
            &13x,i3,/)
            7020  format(//)
            7030  format(4(1x,a8,2x,1pd10.3,:))     
        END IF
    END SUBROUTINE reader

!Writes physical variables and fractional abundances to output file, called every time step.
    SUBROUTINE output
        !1.d-30 stops numbers getting too small for fortran.
        DO l=1,nspec
            IF(abund(l,dstep) .le. 1.d-30) THEN
               abund(l,dstep)= 1.d-30
            ENDIF
        END DO
        !write out cloud properties
        write(1,8020) tage,dens,temp,av(dstep),radfield,zeta,h2form,fc,fo,&
                        &fmg,fhe,dstep
        !and a blank line
        write(1,8000)
        !and then all the abundances for this step
        write(1,8010) (specname(i),abund(i,dstep),i=1,nspec) 
        write(1,8000)
        !If this is the last time step of phase I, write a start file for phase II
        IF (first .eq. 1) THEN
           IF (switch .eq. 0 .and. tage .ge. tfin& 
               &.or. switch .eq. 1 .and.dens .ge. dfin) THEN
               write(7,8020) tage,dens,temp,av(dstep),radfield,zeta,h2form,fc,fo,&
                       &fmg,fhe,dstep
               write(7,8000)
               write(7,8010) (specname(i),abund(i,dstep),i=1,nspec)
               write(7,8000)
           ENDIF
        ENDIF
        8000  format(/)
        8010  format(4(1x,a8,'=',1x,1pd10.3,:))
        8020 format(&
        &'age of cloud             time  = ',1pd11.3,' years',/,&
        &'total hydrogen density   dens  = ',0pf15.4,' cm-3',/,&
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
          write(4,8030) tage,dens,Y(outindx)
          8030  format(1pd11.3,1x,0pf15.4,6(1x,1pd10.3))
        END IF
    END SUBROUTINE output

    SUBROUTINE chem_update
    !Called every time/depth step and updates the abundances of all the species

        !y is at final value of previous depth iteration so set to initial values of this depth with abund
        !reset other variables for good measure        
        y=abund(:,dstep)
        h2form = 1.0d-17*dsqrt(temp)
    
        !evaluate co and h2 column densities for use in rate calculations
        !sum column densities of each point up to dstep. boxlength and dens are pulled out of the sum as common factors  
        IF (dstep.gt.1) THEN
            h2col=(sum(abund(nh2,:dstep-1))+0.5*abund(nh2,dstep))*dens*(size/real(points))
            cocol=(sum(abund(nco,:dstep-1))+0.5*abund(nco,dstep))*dens*(size/real(points))
        ELSE
            h2col=0.5*abund(nh2,dstep)*dens*(size/real(points))
            cocol=0.5*abund(nco,dstep)*dens*(size/real(points))
        ENDIF

        !call the actual ODE integrator
        CALL integrate

        !call evaporation to remove species from grains at certain temperatures
        CALL evaporate

        !Set abundances to output of DLSODE
        abund(:,dstep)=y
    END SUBROUTINE chem_update

    SUBROUTINE integrate
    !This subroutine calls DLSODE (3rd party ODE solver) until it can reach tout with acceptable errors (RTOL/ATOL)

        DO WHILE(t0 .lt. tout)            
            !reset parameters for DLSODE
            ITOL=1
            RTOL=1e-15
            ATOL=1e-20
            MF=22
            ITASK=1
            NEQ=nspec+1
            IF(ISTATE .EQ. 1) THEN
                IOPT=1
                IWORK(6)=MXSTEP
            ENDIF
            IF(MXSTEP .GT. IWORK(6)) THEN
                ISTATE=3
                IOPT=1
                IWORK(6)=MXSTEP
            ENDIF
            write(*,*) IWORK(6),' ',MXSTEP,' ',IOPT,' ',ISTATE
           
            !get reaction rates for this iteration
            CALL reacrates
            !Call the integrator.
            CALL DLSODE(F,NEQ,Y,T0,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
            &             RWORK,LRW,IWORK,LIW,JAC,MF)

            IF(ISTATE.EQ.2) THEN
                write(*,*)'Call to LSODE successful at time: ',(TOUT*year),' years'
                write(*,*)'        Steps: ',IWORK(6)
                IOPT=0
            ELSEIF(ISTATE.EQ.-1) THEN
                write(*,*)'Call to LSODE returned -1 meaning that MXSTEP exceeded'
                write(*,*)'but the integration was successful'
                write(*,*)'Doubling MXSTEP (thanks Tom :) ) from:',MXSTEP,' to:',MXSTEP*2
                MXSTEP=MXSTEP*2
            ELSE
                write(*,*)'ISTATE value: ', ISTATE,'...no trap written just yet!'
                stop
            ENDIF
            ISTATE=2

        END DO
        !DLSODE USES THIS COUNTER TO CHECK IF IT HAS INTEGRATED ODES BEFORE, NEEDS TO FORGET THAT BETWEEN TIME STEPS
        ISTATE=1                    
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
        IF (collapse .eq. 0) THEN
            y(nspec+1)=dens
        ELSE
            dens=y(nspec+1)
        END IF

        !Set D to the gas density for use in the ODEs
        D=dens

        !The ODEs created by MakeRates go here, they are essentially sums of terms that look like k(1,2)*y(1)*y(2)*dens. Each species ODE is made up
        !of the reactions between it and every other species it reacts with.
        INCLUDE 'odes.f90'

        !Sum of abundaces of all mantle species. mantleindx stores the indices of mantle species.
        mantle=sum(y(mgrainlist))
        !updated just in case temp changed
        h2form=1.0d-17*dsqrt(temp) 

        !H2 formation should occur at both steps - however note that here there is no 
        !temperature dependence. y(nh) is hydrogen fractional abundance.      
        ydot(nh)  = ydot(nh) - 2.0*( h2form*dens*y(nh) - h2dis*y(nh2) )
        !                             h2 formation - h2-photodissociation
        ydot(nh2) = ydot(nh2) + h2form*dens*y(nh) - h2dis*y(nh2)
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
    IF (tstep .gt. 0) THEN
        !Viti 04 evap
        IF (evap .eq. 1) THEN
            ! this code only works if you use the 288,000 temp increase from old uclchem
            !if you do y(a)=y(a)+y(b) with a and b as integer arrays. uses corresponding elements
            !eg. y(a[1])=y(a[1])+y(b[1]). So make sure they match by comparing evaplist to species.csv

            !Solid Evap
            IF (solidflag .eq. 1) THEN
                y(colist)=y(colist)+0.35*y(mcolist)
                y(mcolist)=0.65*y(mcolist)
                !Set flag to 2 to stop it being recalled
                solidflag=2
            ENDIF

            !mono evap
            CALL temper

            !Volcanic, int should really be separated out.
            IF (volcflag .eq. 1) THEN
                y(colist)=y(colist)+0.667*y(mcolist)
                y(mcolist)=0.333*y(mcolist)
                y(co2list)=y(co2list)+0.667*y(mco2list)
                y(mco2list)=0.333*y(mco2list)
                y(intlist)=y(intlist)+0.5*y(mintlist)
                y(mintlist)=0.5*y(mintlist)
                y(int2list)=y(int2list)+0.5*y(mint2list)
                y(mint2list)=0.5*y(mint2list)
                !Set flag to 2 to stop it being recalled
                volcflag=2
            ENDIF

            !Co-desorption
            IF (coflag .eq. 1) THEN
                y(grainlist)=y(grainlist)+y(mgrainlist)
                y(mgrainlist)=1d-30
                coflag=2
            ENDIF
        ELSE IF (evap .eq. 2 .and. coflag .ne. 2) THEN
            !Alternative evap. Instaneous evaporation of all grain species
            y(grainlist)=y(grainlist)+y(mgrainlist)
            y(mgrainlist)=1d-30
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
            en=cobindener(i)*kbolt
            expdust=cobindener(i)/temp
            newm = mass(speci)*1.66053e-27
            freq = dsqrt((2*(surden)*en)/((pi**2)*newm))
            kevap=freq*exp(-expdust)
            IF (kevap .ge. 0.99) THEN
                write(*,*)i
                y(speci)=y(speci)+(0.7*y(mcolist(i)))
                y(mcolist(i))=0.3*y(mcolist(i))
                !set to 1d50 so it can't happen again
                cobindener(i)=1d50
            END IF 
        END DO

        !co2 monoevap
        DO i=lbound(co2list,1),ubound(co2list,1)
            speci=co2list(i)
            en=co2bindener(i)*kbolt
            expdust=co2bindener(i)/temp
            newm = mass(speci)*1.66053e-27
            freq = dsqrt((2*(surden)*en)/((pi**2)*newm))
            kevap=freq*exp(-expdust)
            IF (kevap .ge. 0.99) THEN
                write(*,*)i
                y(speci)=y(speci)+(0.7*y(mco2list(i)))
                y(mco2list(i))=0.3*y(mco2list(i))
                co2bindener(i)=1d50
            END IF 
        END DO

        !int mono evap
        DO i=lbound(intlist,1),ubound(intlist,1)
            speci=intlist(i)
            en=intbindener(i)*kbolt
            expdust=intbindener(i)/temp
            newm = mass(speci)*1.66053e-27
            freq = dsqrt((2*(surden)*en)/((pi**2)*newm))
            kevap=freq*exp(-expdust)
            IF (kevap .ge. 0.99) THEN
                write(*,*)i
                y(speci)=y(speci)+(0.1*y(mintlist(i)))
                y(mintlist(i))=0.9*y(mintlist(i))
                intbindener(i)=1d50
            END IF 
        END DO
    END SUBROUTINE temper
!This is a dummy for DLSODE, it has to call it but we do not use it.
    subroutine JAC (NEQ,T,Y,ML,MU,PD,NROWPD)
         INTEGER  NEQ, ML, MU, NROWPD
         DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
    END SUBROUTINE JAC        
END MODULE chem
