! Chemistry module of UCL_CHEM. Contains all the core machinery of the code, not really intended to be altered.
! Use physics module to alter temp/density behaviour etc. This module should solve chemistry for a cloud of gas
MODULE chem
USE physics
IMPLICIT NONE
EXTERNAL dlsode
    !makerates gives these numbers, nspec includes electrons, ngrain is number of mantle species
    integer,parameter :: nreac=2462,nspec=209,ngrain=52,nout=6

    !These integers store the species index of important species, x is for ions    
    integer :: nh,nh2,nc,ncx,no,nn,ns,nhe,nco,nmg,nh2o,nsi,nsix,ncl,nclx,nch3oh

    !These integer arrays store numbers labelling reactions, species and specific reactions of interest
    integer ::reacindx(nreac),specindx(nspec), nrco,mantleindx(ngrain),outindx(nout),writestep

    !loop counters    
    integer :: i,j,l

    !These are variables for reaction rates, alpha/beta/gamas are combined each time step to make rate,the total reaction rate
    double precision :: rate(nreac),alpha(nreac),beta(nreac),gama(nreac),mass(nspec)
    character(LEN=10) :: re1(nreac),re2(nreac),re3(nreac),p1(nreac),p2(nreac),p3(nreac),p4(nreac)
    character(LEN=10) :: specname(nspec)  
    
    !DLSODE variables    
    integer :: ITOL,ITASK,ISTATE,IOPT,MESFLG,LUNIT,NEQ,mf
    integer :: LRW,LIW,MXSTEP,IWORK(500)
    double precision :: RTOL,ATOL,RWORK(100000)

    !initial fractional elemental abudances and arrays to store abundances
    double precision :: fh,fhe,fc,fo,fn,fs,fmg,fsi,fcl,h2col,cocol,junk1,junk2
    double precision :: y(nspec+1),abund(nspec+1,points),ax(nspec+1),mantle
    
    !Variables controlling chemistry
    double precision :: radfield,zeta,fr,omega,grain,radg,cion,h2form,h2dis
    double precision :: ebmaxh2,epsilon,ebmaxcrf,ebmaxcr,phi,ebmaxuvcr,uvy
    double precision :: taud,dopw,radw,xl,fosc

    logical :: startr
    integer :: dimco, dimh2
    double precision :: corates(7,6), y2r(7,6), ncogr(7), nh2gr(6)    

CONTAINS
!This gets called immediately by main so put anything here that you want to happen before the time loop begins, reader is necessary.
    SUBROUTINE initialise
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
       write(*,*) abund(nspec+1,:),dens,d0
        !h2 formation rate initially set
        h2form = 1.0d-17*dsqrt(temp)
        
    END SUBROUTINE initialise

!Reads input reaction and species files as well as the final step of previous run if this is phase 2
    SUBROUTINE reader
        IMPLICIT NONE
        integer i,j,l,m

        !read start file IF not first phase to get all the abundances from the last step of previous phase as well as density and temp
        IF (first .eq. 0) THEN
            DO l=1,points
                read(7,7040)dens
                read(7,7050)temp
                read(7,7030)
                read(7,7020)tage,dstep
                read(7,7010) (specname(i),abund(i,l),i=1,nspec)
                read(7,7000)
            END DO
            7000  format(//)
            7010  format(4(1x,a8,'=',1x,1pd10.3,:))
            7020  format(1x,'time is (yrs) ',1pd11.3,' depth = ',i3,//)
            7030  format(////////)
            7040  format("total hydrogen density   dens  =    ",D9.3," cm-3")
            7050  format("cloud temperature        temp  =    ",F4.2," k")
        END IF

        !read species file and assign specindx for important species to the integers used to store them.
        !format was 8010, mantle is dummy for the abundances
        read(3,*)(specindx(j),specname(j),junk1,mass(j),j=1,nspec-1)
        write(71,8010)(specindx(j),specname(j),mass(j),j=1,nspec-1)
        8010  format(i4,1x,a8,1x,f5.1)
        specname(nspec)='electr'
        j=1
        DO i=1,nspec-1
            IF (specname(i).eq.'H')  nh  = specindx(i)
            IF (specname(i).eq.'H2')  nh2 = specindx(i)
            IF (specname(i).eq.'C')  nc  = specindx(i)
            IF (specname(i).eq.'C+')  ncx = specindx(i)
            IF (specname(i).eq.'O')  no  = specindx(i)
            IF (specname(i).eq.'N')        nn  = specindx(i)
            IF (specname(i).eq.'S+')  ns  = specindx(i)
            IF (specname(i).eq.'HE')  nhe = specindx(i)
            IF (specname(i).eq.'CO')  nco = specindx(i)
            IF (specname(i).eq.'MG')  nmg = specindx(i)
            IF (specname(i).eq.'H2O')  nh2o = specindx(i)
            IF (specname(i).eq.'SI')  nsi = specindx(i)
            IF (specname(i).eq.'SI+')  nsix= specindx(i)
            IF (specname(i).eq.'CL')  ncl = specindx(i)
            IF (specname(i).eq.'CL+')  nclx= specindx(i)
            IF (specname(i).eq.'CH3OH')  nch3oh= specindx(i)

            !Finds all the mantle species and adds their index to a list called mantleindx
            IF (specname(i)(:1).eq.'#')  THEN
                mantleindx(j)=specindx(i)
                write(78,*) mantleindx(j)
                j=j+1
            END IF
        END DO

        !read reac file, alpha, beta and gama are used for working out reaction rate each time step
        DO j=1,nreac
            read(2,*) reacindx(j),re1(j),re2(j),re3(j),p1(j),p2(j),p3(j),&
            &p4(j),alpha(j),beta(j),gama(j),junk1,junk2
            write(72,8020) reacindx(j),re1(j),re2(j),re3(j),p1(j),p2(j),p3(j),&
            &p4(j),alpha(j),beta(j),gama(j)
        END DO    
        8020  format(i4,5(1x,a8),2(1x,a4),1x,1pe8.2,3x,1pe8.2,2x,1pe8.2)
            END SUBROUTINE reader

!Writes physical variables and fractional abundances to output file, called every time step.
    SUBROUTINE output
        DO l=1,nspec
            IF(abund(l,dstep) .le. 1.d-30) THEN
                abund(l,dstep)= 1.d-30
            ENDIF
        END DO
        write(1,8020) tage,dens,temp,av(dstep),radfield,zeta,h2form,fc,fo,&
                        &fmg,fhe,dstep
        write(1,8000)
        write(1,8010) (specname(i),abund(i,dstep),i=1,nspec) 
        write(1,8000)
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
        &'cosmic ray ioniz. rate   zeta  = ',0pf10.2,' (unit = 1.3e-17'&
        &'s-1)',/,&
        &'h2 formation rate coef.        = ',1pe8.2,' cm3 s-1',/,&
        &'c / htot = ',1pe7.1,4x,' o / htot = ',1pe7.1,/&
        &'mg / htot = ',1pe7.1,&
        &'he / htot = ',1pe7.1,&
        &'depth     = ',i3)

        IF ( mod(tstep,writestep) .eq. 0) THEN
          write(4,8030) tage,dens,Y(outindx)
          8030  format(1pd11.3,1x,0pf15.4,6(1x,1pd10.3))
        END IF
    END SUBROUTINE output

    SUBROUTINE evaporate
    END SUBROUTINE evaporate

!Called every time/depth step and updates the abundances of all the species
    SUBROUTINE chem_update
        !y is at final value of previous depth iteration so set to initial values of this depth with abund
        !reset other variables for good measure        
        y=abund(:,dstep)
        dens=abund(nspec+1,dstep)
        ax=y*dens
        h2form = 1.0d-17*dsqrt(temp)

    
        !evaluate co and h2 column densities for use in rate calculations
        !sum column densities of each point up to dstep. boxsep and dens are pulled out of the sum as common factors  
        IF (dstep.gt.1 .or. points .eq. 1) THEN
            h2col=sum(abund(nh2,:dstep))*dens*size*(real(dstep)/real(points))
            cocol=sum(abund(nco,:dstep))*dens*size*(real(dstep)/real(points))
        ELSE
            h2col=0.
            cocol=0.
        ENDIF
        !call the actual ODE integrator
        call integrate
        !Set abundances to output of DLSODE
        abund(:,dstep)=y
    END SUBROUTINE chem_update

!This subroutine calls DLSODE (3rd party ODE solver) until it can reach tout with acceptable errors (RTOL/ATOL)
    SUBROUTINE integrate
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
            write(79,*) Y(nspec)

        END DO
        !DLSODE USES THIS COUNTER TO CHECK IF IT HAS INTEGRATED ODES BEFORE, NEEDS TO FORGET THAT BETWEEN TIME STEPS
        !ISTATE=1                    
    END SUBROUTINE integrate

    include 'rates.f90'

    !DLSODE calls this subroutine to ask it what the RHS of the equations dy/dt=... are    
    SUBROUTINE  F (NEQ, T, Y, YDOT)
        INTEGER :: NEQ
        DOUBLE PRECISION :: T,Y(nspec+1),YDOT(nspec+1)
        DOUBLE PRECISION :: D,loss,prod
        
        !Dens is updated by DLSODE just like abundances so this ensures dens is at correct value for this timestep
        dens=y(nspec+1)
        !Set D to the gas density for use in the ODEs
        D=dens
        !The ODEs created by MakeRates go here, they are essentially sums of terms that look like k(1,2)*y(1)*y(2)*dens. Each species ODE is made up
        !of the reactions between it and every other species it reacts with.
        INCLUDE 'odes.f90'
        !Sum of abundaces of all mantle species. mantleindx stores the indices of mantle species.
        !JH: I should test the value of mantle against Serenas code to check F95 intrinsic functions do what I think they do
        mantle=sum(y(mantleindx))
        h2form=1.0d-17*dsqrt(temp)        
        
        !H2 formation should occur at both steps - however note that here there is no 
        !temperature dependence. y(nh) is hydrogen fractional abundance.      
        ydot(nh)  = ydot(nh) - 2.0*( h2form*dens*y(nh) - h2dis*y(nh2) )
        !                             h2 formation - h2-photodissociation
        ydot(nh2) = ydot(nh2) + h2form*dens*y(nh) - h2dis*y(nh2)
        !                       h2 formation  - h2-photodissociation

        !IF (evap .ne. 2)  THEN
        !   radc=newradc
        !   rad=newrad
        !ENDIF

        ! get density change from physics module to send to DLSODE
        IF (collapse .eq. 1) THEN
            ydot(nspec+1)=densdot()
        ENDIF
        write(79,*) Y(nspec)
        ydot(nspec)=0.0


    END SUBROUTINE F

!integrate calls reacrates to get the reaction rates at every iteration. reacrates calls further functions.
!This file is already long so I've hidden those subroutines in rates.f90

!This is a dummy for DLSODE, it has to call it but we do not use it.
    subroutine JAC (NEQ,T,Y,ML,MU,PD,NROWPD)
         INTEGER  NEQ, ML, MU, NROWPD
         DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
    END SUBROUTINE JAC        
END MODULE chem
