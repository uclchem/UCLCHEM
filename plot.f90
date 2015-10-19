
PROGRAM plot
                 
!This program will pick species of interest and output frac abundance    

character*9 speci(300)
integer i,j,ncount,nspec,depth,dd,k
double precision tdum,dens,rad,zeta,gr,fc,fo,fmg,fhe
double precision xa1,xa2,xa3,xa4,xa5,xa6,tage,b(300,10)

! Here open the output of the PDR/chemical model4
open (1,file='sulfur-phase1-1',status='old')

!Here open the output file(s) to be read in DIPSO
open (2,file='plotsulfur',status='unknown')

! here define ncount=number of time spteps
! nspec=total number of species, dd=number of depth points

ncount=1091;nspec=219;dd=1
!loop over time steps and depth
DO i=1,ncount
  DO k=1,dd
    write(*,*)i
    !read the cloud properties
    read(1,8000) tage,dens,tdum,av,rad,zeta,gr,fc,fo,fmg,fhe
    !skip lines
    read(1,8010)
    !read species in
    read(1,8020) (speci(j),b(j,k),j=1,nspec)
    !skip lines
    read(1,8010)

    !fixed format for cloud properties block
    8000 format(&
    &'age of cloud             time  = ',1pd11.3,' years',/,&
    &'total hydrogen density   dens  = ',0pf15.4,' cm-3',/,&
    &'cloud temperature        temp  = ',0pf8.2,' k',/,&
    &'visual extinction        av    = ',0pf12.4,' mags',/,&
    &'radiation field          rad   = ',0pf10.2,' (habing = 1)',/,&
    &'cosmic ray ioniz. rate   zeta  = ',0pf10.2,' (unit = 1.3e-17'&
    &'s-1)',/,&
    &'h2 formation rate coef.        = ',1pe8.2,' cm3 s-1',/,&
    &'c / htot = ',1pe7.1,4x,' o / htot = ',1pe7.1,/&
    &'mg / htot = ',1pe7.1,' he / htot = ',1pe7.1,' depth     = ',i3)
    8010 format(//)
    !format for species output
!    8020  format(4(1x,a8,'=',1x,1pd10.3,:))
    8020  format(4(1x,a8,'=',1x,d8.3E3,:))

    !output as log10 for convenience
    dens=log10(dens)
    xa1=log10(b(2,k))
    xa2=log10(b(219,k))
    xa3=log10(b(220,k))
    xa4=log10(b(218,k))
    xa5=log10(b(4,k))
    xa6=log10(b(3,k))

    !write this time step then move on
    write(2,33)tage,dens,xa1,xa2,xa3,xa4,xa5
    33 format(f9.4,1x,f9.4,1x,f9.4,1x,5(1x,f8.4))
  END DO
END DO
END PROGRAM PLOT