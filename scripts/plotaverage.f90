!This program plots the density weighted average of fractional abundance over depth points.
!It prints columns of Tage, average density, fractional abundances
!Tage is also averaged but if your age varies much between depth points, there are problems with your output.

!Much of the functionality of this program is already embedded in UCL_CHEM
!You can set outindx in parameters.f90 and the timestep interval for which to print outputs
!Then it plots log age,dens, fractional abundances.

PROGRAM plotaverage
IMPLICIT NONE
!This program will pick species of interest and output frac abundance
!you only need to change the three lines below the variable declarations
! ie. open files 1 and 2 and set outindx.

!variable declarations.
character*10 :: speci
character*92 :: test
integer i,j,k,nspec,depth,points, outindx(6)
double precision tdum,dens,rad,zeta,gr,fc,fo,fmg,fhe
double precision outy(6),tage,denstore,tstore,av
double precision, allocatable :: abund(:)


!Here are the user options
!-----------------------------------------------------------------------------------
! Open the output file from UCL_CHEM
open (1,file='output',status='old')

! Open the output file(s) to be read in DIPSO
open (2,file='plotout',status='unknown')

!identify species here for output
outindx=(/66,77,30,80,147,3/)
!-----------------------------------------------------------------------------------



!First read in header of first time/depth step
read(1,8000) tstore,denstore,tdum,av,rad,zeta,gr,fc,fo,fmg,fhe,depth
!skip empty lines
read(1,8010)

!Read file until no longer reading in species. This is number of lines of species.
!four species per line so nspec= lines * 4
nspec=0
read(1,'(a92)') test
DO WHILE (test(2:2) .NE. ' ')
    nspec=nspec+1
    read(1,'(a92)') test
END DO
nspec=nspec*4

!allocate enough space in arrays for number of species
allocate(abund(nspec))


read(1,'(a92)') test

!No read file until depth=1 again, this gives number of depth points
depth=2
points=0
DO WHILE (depth .ne. 1)
    points=points+1
    read(1,8000) tstore,denstore,tdum,av,rad,zeta,gr,fc,fo,fmg,fhe,depth
    read(1,'(a92)')
    read(1,'(a92)')
    DO i=1,(nspec/4)
        read(1,'(a92)') test
    END DO
    read(1,'(a92)')
    read(1,'(a92)')
END DO

!rewind file to start reading in.
rewind(1)

!quick output for users to check
write(*,*) " "
write(*,*) "Number of species:",nspec," ...assuming all rows are full"
write(*,*) "Number of depth points:",points
write(*,*) "Be aware that the program will exit ungracefully at end of file"
write(*,*) "This is not a problem, your output will be complete."
write(*,*) "Compare the final value of tage to your outputs if necessary"



!loop over time steps and depth
DO i=1,1000
    tage=0.0
    dens=0.0
    outy=0.0
    DO k=1,points
        !read the cloud properties
        read(1,8000) tstore,denstore,tdum,av,rad,zeta,gr,fc,fo,fmg,fhe
        !skip lines
        read(1,8010) 
        !read(1,'(a92)') test
        !read species in
        read(1,8020) (speci,abund(j),j=1,nspec)
        !skip lines
        read(1,8010)
 
        !Format for cloud properties block
        8000 format(&
        &33x,1pd11.3,5x,/&
        &33x,0pf15.4,5x,/,&
        &33x,0pf8.2,2x,/,&
        &33x,0pf12.4,4x,/&
        &33x,0pf10.2,13x,/,&
        &33x,0pf10.2,20x,/,&
        &33x,1pe8.2,8x,/,&
        &11x,1pe7.1,4x,12x,1pe7.1,/&
        &12x,1pe7.1,13x,1pe7.1,&
        &13x,i3)
        !Line skip
        8010  format(/)
        !Format for species
        8020  format(4(1x,a10,2x,D10.3,:))
        tage=tage+tstore
        dens=dens+denstore
        outy=outy+(abund(outindx)*denstore)
    END DO
    tage=log10(tage/points)
    outy=log10(outy/dens)
    dens=log10(dens/points)
    write(2,33)tage,dens,outy
    33 format(f9.4,1x,f9.4,1x,f9.4,1x,5(1x,f8.4))
END DO
END PROGRAM plotaverage