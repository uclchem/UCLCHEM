CALL GET_COMMAND_ARGUMENT(1, paramFile)  
IF (paramFile .ne. "") THEN
    open(fileNum, file=paramFile)
    DO WHILE (ios == 0)
        READ(fileNum, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1

            ! Find the first instance of whitespace.  Split label and data.
            pos = scan(buffer, ' ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)
            SELECT CASE (label)
            CASE('initialTemp')
                READ(buffer,*,iostat=ios) initialTemp
            CASE('maxTemp')
                READ(buffer,*,iostat=ios) maxTemp
            CASE('initialDens')
                READ(buffer,*,iostat=ios) initialDens
            CASE('finalDens')
                READ(buffer,*,iostat=ios) finalDens
            CASE('currentTime')
                READ(buffer,*,iostat=ios) currentTime
            CASE('finalTime')
                READ(buffer,*,iostat=ios) finalTime
            CASE('radfield')
                READ(buffer,*,iostat=ios) radfield
            CASE('zeta')
                READ(buffer,*,iostat=ios) zeta
            CASE('fr')
                READ(buffer,*,iostat=ios) fr
            CASE('rout')
                READ(buffer,*,iostat=ios) rout
            CASE('rin')
                READ(buffer,*,iostat=ios) rin
            CASE('baseAv')
                READ(buffer,*,iostat=ios) baseAv
            CASE('points')
                READ(buffer,*,iostat=ios) points
            CASE('switch')
                READ(buffer,*,iostat=ios) switch
            CASE('collapse')
                READ(buffer,*,iostat=ios) collapse
            CASE('bc')
                READ(buffer,*,iostat=ios) bc
            CASE('readAbunds')
                READ(buffer,*,iostat=ios) readAbunds
            CASE('phase')
                READ(buffer,*,iostat=ios) phase
            CASE('desorb')
                READ(buffer,*,iostat=ios) desorb
            CASE('h2desorb')
                READ(buffer,*,iostat=ios) h2desorb
            CASE('crdesorb')
                READ(buffer,*,iostat=ios) crdesorb
            CASE('uvcr')
                READ(buffer,*,iostat=ios) uvcr
            CASE('instantSublimation')
                READ(buffer,*,iostat=ios) instantSublimation
            CASE('ion')
                READ(buffer,*,iostat=ios) ion
            CASE('tempindx')
                READ(buffer,*,iostat=ios) tempindx
            CASE('fhe')
                READ(buffer,*,iostat=ios) fhe
            CASE('fc')
                READ(buffer,*,iostat=ios) fc
            CASE('fo')
                READ(buffer,*,iostat=ios) fo
            CASE('fn')
                READ(buffer,*,iostat=ios) fn
            CASE('fs')
                READ(buffer,*,iostat=ios) fs
            CASE('fmg')
                READ(buffer,*,iostat=ios) fmg
            CASE('fsi')
                READ(buffer,*,iostat=ios) fsi
            CASE('fcl')
                READ(buffer,*,iostat=ios) fcl
            CASE('fp')
                READ(buffer,*,iostat=ios) fp
            CASE('ff')
                READ(buffer,*,iostat=ios) ff
            CASE('outSpecies')
                pos = scan(buffer, ' ')
                READ(buffer(1:pos),*,iostat=ios) readIndx
                ALLOCATE(outSpecies(readIndx))
                READ(buffer(pos:),*,iostat=ios) outSpecies
            CASE('writeStep')
                READ(buffer,*,iostat=ios) writeStep 
            CASE('abundFile')
                READ(buffer,*,iostat=ios) abundFile
            CASE('outputFile')
                READ(buffer,*,iostat=ios) outputFile
            CASE('columnFile')
                READ(buffer,*,iostat=ios) columnFile 
                columnFileRead=.True.
            CASE('ebmaxh2')
                READ(buffer,*,iostat=ios) ebmaxh2
            CASE('epsilon')
                READ(buffer,*,iostat=ios) epsilon
            CASE('ebmaxcrf')
                READ(buffer,*,iostat=ios) ebmaxcrf
            CASE('uvcreff')
                READ(buffer,*,iostat=ios) uvcreff
            CASE('ebmaxcr')
                READ(buffer,*,iostat=ios) ebmaxcr
            CASE('phi')
                READ(buffer,*,iostat=ios) phi
            CASE('ebmaxuvcr')
                READ(buffer,*,iostat=ios) ebmaxuvcr
            CASE('uv_yield')
                READ(buffer,*,iostat=ios) uv_yield
            CASE('omega')
                READ(buffer,*,iostat=ios) omega
            CASE('vs')
                READ(buffer,*,iostat=ios) vs
            CASE('alpha')
                pos = scan(buffer, ' ')
                READ(buffer(1:pos),*,iostat=ios) readIndx
                READ(buffer(pos+1:),*,iostat=ios) alpha(readIndx)
            CASE DEFAULT
                WRITE(*,*) 'Skipping invalid label at line', line
            END SELECT
            IF (ios.ne.0) WRITE(*,*) "Problem with line ", label(1:pos),buffer
        END IF
    END DO
END IF
open(10,file=outputFile,status='unknown')

!Optionally, user can set species and have their abundances written out in columns with temp,dens and time.
!This checks whether those species were chosen and creates the output. Also warns users if set up isn't complete.
columnFlag=ALLOCATED(outSpecies) 
IF (columnFlag) THEN
    open(11,file=columnFile,status='unknown')
ELSE IF (columnFileRead) THEN
    write(*,*) "columnFile read but no species chosen, see example inputs for how to set outSpecies array."
    write(*,*) "Continuing without columnated output."
END IF

open(7,file=abundFile,status='unknown')