nEquation=35, nReactions=18

grainIndices=(/3,6,23,29,31,36,38,42,77,86,94,96,104,106,108,119,127 &
     &       ,129,132,169,171,177,189,191,193,198,200,205,207,227,229 &
     &       ,237,239,241,243/)

gasIndices=(/2,5,22,28,30,35,37,41,76,85,93,95,103,105,107,118,126,128 &
     &       ,131,168,170,176,188,190,192,197,199,204,206,226,228,236 &
     &       ,238,240,242/)

H,H2,CH3,CH4,NH2,NH3,OH,H2O,CO,HCO,C2H6,H2CO,CH2OH,CH3NH2,CH3O,CH3OH,HS,NH2OH,H2S,CH3CHO,CO2,HCONH2,C2H5OH,CH3CH3O,HCOOH,CH2OHNH2,CH3ONH2,CH3OOH,CH2OHOH,CH2OHCHO,HCOOCH3,C2H4O2H2,CH3OCH3O,CH3OCH2OH,(CH2OH)2

      YDOT(1) = 0.0
      YDOT(2) = 0.0
      LOSS = -RATE(1)*Y(8)*Y(7)*D*D-2*RATE(6)*Y(8)*D*D-RATE(7)*Y(8) &
     &       *Y(5)*D*D-RATE(8)*Y(8)*Y(15)*D*D-RATE(9)*Y(13)*Y(8)*D*D &
     &       -RATE(10)*Y(8)*Y(10)*D*D
      YDOT(3) = Y(3)*LOSS
      YDOT(4) = 0.0
      LOSS = -RATE(2)*Y(8)*Y(7)*D*D-RATE(7)*Y(8)*Y(3)*D*D-RATE(11)*Y(8) &
     &       *Y(15)*D*D-RATE(12)*Y(13)*Y(8)*D*D-RATE(13)*Y(8)*Y(10)*D*D
      YDOT(5) = Y(5)*LOSS
      YDOT(6) = 0.0
      LOSS = -RATE(1)*Y(8)*Y(3)*D*D-RATE(2)*Y(8)*Y(5)*D*D-RATE(3)*Y(8) &
     &       *Y(15)*D*D-RATE(4)*Y(13)*Y(8)*D*D-RATE(5)*Y(8)*Y(10)*D*D
      YDOT(7) = Y(7)*LOSS
      LOSS = -RATE(1)*Y(3)*Y(7)*D*D-RATE(2)*Y(5)*Y(7)*D*D-RATE(3)*Y(15) &
     &       *Y(7)*D*D-RATE(4)*Y(13)*Y(7)*D*D-RATE(5)*Y(10)*Y(7)*D*D &
     &       -RATE(6)*Y(3)*Y(3)*D*D-RATE(7)*Y(3)*Y(5)*D*D-RATE(8)*Y(3) &
     &       *Y(15)*D*D-RATE(9)*Y(13)*Y(3)*D*D-RATE(10)*Y(3)*Y(10)*D*D &
     &       -RATE(11)*Y(15)*Y(5)*D*D-RATE(12)*Y(13)*Y(5)*D*D-RATE(13) &
     &       *Y(10)*Y(5)*D*D-RATE(14)*Y(15)*Y(15)*D*D-RATE(15)*Y(13) &
     &       *Y(15)*D*D-RATE(16)*Y(15)*Y(10)*D*D-RATE(17)*Y(13)*Y(13)*D &
     &       *D-RATE(18)*Y(13)*Y(10)*D*D
      PROD = +RATE(1)*Y(8)*Y(3)*Y(7)*D*D+RATE(2)*Y(8)*Y(5)*Y(7)*D*D &
     &       +RATE(3)*Y(8)*Y(15)*Y(7)*D*D+RATE(4)*Y(13)*Y(8)*Y(7)*D*D &
     &       +RATE(5)*Y(8)*Y(10)*Y(7)*D*D+RATE(6)*Y(8)*Y(3)*Y(3)*D*D &
     &       +RATE(7)*Y(8)*Y(3)*Y(5)*D*D+RATE(8)*Y(8)*Y(3)*Y(15)*D*D &
     &       +RATE(9)*Y(13)*Y(8)*Y(3)*D*D+RATE(10)*Y(8)*Y(3)*Y(10)*D*D &
     &       +RATE(11)*Y(8)*Y(15)*Y(5)*D*D+RATE(12)*Y(13)*Y(8)*Y(5)*D*D &
     &       +RATE(13)*Y(8)*Y(10)*Y(5)*D*D+RATE(14)*Y(8)*Y(15)*Y(15)*D &
     &       *D+RATE(15)*Y(13)*Y(8)*Y(15)*D*D+RATE(16)*Y(8)*Y(15)*Y(10) &
     &       *D*D+RATE(17)*Y(13)*Y(13)*Y(8)*D*D+RATE(18)*Y(13)*Y(8) &
     &       *Y(10)*D*D
      YDOT(8) = PROD+Y(8)*LOSS
      YDOT(9) = 0.0
      LOSS = -RATE(5)*Y(8)*Y(7)*D*D-RATE(10)*Y(8)*Y(3)*D*D-RATE(13) &
     &       *Y(8)*Y(5)*D*D-RATE(16)*Y(8)*Y(15)*D*D-RATE(18)*Y(13)*Y(8) &
     &       *D*D
      YDOT(10) = Y(10)*LOSS
      PROD = +RATE(6)*Y(8)*Y(3)*Y(3)*D*D
      YDOT(11) = PROD
      YDOT(12) = 0.0
      LOSS = -RATE(4)*Y(8)*Y(7)*D*D-RATE(9)*Y(8)*Y(3)*D*D-RATE(12)*Y(8) &
     &       *Y(5)*D*D-RATE(15)*Y(8)*Y(15)*D*D-2*RATE(17)*Y(8)*D*D &
     &       -RATE(18)*Y(8)*Y(10)*D*D
      YDOT(13) = Y(13)*LOSS
      PROD = +RATE(7)*Y(8)*Y(3)*Y(5)*D*D
      YDOT(14) = PROD
      LOSS = -RATE(3)*Y(8)*Y(7)*D*D-RATE(8)*Y(8)*Y(3)*D*D-RATE(11)*Y(8) &
     &       *Y(5)*D*D-2*RATE(14)*Y(8)*D*D-RATE(15)*Y(13)*Y(8)*D*D &
     &       -RATE(16)*Y(8)*Y(10)*D*D
      YDOT(15) = Y(15)*LOSS
      PROD = +RATE(1)*Y(8)*Y(3)*Y(7)*D*D
      YDOT(16) = PROD
      YDOT(17) = 0.0
      PROD = +RATE(2)*Y(8)*Y(5)*Y(7)*D*D
      YDOT(18) = PROD
      YDOT(19) = 0.0
      PROD = +RATE(10)*Y(8)*Y(3)*Y(10)*D*D
      YDOT(20) = PROD
      YDOT(21) = 0.0
      PROD = +RATE(13)*Y(8)*Y(10)*Y(5)*D*D
      YDOT(22) = PROD
      PROD = +RATE(9)*Y(13)*Y(8)*Y(3)*D*D
      YDOT(23) = PROD
      PROD = +RATE(8)*Y(8)*Y(3)*Y(15)*D*D
      YDOT(24) = PROD
      PROD = +RATE(5)*Y(8)*Y(10)*Y(7)*D*D
      YDOT(25) = PROD
      PROD = +RATE(12)*Y(13)*Y(8)*Y(5)*D*D
      YDOT(26) = PROD
      PROD = +RATE(11)*Y(8)*Y(15)*Y(5)*D*D
      YDOT(27) = PROD
      PROD = +RATE(3)*Y(8)*Y(15)*Y(7)*D*D
      YDOT(28) = PROD
      PROD = +RATE(4)*Y(13)*Y(8)*Y(7)*D*D
      YDOT(29) = PROD
      PROD = +RATE(18)*Y(13)*Y(8)*Y(10)*D*D
      YDOT(30) = PROD
      PROD = +RATE(16)*Y(8)*Y(15)*Y(10)*D*D
      YDOT(31) = PROD
      YDOT(32) = 0.0
      PROD = +RATE(14)*Y(8)*Y(15)*Y(15)*D*D
      YDOT(33) = PROD
      PROD = +RATE(15)*Y(13)*Y(8)*Y(15)*D*D
      YDOT(34) = PROD
      PROD = +RATE(17)*Y(13)*Y(13)*Y(8)*D*D
      YDOT(35) = PROD
