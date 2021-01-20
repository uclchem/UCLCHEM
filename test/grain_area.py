radii=[3.57699069369688527E-007,
	6.98799352454772693E-007,
	1.66467328055428966E-006, 
	3.55645670102603961E-006,  
	6.38977929313982604E-006, 
	9.97829390939405377E-006,  
	1.40343774171740631E-005,  
	1.82200037492632748E-005,  
	2.21826037436107715E-005,  
	2.55871898379409544E-005,  
	2.81455972293716657E-005,  
	2.96417164382822754E-005]

weights=[5792151389936155.0,
		2051685967306352.5,
		235554962470170.16,
		28071872396201.828,
		4932081459260.9521,
		1234971510565.4812,
		403598849340.45422,
		160061512934.20764,
		72679832471.918198,
		35689784319.592133,
		17382741166.856663,
		6457454266.3027544]

average_area=0.0
total_weight=0.0
for i,radius in enumerate(radii):
	area=radius*radius*4*3.141592654
	average_area+=weights[i]*area
	total_weight+=weights[i]

smallest_area=4*3.141592654*radii[0]**2.0
biggest_area=4*3.141592654*radii[-1]**2.0
average_area=average_area/total_weight
print(smallest_area,average_area,biggest_area)