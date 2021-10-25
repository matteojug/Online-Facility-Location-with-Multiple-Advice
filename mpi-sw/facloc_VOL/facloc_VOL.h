// facloc_VOL.h
#ifndef _FACLOC_VOL
#define _FACLOC_VOL

/* Adjusted version of Volume algorithm

   UFLP-call interface:
   --------------------
   Please provide the following arrays:	

   open_cost   = double[facilities];
   cost_matrix = double[facilities*cities];

   All data non-negative. The connection cost from city j to facility i
   is denoted in position cost_matrix[(i*cities) + j]
  
   connected   = int[cities];
   
   is the solution. connected[i] contains the index of the facility 
   city i is connected to in the best solution found.

*/



// Call of volume algorithm
bool UNCAP_FACILITY_LOCATION_VOL(const double *open_cost,
				 const double *cost_matrix,
				 const int facilities,
				 const int cities,
				 int *connected,
				 double &cost);

#endif
