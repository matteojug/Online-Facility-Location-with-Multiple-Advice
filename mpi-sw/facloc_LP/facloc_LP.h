// facloc_LP.h
#ifndef _FACLOC_LP
#define _FACLOC_LP

/* LP-Solution - solves the underlying LP fractional or integral with CPLEX

   UFLP-call interface:
   --------------------
   Please provide the following arrays:	

   open_cost   = double[facilities];
   cost_matrix = double[facilities*cities];

   All data non-negative. The connection cost from city j to facility i
   is denoted in position cost_matrix[(i*cities) + j]
  
   connected   = int[cities];
   
   is the solution. connected[i] and contains the index of the facility 
   city i is connected to in the best solution found.
*/
  
bool LP_solution(const double *open_cost,
	       const double *cost_matrix,
	       const int facilities,
	       const int cities,
	       const bool ip,  // true - solves MIP, false - solves LP
		     int *connected,
		     double& cost);      

#endif
