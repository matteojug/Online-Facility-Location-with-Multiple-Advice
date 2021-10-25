#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cplex.h>

// LP-solution 
bool LP_solution(const double *open_cost, const double *cost_matrix, const int facilities, 
		 const int cities, const bool ip, int *connected, double& cost)
{
  // This is the LP-solution for the underlying linear problem for uncapacitated
  // facility location. For further information on CPLEX refer to the CPLEX Manuals.

  int  cnt;
  bool val = false;
  int edges = facilities*cities;

  // Variables: Facilities yi, edges xij - overall facilities + edges
  int variables = edges + facilities;

  // Constraints: For every city a constraint, for every connection a constraint - overall: cities + edges
  int constraints = edges + cities;

  // 1*cities entries for every yi, 2 entries for every xij - overall 3*edges nonzero entries in the matrix
  int entries = 3*edges;

  double *obj = new double [variables];  // Objective Function
  double *rhs = new double [constraints];  // Right-Hand-Side
  char *sense = new char [constraints];   // Constraint type
  int *matbeg = new int [variables]; // Start of entries of column [j] at matval
  int *matcnt = new int [variables]; // Number of entries of column [j]
  int *matind = new int [entries];        // Constraint, in which the entry of matval appears
  double *matval = new double [entries];  // An array of 1's and -1's.
  double *lb = new double [variables]; // Lower bounds
  double *ub = new double [variables]; // Upper bounds
  // rngval is NULL

  // For MIP-optimization
  char *ctype = new char [variables];  // Variable type
  double *values = new double [variables]; // Values for init-solution
  int *indices = new int [variables]; // Indices for variables, that get initialized by init-solution

  // Solution
  double *x = new double [variables];

  // First we set up facilities
  for (cnt = 0; cnt < facilities; cnt++){
    matbeg[cnt] = cnt*cities; // Setting up indicator variables yi
    matcnt[cnt] = cities;
    obj[cnt] = open_cost[cnt];
    lb[cnt] = 0;
    ub[cnt] = CPX_INFBOUND;
  }

  // First the constraints: every city connected to 1 or more facilities
  for (cnt = 0; cnt < cities; cnt++) {
    rhs[cnt] = 1;
    sense[cnt] = 'G';
    connected[cnt] = -1;
  }

  // Now the constraints that every city is connected to an opened facility
  for (cnt = 0; cnt < edges; cnt++) {
    obj[facilities + cnt] = cost_matrix[cnt];
    rhs[cities + cnt] = 0;
    sense[cities + cnt] = 'G';
    matval[cnt] = 1;                         // Coefficients for yi
    matind[cnt] = cities+cnt;

    matbeg[facilities + cnt] = edges+(2*cnt); // Setting up variables
    matcnt[facilities + cnt] = 2;
    matval[edges+(2*cnt)] = 1;                // Coefficients for yij and yi
    matval[edges+(2*cnt)+1] = -1;
    matind[edges+(2*cnt)] = cnt % cities;
    matind[edges+(2*cnt)+1] = cities+cnt;      
    lb[facilities + cnt] = 0;                 // Bounds
    ub[facilities + cnt] = CPX_INFBOUND;
  }

  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int           status = 0;


  // Initialize the CPLEX environment

  env = CPXopenCPLEX (&status);

  if ( env == NULL ) {
    char  errmsg[1024];
    fprintf (stderr, "Could not open CPLEX environment.\n");
    CPXgeterrorstring (env, status, errmsg);
    fprintf (stderr, "%s", errmsg);
    goto TERMINATE_NO_CPX;
  }
  
  /* Turn on output to the screen 
  
  status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
  if ( status ) {
    fprintf (stderr, 
	     "Failure to turn on screen indicator, error %d.\n", status);
    goto TERMINATE;
  }
  
  // Turn on data checking

  status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
  if ( status ) {
    fprintf (stderr, 
	     "Failure to turn on data checking, error %d.\n", status);
    goto TERMINATE;
    }
  */

  // Create the problem. 
  
  lp = CPXcreateprob (env, &status, "lpex1");

  if ( lp == NULL ) {
    fprintf (stderr, "Failed to create LP.\n");
    goto TERMINATE;
  }

  // Set up variables

  status = CPXcopylp(env, lp, variables, constraints, CPX_MIN, obj, rhs, sense,
		     matbeg, matcnt, matind, matval, lb, ub, NULL);
  
  if (status != 0) {
    fprintf(stderr, "Error setting up LP.\n");
    goto TERMINATE;
  }

  if (!ip) {
    // This is the LP-version, factional output

    // write problem to a file
    // CPXwriteprob(env, lp, "lp.txt", "LP"); 

    // Solve problem
    status = CPXdualopt(env, lp);

    if (status != 0) {
      fprintf(stderr, "Error solving LP.\n");
      goto TERMINATE;
    }  

    // Get values
    status = CPXgetobjval(env, lp, &cost);
    status = CPXgetx(env, lp, x, 0, variables-1);

    if (status != 0) {
      fprintf(stderr, "Error retrieving solution.\n");
      goto TERMINATE;
    }  
    
    // If it is integer (enough), then we mark the corresponding connection
    for (cnt = facilities; cnt < variables; cnt++) 
      connected[(cnt-facilities) % cities] = (x[cnt] > 0.95) ? (cnt-facilities) / cities : -1;
    // This - of course - does not give us a complete solution !!

  } // end LP

  else { 
    // Here MIP-solution

    // Set up init solution and variable types
    for (cnt = 0; cnt < variables; cnt++) {
      indices[cnt] = cnt;
      values[cnt] = (cnt < (2*facilities)) ? 1 : 0;
      ctype[cnt] = (cnt < facilities) ? 'B' : 'C';
    }

    // Setting up variable types in the problem
    status = CPXcopyctype(env, lp, ctype);
    
    if (status != 0) {
      fprintf(stderr, "Error setting up variable types.\n");
      goto TERMINATE;
    }

    // Write problem to a file
    // CPXwriteprob(env, lp, "lp.txt", "LP"); 

    // Write init-solution
    status = CPXcopymipstart(env, lp, variables, indices, values);

    if (status != 0) {
      fprintf(stderr, "Error setting up feasible solution.\n");
      goto TERMINATE;
    }

    // Start the algorithm

    status = CPXmipopt(env, lp);

    if (status != 0) {
      fprintf(stderr, "Error solving Dual/MIP\n");
      goto TERMINATE;
    }

    // Get the solution values
    status = CPXgetmipobjval(env, lp, &cost);
    status = CPXgetmipx(env, lp, x, 0, variables-1);
    
    if (status != 0) {
      fprintf(stderr, "Error retrieving solution.\n");
      goto TERMINATE;
    }
    // cout << "Cost  : " << cost << endl;

    // set up solution
    for (cnt = facilities; cnt < variables; cnt++) 
      if (0.5 <= x[cnt]) {
	if (connected[(cnt-facilities) % cities] == -1)
	  connected[(cnt-facilities) % cities] = (cnt-facilities) / cities;
      }
  }

  val = true;

 TERMINATE:
  CPXfreeprob(env, &lp);
  CPXcloseCPLEX(&env);
 TERMINATE_NO_CPX:
  delete[] values;
  delete[] indices;
  delete[] x;
  delete[] lb;
  delete[] ub;
  delete[] matval;
  delete[] matcnt;
  delete[] matind;
  delete[] matbeg;
  delete[] obj;
  delete[] rhs;
  delete[] sense;
  delete[] ctype;
  
  return val;
}
