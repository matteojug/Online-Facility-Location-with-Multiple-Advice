LP-Solution for UFLP
-------------------- 
This function solves the underlying LP fractional or integral with CPLEX
 
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

The additional boolean parameter ip indicated whether the IP should be
solved (true) or the LP-relaxation (false). When ip = false is chosen, there
may not be a valid solution reported in connected.


Compile with CPLEX-include-path (needs standard cplex.h).


Martin Hoefer, 2002