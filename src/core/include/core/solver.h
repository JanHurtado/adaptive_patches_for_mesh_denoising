#ifndef SOLVER_H
#define SOLVER_H

#include <ilcplex/cplex.h>
#include "core/mesh.h"

void quadratic_programming_solver(int num_vars, AKMatrix & linear_term, AKMatrix & quadratic_term, AKMatrix & constraint_coefficients_diag, double single_constraint, vector<double> & solution);

#endif // SOLVER_H