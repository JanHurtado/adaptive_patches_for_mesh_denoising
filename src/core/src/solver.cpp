#include "core/solver.h"

void quadratic_programming_solver(int num_vars, AKMatrix & linear_term, AKMatrix & quadratic_term, AKMatrix & constraint_coefficients_diag, double single_constraint, vector<double> & solution)
{
	if (num_vars > 0)
	{
		char     probname[] = "quadratic programming solver";
		int      numcols = num_vars;
		int      numrows = 1;
		int      objsen = 1;
		AKNumber   *obj = new AKNumber[numcols];
		AKNumber   *rhs = new AKNumber[numrows];
		char     *sense = new char[numrows];
		int      *matbeg = new int[numcols];
		int      *matcnt = new int[numcols];
		int      *matind = new int[numcols];
		AKNumber   *matval = new AKNumber[numcols];
		AKNumber   *lb = new AKNumber[numcols];
		AKNumber   *ub = new AKNumber[numcols];
		int      *qmatbeg = new int[numcols];
		int      *qmatcnt = new int[numcols];
		int      *qmatind = new int[numcols * numcols];
		AKNumber   *qmatval = new AKNumber[numcols * numcols];

		/* Declare and allocate space for the variables and arrays where we
		will store the optimization results including the status, objective
		value, variable values, dual values, row slacks and variable
		reduced costs. */

		int      solstat;
		AKNumber   objval;
		AKNumber * x = new AKNumber[numcols];
		AKNumber * pi = new AKNumber[numrows];
		AKNumber * slack = new AKNumber[numrows];
		AKNumber * dj = new AKNumber[numcols];


		CPXENVptr     env = NULL;
		CPXLPptr      lp = NULL;
		int           status;
		int           cur_numrows, cur_numcols;
		int i, j;


		env = CPXopenCPLEX(&status);

		/* If an error occurs, the status value indicates the reason for
		failure.  A call to CPXcomputeerrorstring will produce the text of
		the error message.  Note that CPXopenCPLEX produces no output,
		so the only way to see the cause of the error is to use
		CPXcomputeerrorstring.  For other CPLEX routines, the errors will
		be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */

		if (env == NULL) {
			char  errmsg[CPXMESSAGEBUFSIZE];
			fprintf(stderr, "Could not open CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
			goto TERMINATE;
		}
		//status = CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
		status = CPXsetintparam(env, CPXPARAM_OptimalityTarget, CPX_OPTIMALITYTARGET_FIRSTORDER);
		if (status) {
			fprintf(stderr,
				"Failure to turn on screen indicator, error %d.\n", status);
			goto TERMINATE;
		}
		/* Fill in the data for the problem.  */

		rhs[0] = single_constraint;
		sense[0] = 'E';
		int k = 0;
		for (i = 0; i < num_vars; i++)
		{
			lb[i] = 0.0f;
			ub[i] = 1.0f;
			obj[i] = linear_term(i, 0);
			matbeg[i] = i;
			matcnt[i] = 1;
			matind[i] = 0;
			matval[i] = constraint_coefficients_diag(i, i);
			qmatbeg[i] = i * num_vars;
			qmatcnt[i] = num_vars;
			for (j = 0; j < num_vars; j++)
			{
				qmatind[k] = j;
				qmatval[k] = quadratic_term(j, i);
				k++;
			}

		}


		/* Create the problem. */

		lp = CPXcreateprob(env, &status, probname);

		/* A returned pointer of NULL may mean that not enough memory
		was available or there was some other problem.  In the case of
		failure, an error message will have been written to the error
		channel from inside CPLEX.  In this example, the setting of
		the parameter CPXPARAM_ScreenOutput causes the error message to
		appear on stdout.  */

		if (lp == NULL) {
			fprintf(stderr, "Failed to create problem.\n");
			goto TERMINATE;
		}

		/* Now copy the LP part of the problem data into the lp */

		status = CPXcopylp(env, lp, numcols, numrows, objsen, obj, rhs,
			sense, matbeg, matcnt, matind, matval,
			lb, ub, NULL);

		if (status) {
			fprintf(stderr, "Failed to copy problem data.\n");
			goto TERMINATE;
		}

		status = CPXcopyquad(env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
		if (status) {
			fprintf(stderr, "Failed to copy quadratic matrix.\n");
			cout << "status " << status << endl;
			goto TERMINATE;
		}


		/* Optimize the problem and obtain solution. */
		status = CPXqpopt(env, lp);
		//status = CPXmipopt(env, lp);
		if (status) {
			fprintf(stderr, "Failed to optimize QP.\n");
			goto TERMINATE;
		}

		status = CPXsolution(env, lp, &solstat, &objval, x, pi, slack, dj);
		if (status) {
			fprintf(stderr, "Failed to obtain solution.\n");
			goto TERMINATE;
		}


		/* Write the output to the screen. */

		//printf("%d ", solstat);
		//printf("\nSolution status = %d\n", solstat);
		//printf("Solution value  = %f\n\n", objval);

		/* The size of the problem should be obtained by asking CPLEX what
		the actual size is, rather than using what was passed to CPXcopylp.
		cur_numrows and cur_numcols store the current number of rows and
		columns, respectively.  */

		cur_numrows = CPXgetnumrows(env, lp);
		cur_numcols = CPXgetnumcols(env, lp);
		for (i = 0; i < cur_numrows; i++) {
			//printf("Row %d:  Slack = %10f  Pi = %10f\n", i, slack[i], pi[i]);
		}

		solution.clear();
		solution.resize(cur_numcols);
		for (j = 0; j < cur_numcols; j++) {
			cur_numcols;
			solution[j] = x[j];
			//_patch_function.push_back(pair<size_t, AKNumber>(_distance_function[j].first, x[j]));
			//printf("Column %d:  Value = %10f  Reduced cost = %10f\n", j, x[j], dj[j]);
		}

		/* Finally, write a copy of the problem to a file. */

	TERMINATE:

		if (lp != NULL) {
			status = CPXfreeprob(env, &lp);
			if (status) {
				fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
			}
		}

		/* Free up the CPLEX environment, if necessary */

		if (env != NULL) {
			status = CPXcloseCPLEX(&env);

			/* Note that CPXcloseCPLEX produces no output,
			so the only way to see the cause of the error is to use
			CPXcomputeerrorstring.  For other CPLEX routines, the errors will
			be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON. */

			if (status) {
				char  errmsg[CPXMESSAGEBUFSIZE];
				fprintf(stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring(env, status, errmsg);
				fprintf(stderr, "%s", errmsg);
			}
		}

		delete[] obj;
		delete[] rhs;
		delete[] sense;
		delete[] matbeg;
		delete[] matcnt;
		delete[] matval;
		delete[] lb;
		delete[] ub;
		delete[] qmatbeg;
		delete[] qmatcnt;
		delete[] qmatind;
		delete[] qmatval;

		delete[] x;
		delete[] pi;
		delete[] slack;
		delete[] dj;
	}
}
