#include "fit.h"
#include <math.h>

double hyperModel(double t, const double *p) {
	return 1.0/(p[0] + p[1] * t);
}

double expModel(double t, const double *p) {
	return 1.0/(p[0] + p[1]*(1 - exp(-t/p[2])));
}

/*
 * Fit HyperModel.
 * m: number of data point;
 * t: x
 * y: y
 */
int fitHyper(int n, double *par, int m, double *t, double *y) {
	lm_control_struct control = lm_control_double;
	lm_status_struct status;
	control.verbosity = 7;

	lmcurve(n, par, m, t, y, hyperModel, &control, &status);

	return 0;
}

int fitExp(int n, double *par, int m, double *t, double *y) {
	lm_control_struct control = lm_control_double;
	lm_status_struct status;
	control.verbosity = 7;

	lmcurve(n, par, m, t, y, expModel, &control, &status);

	return 0;
}