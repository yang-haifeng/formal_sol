#ifndef _UTIL_H
#define _UTIL_H

#include <math.h>
#include "typedef.h"
#include "Eigen/Dense"

using namespace Eigen;

double get_psi(double btheta, double bphi, double theta, double phi);

void to_grain_frame(double btheta, double bphi, double theta, double phi, 
		double &thetap, double &phip);

Matrix4d rotation_Matrix(double psi);

#endif
