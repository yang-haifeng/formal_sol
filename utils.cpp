#include "utils.h"

double get_psi(double btheta, double bphi, double theta, double phi){
	double dphi = phi-bphi;
	double tanpsi = sin(btheta) * sin(dphi)/(cos(theta) * sin(btheta) * cos(dphi) - sin(theta) * cos(btheta));
	return atan(tanpsi);
}

void to_grain_frame(double btheta, double bphi, double theta, double phi, 
		double &thetap, double &phip){
	double dphi = phi-bphi;
	double costheta = sin(btheta) * sin(theta) * cos(dphi) + cos(btheta) * cos(theta);
	thetap = acos(costheta);
	phip = atan2(sin(theta) * sin(phi), cos(btheta) * sin(theta) * cos(dphi) - sin(btheta) * cos(theta) );
}

Matrix4d rotation_Matrix(double psi){
	Matrix4d result;
	result << 1, 0, 0, 0,
	       0, cos(2*psi), -sin(2*psi), 0,
	       0, sin(2*psi),  cos(2*psi), 0,
	       0, 0, 0, 1;
	return result;
}
