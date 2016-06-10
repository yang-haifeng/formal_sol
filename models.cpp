#include "models.h"

Model::Model(){
	r_max = 1e3*AU; 
	Kext=1; Kpol=0; Kcpol=0; 
	lambda=0.1;
	los_phi=0; los_theta=0;
	//Mext << Matrix4d(Vector4d::Constant(Kext).asDiagonal());
	//Vabs << Kext, Kpol, 0, 0;
}

double Model::get_BnuT(double x, double y, double z){
	cout<<"Warning: get_BnuT method not implemented!"<<endl;
	return 0;
}

double Model::get_Rho(double x, double y, double z){
	cout<<"Warning: get_Rho method not implemented!"<<endl;
	return 0;
}

void Model::get_Orientation(double x, double y, double z, double &theta, double &phi){
	cout<<"Warning: get_Orientation method not implemented!"<<endl;
	theta = 0; phi=0;
}

Vector4d Model::Integrate(double x, double y, double z, double n_theta, double n_phi, double step){
	Matrix4d T = Matrix4d(Vector4d::Constant(1).asDiagonal());
	double s = 0;
	Vector4d result = Vector4d::Constant(0);
	double rho, bnuT;
	double xp=x, yp=y, zp=z;
	double dx = - step * sin(n_theta) * cos(n_phi);
	double dy = - step * sin(n_theta) * sin(n_phi);
	double dz = - step * cos(n_theta);
	while (true){
		rho = get_Rho(xp, yp, zp);
		bnuT = get_BnuT(xp, yp, zp);

		result += T * cal_Vabs(xp, yp, zp, n_theta, n_phi) * step;
		T -= T * cal_Mext(xp, yp, zp, n_theta, n_phi) * step;

		xp += dx; yp += dy; zp += dz;
		if ( reachBoundary(xp, yp, zp) ) break;
	}
	return result;
}

void Model::test(){
	//cout<<Vabs<<endl;
	//cout<<Mext<<endl;
	cout<<cal_Mext(0,0,0,0,0)<<endl;
	Kpol=0.1;
	cout<<cal_Mext(0,0,0,0,0)<<endl;
	cout<<cal_Mext(0,0,0,PI/2,0)<<endl;
	Kpol=0;
}

bool Model::reachBoundary(double x, double y, double z){
	if (x*x + y*y + z*z >= r_max*r_max) return true;
	else return false;
}

Vector4d Model::cal_Vabs(double x, double y, double z, double n_theta, double n_phi){
	double cosDTheta;
	double theta, phi;
	get_Orientation(x, y, z, theta, phi);
	double dphi = phi-n_phi;
	//cosDTheta = sin(theta)*cos(phi)*sin(n_theta)*cos(n_phi) +
		//sin(theta)*sin(phi)*sin(n_theta)*sin(n_phi) + cos(theta)*cos(n_theta);
	cosDTheta = sin(theta)*cos(dphi)*sin(n_theta) + cos(theta)*cos(n_theta);
	Vector4d Vabs;
	Vabs << Kext + Kpol * cosDTheta*cosDTheta, 
	     Kpol * (1- cosDTheta*cosDTheta), 0, 0;
	return Vabs;
}

Matrix4d Model::cal_Mext(double x, double y, double z, double n_theta, double n_phi){
	double cosDTheta;
	double theta, phi;
	get_Orientation(x, y, z, theta, phi);
	double dphi = phi-n_phi;
	//cosDTheta = sin(theta)*cos(phi)*sin(n_theta)*cos(n_phi) +
		//sin(theta)*sin(phi)*sin(n_theta)*sin(n_phi) + cos(theta)*cos(n_theta);
	cosDTheta = sin(theta)*cos(dphi)*sin(n_theta) + cos(theta)*cos(n_theta);
	double Cep = Kext + Kpol * cosDTheta*cosDTheta;
	double Cpp = Kpol * (1- cosDTheta*cosDTheta);
	double Ccpp = Kcpol  * (1- cosDTheta*cosDTheta);

	double xp, yp; // The xp,yp after rotation, not the global coordinates.
	xp = sin(theta)*cos(dphi)*cos(n_theta) - cos(theta)*sin(n_theta);
	yp = sin(theta)*sin(dphi);
	double cos2phi = (xp*xp-yp*yp)/(xp*xp+yp*yp);
	double sin2phi = 2*xp*yp/(xp*xp+yp*yp);
	if ((xp*xp+yp*yp)==0) cos2phi=sin2phi=0;

	Matrix4d Mext;
	Mext << Cep, -Cpp*cos2phi, -Cpp*sin2phi, 0,
	     -Cpp*cos2phi, Cep, 0, Ccpp*sin2phi,
	     -Cpp*sin2phi, 0, Cep, -Ccpp*cos2phi,
	     0, -Ccpp*sin2phi, Ccpp*cos2phi, Cep;
	return Mext;
}
