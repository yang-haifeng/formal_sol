#include "models.h"

/////////////////////////////////////////////////////////////////////////////
// Member methods for Model class.

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
	Vector4d Vabs; Matrix4d Mext;
	while (true){
		rho = get_Rho(xp, yp, zp);
		bnuT = get_BnuT(xp, yp, zp);

		cal_VM(xp, yp, zp, n_theta, n_phi, Vabs, Mext);
		result += T * Vabs * step*rho*bnuT;
		T -= T * Mext * step*rho;

		xp += dx; yp += dy; zp += dz;
		if ( reachBoundary(xp, yp, zp) ) break;
	}
	return result;
}

void Model::test(){
	//cout<<Vabs<<endl;
	//cout<<Mext<<endl;
}

bool Model::reachBoundary(double x, double y, double z){
	if (x*x + y*y + z*z >= r_max*r_max) return true;
	else return false;
}

void Model::cal_VM(double x, double y, double z, double n_theta, double n_phi, 
		Vector4d &Vout, Matrix4d &Mout){
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

	Vout << Cep, -Cpp*cos2phi, -Cpp*sin2phi, 0;
	Mout << Cep, -Cpp*cos2phi, -Cpp*sin2phi, 0,
	     -Cpp*cos2phi, Cep, 0, Ccpp*sin2phi,
	     -Cpp*sin2phi, 0, Cep, -Ccpp*cos2phi,
	     0, -Ccpp*sin2phi, Ccpp*cos2phi, Cep;
}

/////////////////////////////////////////////////////////////////////////////
// Member methods for SlabUniform class.

SlabUniform::SlabUniform(){
	r_max = 100*AU; h = 10*AU;
	Kext = 1.0; Kpol = 0.1; Kcpol = 0.01;

	lambda = 0.1; los_theta = 0; los_phi=0; 

	Bnu0 = BnuT(30., con_c/lambda);
	rho0 = 1e-15;
}

SlabUniform::SlabUniform(double Temp, double rho, double height, double rm){
	r_max = rm; h = height;
	Kext = 1.0; Kpol = 0.1; Kcpol = 0.01;

	lambda = 0.1; los_theta = 0; los_phi=0; 

	Bnu0 = BnuT(Temp, con_c/lambda);
	rho0 = rho;
}

double SlabUniform::get_BnuT(double x, double y, double z){
	if (fabs(z)<=h)
		//return Bnu0;
		return 1.; // For testing purpose, BnuT is set as 1;
	else return 0;
}

double SlabUniform::get_Rho(double x, double y, double z){
	if (fabs(z)<=h)
		return rho0;
	else return 0;
}

void SlabUniform::get_Orientation(double x, double y, double z, double &theta, double &phi){
	theta = PI/2; phi=0;
}

bool SlabUniform::reachBoundary(double x, double y, double z){
	if (fabs(z)>1.1*h) return true;
	if (x*x+y*y>r_max*r_max) return true;
	return false;
}

/////////////////////////////////////////////////////////////////////////////
// Member methods for SlabTwisted class.

SlabTwisted::SlabTwisted(){
	r_max = 100*AU; h = 10*AU;
	Kext = 1.0; Kpol = 0.1; Kcpol = 0.01;

	lambda = 0.1; los_theta = 0; los_phi=0; 

	Bnu0 = BnuT(30., con_c/lambda);
	rho0 = 1e-15;

	k = 10.;
}

SlabTwisted::SlabTwisted(double Temp, double rho, double height, double rm, double K){
	r_max = rm; h = height;
	Kext = 1.0; Kpol = 0.1; Kcpol = 0.01;

	lambda = 0.1; los_theta = 0; los_phi=0; 

	Bnu0 = BnuT(Temp, con_c/lambda);
	rho0 = rho;

	k = K;
}

void SlabTwisted::get_Orientation(double x, double y, double z, double &theta, double &phi){
	theta = PI/2; 
	phi = k*z*rho0;
}

void SlabTwisted::test(){
	double z;
	double theta, phi;
	for (int i=0; i<21; i++){
		z = i-10;
		get_Orientation(0,0,z*AU, theta, phi);
		cout<<z<<"\t"<<phi/PI*180<<endl;
	}
}

/////////////////////////////////////////////////////////////////////////////
// Miscellaneous functions

double BnuT(double T, double nu){
	return 2*con_h*pow(nu,3)/(con_c*con_c) / (exp(con_h*nu/con_k/T) - 1);
}

