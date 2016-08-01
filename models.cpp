#include "models.h"
#include "utils.h"

/////////////////////////////////////////////////////////////////////////////
// Member methods for Model class.

Model::Model(){
	r_max = 1e3*AU; 
	Kext=1; Kpol=0; Kcpol=0; 
	lambda=0.1;
	Ksca=0.1;
	//los_phi=0; los_theta=0;
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

Matrix4d Model::get_ZMatrix(double theta_i, double phi_i, double theta_o, double phi_o){
	double dphi = phi_o-phi_i;
	double cosTheta = sin(theta_i) * sin(theta_o) * cos(dphi) + cos(theta_i) * cos(theta_o);
	double factor = 3./16./PI * Ksca;
	Matrix4d M0;
	double cosTheta2 = cosTheta*cosTheta;
	M0 << factor * (1+cosTheta2), -factor * (1-cosTheta2), 0, 0,
	   -factor * (1-cosTheta2), factor * (1+cosTheta2), 0, 0,
	   0, 0, 2*cosTheta, 0,
	   0, 0, 0, 2*cosTheta;

	double cosi1 = -sin(theta_o)*cos(theta_i) * cos(dphi) + sin(theta_i)*cos(theta_o);
	double cosi2 = -sin(theta_i)*cos(theta_o) * cos(dphi) + sin(theta_o)*cos(theta_i);
	double i1 = acos(cosi1);
	double i2 = acos(cosi2);
	return rotation_Matrix(-i2)*M0*rotation_Matrix(i1);
}

Vector4d Model::Integrate(double x, double y, double z, double n_theta, double n_phi, double step){
	Matrix4d T = Matrix4d(Vector4d::Constant(1).asDiagonal());
	//double s = 0;
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

Vector4d Model::Image(double x, double y, double z, double l_theta, double l_phi, double step){
	Matrix4d T = Matrix4d(Vector4d::Constant(1).asDiagonal());
	Vector4d result = Vector4d::Constant(0);
	double rho, bnuT;
	double xp=x, yp=y, zp=z;
	double dx = -step*sin(l_theta)*cos(l_phi);
	double dy = -step*sin(l_theta)*sin(l_phi);
	double dz = -step*cos(l_theta);
	Vector4d Vabs; Matrix4d Mext;
	Matrix4d Z;
	double theta, phi;
	Vector4d Sin, Ssca;
	Matrix4d R1, R2;
	double psi1, psi2;
	double btheta, bphi;
	double theta1, phi1;
	double theta2, phi2;
	while (true){
		//cout<<"***********************************************"<<endl;
		//cout<<xp/AU<<"\t"<<yp/AU<<"\t"<<zp/AU<<endl;
		rho = get_Rho(xp, yp, zp);
		bnuT = get_BnuT(xp, yp, zp);

		cal_VM(xp, yp, zp, l_theta, l_phi, Vabs, Mext);

		Ssca = Vector4d::Constant(0);
		/*
		*/
		get_Orientation(xp, yp, zp, btheta, bphi);
		to_grain_frame(btheta, bphi, l_theta, l_phi, theta2, phi2);
		psi2 = get_psi(btheta, bphi, theta, phi);
		for (int i=0; i<Ntheta; i++){
			theta = PI/Ntheta*i;
			for (int j=0; j<Nphi; j++){
				phi = 2*PI/Nphi*j;
				Sin = Integrate(xp, yp, zp, theta, phi);

				//cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"<<endl;
				//cout<<theta<<"\t"<<phi<<endl;
				//cout<<Sin<<endl;

				psi1 = get_psi(btheta, bphi, theta, phi);
				to_grain_frame(btheta, bphi, theta, phi, theta1, phi1);

				//cout<<psi2<<" rotation Matrix:"<<endl;
				//cout<<rotation_Matrix(-psi2)<<endl;
				//cout<<"Z Matrix:"<<endl;
				//cout<<get_ZMatrix(theta1, phi1, theta2, phi2)<<endl;
				//cout<<psi1<<" rotation Matrix:"<<endl;
				//cout<<rotation_Matrix(psi1)<<endl;

				Ssca += rotation_Matrix(-psi2) * 
					get_ZMatrix(theta1, phi1, theta2, phi2) * 
					rotation_Matrix(psi1) *
					Sin * 
					rho * sin(theta) * PI/Ntheta * 2*PI/Nphi;
				//cout<<Ssca<<endl;
			}
		}

		result += (Vabs * rho * bnuT + Ssca)*step;
		//cout<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<endl;

		T -= T*Mext * step * rho;
		xp += dx; yp += dy; zp += dz;
		if (reachBoundary(xp, yp, zp)) break;
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

	lambda = 0.1; //los_theta = 0; los_phi=0; 

	Bnu0 = BnuT(30., con_c/lambda);
	rho0 = 1e-15;
}

SlabUniform::SlabUniform(double Temp, double rho, double height, double rm){
	r_max = rm; h = height;
	Kext = 1.0; Kpol = 0.1; Kcpol = 0.01;

	lambda = 0.1; //los_theta = 0; los_phi=0; 

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

	lambda = 0.1; //los_theta = 0; los_phi=0; 

	Bnu0 = BnuT(30., con_c/lambda);
	rho0 = 1e-15;

	k = 10.;
}

SlabTwisted::SlabTwisted(double Temp, double rho, double height, double rm, double K){
	r_max = rm; h = height;
	Kext = 1.0; Kpol = 0.1; Kcpol = 0.01;

	lambda = 0.1; //los_theta = 0; los_phi=0; 

	Bnu0 = BnuT(Temp, con_c/lambda);
	rho0 = rho;

	k = K;
}

void SlabTwisted::get_Orientation(double x, double y, double z, double &theta, double &phi){
	theta = PI/2; 
	phi = k*(z+h)*rho0;
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
// Member methods for SlabUniform class.

SlabSphGrain::SlabSphGrain(){
	r_max = 100*AU; h = 10*AU;
	Kext = 1.0; Kpol = 0.0; Kcpol = 0.0;
	Ksca = 0.1;

	lambda = 0.1; //los_theta = 0; los_phi=0; 

	Bnu0 = BnuT(30., con_c/lambda);
	rho0 = 1e-15;
}

SlabSphGrain::SlabSphGrain(double Temp, double rho, double height, double rm){
	r_max = rm; h = height;
	Kext = 1.0; Kpol = 0.0; Kcpol = 0.0;
	Ksca = 0.1;

	lambda = 0.1; //los_theta = 0; los_phi=0; 

	Bnu0 = BnuT(Temp, con_c/lambda);
	rho0 = rho;
}

void SlabSphGrain::get_Orientation(double x, double y, double z, double &theta, double &phi){
	//double epsilon = 0.05;
	theta = PI/2; 
	//phi = 0;
	phi = PI/Nphi;
}

/////////////////////////////////////////////////////////////////////////////
// Miscellaneous functions

double BnuT(double T, double nu){
	return 2*con_h*pow(nu,3)/(con_c*con_c) / (exp(con_h*nu/con_k/T) - 1);
}

