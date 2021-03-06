#include <iostream>
#include <fstream>
#include "typedef.h"
#include "models.h"
#include "utils.h"
using namespace std;

int main(){
	// Testing No.1: I = 1-e^(-tau);
	/*
	Vector4d Vt;
	for (int i=0; i<10; i++){
		SlabUniform M1(30, 3.33e-16*(i+1), 10*AU, 100*AU);
		Vt<<M1.Integrate(0,0,10*AU, 0, 0);
		cout<<Vt[0]<<"\t"<<Vt[1]<<endl;
	}
	*/


	// This simply tests the orientation setup
	//SlabTwisted M1 = SlabTwisted();
	//M1.test();
	
	// This test against Martin 1974
	/*
	Vector4d Vt;
	double k;
	int N=180;
	for (int i=0; i<N; i++){
		k = PI/N*i/(3.33e-16*20*AU);
		SlabTwisted M1(30, 3.33e-16, 10*AU, 100*AU, k);
		Vt<<M1.Integrate(0,0,10*AU, 0, 0);
		//cout<<Vt[0]<<"\t"<<Vt[1]<<endl;
		cout<<180./N*i<<"\t"<<Vt[0]<<"\t"<<Vt[1]<<"\t"<<Vt[2]<<"\t"<<Vt[3]<<endl;
	}
	*/

	// Tests for get_psi
	/*
	double theta;
	cout<<"Test1"<<endl;
	for (int i=0; i<10; i++){
		theta = PI/10*i;
	  	cout<<theta<<" "<<get_psi(0, 0, theta, 0)<<endl;
	}
	cout<<"Test2"<<endl;
	for (int i=0; i<10; i++){
		theta = PI/10*i;
	  	cout<<theta<<" "<<get_psi(theta, 0, PI/2, PI/2)<<endl;
	}
	cout<<"Test3"<<endl;
	double phi;
	for (int i=0; i<10; i++){
		phi = PI/10*i;
	  	cout<<phi<<" "<<get_psi(PI/2, 0, PI/2, phi)<<endl;
	}
	*/

	// Tests for to_grain_frame
	/*
	double btheta = 1., bphi = 0;
	double theta, phi;
	double thetap, phip;
	for (int i=0; i<10; i++){
		theta = i/10.*PI;
	//for (int j=0; j<10; j++){
		//phi = i/10.*PI*2;
		phi = 0.;
		to_grain_frame(btheta, bphi, theta, phi, thetap, phip);
		cout<<theta<<"\t"<<thetap<<"\t"<<phi<<"\t"<<phip<<endl;
	}
	//}
	*/

	// Test for Z_Matrix from Rayleigh Limit.
	/*
	Model M = Model();
	Matrix4d Z;
	for (int i=0; i<180; i++){
		Z = M.get_ZMatrix(0., 0., i/180.*PI, 0.);
		cout << i <<" "<< -Z(0,1)/Z(0,0)<<endl;
	}
	*/

	// The following test is the test for Rayleigh limit with diff incl.
	/*
	SlabSphGrain M1 = SlabSphGrain(30, 1e-16, 10*AU, 200*AU);
	Vector4d result;
	double theta;
	ofstream Fout;
	Fout.open("test/RTv1.txt");
	for (int i=0;i<10;i++){
	  theta = i/10.*PI/2;
	  cout<<"Working on: l_theta = "<<theta<<endl;
	  result = M1.Image(0, 0, 10*AU, theta, 0);
	  //result = M1.Integrate(0, 0, 10*AU, theta, 0);
	  Fout<<theta<<" "<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<endl;
	}
	*/

	// Rayleigh Test with different phi angle at theta = PI/4.
	/*
	SlabSphGrain M1 = SlabSphGrain(30, 1e-15, 10*AU, 200*AU);
	Vector4d result;
	double phi;
	ofstream Fout;
	Fout.open("test/RTv2.txt");
	for (int i=0; i<36; i++){
	  phi = i*10./180.*PI;
	  cout<<"Working on: l_phi = "<<phi<<endl;
	  result = M1.Image(0, 0, 10*AU, PI/4., phi);
	  //result = M1.Integrate(0, 0, 10*AU, theta, 0);
	  Fout<<phi<<" "<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<endl;
	}
	*/

	//SlabSphGrain M1 = SlabSphGrain(30, 1e-15, 10*AU, 200*AU);
	//cout<<M1.Image(0,0,10*AU,PI/4,0)<<endl;

	// Let's generate a full image.
	/*
	SlabSphGrain M = SlabSphGrain(30, 1e-15, 10*AU, 200*AU);
	int Nx = 10, Ny = 10;
	double theta = PI/4.;
	double xmin = -100*AU, xmax = 100*AU;
	double ymin = -100*AU, ymax = 100*AU;
	double dx = (xmax-xmin)/Nx;
	double dy = (ymax-ymin)/Ny;
	double x,y;
	double deltaX = 10*AU * tan(theta);
	ofstream Fout;
	Fout.open("test/Slab45d.txt");
	Vector4d result;
	Fout<<"#x\ty\tI\tQ\tU\tV"<<endl;
	for (int i=0; i<Nx+1; i++){
	for (int j=0; j<Ny+1; j++){
		cout<<"Working on: "<<i<<", "<<j<<endl;
		x = xmin + i*dx + deltaX;
		y = ymin + j*dy;
		result = M.Image(x, y, 10*AU, PI/4, 0.);
	  	Fout<<x/AU<<"\t"<<y/AU<<"\t"<<result[0]<<"\t"<<result[1]<<"\t"<<result[2]<<"\t"<<result[3]<<endl;
	}
	}
	*/

	// This is test for the weird angle problem.
	/*
	SlabSphGrain M = SlabSphGrain(30, 1e-15, 10*AU, 200*AU);
	double y, theta;
	ofstream Fout;
	Fout.open("test/RDisk_Chenliang.txt");
	Vector4d result;
	Fout<<"#theta\ty\tI\tQ\tU\tV"<<endl;

	y = 100*AU; theta = PI/6;
	result = M.Image(y, 0, 10*AU, theta, PI/2.);
	cout<<theta/PI*180<<"\t"<<y/AU<<"\t"<<result[0]<<"\t"<<result[1]<<"\t"<<result[2]<<"\t"<<result[3]<<endl;

	y =  50*AU; theta = PI/4;
	result = M.Image(y, 0, 10*AU, theta, PI/2.);
	Fout<<theta/PI*180<<"\t"<<y/AU<<"\t"<<result[0]<<"\t"<<result[1]<<"\t"<<result[2]<<"\t"<<result[3]<<endl;

	y = 100*AU; theta = PI/6;
	result = M.Image(y, 0, 10*AU, theta, PI/2.);
	Fout<<theta/PI*180<<"\t"<<y/AU<<"\t"<<result[0]<<"\t"<<result[1]<<"\t"<<result[2]<<"\t"<<result[3]<<endl;

	y =  50*AU; theta = PI/6;
	result = M.Image(y, 0, 10*AU, theta, PI/2.);
	Fout<<theta/PI*180<<"\t"<<y/AU<<"\t"<<result[0]<<"\t"<<result[1]<<"\t"<<result[2]<<"\t"<<result[3]<<endl;
	*/

	// Here's the HL Tau model
	/*
	HLTau M = HLTau();
	//M.set_kappa(0.3899+0.418e-2, 0.4178e-2);
	M.set_kappa(0.553521E+00 + 0.469887E+01, 0.469887E+01);
	M.set_adaptive(0.1);
	//cout<<M.Image(10*AU,0,0,0.,0)<<endl;
	cout<<M.Image(10*AU,0,0,0.,0)<<endl;
	*/
	//M.multiply_rho0();
	//M.multiply_rho0(10);
	//cout<<M.Image(10*AU,10*AU,100*AU,PI/4.,0)<<endl;
	//M.get_Image(PI/4, 30, 200*AU, "test/45degree_30x30_200au_doublerho0.out");
	//M.restart_Circle_Image(PI/4., 15, 16, 10*AU, 150*AU, "test/HLTau_vo_a100micron.out", 127);

	// Generating Density plotting here.
	/*
	HLTau M = HLTau();
	ofstream FHR;
	FHR.open("test/HLTau_HR.out");
	for (int i=0; i<200; i++){
		FHR<<i<<"\t"<<M.get_HR(i*AU)/AU<<endl;
	}
	*/

	// This doesn't seems to be useful anymore.
	/*
	HLTau M = HLTau();
	M.set_adaptive(0.005);
	double theta = PI/4;
	int Npx = 30;
	double FoV = 200*AU;
        Vector3d de1, de2, de3;
        de1 << cos(theta), 0, -sin(theta);
        de2 << 0, 1, 0;
        de3 << sin(theta), 0, cos(theta);
        de1 *= FoV/(Npx-1);
        de2 *= FoV/(Npx-1);
        de3 *= AU; 
        Vector3d P0; 
        P0 << 0,0,0;
        P0 -= de1 * (Npx/2-0.5);
        P0 -= de2 * (Npx/2-0.5);
        Vector3d P;
        Vector4d result;
        //for (int i=0;i<Npx; i++){
        //for (int j=0;j<Npx; j++){
	int i=15, j=28;
          cout<<i<<"\t"<<j<<endl;
          P = P0 + de1 *i + de2*j;
          cout<<P(0)/AU<<"\t"<<P(1)/AU<<"\t"<<P(2)/AU<<"\t"<<endl;
          while (!M.reachBoundary(P+de3)){
                P+=de3;
          }
          cout<<P(0)/AU<<"\t"<<P(1)/AU<<"\t"<<P(2)/AU<<"\t"<<endl;
          result = M.Image(P(0), P(1), P(2), theta, 0); 
          cout<<result(0)<<"\t"<<result(1)<<"\t"<<result(2)<<"\t"<<result(3)<<endl;
        //}
        //}
	*/

	// Some test for the ZMatrix
	// This is the case for xOy plane scattering
	/*
	Model M = Model();
	Matrix4d Z;
	for (int i=0; i<360; i++)
	{
		Z = M.get_ZMatrix(PI/2, 0, PI/2, i/180.*PI);
		cout<<i<<"\t"<<Z(0,0)<<"\t"<<Z(0,1)<<"\t"<<Z(0,2)<<"\t"<<Z(1,0)<<"\t"<<Z(1,1)<<"\t"<<Z(1,2)<<"\t"<<Z(2,0)<<"\t"<<Z(2,1)<<"\t"<<Z(2,2)<<"\t"<<Z(3,3)<<endl;
	}
	*/
	// Inclined disk planar scattering
	/*
	Model M = Model();
	Matrix4d Z;
	for (int i=0; i<360; i++)
	{
		Z = M.get_ZMatrix(PI/2, i/180.*PI, PI/4, 0);
		cout<<i<<"\t"<<Z(0,0)<<"\t"<<Z(0,1)<<"\t"<<Z(0,2)<<"\t"<<Z(1,0)<<"\t"<<Z(1,1)<<"\t"<<Z(1,2)<<"\t"<<Z(2,0)<<"\t"<<Z(2,1)<<"\t"<<Z(2,2)<<"\t"<<Z(3,3)<<endl;
	}
	*/

	// Let's do some ConeModel.
	// Start with reachBoundary test and density, temperature.
	/*
	double x, z;
	ConeModel M = ConeModel(PI/4, 1e-15);
	for (int i=0; i<101; i++){
	for (int j=0; j<101; j++){
		x = (-200+4*i)*AU;
		z = (-200+4*j)*AU;
		if (!M.reachBoundary(x, 0., z)){
			cout<<x/AU<<"\t"<<z/AU<<"\t"<<M.get_BnuT(x, 0., z)<<"\t"<<M.get_Rho(x, 0., z)<<endl;
		}
	}
	}
	*/
	// Now onto some business.
	/*
	ConeModel M = ConeModel(PI/6, 1e-14);
	M.set_adaptive(0.1); // This might be too big of an adaptive. 
	M.get_Circle_Image(PI/4, 15, 16, 10*AU, 150*AU, "test/C10x_45d.out");
	*/

	// Need some tests for BnuT.
	/*
	double lambda, nu;
	double T0=5000;
	for (int i=0; i<100; i++){
	  lambda = 0.03*(i+1)/1e4;
	  nu = con_c/lambda;
	  cout<<lambda<<" "<<BnuT(T0,nu)<<endl;
	}
	*/

	// Check for HL Tau Temperature structure
	/*
	HLTau H = HLTau();
	double a = 1.5625*AU;
	cout<<H.get_BnuT(a,a,a)<<endl;
	double r;
	for (int i=0; i<100; i++){
	  r = 2.0*(i+1)*AU;
	  cout<<2*(i+1)<<" "<<H.get_BnuT(r,0,0)/AU<<endl;
	}
	*/

	// Let's use test to calculate polarization degree along minor axis
	// for HL Tau model
	//HLTau M = HLTau();
	//M.set_kappa(7.381826e-01+6.350236, 6.350236);
	//M.set_kappa(5.176868e-01+5.640267e-03, 5.640267e-03);
	//M.set_adaptive(0.1);
	//M.multiply_rho0(0.1);
	//M.get_Image_Minor(PI/4, 10, 10*AU, 150*AU, "test/minor_original.dat");

	// Uniform Slab model to compare with analytical work
	double H = 50*AU;
	SlabSphGrain SSG = SlabSphGrain(30, 1e-14, H, 200*AU);
	SSG.set_Kappa(7.381826e-01+6.350236, 6.350236);
	double inc = PI/4;
	cout<<SSG.Image(H*tan(inc), 0, H, inc, 0);

}



