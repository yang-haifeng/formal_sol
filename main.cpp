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
	double xmin = -50*AU, xmax = 50*AU;
	double ymin = -50*AU, ymax = 50*AU;
	double dx = (xmax-xmin)/Nx;
	double dy = (ymax-ymin)/Ny;
	double x,y;
	double deltaX = 10*AU * tan(theta);
	ofstream Fout;
	Fout.open("test/RDisk10x10.txt");
	Vector4d result;
	Fout<<"#x\ty\tI\tQ\tU\tV"<<endl;
	for (int i=0; i<Nx+1; i++){
	for (int j=0; j<Ny+1; j++){
		cout<<"Working on: "<<i<<", "<<j<<endl;
		x = xmin + i*dx + deltaX;
		y = ymin + j*dy;
		result = M.Image(x, y, 10*AU, PI/4., 0.);
	  	Fout<<x<<y<<"\t"<<result[0]<<"\t"<<result[1]<<"\t"<<result[2]<<"\t"<<result[3]<<endl;
	}
	}
	*/

	// Here's the HL Tau model
	HLTau M = HLTau();
	M.set_adaptive(0.01);
	//cout<<M.Image(10*AU,10*AU,100*AU,PI/4.,0)<<endl;
	M.get_Image(PI/4, 100, 100*AU, "test/45.out");
}




