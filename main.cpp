#include <iostream>
#include <fstream>
#include "typedef.h"
#include "models.h"
#include "utils.h"
using namespace std;

int main(){
	/*
	Vector4d Vt;
	for (int i=0; i<10; i++){
		SlabUniform M1(30, 3.33e-16*(i+1), 10*AU, 100*AU);
		Vt<<M1.Integrate(0,0,10*AU, 0, 0);
		cout<<Vt[0]<<"\t"<<Vt[1]<<endl;
	}
	*/
	//SlabTwisted M1 = SlabTwisted();
	//M1.test();
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

	/*
	*/
	SlabSphGrain M1 = SlabSphGrain();
	Vector4d result;
	double theta;
	ofstream Fout;
	Fout.open("test/out2");
	for (int i=0;i<10;i++){
	  theta = i/10.*PI/2;
	  cout<<"Working on: l_theta = "<<theta<<endl;
	  result = M1.Image(0, 0, 10*AU, theta, 0);
	  Fout<<theta<<" "<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<endl;
	}
	//cout<<get_psi(PI/2, 0, PI/2, 0)<<endl;
	//
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

	//SlabSphGrain M1 = SlabSphGrain();
	//cout<< M1.Image(0, 0, 10*AU, 0., 0.)<<endl;
}
