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

	SlabSphGrain M1 = SlabSphGrain();
	cout<<M1.Image(0,0,10*AU, PI/4, 0);
	//cout<<get_psi(PI/2, 0, PI/2, 0)<<endl;
}
