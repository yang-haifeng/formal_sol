#include <iostream>
#include <fstream>
#include "typedef.h"
#include "models.h"
using namespace std;

int main(){
	Vector4d Vt;
	for (int i=0; i<10; i++){
		SlabUniform M1 = SlabUniform(30, 3.33e-16*(i+1), 10*AU, 100*AU);
		Vt<<M1.Integrate(0,0,10*AU, 0, 0);
		cout<<Vt[0]<<"\t"<<Vt[1]<<endl;
	}
}
