#include <iostream>
#include <fstream>
#include "typedef.h"
#include "models.h"
using namespace std;

int main(){
	Model M1 = Model();
	M1.get_BnuT(10, 10, 10);
	M1.test();
	return 0;
}
