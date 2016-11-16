#include <iostream>
#include <fstream>
#include "typedef.h"
#include "models.h"
#include "utils.h"
using namespace std;

int main(){
  ifstream Fpar;
  Fpar.open("param.inp")
  if(!Fpar.is_open()){
    cout<<"Failed to open parameter file: param.inp"<<endl;
    return 1;
  }
}
