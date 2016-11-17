#include <iostream>
#include <fstream>
#include <string>
#include "typedef.h"
#include "models.h"
#include "utils.h"
using namespace std;

int main(){
  ifstream Fpar;
  Fpar.open("param.inp");
  if(!Fpar.is_open()){
    cout<<"Failed to open parameter file: param.inp"<<endl;
    return 1;
  }
  string line;

  // First line is problem type. Only accept hltau at this point.
  getline(Fpar, line); 
  if( strncmp(line.substr(0,5).c_str(), "hltau", 5) ){
    cout<<"Error: problems other than HLTau are not implemented yet."<<endl;
    return 2;
  }

  // Second line is rho0. Will take -1 to use default value.
  double rho0;
  Fpar>>rho0;
  getline(Fpar, line);

  // Third and forth lines are extinction and scattering opacity, respectively.
  double Cabs, Csca;
  Fpar>>Cabs;
  getline(Fpar, line);
  Fpar>>Csca;
  getline(Fpar, line);

  // Fifth line is viewing angle
  double theta;
  Fpar>>theta;
  getline(Fpar, line);
  theta *= PI/180.;

  // Sixth line is output file name.
  getline(Fpar, line);
  size_t spliter = line.find_first_of(" ;\t");
  if (spliter!=string::npos){
    line = line.substr(0,spliter);
  }
  string Fout=line;
  cout<<"Output file: "<<Fout<<endl;

  HLTau M = HLTau();
  M.set_kappa(Cabs+Csca, Csca);
  if (rho0!=-1) M.set_rho0(rho0);
  else cout<<"Using default density profile."<<endl;
  M.get_Circle_Image(theta, 15, 16, 10*AU, 150*AU, Fout);

  return 0;
}
