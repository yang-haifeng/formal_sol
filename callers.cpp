#include "callers.h"

int call_HLTau(){
  ifstream Fpar;
  Fpar.open("param.inp");
  if(!Fpar.is_open()){
    cout<<"ERROR: failed to open parameter file: param.inp"<<endl;
    return 1;
  }
  string line;

  getline(Fpar, line); 
  if( strncmp(line.substr(0,5).c_str(), "hltau", 5) ){
    cout<<"ERROR: Problems other than HLTau are passed to call_HLTau."<<endl;
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

  // Seventh line is the calculation mode
  int NMode;
  Fpar>>NMode;
  getline(Fpar, line);

  HLTau M = HLTau();
  M.set_kappa(Cabs+Csca, Csca);
  if (rho0!=-1) M.set_rho0(rho0);
  else cout<<"Using default density profile."<<endl;

  switch(NMode){
  case 0:
    cout<<"Doing HLTau calculation in Standard Mode with a grid of 15x16 points."<<endl;
    M.get_Circle_Image(theta, 15, 16, 10*AU, 150*AU, Fout);
    break;
  case 1:
    cout<<"Doing calculation at given radius."<<endl;
    // Eighth line is Number of points in Radial direction
    int nr, np;
    Fpar>>nr; getline(Fpar, line);
    Fpar>>np; getline(Fpar, line);
    double R;
    for (int i=0;i<nr;i++){
      Fpar>>R;
      cout<<"Calculating at radius: "<<R<<endl;
    }
  case 2:
    cout<<"Restart Mode not implemented yet. I'll do nothing here for now."<<endl;
    break;
  default:
    cout<<"ERROR: Calculation Mode "<<NMode<<" is not valid. please use 0 for standard mode, 1 for given radius and 2 for restart."<<endl;
    return 3;
  }

  return 0;
}
