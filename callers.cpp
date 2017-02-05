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

  // Add one line for Hfactor. 
  // Note that previous parameter file will fail because of this. 
  // Need to add one line with only number 1 in it.
  double Hfactor;
  Fpar>>Hfactor;
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
  M.multiply_H0(Hfactor);

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
      M.get_Circle_Image(theta, 1, np, R*AU, R*AU, Fout, true);
    }
    break;
  case 2:
    cout<<"Restart Mode not implemented yet. I'll do nothing here for now."<<endl;
    break;
  case 3:
    cout<<"Doing calculation slice by slice."<<endl;
    // Eighth line is Number of points in Radial direction
    Fpar>>nr; getline(Fpar, line);
    Fpar>>np; getline(Fpar, line);

    // 10th and 11th lines are number of threads and id of this thread
    int Nthreads, IDthread;
    Fpar>>Nthreads; getline(Fpar, line);
    Fpar>>IDthread; getline(Fpar, line);
    cout<<Nthreads<<" threads in total. This is thread #"<<IDthread<<"."<<endl;

    M.get_Slice_Image(theta, nr, np, 10*AU, 150*AU, Fout, Nthreads, IDthread);
    break;
  default:
    cout<<"ERROR: Calculation Mode "<<NMode<<" is not valid. please use 0 for standard mode, 1 for given radius and 2 for restart. 3 is for doing slice calculation."<<endl;
    return 3;
  }

  return 0;
}
