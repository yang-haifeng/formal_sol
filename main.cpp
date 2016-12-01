#include "callers.h"
using namespace std;

int main(){
  ifstream Fpar;
  Fpar.open("param.inp");
  if(!Fpar.is_open()){
    cout<<"ERROR: failed to open parameter file: param.inp"<<endl;
    return 1;
  }
  string line;

  // First line is problem type. Only accept hltau at this point.
  getline(Fpar, line); 
  if( strncmp(line.substr(0,5).c_str(), "hltau", 5) ){
    cout<<"ERROR: problems other than HLTau are not implemented yet."<<endl;
    return 2;
  }
  Fpar.close();

  int Nresult = call_HLTau();
  if( Nresult!=0 ){
    cout<<"Error raised in call_HLTau(). Error code: "<<Nresult<<endl;
  }

  return 0;
}
