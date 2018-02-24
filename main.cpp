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

  // First line is problem type. 
  int Nresult;
  getline(Fpar, line); 
  if( strncmp(line.substr(0,5).c_str(), "hltau", 5) == 0 ){
    Nresult = call_HLTau();
  }
  else if( strncmp(line.substr(0,6).c_str(), "warped", 6) == 0 ){
    Nresult = call_warped();
  }
  else{
    cout<<"ERROR: problem type not understood. The first line reads:"<<endl;
    cout<<line<<endl;
  }
  Fpar.close();

  if( Nresult!=0 ){
    cout<<"Error raised during the calculation. Error code: "<<Nresult<<endl;
  }

  return Nresult;
}
