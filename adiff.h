#ifndef _ADIFF_H
#define _ADIFF_H

#include "models.h"
#include "utils.h"

class ADiff : public HLTau{
  protected:
    double Kabs;
  public:
    ADiff();
    void cal_VM(double x, double y, double z, double n_theta, double n_phi, Vector4d &Vout, Matrix4d &Mout);
    Matrix4d get_ZMatrix(double theta_i, double phi_i, double theta_o, double phi_o);
    Vector4d Image(double x, double y, double z, double l_theta, double l_phi, double step=0.1*AU);
    double get_Kext(double x, double y, double z);
    double get_Ksca(double x, double y, double z);
};

#endif
