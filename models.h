#include <iostream>
#include <math.h>
#include "typedef.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Model{
	protected:
		double r_max; 
		double Kext, Kpol, Kcpol; 
		double lambda; 
		double los_theta, los_phi;
		//Matrix4d Mext;
		//Vector4d Vabs;
	public:
		Model();
		virtual double get_BnuT(double x, double y, double z);
		virtual double get_Rho(double x, double y, double z);
		virtual void get_Orientation(double x, double y, double z, double &theta, double &phi);
		virtual bool reachBoundary(double x, double y, double z);

		Vector4d Integrate(double x, double y, double z, double n_theta, double n_phi, double step=0.1*AU);
		virtual void cal_VM(double x, double y, double z, double n_theta, double n_phi,
				Vector4d &Vout, Matrix4d &Mout);

		virtual void test();
};

class SlabUniform : public Model{
	protected:
		double Bnu0, rho0;
		double h; 
	public:
		SlabUniform();
		SlabUniform(double Temp, double rho, double height, double rm);
		double get_BnuT(double x, double y, double z);
		double get_Rho(double x, double y, double z);
		virtual void get_Orientation(double x, double y, double z, double &theta, double &phi);
		bool reachBoundary(double x, double y, double z);
};

class SlabTwisted : public SlabUniform{
	protected:
		double k;
	public:
		SlabTwisted();
		SlabTwisted(double Temp, double rho, double height, double rm, double K);
		virtual void get_Orientation(double x, double y, double z, double &theta, double &phi);
		void test();
};

double BnuT(double T, double nu);
