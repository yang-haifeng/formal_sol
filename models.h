#ifndef _MODEL_H
#define _MODEL_H

#include <iostream>
#include <fstream>
#include <math.h>
#include "typedef.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

class Model{
	protected:
		double r_max; 
		double Kext, Kpol, Kcpol; 
		double Ksca;
		double lambda; 
		double tau_ad;
		//double los_theta, los_phi;
		//Matrix4d Mext;
		//Vector4d Vabs;
	public:
		Model();
		virtual double get_BnuT(double x, double y, double z);
		virtual double get_Rho(double x, double y, double z);
		virtual void get_Orientation(double x, double y, double z, double &theta, double &phi);
		virtual bool reachBoundary(double x, double y, double z);
		bool reachBoundary(Vector3d P);
		virtual Matrix4d get_ZMatrix(double theta_i, double phi_i, double theta_o, double phi_o);

		Vector4d Integrate(double x, double y, double z, double n_theta, double n_phi, double step=0.1*AU);
		Vector4d Image(double x, double y, double z, double l_theta, double l_phi, double step=0.1*AU);
		void get_Image(double theta, int Npx, double FoV, string fName);
		void get_Circle_Image(double theta, int Nr, int Nph, double Rin, double Rout, string fName);
		virtual void cal_VM(double x, double y, double z, double n_theta, double n_phi,
				Vector4d &Vout, Matrix4d &Mout);

		virtual void test();
		void set_adaptive(double tau);
};

class SlabUniform : public Model{
	protected:
		double Bnu0, rho0;
		double h; 
	public:
		using Model::reachBoundary;
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

class SlabSphGrain : public SlabUniform{
	public:
		SlabSphGrain();
		SlabSphGrain(double Temp, double rho, double height, double rm);
		void get_Orientation(double x, double y, double z, double &theta, double &phi);
};

class HLTau : public Model{
	protected:
		double Rc, T0, H0, p, q, rho0;
	public:
		using Model::reachBoundary;
		HLTau();
		double get_BnuT(double x, double y, double z);
		double get_Rho(double x, double y, double z);
		bool reachBoundary(double x, double y, double z);
		void set_rho0(double rho);
		void multiply_rho0(double ratio=2);
		double get_HR(double R);
};

class ConeModel : public Model{
	protected:
		double theta_cone;
		double rho0;
	public:
		using Model::reachBoundary;
		ConeModel(double theta, double rho);
		double get_BnuT(double x, double y, double z);
		double get_Rho(double x, double y, double z);
		bool reachBoundary(double x, double y, double z);
};

double BnuT(double T, double nu);

#endif
