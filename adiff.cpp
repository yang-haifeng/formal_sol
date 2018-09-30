#include "adiff.h"

ADiff::ADiff(){
        r_max = 200*AU; Rc = 79*AU;
        T0 = 70; H0 = 16.8 * AU;
        //rho0 = 1.964e-15;
        //rho0 = 4.7166961619e-15; // New number from problem_setup.py. Need to find why there was a difference.
        rho0 = 1.12750603885e-15; // WTF man, WTF!
        p = 1.064; q = 0.43;
        Kext = 1.29; Kpol=Kcpol=0;
        Ksca = 0.78;
        lambda = 0.1;
        tau_ad = 0.1;

	Kabs = Kext - Ksca;
}

void ADiff::cal_VM(double x, double y, double z, double n_theta, double n_phi,
                Vector4d &Vout, Matrix4d &Mout){
        double cosDTheta;
        double theta, phi;
        get_Orientation(x, y, z, theta, phi);
        double dphi = phi-n_phi;
        //cosDTheta = sin(theta)*cos(phi)*sin(n_theta)*cos(n_phi) +
                //sin(theta)*sin(phi)*sin(n_theta)*sin(n_phi) + cos(theta)*cos(n_theta);
        cosDTheta = sin(theta)*cos(dphi)*sin(n_theta) + cos(theta)*cos(n_theta);

	double kext = get_Kext(x, y, z);

        double Cep = kext + Kpol * cosDTheta*cosDTheta;
        double Cpp = Kpol * (1- cosDTheta*cosDTheta);
        double Ccpp = Kcpol  * (1- cosDTheta*cosDTheta);

        double xp, yp; // The xp,yp after rotation, not the global coordinates.
        xp = sin(theta)*cos(dphi)*cos(n_theta) - cos(theta)*sin(n_theta);
        yp = sin(theta)*sin(dphi);
        double cos2phi = (xp*xp-yp*yp)/(xp*xp+yp*yp);
        double sin2phi = 2*xp*yp/(xp*xp+yp*yp);
        if ((xp*xp+yp*yp)==0) cos2phi=sin2phi=0;

        //Vout << Cep, -Cpp*cos2phi, -Cpp*sin2phi, 0; // This needs some fixing. Cep v.s. Cap
        Vout << Kabs, 0., 0., 0; // Temporary fixing. Only deals with spherical grains.
        Mout << Cep, -Cpp*cos2phi, -Cpp*sin2phi, 0,
             -Cpp*cos2phi, Cep, 0, Ccpp*sin2phi,
             -Cpp*sin2phi, 0, Cep, -Ccpp*cos2phi,
             0, -Ccpp*sin2phi, Ccpp*cos2phi, Cep;
}

Matrix4d ADiff::get_ZMatrix(double theta_i, double phi_i, double theta_o, double phi_o){
        //double dphi = phi_o-phi_i;
        Vector3d et1, ep1, et2, ep2;
        et1 << cos(theta_i) * cos(phi_i), cos(theta_i) * sin(phi_i), -sin(theta_i);
        ep1 << -sin(phi_i), cos(phi_i), 0;
        et2 << cos(theta_o) * cos(phi_o), cos(theta_o) * sin(phi_o), -sin(theta_o);
        ep2 << -sin(phi_o), cos(phi_o), 0;

        double factor = 3./8./PI;
        double ftt, ftp, fpt, fpp;
        ftt = et2.transpose() * et1;
        ftp = et2.transpose() * ep1;
        fpt = ep2.transpose() * et1;
        fpp = ep2.transpose() * ep1;

        double Z11, Z12, Z13, Z14;
        double Z21, Z22, Z23, Z24;
        double Z31, Z32, Z33, Z34;
        double Z41, Z42, Z43, Z44;

        Z11 = 0.5 * (ftt*ftt + ftp*ftp + fpt*fpt + fpp*fpp);
        Z12 = 0.5 * (ftt*ftt - ftp*ftp + fpt*fpt - fpp*fpp);
        Z13 = ftt*ftp + fpp*fpt;
        Z14 = 0.;
        Z21 = 0.5 * (ftt*ftt + ftp*ftp - fpt*fpt - fpp*fpp);
        Z22 = 0.5 * (ftt*ftt - ftp*ftp - fpt*fpt + fpp*fpp);
        Z23 = ftt*ftp - fpp*fpt;
        Z24 = 0.;
        Z31 = ftt*fpt + fpp*ftp;
        Z32 = ftt*fpt - fpp*ftp;
        Z33 = fpp*ftt + ftp*fpt;
        Z34 = 0.;
        Z41 = 0.;
        Z42 = 0.;
        Z43 = 0.;
        Z44 = fpp*ftt - fpt*ftp;

        Matrix4d M;
        M <<    Z11, Z12, Z13, Z14,
                Z21, Z22, Z23, Z24,
                Z31, Z32, Z33, Z34,
                Z41, Z42, Z43, Z44;
        M*=factor;
        return M;
}

Vector4d ADiff::Image(double x, double y, double z, double l_theta, double l_phi, double step0){
        Matrix4d T = Matrix4d(Vector4d::Constant(1).asDiagonal());
        Vector4d result = Vector4d::Constant(0);
        double rho, bnuT;
        double xp=x, yp=y, zp=z;
        double step = step0;
        double dx = -step*sin(l_theta)*cos(l_phi);
        double dy = -step*sin(l_theta)*sin(l_phi);
        double dz = -step*cos(l_theta);
        Vector4d Vabs; Matrix4d Mext;
        Matrix4d Z;
        double theta, phi;
        Vector4d Sin, Ssca;
        Matrix4d R1, R2;
        double psi1, psi2;
        double btheta, bphi;
        double theta1, phi1;
        double theta2, phi2;
        while (true){
                //cout<<"***********************************************"<<endl;
                //cout<<xp/AU<<"\t"<<yp/AU<<"\t"<<zp/AU<<"\t"<<step/AU<<endl;
                rho = get_Rho(xp, yp, zp);
                bnuT = get_BnuT(xp, yp, zp);

                if (tau_ad>0){
                        if (rho != 0){
                        step = tau_ad / rho / Kext;
                        if (step>AU) step=AU;
                        dx = -step*sin(l_theta)*cos(l_phi);
                        dy = -step*sin(l_theta)*sin(l_phi);
                        dz = -step*cos(l_theta);
                        }
                        else{
                        step = AU;
                        dx = -step*sin(l_theta)*cos(l_phi);
                        dy = -step*sin(l_theta)*sin(l_phi);
                        dz = -step*cos(l_theta);
                        }
                }

                cal_VM(xp, yp, zp, l_theta, l_phi, Vabs, Mext);

                Ssca = Vector4d::Constant(0);
                /*
                */
                get_Orientation(xp, yp, zp, btheta, bphi);
                to_grain_frame(btheta, bphi, l_theta, l_phi, theta2, phi2);
                psi2 = get_psi(btheta, bphi, theta, phi);
                for (int i=0; i<Ntheta; i++){
                        theta = PI/Ntheta*i;
                        for (int j=0; j<Nphi; j++){
                                phi = 2*PI/Nphi*j;
                                Sin = Integrate(xp, yp, zp, theta, phi);

                                //cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"<<endl;
                                //cout<<theta<<"\t"<<phi<<endl;
                                //cout<<Sin<<endl;

                                psi1 = get_psi(btheta, bphi, theta, phi);
                                to_grain_frame(btheta, bphi, theta, phi, theta1, phi1);

                                //cout<<psi2<<" rotation Matrix:"<<endl;
                                //cout<<rotation_Matrix(-psi2)<<endl;
                                //cout<<"Z Matrix:"<<endl;
                                //cout<<get_ZMatrix(theta1, phi1, theta2, phi2)<<endl;
                                //cout<<psi1<<" rotation Matrix:"<<endl;
                                //cout<<rotation_Matrix(psi1)<<endl;

                                Ssca += rotation_Matrix(-psi2) *
                                        get_ZMatrix(theta1, phi1, theta2, phi2) *
					get_Ksca(xp, yp, zp) *
                                        rotation_Matrix(psi1) *
                                        Sin *
                                        rho * sin(theta) * PI/Ntheta * 2*PI/Nphi;
                                //cout<<Ssca<<endl;
                        }
                }

                result += T*(Vabs * rho * bnuT + Ssca)*step;
                cout<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<endl;

                T -= T*Mext * step * rho;
                xp += dx; yp += dy; zp += dz;
                if (reachBoundary(xp, yp, zp)) {
                //cout<<xp/AU<<"\t"<<yp/AU<<"\t"<<zp/AU<<"\tBoundary reached"<<endl;
                break;}
        }
        return result;
}

double ADiff::get_Kext(double x, double y, double z){
  return Kabs + get_Ksca(x,y,z);
}

double ADiff::get_Ksca(double x, double y, double z){
  // Let's do a very simple model with grain size decrease by a factor of 10 
  // linearly in radius as we go from center to r_max
  double R = sqrt(x*x+y*y);
  return Ksca * (1 - 0.9 * (R/r_max));
}



