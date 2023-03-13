// EMDecayWidths includes
#ifndef EMECAYWIDTHS_H
#define EMDECAYWIDTHS_H

#include <string>
#include <vector>

class EMDecayWidths{

public:
  EMDecayWidths();
  virtual ~EMDecayWidths();
  virtual double execute(double ma_val, double mb_val, double mc_val, double ga_val, double sa_val,
			 double la_val, double ja_val, double sl_val, double al_val, double ar_val,
			 int baryon, int excMode, int prodDecay);
private:
  double MA; double MB; double MC;
  int modeExcitation=0;
  double pi_val = 3.1415926536;
  double JA = 0.5;   std::vector<double> mJA;
  double LA = 1.;    std::vector<double> mLA;
  double L1 = 0.;    std::vector<double> mL1;
  double L2 = 0.;    std::vector<double> mL2;  
  double SA = 0.5;   std::vector<double> mSA;
  double SB = 0.5;   std::vector<double> mSB;
  double SC = 0;     std::vector<double> mSC;
  double slight = 0; std::vector<double> m12; 
  std::vector<double> m23;
  double slightf = 1;std::vector<double> m24; 
  double s = 1.0;    std::vector<double> m;   
  double s1 = 0.5;   std::vector<double> m1;  
  double s2 = 0.5;   std::vector<double> m2;  
  double s3 = 0.5;   std::vector<double> m3;  
  double s4 = 0.5;   std::vector<double> m4;  
  double s5 = 0.5;   std::vector<double> m5;

  double LB = 0.0;  std::vector<double> mLB;
  double JB = 0.0;  std::vector<double> mJB;

  int baryonFlag=0;
  int decayProd=0;
  double E = 2.718281828;

  // check later if still needed
  int p_imag = 1;

  virtual double ALPHA_MES(int diagram);

  virtual int KroneckerDelta(double i, double j);
  virtual std::vector<double> getMomentumProjections(double j_angular);
  virtual double ANGULAR_SUM(double alpha_rho, double alpha_lam,
			     double alpha_mes, double k_value);
  virtual double DecayWidth(double decay, double gamma, double fi2_value, double angular_sum_value);

  virtual double EWCC(double MA, double MB, double MC);
  virtual double EB(double MA, double MB, double MC);
  virtual double K(double EB, double MB);
  virtual double FI2(double EB, double EWCC, double MA, double k_value);

  // SPIN-FLIP Integrals
  virtual double SPINFLIP_U1_GS_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double SPINFLIP_U2_GS_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double SPINFLIP_U3_GS_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double SPINFLIP_U1_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U1_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U1_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U2_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U2_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U2_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U3_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U3_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U3_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U1_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U1_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U1_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U2_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U2_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U2_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U3_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U3_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double SPINFLIP_U3_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);

  // ORBIT-SPLIT INTEGRALS
  virtual double ORBITALSPLIT_U1_1l_m1_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m1_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m1_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m0_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m0_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m1m_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m1m_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m1m_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);

  // New integrals
  virtual double ORBITALSPLIT_U3_1l_m0_1l_m0(double k_value, double alpha_lam, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U3_1l_m1m_1l_m1m(double k_value, double alpha_lam, double MB, double ML, double thethak);
  virtual double ORBITALSPLIT_U1_1l_m0_rl_m0(double k_value, double alpha_rho, double alpha_lam, double MB, double ML);
  virtual double ORBITALSPLIT_U2_1l_m0_rl_m0(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double ORBITALSPLIT_U1_1r_m1_rl_m1(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double ORBITALSPLIT_U1_1r_m0_rl_m0(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U1_1r_m1m_rl_m1m(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1r_m1_rl_m1(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1r_m0_rl_m0(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1r_m1m_rl_m1m(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U3_1r_m1_rl_m1(double k_value, double alpha_lam, double MB, double ML);
  virtual double ORBITALSPLIT_U3_1r_m0_rl_m0(double k_value, double alpha_lam, double MB, double ML);
  virtual double ORBITALSPLIT_U3_1r_m1m_rl_m1m(double k_value, double alpha_lam, double MB, double ML);
  virtual double ORBITALSPLIT_U1_1r_m0_ll_m0(double k_value, double alpha_lam, double alpha_rho,
						    double MB, double ML, double thetak, double phik);
  virtual double ORBITALSPLIT_U2_1r_m0_ll_m0(double k_value, double alpha_lam, double alpha_rho,
						    double MB, double ML, double thetak, double phik);

  //(13.03.2023)
  virtual double ORBITALSPLIT_U1_2l_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_2l_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U3_2l_m0_GS(double k_value, double alpha_lam, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U1_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U1_2nl_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1nl_m0_GS_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U3_1nl_m0_GS(double k_value, double alpha_lam, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U1_1nr_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1nr_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
};


//to talk to python
extern "C"{
  double electro_execute(double ma_val, double mb_val, double mc_val, double ga_val, double sa_val,
			double la_val, double ja_val, double sl_val, double al_val, double ar_val,
			int baryon, int excMode, int prodDecay){

    EMDecayWidths *  m_decays = new EMDecayWidths();
    return m_decays->execute(ma_val, mb_val, mc_val, ga_val, sa_val,
			     la_val, ja_val, sl_val, al_val, ar_val,
			     baryon, excMode, prodDecay);
  }
}

#endif //> !EMDECAYWIDTHS_H
