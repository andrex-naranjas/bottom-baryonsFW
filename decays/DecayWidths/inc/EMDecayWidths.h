// EMDecayWidths includes
#ifndef EMECAYWIDTHS_H
#define EMDECAYWIDTHS_H

#include "WignerSymbols.h"

#include <string>
#include <vector>

class EMDecayWidths{

 public:
  EMDecayWidths();
  virtual ~EMDecayWidths();
  virtual double execute(double ma_val, double sa_val, double ja_val, double la_val, double sla_val, double lla_val, double lra_val,
			 double mb_val, double sb_val, double jb_val, double lb_val, double slb_val, double llb_val, double lrb_val,
			 double al_val, double ar_val, double mstrange_val,
			 double mbottom_val, double mlight_val, int baryon, int excMode, int prodDecay);

 private:
  double MA; double MB; double MC; double mstrange; double mbottom; double mlight;
  int modeExcitation=0;
  double pi_val = 3.1415926536;

  double L   = 0.0;  std::vector<double> mL;

  double SA  = 0.0;  std::vector<double> mSA;
  double JA  = 0.0;  std::vector<double> mJA;
  double SlA = 0.0;  std::vector<double> mSlA;
  double LA  = 0.0;  std::vector<double> mLA;
  double LlA = 0.0;  std::vector<double> mLlA;
  double LrA = 0.0;  std::vector<double> mLrA;

  double SB  = 0.0;  std::vector<double> mSB;
  double JB  = 0.0;  std::vector<double> mJB;
  double SlB = 0.0;  std::vector<double> mSlB;
  double LB  = 0.0;  std::vector<double> mLB;
  double LlB = 0.0;  std::vector<double> mLlB;
  double LrB = 0.0;  std::vector<double> mLrB;

  double S1  = 0.5;  std::vector<double> mS1;
  double S2  = 0.5;  std::vector<double> mS2;
  double S3  = 0.5;  std::vector<double> mS3;

  int baryonFlag=0;
  int decayProd=0;

  // check later if still needed
  int p_imag = 1;

  virtual int KroneckerDelta(double i, double j);
  virtual int KroneckerDelta_extended(double i, double j, int excMode);
  virtual std::vector<double> getMomentumProjections(double j_angular);
  virtual double ClebschGordan(WignerSymbols *m_wigner,
			       double l1, double l2, double l3,
			       double m1, double m2, double m3);
  virtual std::vector<double> FlavorVector(double charge=0);
  virtual double ANGULAR_SUM_SQUARED(double alpha_rho, double alpha_lam, double k_value, std::vector<double>, int excMode);
  virtual double ANGULAR_SUM_SQUARED_GROUND(double alpha_rho, double alpha_lam, double k_value);
  virtual double DecayWidth(double decay, double fi2_value, double angular_sum_value);
  virtual double EB(double MB, double K);
  virtual double K(double MA, double MB);
  virtual double FI2(double EB, double MA, double k_value);

  // SUMMARY INTEGRALS
  virtual double U1_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode);
  virtual double T1_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode);
  virtual double U2_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode);
  virtual double T2_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode);
  virtual double U3_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode);
  virtual double T3_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode);

  // SPIN-FLIP Integrals
  virtual double SPINFLIP_U1_GS_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double SPINFLIP_U2_GS_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double SPINFLIP_U3_GS_GS(double k_value, double alpha_lam, double MB, double ML);
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
  virtual double ORBITALSPLIT_U1_1l_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m0_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m1m_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m1m_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U1_1l_m1m_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double thetak);

  // New integrals
  virtual double ORBITALSPLIT_U3_1l_m0_1l_m0(double k_value, double alpha_lam, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U3_1l_m1m_1l_m1m(double k_value, double alpha_lam, double MB, double ML, double thethak);
  virtual double ORBITALSPLIT_U1_1l_m0_1r_m0(double k_value, double alpha_rho, double alpha_lam, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1l_m0_1r_m0(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U1_1r_m1_1r_m1(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U1_1r_m0_1r_m0(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U1_1r_m1m_1r_m1m(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1r_m1_1r_m1(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1r_m0_1r_m0(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1r_m1m_1r_m1m(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U3_1r_m1_1r_m1(double k_value, double alpha_lam, double MB, double ML);
  virtual double ORBITALSPLIT_U3_1r_m0_1r_m0(double k_value, double alpha_lam, double MB, double ML);
  virtual double ORBITALSPLIT_U3_1r_m1m_1r_m1m(double k_value, double alpha_lam, double MB, double ML);
  virtual double ORBITALSPLIT_U1_1r_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,
					     double MB, double ML, double thetak, double phik);
  virtual double ORBITALSPLIT_U2_1r_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,
					     double MB, double ML, double thetak, double phik);


  virtual double ORBITALSPLIT_U2_1l_m1_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1l_m1_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U2_1l_m1_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U2_1l_m0_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U2_1l_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_1l_m0_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U2_1l_m1m_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U2_1l_m1m_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak);
  virtual double ORBITALSPLIT_U2_1l_m1m_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U3_1l_m1_1l_m1(double k_value, double alpha_lam, double MB, double ML, double thetak);

  //(13.03.2023)
  virtual double ORBITALSPLIT_U1_2l_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_2l_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U3_2l_m0_GS(double k_value, double alpha_lam, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U1_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  virtual double ORBITALSPLIT_U2_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak);
  // virtual double ORBITALSPLIT_U1_2nl_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  // virtual double ORBITALSPLIT_U2_1nl_m0_GS_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double ORBITALSPLIT_U3_1nl_m0_GS(double k_value, double alpha_lam, double MB, double ML);
  virtual double ORBITALSPLIT_U1_1nr_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double ORBITALSPLIT_U2_1nr_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double ORBITALSPLIT_U1_1nl_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);
  virtual double ORBITALSPLIT_U2_1nl_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML);

  // Tensors
  //T1l
  virtual double T1l(double k_value, double alpha_lam, double alpha_rho,
		     double MB, double ML, double thetak, double phik, double mLlA);
  virtual double T2l(double k_value, double alpha_lam, double alpha_rho,
			    double MB, double ML, double thetak, double phik, double mLlA);
  virtual double T3l(double k_value, double alpha_lam, double alpha_rho,
			    double MB, double ML, double thetak, double phik, double mLlA);
  virtual double T1r(double k_value, double alpha_lam, double alpha_rho,
			    double MB, double ML, double thetak, double phik, double mLrA);
  virtual double T2r(double k_value, double alpha_lam, double alpha_rho,
			    double MB, double ML, double thetak, double phik, double mLrA);
  virtual double T3r();
};


//to talk to python
extern "C"{
  double electro_execute(double ma_val, double sa_val, double ja_val, double la_val, double sla_val, double lla_val, double lra_val,
			 double mb_val, double sb_val, double jb_val, double lb_val, double slb_val, double llb_val, double lrb_val,
			 double al_val, double ar_val,
			 double mbottom_val, double mlight_val, int baryon, int excMode, int prodDecay){
    EMDecayWidths *m_decays = new EMDecayWidths();
    return m_decays->execute(ma_val, sa_val, ja_val, la_val, sla_val, lla_val, lra_val,
			     mb_val, sb_val, jb_val, lb_val, slb_val, llb_val, lrb_val,
			     al_val, ar_val,
			     mbottom_val, mlight_val, baryon, excMode, prodDecay);
  }
}

#endif //> !EMDECAYWIDTHS_H
