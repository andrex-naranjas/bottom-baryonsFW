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
			 double mb_val,
			 double al_val, double ar_val,
			 double mbottom_val, double mupdown_val, double mstrange_val,
			 int baryon, int excMode, int prodDecay);

 private:
  double MA; double MB; double MC; double mbottom; double mupdown; double mstrange; double mlight;
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

  double flav_q1 = 0.;
  double flav_q2 = 0.;
  double flav_q3 = 0.;

  // check later if still needed
  int p_imag = 1;

  virtual int KroneckerDelta(double i, double j);
  virtual int KroneckerDelta_extended(double i, double j, int excMode);

  virtual std::vector<double> getMomentumProjections(double j_angular, bool onlyPositive=false);
  virtual double ClebschGordan(WignerSymbols *m_wigner,
			       double l1, double l2, double l3,
			       double m1, double m2, double m3);
  virtual double ANGULAR_SUM_SQUARED(double alpha_rho, double alpha_lam, double k_value, int excMode);
  virtual double DecayWidth(double fi2_value, double angular_sum_value);
  virtual double EB(double MB, double K);
  virtual double K(double MA, double MB);
  virtual double FI2(double EB, double MA, double k_value);

  // SUMMARY INTEGRALS
  virtual double U1_rho_lambda(double k_value, double alpha_rho, double alpha_lam, double LA, int MLA, int MLB, int mLrA, int mLlA, int mLrB, int mLlB, int excMode);
  virtual double T1_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int LA, int mLrA, int mLlA, int mLrB, int mLlB, int excMode);
  virtual double U2_rho_lambda(double k_value, double alpha_rho, double alpha_lam, double LA, int MLA, int MLB, int mLrA, int mLlA, int mLrB, int mLlB, int excMode);
  virtual double T2_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int LA, int mLrA, int mLlA,  int mLrB, int mLlB, int excMode);
  virtual double U3_rho_lambda(double k_value, double alpha_rho, double alpha_lam, double LA, int MLA, int MLB, int mLrA, int mLlA, int mLrB, int mLlB, int excMode);
  virtual double T3_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int LA, int mLrA, int mLlA,  int mLrB, int mLlB, int excMode);

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


  // Definitions of the Integrals and Tensor Operators for decays from D-wave lambda to P-wave lambda
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_l2_m1_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U1_l2_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U1_l2_m1m_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_l2_m1_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_l2_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_l2_m1m_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_l2_m1_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_l2_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_l2_m1m_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double T1_l2_m2_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_l2_m1_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_l2_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho,double mbottom, double mlight, double thetak, double phik);
  virtual double T2_l2_m2_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_l2_m1_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_l2_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho,double mbottom, double mlight, double thetak, double phik);
  virtual double T3_l2_m2_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_l2_m1_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_l2_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);

  // Definitions of the Integrals and Tensor Operators for decays from D-wave rho to P-wave rho
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_r2_m1_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U1_r2_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U1_r2_m1m_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r2_m1_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r2_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r2_m1m_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_r2_m1_r1_m1();
  virtual double SPINFLIP_U3_r2_m0_r1_m0();
  virtual double SPINFLIP_U3_r2_m1m_r1_m1m();
  
  // ORBIT-FLIP Tensor operators
  virtual double T1_r2_m2_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r2_m1_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r2_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho,double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r2_m2_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r2_m1_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r2_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho,double mbottom, double mlight, double thetak, double phik);
  virtual double T3_r2_m2_r1_m1();
  virtual double T3_r2_m1_r1_m0();
  virtual double T3_r2_m0_r1_m1m();

  //Definitions of the Integrals and Tensor Operators for decays from D-wave lambda to P-wave rho
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_l2_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_l2_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_l2_m0_r1_m0();
  virtual double T1_l2_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_l2_m1_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_l2_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_l2_m1_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_l2_m0_r1_m1m();
  virtual double T3_l2_m1_r1_m0();

  //Definitions of the Integrals and Tensor Operators for decays from D-wave rho to P-wave lambda
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_r2_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r2_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_r2_m0_l1_m0();

  // ORBIT-FLIP Tensor operators
  virtual double T1_r2_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r2_m1_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r2_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r2_m1_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_r2_m0_l1_m1m();
  virtual double T3_r2_m1_l1_m0();

  //Definitions of the Integrals and Tensor Operators for decays from Mixed to P-wave rho
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_r1_m1_l1_m0_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U1_r1_m0_l1_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U1_r1_m1m_l1_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r1_m1_l1_m0_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r1_m0_l1_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r1_m1m_l1_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_r1_m1_l1_m0_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_r1_m0_l1_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_r1_m1m_l1_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);

  // ORBIT-FLIP Tensor operators
  virtual double T1_r1_m0_l1_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r1_m1_l1_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r1_m1_l1_m1_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r1_m0_l1_m1_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r1_m1_ml_1m_1r_1m_1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m0_l1_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m1_l1_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m1_l1_m1_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m0_l1_m1_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m1_ml_1m_1r_1m_1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_r1_m0_l1_m0_r1_m1m();
  virtual double T3_r1_m1_l1_m0_r1_m0();
  virtual double T3_r1_m1_l1_m1_r1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_r1_m0_l1_m1_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_r1_m1_ml_1m_1r_1m_1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);

  //Definitions of the Integrals and Tensor Operators for decays from Mixed to P-wave lambda
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_r1_m0_l1_m1_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U1_r1_m0_l1_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U1_r1_m0_l1_m1m_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r1_m0_l1_m1_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r1_m0_l1_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r1_m0_l1_m1m_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_r1_m0_l1_m1_l1_m1();
  virtual double SPINFLIP_U3_r1_m0_l1_m0_l1_m0();
  virtual double SPINFLIP_U3_r1_m0_l1_m1m_l1_m1m();

  // ORBIT-FLIP Tensor operators
  virtual double T1_r1_m0_l1_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r1_m0_l1_m1_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r1_m1_l1_m1_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r1_m1_l1_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r1_m1_l1_m1m_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m0_l1_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m0_l1_m1_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m1_l1_m1_l1_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m1_l1_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m1_l1_m1m_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_r1_m0_l1_m0_l1_m1m();
  virtual double T3_r1_m0_l1_m1_l1_m0();
  virtual double T3_r1_m1_l1_m1_l1_m1();
  virtual double T3_r1_m1_l1_m0_l1_m0();
  virtual double T3_r1_m1_l1_m1m_l1_m1m();

  //Definitions of the Integrals and Tensor Operators for decays from D-lambda to ground state
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_l2_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_l2_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_l2_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);

  // ORBIT-FLIP Tensor operators
  virtual double T1_l2_m1_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_l2_m1_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_l2_m1_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);

  //Definitions of the Integrals and Tensor Operators for decays from D-rho to ground state
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_r2_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r2_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_r2_m0_GS();

  // ORBIT-FLIP Tensor operators
  virtual double T1_r2_m1_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r2_m1_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_r2_m1_GS();

  //Definitions of the Integrals and Tensor Operators for decays from Mixed to ground state
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_r1_m0_l1_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_r1_m0_l1_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_r1_m0_l1_m0_GS();

  // ORBIT-FLIP Tensor operators
  virtual double T1_r1_m1_l1_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T1_r1_m0_l1_m1_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m1_l1_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_r1_m0_l1_m1_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_r1_m1_l1_m0_GS();
  virtual double T3_r1_m0_l1_m1_GS();

  //Definitions of the Integrals and Tensor Operators for decays from Radial lambda to ground state
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_nl1_l0_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_nl1_l0_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_nl1_l0_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);

  //Definitions of the Integrals and Tensor Operators for decays from Radial rho to ground state
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_nr1_r0_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_nr1_r0_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_nr1_r0_m0_GS();

  //Definitions of the Integrals and Tensor Operators for decays from Radial lambda to P wave lambda
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_nl1_l0_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_nl1_l0_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_nl1_l0_m0_l1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);

  // ORBIT-FLIP Tensor operators
  virtual double T1_nl1_l0_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_nl1_l0_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_nl1_l0_m0_l1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);

  //Definitions of the Integrals and Tensor Operators for decays from Radial lambda to P wave rho
  // SPIN-FLIP INTEGRALS
  virtual double SPINFLIP_U1_nl1_l0_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U2_nl1_l0_m0_r1_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight);
  virtual double SPINFLIP_U3_nl1_l0_m0_r1_m0();

  // ORBIT-FLIP Tensor operators
  virtual double T1_nl1_l0_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T2_nl1_l0_m0_r1_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak, double phik);
  virtual double T3_nl1_l0_m0_r1_m1m();




  // Sum of integrals
  virtual double U1Dl2Pl2m1m(double k_value, double alpha_lam, double alpha_rho, int MLlA, int MLlB);
  virtual double U2Dl2Pl2m1m(double k_value, double alpha_lam, double alpha_rho, int MLlA, int MLlB);
  virtual double U3Dl2Pl2m1m(double k_value, double alpha_lam, double alpha_rho, int MLlA, int MLlB);
  virtual double T1Dl2Pll2ml1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA, int MLlB);
  virtual double T2Dl2Pll2ml1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA, int MLlB);
  virtual double T3Dl2Pll2ml1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA, int MLlB);

  // Sum integrals 2 // Dwave lambda ->Pwave rho
  virtual double U1Dl2Pr2m1m(double k_value, double alpha_lam, double alpha_rho, int MLlA_val, int MLrB_val);     
  virtual double U2Dl2Pr2m1m(double k_value, double alpha_lam, double alpha_rho, int MLlA_val, int MLrB_val);
  virtual double U3Dl2Pr2m1m(double k_value, double alpha_lam, double alpha_rho, int MLlA_val, int MLrB_val);
  virtual double T1Dl2Prl2mr1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA_val, int MLrB_val);
  virtual double T2Dl2Prl2mr1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA_val, int MLrB_val);
  virtual double T3Dl2Prl2mr1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA_val, int MLrB_val);
  virtual double T1Dr2Plr2ml1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlB_val);
  virtual double T1Dr2Prr2mr1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLrB_val);
  virtual double T2Dr2Plr2ml1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlB_val);
  virtual double T2Dr2Prr2mr1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLrB_val);
  virtual double T3Dr2Plr2ml1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlB_val);
  virtual double T3Dr2Prr2mr1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLrB_val);
  virtual double U1Dr2Pl2m1m(double k_value, double alpha_lam, double alpha_rho, int MLrA_val, int MLlB_val);
  virtual double U1Dr2Pr2m1m(double k_value, double alpha_lam, double alpha_rho, int MLrA_val, int MLrB_val);
  virtual double U2Dr2Pl2m1m(double k_value, double alpha_lam, double alpha_rho, int MLrA_val, int MLlB_val);
  virtual double U2Dr2Pr2m1m(double k_value, double alpha_lam, double alpha_rho, int MLrA_val, int MLrB_val);
  virtual double U3Dr2Pl2m1m(double k_value, double alpha_lam, double alpha_rho, int MLrA_val, int MLlB_val);
  virtual double U3Dr2Pr2m1m(double k_value, double alpha_lam, double alpha_rho, int MLrA_val, int MLrB_val);

  virtual double U1Mix2Pl1m(double k_value, double alpha_lam, double alpha_rho, int MLlA_val, int MLrA_val, int MLlB_val);
  virtual double U2Mix2Pl1m(double k_value, double alpha_lam, double alpha_rho, int MLlA_val, int MLrA_val, int MLlB_val);
  virtual double U3Mix2Pl1m(double k_value, double alpha_lam, double alpha_rho, int MLlA_val, int MLrA_val, int MLlB_val);
  virtual double T1Mix2Pl1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA_val, int MLrA_val, int MLlB_val);
  virtual double T2Mix2Pl1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA_val, int MLrA_val, int MLlB_val);
  virtual double T3Mix2Pl1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA_val, int MLrA_val, int MLlB_val);
  virtual double U1Mix2Pr1m(double k_value, double alpha_lam, double alpha_rho, int MLlA_val, int MLrA_val, int MLrB_val);
  virtual double U2Mix2Pr1m(double k_value, double alpha_lam, double alpha_rho, int MLlA_val, int MLrA_val, int MLrB_val);
  virtual double U3Mix2Pr1m(double k_value, double alpha_lam, double alpha_rho, int MLlA_val, int MLrA_val, int MLrB_val);
  virtual double T1Mix2Pr1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA_val, int MLrA_val, int MLrB_val);
  virtual double T2Mix2Pr1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA_val, int MLrA_val, int MLrB_val);
  virtual double T3Mix2Pr1m(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLlA_val, int MLrA_val, int MLrB_val);

  // Sum of integrals Dwave lambda -> ground state
  virtual double U1Dl2GS2mgs(double k_value, double alpha_lam, double alpha_rho, int MLlA_val);
  virtual double U2Dl2GS2mgs(double k_value, double alpha_lam, double alpha_rho, int MLlA_val);
  virtual double U3Dl2GS2mgs(double k_value, double alpha_lam, double alpha_rho, int MLlA_val);
  virtual double T1Dl2GS2mgs(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlA_val);
  virtual double T2Dl2GS2mgs(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlA_val);
  virtual double T3Dl2GS2mgs(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlA_val);

  // Sum of integrals Dwave rho -> ground state
  virtual double U1Dr2GS2mgs(double k_value, double alpha_lam, double alpha_rho, int MLlA_val);
  virtual double U2Dr2GS2mgs(double k_value, double alpha_lam, double alpha_rho, int MLlA_val);
  virtual double U3Dr2GS2mgs(double k_value, double alpha_lam, double alpha_rho, int MLlA_val);
  virtual double T1Dr2GS2mgs(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlA_val);
  virtual double T2Dr2GS2mgs(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlA_val);
  virtual double T3Dr2GS2mgs(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlA_val);

  // Sum of integrals Mixed -> ground state
  virtual double U1Mix2GS1mgs(double k_value, double alpha_lam, double alpha_rho, int MLrA_val, int MLlA_val);
  virtual double U2Mix2GS1mgs(double k_value, double alpha_lam, double alpha_rho, int MLrA_val, int MLlA_val);
  virtual double U3Mix2GS1mgs(double k_value, double alpha_lam, double alpha_rho, int MLrA_val, int MLlA_val);
  virtual double T1Mix2GS1mgs(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlA_val);
  virtual double T2Mix2GS1mgs(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlA_val);
  virtual double T3Mix2GS1mgs(double k_value, double alpha_lam, double alpha_rho, double thetak, double phik,  int MLrA_val, int MLlA_val);
  
};



//to talk to python
extern "C"{
  double electro_execute(double ma_val, double sa_val, double ja_val, double la_val, double sla_val, double lla_val, double lra_val,
			 double mb_val,
			 double al_val, double ar_val,
			 double mbottom_val,  double mupdown_val, double mstrange_val,
			 int baryon, int excMode, int prodDecay){
    EMDecayWidths *m_decays = new EMDecayWidths();
    return m_decays->execute(ma_val, sa_val, ja_val, la_val, sla_val, lla_val, lra_val,
			     mb_val,
			     al_val, ar_val,
			     mbottom_val, mupdown_val, mstrange_val,
			     baryon, excMode, prodDecay);
  }
}

#endif //> !EMDECAYWIDTHS_H
