// EMDecay Class includes
#ifndef EMDECAYWIDTHS_CXX
#define EMDECAYWIDTHS_CXX

#include "EMDecayWidths.h"
#include "WignerSymbols.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
//#include <complex>


EMDecayWidths::EMDecayWidths()
{
}

EMDecayWidths::~EMDecayWidths(){}

double EMDecayWidths::execute(double ma_val, double mb_val, double sa_val,
			      double la_val, double ja_val, double sl_val, double al_val, double ar_val,
			      double mbottom_val, double mlight_val, int baryon, int excMode, int prodDecay){  
  // decay product masses
  MA = ma_val;
  MB = mb_val;  
  if(MA<MB) return 0.; //energy conservation

  // quark masses
  mbottom = mbottom_val;
  mlight = mlight_val;

  // which baryon, mode, decay product
  baryonFlag = baryon;
  modeExcitation = excMode;
  decayProd  = prodDecay;
  double alpha_rho = 0.,alpha_lam = 0.,flav_coup= 0.;
  double slf_val=0., sb_val=0., lb_val=0., jb_val=0.;
  double l1_val=0., l2_val=0.;

  alpha_rho = ar_val;
  alpha_lam = al_val;

  //fetch quantum numbers and projections
  SA = sa_val;      mSA = getMomentumProjections(SA);
  LA = la_val;      mLA = getMomentumProjections(LA);
  JA = ja_val;      mJA = getMomentumProjections(JA);
  SB = sb_val;      mSB = getMomentumProjections(SB);
  LB = lb_val;      mLB = getMomentumProjections(LB);
  JB = jb_val;      mJB = getMomentumProjections(JB);
  
  slight = sl_val;  m12 = getMomentumProjections(slight);
  slight = sl_val;  m23 = getMomentumProjections(slight);
  slightf= slf_val; m24 = getMomentumProjections(slightf);
  L1 = l1_val;      mL1 = getMomentumProjections(L1);
  L2 = l2_val;      mL2 = getMomentumProjections(L2);

  mSC = getMomentumProjections(SC);  
  //values are the same for all states (at least for now!)
  s  = 1.0;   m  = getMomentumProjections(s);
  s1 = 0.5;   m1 = getMomentumProjections(s1);
  s2 = 0.5;   m2 = getMomentumProjections(s2);
  s3 = 0.5;   m3 = getMomentumProjections(s3);
  s4 = 0.5;   m4 = getMomentumProjections(s4);
  s5 = 0.5;   m5 = getMomentumProjections(s5);

  double k_value; k_value = K(MA, MB);
  double EB_value = EB(MB, k_value);
  
  double sum_value  = 0.; //test
  sum_value = ANGULAR_SUM(alpha_rho, alpha_lam, k_value);

  double fi2_value  = FI2(EB_value, MA, k_value);
  double decayWidth = DecayWidth(flav_coup, fi2_value, sum_value);
  decayWidth = 1;

  double test_integral = 1; //SPINFLIP_U1_GS_GS();

  // test function
  double thetak = 0; double phik = 0; double mLlA=1;
  k_value = 0.249748;
  alpha_lam = 0.524626;
  alpha_rho = 0.413255;
  mbottom = 4.928;
  mlight = 0.382;


  double tensor1 = T1l(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  double tensor2 = T2l(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);

  std::cout<<tensor1<<std::endl;

  return decayWidth * test_integral;
}

double EMDecayWidths::DecayWidth(double flav_coup, double fi2_value, double angular_sum_value){
  double GeV = 1000.;
  double decayWidth = flav_coup  * fi2_value * (1./(2 * JA + 1)) * angular_sum_value;
  return decayWidth * GeV;
}

std::vector<double> EMDecayWidths::getMomentumProjections(double j_angular){
  // Method to obtain the projections "m_projection" for a given angular momentum "j_angular"
  std::vector<double> angularProjections; angularProjections.clear();
  if(j_angular==0.){angularProjections.push_back(0); return angularProjections;}  
  double m_projection = (-1.0)*j_angular;
  do{
    angularProjections.push_back(m_projection);
    m_projection++;
  }while(m_projection<=j_angular);  
  return angularProjections;
}

int EMDecayWidths::KroneckerDelta(double i, double j){
  if(i==j) return 1;
  else return 0;
}

double EMDecayWidths::EB(double MB, double K){
  double value = std::sqrt(std::pow(MB, 2) - std::pow(K, 2));
  return value;
}

double EMDecayWidths::K(double MA, double MB){
  double value = (0.5)*(std::pow(MA, 2) - std::pow(MB, 2))/ MA;
  return value;
}

double EMDecayWidths::FI2(double EB, double MA, double k_value){
  double value = 4.* pi_val * (EB/MA) * std::pow(k_value, 2);
  return value;
}

double EMDecayWidths::ClebshGordan(WignerSymbols *m_wigner,
				   double l1, double l2, double l3,
				   double m1, double m2, double m3){
  double coef = std::pow(-1.0, m2-m1-m3) * std::pow(2*l3 + 1, 0.5);
  double three_j = m_wigner->wigner3j(l1, l2, l3, m1, m2, (-1.0)*m3);
  return coef * three_j;
}

double EMDecayWidths::ANGULAR_SUM(double alpha_rho, double alpha_lam, double k_value){  
  WignerSymbols  *m_wigner = new WignerSymbols();
  return 1.0;
}

// SPIN-FLIP INTEGRALS
double EMDecayWidths::SPINFLIP_U1_GS_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8.;
  double value2 = 3. * std::pow(mbottom, 2) / (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = 1./std::pow(alpha_rho, 2);
  double value = std::exp(value1 * (value2 + value3));
  return value;
}

double EMDecayWidths::SPINFLIP_U2_GS_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8.;
  double value2 = 3. * std::pow(mbottom, 2) / (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = 1./std::pow(alpha_rho, 2);
  double value = std::exp(value1 * (value2 + value3));
  return value;
}

double EMDecayWidths::SPINFLIP_U3_GS_GS(double k_value, double alpha_lam,  double mbottom, double mlight){
  double value1 = (-3.0) * std::pow(k_value, 2) * std::pow(mlight, 2) / 2. * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = std::exp(value1);
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = std::sqrt(6) * p_imag * mbottom * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = std::sqrt(3) * p_imag * mbottom * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-1.0) * std::sqrt(6) * p_imag * mbottom * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = std::sqrt(6) * p_imag * mbottom * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (-1.0) * std::sqrt(3) * p_imag * mbottom * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-1.0) * std::sqrt(6) * p_imag * mbottom * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = p_imag * phik;
  double value2 = (-3.0) * std::pow(mlight, 2) * std::pow(k_value, 2) / 2 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (-1.0) * std::sqrt(6) * p_imag * mlight * k_value * std::exp(value1 + value2) * std::sin(thetak)/ (2 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-3.0) * std::pow(mlight, 2) * std::pow(k_value, 2) / 2 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = std::sqrt(3) * p_imag * mlight * k_value * std::exp(value1) * std::cos(thetak)/ (alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * p_imag * phik;
  double value2 = (-3.0) * std::pow(mlight, 2) * std::pow(k_value, 2) / 2 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = std::sqrt(6) * p_imag * mlight * k_value * std::exp(value1 + value2) * std::sin(thetak)/ (2 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = std::sqrt(2) * p_imag * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value =  (-1.0) * p_imag * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-1.0) * std::sqrt(2) * p_imag * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = (-1.0) * std::sqrt(2) * p_imag * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value =  p_imag * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = std::sqrt(2) * p_imag * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  k_value=0; alpha_lam=0; alpha_rho=0; mbottom=0; mlight=0; phik=0;  thetak=0;
  double value = 0;
  return value * k_value * alpha_lam * alpha_rho * mbottom * mlight * phik * thetak;
}

double EMDecayWidths::SPINFLIP_U3_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  k_value=0; alpha_lam=0; alpha_rho=0; mbottom=0; mlight=0; phik=0;  thetak=0;
  double value = 0;
  return value * k_value * alpha_lam * alpha_rho * mbottom * mlight * phik * thetak;
}

double EMDecayWidths::SPINFLIP_U3_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  k_value=0; alpha_lam=0; alpha_rho=0; mbottom=0; mlight=0; phik=0;  thetak=0;
  double value = 0;
  return value * k_value * alpha_lam * alpha_rho * mbottom * mlight * phik * thetak;
}

// ORBIT-SPLIT INTEGRALS
// U1_1lambda-1lambda
double EMDecayWidths::ORBITALSPLIT_U1_1l_m1_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value = std::exp(value1 + value2) * (1 - value2 * std::pow(std::sin(thetak), 2)) ;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = 3. * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/ 16 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = 2 * p_imag * phik;
  double value = 0.375 * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m0_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = 3. * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/ 16 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (((-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 4 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2))) * std::pow(std::cos(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m0_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = (-1.0) * 0.1875 * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/(std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
} 

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1m_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = 2 * p_imag * phik;
  double value = 3 * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1m_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-3.0) * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/16 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1m_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (value2 * std::pow(std::sin(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

// U2_1lambda-1lambda
double EMDecayWidths::ORBITALSPLIT_U2_1l_m1_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (value2 * std::pow(std::sin(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = 3. * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/ 16 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = 2 * p_imag * phik;
  double value = 0.375 * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = 3. * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/ 16 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (((-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 4 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2))) * std::pow(std::cos(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = (-1.0) * 0.1875 * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/(std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
} 

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1m_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = 2 * p_imag * phik;
  double value = 3 * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1m_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-3.0) * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/16 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1m_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (value2 * std::pow(std::sin(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}


// New integrals
//U3 1lam->1lam;
double EMDecayWidths::ORBITALSPLIT_U3_1l_m1_1l_m1(double k_value, double alpha_lam, double mbottom, double mlight, double thetak){

  double value1 = (-3.) * std::pow(k_value,2) * std::pow(mlight,2) + 2 * std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2)
    + 3 * std::pow(k_value,2) * std::pow(mlight,2) * std::pow(std::cos(thetak),2);
  double value2 =  (2.) * std::pow(alpha_lam,2) * std::exp((3 * std::pow(k_value,2) * std::pow(mlight,2))/(2. * std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))) * std::pow(mbottom + 2 * mlight,2);
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U3_1l_m0_1l_m0(double k_value, double alpha_lam, double mbottom, double mlight, double thetak){
  double value1 = std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2) - 3 * std::pow(k_value,2) * std::pow(mlight,2) * std::pow(std::cos(thetak),2);
  double value2 = std::pow(alpha_lam,2) * std::exp((3 * std::pow(k_value,2) * std::pow(mlight,2)) / (2. * std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))) * std::pow(mbottom + 2 * mlight,2);
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U3_1l_m1m_1l_m1m(double k_value, double alpha_lam, double mbottom, double mlight, double thetak){

  double value1 = -3 * std::pow(k_value,2) * std::pow(mlight,2) + 2 * std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2) + 3 * std::pow(k_value,2) * std::pow(mlight,2) * std::pow(std::cos(thetak),2);
  double value2 = 2. * std::pow(alpha_lam,2) * std::exp((3 * std::pow(k_value,2) * std::pow(mlight,2))/(2. * std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))) * std::pow(mbottom + 2 * mlight,2);
  double value = value1/value2;
  return value;
}

//U1 1lam->1rho
double EMDecayWidths::ORBITALSPLIT_U1_1l_m0_1r_m0(double k_value, double alpha_rho, double alpha_lam, double mbottom, double mlight, double thetak){
  double value1 = std::pow(k_value,2) * mbottom * std::pow(std::cos(thetak),2);
  double value2 = std::exp((std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(mbottom,2))/(std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))))/8.);
  double value2p = (4 * alpha_lam * alpha_rho * mbottom + 8 * alpha_lam * alpha_rho * mlight);
  value2 = value2 * value2p;
  double value = value1/value2;
  return (-1.0) * std::sqrt(3) * value;
}

//U2 1lam->1rho
double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1r_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = std::pow(k_value,2) * mbottom * std::pow(std::cos(thetak),2);
  double value2 = std::exp((std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(mbottom,2))/(std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))))/8.);
  double value2p = 4 * alpha_lam * alpha_rho * mbottom + 8 * alpha_lam * alpha_rho * mlight;
  value2 = value2 * value2p;
  double value = value1/value2;
  return (-1.0) * std::sqrt(3) * value;
}


//U1 1rho->1rho
double EMDecayWidths::ORBITALSPLIT_U1_1r_m1_1r_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = 8 * std::pow(alpha_rho,2) - std::pow(k_value,2) + std::pow(k_value,2) * std::pow(std::cos(thetak),2);
  double value2 = 8. * std::pow(alpha_rho,2);
  double value2p = std::exp((std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(mbottom,2))/(std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1r_m0_1r_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = 4 * std::pow(alpha_rho,2) - std::pow(k_value,2) * std::pow(std::cos(thetak),2);
  double value2 = 4. * std::pow(alpha_rho,2);
  double value2p = std::exp((std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(mbottom,2))/(std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1r_m1m_1r_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = 8 * std::pow(alpha_rho,2) - std::pow(k_value,2) + std::pow(k_value,2) * std::pow(std::cos(thetak),2);
  double value2 = 8. * std::pow(alpha_rho,2);
  double value2p = std::exp((std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(mbottom,2))/(std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}

//U2 1rho->1rho
double EMDecayWidths::ORBITALSPLIT_U2_1r_m1_1r_m1(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 =(8 * std::pow(alpha_rho,2) - std::pow(k_value,2) + std::pow(k_value,2) * std::pow(std::cos(thetak),2));
  double value2 = 8. * std::pow(alpha_rho,2);
  double value2p = std::exp((std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(mbottom,2))/(std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1r_m0_1r_m0(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = (4 * std::pow(alpha_rho,2) - std::pow(k_value,2) * std::pow(std::cos(thetak),2));
  double value2 = 4. * std::pow(alpha_rho,2);
  double value2p = std::exp((std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(mbottom,2))/(std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1r_m1m_1r_m1m(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = 8 * std::pow(alpha_rho,2) - std::pow(k_value,2) + std::pow(k_value,2) * std::pow(std::cos(thetak),2);
  double value2 = 8 * std::pow(alpha_rho,2);
  double value2p = std::exp((std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(mbottom,2))/(std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}

//U3 1rho->1rho
double EMDecayWidths::ORBITALSPLIT_U3_1r_m1_1r_m1(double k_value, double alpha_lam, double mbottom, double mlight){
  double value = std::exp((-3 * std::pow(k_value,2) * std::pow(mlight,2))/(2. * std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2)));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U3_1r_m0_1r_m0(double k_value, double alpha_lam, double mbottom, double mlight){
  double value = std::exp((-3 * std::pow(k_value,2) * std::pow(mlight,2))/(2. * std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2)));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U3_1r_m1m_1r_m1m(double k_value, double alpha_lam, double mbottom, double mlight){
  double value = std::exp((-3 * std::pow(k_value,2) * std::pow(mlight,2))/(2. * std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2)));
  return value;
}

//U1 1rho->1lam
double EMDecayWidths::ORBITALSPLIT_U1_1r_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,
						  double mbottom, double mlight, double thetak, double phik){
  double value1 = std::exp(-std::pow(k_value,2)/(8. * std::pow(alpha_rho,2)) -
			   (3 * std::pow(k_value,2) * std::pow(mbottom,2))/(8. * std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2)) + p_imag * phik);
  double value1p = std::pow(k_value,2) * mbottom * std::cos(thetak) * std::sin(thetak);
  value1 = value1 * value1p;
  double value2 = 4 * alpha_lam * alpha_rho * mbottom + 8 * alpha_lam * alpha_rho * mlight;
  double value = value1/value2;
  return std::sqrt(1.5) * value;
}

//U2 1rho->1lam
double EMDecayWidths::ORBITALSPLIT_U2_1r_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,
						  double mbottom, double mlight, double thetak, double phik){
  double value1 = std::exp(-(std::pow(k_value,2) * (std::pow(alpha_rho,-2) +
						    (3 * std::pow(mbottom,2))/(std::pow(alpha_lam,2) * std::pow(mbottom + 2 * mlight,2))))/8. + p_imag * phik);
  double value1p = std::pow(k_value,2) * mbottom * std::cos(thetak) * std::sin(thetak);
  value1 = value1 * value1p;
  double value2 = 4 * alpha_lam * alpha_rho * mbottom + 8 * alpha_lam * alpha_rho * mlight;
  double value = value1/value2;
  return -std::sqrt(1.5) * value;
}

//U1 2lam->gs (13.03.2023)
double EMDecayWidths::ORBITALSPLIT_U1_2l_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = std::pow(k_value,2)*std::pow(mbottom,2)*(1 + 3*std::cos(2*thetak));///
  double value2 = 16.*std::pow(alpha_lam,2)*std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
  double value2p = std::pow(mbottom + 2*mlight,2);  
  value2 = value2 * value2p;
  double value = value1/value2;
  return (-1.0) * std::sqrt(3) * value;
}

//U2 2lam->gs
double EMDecayWidths::ORBITALSPLIT_U2_2l_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = std::pow(k_value,2)*std::pow(mbottom,2)*(1 + 3*std::cos(2*thetak));///
  double value2 = 16.*std::pow(alpha_lam,2)*std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
  double value2p = std::pow(mbottom + 2*mlight,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return (-1.0) * std::sqrt(3) * value;
}

// U3 2lam->gs
double EMDecayWidths::ORBITALSPLIT_U3_2l_m0_GS(double k_value, double alpha_lam, double mbottom, double mlight, double thetak){
  double value1 = std::pow(k_value,2)*std::pow(mlight,2)*(1 + 3*std::cos(2*thetak));
  double value2 = 4.*std::pow(alpha_lam,2)*std::exp((3*std::pow(k_value,2)*std::pow(mlight,2))/(2.*std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2)));
  double value2p = std::pow(mbottom + 2*mlight,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return (-1.0) * std::sqrt(3) * value;  
}

// U1 2rho->gs
double EMDecayWidths::ORBITALSPLIT_U1_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = std::pow(k_value,2)*(1 - 3*std::pow(std::cos(thetak),2));///
  double value2 = 8.*std::sqrt(3)*std::pow(alpha_rho,2)*std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
  double value = value1/value2;
  return value;
}

// U2 2rho->gs
double EMDecayWidths::ORBITALSPLIT_U2_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight, double thetak){
  double value1 = std::pow(k_value,2)*(1 - 3*std::pow(std::cos(thetak),2));///
  double value2 = 8.*std::sqrt(3)*std::pow(alpha_rho,2)*std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
  double value = value1/value2;
  return value;
}
  
// U1 1nlam->gs
double EMDecayWidths::ORBITALSPLIT_U1_1nl_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight){
  double value1 = std::pow(k_value,2)*std::pow(mbottom,2);///
  double value2 = 4.*std::pow(alpha_lam,2) * std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
  double value2p = std::pow(mbottom + 2*mlight,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return std::sqrt(1.5)*value;
}

// U2 1nlam->gs
double EMDecayWidths::ORBITALSPLIT_U2_1nl_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight){
  double value1 = std::pow(k_value,2)*std::pow(mbottom,2);///  
  double value2 = 4.*std::pow(alpha_lam,2) * std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
  double value2p = std::pow(mbottom + 2*mlight,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return std::sqrt(1.5)*value;
}

// U3 1nlam->gs
double EMDecayWidths::ORBITALSPLIT_U3_1nl_m0_GS(double k_value, double alpha_lam, double mbottom, double mlight){
  double value1 = std::pow(k_value,2)*std::pow(mlight,2);///
  double value2 = std::pow(alpha_lam,2)*std::exp((3*std::pow(k_value,2)*std::pow(mlight,2))/(2.*std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2)));
  double value2p = std::pow(mbottom + 2*mlight,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return std::sqrt(1.5)*value;
}


// U1 1nrho->gs
double EMDecayWidths::ORBITALSPLIT_U1_1nr_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight){
  double value1 = std::pow(k_value,2);///
  double value2 = 4.*std::sqrt(6)*std::pow(alpha_rho,2)*std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
  double value = value1/value2;
  return value;
}

// U2 1nrho->gs
double EMDecayWidths::ORBITALSPLIT_U2_1nr_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight){
  double value1 = std::pow(k_value,2);///
  double value2 = 4.*std::sqrt(6)*std::pow(alpha_rho,2)*std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
  double value = value1/value2;
  return value;
}

//Tensor operators
//T1l
double EMDecayWidths::T1l(double k_value, double alpha_lam, double alpha_rho,
			  double mbottom, double mlight, double thetak, double phik, double mLlA){
  double value1 = ORBITALSPLIT_U1_1l_m1_1l_m1(k_value, alpha_lam, alpha_rho, mbottom, mlight, phik, thetak);
  double value2 = ORBITALSPLIT_U1_2l_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value3 = SPINFLIP_U1_GS_GS(k_value, alpha_lam, alpha_rho,  mbottom, mlight);
  double value4 = ORBITALSPLIT_U1_1nl_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value = value1 + std::sqrt(1/3) * value2 +  value3 +  std::sqrt(2/3) * value4;
  std::cout<<value1<<"   "<<value2<<"   "<<value3<<"   "<<value4<<"   "<<"   "<<value<<std::endl;
  return p_imag * std::sqrt(1/6) * alpha_lam * value;
}

//T2l
double EMDecayWidths::T2l(double k_value, double alpha_lam, double alpha_rho,
			  double mbottom, double mlight, double thetak, double phik, double mLlA){
  double value1 = ORBITALSPLIT_U2_1l_m1_1l_m1(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value2 = ORBITALSPLIT_U2_2l_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value3 = SPINFLIP_U2_GS_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value4 = ORBITALSPLIT_U2_1nl_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value = value1 + std::sqrt(1/3) * value2 +  value3 +  std::sqrt(2/3) * value4;
  std::cout<<value1<<"   "<<value2<<"   "<<value3<<"   "<<value4<<"   "<<"   "<<value<<std::endl;
  return p_imag * std::sqrt(1/6) * alpha_lam * value;
}

//T3l
double EMDecayWidths::T3l(double k_value, double alpha_lam, double alpha_rho,
			  double mbottom, double mlight, double thetak, double phik, double mLlA){
  double value1 = ORBITALSPLIT_U3_1l_m1_1l_m1(k_value, alpha_lam, mbottom, mlight, thetak);
  double value2 = ORBITALSPLIT_U3_2l_m0_GS(k_value, alpha_lam, mbottom, mlight, thetak);
  double value3 = SPINFLIP_U3_GS_GS(k_value, alpha_lam, mbottom, mlight);
  double value4 = ORBITALSPLIT_U3_1nl_m0_GS(k_value, alpha_lam, mbottom, mlight);
  double value = value1 + std::sqrt(1/3) * value2 +  value3 +  std::sqrt(2/3) * value4;
  return (-1.0) * p_imag * std::sqrt(2/3) * alpha_lam * value;
}

//T1r
double EMDecayWidths::T1r(double k_value, double alpha_lam, double alpha_rho,
			  double mbottom, double mlight, double thetak, double phik, double mLrA){
  double value1 = ORBITALSPLIT_U1_1r_m1_1r_m1(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value2 = ORBITALSPLIT_U1_2r_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value3 = SPINFLIP_U1_GS_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value4 = ORBITALSPLIT_U1_1nr_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value = value1 + std::sqrt(1/3) * value2 +  value3 +  std::sqrt(2/3) * value4;
  return p_imag * std::sqrt(1/2) * alpha_rho * value;
}

//T2r
double EMDecayWidths::T2r(double k_value, double alpha_lam, double alpha_rho,
			  double mbottom, double mlight, double thetak, double phik, double mLrA){
  double value1 = ORBITALSPLIT_U2_1r_m1_1r_m1(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value2 = ORBITALSPLIT_U2_2r_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value3 = SPINFLIP_U2_GS_GS(k_value, alpha_lam, alpha_rho,  mbottom, mlight);
  double value4 = ORBITALSPLIT_U2_1nr_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value = value1 + std::sqrt(1/3) * value2 +  value3 +  std::sqrt(2/3) * value4;
  return (-1.0) * p_imag * std::sqrt(1/2) * alpha_rho * value;
}

//T3r
double EMDecayWidths::T3r(){
  double value = 0;
  return value;
}

#endif
