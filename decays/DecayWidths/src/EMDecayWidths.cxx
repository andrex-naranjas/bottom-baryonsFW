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

double EMDecayWidths::execute(double ma_val, double mb_val, double mc_val, double ga_val, double sa_val,
				 double la_val, double ja_val, double sl_val, double al_val, double ar_val,
				 int baryon, int excMode, int prodDecay){  
  // decay product masses
  MA = ma_val;
  MB = mb_val;
  MC = mc_val;
  if(MA<MB+MC) return 0.; //energy conservation

  // which baryon, mode, decay product
  baryonFlag = baryon;
  modeExcitation = excMode;
  decayProd  = prodDecay;
  double alpha_rho = 0.,alpha_lam = 0.,alpha_mes = 0.,flav_coup= 0.;
  double slf_val=0., sb_val=0., lb_val=0., jb_val=0.;
  double l1_val=0., l2_val=0.;

  double gamma = ga_val;

  alpha_rho = ar_val;
  alpha_lam = al_val;
  alpha_mes = ALPHA_MES(diagram);

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
  s  = 1.0;   m   = getMomentumProjections(s);
  s1 = 0.5;   m1  = getMomentumProjections(s1);
  s2 = 0.5;   m2  = getMomentumProjections(s2);
  s3 = 0.5;   m3  = getMomentumProjections(s3);
  s4 = 0.5;   m4  = getMomentumProjections(s4);
  s5 = 0.5;   m5  = getMomentumProjections(s5);

  double EB_value = EB(MA,MB,MC);
  double k_value; k_value = K(EB_value, MB);
  double EWCC_value = EWCC(MA, MB, MC);
  
  double sum_value  = 0.; //test
  sum_value = ANGULAR_SUM(alpha_rho, alpha_lam, alpha_mes, k_value);

  double fi2_value  = FI2(EB_value, EWCC_value, MA, k_value);
  double decayWidth = DecayWidth(flav_coup, gamma, fi2_value, sum_value);
  decayWidth = 1;

  test_integral = SPINFLIP_U1_GS_GS();

  return decayWidth * test_integral;
}

double EMDecayWidths::DecayWidth(double flav_coup, double gamma, double fi2_value, double angular_sum_value){
  double GeV = 1000.;
  double decayWidth = flav_coup * std::pow(gamma, 2) * fi2_value * (1./(2 * JA + 1)) * angular_sum_value;
  return decayWidth * GeV;
}


double EMDecayWidths::ANGULAR_SUM(double alpha_rho, double alpha_lam,
				  double alpha_mes, double k_value){
  
  WignerSymbols  * m_wigner = new WignerSymbols();
  return 1.0;
}


// SPIN-FLIP INTEGRALS
double EMDecayWidths::SPINFLIP_U1_GS_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8.;
  double value2 = 3. * std::pow(MB, 2) / (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = 1./std::pow(alpha_rho, 2);
  double value = std::exp(value1 * (value2 + value3));
  return value;
}

double EMDecayWidths::SPINFLIP_U2_GS_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8.;
  double value2 = 3. * std::pow(MB, 2) / (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = 1./std::pow(alpha_rho, 2);
  double value = std::exp(value1 * (value2 + value3));
  return value;
}

double EMDecayWidths::SPINFLIP_U3_GS_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML){
  double value1 = (-3.0) * std::pow(k_value, 2) * std::pow(ML, 2) / 2. * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = std::exp(value1);
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = i * phik;
  double value = std::sqrt(6) * i * MB * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (MB + 2 * ML) );
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = std::sqrt(3) * i * MB * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_lam  * (MB + 2 * ML) );
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = (-1.0) * i * phik;
  double value = (-1.0) * std::sqrt(6) * i * MB * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (MB + 2 * ML) );
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = i * phik;
  double value = std::sqrt(6) * i * MB * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (MB + 2 * ML) );
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = (-1.0) * std::sqrt(3) * i * MB * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_lam  * (MB + 2 * ML) );
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = (-1.0) * i * phik;
  double value = (-1.0) * std::sqrt(6) * i * MB * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (MB + 2 * ML) );
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = i * phik;
  double value2 = (-3.0) * std::pow(ML, 2) * std::pow(k_value, 2) / 2 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = (-1.0) * std::sqrt(6) * i * ML * k_value * std::exp(value1 + value2) * std::sin(thetak)/ (2 * alpha_lam  * (MB + 2 * ML) );
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-3.0) * std::pow(ML, 2) * std::pow(k_value, 2) / 2 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = std::sqrt(3) * i * ML * k_value * std::exp(value1) * std::cos(thetak)/ (alpha_lam  * (MB + 2 * ML) );
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * i * phik;
  double value2 = (-3.0) * std::pow(ML, 2) * std::pow(k_value, 2) / 2 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = std::sqrt(6) * i * ML * k_value * std::exp(value1 + value2) * std::sin(thetak)/ (2 * alpha_lam  * (MB + 2 * ML) );
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = i * phik;
  double value = std::sqrt(2) * i * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value =  (-1.0) * i * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = (-1.0) * i * phik;
  double value = (-1.0) * std::sqrt(2) * i * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = i * phik;
  double value = (-1.0) * std::sqrt(2) * i * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value =  i * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = (-1.0) * i * phik;
  double value = std::sqrt(2) * i * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value = 0;
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value = 0;
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value = 0;
  return value;
}

// ORBIT-SPLIT INTEGRALS
// U1_1lambda-1lambda
double EMDecayWidths::ORBITALSPLIT_U1_1l_m1_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = (value2 * std::pow(std::sin(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = i * phik;
  double value = 3. * std::sqrt(2) * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/ 16 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = 2 * i * phik;
  double value = 0.375 * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m0_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = (-1.0) * i * phik;
  double value = 3. * std::sqrt(2) * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/ 16 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = (((-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 4 * (std::pow(alpha_lam * (MB + 2. * ML), 2))) * std::pow(std::cos(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m0_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = i * phik;
  double value = (-1.0) * 0.1875 * std::sqrt(2) * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/(std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
} 

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1m_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = 2 * i * phik;
  double value = 3 * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1m_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = (-1.0) * i * phik;
  double value = (-3.0) * std::sqrt(2) * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/16 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1m_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = (value2 * std::pow(std::sin(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

// U2_1lambda-1lambda
double EMDecayWidths::ORBITALSPLIT_U2_1l_m1_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = (value2 * std::pow(std::sin(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = i * phik;
  double value = 3. * std::sqrt(2) * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/ 16 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = 2 * i * phik;
  double value = 0.375 * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = (-1.0) * i * phik;
  double value = 3. * std::sqrt(2) * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/ 16 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = (((-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 4 * (std::pow(alpha_lam * (MB + 2. * ML), 2))) * std::pow(std::cos(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = i * phik;
  double value = (-1.0) * 0.1875 * std::sqrt(2) * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/(std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
} 

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1m_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = 2 * i * phik;
  double value = 3 * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1m_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value3 = (-1.0) * i * phik;
  double value = (-3.0) * std::sqrt(2) * std::pow(MB, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/16 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1m_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double MB, double ML, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(MB, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (MB + 2. * ML), 2));
  double value = (value2 * std::pow(std::sin(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}


// New integrals
//U3 1lam->1lam;
double EMDecayWidths::ORBITALSPLIT_U3_1l_m1_1l_m1(double k_value, double alpha_lam, double MB, double ML, double thetak){

  double value1 = (-3.) * std::pow(k_value,2) * std::pow(ML,2) + 2 * std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2)
    + 3 * std::pow(k_value,2) * std::pow(ML,2) * std::pow(std::cos(thetak),2);
  double value2 =  (2.) * std::pow(alpha_lam,2) * std::pow(E,(3 * std::pow(k_value,2) * std::pow(ML,2))/(2. * std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))) * std::pow(MB + 2 * ML,2);
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U3_1l_m0_1l_m0(double k_value, double alpha_lam, double MB, double ML, double thetak){
  double value1 = std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2) - 3 * std::pow(k_value,2) * std::pow(ML,2) * std::pow(std::cos(thetak),2);
  double value2 = std::pow(alpha_lam,2) * std::pow(E,(3 * std::pow(k_value,2) * std::pow(ML,2)) / (2. * std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))) * std::pow(MB + 2 * ML,2);
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U3_1l_m1m_1l_m1m(double k_value, double alpha_lam, double MB, double ML, double thethak){

  double value1 = -3 * std::pow(k_value,2) * std::pow(ML,2) + 2 * std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2) + 3 * std::pow(k_value,2) * std::pow(ML,2) * std::pow(std::cos(thetak),2);
  double value2 = 2. * std::pow(alpha_lam,2) * std::pow(E,(3 * std::pow(k_value,2) * std::pow(ML,2))/(2. * std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))) * std::pow(MB + 2 * ML,2);
  double value = value1/value2;
  return value;
}

//U1 1lam->1rho
double EMDecayWidths::ORBITALSPLIT_U1_1l_m0_rl_m0(double k_value, double alpha_rho, double alpha_lam, double MB, double ML){
  double value1 = std::pow(k_value,2) * MB * std::pow(std::cos(thetak),2);
  double value2 = std::pow(E,(std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(MB,2))/(std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))))/8.);
  double value2p = (4 * alpha_lam * alpha_rho * MB + 8 * alpha_lam * alpha_rho * ML);
  value2 = value2 * value2p;
  double value = value1/value2;
  return -Sqrt(3) * value;
}

//U2 1lam->1rho
double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_rl_m0(double k_value, double alpha_lam, double alpha_rho, double MB, double ML){
  double value1 = std::pow(k_value,2) * MB * std::pow(std::cos(thetak),2);
  double value2 = std::pow(E,(std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(MB,2))/(std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))))/8.);
  double value2p = 4 * alpha_lam * alpha_rho * MB + 8 * alpha_lam * alpha_rho * ML;
  value2 = value2 * value2p;
  double value = value1/value2;
return -Sqrt(3) * value;
}


//U1 1rho->1rho
double EMDecayWidths::ORBITALSPLIT_U1_1r_m1_rl_m1(double k_value, double alpha_lam, double alpha_rho, double MB, double ML){
  double value1 = 8 * std::pow(alpha_rho,2) - std::pow(k_value,2) + std::pow(k_value,2) * std::pow(std::cos(thetak),2);
  double value2 = 8. * std::pow(alpha_rho,2);
  double value2p = std::pow(E,(std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(MB,2))/(std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1r_m0_rl_m0(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = 4 * std::pow(alpha_rho,2) - std::pow(k_value,2) * std::pow(std::cos(thetak),2);
  double value2 = 4. * std::pow(alpha_rho,2);
  double value2p = std::pow(E,(std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(MB,2))/(std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1r_m1m_rl_m1m(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = 8 * std::pow(alpha_rho,2) - std::pow(k_value,2) + std::pow(k_value,2) * std::pow(std::cos(thetak),2);
  double value2 = 8. * std::pow(alpha_rho,2)
    double value2p = std::pow(E,(std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(MB,2))/(std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}


//U2 1rho->1rho
double EMDecayWidths::ORBITALSPLIT_U2_1r_m1_rl_m1(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 =(8 * std::pow(alpha_rho,2) - std::pow(k_value,2) + std::pow(k_value,2) * std::pow(std::cos(thetak),2));
  double value2 = 8. * std::pow(alpha_rho,2);
  double value2p = std::pow(E,(std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(MB,2))/(std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1r_m0_rl_m0(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = (4 * std::pow(alpha_rho,2) - std::pow(k_value,2) * std::pow(std::cos(thetak),2));
  double value2 = 4. * std::pow(alpha_rho,2);
  double value2p = std::pow(E,(std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(MB,2))/(std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}


double EMDecayWidths::ORBITALSPLIT_U2_1r_m1m_rl_m1m(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = 8 * std::pow(alpha_rho,2) - std::pow(k_value,2) + std::pow(k_value,2) * std::pow(std::cos(thetak),2);
  double value2 = 8 * std::pow(alpha_rho,2);
  double value2p = std::pow(E,(std::pow(k_value,2) * (std::pow(alpha_rho,-2) + (3 * std::pow(MB,2))/(std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))))/8.);
  value2 = value2 * value2p;
  double value = value1/value2;
  return value;
}


//U3 1rho->1rho
double EMDecayWidths::ORBITALSPLIT_U3_1r_m1_rl_m1(double k_value, double alpha_lam, double MB, double ML){
  double value = std::pow(E,(-3 * std::pow(k_value,2) * std::pow(ML,2))/(2. * std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2)));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U3_1r_m0_rl_m0(double k_value, double alpha_lam, double MB, double ML){
  double value = std::pow(E,(-3 * std::pow(k_value,2) * std::pow(ML,2))/(2. * std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2)));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U3_1r_m1m_rl_m1m(double k_value, double alpha_lam, double MB, double ML){
  double value = std::pow(E,(-3 * std::pow(k_value,2) * std::pow(ML,2))/(2. * std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2)));
  return value;
}


//U1 1rho->1lam
double EMDecayWidths::ORBITALSPLIT_U1_1r_m0_ll_m0(double k_value, double alpha_lam, double alpha_rho,
						  double MB, double ML, double thetak, double phik){

  double value1 = std::pow(E,-std::pow(k_value,2)/(8. * std::pow(alpha_rho,2)) -
			   (3 * std::pow(k_value,2) * std::pow(MB,2))/(8. * std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2)) + Complex(0,1) * phik);
  double value1p = std::pow(k_value,2) * MB * std::cos(thetak) * std::sin(thetak);
  value1 = value1 * value1p;
  double value2 = 4 * alpha_lam * alpha_rho * MB + 8 * alpha_lam * alpha_rho * ML;
  double value = value1/value2;
  return Sqrt(1.5) * value;
}


//U2 1rho->1lam
double EMDecayWidths::ORBITALSPLIT_U2_1r_m0_ll_m0(double k_value, double alpha_lam, double alpha_rho,
						  double MB, double ML, double thetak, double phik){

  double value1 = std::pow(E,-(std::pow(k_value,2) * (std::pow(alpha_rho,-2) +
						   (3 * std::pow(MB,2))/(std::pow(alpha_lam,2) * std::pow(MB + 2 * ML,2))))/8. + Complex(0,1) * phik);
  double value1p = std::pow(k_value,2) * MB * std::cos(thetak) * std::sin(thetak);
  double value2 = 4 * alpha_lam * alpha_rho * MB + 8 * alpha_lam * alpha_rho * ML;
  double value = value1/value2;
  return -Sqrt(1.5) * value;
}


//U1 2lam->gs (13.03.2023)
double EMDecayWidths::ORBITALSPLIT_U1_2l_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2)*std::pow(MB,2)*(1 + 3*std::cos(2*thetak));///
  double value2 = 16.*std::pow(alpha_lam,2)*std::pow(E,(std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(MB,2))/(std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2))))/8.);
  double value2p = std::pow(MB + 2*ML,2);  
  value2 = value2 * value2p;
  double value = value1/value2;
  return (-1.0) * Sqrt(3) * value;
}

//U2 2lam->gs
double EMDecayWidths::ORBITALSPLIT_U2_2l_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2)*std::pow(MB,2)*(1 + 3*std::cos(2*thetak));///
  double value2 = 16.*std::pow(alpha_lam,2)*std::pow(E,(std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(MB,2))/(std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2))))/8.);
  double value2p = std::pow(MB + 2*ML,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return (-1.0) * Sqrt(3) * value;
}

// U3 2lam->gs
double EMDecayWidths::ORBITALSPLIT_U3_2l_m0_GS(double k_value, double alpha_lam, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2)*std::pow(ML,2)*(1 + 3*std::cos(2*thetak));
  double value2 = 4.*std::pow(alpha_lam,2)*std::pow(E,(3*std::pow(k_value,2)*std::pow(ML,2))/(2.*std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2)));
  double value2p = std::pow(MB + 2*ML,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return (-1.0) * Sqrt(3) * value;  
}

// U1 2rho->gs
double EMDecayWidths::ORBITALSPLIT_U1_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2)*(1 - 3*std::pow(std::cos(thetak),2));///
  double value2 = 8.*Sqrt(3)*std::pow(alpha_rho,2)*std::pow(E,(std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(MB,2))/(std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2))))/8.);
  double value = value1/value2;
  return value;
}

// U2 2rho->gs
double EMDecayWidths::ORBITALSPLIT_U2_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2)*(1 - 3*std::pow(std::cos(thetak),2));///
  double value2 = 8.*Sqrt(3)*std::pow(alpha_rho,2)*std::pow(E,(std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(MB,2))/(std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2))))/8.);
  double value = value1/value2;
  return value;
}
  
// U1 1nlam->gs
double EMDecayWidths::ORBITALSPLIT_U1_2nl_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2)*std::pow(MB,2);///
  double value2 = 4.*std::pow(alpha_lam,2) * std::pow(E,(std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(MB,2))/(std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2))))/8.);
  double value2p = std::pow(MB + 2*ML,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return Sqrt(1.5)*value;
}

// U2 1nlam->gs
double EMDecayWidths::ORBITALSPLIT_U2_1nl_m0_GS_2r_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2)*std::pow(MB,2);///  
  double value2 = 4.*std::pow(alpha_lam,2) * std::pow(E,(std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(MB,2))/(std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2))))/8.);
  double value2p = std::pow(MB + 2*ML,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return Sqrt(1.5)*value;
}

// U3 1nlam->gs
double EMDecayWidths::ORBITALSPLIT_U3_1nl_m0_GS(double k_value, double alpha_lam, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2)*std::pow(ML,2);///
  double value2 = std::pow(alpha_lam,2)*std::pow(E,(3*std::pow(k_value,2)*std::pow(ML,2))/(2.*std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2)));
  double value2p = std::pow(MB + 2*ML,2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return Sqrt(1.5)*value;
}


// U1 1nrho->gs
double EMDecayWidths::ORBITALSPLIT_U1_1nr_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2);///
  double value2 = 4.*Sqrt(6)*std::pow(alpha_rho,2)*std::pow(E,(std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(MB,2))/(std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2))))/8.);
  double value = value1/value2;
  return value;
}

// U2 1nrho->gs
double EMDecayWidths::ORBITALSPLIT_U2_1nr_m0_GS(double k_value, double alpha_lam, double alpha_rho, double MB, double ML, double thetak){
  double value1 = std::pow(k_value,2);///
  double value2 = 4.*Sqrt(6)*std::pow(alpha_rho,2)*std::pow(E,(std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(MB,2))/(std::pow(alpha_lam,2)*std::pow(MB + 2*ML,2))))/8.);
  double value = value1/value2;
  return value;
}


//Tensor operators
//T1l
double EMDecayWidths::T1l(double k_value, double alpha_lam, double alpha_rho,
			  double MB, double ML, double thetak, double phik, double mLlA){

  double value1 = ORBITALSPLIT_U1_1l_m1_1l_m1(k_value, alpha_lam, alpha_rho, MB, ML, phik, thetak);
  double value2 = ORBITALSPLIT_U1_1l_m1_1l_m1(k_value, alpha_lam, alpha_rho, MB, ML, phik, thetak);
  double value = value1/value2;
  return (-1.0)*Sqrt(1.5) * value;
}

