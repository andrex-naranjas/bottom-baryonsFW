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
`
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
  double decayWidth = flav_coup * std::pow(gamma, 2) * fi2_value * (1./(2*JA + 1)) * angular_sum_value;
  return decayWidth*GeV;
}


double EMDecayWidths::ANGULAR_SUM(double alpha_rho, double alpha_lam,
				  double alpha_mes, double k_value){
  
  WignerSymbols *m_wigner = new WignerSymbols();
  return 1.0;
}


// SPIN-FLIP INTEGRALS
double EMDecayWidths::SPINFLIP_U1_GS_GS(double k_value, double alpha_lam,  double MB, double ML){
  double value1 = (-1.0)*std::pow(k_value, 2) / 8.;
  double value2 = 3.*std::pow(MB, 2) / (std::pow(alpha_lam * (MB + 2.*ML), 2));
  double value3 = 1./std::pow(alpha_lam, 2);
  double value = std::exp(value1 * (value2 + value3));
  return value;
}
