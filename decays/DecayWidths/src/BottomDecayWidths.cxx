//BottomDecay Class includes
#ifndef BOTTOMDECAYWIDTHS_CXX
#define BOTTOMDECAYWIDTHS_CXX

#include "BottomDecayWidths.h"
#include "WignerSymbols.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
//#include <complex>


BottomDecayWidths::BottomDecayWidths()
{
}

BottomDecayWidths::~BottomDecayWidths(){}

double BottomDecayWidths::execute(double ma_val, double mb_val, double mc_val, double ga_val, double sa_val,
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
  int diagram = 0;
  
  // options according baryons and decay products
  if(baryonFlag==1){// omegas
    if(decayProd == 1)     {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.; SC=0; diagram = 1;} //Xi+K
    else if(decayProd==2)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.; SC=0; diagram = 1;} //Xi'+K
    else if(decayProd==3)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.; SC=0; diagram = 1;} //Xi*+K
    else if(decayProd==4)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.; SC=1; diagram = 1;} //Xi+K*
    else if(decayProd==5)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.; SC=1; diagram = 1;} //Xi'+K*
    else if(decayProd==6)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.; SC=1; diagram = 1;} //Xi'+K*
    else if(decayProd==7)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 2./9.; SC=0; diagram = 1;} //Omega+eta
    else if(decayProd==8)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 2./9.; SC=0; diagram = 1;} //Omega*+eta
    else if(decayProd==9)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./9.; SC=1; diagram = 1;} //Omega+phi
    else if(decayProd==10) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./9.; SC=1; diagram = 1;} //Omega*+phi
    else if(decayProd==11) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.; SC=0; diagram = 1;} //Omega+eta'
    else if(decayProd==12) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.; SC=0; diagram = 1;} //Omega*+eta'
    // SECOND DIAGRAM
    else if(decayProd==13) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 2./9.; SC=0; diagram = 2;} //cascade_8 D
    else if(decayProd==14) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=1.5; flav_coup = 2./9.; SC=0; diagram = 2;} //cascade_10 D

  }else if(baryonFlag==2 or baryonFlag==5){//cascades (6-plet, 3-plet)
    if(decayProd == 1)     {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./12.; SC=0; diagram = 1;} //lambda+K
    else if(decayProd==2)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./8.;  SC=0; diagram = 1;} //Xi+Pi
    else if(decayProd==3)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./8.;  SC=0; diagram = 1;} //Xi'(2578)+Pi comp. w/mathematica
    else if(decayProd==4)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./8.;  SC=0; diagram = 1;} //Xi*+pi
    else if(decayProd==5)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./4.;  SC=0; diagram = 1;} //Sigma+K comp. w/ mathematica
    else if(decayProd==6)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./4.;  SC=0; diagram = 1;} //Sigma*+K
    else if(decayProd==7)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./72.; SC=0; diagram = 1;} //Xi+eta    
    else if(decayProd==8)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./12.; SC=1; diagram = 1;} //Lambda+K*
    else if(decayProd==9)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./8.;  SC=1; diagram = 1;} //Xi+rho
    else if(decayProd==10) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./8.;  SC=1; diagram = 1;} //Xi'+rho
    else if(decayProd==11) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./8.;  SC=1; diagram = 1;} //Xi*+rho
    else if(decayProd==12) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./4.;  SC=1; diagram = 1;} //Sigma+K*
    else if(decayProd==13) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./4.;  SC=1; diagram = 1;} //Sigma*+K*
    else if(decayProd==14) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./72.; SC=0; diagram = 1;} //Xi'+eta
    else if(decayProd==15) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./72.; SC=0; diagram = 1;} //Xi*+eta
    else if(decayProd==16) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./9.;  SC=0; diagram = 1;} //Xi+eta'
    else if(decayProd==17) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./9.;  SC=0; diagram = 1;} //Xi'+eta'
    else if(decayProd==18) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./9.;  SC=0; diagram = 1;} //Xi*+eta'
    else if(decayProd==19) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./24.; SC=1; diagram = 1;} //Xi+omega
    else if(decayProd==20) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./24.; SC=1; diagram = 1;} //Xi'+omega
    else if(decayProd==21) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./24.; SC=1; diagram = 1;} //Xi*+omega
    else if(decayProd==22) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./12.; SC=1; diagram = 1;} //Xi+phi
    else if(decayProd==23) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./12.; SC=1; diagram = 1;} //Xi'+phi
    else if(decayProd==24) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./12.; SC=1; diagram = 1;} //Xi*+phi
    // SECOND DIAGRAM FOR CASCADES PRIME    
    else if(decayProd==25 && baryonFlag==2) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./12.; SC=0; diagram = 2;} //Sigma_8_lam+D
    else if(decayProd==26 && baryonFlag==2) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./18.; SC=0; diagram = 2;} //Xi_8_lam+Ds
    else if(decayProd==27 && baryonFlag==2) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./12.; SC=1; diagram = 2;} //Sigma_8_lam+D*
    else if(decayProd==28 && baryonFlag==2) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=1.5; flav_coup = 1./3. ; SC=0; diagram = 2;} //Sigma_10_lam+D
    // SECOND DIAGRAM FOR CASCADES-ANTITRIPLET
    else if(decayProd==25 && baryonFlag==5) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./36.; SC=0; diagram = 2;} //Lambda_8_rho+D
    else if(decayProd==26 && baryonFlag==5) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./36.; SC=1; diagram = 2;} //Lambda_8_rho+D*
    else if(decayProd==27 && baryonFlag==5) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./4.;  SC=0; diagram = 2;} //Sigma_8_rho+D
    else if(decayProd==28 && baryonFlag==5) {slf_val=0.; sb_val=0.5; lb_val=1.0; jb_val=1.5; flav_coup = 1./72.; SC=0; diagram = 2;} //Lambda*_8_rho+D
  }else if(baryonFlag==3){// sigmas
    if(decayProd == 1)     {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.;  SC=0; diagram = 1;} //Sigma+Pi
    else if(decayProd==2)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.;  SC=0; diagram = 1;} //Sigma*+pi
    else if(decayProd==3)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./2.;  SC=0; diagram = 1;} //Lambda+Pi
    else if(decayProd==4)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./18.; SC=0; diagram = 1;} //Sigma+eta
    else if(decayProd==5)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.;  SC=0; diagram = 1;} //Xi+K    
    else if(decayProd==6)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.;  SC=1; diagram = 1;} //Sigma+rho
    else if(decayProd==7)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.;  SC=1; diagram = 1;} //Sigma*+rho
    else if(decayProd==8)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./3.;  SC=1; diagram = 1;} //Lambda+rho
    else if(decayProd==9)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./18.; SC=0; diagram = 1;} //Sigma*+eta
    else if(decayProd==10) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./9.;  SC=0; diagram = 1;} //Sigma+eta'
    else if(decayProd==11) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./9.;  SC=0; diagram = 1;} //Sigma*+eta'
    else if(decayProd==12) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.;  SC=0; diagram = 1;} //Xi'+K
    else if(decayProd==13) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.;  SC=0; diagram = 1;} //Xi*+K
    else if(decayProd==14) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.;  SC=1; diagram = 1;} //Xi+K*
    else if(decayProd==15) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.;  SC=1; diagram = 1;} //Xi'+K*
    else if(decayProd==16) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.;  SC=1; diagram = 1;} //Xi*+K*
    else if(decayProd==17) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.;  SC=1; diagram = 1;} //Sigma+omega
    else if(decayProd==18) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.;  SC=1; diagram = 1;} //Sigma+omega
    // SECOND DIAGRAM
    else if(decayProd==19) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./9.;  SC=0; diagram = 2;} //N+D
    else if(decayProd==20) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./9.;  SC=0; diagram = 2;} //Sigma_8+Ds
    else if(decayProd==21) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./9.;  SC=1; diagram = 2;} //N+D*
    else if(decayProd==22) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=1.5; flav_coup = 4./9.;  SC=0; diagram = 2;} //Delta+D    
    else if(decayProd==23) {slf_val=1.; sb_val=0.5; lb_val=1.0; jb_val=1.5; flav_coup = 1./18.; SC=0; diagram = 2;} //N*(1520)+D
    else if(decayProd==24) {slf_val=1.; sb_val=0.5; lb_val=1.0; jb_val=0.5; flav_coup = 1./18.; SC=0; diagram = 2;} //N*(1535)+D    
    else if(decayProd==25) {slf_val=1.; sb_val=0.5; lb_val=2.0; jb_val=2.5; flav_coup = 1./9.;  SC=0; diagram = 2;} //N*(1680)+D
    else if(decayProd==26) {slf_val=1.; sb_val=0.5; lb_val=2.0; jb_val=1.5; flav_coup = 1./9.;  SC=0; diagram = 2;} //N*(1720)+D
  }else if(baryonFlag==4){// lambdas
    if(decayProd == 1)     {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./2.; SC=0; diagram = 1;} //Sigma+Pi
    else if(decayProd==2)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./2.; SC=0; diagram = 1;} //Sigma*+pi
    else if(decayProd==3)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./18.;SC=0; diagram = 1;} //Lambda+eta    
    else if(decayProd==4)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./2.; SC=1; diagram = 1;} //Sigma+rho
    else if(decayProd==5)  {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./2.; SC=1; diagram = 1;} //Sigma*+rho
    else if(decayProd==6)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./9.; SC=0; diagram = 1;} //Lambda+eta'
    else if(decayProd==7)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.; SC=1; diagram = 1;} //Lambda+omega
    else if(decayProd==8)  {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.; SC=0; diagram = 1;} //Xi+K
    else if(decayProd==9)  {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.; SC=0; diagram = 1;} //Xi'+K
    else if(decayProd==10) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.; SC=0; diagram = 1;} //Xi*+K
    else if(decayProd==11) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.; SC=1; diagram = 1;} //Xi+K*
    else if(decayProd==12) {slf_val=1.; sb_val=0.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.; SC=1; diagram = 1;} //Xi'+K*
    else if(decayProd==13) {slf_val=1.; sb_val=1.5; lb_val=0.0; jb_val=0.0; flav_coup = 1./6.; SC=1; diagram = 1;} //Xi*+K*
    // SECOND DIAGRAM    
    else if(decayProd==14) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./3.; SC=0; diagram = 2;} //N+D
    //else if(decayProd==15) {slf_val=0.; sb_val=0.5; lb_val=0.0; jb_val=0.5; flav_coup = 1./3.; SC=1; diagram = 2;} //N+D*
  }

  alpha_rho = ar_val;
  alpha_lam = al_val;
  alpha_mes = ALPHA_MES(diagram);
  if(modeExcitation==5){l1_val=1.; l2_val=1.;}//for mixed states
  else{l1_val=0.; l2_val=0.;}
  
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
  
  double sum_value  = 0.;
  if(diagram==1) sum_value  = ANGULAR_SUM(alpha_rho, alpha_lam, alpha_mes, k_value);
  else if(diagram==2) sum_value = ANGULAR_SUM_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);

  double fi2_value  = FI2(EB_value, EWCC_value, MA, k_value);
  double decayWidth = DecayWidth(flav_coup, gamma, fi2_value, sum_value);

  return decayWidth;
}

double BottomDecayWidths::DecayWidth(double flav_coup, double gamma, double fi2_value, double angular_sum_value){
  double GeV = 1000.;
  double decayWidth = flav_coup * std::pow(gamma, 2) * fi2_value * (1./(2*JA + 1)) * angular_sum_value;
  return decayWidth*GeV;
}

double BottomDecayWidths::ALPHA_MES(int diagram){
  if(diagram == 1) return 0.46;
  else return 0.6;
}

double BottomDecayWidths::ANGULAR_SUM(double alpha_rho, double alpha_lam,
 				     double alpha_mes, double k_value){
  
  WignerSymbols *m_wigner = new WignerSymbols();
  double outerSum = 0;
  double finalIntegral1=0., finalIntegral2=0.; //, finalIntegral3=0.;
  double dummy=0;
  double delta1=0, delta2=0, delta_mix=1.;;
  double cmix_d1_a=0,cmix_d1_b=0;
  double cmix_d2_a=0,cmix_d2_b=0,cmix_d2_c=0;

  if(modeExcitation == 0){//GROUND
    finalIntegral1=I01B0TOT_GROUND(alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 1 && LA==1){//P-WAVE, LAMBDA
    finalIntegral1=I01B0TOT_PWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT_PWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 2 && LA==1){//P-WAVE, RHO
    finalIntegral1=I01B0TOT_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 1 && LA==2){//D-WAVE, LAMBDA
    finalIntegral1=I01B0TOT_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 2 && LA==2){//D-WAVE, RHO
    finalIntegral1=I01B0TOT_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 3){//RADIAL, LAMBDA
    finalIntegral1=I01B0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 4){//RADIAL, RHO
    finalIntegral1=I01B0TOT_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 5){//MIXED
    finalIntegral1=I01B0TOT_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value);
    cmix_d1_a = std::pow(6,0.5)/(2*std::pow(ARO0_MIXED(alpha_lam, alpha_mes),2));
    cmix_d1_b = 3.*(ARO5_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value) * ARO6_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value));
    cmix_d2_a = ARO5_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value) * ARO4_MIXED(alpha_lam, alpha_mes);
    cmix_d2_b = ARO6_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value) * ARO1_MIXED(alpha_lam, alpha_mes);
    cmix_d2_c = (BRO1_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_MIXED(alpha_rho, alpha_lam, alpha_mes), 2)) *
      ARO1_MIXED(alpha_lam, alpha_mes)* ARO4_MIXED(alpha_lam, alpha_mes);
  }
  
  for(int iMJA = 0; iMJA<(int)mJA.size(); iMJA++){
    
    double innerSum = 0;
    for(int iMLA = 0; iMLA<(int)mLA.size(); iMLA++)
      for(int iML1 = 0; iML1<(int)mL1.size(); iML1++)
	for(int iML2 = 0; iML2<(int)mL2.size(); iML2++)
	  for(int iMSA = 0; iMSA<(int)mSA.size(); iMSA++)
	    for(int iM = 0; iM<(int)m.size(); iM++)
	      for(int iM24 = 0; iM24<(int)m24.size(); iM24++)
		for(int iM1 = 0; iM1<(int)m1.size(); iM1++)
		  for(int iMSB = 0; iMSB<(int)mSB.size(); iMSB++)
		    for(int iM3 = 0; iM3<(int)m3.size(); iM3++)
		      for(int iM5 = 0; iM5<(int)m5.size(); iM5++)
			for(int iMSC = 0; iMSC<(int)mSC.size(); iMSC++)
			  for(int iM23 = 0; iM23<(int)m23.size(); iM23++)
			    for(int iM4 = 0; iM4<(int)m4.size(); iM4++)
			      for(int iM2 = 0; iM2<(int)m2.size(); iM2++){
				delta1=0.; delta2=0.; delta_mix=1.0;
				if(modeExcitation == 0){//GROUND
				  delta1 = KroneckerDelta(0, m.at(iM));
				}else if(modeExcitation == 1 && LA == 1){//P-WAVE, LAMBDA
				  delta1 = KroneckerDelta(m.at(iM), 0)*KroneckerDelta(mLA.at(iMLA), 0);
				  delta2 = KroneckerDelta(m.at(iM), mLA.at(iMLA));
				}else if(modeExcitation == 2 && LA ==1){//P-WAVE, RHO
				  delta1 = KroneckerDelta(m.at(iM), 0)*KroneckerDelta(mLA.at(iMLA), 0);
				  delta2 = KroneckerDelta(m.at(iM), mLA.at(iMLA));
				}else if((modeExcitation == 1 || modeExcitation == 2)  && LA == 2){//D-WAVE,same for LAMBDA & RHO
				  delta1 = KroneckerDelta(m.at(iM), 0) * KroneckerDelta(mLA.at(iMLA), 0);
				  delta2 = (m_wigner->wigner3j(1, 1, LA, m.at(iM), 0, (-1.0)*m.at(iM))) * std::pow(-1.0, (-1.0)*m.at(iM)) *
				    std::pow((2*LA+1),0.5) * KroneckerDelta(m.at(iM), mLA.at(iMLA));
				}else if(modeExcitation == 3 || modeExcitation == 4){//RADIAL,same for LAMBDA & RHO
				  delta1 = KroneckerDelta(m.at(iM), 0);
				  delta2 = KroneckerDelta(m.at(iM), 0);}
				else if(modeExcitation == 5){// MIXED states
				  delta1 = (cmix_d1_a*KroneckerDelta(m.at(iM), mL2.at(iML2)) + cmix_d1_b*KroneckerDelta(m.at(iM),0)*KroneckerDelta(0, mL2.at(iML2)))*KroneckerDelta(mL1.at(iML1), 0);
				  delta2 = cmix_d2_a*KroneckerDelta(m.at(iM), 0)*KroneckerDelta(mL1.at(iML1), mL2.at(iML2)) +
				    cmix_d2_b*KroneckerDelta(0, mL2.at(iML2))*KroneckerDelta(mL1.at(iML1), m.at(iM)) + cmix_d2_c*KroneckerDelta(mL1.at(iML1), 0)*KroneckerDelta(m.at(iM), mL2.at(iML2));
				  delta_mix = m_wigner->wigner3j(1, 1, LA, mL1.at(iML1), mL2.at(iML2), (-1.0)*mLA.at(iMLA)) * std::pow(-1.0, (-1.0)*mLA.at(iMLA) ) * std::pow(2*LA+1,0.5);
				}
				
				if(delta1!=0 || delta2!=0){
				  dummy = (finalIntegral1*delta1 + finalIntegral2*delta2)*delta_mix*
				    m_wigner->wigner3j(LA, SA, JA, mLA.at(iMLA), mSA.at(iMSA), (-1.0)*mJA.at(iMJA))*
				    m_wigner->wigner3j(1, 1, 0, m.at(iM), (-1.0)*m.at(iM), 0)*
				    m_wigner->wigner3j(slightf, 0.5, SB, m24.at(iM24), m1.at(iM1), (-1.0)*mSB.at(iMSB))*
				    m_wigner->wigner3j(0.5, 0.5, SC, m3.at(iM3), m5.at(iM5), (-1.0)*mSC.at(iMSC))*
				    m_wigner->wigner3j(slight, 0.5, SA, m23.at(iM23), m1.at(iM1), (-1.0)*mSA.at(iMSA))*
				    m_wigner->wigner3j(0.5, 0.5, 1, m4.at(iM4), m5.at(iM5), m.at(iM))*
				    m_wigner->wigner3j(0.5, 0.5, slightf, m2.at(iM2), m4.at(iM4), (-1.0)*m24.at(iM24))*
				    m_wigner->wigner3j(0.5,0.5,slight, m2.at(iM2), m3.at(iM3), (-1.0)*m23.at(iM23))*
				    std::pow(3 * (2*JA+1) * (2*slight+1) * (2*slightf+1) * (2*SA+1) * (2*SB+1) * (2*SC+1),0.5)*
				    std::pow(-1.0,SA-LA-mJA.at(iMJA))*
				    std::pow(-1.0,1+m.at(iM)-mSA.at(iMSA)-mSB.at(iMSB)-slight-slightf-m23.at(iM23)-m24.at(iM24)-mSC.at(iMSC));
				}else dummy = 0;
				
				innerSum+=dummy;
				dummy = 0;
			      }
    outerSum += std::pow(innerSum,2);        
  }
  delete m_wigner;
  return outerSum;
}

double BottomDecayWidths::ANGULAR_SUM_SECOND(double alpha_rho, double alpha_lam,
					    double alpha_mes, double k_value){
  
  if(modeExcitation == 0 ||  // GROUND, RADIAL RHO do not contribute
     modeExcitation == 4 ) return 0.;
  
  WignerSymbols *m_wigner = new WignerSymbols();
  double outerSum = 0;
  double finalIntegral1=0., finalIntegral2=0.;
  double dummy=0;
  double delta1=0, delta2=0, factor=0;

  if(modeExcitation == 1 && LA==1){//P-WAVE, LAMBDA
    if(baryonFlag==1 || baryonFlag==2 || baryonFlag==4 || baryonFlag==5) return 0.;
    if(baryonFlag==3){
      if(decayProd>=19 && decayProd<=22){ // sigmas decay only through
	finalIntegral1=I01B0TOT_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
	finalIntegral2=I02B0TOT_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
    }else return 0.;}
    
  }else if(modeExcitation == 1 && LA==2){//D-WAVE, LAMBDA    
    if(baryonFlag==2){
      if(decayProd>=25 && decayProd<=28){
	finalIntegral1=I01B0TOT_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
	finalIntegral2=I02B0TOT_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
      }
    }else if(baryonFlag==5){
      if(decayProd>=25 && decayProd<=28){
	finalIntegral1=I01B0TOT_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
	finalIntegral2=I02B0TOT_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
      }
    }else{
      finalIntegral1=I01B0TOT_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
      finalIntegral2=I02B0TOT_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
    }
    
  }else if(modeExcitation == 2 && LA==2){//D-WAVE, RHO
    if(baryonFlag==1 || baryonFlag==2 || baryonFlag==4 || baryonFlag==5) return 0.;
    if(baryonFlag==3){
      if(decayProd>=25 && decayProd<=26){ // sigmas decay only through 
      finalIntegral1=I01B0TOT_DWAVE_RHO_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
      finalIntegral2=0;
    }else return 0;}
        
  }else if(modeExcitation == 3){//RADIAL, LAMBDA
    if(baryonFlag==2){
      if(decayProd>=25 && decayProd<=28){
	finalIntegral1=I01B0TOT_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
	finalIntegral2=I02B0TOT_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
      }
    }else if(baryonFlag==5){
      if(decayProd>=25 && decayProd<=28){
	finalIntegral1=I01B0TOT_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
	finalIntegral2=I02B0TOT_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
      }
    }else{
      finalIntegral1=I01B0TOT_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
      finalIntegral2=I02B0TOT_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
    }

  }else if(modeExcitation == 5){//MIXED
    
    if(baryonFlag==1 || baryonFlag==2 || baryonFlag==4) return 0.;

    if(baryonFlag==3){
      if(decayProd>=23 && decayProd<=24){ // sigmas decay only through
	finalIntegral1=I01B0TOT_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
	finalIntegral2=I02B0TOT_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
      }else return 0;}
    
    if(baryonFlag==5){
      if(decayProd==32){
	finalIntegral1=I01B0TOT_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
	finalIntegral2=I02B0TOT_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
    }else return 0;}
    
  }
  
  for(int iMJA = 0; iMJA<(int)mJA.size(); iMJA++){
    for(int iMJB = 0; iMJB<(int)mJB.size(); iMJB++){
    
      double innerSum = 0;
      for(int iMLA = 0; iMLA<(int)mLA.size(); iMLA++)
	for(int iMSA = 0; iMSA<(int)mSA.size(); iMSA++)
	  for(int iM = 0; iM<(int)m.size(); iM++)
	    for(int iMLB = 0; iMLB<(int)mLB.size(); iMLB++)
	      for(int iMSB = 0; iMSB<(int)mSB.size(); iMSB++)
		for(int iMSC = 0; iMSC<(int)mSC.size(); iMSC++)
		  for(int iM3 = 0; iM3<(int)m3.size(); iM3++)
		    for(int iM4 = 0; iM4<(int)m4.size(); iM4++)
		      for(int iM5 = 0; iM5<(int)m5.size(); iM5++)
			for(int iM12 = 0; iM12<(int)m12.size(); iM12++){
			  delta1=0.; delta2=0.;				
			  if(modeExcitation == 1 && LA == 1){//P-WAVE, LAMBDA
			    delta1 = KroneckerDelta(m.at(iM), 0)*KroneckerDelta(mLA.at(iMLA), 0);
			    delta2 = KroneckerDelta(m.at(iM), mLA.at(iMLA));
			  }else if((modeExcitation == 1 || modeExcitation == 2)  && LA == 2){//D-WAVE,same for LAMBDA & RHO
			    delta1 = KroneckerDelta(m.at(iM), 0) * KroneckerDelta(mLA.at(iMLA), mLB.at(iMLB));
			    delta2 = (m_wigner->wigner3j(1, 1, LA, m.at(iM), 0, (-1.0)*m.at(iM))) * std::pow(-1.0, (-1.0)*m.at(iM)) *
			      std::pow((2*LA+1),0.5) * KroneckerDelta(m.at(iM), mLA.at(iMLA));
			  }else if(modeExcitation == 2 && LA==2){//D-WAVE, RHO
			    delta1 = KroneckerDelta(m.at(iM), 0) * KroneckerDelta(mLA.at(iMLA), mLB.at(iMLB));
			    delta2=0;
			  }else if(modeExcitation == 3 || modeExcitation == 4){//RADIAL,same for LAMBDA & RHO
			    delta1 = KroneckerDelta(m.at(iM), 0);
			    delta2 = KroneckerDelta(m.at(iM), 0);
			  }else if(modeExcitation == 5 ){//MIXED stated
			    factor = m_wigner->wigner3j(1, 1, LA, m.at(iM), mLB.at(iMLB), (-1.0)*mLA.at(iMLA) ) * std::pow(-1, (-1.0)*mLA.at(iMLA)) * std::pow(2*LA + 1, 0.5);
			    delta1 = KroneckerDelta(m.at(iM), 0) * factor;
			    delta2 = factor;}
		    
			  if(delta1!=0 || delta2!=0){
			    dummy = (finalIntegral1*delta1 + finalIntegral2*delta2)*
			      m_wigner->wigner3j(LA, SA, JA, mLA.at(iMLA), mSA.at(iMSA), (-1.0)*mJA.at(iMJA))*
			      m_wigner->wigner3j(LB, SB, JB, mLB.at(iMLB), mSB.at(iMSB), (-1.0)*mJB.at(iMJB))*
			      m_wigner->wigner3j(1, 1, 0, m.at(iM), (-1.0)*m.at(iM), 0)*
			      m_wigner->wigner3j(slight, 0.5, SB, m12.at(iM12), m4.at(iM4), (-1.0)*mSB.at(iMSB))*
			      m_wigner->wigner3j(0.5, 0.5, SC, m3.at(iM3), m5.at(iM5), (-1.0)*mSC.at(iMSC))*
			      m_wigner->wigner3j(slight, 0.5, SA, m12.at(iM12), m3.at(iM3), (-1.0)*mSA.at(iMSA))*
			      m_wigner->wigner3j(0.5, 0.5, 1, m4.at(iM4), m5.at(iM5), m.at(iM))*				   
			      std::pow(3 * (2*JA+1) * (2*SA+1) * (2*JB+1)* (2*SB+1) * (2*SC+1), 0.5)*
			      std::pow(-1.0,SA-LA-mJA.at(iMJA))*
			      std::pow(-1.0,SB-LB-mJB.at(iMJB))*
			      std::pow(-1.0,1+m.at(iM)-mSA.at(iMSA)-mSB.at(iMSB)-mSC.at(iMSC));
			  }else dummy = 0;
		      
			  innerSum+=dummy;
			  dummy = 0;
			}
      outerSum += std::pow(innerSum,2); 
    }
  }
  delete m_wigner;
  return outerSum;
}

std::vector<double> BottomDecayWidths::getMomentumProjections(double j_angular){
  //gets the m projections "m_projection" for a given angular momentum "j_angular"
  std::vector<double> angularProjections; angularProjections.clear();
  if(j_angular==0.){ angularProjections.push_back(0); return angularProjections;}
  
  double m_projection = (-1.0)*j_angular;
  do{
    angularProjections.push_back(m_projection);
    m_projection++;
  }while(m_projection<=j_angular);
  
  return angularProjections;
}

int BottomDecayWidths::KroneckerDelta(double i, double j){
  if(i==j) return 1;
  else return 0;
}

double BottomDecayWidths::EWCC(double MA, double MB, double MC){
  double value = (0.5*(std::pow(MA,2) + std::pow(MC,2) - std::pow(MB,2)) )/ MA;
  return value;
}

double BottomDecayWidths::EB(double MA, double MB, double MC){
  double value = ( 0.5*(std::pow(MA,2) - std::pow(MC,2) + std::pow(MB,2)) )/ MA;
  return value;
}

double BottomDecayWidths::K(double EB, double MB){
  double value = std::sqrt( std::pow(EB,2) - std::pow(MB,2) );
  return value;
}

double BottomDecayWidths::FI2(double EB, double EWCC, double MA, double k_value){
  double value = (2*(pi_val)*k_value*(EB * EWCC)) / MA ; 
  return value;
}

double BottomDecayWidths::CBARIN(double alpha_rho, double alpha_lam){
  double value1 = std::pow(3,0.75);
  double value2 = std::pow(1.0/( pi_val * std::pow(alpha_rho, 2) ), 0.75);
  double value3 = std::pow(1.0/( pi_val * std::pow(alpha_lam, 2) ), 0.75);
  return value1 * value2 * value3 ;
}

double BottomDecayWidths::CMESON(double alpha_mes){
  double value = std::pow(1.0/( pi_val * std::pow(alpha_mes, 2) ), 0.75);
  return value;
}

//GROUND STATES
//CBARIN_GROUND --> same as P-wave (CBARIN)
//CMESON_GROUND --> same as P-wave (CMESON)
double BottomDecayWidths::C0_GROUND(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3*std::pow(3,0.5));
  double value2 = CBARIN(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double BottomDecayWidths::ARO0_GROUND(double alpha_lam, double alpha_mes){
  double value1 = 1./std::pow(alpha_lam,2);
  double value2 = 1./(12.*std::pow(alpha_mes,2));
  return std::pow(value1+value2,0.5);
}

double BottomDecayWidths::F0_GROUND(double alpha_rho, double alpha_lam, double k_value){
  double value1 = (-1./3.)/(alpha_rho*alpha_lam);
  return k_value*k_value*value1;
}

double BottomDecayWidths::F01_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1./(2*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(6,0.5)*alpha_mes*alpha_mes);  
  double value1pp = value1p*value1p;
  double value1ppp= 1./std::pow(ARO0_GROUND(alpha_lam, alpha_mes),2);
  double value1 = value1pp*value1ppp;
  double value2 = 1./(8*std::pow(alpha_mes,2));  
  return k_value*k_value*(value1-value2);
}

double BottomDecayWidths::BRO2_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1./(2*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(6,0.5)*alpha_mes*alpha_mes);  
  double value1pp= 1./std::pow(ARO0_GROUND(alpha_lam, alpha_mes),2);
  double value1 = (std::pow(2./3.,0.5))*value1p*value1pp;
  double value2 = 2.;  
  return k_value*(value1-value2);  
}

double BottomDecayWidths::ARO1_GROUND(double alpha_lam, double alpha_mes){
  double value1 = std::pow(2,0.5);
  double value2 = (std::pow(2./3.,0.5))/(4 * std::pow(3,0.5) * std::pow(ARO0_GROUND(alpha_lam, alpha_mes),2) * std::pow(alpha_mes,2) );
  return -value1+value2;
}

double BottomDecayWidths::ARO2_GROUND(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4./std::pow(alpha_rho,2);
  double value1pp = 1./std::pow(alpha_mes,2);
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1./(48*std::pow(alpha_mes,4)*std::pow(ARO0_GROUND(alpha_lam, alpha_mes),2));
  return std::pow(value1-value2, 0.5);
}

double BottomDecayWidths::BRO1_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1  = k_value / ( 4*std::pow(2,0.5) * std::pow(alpha_mes,2) );
  double value2  = k_value / ( 2*std::pow(2,0.5) * alpha_rho * alpha_lam );
  double value3p = k_value / ( 4*std::pow(3,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0_GROUND(alpha_lam, alpha_mes), 2) );
  double value3pp= 1.0/( 2*std::pow(6,0.5)*alpha_lam*alpha_rho ) + 1.0/(4*std::pow(6,0.5)*alpha_mes*alpha_mes);
  double value3  = value3p*value3pp;
  return (-1.0)*(value1 + value2 + value3);
}

double BottomDecayWidths::F0TOT_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1 = F01_GROUND(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = BRO1_GROUND(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO2_GROUND(alpha_rho, alpha_lam, alpha_mes) ;
  value2 = std::pow(value2,2);  
  double value3 = F0_GROUND(alpha_rho, alpha_lam, k_value);  
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1   = C0_GROUND(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_GROUND(alpha_lam, alpha_mes), 3) );
  double value2p  = BRO1_GROUND(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_GROUND(alpha_rho, alpha_lam, alpha_mes), 2);
  double value2   = value2p*ARO1_GROUND(alpha_lam, alpha_mes) + BRO2_GROUND(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3   = std::pow(0.75,0.5) / (1.0 * std::pow( ARO2_GROUND(alpha_rho, alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*std::pow(pi_val,2.5);
}

double BottomDecayWidths::I01B0TOT_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_GROUND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_GROUND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

//P-WAVE-LAMBDA STATES
double BottomDecayWidths::CBARFIN_PWAVE(double alpha_rho, double alpha_lam){//backup
  double value1 = std::pow(3,0.75) * std::pow( 8.0/(3.0*std::pow(pi_val, 0.5)), 0.5);
  double value2 = std::pow(1.0/( std::pow(alpha_lam, 2) ), 1.25);
  double value3 = std::pow(1.0/( pi_val * std::pow(alpha_rho, 2) ), 0.75);
  return value1 * value2 * value3;// * mycomplex;//define complex if needed
}

double BottomDecayWidths::C0_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3*std::pow(3,0.5));
  double value2 = CBARFIN_PWAVE(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double BottomDecayWidths::ARO0_PWAVE(double alpha_rho, double alpha_mes){
  double value1 = (1.0/1.0) * (1./std::pow(alpha_rho,2));
  double value2 = (1.0/4.0) * (1./std::pow(alpha_mes,2));
  return std::pow(value1+value2, 0.5);
}

double BottomDecayWidths::F0_PWAVE(double alpha_rho, double alpha_lam, double k_value){
  double value1 = std::pow(k_value, 2);
  double value2 =1./(alpha_rho*alpha_lam);
  return (-1./3.)*(value1*value2);  
}

double BottomDecayWidths::F01_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1.0/(2.0*std::pow(2,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(2,0.5)*alpha_mes*alpha_mes);
  double value1pp = value1p*value1p;
  double value1ppp= 1./(std::pow(ARO0_PWAVE(alpha_rho, alpha_mes),2));  
  double value1   = value1pp*value1ppp;
  double value2   = 1.0/(8.0*std::pow(alpha_mes,2));  
  return std::pow(k_value,2)*(value1-value2);
}

double BottomDecayWidths::BRO2_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = 1.0/(2.0*alpha_lam*alpha_rho) + 1.0/(4.0*alpha_mes*alpha_mes);
  double value2 = 1.0/(std::pow(ARO0_PWAVE(alpha_rho, alpha_mes),2));
  return k_value*(value1*value2 - 2.0);
}

double BottomDecayWidths::ARO1_PWAVE(double alpha_rho, double alpha_mes){
  double value1 = std::pow(2./3.,0.5);
  double value2 = 1./(2.0 * std::pow(6,0.5) * std::pow(ARO0_PWAVE(alpha_rho, alpha_mes),2) * std::pow(alpha_mes,2) );
  return value1 - value2;
}

double BottomDecayWidths::ARO2_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4.0/(std::pow(alpha_lam,2));
  double value1pp = 1.0/(3.0*std::pow(alpha_mes,2));
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1.0/(48.0 * std::pow(alpha_mes,4) * std::pow(ARO0_PWAVE(alpha_rho, alpha_mes),2) );
  return std::pow(value1-value2, 0.5);
}

double BottomDecayWidths::BRO1_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1  = k_value / ( 4*std::pow(6,0.5) * std::pow(alpha_mes,2) );
  double value2  = k_value / ( 2*std::pow(6,0.5) * alpha_rho * alpha_lam );
  double value3p = k_value / ( 4*std::pow(3,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0_PWAVE(alpha_rho, alpha_mes), 2) );
  double value3pp= 1.0/( 2*std::pow(2,0.5)*alpha_lam*alpha_rho ) + 1.0/(4*std::pow(2,0.5)*alpha_mes*alpha_mes);
  double value3  = value3p*value3pp;
  return value1 + value2 - value3;
}

double BottomDecayWidths::F0TOT_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_PWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = BRO1_PWAVE(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO2_PWAVE(alpha_rho, alpha_lam, alpha_mes) ;
  value2 = std::pow(value2,2);
  double value3 = F0_PWAVE(alpha_rho, alpha_lam, k_value);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_PWAVE(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_PWAVE(alpha_rho, alpha_mes), 3) );
  double value1pp = BRO1_PWAVE(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_PWAVE(alpha_rho, alpha_lam, alpha_mes), 2);
  double value1   = value1p*value1pp;
  double value2   = value1pp*ARO1_PWAVE(alpha_rho, alpha_mes) + BRO2_PWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3   = 3.0 / (4.0 * std::pow( ARO2_PWAVE(alpha_rho, alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*pi_val*pi_val;
}

double BottomDecayWidths::I02B0_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = C0_PWAVE(alpha_rho, alpha_lam, alpha_mes) / std::pow(ARO0_PWAVE(alpha_rho, alpha_mes), 3) ;
  double value2 = (3. * ARO1_PWAVE(alpha_rho, alpha_mes)) / (8.* std::pow( ARO2_PWAVE(alpha_rho, alpha_lam, alpha_mes), 5));
  return value1*value2*pi_val*pi_val;
}

double BottomDecayWidths::I01B0TOT_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_PWAVE(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_PWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_PWAVE(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_PWAVE(alpha_rho, alpha_lam, alpha_mes);
  return value1*value2;
}


// P-WAVE-RHO STATES
double BottomDecayWidths::CBARFIN_PWAVE_RHO(double alpha_rho, double alpha_lam){
  double value1 = std::pow(3,0.75) * std::pow( 8.0/(3.0*std::pow(pi_val, 0.5)), 0.5);
  double value2 = std::pow(1.0/(pi_val* std::pow(alpha_lam,2) ), 0.75);
  double value3 = std::pow(1.0/( std::pow(alpha_rho, 2) ), 1.25);
  return value1 * value2 * value3;// * mycomplex;//define complex if needed
}

double BottomDecayWidths::C0_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3*std::pow(3,0.5));
  double value2 = CBARFIN_PWAVE_RHO(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double BottomDecayWidths::ARO0_PWAVE_RHO(double alpha_lam, double alpha_mes){
  double value1 = 1./std::pow(alpha_lam,2);
  double value2 = 1./(12*std::pow(alpha_mes, 2));  
  return std::pow(value1+value2,0.5);
}

double BottomDecayWidths::F01_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1./(2*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(6,0.5)*alpha_mes*alpha_mes);  
  double value1pp = value1p*value1p;
  double value1ppp= 1./std::pow(ARO0_PWAVE_RHO(alpha_lam, alpha_mes),2);
  double value1 = value1pp*value1ppp;
  double value2 = 1./(8*std::pow(alpha_mes,2));  
  return k_value*k_value*(value1-value2);
}

double BottomDecayWidths::BRO1_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1./(2*alpha_lam*alpha_rho) + 1./(4*alpha_mes*alpha_mes);
  double value1pp = 1./(3*std::pow(ARO0_PWAVE_RHO(alpha_lam, alpha_mes), 2));
  double value1 = value1p * value1pp;
  double value2 = 2.;
  return k_value*(value1-value2);
}

double BottomDecayWidths::ARO1_PWAVE_RHO(double alpha_lam, double alpha_mes){
  double value1 = std::pow(2,0.5);
  double value2 = 1./(6 * std::pow(2,0.5) * std::pow(ARO0_PWAVE_RHO(alpha_lam, alpha_mes),2) * std::pow(alpha_mes,2) );
  return -value1+value2;
}

double BottomDecayWidths::ARO2_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4./std::pow(alpha_rho,2);
  double value1pp = 1./std::pow(alpha_mes,2);
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1./(48*std::pow(alpha_mes,4)*std::pow(ARO0_PWAVE_RHO(alpha_lam, alpha_mes),2));
  return std::pow(value1-value2, 0.5);
}

double BottomDecayWidths::BRO2_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1    = k_value/(4*std::pow(2,0.5)*std::pow(alpha_mes,2));
  double value2    = k_value/(2*std::pow(2,0.5)*alpha_rho*alpha_lam);
  double value3p   = k_value/(4*std::pow(3,0.5)*std::pow(alpha_mes,2));
  double value3ppp =      1./(2*std::pow(6,0.5)*alpha_rho*alpha_lam) + 1./(4*std::pow(6,0.5)*std::pow(alpha_mes,2));
  double value3pp  = 1./std::pow(ARO0_PWAVE_RHO(alpha_lam, alpha_mes),2);
  double value3 = value3p * value3pp * value3ppp;
  return -value1 - value2 + value3;
}

double BottomDecayWidths::F0TOT_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = std::pow(BRO2_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value)/ARO2_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 2);
  double value3 = (1./3.)*(std::pow(k_value,2)*(1./(alpha_rho*alpha_lam)));
  return value1 + value2 - value3;
}

double BottomDecayWidths::I01B0_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes) / std::pow(ARO0_PWAVE_RHO(alpha_lam, alpha_mes), 3) ;
  double value1pp = BRO2_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 2);
  double value1 = value1p*value1pp;
  double value2 = value1pp*ARO1_PWAVE_RHO(alpha_lam, alpha_mes) + BRO1_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3 = 3./(4.* std::pow( ARO2_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 3));
  return value1*value2*value3*pi_val*pi_val;
}

double BottomDecayWidths::I02B0_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = C0_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes) / std::pow(ARO0_PWAVE_RHO(alpha_lam, alpha_mes), 3);
  double value2 = (3. * ARO1_PWAVE_RHO(alpha_lam, alpha_mes) ) / (8.* std::pow(  ARO2_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 5) );
  return pi_val*pi_val*value1*value2;
}

double BottomDecayWidths::I01B0TOT_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_PWAVE_RHO(alpha_rho, alpha_lam, alpha_mes);
  return value1*value2;
}


// D-wave excitations
// lambda wave functions
double BottomDecayWidths::CBARFIN_DWAVE(double alpha_rho, double alpha_lam){
  double value1 = (std::pow(3,0.75) * std::pow(16,0.5))/(std::pow(15,0.5) * (std::pow(pi_val,0.25)) );
  double value2 = std::pow(( 1.0/alpha_lam ), 3.5);
  double value3 = std::pow(1.0/( pi_val * std::pow(alpha_rho, 2) ), 0.75);
  return value1 * value2 * value3;// * mycomplex;//define complex if needed
}

//CBARIN_DWAVE --> same as P-wave (CBARIN)
//CMESON_DWAVE --> same as P-wave (CMESON)

double BottomDecayWidths::C0_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3.0*std::pow(3,0.5));
  double value2 = CBARFIN_DWAVE(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double BottomDecayWidths::ARO0_DWAVE(double alpha_rho, double alpha_mes){
  double value1 = 1.0/std::pow(alpha_rho,2);  
  double value2 = 1.0/(4.0*std::pow(alpha_mes,2));
  return std::pow(value1+value2,0.5);
}

double BottomDecayWidths::F0_DWAVE(double k_value, double alpha_rho, double alpha_lam){
  double value1 = std::pow(k_value, 2) /(alpha_rho * alpha_lam) ;
  return (-1.0/3.0) * value1;
}

double BottomDecayWidths::F01_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1p  = 1.0/(2.0*std::pow(2,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(2,0.5)*alpha_mes*alpha_mes);
  double value1pp = value1p*value1p;
  double value1ppp= 1./(std::pow(ARO0_DWAVE(alpha_rho, alpha_mes),2));  
  double value1   = value1pp*value1ppp;
  double value2   = 1.0/(8.0*std::pow(alpha_mes,2));  
  return std::pow(k_value,2)*(value1-value2);
}

double BottomDecayWidths::BRO2_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = 1.0/(2.0*alpha_lam*alpha_rho) + 1.0/(4.0*alpha_mes*alpha_mes);
  double value2 = 1.0/(std::pow(ARO0_DWAVE(alpha_rho, alpha_mes),2));
  return k_value*(value1*value2 - 2.0);
}

double BottomDecayWidths::ARO1_DWAVE(double alpha_rho, double alpha_mes){
  double value1 = std::pow(2./3.,0.5);
  double value2 = 1./(2.0 * std::pow(6,0.5) * std::pow(ARO0_DWAVE(alpha_rho, alpha_mes),2) * std::pow(alpha_mes,2) );
  return value1 - value2;
}

double BottomDecayWidths::ARO2_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4.0/(std::pow(alpha_lam,2));
  double value1pp = 1.0/(3.0*std::pow(alpha_mes,2));
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1.0/(48.0 * std::pow(alpha_mes,4) * std::pow(ARO0_DWAVE(alpha_rho, alpha_mes),2) );
  return std::pow(value1-value2, 0.5);
}

double BottomDecayWidths::BRO1_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1  = k_value / ( 4*std::pow(6,0.5) * std::pow(alpha_mes,2) );
  double value2  = k_value / ( 2*std::pow(6,0.5) * alpha_rho * alpha_lam );
  double value3p = k_value / ( 4*std::pow(3,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0_DWAVE(alpha_rho, alpha_mes), 2) );
  double value3pp= 1.0/( 2*std::pow(2,0.5)*alpha_lam*alpha_rho ) + 1.0/(4*std::pow(2,0.5)*alpha_mes*alpha_mes);
  double value3  = value3p*value3pp;
  return value1 + value2 - value3;
}

double BottomDecayWidths::F0TOT_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes) ;
  value2 = std::pow(value2,2);
  double value3 = F0_DWAVE(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_DWAVE(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_DWAVE(alpha_rho, alpha_mes), 3) );
  double value1pp = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 2);
  double value1   = value1p*std::pow(value1pp, 2);
  double value2   = value1pp*ARO1_DWAVE(alpha_rho, alpha_mes) + BRO2_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3   = std::pow(15,0.5) / (4.0 * std::pow( ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*pi_val*pi_val;
}

double BottomDecayWidths::I02B0_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_DWAVE(alpha_rho, alpha_lam, alpha_mes) / std::pow(ARO0_DWAVE(alpha_rho, alpha_mes), 3) ;
  double value2 = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value) / std::pow(ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 2);
  double value3 = (3*std::pow(10.,0.5) * ARO1_DWAVE(alpha_rho, alpha_mes)) / (8.* std::pow( ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 5));
  return value1*value2*value3*pi_val*pi_val;
}

double BottomDecayWidths::I01B0TOT_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

// D-wave excitations
// rho wave functions
double BottomDecayWidths::CBARFIN_DWAVE_RHO(double alpha_rho, double alpha_lam){
  double value1 = (std::pow(3,0.75) * std::pow(16,0.5))/(std::pow(15,0.5) * (std::pow(pi_val,0.25)) );
  double value2 = std::pow(( 1.0/alpha_rho ), 3.5);
  double value3 = std::pow(1.0/( pi_val * std::pow(alpha_lam, 2)  ), 0.75);
  return value1 * value2 * value3;// * mycomplex;//define complex if needed
}

//CBARIN_DWAVE --> same as P-wave (CBARIN)
//CMESON_DWAVE --> same as P-wave (CMESON)

double BottomDecayWidths::C0_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3.0*std::pow(3,0.5));
  double value2 = CBARFIN_DWAVE_RHO(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double BottomDecayWidths::ARO0_DWAVE_RHO(double alpha_lam, double alpha_mes){
  double value1 = 1.0/std::pow(alpha_lam,2);  
  double value2 = 1.0/(12.0*std::pow(alpha_mes,2));
  return std::pow(value1+value2,0.5);
}

double BottomDecayWidths::F0_DWAVE_RHO(double k_value, double alpha_rho, double alpha_lam){
  double value1 = std::pow(k_value, 2) /(alpha_rho * alpha_lam) ;
  return (-1.0/3.0) * value1;
}

double BottomDecayWidths::F01_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1p  = 1.0/(2.0*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(6,0.5)*alpha_mes*alpha_mes);
  double value1pp = value1p*value1p;
  double value1ppp= 1./(std::pow(ARO0_DWAVE_RHO(alpha_lam, alpha_mes),2));  
  double value1   = value1pp*value1ppp;
  double value2   = 1.0/(8.0*std::pow(alpha_mes,2));  
  return std::pow(k_value,2)*(value1-value2);
}

double BottomDecayWidths::BRO2_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = 1.0/(2.0*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1.0/(4.0*std::pow(6,0.5)*alpha_mes*alpha_mes);
  double value2 = 1.0/(std::pow(ARO0_DWAVE_RHO(alpha_lam, alpha_mes),2));
  return k_value*(value1*value2*(std::pow(2./3., 0.5)) - 2.0);
}

double BottomDecayWidths::ARO1_DWAVE_RHO(double alpha_lam, double alpha_mes){
  double value1 = std::pow(2,0.5);
  double value2 = 1./(4.0 * std::pow(3,0.5) * std::pow(ARO0_DWAVE_RHO(alpha_lam, alpha_mes),2) * std::pow(alpha_mes,2) );
  return value2*(std::pow(2./3.,0.5)) - value1;
}

double BottomDecayWidths::ARO2_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4.0/(std::pow(alpha_rho,2));
  double value1pp = 1.0/(std::pow(alpha_mes,2));
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1.0/(48.0 * std::pow(alpha_mes,4) * std::pow(ARO0_DWAVE_RHO(alpha_lam, alpha_mes),2) );
  return std::pow(value1-value2, 0.5);
}

double BottomDecayWidths::BRO1_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1  = k_value / ( 4*std::pow(2,0.5) * std::pow(alpha_mes,2) );
  double value2  = k_value / ( 2*std::pow(2,0.5) * alpha_rho * alpha_lam );
  double value3p = k_value / ( 4*std::pow(3,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0_DWAVE_RHO(alpha_lam, alpha_mes), 2) );
  double value3pp= 1.0/( 2*std::pow(6,0.5)*alpha_lam*alpha_rho ) + 1.0/(4*std::pow(6,0.5)*alpha_mes*alpha_mes);
  double value3  = value3p*value3pp;
  return (-1.0)*(value1 + value2 - value3);
}

double BottomDecayWidths::F0TOT_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = BRO1_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO2_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes) ;
  value2 = std::pow(value2,2);
  double value3 = F0_DWAVE_RHO(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_DWAVE_RHO(alpha_lam, alpha_mes), 3) );
  double value1pp = BRO1_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 2);
  double value1   = value1p*std::pow(value1pp, 2);
  double value2   = value1pp*ARO1_DWAVE_RHO(alpha_lam, alpha_mes) + BRO2_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3   = std::pow(15,0.5) / (4.0 * std::pow( ARO2_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*pi_val*pi_val;
}

double BottomDecayWidths::I02B0_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes) / std::pow(ARO0_DWAVE_RHO(alpha_lam, alpha_mes), 3) ;
  double value2 = BRO1_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value) / std::pow(ARO2_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 2);
  double value3 = (3*std::pow(10.,0.5) * ARO1_DWAVE_RHO(alpha_lam, alpha_mes)) / (8.* std::pow( ARO2_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 5));
  return value1*value2*value3*pi_val*pi_val;
}

double BottomDecayWidths::I01B0TOT_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}


// RADIAL excitations
// lambda wave functions
double BottomDecayWidths::CBARFIN_RADIAL(double alpha_rho, double alpha_lam){
  double value1 = std::pow(3,0.75) * std::pow(2./3.,0.5);
  double value2 = std::pow(1.0/( std::pow(pi_val, 4) * std::pow(alpha_lam, 2) *std::pow(alpha_rho, 2) ), 0.75);
  return value1 * value2 ;// * mycomplex;//define complex if needed
}

//CBARIN_RADIAL --> same as P-wave (CBARIN)
//CMESON_RADIAL --> same as P-wave (CMESON)

double BottomDecayWidths::C0_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3.0*std::pow(3,0.5));
  double value2 = CBARFIN_RADIAL(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double BottomDecayWidths::ARO0_RADIAL(double alpha_rho, double alpha_mes){
  return ARO0_DWAVE(alpha_rho, alpha_mes);}

double BottomDecayWidths::F0_RADIAL(double k_value, double alpha_rho, double alpha_lam){
  return F0_DWAVE(k_value, alpha_rho, alpha_lam);
}

double BottomDecayWidths::F01_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1.0/(2.0*std::pow(2,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(2,0.5)*alpha_mes*alpha_mes);  
  double value1pp = value1p*value1p;  
  double value1ppp= 1./(std::pow(ARO0_RADIAL(alpha_rho, alpha_mes),2));  
  double value1 = value1pp*value1ppp;
  double value2 = 1.0/(8.0*std::pow(alpha_mes,2));
  return std::pow(k_value,2)*(value1-value2);
}

double BottomDecayWidths::BRO2_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = 1.0/(2.0*alpha_lam*alpha_rho) + 1.0/(4.0*alpha_mes*alpha_mes);
  double value2 = 1.0/(std::pow(ARO0_RADIAL(alpha_rho, alpha_mes),2));
  return k_value*(value1*value2 - 2.0);
}

double BottomDecayWidths::ARO1_RADIAL(double alpha_rho, double alpha_mes){
  double value1 = std::pow(2./3., 0.5);
  double value2 = 1./(2.0 * std::pow(6,0.5) * std::pow(ARO0_RADIAL(alpha_rho, alpha_mes),2) * std::pow(alpha_mes,2) );
  return value1 - value2;
}

double BottomDecayWidths::ARO2_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4.0/(std::pow(alpha_lam,2));
  double value1pp = 1.0/(3.0*std::pow(alpha_mes,2));
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1.0/(48.0 * std::pow(alpha_mes,4) * std::pow(ARO0_RADIAL(alpha_rho, alpha_mes),2) );
  return std::pow(value1-value2, 0.5);
}

double BottomDecayWidths::BRO1_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1  = k_value / ( 4*std::pow(6,0.5) * std::pow(alpha_mes,2) );
  double value2  = k_value / ( 2*std::pow(6,0.5) * alpha_rho * alpha_lam );
  double value3p = k_value / ( 4*std::pow(3,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0_RADIAL(alpha_rho, alpha_mes), 2) );
  double value3pp= 1.0/( 2*std::pow(2,0.5)*alpha_lam*alpha_rho ) + 1.0/(4*std::pow(2,0.5)*alpha_mes*alpha_mes);
  double value3  = value3p*value3pp;
  return value1 + value2 - value3;
}

double BottomDecayWidths::F0TOT_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = BRO1_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO2_RADIAL(alpha_rho, alpha_lam, alpha_mes) ;
  value2 = std::pow(value2,2);
  double value3 = F0_RADIAL(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_RADIAL(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_RADIAL(alpha_rho, alpha_mes), 3) );
  double value2p = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 2);
  double value2 = value2p*ARO1_RADIAL(alpha_rho, alpha_mes) + BRO2_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3 = std::pow(12.,0.5)/(4 * std::pow(ARO2_RADIAL(alpha_rho, alpha_lam, alpha_mes) , 3));
  double value4 = 1.5 - std::pow(value2p,2)/std::pow(alpha_lam,2);
  return value1 * value2 * value3 * value4 * std::pow(pi_val, 2.5);
}

double BottomDecayWidths::I02B0_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_RADIAL(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_RADIAL(alpha_rho, alpha_mes), 3) );  
  double value2p = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 2);  
  double value2pp  = value2p*ARO1_RADIAL(alpha_rho, alpha_mes) + BRO2_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2ppp = (3.*std::pow(12.,0.5))/(8 * std::pow(alpha_lam,2) * std::pow(ARO2_RADIAL(alpha_rho, alpha_lam, alpha_mes) , 5) );
  double value2 = (-1.0)*value2pp * value2ppp;
  double value3 = (-0.5*std::pow(3.,0.5)*value2p*ARO1_RADIAL(alpha_rho, alpha_mes))/(std::pow(alpha_lam,2) * std::pow(ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 5));  
  return value1 * (value2 + value3 )* std::pow(pi_val, 2.5);
}

double BottomDecayWidths::I01B0TOT_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

// RADIAL excitations
// rho wave functions
double BottomDecayWidths::CBARFIN_RADIAL_RHO(double alpha_rho, double alpha_lam){
  double value1 = std::pow(3,0.75) * std::pow(2./3.,0.5);
  double value2 = std::pow(1.0/( std::pow(pi_val, 4) * std::pow(alpha_lam, 2) *std::pow(alpha_rho, 2) ), 0.75);
  return value1 * value2 ;// * mycomplex;//define complex if needed
}

//CBARIN_RADIAL_RHO --> same as P-wave (CBARIN)
//CMESON_RADIAL_RHO --> same as P-wave (CMESON)

double BottomDecayWidths::C0_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3.0*std::pow(3,0.5));
  double value2 = CBARFIN_RADIAL_RHO(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double BottomDecayWidths::ARO0_RADIAL_RHO(double alpha_lam, double alpha_mes){
  return ARO0_DWAVE_RHO(alpha_lam, alpha_mes);}

double BottomDecayWidths::F0_RADIAL_RHO(double k_value, double alpha_rho, double alpha_lam){
  return F0_DWAVE_RHO(k_value, alpha_rho, alpha_lam);
}

double BottomDecayWidths::F01_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1.0/(2.0*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(6,0.5)*alpha_mes*alpha_mes);  
  double value1pp = value1p*value1p;  
  double value1ppp= 1./(std::pow(ARO0_RADIAL_RHO(alpha_lam, alpha_mes),2));  
  double value1 = value1pp*value1ppp;
  double value2 = 1.0/(8.0*std::pow(alpha_mes,2));
  return std::pow(k_value,2)*(value1-value2);
}

double BottomDecayWidths::BRO2_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = (std::pow(2./3., 0.5))/(2.0*std::pow(6.,0.5)*alpha_lam*alpha_rho) + 1.0/(4.0*std::pow(6.,0.5)*alpha_mes*alpha_mes);
  double value2 = 1.0/(std::pow(ARO0_RADIAL_RHO(alpha_lam, alpha_mes),2));
  return k_value*(value1*value2 - 2.0);
}

double BottomDecayWidths::ARO1_RADIAL_RHO(double alpha_lam, double alpha_mes){
  double value1 = std::pow(2., 0.5);
  double value2 = (std::pow(2./3., 0.5))/(4.0 * std::pow(3,0.5) * std::pow(ARO0_RADIAL_RHO(alpha_lam, alpha_mes),2) * std::pow(alpha_mes,2) );
  return value2 - value1;
}

double BottomDecayWidths::ARO2_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4.0/(std::pow(alpha_rho,2));
  double value1pp = 1.0/(std::pow(alpha_mes,2));
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1.0/(48.0 * std::pow(alpha_mes,4) * std::pow(ARO0_RADIAL_RHO(alpha_lam, alpha_mes),2) );
  return std::pow(value1-value2, 0.5);
}

double BottomDecayWidths::BRO1_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1  = k_value / ( 4*std::pow(2,0.5) * std::pow(alpha_mes,2) );
  double value2  = k_value / ( 2*std::pow(2,0.5) * alpha_rho * alpha_lam );
  double value3p = k_value / ( 4*std::pow(3,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0_RADIAL_RHO(alpha_lam, alpha_mes), 2) );
  double value3pp= 1.0/( 2*std::pow(6,0.5)*alpha_lam*alpha_rho ) + 1.0/(4*std::pow(6,0.5)*alpha_mes*alpha_mes);
  double value3  = value3p*value3pp;
  return (-1.0)*(value1 + value2 - value3);
}

double BottomDecayWidths::F0TOT_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = BRO1_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO2_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes) ;
  value2 = std::pow(value2,2);
  double value3 = F0_RADIAL_RHO(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3;
}

double BottomDecayWidths::I01B0_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_RADIAL_RHO(alpha_lam, alpha_mes), 3) );  
  double value2p = BRO1_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes), 2);  
  double value2 = value2p*ARO1_RADIAL_RHO(alpha_lam, alpha_mes) + BRO2_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3 = std::pow(12.,0.5)/(4 * std::pow(ARO2_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes) , 3));
  double value4 = 1.5 - std::pow(value2p,2)/std::pow(alpha_rho,2);
  return value1 * value2 * value3 * value4 * std::pow(pi_val, 2.5);
}

double BottomDecayWidths::I02B0_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_RADIAL_RHO(alpha_lam, alpha_mes), 3) );  
  double value2p = BRO1_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 2);  
  double value2pp  = value2p*ARO1_RADIAL_RHO(alpha_lam, alpha_mes) + BRO2_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2ppp = (3.*std::pow(12.,0.5))/(8 * std::pow(alpha_rho,2) * std::pow(ARO2_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes) , 5) );
  double value2 = (-1.0)*value2pp * value2ppp;
  double value3 = (-std::pow(3.,0.5)*value2p*ARO1_RADIAL_RHO(alpha_lam, alpha_mes))/(2.*std::pow(alpha_rho,2) * std::pow(ARO2_DWAVE_RHO(alpha_rho, alpha_lam, alpha_mes), 5));
  return value1 * (value2 + value3 )* std::pow(pi_val, 2.5);
}

double BottomDecayWidths::I01B0TOT_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_RADIAL_RHO(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

// Mixed excitations
double BottomDecayWidths::CBARFIN_MIXED(double alpha_rho, double alpha_lam){
  double value1 = (std::pow(3,0.75) * 8)/(3 * std::pow(pi_val,0.5) );
  double value2 = std::pow((1.0/alpha_rho ), 2.5);
  double value3 = std::pow((1.0/alpha_lam ), 2.5);
  return value1 * value2 * value3;// * mycomplex;//define complex if needed
}

//CBARIN_DWAVE --> same as P-wave (CBARIN)
//CMESON_DWAVE --> same as P-wave (CMESON)

double BottomDecayWidths::C0_MIXED(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3.0*std::pow(3,0.5));
  double value2 = CBARFIN_MIXED(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double BottomDecayWidths::ARO0_MIXED(double alpha_lam, double alpha_mes){
  double value1 = 1.0/std::pow(alpha_lam,2);  
  double value2 = 1.0/(12.0*std::pow(alpha_mes,2));
  return std::pow(value1+value2,0.5);
}

double BottomDecayWidths::F0_MIXED(double k_value, double alpha_rho, double alpha_lam){
  double value1 = std::pow(k_value, 2) / (alpha_rho * alpha_lam) ;
  return (-1.0/3.0) * value1;
}

double BottomDecayWidths::F01_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1p  = 1.0/(2.0*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(6,0.5)*alpha_mes*alpha_mes);
  double value1pp = value1p*value1p;
  double value1ppp= 1./(std::pow(ARO0_MIXED(alpha_lam, alpha_mes),2));  
  double value1   = value1pp*value1ppp;
  double value2   = 1.0/(8.0*std::pow(alpha_mes,2));  
  return std::pow(k_value,2)*(value1-value2);
}

double BottomDecayWidths::BRO2_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = 1.0/(2.0*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1.0/(4.0*std::pow(6,0.5)*alpha_mes*alpha_mes);
  double value2 = 1.0/(std::pow(ARO0_MIXED(alpha_lam, alpha_mes),2));
  return k_value*(value1*value2*(std::pow(2./3., 0.5)) - 2.0);
}

double BottomDecayWidths::ARO1_MIXED(double alpha_lam, double alpha_mes){
  double value1 = std::pow(2,0.5);
  double value2 = (std::pow(2./3.,0.5))/(4.0 * std::pow(3,0.5) * std::pow(ARO0_MIXED(alpha_lam, alpha_mes),2) * std::pow(alpha_mes,2) );
  return value2 - value1;
}

double BottomDecayWidths::ARO2_MIXED(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4.0/(std::pow(alpha_rho,2));
  double value1pp = 1.0/(std::pow(alpha_mes,2));
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1.0/(48.0 * std::pow(alpha_mes,4) * std::pow(ARO0_MIXED(alpha_lam, alpha_mes),2) );
  return std::pow(value1-value2, 0.5);
}

double BottomDecayWidths::BRO1_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1  = k_value / ( 4*std::pow(2,0.5) * std::pow(alpha_mes,2) );
  double value2  = k_value / ( 2*std::pow(2,0.5) * alpha_rho * alpha_lam );
  double value3p = k_value / ( 4*std::pow(3,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0_MIXED(alpha_lam, alpha_mes), 2) );
  double value3pp= 1.0/( 2*std::pow(6,0.5)*alpha_lam*alpha_rho ) + 1.0/(4*std::pow(6,0.5)*alpha_mes*alpha_mes);
  double value3  = value3p*value3pp;
  return (-1.0)*(value1 + value2 + value3);
}

double BottomDecayWidths::ARO4_MIXED(double alpha_lam, double alpha_mes){
  double value1 = 1./(4.0 * std::pow(3,0.5) * std::pow(ARO0_MIXED(alpha_lam, alpha_mes),2) * std::pow(alpha_mes,2) );
  return value1;
}

double BottomDecayWidths::BRO4_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = 1.0/(2.0*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1.0/(4.0*std::pow(6,0.5)*alpha_mes*alpha_mes);
  double value2 = 1.0/(std::pow(ARO0_MIXED(alpha_lam, alpha_mes),2));
  return k_value*(value1*value2);
}

double BottomDecayWidths::ARO5_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p = BRO1_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value) / std::pow(ARO2_MIXED(alpha_rho, alpha_lam, alpha_mes), 2);
  double value1 = value1p*ARO1_MIXED(alpha_lam, alpha_mes);
  double value2 = BRO2_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1 + value2;
}

double BottomDecayWidths::ARO6_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p = BRO1_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value) / std::pow(ARO2_MIXED(alpha_rho, alpha_lam, alpha_mes), 2);
  double value1 = value1p*ARO4_MIXED(alpha_lam, alpha_mes);
  double value2 = BRO4_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1 + value2;
}

double BottomDecayWidths::F0TOT_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = BRO1_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO2_MIXED(alpha_rho, alpha_lam, alpha_mes) ;
  value2 = std::pow(value2,2);
  double value3 = F0_MIXED(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_MIXED(alpha_rho, alpha_lam, alpha_mes) / (4.0*std::pow(ARO0_MIXED(alpha_lam, alpha_mes), 3) );
  double value2 = BRO1_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_MIXED(alpha_rho, alpha_lam, alpha_mes), 2);
  double value3 = std::pow(3./4.,0.5) / (std::pow( ARO2_MIXED(alpha_rho, alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*std::pow(pi_val,1.5);
}

double BottomDecayWidths::I02B0_MIXED(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = C0_MIXED(alpha_rho, alpha_lam, alpha_mes) / (1.0*std::pow(ARO0_MIXED(alpha_lam, alpha_mes), 3) );
  double value2 = (3.0*std::pow(3./4.,0.5)) / (8.0*std::pow( ARO2_MIXED(alpha_rho, alpha_lam, alpha_mes), 5) );
  return value1*value2*std::pow(pi_val,1.5);
}

double BottomDecayWidths::I01B0TOT_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_MIXED(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_MIXED(alpha_rho, alpha_lam, alpha_mes);
  return value1*value2;
}

// SECOND DIAGRAM
// D-WAVE LAMBDA
// same as   C0_DWAVE
double BottomDecayWidths::ARO0_DWAVE_SECOND(double alpha_rho){
  return 1./alpha_rho;
}

double BottomDecayWidths::F0_DWAVE_SECOND(double k_value, double alpha_rho, double alpha_lam){
  double value1 = std::pow(k_value, 2) /(alpha_rho * alpha_lam);
  return (-1.0/3.0) * value1;
}

double BottomDecayWidths::F01_DWAVE_SECOND(double alpha_mes, double k_value){
  double value1 = std::pow(k_value, 2);
  double value2 = (-1.0)/(8 * std::pow(alpha_mes, 2));
  return value1 * value2;
}

double BottomDecayWidths::BRO1_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = (-1.0)/(2.*std::pow(6, 0.5)*alpha_mes*alpha_mes);
  double value2 = (-1.0)/(std::pow(6, 0.5)*alpha_rho*alpha_lam);
  return k_value*(value1 + value2);
}

double BottomDecayWidths::ARO1_DWAVE_SECOND(double alpha_lam, double alpha_mes){
  double value1 = 1.0/(alpha_lam*alpha_lam);
  double value2 = 1.0/(3*alpha_mes*alpha_mes);
  return std::pow(value1 + value2, 0.5);
}

double BottomDecayWidths::F0TOT_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_DWAVE_SECOND(alpha_mes, k_value);
  double value2 = BRO1_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO1_DWAVE_SECOND(alpha_lam, alpha_mes) ;
  value2 = std::pow(value2, 2);
  double value3 = F0_DWAVE_SECOND(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_DWAVE(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_DWAVE_SECOND(alpha_rho), 3) );
  double value1pp = BRO1_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO1_DWAVE_SECOND(alpha_lam, alpha_mes), 2);
  double value1   = value1p*std::pow(value1pp, 2);
  double value2   = value1pp*((2*std::pow(6, 0.5))/3.) +  2*k_value;
  double value3   = std::pow(15, 0.5) / (4. * std::pow( ARO1_DWAVE_SECOND(alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*pi_val*pi_val;
}

double BottomDecayWidths::I02B0_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_DWAVE(alpha_rho, alpha_lam, alpha_mes) / std::pow(ARO0_DWAVE_SECOND(alpha_rho), 3) ;
  double value2 = BRO1_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value) / std::pow(ARO1_DWAVE_SECOND(alpha_lam, alpha_mes), 2);
  double value3 = std::pow(15.,0.5) / (2.* std::pow( ARO1_DWAVE_SECOND(alpha_lam, alpha_mes), 5));
  return value1*value2*value3*pi_val*pi_val;
}

double BottomDecayWidths::I01B0TOT_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_DWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}


// SECOND DIAGRAM
// D-WAVE RHO
// C0 same as C0_DWAVE_RHO
double BottomDecayWidths::C0_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = std::pow(CBARFIN_DWAVE_RHO(alpha_rho, alpha_lam), 2);
  double value2 = CMESON(alpha_mes);
  return (1/(3.*std::pow(3, 0.5)))*value1*value2;
}

double BottomDecayWidths::ARO0_DWAVE_RHO_SECOND(double alpha_rho){
  return 1./alpha_rho;
}

double BottomDecayWidths::F0_DWAVE_RHO_SECOND(double k_value, double alpha_rho, double alpha_lam){
  double value1 = std::pow(k_value, 2) /(alpha_rho * alpha_lam);
  return (-1.0/3.0) * value1;
}

double BottomDecayWidths::F01_DWAVE_RHO_SECOND(double alpha_mes, double k_value){
  double value1 = std::pow(k_value, 2);
  double value2 = (-1.0)/(8 * std::pow(alpha_mes, 2));
  return value1 * value2;
}

double BottomDecayWidths::BRO1_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = (-1.0)/(2.*std::pow(6, 0.5)*alpha_mes*alpha_mes);
  double value2 = (-1.0)/(std::pow(6, 0.5)*alpha_rho*alpha_lam);
  return k_value*(value1 + value2);
}

double BottomDecayWidths::ARO1_DWAVE_RHO_SECOND(double alpha_lam, double alpha_mes){
  double value1 = 1.0/(alpha_lam*alpha_lam);
  double value2 = 1.0/(3*alpha_mes*alpha_mes);
  return std::pow(value1 + value2, 0.5);
}

double BottomDecayWidths::F0TOT_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_DWAVE_RHO_SECOND(alpha_mes, k_value);
  double value2 = BRO1_DWAVE_RHO_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO1_DWAVE_RHO_SECOND(alpha_lam, alpha_mes) ;
  value2 = std::pow(value2, 2);
  double value3 = F0_DWAVE_RHO_SECOND(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_DWAVE_RHO_SECOND(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_DWAVE_RHO_SECOND(alpha_rho), 7) );
  double value1pp = BRO1_DWAVE_RHO_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO1_DWAVE_RHO_SECOND(alpha_lam, alpha_mes), 2);
  double value1   = value1p * (15./16.);
  double value2   = value1pp*((2*std::pow(6, 0.5))/3.) + 2*k_value;
  double value3   = std::pow(0.75, 0.5) / (std::pow( ARO1_DWAVE_RHO_SECOND(alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*std::pow(pi_val, 1.5);
}

double BottomDecayWidths::I01B0TOT_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_DWAVE_RHO_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_DWAVE_RHO_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

// SECOND DIAGRAM
// P-WAVE lambda
// CBARIN the same
double BottomDecayWidths::C0_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = CBARFIN_PWAVE(alpha_rho, alpha_lam);
  double value2 = CBARIN(alpha_rho, alpha_lam);
  double value3 = CMESON(alpha_mes);
  return (1./(3.*std::pow(3, 0.5)))*value1*value2*value3;
}

double BottomDecayWidths::ARO0_PWAVE_SECOND(double alpha_rho){
  return 1./alpha_rho;
}

double BottomDecayWidths::F0_PWAVE_SECOND(double k_value, double alpha_rho, double alpha_lam){
  double value1 = std::pow(k_value, 2) /(alpha_rho * alpha_lam);
  return (-1.0/3.0) * value1;
}

double BottomDecayWidths::F01_PWAVE_SECOND(double alpha_mes, double k_value){
  double value1 = std::pow(k_value, 2);
  double value2 = (-1.0)/(8 * std::pow(alpha_mes, 2));
  return value1 * value2;
}

double BottomDecayWidths::BRO1_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = (-1.0)/(2.*std::pow(6, 0.5)*alpha_mes*alpha_mes);
  double value2 = (-1.0)/(std::pow(6, 0.5)*alpha_rho*alpha_lam);
  return k_value*(value1 + value2);
}

double BottomDecayWidths::ARO1_PWAVE_SECOND(double alpha_lam, double alpha_mes){
  double value1 = 1.0/(alpha_lam*alpha_lam);
  double value2 = 1.0/(3*alpha_mes*alpha_mes);
  return std::pow(value1 + value2, 0.5);
}

double BottomDecayWidths::F0TOT_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_PWAVE_SECOND(alpha_mes, k_value);
  double value2 = BRO1_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO1_PWAVE_SECOND(alpha_lam, alpha_mes) ;
  value2 = std::pow(value2, 2);
  double value3 = F0_PWAVE_SECOND(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_PWAVE_SECOND(alpha_rho), 3) );
  double value1pp = BRO1_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO1_PWAVE_SECOND(alpha_lam, alpha_mes), 2);
  double value1   = value1p * value1pp;
  double value2   = value1pp*((2*std::pow(6, 0.5))/3.) + 2*k_value;
  double value3   = std::pow(0.75, 1.0) / (std::pow( ARO1_PWAVE_SECOND(alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*std::pow(pi_val, 2);
}

double BottomDecayWidths::I02B0_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_PWAVE_SECOND(alpha_rho), 3) );
  double value1   = value1p * ((2*std::pow(6, 0.5))/3.);
  double value2   = (3./8.)/ (std::pow( ARO1_PWAVE_SECOND(alpha_lam, alpha_mes), 5) );
  return value1*value2*std::pow(pi_val, 2);
}

double BottomDecayWidths::I01B0TOT_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_PWAVE_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

// SECOND DIAGRAM
// MIXED
// CMESON the same
double BottomDecayWidths::C0_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = CBARFIN_MIXED(alpha_rho, alpha_lam);
  double value2 = CBARFIN_PWAVE_RHO(alpha_rho, alpha_lam);
  double value3 = CMESON(alpha_mes);
  return (1/(3.*std::pow(3, 0.5)))*value1*value2*value3;
}

double BottomDecayWidths::ARO0_MIXED_SECOND(double alpha_rho){
  return 1./alpha_rho;
}

double BottomDecayWidths::F0_MIXED_SECOND(double k_value, double alpha_rho, double alpha_lam){
  double value1 = std::pow(k_value, 2) /(alpha_rho * alpha_lam);
  return (-1.0/3.0) * value1;
}

double BottomDecayWidths::F01_MIXED_SECOND(double alpha_mes, double k_value){
  double value1 = std::pow(k_value, 2);
  double value2 = (-1.0)/(8 * std::pow(alpha_mes, 2));
  return value1 * value2;
}

double BottomDecayWidths::BRO1_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = (-1.0)/(2.*std::pow(6, 0.5)*alpha_mes*alpha_mes);
  double value2 = (-1.0)/(std::pow(6, 0.5)*alpha_rho*alpha_lam);
  return k_value*(value1 + value2);
}

double BottomDecayWidths::ARO1_MIXED_SECOND(double alpha_lam, double alpha_mes){
  double value1 = 1.0/(alpha_lam*alpha_lam);
  double value2 = 1.0/(3*alpha_mes*alpha_mes);
  return std::pow(value1 + value2, 0.5);
}

double BottomDecayWidths::F0TOT_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_MIXED_SECOND(alpha_mes, k_value);
  double value2 = BRO1_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO1_MIXED_SECOND(alpha_lam, alpha_mes) ;
  value2 = std::pow(value2, 2);
  double value3 = F0_MIXED_SECOND(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_MIXED_SECOND(alpha_rho), 5) );
  double value1pp = BRO1_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO1_MIXED_SECOND(alpha_lam, alpha_mes), 2);
  double value1   = value1p * value1pp * (3./8.);
  double value2   = value1pp*((2*std::pow(6, 0.5))/3.) + 2*k_value;
  double value3   = std::pow(0.75, 1.0) / (std::pow( ARO1_MIXED_SECOND(alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*std::pow(pi_val, 1);
}

double BottomDecayWidths::I02B0_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_MIXED_SECOND(alpha_rho), 5) );
  double value1   = value1p * (3./8.);
  double value2   = ((std::pow(6, 0.5))/4.)/ (std::pow( ARO1_MIXED_SECOND(alpha_lam, alpha_mes), 5) );
  return value1*value2*std::pow(pi_val, 1);
}

double BottomDecayWidths::I01B0TOT_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_MIXED_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

//SECOND DIAGRAM
//RADIAL excitations
//lambda wave functions
//CBARFIN_RADIAL_SECOND--> same as D-wave 
//CBARIN_RADIAL_SECOND --> same as P-wave (CBARIN)
//CMESON_RADIAL_SECOND --> same as P-wave (CMESON)

double BottomDecayWidths::C0_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3.0*std::pow(3,0.5));
  double value2 = CBARFIN_RADIAL(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double BottomDecayWidths::ARO0_RADIAL_SECOND(double alpha_rho){
  return 1./alpha_rho;
}

double BottomDecayWidths::F0_RADIAL_SECOND(double k_value, double alpha_rho, double alpha_lam){
  double value1 = std::pow(k_value, 2) /(alpha_rho * alpha_lam);
  return (-1.0/3.0) * value1;
}

double BottomDecayWidths::F01_RADIAL_SECOND(double alpha_mes, double k_value){
  double value1 = std::pow(k_value, 2);
  double value2 = (-1.0)/(8 * std::pow(alpha_mes, 2));
  return value1 * value2;
}

double BottomDecayWidths::BRO1_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = (-1.0)/(2.*std::pow(6, 0.5)*alpha_mes*alpha_mes);
  double value2 = (-1.0)/(std::pow(6, 0.5)*alpha_rho*alpha_lam);
  return k_value*(value1 + value2);
}

double BottomDecayWidths::ARO1_RADIAL_SECOND(double alpha_lam, double alpha_mes){
  double value1 = 1.0/(alpha_lam*alpha_lam);
  double value2 = 1.0/(3*alpha_mes*alpha_mes);
  return std::pow(value1 + value2, 0.5);
}

double BottomDecayWidths::F0TOT_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_RADIAL_SECOND(alpha_mes, k_value);
  double value2 = BRO1_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO1_RADIAL_SECOND(alpha_lam, alpha_mes) ;
  value2 = std::pow(value2, 2);
  double value3 = F0_RADIAL_SECOND(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double BottomDecayWidths::I01B0_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1   = C0_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_RADIAL_SECOND(alpha_rho), 3) );
  double value1pp = BRO1_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO1_RADIAL_SECOND(alpha_lam, alpha_mes), 2);
  double value2   = 1.5 - std::pow(value1pp/alpha_lam, 2);
  double value3   = value1pp*((2*std::pow(6, 0.5))/3.) + 2*k_value;
  double value4   = std::pow(0.75, 0.5) / (std::pow( ARO1_RADIAL_SECOND(alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*value4*std::pow(pi_val, 2.5);
}

double BottomDecayWidths::I02B0_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1   = C0_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_RADIAL_SECOND(alpha_rho), 3) );
  double value1pp = BRO1_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO1_RADIAL_SECOND(alpha_lam, alpha_mes), 2);
  double value2p  = ( (2*std::pow(3, 0.5)) * ( value1pp*((2*std::pow(6, 0.5))/3.) + 2*k_value) ) / std::pow(alpha_lam, 2);  
  double value2pp = ( 8*std::pow(6, 0.5) * value1pp ) / (3. * std::pow(3, 0.5) * std::pow(alpha_lam, 2) );
  double value2   = (-1.0)*(value2p+value2pp);
  double value3   = (3./8.)/ (std::pow( ARO1_MIXED_SECOND(alpha_lam, alpha_mes), 5) );
  return value1*value2*value3*std::pow(pi_val, 2.5);
}

double BottomDecayWidths::I01B0TOT_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double BottomDecayWidths::I02B0TOT_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_RADIAL_SECOND(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

#endif
