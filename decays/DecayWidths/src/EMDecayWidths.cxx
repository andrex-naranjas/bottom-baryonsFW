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

double EMDecayWidths::execute(double ma_val, double sa_val, double ja_val, double la_val, double sla_val, double lla_val, double lra_val,
			      double mb_val, double sb_val, double jb_val, double lb_val, double slb_val, double llb_val, double lrb_val,
			      double al_val, double ar_val,
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

  alpha_rho = ar_val;
  alpha_lam = al_val;

  //fetch quantum numbers and projections
  L   = 0;  mL = getMomentumProjections(L); //check if this is always the case

  // sa_val = 3./2.;
  // ja_val = 3./2.;
  // sla_val = 1.;

  // sb_val = 1./2.;
  // jb_val = 1./2.;
  // slb_val = 0.;

  // sa_val    = 0.5;
  // sla_val   = 1.;  
  // la_val    = 1.;  
  // lla_val   = 1.;  
  // lra_val   = 0.;  
  // ja_val    = 0.5;
  

  // sb_val     = 0.5;
  // slb_val    = 1.;  
  // lb_val     = 0.;  
  // llb_val    = 0.;  
  // lrb_val    = 0.;  
  // jb_val     = 0.5;
  

  SA  = sa_val;    mSA  = getMomentumProjections(SA);
  JA  = ja_val;    mJA  = getMomentumProjections(JA);
  SlA = sla_val;   mSlA = getMomentumProjections(SlA);
  LA  = la_val;    mLA  = getMomentumProjections(LA);
  LlA = lla_val;   mLlA = getMomentumProjections(LlA);
  LrA = lra_val;   mLrA = getMomentumProjections(LrA);

  SB  = sb_val;    mSB  = getMomentumProjections(SB);
  JB  = jb_val;    mJB  = getMomentumProjections(JB);
  SlB = slb_val;   mSlB = getMomentumProjections(SlB);
  LB  = lb_val;    mLB  = getMomentumProjections(LB);
  LlB = llb_val;   mLlB = getMomentumProjections(LlB);
  LrB = lrb_val;   mLrB = getMomentumProjections(LrB);

  S1  = 0.5;       mS1  = getMomentumProjections(S1);
  S2  = 0.5;       mS2  = getMomentumProjections(S2);
  S3  = 0.5;       mS3  = getMomentumProjections(S3);

  double k_value; k_value = K(MA, MB);
  double EB_value = EB(MB, k_value);
  
  double fi2_value  = FI2(EB_value, MA, k_value);
  double decayWidth = 0;//DecayWidth(flav_coup, fi2_value, sum_value);
  decayWidth = 1;

  double test_integral = 1; //SPINFLIP_U1_GS_GS();
  std::vector<double> flavor_vector = FlavorVector(1.);

  double test_sum = ANGULAR_SUM_SQUARED(alpha_rho, alpha_lam, k_value, flavor_vector, modeExcitation);
  std::cout<<test_sum<<"    test sum"<<std::endl;


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

double EMDecayWidths::ClebschGordan(WignerSymbols *m_wigner,
				    double l1, double l2, double l3,
				    double m1, double m2, double m3){
  double coef = std::pow(-1.0, l2-l1-m3) * std::pow(2.0*l3 + 1.0, 0.5);
  double three_j = m_wigner->wigner3j(l1, l2, l3, m1, m2, (-1.0)*m3);
  return coef * three_j;
}

std::vector<double> EMDecayWidths::FlavorVector(double charge){
  double mu_qu = (+1.) * ((2./3.) * std::sqrt(1./137.)/(2. * mlight));
  double mu_qd = (-1.) * ((1./3.) * std::sqrt(1./137.)/(2. * mlight));
  double mu_qs = (-1.) * ((1./3.) * std::sqrt(1./137.)/(2. * mlight));
  double mu_qb = (-1.) * ((1./3.) * std::sqrt(1./137.)/(2. * mbottom));
  std::vector<double> flavor_vector; flavor_vector.clear();
  flavor_vector.push_back(mu_qb);
  if(charge==0){
    flavor_vector.push_back( 0.5*(mu_qd + mu_qs) );
    flavor_vector.push_back( 0.5*(mu_qd + mu_qs) );
  }else if(charge==1){
    flavor_vector.push_back( 0.5*(mu_qd + mu_qs) );
    flavor_vector.push_back( 0.5*(mu_qd + mu_qs) );
  }
  return flavor_vector;    
}

double EMDecayWidths::ANGULAR_SUM_SQUARED(double alpha_rho, double alpha_lam, double k_value, std::vector<double> flavor_vector, int excMode){

  WignerSymbols *m_wigner = new WignerSymbols();
  double dummy = 0;
  double AMP1_1 = 0.; double AMP1_2 = 0.; double AMP1 = 0.;
  double AMP2_1 = 0.; double AMP2_2 = 0.; double AMP2 = 0.; 
  double AMP3_1 = 0.; double AMP3_2 = 0.; double AMP3 = 0.;
  double TOT_AMP = 0.; double SUM_SQUARED_AMP = 0.; 

  for(int iMJA = 0;  iMJA<(int)mJA.size(); iMJA++){//SUM SQUARED
    TOT_AMP = 0.;
    AMP2_1 = 0.; AMP2_2 = 0.; AMP2 = 0.;
    AMP1_1 = 0.; AMP1_2 = 0.; AMP1 = 0.;
    AMP3_1 = 0.; AMP3_2 = 0.; AMP3 = 0.;    
    for(int iMSA = 0;  iMSA<(int)mSA.size();  iMSA++)// AMP1
      for(int iMSlA = 0; iMSlA<(int)mSlA.size(); iMSlA++)
	for(int iMLA = 0; iMLA<(int)mLA.size(); iMLA++)
	  for(int iMLlA = 0; iMLlA<(int)mLlA.size(); iMLlA++)
	    for(int iMLrA = 0; iMLrA<(int)mLrA.size(); iMLrA++)
	      for(int iMSB = 0;  iMSB<(int)mSB.size(); iMSB++)
		for(int iMSlB = 0;  iMSlB<(int)mSlB.size(); iMSlB++)
		  for(int iMLB = 0;  iMLB<(int)mLB.size(); iMLB++)
		    for(int iMLlB = 0; iMLlB<(int)mLlB.size(); iMLlB++)
		      for(int iMLrB = 0; iMLrB<(int)mLrB.size(); iMLrB++)
			for(int iMJB = 0;  iMJB<(int)mJB.size(); iMJB++)
			  for(int iMS1 = 0; iMS1 <(int)mS1.size(); iMS1++)
			    for(int iMS2 = 0; iMS2 <(int)mS2.size(); iMS2++)
			      for(int iMS3 = 0; iMS3 <(int)mS3.size(); iMS3++){
				dummy = U1_rho_lambda(k_value, alpha_rho, alpha_lam, mLrA.at(iMLrA), mLlA.at(iMLlA), excMode)*
				  std::sqrt((S1 + mS1.at(iMS1)) * (S1 - mS1.at(iMS1) + 1))*
				  ClebschGordan(m_wigner, LB,  SB,  JB,  mLB.at(iMLB),   mSB.at(iMSB),   mJB.at(iMJB))*
				  ClebschGordan(m_wigner, LA,  SA,  JA,  mLA.at(iMLA),   mSA.at(iMSA),   0.5/*mJA.at(iMJA)*/)*
				  ClebschGordan(m_wigner, LlA, LrA, LA,  mLlA.at(iMLlA), mLrA.at(iMLrA), mLA.at(iMLA))*
				  ClebschGordan(m_wigner, LlB, LrB, LB,  mLlB.at(iMLlB), mLrB.at(iMLrB), mLB.at(iMLB))*
				  ClebschGordan(m_wigner, SlB, S3,  SB,  mSlB.at(iMSlB), mS3.at(iMS3),   mSB.at(iMSB))*
				  ClebschGordan(m_wigner, S1,  S2,  SlB, mS1.at(iMS1)-1, mS2.at(iMS2),   mSlB.at(iMSlB))*
				  ClebschGordan(m_wigner, SlA, S3,  SA,  mSlA.at(iMSlA), mS3.at(iMS3),   mSA.at(iMSA))*
				  ClebschGordan(m_wigner, S1,  S2,  SlA, mS1.at(iMS1),   mS2.at(iMS2),   mSlA.at(iMSlA));
				AMP1_1+=dummy;
			      }
    AMP1_1 *= flavor_vector.at(2) * (2.*std::sqrt(pi_val * k_value));

    for(int iMSA = 0; iMSA<(int)mSA.size(); iMSA++)
      for(int iMLA = 0; iMLA<(int)mLA.size(); iMLA++)
	for(int iMLlA = 0; iMLlA<(int)mLlA.size(); iMLlA++)
	  for(int iMLrA = 0; iMLrA<(int)mLrA.size(); iMLrA++)                         
	    for(int iMJB = 0; iMJB<(int)mJB.size(); iMJB++)
	      for(int iMSB = 0; iMSB<(int)mSB.size(); iMSB++)
		for(int iMLB = 0; iMLB<(int)mLB.size(); iMLB++)
		  for(int iMLlB = 0; iMLlB<(int)mLlB.size(); iMLlB++)
		    for(int iMLrB = 0; iMLrB<(int)mLrB.size(); iMLrB++){
		      dummy = T1_rho_lambda(k_value, alpha_rho, alpha_lam, mLrA.at(iMLrA), mLlA.at(iMLlA), excMode)*
			KroneckerDelta_extended(mLrA.at(iMLrA), mLlA.at(iMLlA), excMode)*
			KroneckerDelta(mSB.at(iMSB), mSA.at(iMSA)) *
			ClebschGordan(m_wigner, LB,  SB,  JB, mLB.at(iMLB),   mSB.at(iMSB),    mJB.at(iMJB))*
			ClebschGordan(m_wigner, LlB, LrB, LB, mLlB.at(iMLlB), mLrB.at(iMLrB),  mLB.at(iMLB))*
			ClebschGordan(m_wigner, LlA, LrA, LA, mLlA.at(iMLlA), mLrA.at(iMLrA),  mLA.at(iMLA))*
			ClebschGordan(m_wigner, LA,  SA,  JA, mLA.at(iMLA),   mSA.at(iMSA),    0.5/*mJA.at(iMJA)*/);
		      AMP1_2+=dummy;
		    }
    AMP1_2 *= KroneckerDelta(SlA, SlB) * KroneckerDelta(SA, SB) * (1.*std::sqrt(pi_val / k_value));
    AMP1 = AMP1_1 - AMP1_2;

    //AmpSquared AMP2
    for(int iMSA = 0;  iMSA<(int)mSA.size();  iMSA++)
      for(int iMSlA = 0; iMSlA<(int)mSlA.size(); iMSlA++)
	for(int iMLA = 0; iMLA<(int)mLA.size(); iMLA++)
	  for(int iMLlA = 0; iMLlA<(int)mLlA.size(); iMLlA++)
	    for(int iMLrA = 0; iMLrA<(int)mLrA.size(); iMLrA++)
	      for(int iMSB = 0;  iMSB<(int)mSB.size(); iMSB++)
		for(int iMSlB = 0;  iMSlB<(int)mSlB.size(); iMSlB++)
		  for(int iMLB = 0;  iMLB<(int)mLB.size(); iMLB++)
		    for(int iMLlB = 0; iMLlB<(int)mLlB.size(); iMLlB++)
		      for(int iMLrB = 0; iMLrB<(int)mLrB.size(); iMLrB++)
			for(int iMJB = 0;  iMJB<(int)mJB.size(); iMJB++)
			  for(int iMS1 = 0; iMS1 <(int)mS1.size(); iMS1++)
			    for(int iMS2 = 0; iMS2 <(int)mS2.size(); iMS2++)
			      for(int iMS3 = 0; iMS3 <(int)mS3.size(); iMS3++){
				dummy = U2_rho_lambda(k_value, alpha_rho, alpha_lam, mLrA.at(iMLrA), mLlA.at(iMLlA), excMode)*
				  std::sqrt((S2 + mS2.at(iMS2)) * (S2 - mS2.at(iMS2) + 1))*
				  ClebschGordan(m_wigner, LB,  SB,  JB,   mLB.at(iMLB),    mSB.at(iMSB),     mJB.at(iMJB))*
				  ClebschGordan(m_wigner, LA,  SA,  JA,   mLA.at(iMLA),    mSA.at(iMSA),     0.5/*mJA.at(iMJA)*/)*
				  ClebschGordan(m_wigner, LlA, LrA, LA,   mLlA.at(iMLlA),  mLrA.at(iMLrA),   mLA.at(iMLA))*
				  ClebschGordan(m_wigner, LlB, LrB, LB,   mLlB.at(iMLlB),  mLrB.at(iMLrB),   mLB.at(iMLB))*
				  ClebschGordan(m_wigner, SlB, S3,  SB,   mSlB.at(iMSlB),  mS3.at(iMS3),     mSB.at(iMSB))*
				  ClebschGordan(m_wigner, S1,  S2,  SlB,  mS1.at(iMS1),    mS2.at(iMS2) - 1, mSlB.at(iMSlB))*
				  ClebschGordan(m_wigner, SlA, S3,  SA,   mSlA.at(iMSlA),  mS3.at(iMS3),     mSA.at(iMSA))*
				  ClebschGordan(m_wigner, S1,  S2,  SlA,  mS1.at(iMS1),    mS2.at(iMS2),     mSlA.at(iMSlA));			     
				AMP2_1+=dummy;
			      }
    AMP2_1 *= flavor_vector.at(1) * (2.*std::sqrt(pi_val * k_value));

    for(int iMSA = 0; iMSA<(int)mSA.size(); iMSA++)
      for(int iMLA = 0; iMLA<(int)mLA.size(); iMLA++)
	for(int iMLlA = 0; iMLlA<(int)mLlA.size(); iMLlA++)
	  for(int iMLrA = 0; iMLrA<(int)mLrA.size(); iMLrA++)                    
	    for(int iMJB = 0; iMJB<(int)mJB.size(); iMJB++)
	      for(int iMSB = 0; iMSB<(int)mSB.size(); iMSB++)
		for(int iMLB = 0; iMLB<(int)mLB.size(); iMLB++)
		  for(int iMLlB = 0; iMLlB<(int)mLlB.size(); iMLlB++)
		    for(int iMLrB = 0; iMLrB<(int)mLrB.size(); iMLrB++){
		      dummy = T2_rho_lambda(k_value, alpha_rho, alpha_lam, mLrA.at(iMLrA), mLlA.at(iMLlA), excMode)*
			KroneckerDelta_extended(mLrA.at(iMLrA), mLlA.at(iMLlA), excMode)*
			KroneckerDelta(mSB.at(iMSB), mSA.at(iMSA))*
			ClebschGordan(m_wigner, LB,   SB,  JB,  mLB.at(iMLB),    mSB.at(iMSB),    mJB.at(iMJB))*
			ClebschGordan(m_wigner, LlB,  LrB, LB,  mLlB.at(iMLlB),  mLrB.at(iMLrB),  mLB.at(iMLB))*
			ClebschGordan(m_wigner, LlA,  LrA, LA,  mLlA.at(iMLlA),  mLrA.at(iMLrA),  mLA.at(iMLA))*
			ClebschGordan(m_wigner, LA,   SA,  JA,  mLA.at(iMLA),    mSA.at(iMSA),    0.5/*mJA.at(iMJA)*/);
		      AMP2_2+=dummy;
		    }		
    AMP2_2 *= KroneckerDelta(SA, SB) * KroneckerDelta(SlA, SlB) * (1.*std::sqrt(pi_val / k_value));
    AMP2 = AMP2_1 - AMP2_2;
    
    // AmpSquared AMP3
    for(int iMSA = 0;  iMSA<(int)mSA.size();  iMSA++)
      for(int iMSlA = 0; iMSlA<(int)mSlA.size(); iMSlA++)
	for(int iMLA = 0; iMLA<(int)mLA.size(); iMLA++)
	  for(int iMLlA = 0; iMLlA<(int)mLlA.size(); iMLlA++)
	    for(int iMLrA = 0; iMLrA<(int)mLrA.size(); iMLrA++)
	      for(int iMSB = 0;  iMSB<(int)mSB.size(); iMSB++)
		for(int iMSlB = 0;  iMSlB<(int)mSlB.size(); iMSlB++)
		  for(int iMLB = 0;  iMLB<(int)mLB.size(); iMLB++)
		    for(int iMLlB = 0; iMLlB<(int)mLlB.size(); iMLlB++)
		      for(int iMLrB = 0; iMLrB<(int)mLrB.size(); iMLrB++)
			for(int iMJB = 0;  iMJB<(int)mJB.size(); iMJB++)
			  for(int iMS1 = 0; iMS1 <(int)mS1.size(); iMS1++)
			    for(int iMS2 = 0; iMS2 <(int)mS2.size(); iMS2++)
			      for(int iMS3 = 0; iMS3 <(int)mS3.size(); iMS3++){
				dummy = U3_rho_lambda(k_value, alpha_rho, alpha_lam, mLrA.at(iMLrA), mLlA.at(iMLlA), excMode)*
				  std::sqrt((S3 + mS3.at(iMS3))*(S3 - mS3.at(iMS3) + 1))*
				  ClebschGordan(m_wigner, LB,   SB,   JB,   mLB.at(iMLB),    mSB.at(iMSB),     mJB.at(iMJB))*
				  ClebschGordan(m_wigner, LA,   SA,   JA,   mLA.at(iMLA),    mSA.at(iMSA),     0.5/*mJA.at(iMJA)*/)*
				  ClebschGordan(m_wigner, LlA,  LrA,  LA,   mLlA.at(iMLlA),  mLrA.at(iMLrA),   mLA.at(iMLA))*
				  ClebschGordan(m_wigner, LlB,  LrB,  LB,   mLlB.at(iMLlB),  mLrB.at(iMLrB),   mLB.at(iMLB))*
				  ClebschGordan(m_wigner, SlB,  S3,   SB,   mSlB.at(iMSlB),  mS3.at(iMS3) - 1, mSB.at(iMSB))*
				  ClebschGordan(m_wigner, S1,   S2,   SlB,  mS1.at(iMS1),    mS2.at(iMS2),     mSlB.at(iMSlB))*
				  ClebschGordan(m_wigner, SlA,  S3,   SA,   mSlA.at(iMSlA),  mS3.at(iMS3),     mSA.at(iMSA))*
				  ClebschGordan(m_wigner, S1,   S2,   SlA,  mS1.at(iMS1),    mS2.at(iMS2),     mSlA.at(iMSlA));
				AMP3_1+=dummy;
			      }
    AMP3_1 *= flavor_vector.at(0) * (2.*std::sqrt(pi_val * k_value));

    for(int iMSA = 0; iMSA<(int)mSA.size(); iMSA++)
      for(int iMLA = 0; iMLA<(int)mLA.size(); iMLA++)
	for(int iMLlA = 0; iMLlA<(int)mLlA.size(); iMLlA++)
	  for(int iMLrA = 0; iMLrA<(int)mLrA.size(); iMLrA++)
	    for(int iMJB = 0; iMJB<(int)mJB.size(); iMJB++)
	      for(int iMSB = 0; iMSB<(int)mSB.size(); iMSB++)
		for(int iMLB = 0; iMLB<(int)mLB.size(); iMLB++)
		  for(int iMLlB = 0; iMLlB<(int)mLlB.size(); iMLlB++)
		    for(int iMLrB = 0; iMLrB<(int)mLrB.size(); iMLrB++){
		      dummy = T3_rho_lambda(k_value, alpha_rho, alpha_lam, mLrA.at(iMLrA), mLlA.at(iMLlA), excMode)*
			KroneckerDelta_extended(mLrA.at(iMLrA), mLlA.at(iMLlA), excMode)*
			KroneckerDelta(mSB.at(iMSB), mSA.at(iMSA))*
			ClebschGordan(m_wigner, LB,   SB,  JB, mLB.at(iMLB),   mSB.at(iMSB),   mJB.at(iMJB))*
			ClebschGordan(m_wigner, LlB,  LrB, LB, mLlB.at(iMLlB), mLrB.at(iMLrB), mLB.at(iMLB))*
			ClebschGordan(m_wigner, LlA,  LrA, LA, mLlA.at(iMLlA), mLrA.at(iMLrA), mLA.at(iMLA))*
			ClebschGordan(m_wigner, LA,   SA,  JA, mLA.at(iMLA),   mSA.at(iMSA),   0.5/*mJA.at(iMJA)*/);
		      AMP3_2+=dummy;
		    }
    AMP3_2 *= KroneckerDelta(SlA, SlB) * KroneckerDelta(SA, SB) * (1.*std::sqrt(pi_val / k_value));
    AMP3 = AMP3_1 - AMP3_2;
    // sum quark amplitudes, squared them and get squared the total
    TOT_AMP = AMP1 + AMP2 + AMP3;
    SUM_SQUARED_AMP += TOT_AMP * TOT_AMP;
  }
  return SUM_SQUARED_AMP;
}

double EMDecayWidths::T1_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode){
  double thetak=0., phik=0.;
  if(excMode==0)
    return T1r(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  else if(excMode==1)
    return T1l(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  return 0.;
}

double EMDecayWidths::T2_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode){
  double thetak=0., phik=0.;
  if(excMode==0)
    return T2r(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  else if(excMode==1)
    return T2l(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  return 0.;
}

double EMDecayWidths::T3_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode){
  double thetak=0., phik=0.;
  if(excMode==0)
    return T3r();
  else if(excMode==1)
    return T3l(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  return 0.;
}

double EMDecayWidths::U1_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode){return 1.;}
double EMDecayWidths::U2_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode){return 1.;}
double EMDecayWidths::U3_rho_lambda(double k_value, double alpha_rho, double alpha_lam, int mLrA, int mLlA, int excMode){return 1.;}


int EMDecayWidths::KroneckerDelta_extended(double mLrA, double mLlA, int excMode){
  if(excMode==0)
    return KroneckerDelta(mLrA, 1);
  else if(excMode==1)
    return KroneckerDelta(mLlA, 1);
  return 0;
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
  double value2 = 3. * std::pow(mbottom, 2) / std::pow(alpha_lam * (mbottom + 2. * mlight), 2);
  double value3 = 1./std::pow(alpha_rho, 2);
  double value = std::exp(value1 * (value2 + value3));
  return value;
}

double EMDecayWidths::SPINFLIP_U3_GS_GS(double k_value, double alpha_lam,  double mbottom, double mlight){
  double value1 = (-3.0) * std::pow(k_value, 2) * std::pow(mlight, 2) / (2. * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value = std::exp(value1);
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = std::sqrt(6.) * p_imag * mbottom * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = std::sqrt(3.) * p_imag * mbottom * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U1_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-1.0) * std::sqrt(6.) * p_imag * mbottom * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = p_imag * phik;
  double value = std::sqrt(6.) * p_imag * mbottom * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1l_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (-1.0) * std::sqrt(3.) * p_imag * mbottom * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1l_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-1.0) * std::sqrt(6.) * p_imag * mbottom * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_lam  * (mbottom + 2 * mlight) );
  return value;
}

double EMDecayWidths::SPINFLIP_U3_1l_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = p_imag * phik;
  double value2 = (-3.0) * std::pow(mlight, 2) * std::pow(k_value, 2) / 2 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (-1.0) * std::sqrt(6.) * p_imag * mlight * k_value * std::exp(value1 + value2) * std::sin(thetak)/ (2 * alpha_lam  * (mbottom + 2 * mlight) );
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
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
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
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-1.0) * std::sqrt(2) * p_imag * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1r_m1_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value3 = p_imag * phik;
  double value = (-1.0) * std::sqrt(2) * p_imag * k_value * std::exp(value1 + value2 + value3) * std::sin(thetak)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1r_m0_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value =  p_imag * k_value * std::exp(value1 + value2) * std::cos(thetak)/ (2 * alpha_rho);
  return value;
}

double EMDecayWidths::SPINFLIP_U2_1r_m1m_GS(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
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
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
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
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value3 = p_imag * phik;
  double value = (-1.0) * 0.1875 * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/(std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
} 

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1m_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = 2 * p_imag * phik;
  double value = 3 * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1m_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / 8 * std::pow(alpha_rho, 2);
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-3.0) * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/16 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U1_1l_m1m_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value = (value2 * std::pow(std::sin(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

// U2_1lambda-1lambda
double EMDecayWidths::ORBITALSPLIT_U2_1l_m1_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value = std::exp(value1 + value2) * (1 - value2 * std::pow(std::sin(thetak), 2)) ;
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
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = 2 * p_imag * phik;
  double value = 0.375 * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value3 = (-1.0) * p_imag * phik;
  double value = 3. * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/ 16 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  double value = (((-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / 4 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2))) * std::pow(std::cos(thetak), 2) + 1) * k_value * std::exp(value1 + value2)/ (4 * alpha_rho);
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m0_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value3 = p_imag * phik;
  double value = (-1.0) * 0.1875 * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/(std::pow(alpha_lam * (mbottom + 2. * mlight), 2));
  return value;
} 

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1m_1l_m1(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value3 = 2 * p_imag * phik;
  double value = 3 * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::pow(std::sin(thetak), 2)/ (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1m_1l_m0(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double phik, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
  double value2 = (-3.0) * std::pow(mbottom, 2) * std::pow(k_value, 2) / (8 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  double value3 = (-1.0) * p_imag * phik;
  double value = (-3.0) * std::sqrt(2) * std::pow(mbottom, 2) * std::pow(k_value, 2) * std::exp(value1 + value2 + value3) * std::sin(2 * thetak)/(16 * (std::pow(alpha_lam * (mbottom + 2. * mlight), 2)));
  return value;
}

double EMDecayWidths::ORBITALSPLIT_U2_1l_m1m_1l_m1m(double k_value, double alpha_lam, double alpha_rho,  double mbottom, double mlight, double thetak){
  double value1 = (-1.0) * std::pow(k_value, 2) / (8 * std::pow(alpha_rho, 2));
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
  double value1 = std::pow(k_value, 2)*std::pow(mlight, 2)*(1 + 3.*std::cos(2.*thetak) );
  double value2 = 4.*std::pow(alpha_lam, 2)*std::exp((3.*std::pow(k_value,2)*std::pow(mlight, 2))/(2.*std::pow(alpha_lam, 2)*std::pow(mbottom + 2.*mlight, 2)));
  double value2p = std::pow(mbottom + 2.*mlight, 2);
  value2 = value2 * value2p;
  double value = value1/value2;
  return (-1.0) * std::sqrt(3.) * value;  
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
  double value2 = 4.*std::sqrt(6.)*std::pow(alpha_rho,2)*std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
  double value = value1/value2;
  return value;
}

// U2 1nrho->gs
double EMDecayWidths::ORBITALSPLIT_U2_1nr_m0_GS(double k_value, double alpha_lam, double alpha_rho, double mbottom, double mlight){
  double value1 = std::pow(k_value,2);///
  double value2 = 4.*std::sqrt(6.)*std::pow(alpha_rho,2)*std::exp((std::pow(k_value,2)*(std::pow(alpha_rho,-2) + (3*std::pow(mbottom,2))/(std::pow(alpha_lam,2)*std::pow(mbottom + 2*mlight,2))))/8.);
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
  double value = value1 + std::sqrt(1./3.) * value2 +  value3 +  std::sqrt(2./3.) * value4;
  return p_imag * std::sqrt(1./6.) * alpha_lam * value;
}

//T2l
double EMDecayWidths::T2l(double k_value, double alpha_lam, double alpha_rho,
			  double mbottom, double mlight, double thetak, double phik, double mLlA){
  double value1 = ORBITALSPLIT_U2_1l_m1_1l_m1(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value2 = ORBITALSPLIT_U2_2l_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value3 = SPINFLIP_U2_GS_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value4 = ORBITALSPLIT_U2_1nl_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value = value1 + std::sqrt(1/3) * value2 +  value3 +  std::sqrt(2/3) * value4;
  return p_imag * std::sqrt(1./6.) * alpha_lam * value;
}

//T3l
double EMDecayWidths::T3l(double k_value, double alpha_lam, double alpha_rho,
			  double mbottom, double mlight, double thetak, double phik, double mLlA){
  double value1 = ORBITALSPLIT_U3_1l_m1_1l_m1(k_value, alpha_lam, mbottom, mlight, thetak);
  double value2 = ORBITALSPLIT_U3_2l_m0_GS(k_value, alpha_lam, mbottom, mlight, thetak);
  double value3 = SPINFLIP_U3_GS_GS(k_value, alpha_lam, mbottom, mlight);
  double value4 = ORBITALSPLIT_U3_1nl_m0_GS(k_value, alpha_lam, mbottom, mlight);
  double value = value1 + std::sqrt(1./3.) * value2 +  value3 +  std::sqrt(2./3.) * value4;
  return (-1.0) * p_imag * std::sqrt(2./3.) * alpha_lam * value;
}

//T1r
double EMDecayWidths::T1r(double k_value, double alpha_lam, double alpha_rho,
			  double mbottom, double mlight, double thetak, double phik, double mLrA){
  double value1 = ORBITALSPLIT_U1_1r_m1_1r_m1(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value2 = ORBITALSPLIT_U1_2r_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value3 = SPINFLIP_U1_GS_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value4 = ORBITALSPLIT_U1_1nr_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value = value1 + std::sqrt(1./3.) * value2 +  value3 +  std::sqrt(2./3.) * value4;
  return p_imag * std::sqrt(1./2.) * alpha_rho * value;
}

//T2r
double EMDecayWidths::T2r(double k_value, double alpha_lam, double alpha_rho,
			  double mbottom, double mlight, double thetak, double phik, double mLrA){
  double value1 = ORBITALSPLIT_U2_1r_m1_1r_m1(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value2 = ORBITALSPLIT_U2_2r_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak);
  double value3 = SPINFLIP_U2_GS_GS(k_value, alpha_lam, alpha_rho,  mbottom, mlight);
  double value4 = ORBITALSPLIT_U2_1nr_m0_GS(k_value, alpha_lam, alpha_rho, mbottom, mlight);
  double value = value1 + std::sqrt(1./3.) * value2 +  value3 +  std::sqrt(2./3.) * value4;
  return (-1.0) * p_imag * std::sqrt(1./2.) * alpha_rho * value;
}

//T3r
double EMDecayWidths::T3r(){
  double value = 0;
  return value;
}


double EMDecayWidths::ANGULAR_SUM_SQUARED_GROUND(double alpha_rho, double alpha_lam, double k_value){

  WignerSymbols  *m_wigner = new WignerSymbols();
  double dummy = 0;

  double innerSum = 0.;

  //for(int iMJA = 0; iMJA<(int)mJA.size(); iMJA++)
  for(int iML = 0; iML<(int)mL.size(); iML++)
    for(int iMJB = 0; iMJB<(int)mJB.size(); iMJB++)
      for(int iMSA = 0; iMSA<(int)mSA.size(); iMSA++)
 	for(int iMSlB = 0; iMSlB<(int)mSlB.size(); iMSlB++)
 	  for(int iMSB = 0; iMSB<(int)mSB.size(); iMSB++)
 	    for(int iMS1 = 0; iMS1<(int)mS1.size(); iMS1++)
 	      for(int iMS2 = 0; iMS2<(int)mS2.size(); iMS2++)
 		for(int iMSlA = 0; iMSlA<(int)mSlA.size(); iMSlA++)
 		  for(int iMS3 = 0; iMS3<(int)mS3.size(); iMS3++){
 		    dummy = std::sqrt((S1 + mS1.at(iMS1)) * (S1 - mS1.at(iMS1) + 1)) *
 		      ClebschGordan(m_wigner, L, SB, JB, mL.at(iML), mSB.at(iMSB), mJB.at(iMJB))*
 		      ClebschGordan(m_wigner, L, SA, JA, mL.at(iML), mSA.at(iMSA), 0.5/*mJA.at(iMJA)*/)*
 		      ClebschGordan(m_wigner, SlB, S3, SB, mSlB.at(iMSlB), mS3.at(iMS3), mSB.at(iMSB))*
 		      ClebschGordan(m_wigner, S1, S2, SlB, mS1.at(iMS1) - 1, mS2.at(iMS2), mSlB.at(iMSlB))*
 		      ClebschGordan(m_wigner, SlA, S3, SA, mSlA.at(iMSlA), mS3.at(iMS3), mSA.at(iMSA))*
 		      ClebschGordan(m_wigner, S1, S2, SlA, mS1.at(iMS1), mS2.at(iMS2), mSlA.at(iMSlA));
 		    innerSum+=dummy;
 		  }

  innerSum = 0.;

  //for(int iMJA = 0; iMJA<(int)mJA.size(); iMJA++)
  for(int iML = 0; iML<(int)mL.size(); iML++)
    for(int iMJB = 0; iMJB<(int)mJB.size(); iMJB++)
      for(int iMSA = 0; iMSA<(int)mSA.size(); iMSA++)	
 	for(int iMSlB = 0; iMSlB<(int)mSlB.size(); iMSlB++)
 	  for(int iMSB = 0; iMSB<(int)mSB.size(); iMSB++)
 	    for(int iMS1 = 0; iMS1<(int)mS1.size(); iMS1++)
 	      for(int iMS2 = 0; iMS2<(int)mS2.size(); iMS2++)
 		for(int iMSlA = 0; iMSlA<(int)mSlA.size(); iMSlA++)
 		  for(int iMS3 = 0; iMS3<(int)mS3.size(); iMS3++){
 		    dummy = std::sqrt((S2 + mS2.at(iMS2)) *(S2 - mS2.at(iMS2) + 1)) *
 		      ClebschGordan(m_wigner, L, SB, JB, mL.at(iML), mSB.at(iMSB), mJB.at(iMJB))*
 		      ClebschGordan(m_wigner, L, SA, JA, mL.at(iML), mSA.at(iMSA), 0.5/*mJA.at(iMJA)*/)*
 		      ClebschGordan(m_wigner, SlB, S3, SB, mSlB.at(iMSlB), mS3.at(iMS3), mSB.at(iMSB))*
 		      ClebschGordan(m_wigner, S1, S2, SlB, mS1.at(iMS1), mS2.at(iMS2) - 1, mSlB.at(iMSlB))*
 		      ClebschGordan(m_wigner, SlA, S3, SA, mSlA.at(iMSlA), mS3.at(iMS3), mSA.at(iMSA))* 
 		      ClebschGordan(m_wigner, S1, S2, SlA, mS1.at(iMS1), mS2.at(iMS2), mSlA.at(iMSlA));
 		    innerSum+=dummy;
 		  }

 innerSum = 0.;

  //for(int iMJA = 0; iMJA<(int)mJA.size(); iMJA++)
  for(int iML = 0; iML<(int)mL.size(); iML++)
    for(int iMJB = 0; iMJB<(int)mJB.size(); iMJB++)
      for(int iMSA = 0; iMSA<(int)mSA.size(); iMSA++)	
  	for(int iMSlB = 0; iMSlB<(int)mSlB.size(); iMSlB++)
  	  for(int iMSB = 0; iMSB<(int)mSB.size(); iMSB++)
  	    for(int iMS1 = 0; iMS1<(int)mS1.size(); iMS1++)
  	      for(int iMS2 = 0; iMS2<(int)mS2.size(); iMS2++)
  		for(int iMSlA = 0; iMSlA<(int)mSlA.size(); iMSlA++)
  		  for(int iMS3 = 0; iMS3<(int)mS3.size(); iMS3++){
  		    dummy =  std::sqrt((S3 + mS3.at(iMS3)) * (S3 - mS3.at(iMS3) + 1)) *
  		      ClebschGordan(m_wigner, L, SB, JB, mL.at(iML), mSB.at(iMSB), mJB.at(iMJB))*
  		      ClebschGordan(m_wigner, L, SA, JA,  mL.at(iML), mSA.at(iMSA), 0.5/*mJA.at(iMJA)*/)*
  		      ClebschGordan(m_wigner, SlB, S3, SB, mSlB.at(iMSlB), mS3.at(iMS3) - 1, mSB.at(iMSB))*
  		      ClebschGordan(m_wigner, S1, S2, SlB,  mS1.at(iMS1), mS2.at(iMS2), mSlB.at(iMSlB))*
  		      ClebschGordan(m_wigner, SlA, S3, SA,  mSlA.at(iMSlA), mS3.at(iMS3), mSA.at(iMSA))*
  		      ClebschGordan(m_wigner, S1, S2, SlA, mS1.at(iMS1), mS2.at(iMS2), mSlA.at(iMSlA));
  		    innerSum+=dummy;
  		  }


  return innerSum;
}

#endif


  // // test function
  // double thetak = 0; double phik = 0; double mLlA=1;
  // double tensor1 = T1l(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  // double tensor2 = T2l(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  // double tensor3 = T3l(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  // double tensor4 = T1r(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  // double tensor5 = T2r(k_value, alpha_lam, alpha_rho, mbottom, mlight, thetak, phik, mLlA);
  // std::cout<<tensor1<<"    "<<tensor2<<"   "<<tensor3<<"   "<<tensor4<<"   "<<tensor5<<std::endl;
