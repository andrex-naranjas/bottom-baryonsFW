//BottomDecayWidths includes
#ifndef BOTTOMECAYWIDTHS_H
#define BOTTOMDECAYWIDTHS_H

#include <string>
#include <vector>

class BottomDecayWidths{

public:
  BottomDecayWidths();
  virtual ~BottomDecayWidths();
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

  // FOR ALL STATES
  virtual double CBARIN(double alpha_rho, double alpha_lam);
  virtual double CMESON(double alpha_mes);
  
  //GROUND STATES
  virtual double C0_GROUND(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_GROUND(double alpha_lam, double alpha_mes);
  virtual double F0_GROUND(double alpha_rho, double alpha_lam, double k_value);
  virtual double F01_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double BRO2_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_GROUND(double alpha_lam, double alpha_mes);
  virtual double ARO2_GROUND(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double BRO1_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double F0TOT_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_GROUND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);

  // P-WAVE-lambda
  virtual double CBARFIN_PWAVE(double alpha_rho, double alpha_lam);
  virtual double C0_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes);  
  virtual double ARO0_PWAVE(double alpha_rho, double alpha_mes);
  virtual double F0_PWAVE(double alpha_rho, double alpha_lam, double k_value);
  virtual double F01_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double BRO2_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_PWAVE(double alpha_rho, double alpha_mes);
  virtual double ARO2_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double BRO1_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double F0TOT_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double I01B0TOT_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0TOT_PWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  
  //P-WAVE-rho mode
  virtual double CBARFIN_PWAVE_RHO(double alpha_rho, double alpha_lam);
  virtual double C0_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes);  
  virtual double ARO0_PWAVE_RHO(double alpha_lam, double alpha_mes);
  virtual double F01_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double BRO1_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_PWAVE_RHO(double alpha_lam, double alpha_mes);
  virtual double ARO2_PWAVE_RHO(double alhpha_rho, double alpha_lam, double alpha_mes);
  virtual double BRO2_PWAVE_RHO(double alhpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double F0TOT_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double I01B0TOT_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0TOT_PWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  
  //D-wave-lambda
  virtual double CBARFIN_DWAVE(double alpha_rho, double alpha_lam);
  virtual double C0_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_DWAVE(double alpha_rho, double alpha_mes);
  virtual double F0_DWAVE(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double BRO2_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_DWAVE(double alpha_rho, double alpha_mes);
  virtual double ARO2_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double BRO1_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double F0TOT_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_DWAVE(double alpha_rho, double alpha_lam,double alpha_mes,double k_value);
  virtual double I02B0TOT_DWAVE(double alpha_rho, double alpha_lam,double alpha_mes,double k_value);

  //D-wave-rho
  virtual double CBARFIN_DWAVE_RHO(double alpha_rho, double alpha_lam);
  virtual double C0_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_DWAVE_RHO(double alpha_lam, double alpha_mes);
  virtual double F0_DWAVE_RHO(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double BRO2_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_DWAVE_RHO(double alpha_lam, double alpha_mes);
  virtual double ARO2_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double BRO1_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double F0TOT_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0TOT_DWAVE_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  
  //RADIAL-lambda
  virtual double CBARFIN_RADIAL(double alpha_rho, double alpha_lam);
  virtual double C0_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_RADIAL(double alpha_rho, double alpha_mes);
  virtual double F0_RADIAL(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double BRO2_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_RADIAL(double alpha_rho, double alpha_mes);
  virtual double ARO2_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double BRO1_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double F0TOT_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_RADIAL(double alpha_rho, double alpha_lam,double alpha_mes,double k_value);
  virtual double I02B0TOT_RADIAL(double alpha_rho, double alpha_lam,double alpha_mes,double k_value);
  
  //RADIAL-rho
  virtual double CBARFIN_RADIAL_RHO(double alpha_rho, double alpha_lam);
  virtual double C0_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_RADIAL_RHO(double alpha_lam, double alpha_m);
  virtual double F0_RADIAL_RHO(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double BRO2_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_RADIAL_RHO(double alpha_lam, double alpha_mes);
  virtual double ARO2_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double BRO1_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double F0TOT_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0TOT_RADIAL_RHO(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);

  //MIXED
  virtual double CBARFIN_MIXED(double alpha_rho, double alpha_lam);
  virtual double C0_MIXED(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_MIXED(double alpha_lam, double alpha_mes);
  virtual double F0_MIXED(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double BRO2_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);  
  virtual double ARO1_MIXED(double alpha_lam, double alpha_mes);
  virtual double ARO2_MIXED(double alpha_rho, double alpha_lam, double alpha_mes);  
  virtual double BRO1_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);    
  virtual double ARO4_MIXED(double alpha_lam, double alpha_mes);  
  virtual double BRO4_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);    
  virtual double ARO5_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO6_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);  
  virtual double F0TOT_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);  
  virtual double I01B0_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_MIXED(double alpha_rho, double alpha_lam, double alpha_mes);  
  virtual double I01B0TOT_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);  
  virtual double I02B0TOT_MIXED(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);  

  //SECOND DIAGRAM
  //D-WAVE LAMBDA
  virtual double ANGULAR_SUM_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO0_DWAVE_SECOND(double alpha_rho);
  virtual double F0_DWAVE_SECOND(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_DWAVE_SECOND(double alpha_mes, double k_value);
  virtual double BRO1_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_DWAVE_SECOND(double alpha_lam, double alpha_mes);
  virtual double F0TOT_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0TOT_DWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);


  //SECOND DIAGRAM
  //D-WAVE RHO
  virtual double C0_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_DWAVE_RHO_SECOND(double alpha_rho);
  virtual double F0_DWAVE_RHO_SECOND(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_DWAVE_RHO_SECOND(double alpha_mes, double k_value);
  virtual double BRO1_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_DWAVE_RHO_SECOND(double alpha_lam, double alpha_mes);
  virtual double F0TOT_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_DWAVE_RHO_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);

  //SECOND DIAGRAM
  //P-WAVE LAMBDA
  virtual double C0_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_PWAVE_SECOND(double alpha_rho);
  virtual double F0_PWAVE_SECOND(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_PWAVE_SECOND(double alpha_mes, double k_value);
  virtual double BRO1_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_PWAVE_SECOND(double alpha_lam, double alpha_mes);
  virtual double F0TOT_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0TOT_PWAVE_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);

  //SECOND DIAGRAM
  //MIXED
  virtual double C0_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_MIXED_SECOND(double alpha_rho);
  virtual double F0_MIXED_SECOND(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_MIXED_SECOND(double alpha_mes, double k_value);
  virtual double BRO1_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_MIXED_SECOND(double alpha_lam, double alpha_mes);
  virtual double F0TOT_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0TOT_MIXED_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);


  //SECOND DIAGRAAM
  //RADIAL LAMBDA
  virtual double C0_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ARO0_RADIAL_SECOND(double alpha_rho);
  virtual double F0_RADIAL_SECOND(double k_value, double alpha_rho, double alpha_lam);
  virtual double F01_RADIAL_SECOND(double alpha_mes, double k_value);
  virtual double BRO1_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1_RADIAL_SECOND(double alpha_lam, double alpha_mes);
  virtual double F0TOT_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0TOT_RADIAL_SECOND(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);					      
      
};


//to talk to python
extern "C"{
  double bottom_execute(double ma_val, double mb_val, double mc_val, double ga_val, double sa_val,
			double la_val, double ja_val, double sl_val, double al_val, double ar_val,
			int baryon, int excMode, int prodDecay){

    BottomDecayWidths *  m_decays = new BottomDecayWidths();
    return m_decays->execute(ma_val, mb_val, mc_val, ga_val, sa_val,
			     la_val, ja_val, sl_val, al_val, ar_val,
			     baryon, excMode, prodDecay);
  }
}

#endif //> !CHARMDECAYWIDTHS_H
