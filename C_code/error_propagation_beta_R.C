//This code determines the resolution as function of energy per nucleon of the HELIX experiment components.
// R=p/Z
// uncertainty in p from uncertainty in R and uncertainty in velocity from uncertainty in tof, ck, and rich 
// Process is to pick a charge and a mass (z=4, m=10) and then stepping over energy which will give resolutions in p and beta gamma. 

#include <iostream>
#include <fstream>
#include "TSystem.h"
using namespace std;
double  calculate_beta(double gamma);
//double  calculate_R(double gamma);

const  double mass=10.0;  // mass (number of nucleons)
const  double charge=4.0; // charge
const  double Rm=800.0; // maximum detectable rigidity GV



void error_propagation_beta_R(){
// chosen constants
  gSystem->SetFPEMask();
  double step_size=0.001; // divides up the energy into bins
  double gamma_low=1.0;  // GeV/nucleon
  double gamma_high=1.06; 
  double index_of_refraction_silica=1.16; // of silica aerogel
  double index_of_refraction_water=1.33; // of water
  double norm_1=100.0;  //for the cereknov water counter test
  double res_index_of_refraction_silica=0.01;
  double timing_factor=(100.0/7700.0); // this factor is the delta t from Jim's calculation (100ps) over the time it takes for light in vacuum to traverse from top to bottom TOF (7.7 ns)
  double index_water_squared=index_of_refraction_water*index_of_refraction_water;
  TRandom * fRand = new TRandom(0);
//TH2D histogram for delta_m vs kinetic energy
  TProfile *h_delta_m_vs_T = new TProfile("Mass resolution vs kinetic energy per nucleon", "Mass resolution vs kinetic energy per nucelon", 1000, 0 , 9);
// Histogram of beta and R variables. 
  TH1D *h_beta= new TH1D("velocity histogram", "velocity histogram", 100000, 0, 0.37);
  TH1D *h_beta_log= new TH1D("Log velocity histogram", "Log velocity histogram", 1000, -3, 0);
//  TH1D *h_R= new TH1D("Rigidity histogram", "Rigidity histogram", 100, 0, 1);

// file to save hists
  TFile *file_out0 = new TFile("error_prop_m.root","RECREATE");
  TFile *file_out1 = new TFile("error_prop_beta.root","RECREATE");
//  TFile *file_out2 = new TFile("error_prop_R.root","RECREATE");
//  loop for the gamma's to run through (really the energies to run through)
  double gamma = gamma_low;
  while(gamma <= gamma_high){
  //increase gamma
    gamma+=step_size;

    // calculate things needed, like from gamma, mass and charge calculate beta, R, etc. 

    // we invent a distribution for the velocity measurement with given mean from gamma and standard deviation found through error propagation of the detection method. If the detection method is the TOF, then the error on velocity is from the errors on the length measurement and the timing resolution. The sigma on the velocity is found to be a function of the velocity.
    // So, we sample velocity from a distribution of velocity with mean the true velocity, and sigma related to timing resolution and the velocity mean, true value. 
    double beta_mean=calculate_beta(gamma);
    double beta=beta_mean;

    double delta_beta_over_beta=beta*timing_factor;

    if(beta>1./1.33){
      delta_beta_over_beta=index_water_squared/(2.0*TMath::Sqrt(charge*charge*norm_1))*TMath::Sqrt(1.0-(1.0/index_water_squared))*(beta*beta)*(TMath::Sqrt(1.0-(1.0/(index_water_squared*beta*beta)))); // this huge line has uncertainty in index of refraction for RICH (times gamma squared in the res_pk variable)
    }

    double beta_sigma=beta_mean*beta_mean*(timing_factor); // this starts with uncertainty in TOF, given by notes with Jim's help
//    double R_mean=calculate_R(gamma);
//    double R_sigma=;

    for(int j=0; j<1000; j++){
      double beta_measure=fRand->Gaus(beta_mean,beta_sigma);
      h_beta->Fill(beta_measure);
// fill a logged histogram
      h_beta_log->Fill(TMath::Log(beta_measure));
    }
//    double res_p=beta*gamma*mass/Rm;  //resolution of momentum, p is deltap over p and is equal to delta R over R which is R over MDR (Rm). now for a given mass and gamma, we know delta_p over p or res_p
    double res_p=0;  //turn off momentum resolution
// velocity res needs fixing since at low gamma, doesnt cerenk.
// need to edit to skip over this resolution if not applicable (beta too small).
    double res_pk=gamma*gamma*(delta_beta_over_beta); // deltabeta over beta from cerenkov counters given from Jim: sqrt(1.-(1./n^2)) *n*beta*sqrt((n*beta)^2-1)  / (2.0*charge sqrt(N1))
    //  double res_velocity=res_pk + res_rich;   //resolution of velocity from three different detectors: 1)TOF (function of beta * timing 100ps/charge) 2) Cerenkov Counter (charge*charge*(N1)*[1-1/betasquared nsquared]) 3)RICH = constant
    //  double res_velocity=res_tof + res_pk + res_rich;   //resolution of velocity from three different detectors: 1)TOF (function of beta * timing 100ps/charge) 2) Cerenkov Counter (charge*charge*(N1)*[1-1/betasquared nsquared]) 3)RICH = constant
    //leaving out TOF since I don't understand it yet.
    //  double res_tof=func_beta*100.0/charge; // 100 picoseconds and func_beta is for time of flight processes this needs to be delta betagamma over beta gamma
    //  double res_rich=gamma*gamma*(0.01); // from deltan over n which we know n to within 1 percent

  // calculate delta_m
    // delta m / m squared = delta p over p squared + delta betagamma / betagamme squared.
    // delta p over p is Z*R/Rm where Rm is from the proposal 800 GV and z=4 so delta p over p is done
    // delta betagamma/betagamma is from 3 different things so each of these has a resolution.

    double res_mass=TMath::Sqrt(TMath::Power(res_p,2.0)+TMath::Power(res_pk,2.0));;

    double T=gamma-1.0;

  // fill histogram with calculation (kinetic energy per nucleon is (E-m)/m which is gamma-1
//    h_delta_m_vs_T->Fill(T,res_mass*mass);
    h_delta_m_vs_T->Fill(T,res_mass);

  }



// save graphs to root file:
  file_out1->cd();
  h_beta->Write("beta_measured");
  h_beta_log->Write("beta_log");
  file_out1->Close();

//  file_out2->cd();
//  h_R->Write("Rigidity_measured");
//  file_out2->Close();

  file_out0->cd();
  h_delta_m_vs_T->Write("delta_m_vs_T");
  file_out0->Close();
}


double calculate_beta(double gamma){

  return TMath::Sqrt(1.0-(1.0/TMath::Power(gamma,2.0)));
}
//double calculate_R(double beta){
//  double beta=calculate_beta(gamma);
//  return TMath::Sqrt(1.0-(1.0/TMath::Power(gamma,2.0)));
//}
