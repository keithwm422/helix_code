//This code determines the resolution as function of energy per nucleon of the HELIX experiment components.
// R=p/Z
// uncertainty in p from uncertainty in R and uncertainty in velocity from uncertainty in tof, ck, and rich 
// Process is to pick a charge and a mass (z=4, m=10) and then stepping over energy which will give resolutions in p and beta gamma. 

#include <iostream>
#include <fstream>
#include "TSystem.h"
using namespace std;
double  calculate_beta(double gamma);

void mass_jims_way(){
// chosen constants
  gSystem->SetFPEMask();
  double mass=10.0;  // mass (number of nucleons)
  double charge=4.0; // charge
  double Rm=800.0; // maximum detectable rigidity GV
  double step_size=0.001; // divides up the energy into bins
  double gamma_low=1.0;  // GeV/nucleon
  double gamma_high=10.0; 
  double index_of_refraction_silica=1.16; // of silica aerogel
  double index_of_refraction_water=1.33; // of water
  double norm_1=100.0;  //for the cereknov water counter test
  double res_index_of_refraction_silica=0.01;
  double timing_factor=(100.0/7700.0); // this factor is the delta t from Jim's calculation (100ps) over the time it takes for light in vacuum to traverse from top to bottom TOF (7.7 ns)
  double index_water_squared=index_of_refraction_water*index_of_refraction_water;

//TH2D histogram for delta_m vs kinetic energy
  TProfile *h_delta_m_vs_T = new TProfile("Mass resolution vs kinetic energy per nucleon", "Mass resolution vs kinetic energy per nucelon", 1000, 0 , 9);
// file to save hists
  TFile *file_out = new TFile("mass_jims_way.root","RECREATE");
//  loop for the gamma's to run through (really the energies to run through)
  double gamma = gamma_low;
  while(gamma <= gamma_high){
  //increase gamma
    gamma+=step_size;

    // calculate things needed, like from gamma, mass and charge calculate beta, R, etc. 

    double beta=calculate_beta(gamma);

//    double res_p=beta*gamma*mass/Rm;  //resolution of momentum, p is deltap over p and is equal to delta R over R which is R over MDR (Rm). now for a given mass and gamma, we know delta_p over p or res_p
    double res_p=0;  //turn off momentum resolution
// velocity res needs fixing since at low gamma, doesnt cerenk.
    double delta_beta_over_beta=beta*(timing_factor); // this starts with uncertainty in TOF, given by notes with Jim's help
// need to edit to skip over this resolution if not applicable (beta too small).
    if(beta>1./1.33){
      delta_beta_over_beta=index_water_squared/(2.0*TMath::Sqrt(charge*charge*norm_1))*TMath::Sqrt(1.0-(1.0/index_water_squared))*(beta*beta)*(TMath::Sqrt(1.0-(1.0/(index_water_squared*beta*beta)))); // this huge line has uncertainty in index of refraction for RICH (times gamma squared in the res_pk variable)
    }
    

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
  h_delta_m_vs_T->Write("delta_m_vs_T");

  file_out->Close();
}


double calculate_beta(double gamma){

  return TMath::Sqrt(1.0-(1.0/TMath::Power(gamma,2.0)));
}
