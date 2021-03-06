//This is a code to simulate the resolution of the mass measurement of HELIX instrument
// R=pc/Ze
// quantities are: R, p=beta gamma m , Z, 
// Process is to generate from distrubition R's and betas to then find m's. Arguments for generating the R's and beta's will be the resolution of the detector somehow. 

#include <iostream>
#include <fstream>

using namespace std;
double calculate_gamma(double beta);
//void  generate_cosmic_rays();

void mass_resolution_cardwell(){
// stopwatch for helping with higher statistics and debugging
  TStopwatch timer;
  timer.Start();
// random numbers for sampling
  TRandom * rand = new TRandom(0);
// these randoms will be used for sampling from distr. the observed quantities given an actual quantity (gaussian distr with error bars related to sigma).
//  TRandom3 *fRand = new TRandom3(seed);

//  fRand->Gaus(0,sigma);

// number of events
  int events=1000;
// declare histograms of true CRs
//  TH1D *hRigid_true= new TH1D("Rigidity filled Histogram of True CRs","Rigidity filled Histogram" , 100,-10,10);
//  TH1D *hBeta_true= new TH1D("Beta filled Histogram of True CRs","Beta filled Histogram", 100, -1, 1);
  TH1D *hMass_true= new TH1D("Mass Identification Histogram of True CRs", "Mass Identification Histogram", 20, 0, 20);
//  TGraph *g_Rigid_Beta_true;
// Open a file to save the data
  TFile *file_out = new TFile("mass_resolution.root","RECREATE");

// declare histograms of observed quantities
//  TH1D *hRigid_observed= new TH1D("Rigidity filled Histogram","Rigidity filled Histogram" , 100,-10,10);
//  TH1D *hBeta_observed= new TH1D("Beta filled Histogram","Beta filled Histogram", 100, -1, 1);
  TH1D *hMass_observed= new TH1D("Mass Identification Histogram", "Mass Identification Histogram", 20, 0, 20);
//  TGraph *g_Rigid_Beta_observed;
// energy range of HELIX is [1GeV/n, 5GeV/n] so what do I need to generate a real physical CR. Flux of Carbon is from pdg or other resource
  // first an energy in this range
  // then, a factor for nuclei number (mass) ratios of abundances in CR or in others like universe/etc.
  // then a charge for that nuclei, a proton number, BeZ=4
  // for R, need Beta perp but since energy and mass are determined, gamma and beta squared are determined, so if beta perp is also constrained then this is good. 

//  generate_cosmic_rays();


// Now, to plot the HELIX observations, need to sample from a distrubition that reflects HELIX's expected rigidity and beta resolution, to see the mass resolution. 

// declare variables of physical quantities, including arrays, that are measured (have to make them physical). In other words, make a real nuclei, mass, charge, energy, then sample from a dsitrubition of possible values measured. Pick a composition is the ANSWER!
//  double rigid=0.0;
//  double rigid_A[events];
//  double beta=0.0;
//  double beta_A[events];
// declare variables of quantities that will be found
  double mass=9.0;
//  double mass_A[events];

//  R=sample_Rigid(R, seed);
  for (int i=0;i<events;i++){
    // Sample variables
    // per each event, find gamma
//    double gamma=calculate_gamma(beta);

    // fill hists true
    hMass_observed->Fill(mass);
    hMass_true->Fill(mass);
    // fill arrays
//    rigid_A[i]=rigid;
//    beta_A[i]=beta;

  }


// save graphs to root file:
  file_out->cd();
//  hRigid->Write("Rigidity");
//  hBeta->Write("Beta");
  hMass_observed->Write("Massobserved");
  hMass_true->Write("Masstrue");
//  g_Rigid_Beta->Write("Rigid_Beta");

  file_out->Close();
//  return 0;
}

//double calculate_gamma(beta){
//  double gamma =1./sqrt(1-(beta*beta));
//  return gamma;
//}

//void generate_cosmic_rays(){
//  double cosmic_rays_energy[events];
//  

//}

/*double sample_Rigid(R, seed){
  
}
*/
