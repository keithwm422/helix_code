//This is a code to simulate the resolution of the mass measurement of HELIX instrument
// R=pc/Ze
// quantities are: R, p=beta gamma m , Z, 
// Process is to generate from distrubition R's and betas to then find m's. Arguments for generating the R's and beta's will be the resolution of the detector somehow. 

#include <iostream>
#include <fstream>

using namespace std;
double  generate_resolution_factor();
double calculate_gamma(double beta);

void Beta_vs_R(){
// stopwatch for helping with higher statistics and debugging
//  TStopwatch timer;
//  timer.Start();
// random numbers for sampling
  TRandom * fRand = new TRandom(0);
// canvas for drawing tests
  TCanvas *c1= new TCanvas("c1" , "c1" ,800, 800);
// these randoms will be used for sampling from distr. the observed quantities given an actual quantity (gaussian distr with error bars related to sigma).
//  TRandom3 *fRand = new TRandom3(seed);


// number of events
  int events_m9=10000;
  int events_m10=events_m9/3;
//  int events_m10=events_m9*(fRand->Uniform(0,1));
// declare histograms of Beta versus Rigidity
//  TH2D *hBeta_vs_R= new TH2D("Events of velocity and rigidity dependence on energy","Rigidity filled Histogram" , 100,-10,10);
  TH1D *hMass_true= new TH1D("Mass Identification Histogram of True CRs", "Mass Identification Histogram", 200, 6, 14);
// Open a file to save the data
  TFile *file_out = new TFile("beta_r.root","RECREATE");

// declare histograms of observed quantities
//  TH1D *hRigid_observed= new TH1D("Rigidity filled Histogram","Rigidity filled Histogram" , 100,-10,10);
//  TH1D *hBeta_observed= new TH1D("Beta filled Histogram","Beta filled Histogram", 100, -1, 1);
  TH1D *hMass_observed= new TH1D("Mass Identification Histogram", "Mass Identification Histogram", 200, 6, 14);
//  TGraph *g_Rigid_Beta_observed;
// energy range of HELIX is [1GeV/n, 5GeV/n] so what do I need to generate a real physical CR. Flux of Carbon is from pdg or other resource
  // first an energy in this range
  // then, a factor for nuclei number (mass) ratios of abundances in CR or in others like universe/etc.
  // then a charge for that nuclei, a proton number, BeZ=4
  // for R, need Beta perp but since energy and mass are determined, gamma and beta squared are determined, so if beta perp is also constrained then this is good. 

//  generate_cosmic_rays();


// Now, to plot the HELIX observations, need to sample from a distrubition that reflects HELIX's expected rigidity and beta resolution, to see the mass resolution. 

// declare variables of physical quantities, including arrays, that are measured (have to make them physical). In other words, make a real nuclei, mass, charge, energy, then sample from a dsitrubition of possible values measured. Pick a composition is the ANSWER!
// declare variables used for sampling the distribution
  double mass=9.0;

// function for calculating deltam from beta and r dependencies.
  double delta_m=generate_resolution_factor();

  double sigma =0.84932*delta_m;// formula for resolution related to sigma of gaussian is sigma=[1/(sqrt(2*ln(2))]*(deltam)
  double sample_mass=0.0;
  cout << "delta_m is " << delta_m << endl;
// fill the mass 9 CRs
  for (int i=0;i<events_m9;i++){
    // Sample variables
    sample_mass=fRand->Gaus(mass,sigma);
    // fill hists true
    hMass_observed->Fill(sample_mass);
    hMass_true->Fill(mass);

  }
// fill the mass 10 CRs
  mass=10.0;
  for (int i=0;i<events_m10;i++){
    // Sample variables
    sample_mass=fRand->Gaus(mass,sigma);
    // fill hists true
    hMass_observed->Fill(sample_mass);
    hMass_true->Fill(mass);

  }
  file_out->cd();
  hMass_observed->Write("MassObserved");
  hMass_true->Write("MassTrue");

// Make fits
//  TFunction *f1= new TFunction();
// log-likilhood goes here?

// save graphs to root file:
//  fgauss->Draw();
// normalise histograms to plot all
/*  cout << "normalising" << endl;
  double norm = 1.0;
  cout << "declaring scale" << endl;
  double integral_n=hMass_observed->Integral();
  double scale = norm/integral_n;
  cout << "scaling" << endl;

  hMass_observed->Scale(scale);
*/
// draw the histogram before doing a fit (required)
  hMass_observed->Draw();
  c1->Update();


  TF1 *fgauss = new TF1("fgauss","[2]*TMath::Gaus(x,[0],[1])",6,14);
  fgauss->SetParameters(9,sigma,1000);
  fgauss->SetParLimits(0,9,9);

  hMass_observed->Fit("fgauss","LB","same",6,14);
  c1->Update();


//  cout << "drawing" << endl;
//  hMass_observed->Draw();
//  hMass_observed->Draw();
//  hMass_observed->Fit(fgauss,"L","same",6,14);
//  hMass_observed->Draw();
//  fgauss->Draw("same");
//  TF1 *fitgauss= hMass_observed->GetFunction("fgauss");
//  cout << fitgauss->GetChisquare() << " is the chi squared value" << endl;

/*
  for (int j=0;j<events_m9+events_m10;j++){
    


  }
*/
//  timer.Stop();
//  cout << "took ";
//  timer.Print();
//  cout  << " seconds" << endl;
  c1->Print("fit_test.png");

  file_out->Close();

}


double generate_resolution_factor(){
  double prefactor=0;
  //generate resolution of rigidity
  double rigidity_frac=0;
  //include here rigidity distributions of particles that would be observed (weighted by flux and such) and resolutions for the detector (depends on rigidity).
  //generate resolution of beta
  double beta=0;
  double beta_frac=0;
  // for chosen beta, calculate gamma
  double gamma = calculate_gamma(beta);
  // calculate prefactor
  prefactor=sqrt((rigidity_frac*rigidity_frac)+(beta_frac*beta_frac)*(gamma*gamma*gamma*gamma));
//  return prefactor;
  return 0.25;
}

double calculate_gamma(double beta){
  double gamma =1./sqrt(1.0-(beta*beta));
  return gamma;
}


//void generate_cosmic_rays(){
//  double cosmic_rays_energy[events];
//  

//}

/*double sample_Rigid(R, seed){
  
}
*/
