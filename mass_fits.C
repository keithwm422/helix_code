//This is a code to simulate the resolution of the mass measurement of HELIX instrument
// R=pc/Ze
// quantities are: R, p=beta gamma m , Z, 
// Process is to generate cosmic ray masses and then try fits to understand better the fitting necessary to rule out a ratio or hypothesis.

#include <iostream>
#include <fstream>

using namespace std;
double  generate_resolution_factor();
double calculate_gamma(double beta);

void mass_fits(){
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
  TFile *file_out = new TFile("mass_fits.root","RECREATE");

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

// make the functions and parameters to fit histogram with
  TF1 *fgauss1 = new TF1("fgauss1","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])",6,14);

//  TF1 *fgauss2 = new TF1("fgauss2","[2]*TMath::Gaus(x,[0],[1])",6,14);
  fgauss1->SetParameters(1000,9,sigma,5000,10,sigma);
  fgauss1->SetParLimits(0,500,1500);
  fgauss1->SetParLimits(1,8,9);
  fgauss1->SetParLimits(3,100,1100);
  fgauss1->SetParLimits(4,9,10);
//  fgauss2->SetParameters(10,sigma,1000);
//  fgauss2->SetParLimits(0,9.9,10.1);
// do the fits
  hMass_observed->Fit("fgauss1","LB","same",6,14);
  c1->Update();
  c1->Print("fit_gauss.png");
//  hMass_observed->Fit("fgauss2","LB","same",6,14);
//  c1->Update();
//std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << std::endl;

//  c1->Print("fit_gauss10.png");

  cout << "declared out func 1" << endl;
  fgauss1->Draw();
  cout << "drawing out func 1" << endl;
  c1->Update();
  c1->Print("fit_gauss_function1.png");



//  timer.Stop();
//  cout << "took ";
//  timer.Print();
//  cout  << " seconds" << endl;

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
  return 1;
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
