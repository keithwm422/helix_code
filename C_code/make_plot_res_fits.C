// this code loads mass histograms and make plots with fits, specifically for group meeting.
// to-do, output the fit results and show as well. 

#include <iostream>
#include <fstream>

using namespace std;
double  generate_resolution_factor();
double calculate_gamma(double beta);

void make_plot_res_fits(){
// canvas for drawing tests
  ofstream outfile;
  outfile.open("fit_params.txt");
  outfile << "fitres0.025" << " " << "fitres0.1" << endl;
  TCanvas *c1= new TCanvas("c1" , "c1" ,800, 800);
  double sigma1 =0.84932*0.25;
  double sigma2 =0.84932*1;
// 1 is res 0.025
  TFile *file1= TFile::Open("/home/keith/Desktop/HELIX/helix_code/mass_fits_1.root");
// 2 is res 0.1
  TFile *file2= TFile::Open("/home/keith/Desktop/HELIX/helix_code/mass_fits_2.root");
// get histograms
  TH1D* mass1 = 0;
  TH1D* mass2 = 0;
  file1->cd();
  file1->GetObject("MassObserved", mass1);
  file2->cd();
  file2->GetObject("MassObserved", mass2);
//make lines thick and colors different
  mass1->SetLineWidth(2);
  mass2->SetLineWidth(2);
  mass1->SetLineColor(2);
  mass2->SetLineColor(4);
// draw before fitting
  gStyle->SetOptStat(0);
  mass1->Draw();
  c1->Update();
  mass2->Draw("same");
  c1->Update();
// Make fit functions (2 since they will be for different resolutions)
  TF1 *fgauss1 = new TF1("fgauss1","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])",6,14);
  fgauss1->SetParameters(1000,9,sigma1,500,10,sigma1);
  fgauss1->SetParLimits(0,100,1500);
  fgauss1->SetParLimits(1,8,9);
  fgauss1->SetParLimits(3,100,1100);
  fgauss1->SetParLimits(4,9,10);
  fgauss1->SetLineColor(1);
  TF1 *fgauss2 = new TF1("fgauss2","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])",6,14);
  fgauss2->SetParameters(1000,9,sigma2,500,10,sigma2);
  fgauss2->SetParLimits(0,100,1500);
  fgauss2->SetParLimits(1,8,9);
  fgauss2->SetParLimits(3,100,1100);
  fgauss2->SetParLimits(4,9,10);
  fgauss2->SetLineColor(6);
// fit the functions to the histograms
  c1->ResetDrawn();
  mass1->Fit("fgauss1","LB","same",6,14);
  c1->Update();
  c1->Print("fit_1.png");
  mass2->Fit("fgauss2","LB","same",6,14);
  c1->Update();
  TLegend *leg= new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(mass1,"res=0.025");
  leg->AddEntry(mass2,"res=0.1");
  leg->Draw("same");
  c1->Update();
//  TFunction *f1= new TFunction();
// log-likilhood goes here?
  outfile << fgauss1->GetChisquare() << " " << fgauss2->GetChisquare() << endl;

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

// make the functions and parameters to fit histogram with

//  TF1 *fgauss2 = new TF1("fgauss2","[2]*TMath::Gaus(x,[0],[1])",6,14);
//  fgauss2->SetParameters(10,sigma,1000);
//  fgauss2->SetParLimits(0,9.9,10.1);
// do the fits
  c1->Print("fit_1_2.png");
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

  file1->Close();
  file2->Close();
  outfile.close();
}
