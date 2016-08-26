//Define a class with the following methods:
// - A method to convolve a vector of values with a Gaussian (width as argument)
// - A method to multiply every entry in a TTree by a factor, then return the
//   scaled data in a vector.

#ifndef __ProcessSim__
#define __ProcessSim__

#include <vector>
#include "TF1.h"
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"

using namespace std;

class ProcessSim{

public:
  ProcessSim(){};
  ~ProcessSim(){};

  void ConvGaus(vector<double>&, vector<double>&, double, vector<double>&);
  void MultiplySpec(TTree*, double, TH1D*, vector<double>&);

  ClassDef(ProcessSim, 1);
};


void ProcessSim::ConvGaus(vector<double>& SimEnergy, vector<double>& SimData,
			  //void ConvGaus(vector<double>& SimEnergy, vector<double>& SimData,
			  double GausWidth, vector<double>& Result){
  
  vector<double> ConvResult(SimEnergy.size(), 0);
  
  if(SimEnergy.size()!=SimData.size()){
    cout << "Energy and Spectrum vectors aren't the same size!" << endl;
    return;
  }
  
  for(size_t n = 0; n<SimData.size(); n++){
    //Do the convolution
    for(int m = -SimData.size(); m<int(SimData.size());
	m++){
      double GaussVal = (1/(sqrt(2*3.1415)*GausWidth))*
	exp(-0.5*pow(((m)/GausWidth),2));
      if(((n-m)>0)&&((n-m)<SimData.size())){
	ConvResult.at(n) += SimData.at(n-m)*GaussVal;
      }
    }
  }
  
  //Normalise
  double DataSum = 0;
  double ConvSum = 0;
  for(size_t i = 0; i<SimData.size(); i++){
    if(ConvResult.at(i)>0){
      DataSum+=SimData.at(i);
      ConvSum+=ConvResult.at(i);
    }
  }
  
  for(size_t j = 0; j<SimData.size(); j++)
    {ConvResult.at(j) = (DataSum/ConvSum)*ConvResult.at(j);}
  
  Result.swap(ConvResult);
  return;
}

//vector<double>& ProcessSim::MultiplySpec(TTree* theTree, double theFactor,
void ProcessSim::MultiplySpec(TTree* theTree, double theFactor,
			      TH1D* ScaledHist, vector<double>& ScaledVec){
//vector<double>& MultiplySpec(double theFactor){

  double theEdep = 0;
  theTree->SetBranchAddress("Edep", &theEdep);
  Long64_t nentries = theTree->GetEntries();
  
  //TH1D* ScaledHist = new TH1D("ScaledHist", "Scaled Histogram", 2000, 0., 2.);
  
  //cout << "Getting entries from the tree..." << endl;
  for(Long64_t i = 0; i<nentries; i++){
    theTree->GetEntry(i);
    ScaledHist->Fill(theEdep*theFactor);
  }
  //cout << "Done!" << endl;  

  ScaledVec.clear();
  //vector<double> ScaledVec;
  //cout << "Setting ScaledVec..." << endl;
  for(int i = 1; i<2001; i++){
    ScaledVec.push_back(ScaledHist->GetBinContent(i));
  }
  //cout << "Done!" << endl;

  //cout << "&ScaledVec[0] = " << &ScaledVec[0] << endl
  //   << "ScaledHist = " << ScaledHist << endl
  //   << "theTree = " << theTree << endl;


  //ScaledHist->Delete("");
  //delete ScaledHist;

  //cout << "ScaledVec.at(0) = " << ScaledVec.at(0) << endl; 

  return;// ScaledVec;
}

#endif

//This will generate a new tree that will have branches to store the (gaussian
//smeared) spectrum (as vector<double>), the xvals (as vector<double>), and the
//multiplicative factor (as double).
void SynthesiseData(){
  
  cout << "Opening data file.. ";
  TFile *f = TFile::Open("ComptonScatterSim_WithMoreScatter.root");
  cout << "Done!" << endl;
  TTree* InTree = (TTree*)f->Get("Results");

  TH1D* ScaledHist = new TH1D("ScaledHist", "Scaled Histogram", 2000, 0., 2.);

  vector<double> Spectrum;
  vector<double> EnVals;
  double Factor = 1;
  double smearing = 36;//*2.67*1.4142; 
  TFile* fout = new TFile("TestOut.root", "RECREATE");
  TTree* outTree = new TTree("Results", "The results");
  outTree->Branch("Spectrum", "vector<double>", &Spectrum);
  outTree->Branch("Energy", "vector<double>", &EnVals);
  outTree->Branch("Factor", &Factor, "Factor/D");

  ProcessSim* ps = new ProcessSim();

  vector<double> ScaledSpec;

  for(Factor=0.8; Factor<1.101; Factor+=0.02){
  //for(smearing = 66; smearing<110; smearing+=17){  
  //Factor = 1;
    cout << "Factor = " << Factor << ", smearing = " << smearing << endl;
    ScaledSpec.clear();
    //cout << "Running MultiplySpec...";
    ps->MultiplySpec(InTree, Factor, ScaledHist, ScaledSpec);
    ScaledHist->Reset("ICES");
    //cout << "Done!" << endl;
    vector<double> ConvScaledSpec;
    vector<double> Energy;
    for(double i = 0.001; i<2; i+= 0.001){Energy.push_back(i);}
    if(Energy.size()!=ScaledSpec.size()){
      cout << "Warning! Energy size = " << Energy.size() 
	   << " and ScaledSpec size = " << ScaledSpec.size() 
	   << " trying to set them equal..." << endl;
      if(Energy.size()>ScaledSpec.size()){
	while(Energy.size()!=ScaledSpec.size()){
	  Energy.pop_back();
	}
      }
      else{ return;}
    }
    //cout << "Running ConvGaus...";    
    ps->ConvGaus(Energy, ScaledSpec, smearing/sqrt((1/Factor)), ConvScaledSpec);
    //cout << "Done!" << endl;
    Spectrum.swap(ConvScaledSpec);
    EnVals.swap(Energy);
    
    outTree->Fill();
  }
  outTree->Write();
  fout->Close();
  return;
}
