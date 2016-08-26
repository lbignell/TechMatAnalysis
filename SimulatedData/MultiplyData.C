#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include <string>
#include <iostream>
#include "TCanvas.h"

using namespace std;

void MultiplyData(string fname, double theFactor){
  TFile *f = TFile::Open(fname.c_str());
  TTree* theTree = (TTree*)f->Get("Results");
  double theEdep = 0;
  theTree->SetBranchAddress("Edep", &theEdep);
  Long64_t nentries = theTree->GetEntries();

  TFile* fout = new TFile("TestMultiplyOut.root", "RECREATE");
  TTree* ScaledTree = new TTree("ScaledTree", "OrigData Scaled by a Factor");
  double ScaledEdep = 0;
  ScaledTree->Branch("Edep", &ScaledEdep, "Edep/D");
  cout << "Doing scaling... ";
  for(Long64_t i = 0; i<nentries; i++){
    theTree->GetEntry(i);
    ScaledEdep = theEdep*theFactor;
    ScaledTree->Fill();
  }
  cout << "Done! Drawing result." << endl;
  TCanvas c1;
  theTree->Draw("Edep>>hOrig", "Edep!=0", "");
  ScaledTree->Draw("Edep>>hScaled", "Edep!=0", "SAME");
  ScaledTree->Write();
  fout->Close();
}
