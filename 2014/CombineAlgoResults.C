//This macro will collect up all of the results in an algorithm results file
//and get the mean + std dev of each sample. It will write these to a new tree.
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include "tinydir.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include "TBranch.h"
#include "TString.h"
#include "TH1D.h"
#include <fstream>
#include "TLeaf.h"
#include <math.h>
#include "TMath.h"
#include "TCanvas.h"

using namespace std;

void* GetPtrToVal(TBranch* theBranch, int entry){
  const char* name = theBranch->GetName();
  theBranch->GetEntry(entry);
  TLeaf* theLeaf = theBranch->GetLeaf(name);
  return theLeaf->GetValuePointer();
}

void CombineAlgoResults(string fname){

  //Open the file.
  TFile* fdata = TFile::Open(fname.c_str());
  TTree* theTree = (TTree*)fdata->Get("Results");
  TBranch* TimeBranch = (TBranch*)theTree->GetBranch("TimeStamp");
  TBranch* FactBranch = (TBranch*)theTree->GetBranch("Factor");
  TBranch* CntTimeBranch = (TBranch*)theTree->GetBranch("CntTime");
  TBranch* RefEndBranch = (TBranch*)theTree->GetBranch("RefEnd");
  TBranch* TotCntsBranch = (TBranch*)theTree->GetBranch("TotCnts");
  TBranch* NameBranch = (TBranch*)theTree->GetBranch("SampleName");
  TBranch* IDnumBranch = (TBranch*)theTree->GetBranch("IDnum");
  TBranch* DoseBranch = (TBranch*)theTree->GetBranch("Dose");

  //Get a vector of the sample names within the file.
  vector< string > vNames;

  Long64_t numEntries = theTree->GetEntries();
  string thisName = "";
  bool isRpt = false;

  for(int i = 0; i<numEntries; i++){
    isRpt = false;
    thisName = *(string*)GetPtrToVal(NameBranch, i);
    for(int j = 0; j<vNames.size(); j++){
      if(thisName == vNames.at(j)){isRpt = true;}
    }
    if(!isRpt){vNames.push_back(thisName);}
  }

  cout << "FILLED THE VECTOR" << endl;

  //Create the file to be written to and the new tree
  TFile* fOut = new TFile("CombinedResults.root", "RECREATE");
  TTree* outTree = new TTree("Results", "Combined rad damage results");
  float theRelChange = 1;
  float theStdDev = 0;
  float theDose = 0;
  string theSampleName = "";
  int theSampleNum = 0;
  outTree->Branch("RelChange", &theRelChange, "RelChange/F");
  outTree->Branch("StdDev", &theStdDev, "StdDev/F");
  outTree->Branch("Dose", &theDose, "Dose/F");
  outTree->Branch("SampleName", "string", &theSampleName);
  outTree->Branch("SampleNum", &theSampleNum, "SampleNum/I");

  //Create a histogram of (1 - Factor/RefEnd) for each sample name.
  //Get the mean/std dev of each.
  TH1F* hRelChange;
  string command;
  TCanvas* theCanvas = new TCanvas();
  Long64_t numRpts;
  TH1F* hDose;
  for(int i=0; i<vNames.size(); i++){
    command = "SampleName==\"" + vNames.at(i) + "\"";
    numRpts = theTree->Draw("(1-Factor/RefEnd)>>hRelChange", command.c_str());
    hRelChange = (TH1F*)gDirectory->Get("hRelChange");
    theRelChange = hRelChange->GetMean();
    theStdDev = hRelChange->GetRMS();
    theSampleName = vNames.at(i);

    theTree->Draw("Dose>>hDose", command.c_str());
    hDose = (TH1F*)gDirectory->Get("hDose");
    theDose = hDose->GetMean();

    outTree->Fill();

  }

  //Write to file
  outTree->Write();
  //Close the files.
  fdata->Close();
  fOut->Close();

}
