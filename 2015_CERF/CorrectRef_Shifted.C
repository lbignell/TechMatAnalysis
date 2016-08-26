//This program will:
// - Loop through the REF vials and get Ref_orig/Sample_orig vals.
// - Loop through every entry name and for each TestSampleName entry:
//   - Normalise to SampleOrig.
//   - Determine the nearest (in time) reference measurement. Get Ref_this.
//   - Correct for (Ref_Orig/Ref_this);
#include <string>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TBranch.h"
#include "TString.h"
#include "TH1D.h"
#include <fstream>
#include <iostream>
#include "TLeaf.h"
#include <math.h>
#include "TMath.h"
#include <sstream>

using namespace std;

//The argument is the file containing a TChain of shifting results.
void CorrectRef(string fname, string TestSampleName){

  TFile* InFile = TFile::Open(fname.c_str());
  TChain* InChain = (TChain*)InFile->Get("ChainResults");

  vector<double>* theChannel = NULL;
  vector<double>* theCounts = NULL;
  Long64_t theTime = -1; 
  float theCntTime = -1;
  int theTotCnts = -1;
  string* theName = NULL;
  int theIDnum = -2;
  bool theisRef = false;
  int thePosition = -1;
  int theSampleNum = 1;
  double theFactor = DBL_MAX;
  double theTestEnd = -1;
  double theRefEnd = -1;
  float theDose = -1;
  int theAbsSampleNum = -1;
  string* theFormulation = NULL;
  int theAbsPos = -1;

  InChain->SetBranchAddress("TestCounts", &theCounts);
  InChain->SetBranchAddress("Channel", &theChannel);
  InChain->SetBranchAddress("TimeStamp", &theTime);
  InChain->SetBranchAddress("CntTime", &theCntTime);
  InChain->SetBranchAddress("TotCnts", &theTotCnts);
  InChain->SetBranchAddress("SampleName", &theName);
  InChain->SetBranchAddress("IDnum", &theIDnum);
  InChain->SetBranchAddress("isRef", &theisRef);
  InChain->SetBranchAddress("Position", &thePosition);
  InChain->SetBranchAddress("SampleNum", &theSampleNum);
  InChain->SetBranchAddress("Factor", &theFactor);
  InChain->SetBranchAddress("TestEnd", &theTestEnd);
  InChain->SetBranchAddress("RefEnd", &theRefEnd);
  InChain->SetBranchAddress("Dose", &theDose);
  InChain->SetBranchAddress("Formulation", &theFormulation);
  InChain->SetBranchAddress("AbsPos", &theAbsPos);

  Long64_t nentries = InChain->GetEntries();
  double theRefOrig;

  cout << "Test sample name = " << TestSampleName << endl;
  
  //Create a tree to write to.
  TFile* OutFile = new TFile("CorrectedShiftedResults.root", "RECREATE");
  TTree* OutTree =
    new TTree("FinalResults",
	      "The shifted results, corrected for instrument drift");

  string OutName = "";
  Long64_t OutTime = 0;
  //vector<double> OutChannel(4096,0);
  //vector<double> OutCounts(4096, 0);
  //  int OutPosition = 0;
  //int OutSampleNum = 0;
  float OutCntTime = 0;
  float OutHnum = 0;
  int OutTotCnts = 0;
  bool OutisRef = false;
  double OutFactor = 0;
  double OutTestEnd = 0;
  double OutRefEnd = 0;
  double OutLY = 0;
  double OutCorrLY = 0;
  //OutTree->Branch("Channel", "vector<double>", &OutChannel);
  //OutTree->Branch("Counts",
  OutTree->Branch("TimeStamp", &OutTime, "TimeStamp/L");
  OutTree->Branch("CntTime", &OutCntTime, "CntTime/F");
  OutTree->Branch("Hnum", &OutHnum, "Hnum/F");
  OutTree->Branch("TotCnts", &OutTotCnts, "TotCnts/F");
  OutTree->Branch("SampleName", "string", &OutName);
  OutTree->Branch("IDnum", &OutIDnum, "IDnum/I");
  OutTree->Branch("isRef", &OutisRef, "isRef/O");
  OutTree->Branch("Factor", &OutFactor, "Factor/D");
  OutTree->Branch("TestEnd", &OutTestEnd, "TestEnd/D");
  OutTree->Branch("RefEnd", &OutRefEnd, "RefEnd/D");
  OutTree->Branch("LY", &OutLY, "LY/D");
  OutTree->Branch("CorrLY", &OutCorrLY, "CorrLY/D");


  Long64_t PreRefMeasTime = 0;
  Long64_t PostRefMeasTime = LONG_MAX;
  double PreRefLY = 0;
  double PostRefLY = 0;
  for(int i = 0; i<nentries; i++){
    InChain->GetEntry(i);
    //See if it is a ref for preref purposes...
    if((theisRef)&&(theIDnum==OutIDnum)){
      //Found a ref; store it as a PreRefMeas.
      PreRefMeasTime = theTime;
      PreRefLY = (1 - theFactor/theTestEnd);
    }
    else if((*theName)==TestSampleName){
      //Find the nearest (future) reference.
      int offset = 1;
      while(i + offset<nentries){
	InChain->GetEntry(i+offset);
	if((theisRef)&&(theIDnum==OutIDnum)){
	  //Found next ref; store it!
	  PostRefMeasTime = theTime;
	  PostRefLY = (1 - theFactor/theTestEnd);
	}
	offset++;
      }

      //Back to our entry
      InChain->GetEntry(i);

      //Determine the smaller time difference
      double ClosestRefLY;
      if((theTime - PreRefMeasTime)<(PostRefMeasTime - theTime)){
	//the previous one came closest.
	ClosestRefLY = PreRefLY;
      }
      else{ClosestRefLY = PostRefLY;}

      //Normalise to the first measurement, then normalise to the instrument.
      OutLY = (1 - theFactor/theTestEnd);
      OutCorrLY = OutLY/RefTestSampleOrig;
      OutCorrLY = OutCorrLY*(theRefOrig/ClosestRefLY);

      //Set the branches.
      OutTime = theTime;
      OutCntTime = theCntTime;
      OutHnum = theHnum;
      OutTotCnts = theTotCnts;
      OutisRef = theisRef;
      OutFactor = theFactor;
      OutTestEnd = theTestEnd;
      OutRefEnd = theRefEnd;
      OutTree->Fill();
    }
    //PostRefMeasTime = LONG_MAX;
  }
  OutTree->Write();
  OutFile->Close();
}
