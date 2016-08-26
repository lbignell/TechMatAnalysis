//This program will:
// - Create a map of the event numbers that correspond to reference vial 
//   measurements, sorted by timestamp.
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
#include <map>

using namespace std;

//The argument is the file containing a TChain of shifting results.
void CorrectRef_Scaled(string fname, string TestSampleName){

  if((TestSampleName.find("_REF")!=string::npos)||
     (TestSampleName.find("_Leach")!=string::npos)||
     (TestSampleName.find("PMT")!=string::npos)||
     (TestSampleName.find("STD")!=string::npos)||
     (TestSampleName.find("EMPTY")!=string::npos)){
    cout << "Skipping file: " << TestSampleName << endl;
    return;
  }

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

  if(TestSampleName.find("\r")!=string::npos){
    cout << "TestSample is terminated with \"\r\", trimming this!" << endl;
    TestSampleName.replace(TestSampleName.find("\r"), string::npos, "");
  }
  
  //Create a tree to write to.
  TFile* OutFile = new TFile("CorrectedShiftedResults.root", "RECREATE");
  TTree* OutTree =
    new TTree("FinalResults",
	      "The shifted results, corrected for instrument drift");

  Long64_t OutTime = -1;
  float OutCntTime = -1;
  int OutTotCnts = -1;
  string OutName = "";
  int OutIDnum = -2;
  bool OutisRef = false;
  int OutPosition = -1;
  int OutSampleNum = -1;
  double OutFactor = DBL_MAX;
  double OutTestEnd = -1;
  double OutRefEnd = -1;
  float OutDose = -1;
  int OutAbsSampleNum = -1;
  string OutFormulation = "";
  int OutAbsPos = -1;
  double OutLY = -1;
  double OutCorrLY = -1;
  double OutRefFactor = -1;

  OutTree->Branch("TimeStamp", &OutTime, "TimeStamp/L");
  OutTree->Branch("CntTime", &OutCntTime, "CntTime/F");
  OutTree->Branch("TotCnts", &OutTotCnts, "TotCnts/I");
  OutTree->Branch("SampleName", "string", &OutName);
  OutTree->Branch("IDnum", &OutIDnum, "IDnum/I");
  OutTree->Branch("isRef", &OutisRef, "isRef/O");
  OutTree->Branch("Position", &OutPosition, "Position/I");
  OutTree->Branch("SampleNum", &OutSampleNum, "SampleNum/I");
  OutTree->Branch("Factor", &OutFactor, "Factor/D");
  OutTree->Branch("TestEnd", &OutTestEnd, "TestEnd/D");
  OutTree->Branch("RefEnd", &OutRefEnd, "RefEnd/D");
  OutTree->Branch("Dose", &OutDose, "Dose/F");
  OutTree->Branch("Formulation", "string", &OutFormulation);
  OutTree->Branch("AbsPos", &OutAbsPos, "AbsPos/I");
  OutTree->Branch("LY", &OutLY, "LY/D");
  OutTree->Branch("CorrLY", &OutCorrLY, "CorrLY/D");
  OutTree->Branch("RefFactor", &OutRefFactor, "RefFactor/D");

  Long64_t PreRefMeasTime = 0;
  Long64_t PostRefMeasTime = LONG_MAX;
  double PreRefLY = 0;
  double PostRefLY = 0;

  map<Long64_t, int> RefEventList;//A map storing the reference events sorted by
  //timestamp. Loop through and fill this...
  for(int i = 0; i<nentries; i++){
    InChain->GetEntry(i);
    if(*theName==(TestSampleName + "_REF")){
      RefEventList[theTime] = i;
    }
  }

  cout<< "Number of reference samples = " << RefEventList.size() << endl;
  cout << "(I looked for " << (TestSampleName + "_REF") << ")." << endl;

  //This loop will go through every measurement. If the correct SampleName is
  //found, get all of the relevant info (looking at the map to find the factor
  //for the next measurement) and fill the tree.

  //I can create an iterator here to return the event number of the next 
  //reference event.
  map<Long64_t, int>::iterator RefIt = RefEventList.begin();

  //I assume that the data are stored in the tree sorted by ascending timestamp.
  for(int i = 0; i<nentries; i++){
    InChain->GetEntry(i);
    //Do something only if it is the sample we want
    if((*theName==TestSampleName)&&(!theisRef)&&(theIDnum!=4)&&(theIDnum!=0)){
      //cout << "Getting Sample Entry # " << i << endl;
      //Set the branches, except for the corrected LY.
      OutTime = theTime;
      OutCntTime = theCntTime;
      OutTotCnts = theTotCnts;
      OutName = *theName;
      OutIDnum = theIDnum;
      OutisRef = theisRef;
      OutPosition = thePosition;
      OutSampleNum = theSampleNum;
      OutFactor = theFactor;
      OutTestEnd = theTestEnd;
      OutRefEnd = theRefEnd;
      OutDose = theDose;
      OutAbsSampleNum = theAbsSampleNum;
      OutFormulation = *theFormulation;
      OutAbsPos = theAbsPos;
      OutLY = 1/(theFactor);

      //cout << "Filled basic info; searching for Ref measurement... " << endl;

      //if(RefIt!=RefEventList.end()){
      //cout << "RefIt->first = " << RefIt->first << ", "
      //     << "RefIt->second = " << RefIt->second << endl;
      //}

      //Now I need to grab the next event according to the map.
      //If the next event is the end, just use the last event.
      //If the time difference between the reference and the sample is more than
      //3600 (an hour), I've done something wrong. I'll check the previous event

      //as the most likely cause of this is a run that stopped before the sample
      //was measured. If that doesn't work I'll spit the dummy and throw an
      //error.
      if(RefIt==RefEventList.end()){RefIt--;}
      else if(abs(OutTime - RefIt->first)>3600){
	//check the previous event;
	RefIt--;
	if(abs(OutTime - RefIt->first)>3600){
	  //Throw an error as something went wrong.
	  cout << "Error: Couldn't find the correct reference vial for sample "
	       << TestSampleName << ", event number " << i << "!" << endl
	       << "Sample timestamp = " << theTime << endl
	       << "Previous timestamp = " << RefIt->first << endl;
	  RefIt++;
	  if(RefIt!=RefEventList.end()){
	    cout << "Next timestamp = " << RefIt->first << endl;
	  }
	  else{cout << "(that was the last reference measurement)" << endl;}
	  return;
	}
      }
      InChain->GetEntry(RefIt->second);      
      OutCorrLY = theFactor/OutFactor;
      OutRefFactor = theFactor;
    
      //Fill the tree and increment the iterator.
      OutTree->Fill();
      RefIt++;
    }
  }

  OutTree->Write();
  OutFile->Close();
}
