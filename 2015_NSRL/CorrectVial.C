//This program will:
// - Open a file containing reference-corrected irradiated scintillator results.
// - Open a file containing vial damage results.
// (I assume that I've measured the damaged vials and the scintillators at the
// same time).
// - I will create a map of event numbers that correspond to the vial damage
//   measurements of the appropriate vial, sorted by timestamp.
// - I will loop through every entry in the scintillator results, getting the
//   results for the vial that I want.
//   When the correct sample name pops up, determine the nearest (in time)
//   reference vial measurement. Check that it is good (factor!=DBL_MAX).
//   If it is good, make the correction; repeat until nearest good one is found.

//Note: This correction only needs to be made AFTER the irradiation has occurred
//so I need an additional argument specifying a timestamp after which to apply
//the correction.

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
void CorrectVial(string Scintfname, string Vialfname,
		       string TestSampleName, Long64_t MinTimeStamp){

  if((TestSampleName.find("_REF")!=string::npos)||
     (TestSampleName.find("_Leach")!=string::npos)||
     (TestSampleName.find("STD")!=string::npos)||
     (TestSampleName.find("EMPTY")!=string::npos)||
     (TestSampleName.find("HDPE")!=string::npos)||
     (TestSampleName.find("PP")!=string::npos)    
     ){
    cout << "Skipping file: " << TestSampleName << endl;
    return;
  }

  cout << "Opening " << Scintfname << endl;

  TFile* ScintInFile = TFile::Open(Scintfname.c_str());
  //TTree* ScintTree = (TTree*)ScintInFile->Get("FinalResults");
  TChain* ScintTree = (TChain*)ScintInFile->Get("ChainFinalResults");

  cout << "Opening " << Vialfname << endl;

  TFile* VialInFile = TFile::Open(Vialfname.c_str());
  TChain* VialChain = (TChain*)VialInFile->Get("ChainResults");


  //Get the Scint dmg tree.
  //  vector<double>* theScintChannel = NULL;
  //vector<double>* theScintCounts = NULL;
  Long64_t theScintTime = -1; 
  float theScintCntTime = -1;
  int theScintTotCnts = -1;
  string* theScintName = NULL;
  int theScintIDnum = -2;
  bool theScintisRef = false;
  int theScintPosition = -1;
  int theScintSampleNum = 1;
  double theScintFactor = DBL_MAX;
  double theScintTestEnd = -1;
  double theScintRefEnd = -1;
  float theScintDose = -1;
  int theScintAbsSampleNum = -1;
  string* theScintFormulation = NULL;
  int theScintAbsPos = -1;
  double theScintCorrLY = DBL_MAX;

  //ScintTree->SetBranchAddress("TestCounts", &theScintCounts);
  //ScintTree->SetBranchAddress("Channel", &theScintChannel);
  ScintTree->SetBranchAddress("TimeStamp", &theScintTime);
  ScintTree->SetBranchAddress("CntTime", &theScintCntTime);
  ScintTree->SetBranchAddress("TotCnts", &theScintTotCnts);
  ScintTree->SetBranchAddress("SampleName", &theScintName);
  ScintTree->SetBranchAddress("IDnum", &theScintIDnum);
  ScintTree->SetBranchAddress("isRef", &theScintisRef);
  ScintTree->SetBranchAddress("Position", &theScintPosition);
  ScintTree->SetBranchAddress("SampleNum", &theScintSampleNum);
  ScintTree->SetBranchAddress("Factor", &theScintFactor);
  ScintTree->SetBranchAddress("TestEnd", &theScintTestEnd);
  ScintTree->SetBranchAddress("RefEnd", &theScintRefEnd);
  ScintTree->SetBranchAddress("Dose", &theScintDose);
  ScintTree->SetBranchAddress("Formulation", &theScintFormulation);
  ScintTree->SetBranchAddress("AbsPos", &theScintAbsPos);
  ScintTree->SetBranchAddress("CorrLY", &theScintCorrLY);

  Long64_t Scintnentries = ScintTree->GetEntries();
  double theScintRefOrig;

  cout << "Test sample name = " << TestSampleName << endl;

  if(TestSampleName.find("\r")!=string::npos){
    cout << "TestSample is terminated with \"\r\", trimming this!" << endl;
    TestSampleName.replace(TestSampleName.find("\r"), string::npos, "");
  }
  
  //Get the Vial dmg chain.
  //Parameters I'll want: SampleName, Factor, TimeStamp.
  Long64_t theVialTime = -1; 
  string* theVialName = NULL;
  double theVialFactor = DBL_MAX;
  VialChain->SetBranchAddress("TimeStamp", &theVialTime);
  VialChain->SetBranchAddress("SampleName", &theVialName);
  VialChain->SetBranchAddress("Factor", &theVialFactor);

  Long64_t Vialnentries = VialChain->GetEntries();


  //Create a tree to write to.
  TFile* OutFile = new TFile("VialCorrectedResults.root", "RECREATE");
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

  ///////////////
  //I need to determine here the appropriate name of the vial to be compared to.
  ///////////////

  string VialSampleName = "";

  if((TestSampleName.find("005A_")!=string::npos)||
     (TestSampleName.find("010A_")!=string::npos)||
     (TestSampleName.find("005A_")!=string::npos)){
    //A WbLS sample, check if it is a PMT sample.
    if((TestSampleName.find("PMT")!=string::npos)){
      //A PMT sample
      VialSampleName = "HDPE_005_High";     
    }
    else{
      //Not a PMT sample, check if it is a front sample
      if((TestSampleName.find("_F")!=string::npos)){
	//A front sample
	VialSampleName = "PP_005_High";
      }
      else{
	//A rear sample
	VialSampleName = "PP_100_High";
      }
    }
  }
  else{
    //A Pure LS sample. Check if it is a PMT sample.
    if((TestSampleName.find("PMT")!=string::npos)){
      //A PMT sample
      VialSampleName = "HDPE_100_High";     
    }
    else{
      //Not a PMT sample, check if it is a front sample
      if((TestSampleName.find("_F")!=string::npos)){
	//A front sample
	VialSampleName = "PP_100_High";
      }
      else{
	//A rear sample
	VialSampleName = "PP_100_High";
      }
    }
  }

  map<Long64_t, int> VialEventList;//A map storing the vial events sorted by
  //timestamp. Loop through and fill this...
  for(int i = 0; i<Vialnentries; i++){
    VialChain->GetEntry(i);
    if(*theVialName==(VialSampleName)){
      VialEventList[theVialTime] = i;
    }
  }


  cout<< "Number of reference samples = " << VialEventList.size() << endl;
  cout << "(I looked for " << VialSampleName << ")." << endl;

  //This loop will go through every measurement. If the correct SampleName is
  //found, get all of the relevant info (looking at the map to find the factor
  //for the next measurement) and fill the tree.

  //I can create an iterator here to return the event number of the next 
  //reference event.
  map<Long64_t, int>::iterator RefIt = VialEventList.begin();

  //I assume that the data are stored in the tree sorted by ascending timestamp.
  for(int i = 0; i<Scintnentries; i++){
    ScintTree->GetEntry(i);
    //Do something only if it is the sample we want
    if((*theScintName==TestSampleName)&&(!theScintisRef)&&(theScintIDnum!=4)
       &&(theScintIDnum!=0)){
      //cout << "Getting Sample Entry # " << i << endl;
      //Set the branches, except for the corrected LY.
      OutTime = theScintTime;
      OutCntTime = theScintCntTime;
      OutTotCnts = theScintTotCnts;
      OutName = *theScintName;
      OutIDnum = theScintIDnum;
      OutisRef = theScintisRef;
      OutPosition = theScintPosition;
      OutSampleNum = theScintSampleNum;
      //OutFactor = theScintFactor;
      OutTestEnd = theScintTestEnd;
      OutRefEnd = theScintRefEnd;
      OutDose = theScintDose;
      OutAbsSampleNum = theScintAbsSampleNum;
      OutFormulation = *theScintFormulation;
      OutAbsPos = theScintAbsPos;
      OutLY = theScintCorrLY;

      if(theScintTime>MinTimeStamp){
	
	//if(RefIt!=VialEventList.end()){
	//  cout << "RefIt->first = " << RefIt->first << ", "
	//     << "RefIt->second = " << RefIt->second << endl;
	//}
	
	//Now I need to grab the next event according to the map.
	//If the next event is the end, just use the last event.
	//If the time difference between the reference and the sample is >10 hr
	//I've done something wrong. I'll check the previous event as the
	//most likely cause of this is a run that stopped before the sample
	//was measured. If that doesn't work I'll spit the dummy and throw an
	//error.
	//if(RefIt==VialEventList.end()){RefIt--;}
	//else if(abs(OutTime - RefIt->first)>36000){
	  //check the previous event;
	  //RefIt--;
	  //if(abs(OutTime - RefIt->first)>36000){
	    //Throw an error as something went wrong.
	    //cout << "Error: Couldn't find the correct reference vial for sample "
	//	 << TestSampleName << ", event number " << i << "!" << endl
	//	 << "Sample timestamp = " << OutTime << endl
	//	 << "Previous timestamp = " << RefIt->first << endl;
	//  RefIt++;
	//  if(RefIt!=VialEventList.end()){
	//    cout << "Next timestamp = " << RefIt->first << endl;
	//  }
	//  else{cout << "(that was the last reference measurement)" << endl;}
	//  return;
	//}
	//}
	VialChain->GetEntry(RefIt->second);      
	OutCorrLY = theVialFactor*OutLY;
	OutRefFactor = theVialFactor;
	RefIt++;
      }
      else{
	//Prior to irrad time, don't apply a vial correction.
	OutCorrLY = OutLY;
	OutRefFactor = 0;
      }
      //Fill the tree and increment the iterator.
      OutTree->Fill();
     
    }
  }

  OutTree->Write();
  OutFile->Close();
}
