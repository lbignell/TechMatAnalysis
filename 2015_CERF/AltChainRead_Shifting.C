//This function just generates a file with a list of unique sample names, in
//a file called samplesfname.
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

void MakeListOfSamples(string fname, string samplesfname){

  //Open the file.
  TFile* fdata = TFile::Open(fname.c_str());
  TTree* theTree = (TTree*)fdata->Get("LSSpec");

  Long64_t nentries = theTree->GetEntries();
  string* theName = NULL;
  Long64_t theTime = 0;
  vector<double>* theChannel = NULL;
  vector<double>* theCounts = NULL;
  int thePosition = 0;
  int theSampleNum = 0;
  float theCntTime = 0;
  float theHnum = 0;
  int theTotCnts = 0;
  int theIDnum = 0;
  bool theisRef = false;

  cout << "Setting Branch Addresses...";
  theTree->SetBranchAddress("SampleName", &theName);
  theTree->SetBranchAddress("TimeStamp", &theTime);
  theTree->SetBranchAddress("Channel", &theChannel);
  theTree->SetBranchAddress("Counts", &theCounts);
  theTree->SetBranchAddress("Position", &thePosition);
  theTree->SetBranchAddress("SampleNum", &theSampleNum);
  theTree->SetBranchAddress("CntTime", &theCntTime);
  theTree->SetBranchAddress("Hnum", &theHnum);
  theTree->SetBranchAddress("TotCnts", &theTotCnts);
  theTree->SetBranchAddress("IDnum", &theIDnum);
  theTree->SetBranchAddress("isRef", &theisRef);
  cout << "Done!" << endl;

  TFile* OutFile = new TFile("EverySample_Shifted.root", "RECREATE");
  TTree* OutTree = new TTree("Results", "The shifted results; all samples");
  string OutName = "";
  Long64_t OutTime = 0;
  int OutAbsMeasNum = 0;
  vector<double> OutChannel(4096, 0);
  vector<double> OutCounts(4096, 0);
  int OutPosition = 0;
  int OutSampleNum = 0;
  float OutCntTime = 0;
  float OutHnum = 0;
  int OutTotCnts = 0;
  int OutIDnum = 0;
  bool OutisRef = false;

  cout << "Creating OutFile branches";
  //Create all of the branches.
  OutTree->Branch("SampleName", "string", &OutName);
  OutTree->Branch("TimeStamp", &OutTime, "TimeStamp/L");
  OutTree->Branch("AbsMeasNum", &OutAbsMeasNum, "AbsMeasNum/I");
  OutTree->Branch("Channel", "vector<double>", &OutChannel);
  OutTree->Branch("Counts", "vector<double>", &OutCounts);
  OutTree->Branch("Position", &OutPosition, "Position/I");
  OutTree->Branch("SampleNum", &OutSampleNum, "SampleNum/I");
  OutTree->Branch("CntTime", &OutCntTime, "CntTime/F");
  OutTree->Branch("Hnum", &OutHnum, "Hnum/F");
  OutTree->Branch("TotCnts", &OutTotCnts, "TotCnts/I");
  OutTree->Branch("IDnum", &OutIDnum, "IDnum/F");
  OutTree->Branch("isRef", &OutisRef, "isRef/O");
  cout << "Done!" << endl;

  //Get a vector of the sample names within the file.
  vector< string > vNames;

  Long64_t numEntries = theTree->GetEntries();
  string thisName = "";
  bool isRpt = false;

  ofstream ofs("SampleNames.txt", ios::out|ios::trunc);

  //cout << "Entering loop..." << endl;
  for(Long64_t i = 0; i<nentries; i++){
    //cout << "Getting entry " << i << endl;
    theTree->GetEntry(i);
    OutAbsMeasNum++;
    //cout << "*theName = " << *theName << endl;
    //cout << "theChannel = " << theChannel << ", *theChannel.size() = "
    //	 << (*theChannel).size() << endl;
    //cout << "theTime = " << theTime << endl;
    //cout << "thePosition = " << thePosition << endl;
    //cout << "theHnum = " << theHnum << endl;
    //cout << "theIDnum = " << theIDnum << endl;
    //cout << "theisRef = " << theisRef << endl;
    if(theisRef){cout << "theisRef = " << theisRef << endl;}

    //    for(int j = 0; j<10; j++){
    //cout << "(*theCounts).at(" << j << ") = " << (*theCounts).at(j)
    //	   << endl;
    //}
    //cout << "Channel...";
    OutChannel.swap(*theChannel);
    //cout << "counts...";
    OutCounts.swap(*theCounts);
    OutPosition = thePosition;
    OutSampleNum = theSampleNum;
    OutCntTime = theCntTime;
    OutHnum = theHnum;
    OutTotCnts = theTotCnts;
    OutIDnum = theIDnum;
    OutisRef = theisRef;
    OutTime = theTime;
    OutName = *theName;
    //if(theisRef){cout << "theOutRef = " << theOutRef << endl;}

    isRpt = false;

    if(vNames.size()!=0){
      for(int i = 0; i<vNames.size(); i++){
	if((*theName)==vNames.at(i)){isRpt = true;}
      }
      if(!isRpt){
	ofs << *theName;
	ofs << "\n";
	vNames.push_back(*theName);
      }
    }
    else{
      ofs << *theName;
      ofs << "\n";
      vNames.push_back(*theName);
    }
    
    //Fill the tree, if the timestamp is OK.
    //if(theTime>0){
      OutTree->Fill();
      //}

  }

  ofs.close();
  OutTree->Write();
  fdata->Close();
  OutFile->Close();

}
